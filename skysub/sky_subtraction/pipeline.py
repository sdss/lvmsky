"""One exposure: Gaia-clean medians, ML prediction, per-fiber LSFs and the compact SkyML SFrame."""

import importlib
import json
import logging
import multiprocessing as mp
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from skysub.medians_computation.src.lvm_medians.metadata import median_meta
from skysub.medians_computation.src.lvm_medians.stack import fiber_roles
from skysub.sky_decomp import reliability as rel
from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf

from . import decomp_adapter as adapter
from .fiber_lsf import FIT_CHANNELS, FiberFit, fit_fiber, sky_designs
from .io import BUNDLED_CALIBRATION, TELESCOPES, median_spectrum, read_cframe, resolve_gaia, select_fibers
from .ml import (
    E_SCALE,
    NEBULAR_LINES_RED,
    UNMEASURED_VELOCITY_KM_S,
    SkyArm,
    data,
    inference,
    line_templates,
    oh_group_labels,
    protection_velocities,
    serialization,
    sky_arm_residual_correction,
    sky_line_scaling_correction,
)
from .product import MLCoefficients, ProductInputs, _json, build_product, sha256, write_product
from .rendering import (
    FACTOR,
    FAMILIES,
    arm_sky,
    channel_slices,
    exposure_line_masks,
    fwhm_grid,
    render,
    residual_protection,
    sky_parts,
)
from .settings import SETTINGS

# Names also used by skysub.sky_subtraction_block and analysis scripts.
__all__ = ['run', 'render', 'fwhm_grid', 'sha256', '_json', 'FACTOR', 'FAMILIES']

LOG = logging.getLogger(__name__)
REPOSITORY = Path(__file__).resolve().parents[2]
SKY_ARMS = ('SkyE', 'SkyW')
# Modules besides this package whose code determines the product.
CODE_MODULES = (
    'skysub.decompose_parallel',
    'skysub.sky_decomp.fit',
    'skysub.sky_decomp.lsf_reference',
    'skysub.sky_decomp.lsf_spline2d',
    'skysub.sky_decomp.lsf_surface_iterative',
    'skysub.sky_decomp.result_io',
    'skysub.sky_decomp.solar_lsf',
    'skysub.sky_decomp.telluric_corrected_lines',
    'skysub.medians_computation.src.lvm_medians.gaia',
    'skysub.medians_computation.src.lvm_medians.metadata',
    'skysub.medians_computation.src.lvm_medians.stack',
    'mlp_predictor.sky_arm_correction',
    'mlp_predictor.data',
    'mlp_predictor.inference',
    'mlp_predictor.serialization',
    'mlp_predictor.sky_line_scaling',
)
_WORK = None


@dataclass
class FiberWork:
    """Read-only inputs of the fiber workers; module-level so that forked workers inherit them."""

    model: object
    wave: np.ndarray
    flux: np.ndarray
    ivar: np.ndarray
    mask: np.ndarray
    fiberids: np.ndarray
    references: dict
    comparison: object
    designs: dict
    science_mask: np.ndarray
    coef: np.ndarray
    oh_coef: np.ndarray
    arms: list
    survey_only: bool
    settings: object
    lines: object = None


@dataclass
class FiberResult:
    index: int
    sky: np.ndarray
    fit: FiberFit


@dataclass
class Corrections:
    matrices: dict
    arms: list
    science_mask: np.ndarray
    oh_coef: np.ndarray
    residual_info: dict
    oh_info: dict
    lines: object = None


def resolve_paths(input_path, output):
    """Absolute input and output paths; protects the CFrame and the standard DRP SFrame name."""
    input_path = Path(input_path).resolve()
    if output:
        output = Path(output).resolve()
    else:
        name = input_path.name.replace('lvmCFrame-', 'lvmSFrame-skyml-').removesuffix('.gz')
        output = input_path.with_name(name)
    standard = output.name.startswith('lvmSFrame-') and not output.name.startswith('lvmSFrame-skyml-')
    if output == input_path or standard:
        raise ValueError('Output must protect the CFrame and standard DRP SFrame names')
    if output.exists():
        raise FileExistsError(output)
    return input_path, output


def code_provenance():
    files = {*Path(__file__).parent.glob('*.py')}
    files.update(Path(importlib.import_module(name).__file__) for name in CODE_MODULES)
    hashes = {str(p.resolve().relative_to(REPOSITORY)): sha256(p) for p in sorted(files)}
    try:
        git = ['git', '-C', str(REPOSITORY)]
        commit = subprocess.run([*git, 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout
        status = subprocess.run(
            [*git, 'status', '--porcelain', '--untracked-files=no'], capture_output=True, text=True, check=True
        ).stdout
        revision = dict(commit=commit.strip(), tracked_changes=bool(status.strip()))
    except (OSError, subprocess.CalledProcessError):
        revision = None
    return hashes, revision


def load_fiber_references(frame, calibration):
    """Survey LSF of every usable Sci/SkyE/SkyW fiber, keyed by CFrame row."""
    slit = frame.slit
    fingerprint = adapter.wave_fingerprint(frame.wave)
    references = {}
    for i, fiberid in enumerate(slit['fiberid']):
        if int(slit['fibstatus'][i]) == 0 and str(slit['telescope'][i]).strip() in TELESCOPES:
            reference = load_survey_ascii_lsf(calibration, int(fiberid), allow_warning=True)
            if fingerprint != (reference.wave_n, reference.wave_min, reference.wave_max, reference.wave_sha256):
                raise ValueError('CFrame grid differs from the empirical calibration grid')
            references[i] = reference
    return references


def mean_reference(states):
    """Approximation: the equal mean of the selected fibers' profiles stands in for the LSF of their median.

    Revisit when profiles or wavelength coverage vary strongly between contributors.
    """
    result = deepcopy(states[0])
    for ch in 'BRZ':
        values = [s.coefficients[ch] for s in states if s.metrics[ch].get('quality') == 'measured']
        if not values:
            raise ValueError(f'No usable {ch} survey reference in selected fibers')
        result.coefficients[ch] = np.mean(values, axis=0)
        result.metrics[ch] = dict(status='Solved', model='selected_fiber_equal_mean', n_fibers=len(values))
    return result


def build_stack(frame, path, medians, references, selected, settings=SETTINGS):
    """One-row median stack in the corpus format read by decompose_parallel."""
    header, wave = frame.header, frame.wave
    positions = {label: adapter.telescope_coords(header, label) for label in selected}
    if not np.isfinite(list(positions.values())).all():
        raise ValueError('Science and sky telescope coordinates must be finite')
    separations = {label: adapter.separation(positions['Sci'], positions[label]) for label in SKY_ARMS}
    near = min(separations, key=separations.get)
    far = 'SkyW' if near == 'SkyE' else 'SkyE'
    counts = {f'{label}_{kind}': len(selected[label]) for label in selected for kind in ('good', 'used')}
    meta = median_meta(
        header=header,
        path=path,
        input_index=0,
        expnum=int(header['EXPOSURE']),
        positions=positions,
        separations=separations,
        near=near,
        far=far,
        counts=counts,
        sci_percentile=np.nan,
        sky_percentile=np.nan,
        gaia_sigma=None,
        gaia_ratio_threshold=settings.gaia_ratio_threshold,
    )
    hdus = [
        fits.PrimaryHDU(),
        fits.ImageHDU(wave, name='WAVE'),
        fits.BinTableHDU(Table(rows=[meta]), name='META'),
        fits.ImageHDU(fiber_roles(frame.slit['fiberid'], selected, near, far)[None, :], name='FIBER_ROLE'),
    ]
    for label, suffix in (('Sci', 'SCI'), (near, 'SKY_NEAR'), (far, 'SKY_FAR')):
        median = medians[label]
        for name, values in (
            ('FLUX', median.flux),
            ('IVAR', median.ivar),
            ('MASK', (median.count == 0).astype(np.uint8)),
            ('LSF', fwhm_grid(references[label], wave)),
        ):
            hdus.append(fits.ImageHDU(values[None, :], name=f'{name}_{suffix}'))
    return fits.HDUList(hdus), near, far


def arm_slots(failed):
    """Which telescope's fit fills each sky slot of the prediction: a failed arm is replaced by the other."""
    good = [label for label in SKY_ARMS if label not in failed]
    if not good:
        raise ValueError('No usable sky arm')
    return {label: (label if label in good else good[0]) for label in SKY_ARMS}


def predict_ifu(ensemble_path, frame, median_fits):
    """One ML coefficient vector at the Science IFU center; the ensemble is released on return."""
    ensemble = serialization.load_ensemble(ensemble_path)
    model = median_fits.science_model
    if list(model.design_names) != list(ensemble['coef_names']):
        raise ValueError('Decomposition coefficient names/order disagree with trained ensemble')
    if ensemble.get('zodi_correction', 'none') != adapter.zodi_correction:
        raise ValueError('Ensemble zodiacal convention differs from decomposition')
    header, wave = frame.header, frame.wave
    fits_by_label = median_fits.results
    pos = {k: adapter.telescope_coords(header, k) for k in TELESCOPES}
    lsf = {k: fwhm_grid(fits_by_label[k].lsf_state, wave) for k in TELESCOPES}
    slot = arm_slots(median_fits.failed_arms)
    if median_fits.failed_arms:
        LOG.warning('Sky arm %s is unusable (misfit or Moon glare); %s stands in for it in the prediction',
                    ', '.join(median_fits.failed_arms), slot[median_fits.failed_arms[0]])
    date = str(header['OBSTIME']).strip()
    prediction = inference.predict_sky_from_minimal_inputs(
        ensemble,
        obstime_mjd=[Time(date, scale='utc').mjd],
        sci_ra=[pos['Sci'][0]],
        sci_dec=[pos['Sci'][1]],
        sky_e_ra=[pos[slot['SkyE']][0]],
        sky_e_dec=[pos[slot['SkyE']][1]],
        sky_w_ra=[pos[slot['SkyW']][0]],
        sky_w_dec=[pos[slot['SkyW']][1]],
        coef_e=fits_by_label[slot['SkyE']].coef[None, :],
        coef_w=fits_by_label[slot['SkyW']].coef[None, :],
        wave=wave,
        lsf_e=lsf[slot['SkyE']],
        lsf_w=lsf[slot['SkyW']],
        lsf_sci=lsf['Sci'],
        date_obs=date,
        expnum=int(header['EXPOSURE']),
        exposure_seconds=float(header['EXPTIME']),
    )
    coef = np.asarray(prediction['coef'][0], float)
    if not np.isfinite(coef).all():
        raise ValueError('ML prediction contains nonfinite coefficients')
    ml = MLCoefficients(
        names=list(model.design_names),
        coef=coef,
        coef_std=np.asarray(prediction['coef_std'][0], float),
        oh_coef=coef.copy(),
    )
    diagnostics = {k: v for k, v in prediction.items() if k not in ('coef', 'coef_std')}
    # -1 means the prediction's own tests did not run; it stays -1 rather than reading as clean.
    bits = np.asarray(diagnostics['reliability'], np.int32)
    extra = rel.fallback_bits(median_fits.failed_arms, median_fits.labels['sky1'])
    diagnostics['reliability'] = np.where(bits < 0, bits, bits | extra).astype(np.int32)
    return ml, diagnostics


def corrected_coefficients(model, matrices, coefficient, parameters):
    """Apply median OH band scales and upper-level-energy tilts to the OH coefficients.

    Mirrors the band merging and tilt definition of ``sky_line_scaling.line_templates``,
    so that rendering the result adds the fitted template corrections exactly.
    """
    output = coefficient.copy()
    sl = model._component_slices(matrices)
    v, e = oh_group_labels(model)
    c = coefficient[sl['oh']]
    flux = np.maximum(c, 0) * matrices['oh'].sum(axis=1)
    bands = sorted(set(v.tolist()))
    keep = [b for b in bands if flux.sum() > 0 and flux[v == b].sum() >= 0.01 * flux.sum()] or bands[:1]
    assign = {b: min(keep, key=lambda k: (abs(k - b), k)) for b in bands}
    for k in keep:
        idx = np.flatnonzero(np.isin(v, [b for b in bands if assign[b] == k]))
        weight = flux[idx].sum()
        mean_energy = float(np.sum(flux[idx] * e[idx]) / weight) if weight > 0 else float(np.mean(e[idx]))
        scale = parameters.get(f'OH_v{k}', 1.0)
        tilt = parameters.get(f'OH_v{k}_tilt', 1.0) - 1.0
        output[sl['oh'].start + idx] *= scale + tilt * (e[idx] - mean_energy) / E_SCALE
    if np.any(output[sl['oh']] < 0) or not np.isfinite(output).all():
        raise ValueError('Median OH correction produces invalid line populations')
    return output


def _sky_arm(median, fit, wave):
    variance = np.divide(1, median.ivar, out=np.full_like(wave, np.inf), where=median.ivar > 0)
    continuum = (fit.components['moon'] + fit.components['zodi']) / FACTOR
    return SkyArm(median.flux, fit.bestfit_lsf / FACTOR, continuum, variance)


def estimate_corrections(frame, median_fits, ml, medians, stack, near, settings=SETTINGS):
    """Blue sky-arm residual transfer and Science-median OH band scales, estimated once."""
    wave, model = frame.wave, median_fits.science_model
    fits_by_label = median_fits.results
    raw, components, matrices = render(model, ml.coef, fits_by_label['Sci'].lsf_state)
    arms = [_sky_arm(medians[k], fits_by_label[k], wave) for k in SKY_ARMS if k not in median_fits.failed_arms]
    velocities = protection_velocities(*adapter.telescope_coords(frame.header, 'Sci'))
    lines = exposure_line_masks(wave, fits_by_label['Sci'].lsf_state, velocities, settings)
    correction, residual_info = sky_arm_residual_correction(
        wave,
        raw / FACTOR,
        (components['moon'] + components['zodi']) / FACTOR,
        arms,
        max_wavelength=settings.residual_max_angstrom,
        **residual_protection(settings, lines, components['oh'] / FACTOR),
        return_info=True,
    )
    # Reconstruction rebuilds the line masks from the Science median LSF and these velocities.
    residual_info['line_protection_velocities_km_s'] = list(velocities)
    science = medians['Sci']
    science_mask = data.science_line_mask_rows(
        stack,
        wave,
        science.flux,
        medians[near].flux,
        extra_lines=NEBULAR_LINES_RED,
        widen_if_unmeasured_km_s=UNMEASURED_VELOCITY_KM_S,
    )
    templates = line_templates(model, matrices, ml.coef, atoms=False, o2=False)
    variance = np.divide(1, science.ivar, out=np.full_like(wave, np.inf), where=science.ivar > 0)
    _, oh_info = sky_line_scaling_correction(
        wave,
        science.flux,
        raw / FACTOR + correction,
        {k: v / FACTOR for k, v in templates.items()},
        variance=variance,
        mask=science_mask,
        return_info=True,
    )
    oh_coef = ml.coef.copy()
    if oh_info['accepted']:
        try:
            oh_coef = corrected_coefficients(model, matrices, ml.coef, oh_info['scales'])
        except ValueError as error:
            oh_info.update(accepted=False, reason=str(error))
            LOG.warning('OH correction reverted to prediction: %s', error)
    return Corrections(matrices, arms, science_mask, oh_coef, residual_info, oh_info, lines)


def _survey_fit(reference, reason):
    qa = [dict(channel=ch, accepted=False, reason=reason) for ch in FIT_CHANNELS]
    return FiberFit(deepcopy(reference), deepcopy(reference), qa, [])


def _fiber_task(index):
    w = _WORK
    reference = w.references[index]
    fiberid = int(w.fiberids[index])
    flux, ivar = w.flux[index] * FACTOR, w.ivar[index] / FACTOR**2
    valid_ivar = np.where(w.mask[index] == 0, ivar, 0.0)
    if w.survey_only:
        fit = _survey_fit(reference, 'survey_only')
    else:
        try:
            fit = fit_fiber(w.wave, flux, valid_ivar, reference, w.comparison, w.designs, w.science_mask, w.settings)
        except Exception as error:
            # One failing fiber keeps its survey LSF instead of stopping the exposure.
            LOG.warning('Fiber %s: LSF fit failed, survey LSF kept: %r', fiberid, error)
            fit = _survey_fit(reference, f'fiber_error:{type(error).__name__}')
    sky = sky_parts(w.model, w.coef, w.oh_coef, fit.adopted, w.arms, w.settings, w.lines).total
    for record in fit.qa:
        record.update(fiberid=fiberid, input_row=index)
        channel = record['channel']
        if not record['accepted'] and reference.metrics[channel].get('quality') != 'measured':
            sky[channel_slices(w.wave, reference, (channel,))[channel]] = np.nan
            record['reason'] = 'unusable_survey_fallback'
    fit.history = [
        dict(input_row=int(index), fiberid=fiberid, **{**row, 'coefficients': row['coefficients'].ravel()})
        for row in fit.history
    ]
    return FiberResult(index, sky.astype('f4'), fit)


def fit_and_render_fibers(work, rows, jobs, started):
    """Fit and render every processed Science row; returns the sky image and per-row fits."""
    global _WORK
    sky = np.zeros(work.flux.shape, dtype='f4')
    fits_by_row = []
    _WORK = work
    pool = ProcessPoolExecutor(jobs, mp_context=mp.get_context('fork')) if jobs > 1 else None
    try:
        results = pool.map(_fiber_task, rows, chunksize=4) if pool else map(_fiber_task, rows)
        for n, result in enumerate(results, 1):
            sky[result.index] = result.sky
            fits_by_row.append(result.fit)
            if n % 100 == 0:
                LOG.info('Science fibers %s/%s; elapsed %.1f s', n, len(rows), time.perf_counter() - started)
    finally:
        if pool:
            pool.shutdown(cancel_futures=True)
        _WORK = None
    return sky, fits_by_row


def sky_arm_rows(slit, references):
    """Usable SkyE/SkyW rows of the CFrame."""
    return np.array([i for i in references if str(slit['telescope'][i]).strip() in SKY_ARMS], int)


def render_sky_arm_fibers(sky, slit, rows, median_fits):
    """Write each Sky fiber's telescope median model into ``sky``; returns the median LSF state used per row."""
    labels = [str(slit['telescope'][i]).strip() for i in rows]
    spectra = {
        label: arm_sky(median_fits.models[label], median_fits.results[label].coef, median_fits.results[label].lsf_state)
        for label in set(labels)
    }
    for row, label in zip(rows, labels, strict=True):
        sky[row] = spectra[label]
    return [median_fits.results[label].lsf_state for label in labels]


def run(
    input_path,
    ensemble_path,
    *,
    gaia_dir=None,
    lsf_calibration=None,
    output=None,
    jobs=1,
    survey_only=False,
    render_sky=True,
    return_hdus=False,
):
    """Process one lvmCFrame and publish ``lvmSFrame-skyml-<exposure>.fits``; returns its path.

    With ``return_hdus`` nothing is written and the product's HDUList is returned (``output`` only
    names the Gaia cache folder); the caller owns it.

    Workers are forked after the ML ensemble is released. Callers outside the CLI should
    limit native threads (OPENBLAS/OMP/MKL_NUM_THREADS=1) before importing numpy.
    """
    started = time.perf_counter()
    settings = SETTINGS
    code_sha256, revision = code_provenance()
    input_path, output = resolve_paths(input_path, output)
    if jobs < 1:
        raise ValueError('jobs must be positive')
    frame = read_cframe(input_path)
    expnum = int(frame.header['EXPOSURE'])
    calibration = Path(lsf_calibration).resolve() if lsf_calibration else BUNDLED_CALIBRATION
    references = load_fiber_references(frame, calibration)

    LOG.info('Exposure %s: Gaia photometry and masked medians', expnum)
    gaia = resolve_gaia(input_path, frame, gaia_dir, output.parent / '.skyml-gaia')
    selection = select_fibers(frame.slit, gaia.fibers, settings)
    medians = {k: median_spectrum(frame.flux, frame.ivar, frame.mask, v) for k, v in selection.indices.items()}
    median_references = {k: mean_reference([references[int(i)] for i in v]) for k, v in selection.indices.items()}
    stack, near, far = build_stack(frame, input_path, medians, median_references, selection.indices, settings)

    LOG.info('Selected fibers %s; decomposing medians', {k: len(v) for k, v in selection.indices.items()})
    labels = dict(sci='Sci', sky1=near, sky2=far)
    median_fits = adapter.fit_medians(stack, {kind: median_references[label] for kind, label in labels.items()}, labels)
    LOG.info('Predicting one Science IFU coefficient vector')
    ml, prediction_diagnostics = predict_ifu(ensemble_path, frame, median_fits)
    corrections = estimate_corrections(frame, median_fits, ml, medians, stack, near, settings)
    ml.oh_coef = corrections.oh_coef

    model = median_fits.science_model
    work = FiberWork(
        model=model,
        wave=frame.wave,
        flux=frame.flux,
        ivar=frame.ivar,
        mask=frame.mask,
        fiberids=np.asarray(frame.slit['fiberid']),
        references=references,
        comparison=median_fits.results['Sci'].lsf_state,
        designs=sky_designs(model, ml.oh_coef, median_references['Sci'], corrections.matrices),
        science_mask=corrections.science_mask,
        coef=ml.coef,
        oh_coef=ml.oh_coef,
        arms=corrections.arms,
        survey_only=survey_only,
        settings=settings,
        lines=corrections.lines,
    )
    rows = np.array([i for i in references if str(frame.slit['telescope'][i]).strip() == 'Sci'], int)
    LOG.info('Fitting and rendering %s Science fibers using %s processes', len(rows), jobs)
    if render_sky:
        sky, fiber_fits = fit_and_render_fibers(work, rows, jobs, started)
        arm_rows = sky_arm_rows(frame.slit, references)
        arm_states = render_sky_arm_fibers(sky, frame.slit, arm_rows, median_fits)
    else:
        # Preparation for a caller that renders its own sky: SKY stays NaN, so the product's rows are flagged invalid.
        sky = np.full(frame.flux.shape, np.nan, dtype='f4')
        fiber_fits = [_survey_fit(references[i], 'survey_only') for i in rows]
        arm_rows, arm_states = np.empty(0, int), []

    qa = [record for fit in fiber_fits for record in fit.qa]
    inputs = ProductInputs(
        sky=sky,
        rows=rows,
        adopted=[f.adopted for f in fiber_fits],
        candidates=[f.candidate for f in fiber_fits],
        qa=qa,
        history=[row for fit in fiber_fits for row in fit.history],
        photometry=selection.annotate(gaia.fibers),
        sources=gaia.sources,
        medians=medians,
        median_fits=median_fits,
        initial_states=[median_references[k] for k in medians],
        stack_meta=Table(stack['META'].data),
        ml=ml,
        fit_mask=corrections.science_mask,
        arm_rows=arm_rows,
        arm_states=arm_states,
    )
    provenance = dict(
        input=str(input_path),
        input_sha256=sha256(input_path),
        ensemble=str(Path(ensemble_path).resolve()),
        ensemble_sha256=sha256(ensemble_path),
        calibration=str(calibration),
        calibration_manifest=json.loads((calibration / 'manifest.json').read_text()),
        gaia=gaia.provenance,
        prediction_diagnostics=prediction_diagnostics,
        code_sha256=code_sha256,
        code_revision=revision,
        residual_correction=corrections.residual_info,
        oh_correction=corrections.oh_info,
        median_fits={
            label: dict(
                status=fit.fit_status,
                summary=fit.fit_summary,
                reliability=median_fits.reliability[label],
                r2=float(fit.r2),
                failed=label in median_fits.failed_arms,
            )
            for label, fit in median_fits.results.items()
        },
        selected_fiberids={k: work.fiberids[v] for k, v in selection.indices.items()},
        median_lsf_rule='equal mean of selected valid fiber profiles; median approximation',
        sky_fiber_model='SkyE/SkyW rows: telescope median decomposition rendered with the median LSF; '
        'FLUX is a residual of the fibers that build that median, not an independent test',
        median_variance='heteroscedastic Gaussian median asymptote; no spatial covariance',
        uncertainty='IVAR: observed only; full sky variance unknown; '
        'reconstructed SKY_EPI: diagonal coefficient ensemble approximation',
        settings=dict(
            fiber_fit_model='full_sky_decomp_plus_nonnegative_spline',
            ml_scope='IFU center',
            survey_only=survey_only,
            render_sky=render_sky,
            jobs=jobs,
            **settings.asdict(),
        ),
        elapsed_seconds=time.perf_counter() - started,
    )
    LOG.info(
        '%s; insitu accepted %s/%s channel fits',
        'Built in memory' if return_hdus else f'Writing {output}',
        sum(bool(r.get('accepted')) for r in qa),
        len(qa),
    )
    if return_hdus:
        return build_product(input_path, inputs, provenance)
    write_product(input_path, output, inputs, provenance)
    return output
