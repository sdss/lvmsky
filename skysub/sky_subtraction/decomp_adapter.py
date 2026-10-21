"""The only contact with private helpers of decompose_parallel, sky_decomp and lvm_medians.

The functions here are thin: they name a private call and fix its arguments, while the
logic stays in the source modules. Model methods of the decomposer (``model._...``) are
used by rendering.py and fiber_lsf.py.
"""

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from skysub import decompose_parallel as dp
from skysub.sky_decomp import reliability as rel
from skysub.medians_computation.src.lvm_medians import gaia
from skysub.medians_computation.src.lvm_medians.stack import _coords, _separation
from skysub.sky_decomp.lsf_spline2d import _line_design
from skysub.sky_decomp.lsf_surface_iterative import _fit_bspline_channel, _wave_fingerprint

from .rendering import FACTOR

KINDS = ('sci', 'sky1', 'sky2')
SKY_LABELS = ('SkyE', 'SkyW')
GAIA_TAP_TIMEOUT_SECONDS = 120
GAIA_TAP_MAX_RECORDS = 1_000_000

line_design = _line_design
wave_fingerprint = _wave_fingerprint
zodi_correction = dp.SPLIT_ZODI_ZODI_CORRECTION


def telescope_coords(header, label):
    return _coords(header, label)


def separation(a, b):
    return _separation(a, b)


@dataclass
class MedianFits:
    """Median decompositions plus everything that depends on the decomposition worker state."""

    labels: dict  # kind -> telescope label: sci/sky1/sky2 -> Sci/near/far
    results: dict  # telescope label -> decomposition result
    reliability: dict  # telescope label -> reliability columns
    science_model: object  # decomposer for the Science median row
    fit_mask: np.ndarray  # science-line mask of the median fits (MED_FIT_MASK)
    render_context: dict  # constructor settings and per-median telluric parameters
    transmission: np.ndarray  # DRP telluric transmission for Sci, near, far
    models: dict = field(default_factory=dict)  # telescope label -> decomposer with its median O2 temperature
    failed_arms: tuple = ()  # sky telescopes whose median fit is a misfit; the other arm stands in


def fit_medians(stack, states, labels):
    """Decompose the three median rows with fixed Blue and R/Z seeded from ``states``.

    ``decompose_in_process`` initializes module globals of ``decompose_parallel``; the
    decomposers, masks and constructor settings derived from them are captured here once.
    Its ``initial_lsf_states`` argument replaces the Blue reference, so B stays fixed.
    """
    fitted = dp.decompose_in_process(stack, initial_lsf_states=states)
    if dp._WORKER_FACTOR != FACTOR:
        raise ValueError('Decomposition flux factor differs from the SkyML renderer')
    results = {labels[kind]: fitted[(kind, 0)][0] for kind in KINDS}
    reliability = {labels[kind]: fitted[(kind, 0)][1] for kind in KINDS}
    for kind in KINDS:
        label, fit = labels[kind], results[labels[kind]]
        if fit.fit_status not in ('Solved', 'AlmostSolved') or not np.isfinite(fit.coef).all():
            raise RuntimeError(f'{label} median decomposition: {fit.fit_status} {fit.fit_summary}')
        if not np.array_equal(fit.lsf_state.coefficients['B'], states[kind].coefficients['B']):
            raise RuntimeError('Blue LSF changed during median fit')
    # A solved fit can still describe nothing, and a good fit can describe Moon glare instead of sky
    # (exposures 13190, 13248).  Science fits are exempt: their R2 is low whenever the target is bright.
    failed = tuple(
        label for label in SKY_LABELS if int(reliability[label]['reliability']) & rel.ARM_FAILURE_MASK
    )
    if len(failed) == len(SKY_LABELS):
        raise RuntimeError(
            'Both sky-arm median decompositions are unusable: '
            + '; '.join(f'{label} {results[label].fit_summary}' for label in failed)
        )
    settings = dict(dp._WORKER_DECOMPOSER_KWARGS)
    settings['config'] = asdict(settings['config'])
    settings['base_dir'] = str(settings['base_dir'])
    models, transmission = [], []
    by_label = {}
    science_model = None
    for kind in KINDS:
        label = labels[kind]
        model = dp._telluric_decomposer(kind, 0)
        model.t_o2 = results[label].t_o2
        by_label[label] = model
        models.append(
            dict(
                telescope=label,
                pwv_mm=model.pwv_mm,
                source_airmass=model.source_airmass,
                t_o2=float(results[label].t_o2),
            )
        )
        transmission.append(model.drp_transmission)
        if kind == 'sci':
            science_model = model
    science_model.t_o2 = results['Sci'].t_o2
    return MedianFits(
        labels=dict(labels),
        results=results,
        reliability=reliability,
        science_model=science_model,
        models=by_label,
        fit_mask=dp._science_line_mask_for_row(0).astype('u1'),
        render_context=dict(settings=settings, models=models),
        transmission=np.asarray(transmission),
        failed_arms=failed,
    )


def fit_lsf_surface(wave, source, target, weight, initial, knots, degree, kernel, background, settings):
    """Constrained LSF solve with a nonnegative background; returns ``(coefficient, metric)``."""
    _, coefficient, _, metric = _fit_bspline_channel(
        wave,
        source,
        target,
        weight,
        np.median(initial, axis=1),
        n_basis=initial.shape[1],
        degree=degree,
        knot_vector=knots,
        kernel_design=kernel,
        background_design=background,
        background_nonnegative=True,
        return_background=True,
        center_offsets=True,
        roughness_fraction=0,
        offset_roughness_fraction=settings.offset_roughness_fraction,
        fallback_prior_fraction=settings.fallback_prior_fraction,
        numerical_ridge=0,
    )
    return coefficient, metric


def query_gaia_sources(cframe, destination):
    """Run the lvm_medians ARI TAP query and cache its source table at ``destination``."""
    service = gaia.TAP_SERVICES['ari']
    table = gaia._download(
        Path(cframe),
        Path(destination),
        service,
        timeout=GAIA_TAP_TIMEOUT_SECONDS,
        token=None,
        maxrec=GAIA_TAP_MAX_RECORDS,
    )
    provenance = dict(mode='ARI TAP', query=gaia.QUERY, service=service, aperture_arcsec=gaia.APERTURE_ARCSEC)
    return table, provenance


def gaia_fiber_table(cframe, sources):
    table, _ = gaia.derive_fiber_table(Path(cframe), sources, observed_total=True)
    return table


def synthetic_g(wave, flux, ivar, mask):
    """Observed (source + sky) Gaia G flux, error and valid fraction for every fiber."""
    passwave, transmission, _, _ = gaia._passband(None)
    weights = np.interp(wave, passwave, transmission, left=0, right=0) * wave * gaia._trapezoid_weights(wave)
    variance = np.divide(1, ivar, out=np.full_like(ivar, np.nan), where=np.isfinite(ivar) & (ivar > 0))
    return gaia._synthetic(flux, variance, weights, mask == 0)
