"""Per-fiber R/Z LSF recovery: full sky-model kernel plus a nonnegative source-continuum spline."""

from copy import deepcopy
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import lsq_linear

from skysub.sky_decomp.lsf_reference import validate_coefficients
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics, native_pixel_edges
from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis

from .decomp_adapter import fit_lsf_surface, line_design
from .lsf_information import line_information
from .rendering import channel_slices, fit_pixels, fit_variance, huber_factor
from .settings import SETTINGS

FIT_CHANNELS = ('R', 'Z')
WIDTH_METRICS = ('w50_angstrom', 'w90_angstrom')
SOLVED = ('Solved', 'AlmostSolved')


def channel_design(model, coefficient, state, channel):
    """Linear map from the channel's LSF coefficients to the rendered sky-line spectrum."""
    indices = model._line_indices[channel]
    line_flux = model._line_weights()[indices] * coefficient[model._line_group[indices]]
    knots, degree = state.knot_vectors[channel], state.degrees[channel]
    coordinates = np.clip(model._line_wave[indices], knots[degree], knots[-degree - 1])
    design = line_design(
        model._line_components[channel],
        line_flux,
        evaluate_bspline_basis(coordinates, knots, degree),
        model.offset_basis_count,
    )
    take = channel_slices(model.wave, state, (channel,))[channel]
    return take, design[take].toarray()


def sky_designs(model, coefficient, state, matrices):
    """Exact linear LSF tensors ``{channel: (take, dict(lines, solar, fixed))}`` of the sky renderer."""
    slices = model._component_slices(matrices)
    line_coef = np.concatenate([coefficient[slices[k]] for k in model._LINE_KEYS])
    solar_source = (
        model.matrix_moon_hr.T @ coefficient[slices['moon']] + model.matrix_zodi_hr.T @ coefficient[slices['zodi']]
    )
    widths = np.diff(native_pixel_edges(model.wave))
    expected_solar = model._build_continuum_operator(state) @ solar_source
    diffuse = matrices['diffuse'].T @ coefficient[slices['diffuse']]
    designs = {}
    for channel in FIT_CHANNELS:
        take, lines = channel_design(model, line_coef, state, channel)
        knots, degree = state.knot_vectors[channel], state.degrees[channel]
        basis = evaluate_bspline_basis(np.clip(model.wave[take], knots[degree], knots[-degree - 1]), knots, degree)
        solar = line_design(
            model._continuum_components[channel], solar_source[take] * widths[take], basis, model.offset_basis_count
        )[take].toarray()
        # Continuum rendering has small cross-channel edge contributions from neighboring source bins.
        edge = expected_solar[take] - solar @ state.coefficients[channel].ravel()
        designs[channel] = (take, dict(lines=lines, solar=solar, fixed=diffuse[take] + edge))
    return designs


def continuum_knots(wave, spacing):
    return np.r_[np.repeat(wave[0], 4), np.arange(wave[0] + spacing, wave[-1], spacing), np.repeat(wave[-1], 4)]


@dataclass
class ChannelFit:
    adopted: np.ndarray  # coefficients applied to the sky: in-situ if accepted, otherwise survey
    candidate: object  # reference state carrying the last in-situ solution of this channel
    qa: dict
    history: list


@dataclass
class FiberFit:
    adopted: object
    candidate: object
    qa: list
    history: list


def _width_change(current, reference):
    return max(float(np.max(np.abs(current[k] / reference[k] - 1))) for k in WIDTH_METRICS)


def _centroid_change(current, reference):
    return float(np.max(np.abs(current['centroid_angstrom'] - reference['centroid_angstrom'])))


def _initial_qa(channel, active, covered, settings):
    n_bins = settings.line_bins
    return dict(
        channel=channel,
        accepted=False,
        reason='insufficient_lines',
        n_line_pixels=int(active.sum()),
        n_supported_basis=int(covered.sum()),
        iterations=0,
        chi2_red=np.nan,
        width_change=np.nan,
        centroid_change=np.nan,
        line_information_pass=False,
        quality_bin_counts=np.zeros(n_bins, dtype=int),
        quality_bin_contrast=np.full(n_bins, np.nan),
        quality_bin_snr=np.full(n_bins, np.nan),
    )


def unfitted_qa(channel, reason, settings=SETTINGS):
    """LSF_QA record, with the schema of a fitted channel, for a profile that no per-fiber fit produced."""
    qa = _initial_qa(channel, np.zeros(0, bool), np.zeros(0, bool), settings)
    qa.update(accepted=True, reason=reason)
    qa.update(
        line_amplitude=np.nan, solar_amplitude=np.nan, downweighted=np.nan, converged=True, median_width_change=np.nan
    )
    return qa


def _solve_amplitudes_and_continuum(line_model, solar_model, continuum, target, weight, good, settings):
    """Nonnegative line/solar amplitudes with a unit prior and source continuum; None if BVLS fails.

    A positive source continuum cannot cancel a wrong sky-line profile.
    """
    templates = np.column_stack([line_model, solar_model, continuum])
    root = np.sqrt(weight[good])
    a = templates[good] * root[:, None]
    y = target[good] * root
    prior = np.zeros((2, a.shape[1]))
    prior[:, :2] = np.eye(2) / settings.amplitude_prior_sigma
    fit = lsq_linear(
        np.vstack([a, prior]),
        np.r_[y, np.ones(2) / settings.amplitude_prior_sigma],
        bounds=(0, np.inf),
        method='bvls',
        tol=1e-7,
    )
    if not fit.success:
        return None
    return fit.x[:2], fit.x[2:].copy()


def _rejection(coefficient, current, changes, information, converged, settings):
    """Name of the first failed acceptance check, or None when the in-situ profile is adopted."""
    try:
        validate_coefficients(coefficient)
        valid = True
    except ValueError:
        valid = False
    low, high = settings.lsf_fwhm_range_angstrom
    fwhm = current['fwhm_angstrom']
    guards = (
        np.all((fwhm >= low) & (fwhm <= high))
        and changes['width_change'] <= settings.lsf_max_relative_width_change
        and changes['centroid_change'] <= settings.lsf_max_centroid_change_angstrom
        and changes['median_width_change'] <= settings.lsf_max_median_width_change
    )
    checks = (
        ('invalid_profile', valid),
        ('line_information', information['line_information_pass']),
        ('not_converged', converged),
        ('profile_guard', guards),
    )
    return next((reason for reason, passed in checks if not passed), None)


def fit_channel(wave, flux, ivar, design, reference, channel, science_mask, comparison, settings=SETTINGS):
    """Alternate nonnegative amplitudes/continuum with a constrained full-sky LSF solve.

    ``design`` holds the ``lines``, ``solar`` and ``fixed`` arrays of ``sky_designs``. The
    survey ``reference`` is adopted unless the fit converges and passes every guard;
    ``comparison`` is the Science-median LSF.
    """
    initial = reference.coefficients[channel]
    knots, degree = reference.knot_vectors[channel], reference.degrees[channel]
    lines, solar, fixed = design['lines'], design['solar'], design['fixed']
    spline_knots = continuum_knots(wave, settings.continuum_spacing_angstrom)
    continuum = BSpline.design_matrix(wave, spline_knots, 3).toarray()
    good = fit_pixels(flux, ivar, science_mask)
    line = lines @ initial.ravel()
    active = good & (np.abs(line) > max(np.percentile(np.abs(line), 99) * settings.active_line_fraction, 1e-30))
    covered = evaluate_bspline_basis(wave, knots, degree)[active].sum(axis=0) > settings.min_pixels_per_basis
    qa = _initial_qa(channel, active, covered, settings)
    candidate = deepcopy(reference)
    history = []
    if active.sum() < settings.min_line_pixels_per_coefficient * initial.size or not covered.all():
        return ChannelFit(initial.copy(), candidate, qa, history)

    target = flux - fixed
    weight = np.where(good, ivar, 0.0)
    coefficient = initial.copy()
    probes = np.linspace(wave[0], wave[-1], settings.n_probes)
    prior_profile = previous = evaluate_lsf_diagnostics(reference, probes)
    dof = max(good.sum() - initial.size - continuum.shape[1] - 2, 1)
    converged, failure = False, None
    for iteration in range(1, settings.max_iterations + 1):
        solution = _solve_amplitudes_and_continuum(
            lines @ coefficient.ravel(), solar @ coefficient.ravel(), continuum, target, weight, good, settings
        )
        if solution is None:
            failure = 'continuum_solver'
            break
        amplitudes, variance_continuum_coef = solution
        # The robust variance uses the source continuum from before this LSF update.
        variance = fit_variance(ivar, good, continuum @ variance_continuum_coef, settings)
        kernel = amplitudes[0] * lines + amplitudes[1] * solar
        solved, metric = fit_lsf_surface(
            wave, line, target, weight, initial, knots, degree, kernel, continuum, settings
        )
        if metric['status'] not in SOLVED:
            failure = 'lsf_solver_' + metric['reason']
            break
        coefficient = solved
        continuum_coef = np.maximum(metric.pop('background_coef'), 0.0)
        stellar = continuum @ continuum_coef
        sky = kernel @ coefficient.ravel() + fixed
        residual = flux - sky - stellar
        z = residual * np.sqrt(np.divide(1, variance, out=np.zeros_like(weight), where=good))
        scale = max(1.0, 1.4826 * np.median(np.abs(z[active] - np.median(z[active]))))
        huber = huber_factor(z, scale, settings)
        weight = np.where(good, huber / variance, 0.0)
        candidate.coefficients[channel] = coefficient
        current = evaluate_lsf_diagnostics(candidate, probes)
        width_step = _width_change(current, previous)
        centroid_step = _centroid_change(current, previous)
        chi2 = float(np.sum((residual[good] ** 2) * weight[good]) / dof)
        history.append(
            dict(
                iteration=iteration,
                chi2_red=chi2,
                line_amplitude=float(amplitudes[0]),
                solar_amplitude=float(amplitudes[1]),
                relative_width_step=width_step,
                centroid_step=centroid_step,
                coefficients=coefficient.copy(),
                probe_wave=probes.copy(),
                w50=current['w50_angstrom'].copy(),
                w90=current['w90_angstrom'].copy(),
                continuum_knots=spline_knots.copy(),
                continuum_coef=continuum_coef.copy(),
                variance_continuum_coef=variance_continuum_coef,
                huber_scale=float(scale),
            )
        )
        line_model = amplitudes[0] * lines @ coefficient.ravel()
        background = amplitudes[1] * solar @ coefficient.ravel() + fixed + stellar
        qa.update(
            iterations=iteration,
            chi2_red=chi2,
            line_amplitude=float(amplitudes[0]),
            solar_amplitude=float(amplitudes[1]),
            downweighted=float(np.mean(huber[active] < 1)),
        )
        if (
            iteration > 1
            and width_step < settings.convergence_width
            and centroid_step < settings.convergence_centroid_angstrom
        ):
            converged = True
            break
        previous = current

    if failure is not None:
        # The survey profile is kept; the candidate holds the last successful iteration, if any.
        qa['reason'] = failure
        return ChannelFit(initial.copy(), candidate, qa, history)
    changes = dict(
        width_change=_width_change(current, prior_profile),
        centroid_change=_centroid_change(current, prior_profile),
        median_width_change=_width_change(current, evaluate_lsf_diagnostics(comparison, probes)),
    )
    information = line_information(wave, flux, ivar, line_model, background, good, settings)
    qa.update(information)
    reason = _rejection(coefficient, current, changes, information, converged, settings)
    accepted = reason is None
    qa.update(accepted=accepted, converged=converged, reason=reason or 'insitu', **changes)
    return ChannelFit(coefficient if accepted else initial.copy(), candidate, qa, history)


def fit_fiber(wave, flux, ivar, reference, comparison, designs, science_mask, settings=SETTINGS):
    """Fit R and Z of one fiber; ``ivar`` is zero on masked input pixels."""
    adopted, candidate, records, history = deepcopy(reference), deepcopy(reference), [], []
    for channel in FIT_CHANNELS:
        take, design = designs[channel]
        fit = fit_channel(
            wave[take], flux[take], ivar[take], design, reference, channel, science_mask[take], comparison, settings
        )
        history.extend(dict(channel=channel, **row) for row in fit.history)
        adopted.coefficients[channel] = fit.adopted
        candidate.coefficients[channel] = fit.candidate.coefficients[channel]
        adopted.metrics[channel] = fit.qa.copy()
        candidate.metrics[channel] = fit.qa.copy()
        records.append(fit.qa)
    adopted.fit_status = 'fiber_insitu_or_survey'
    candidate.fit_status = 'fiber_candidate'
    adopted.completed_cycles = candidate.completed_cycles = max(r['iterations'] for r in records)
    adopted.requested_cycles = candidate.requested_cycles = settings.max_iterations
    return FiberFit(adopted, candidate, records, history)
