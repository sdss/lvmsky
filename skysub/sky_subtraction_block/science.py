"""Science-fiber sky for the block LSF, evaluated for all fibers at once.

The sky renderer is exactly affine in the LSF coefficients (it is a convolution of fixed matrices), and every
Science fiber shares the ML coefficients, so ``outputs = offset + h @ D`` holds with ``h`` the stacked B/R/Z
coefficients of a fiber. ``offset`` and ``D`` come from 1 + 3*n_offset*n_wave renders of the production renderer.
"""
from copy import copy, deepcopy
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import BSpline
from skysub.sky_decomp.lsf_spline2d import mspline_basis
from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
from skysub.sky_subtraction import pipeline as common
from skysub.sky_subtraction.rendering import channel_slices,residual_protection
from mlp_predictor.sky_arm_correction import sky_arm_residual_correction
from .lsf import FIT_CHANNELS, evaluate_block

CHANNELS = 'BRZ'


def _outputs(r, state):
    """Raw ML sky, OH correction, Moon + zodiacal light and raw OH of the Science model, in physical units."""
    raw, components, matrices = common.render(r.model_sci, r.coef, state)
    oh = r.model_sci._components_from_coef(r.oh_coef, matrices)['oh'] - components['oh']
    return np.stack([raw, oh, components['moon'] + components['zodi'], components['oh']]) / common.FACTOR


def affine_sky(r, template):
    """``(offset, D)``: outputs of a stacked coefficient vector ``h`` are ``offset + h @ D.reshape(len(h), -1)``."""
    state = deepcopy(template)
    for c in CHANNELS:
        state.coefficients[c] = np.zeros_like(state.coefficients[c])
    offset = _outputs(r, state)
    rows = []
    for c in CHANNELS:
        flat = state.coefficients[c].reshape(-1)  # a view: writing a unit entry changes the state
        for k in range(flat.size):
            flat[k] = 1.
            rows.append(_outputs(r, state) - offset)
            flat[k] = 0.
    return offset, np.stack(rows)


def fwhm_batch(masses, half_width):
    """FWHM of unit-grid profiles ``masses`` (n, n_offset); same definition as ``evaluate_lsf_diagnostics``."""
    delta = np.linspace(-half_width, half_width, 1201)
    density = masses @ mspline_basis(delta, n_basis=masses.shape[1], half_width=half_width).T
    n, m = density.shape
    index = np.arange(m)
    peak = density.argmax(axis=1)
    half = 0.5 * density[np.arange(n), peak]
    below = density <= half[:, None]
    lo = np.where(below & (index <= peak[:, None]), index, -1).max(axis=1)
    hi = np.where(below & (index >= peak[:, None]), index, m).min(axis=1)
    ok = (lo >= 0) & (hi < m)
    lo, hi = np.clip(lo, 0, m - 2), np.clip(hi, 1, m - 1)
    row = np.arange(n)
    p0, p1, p2, p3 = (density[row, i] for i in (lo, lo + 1, hi - 1, hi))
    x0 = delta[lo] + (half - p0) / (p1 - p0) * (delta[lo + 1] - delta[lo])
    x1 = delta[hi] + (half - p3) / (p2 - p3) * (delta[hi - 1] - delta[hi])
    return np.where(ok, x1 - x0, np.nan)


def _fwhm_rows(state, coefficients, wave, channel, probes_count):
    """Per-fiber FWHM over the channel's pixels; ``coefficients`` is (n_fiber, n_offset, n_wave_basis)."""
    take = channel_slices(wave, state, (channel,))[channel]
    probes = np.linspace(wave[take][0], wave[take][-1], probes_count)
    knots, degree = state.knot_vectors[channel], state.degrees[channel]
    basis = evaluate_bspline_basis(np.clip(probes, knots[degree], knots[-degree - 1]), knots, degree)
    basis = basis.toarray() if sp.issparse(basis) else np.asarray(basis)
    masses = np.einsum('pq,fkq->fpk', basis, coefficients)
    values = fwhm_batch(masses.reshape(-1, masses.shape[-1]), state.config.get('offset_half_width_angstrom', 3.0))
    values = values.reshape(len(coefficients), probes_count)
    # np.interp reproduces fwhm_grid, including its handling of undefined probes.
    return take, np.stack([np.interp(wave[take], probes, v) for v in values])


def uniform_metadata(r):
    """True when every Science fiber's adopted B-channel LSF metadata equals the first fiber's (only coefficients differ)."""
    meta = r.h['LSF_META'].data
    n = int(meta['spectrum_index'].max()) + 1
    for name in meta.columns.names:
        if name == 'spectrum_index':
            continue
        values = np.asarray(meta[name])
        values = values.reshape(n, len(CHANNELS), *values.shape[1:])
        if not np.all((values == values[:1]) | ((values != values) & (values[:1] != values[:1]))):
            return False
    return True


def evaluate_science(r, blocks, slit, science, lsf_row, chunk=128):
    """``(sky, states, lsf_rows)`` of the Science rows; replaces one full render per fiber by a matrix product."""
    first = int(science[0])
    template = evaluate_block(r.state('MED_', r.labels.index('Sci')), blocks, int(slit['spectrographid'][first]), float(slit['fiberid'][first]))
    offset, D = affine_sky(r, template)
    D = D.reshape(len(D), -1)
    shapes = {c: template.coefficients[c].shape for c in CHANNELS}
    cube = r.h['LSF_COEF'].data
    position = np.asarray(slit['fiberid'])[science].astype(float)
    spectrograph = np.asarray(slit['spectrographid'])[science]
    cube_rows = np.array([lsf_row[int(i)] for i in science])
    coefficients = {c: np.empty((len(science),) + shapes[c]) for c in CHANNELS}
    coefficients['B'][:] = cube[cube_rows, 0, :shapes['B'][0], :shapes['B'][1]]
    for sp_id in np.unique(spectrograph):
        take = spectrograph == sp_id
        for c in FIT_CHANNELS:
            block = blocks[(int(sp_id), c)]
            basis = BSpline.design_matrix(position[take], block['slit_knots'], 3).toarray()
            coefficients[c][take] = np.einsum('kqs,fs->fkq', block['coef'], basis)
    h = np.concatenate([coefficients[c].reshape(len(science), -1) for c in CHANNELS], axis=1)

    sky = np.zeros((len(slit), len(r.wave)), dtype='f4')
    for start in range(0, len(science), chunk):
        part = slice(start, start + chunk)
        raw, oh, solar, raw_oh = (h[part] @ D).reshape(-1, 4, len(r.wave)).transpose(1, 0, 2) + offset[:, None, :]
        residual = sky_arm_residual_correction(r.wave, raw, solar, r.arms, max_wavelength=r.settings.residual_max_angstrom,
            **residual_protection(r.settings, r.lines, raw_oh))
        sky[science[part]] = (raw + oh + residual).astype('f4')

    # One LSF state per fiber for the product; B metadata is common to all fibers (checked), coefficients differ.
    if not uniform_metadata(r):
        raise ValueError('Preparation LSF metadata differs between fibers; fast state assembly does not apply')
    base = {int(s): evaluate_block(r.state('', int(cube_rows[spectrograph == s][0])), blocks, s, position[spectrograph == s][0])
            for s in np.unique(spectrograph)}
    states = []
    for i in range(len(science)):
        state = copy(base[int(spectrograph[i])])
        state.coefficients = {c: coefficients[c][i] for c in CHANNELS}
        states.append(state)

    # B widths stay those of the preparation; R and Z come from the block profiles.
    widths = np.array(r.h['LSF'].data[science], dtype=float)
    for c in FIT_CHANNELS:
        for s in np.unique(spectrograph):
            take = spectrograph == s
            pixels, values = _fwhm_rows(base[int(s)], coefficients[c][take], r.wave, c, r.settings.fwhm_probes)
            widths[np.flatnonzero(take)[:, None], np.flatnonzero(pixels)[None, :]] = values
    return sky, states, widths
