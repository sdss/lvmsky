"""Sky rendering and robust fit weights shared by the exposure fit and reconstruction."""

from dataclasses import dataclass

import numpy as np

from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics

from .ml import SkyArm, sky_arm_residual_correction
from .settings import SETTINGS

# Decomposition coefficients and matrices work in flux units of 1e-14 erg/s/cm2/A.
FACTOR = 1e14
FAMILIES = ('oh', 'moon', 'zodi', 'diffuse', 'atom', 'orc', 'o2')

__all__ = [
    'FACTOR',
    'FAMILIES',
    'SkyArm',
    'SkyParts',
    'channel_slices',
    'fwhm_grid',
    'render',
    'sky_parts',
    'fit_variance',
    'huber_factor',
    'fit_pixels',
]


def channel_slices(wave, state, channels=None):
    """Boolean pixel selections from the LSF state's channel bounds; open bounds are unlimited."""
    out = {}
    for channel, (lo, hi) in state.channel_bounds.items():
        if channels is None or channel in channels:
            out[channel] = (wave >= (-np.inf if lo is None else lo)) & (wave < (np.inf if hi is None else hi))
    return out


def fwhm_grid(state, wave, settings=SETTINGS):
    """FWHM image row interpolated from profile diagnostics (not from LSF coefficients)."""
    out = np.empty_like(wave)
    for take in channel_slices(wave, state).values():
        probes = np.linspace(wave[take][0], wave[take][-1], settings.fwhm_probes)
        out[take] = np.interp(wave[take], probes, evaluate_lsf_diagnostics(state, probes)['fwhm_angstrom'])
    return out


def render(model, coefficient, state):
    """``(total, components, matrices)`` of the decomposition model with LSF ``state``."""
    model._set_lsf_state(state)
    matrices = model._assemble_refined_matrices()
    components = model._components_from_coef(coefficient, matrices)
    return sum((components[k] for k in FAMILIES), np.zeros_like(model.wave)), components, matrices


@dataclass
class SkyParts:
    """Physical-unit sky: ML prediction, OH correction and Blue sky-arm residual."""

    ml: np.ndarray
    oh: np.ndarray
    residual: np.ndarray
    matrices: dict

    @property
    def total(self):
        return self.ml + self.oh + self.residual


def sky_parts(model, coefficient, oh_coefficient, state, arms, settings=SETTINGS):
    """The subtracted sky of one fiber; identical for the pipeline and reconstruction."""
    raw, components, matrices = render(model, coefficient, state)
    corrected = model._components_from_coef(oh_coefficient, matrices)
    continuum = (components['moon'] + components['zodi']) / FACTOR
    residual = sky_arm_residual_correction(
        model.wave, raw / FACTOR, continuum, arms, max_wavelength=settings.residual_max_angstrom
    )
    return SkyParts(raw / FACTOR, (corrected['oh'] - components['oh']) / FACTOR, residual, matrices)


def fit_pixels(flux, ivar, science_mask):
    """Pixels entering the fiber LSF fit; ``ivar`` already carries the input mask."""
    return np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & ~science_mask


def fit_variance(ivar, good, continuum, settings=SETTINGS):
    """Observed variance plus a fractional floor on the source continuum."""
    return np.divide(1, ivar, out=np.full_like(ivar, np.inf), where=good) + (settings.variance_floor * continuum) ** 2


def huber_factor(z, scale, settings=SETTINGS):
    return np.minimum(1.0, settings.huber_threshold * scale / np.maximum(np.abs(z), 1e-30))
