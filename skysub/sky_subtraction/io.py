"""CFrame input, Gaia photometry, fiber selection and masked median spectra."""

import logging
import time
import warnings
from pathlib import Path
from typing import NamedTuple

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from skysub.medians_computation.src.lvm_medians.stack import FLUX_UNIT, gaia_ratio_clean_mask

from . import decomp_adapter as adapter
from .settings import SETTINGS

LOG = logging.getLogger(__name__)
TELESCOPES = ('Sci', 'SkyE', 'SkyW')
GAIA_RETRIES = 3
BUNDLED_CALIBRATION = Path(__file__).resolve().parents[1] / 'sky_decomp/data/lsf'


class CFrame(NamedTuple):
    header: fits.Header
    slit: Table
    wave: np.ndarray
    flux: np.ndarray  # FLUX_UNIT
    ivar: np.ndarray
    mask: np.ndarray
    lsf: np.ndarray


class Median(NamedTuple):
    flux: np.ndarray
    ivar: np.ndarray
    count: np.ndarray


class GaiaPhotometry(NamedTuple):
    fibers: Table  # one row per SLITMAP row, in SLITMAP order
    sources: Table | None  # TAP source matches; None for an aggregate cache
    provenance: dict


class Selection(NamedTuple):
    indices: dict  # telescope -> selected row indices
    eligible: np.ndarray
    reasons: np.ndarray  # per row: gaia_clean, lowest_ratio_50, reduced_sample or not_selected

    def annotate(self, photometry):
        """Copy of the photometry table with the selection columns of the GAIA extension."""
        table = photometry.copy()
        table['eligible'] = self.eligible
        table['selected'] = self.reasons != 'not_selected'
        table['selection_reason'] = self.reasons
        return table


def flux_scale(header):
    return u.Unit(header['BUNIT']).to(FLUX_UNIT)


def read_cframe(path):
    if not Path(path).name.startswith('lvmCFrame-'):
        raise ValueError('Input must be an lvmCFrame')
    with fits.open(path, memmap=False) as h:
        scale = flux_scale(h['FLUX'].header)
        frame = CFrame(
            header=h[0].header.copy(),
            slit=Table(h['SLITMAP'].data),
            wave=np.asarray(h['WAVE'].data, dtype=float),
            flux=np.asarray(h['FLUX'].data, dtype=float) * scale,
            ivar=np.asarray(h['IVAR'].data, dtype=float) / scale**2,
            mask=np.asarray(h['MASK'].data),
            lsf=np.asarray(h['LSF'].data, dtype=float),
        )
    wave, slit = frame.wave, frame.slit
    if wave.ndim != 1 or not np.isfinite(wave).all() or np.any(np.diff(wave) <= 0):
        raise ValueError('WAVE must be a finite increasing one-dimensional grid')
    if any(a.shape != (len(slit), len(wave)) for a in (frame.flux, frame.ivar, frame.mask, frame.lsf)):
        raise ValueError('Spectral arrays, WAVE, and SLITMAP dimensions disagree')
    if len(set(slit['fiberid'])) != len(slit):
        raise ValueError('SLITMAP physical fiber IDs must be unique')
    return frame


# Earlier name; the result still unpacks as (header, slit, wave, flux, ivar, mask, lsf).
read_frame = read_cframe


def _cached_fibers(path, slit, expnum):
    """Aggregate cache rows in SLITMAP order, or ``(None, reason)`` when the cache does not match."""
    with fits.open(path) as h:
        cached = Table(h['FIBERS'].data)
        if h[0].header.get('EXPNUM') != expnum:
            return None, 'EXPNUM differs'
    if len(set(cached['fiberid'])) != len(cached) or set(cached['fiberid']) != set(slit['fiberid']):
        return None, 'fiber IDs differ from SLITMAP'
    by_id = {int(v): i for i, v in enumerate(cached['fiberid'])}
    cached = cached[[by_id[int(v)] for v in slit['fiberid']]]
    same_position = np.allclose(cached['fiber_ra'], slit['ra'], atol=1e-7, rtol=0, equal_nan=True) and np.allclose(
        cached['fiber_dec'], slit['dec'], atol=1e-7, rtol=0, equal_nan=True
    )
    if not same_position:
        return None, 'fiber coordinates differ from SLITMAP'
    return cached, None


def _query_sources(path, destination):
    """Three attempts cover transient TAP failures; the last error propagates."""
    for attempt in range(GAIA_RETRIES):
        try:
            return adapter.query_gaia_sources(path, destination)
        except Exception as error:
            if attempt == GAIA_RETRIES - 1:
                raise
            LOG.warning('Gaia TAP attempt %s failed: %s', attempt + 1, error)
            time.sleep(2**attempt)


def resolve_gaia(path, frame, directory, cache):
    """Per-fiber Gaia photometry from a matching aggregate cache, otherwise from ARI TAP.

    The observed (source + sky) G flux of every fiber is always recomputed from the CFrame.
    """
    expnum = int(frame.header['EXPOSURE'])
    fibers, sources, provenance, miss = None, None, None, 'no cache directory'
    candidate = None if directory is None else Path(directory) / f'lvmGAIA-fibers-{expnum:08d}.fits'
    if candidate is not None:
        miss = 'no cache file'
        if candidate.is_file():
            fibers, miss = _cached_fibers(candidate, frame.slit, expnum)
    if fibers is None:
        LOG.info('Gaia cache not used (%s); querying ARI TAP', miss)
        cache = Path(cache)
        cache.mkdir(parents=True, exist_ok=True)
        sources, provenance = _query_sources(path, cache / f'lvmGAIA-sources-{expnum:08d}.fits')
        fibers = adapter.gaia_fiber_table(path, sources)
        provenance['cache_miss'] = miss
    else:
        provenance = dict(mode='aggregate cache', file=str(candidate.resolve()), source_matches='aggregate only')
    observed, error, fraction = adapter.synthetic_g(frame.wave, frame.flux, frame.ivar, frame.mask)
    fibers['lvm_flux_plus_sky_g'] = observed
    fibers['lvm_flux_plus_sky_g_err'] = error
    fibers['lvm_flux_plus_sky_g_valid_fraction'] = fraction
    fibers['ratio_gaia_to_lvm_flux_plus_sky_g'] = np.divide(
        fibers['gaia_g_flux'], observed, out=np.full(len(frame.slit), np.nan), where=observed > 0
    )
    return GaiaPhotometry(fibers, sources, provenance)


def select_fibers(slit, photometry, settings=SETTINGS):
    """Gaia-clean usable fibers per telescope; Science falls back to its lowest-ratio fibers."""
    clean = gaia_ratio_clean_mask(photometry, settings.gaia_ratio_threshold)
    ratio = np.asarray(photometry['ratio_gaia_to_lvm_flux_plus_sky_g'], float)
    n, valid = np.asarray(photometry['n_gaia_sources']), np.asarray(photometry['gaia_g_n_valid_sources'])
    complete = (n == 0) | ((n == valid) & (valid > 0) & np.isfinite(ratio) & (ratio >= 0))
    eligible = (
        (np.asarray(slit['fibstatus']) == 0)
        & complete
        & np.isfinite(ratio)
        & (np.asarray(photometry['lvm_flux_plus_sky_g']) > 0)
        & (np.asarray(photometry['lvm_flux_plus_sky_g_valid_fraction']) >= settings.min_valid_fraction)
    )
    telescopes = np.char.strip(np.asarray(slit['telescope']).astype(str))
    minimum = settings.science_min_fibers
    selected, reasons = {}, np.full(len(slit), 'not_selected', dtype='U24')
    for label in TELESCOPES:
        candidates = np.flatnonzero(eligible & (telescopes == label))
        indices = np.flatnonzero(eligible & clean & (telescopes == label))
        reason = 'gaia_clean'
        if label == 'Sci' and len(indices) < minimum:
            order = np.lexsort((np.asarray(slit['fiberid'])[candidates], ratio[candidates]))
            indices = candidates[order[:minimum]]
            reason = f'lowest_ratio_{minimum}' if len(indices) == minimum else 'reduced_sample'
        if not len(indices):
            raise ValueError(f'{label}: no eligible Gaia-selected fibers')
        selected[label] = indices
        reasons[indices] = reason
    return Selection(selected, eligible, reasons)


def median_spectrum(flux, ivar, mask, indices):
    """Masked median, its Gaussian-asymptote inverse variance and contributor counts.

    Approximation: independent Gaussian fibers sharing one signal. Intrinsic spatial
    variation would need resampling or a covariance model. The heteroscedastic median
    variance is n / (4 (sum_i f_i(mu))^2).
    """
    f, v = flux[indices], ivar[indices]
    good = (mask[indices] == 0) & np.isfinite(f) & np.isfinite(v) & (v > 0)
    count = good.sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        med = np.nanmedian(np.where(good, f, np.nan), axis=0)
    density = np.sum(np.where(good, np.sqrt(np.maximum(v, 0) / (2 * np.pi)), 0), axis=0)
    inverse_variance = np.divide(4 * density**2, count, out=np.zeros_like(density), where=count > 0)
    single = count == 1
    inverse_variance[single] = np.sum(np.where(good, v, 0), axis=0)[single]
    inverse_variance[count == 0] = 0
    return Median(med, inverse_variance, count)
