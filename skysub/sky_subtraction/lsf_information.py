"""Conservative observed-line information gate across a complete R/Z channel."""

import numpy as np
from scipy.signal import find_peaks

from .settings import SETTINGS


def line_information(wave, observed, ivar, lines, continuum, good, settings=SETTINGS):
    """Require independently detected, high-contrast sky features in every wavelength bin.

    The cuts are conservative operating points, not a calibrated recovery probability;
    calibrating them needs held-out real or injected spectra.
    """
    n_bins = settings.line_bins
    step = np.median(np.diff(wave))
    edges = np.linspace(wave[0], wave[-1] + step, n_bins + 1)
    counts = np.zeros(n_bins, dtype=int)
    contrasts = [[] for _ in range(n_bins)]
    snrs = [[] for _ in range(n_bins)]
    distance = max(1, int(np.ceil(settings.line_min_separation_angstrom / step)))
    peaks = find_peaks(np.where(np.isfinite(lines), lines, 0), distance=distance)[0]
    valid = good & np.isfinite(observed) & np.isfinite(lines) & np.isfinite(continuum) & np.isfinite(ivar) & (ivar > 0)
    inner, outer = settings.line_side_angstrom
    for i in peaks:
        if not valid[i] or lines[i] <= 0:
            continue
        offset = np.abs(wave - wave[i])
        core = offset <= settings.line_core_angstrom
        used = core & valid
        side = (offset >= inner) & (offset <= outer) & valid
        if (
            used.sum() < settings.line_min_pixels
            or used.sum() < settings.line_min_core_fraction * core.sum()
            or side.sum() < settings.line_min_pixels
        ):
            continue
        # Project out a local continuum and slope before measuring the observed line.
        design = np.column_stack([lines[used] / lines[i], np.ones(used.sum()), wave[used] - wave[i]])
        weighted = design * np.sqrt(ivar[used, None])
        if np.linalg.matrix_rank(weighted) < 3:
            continue
        covariance = np.linalg.inv(weighted.T @ weighted)
        coefficient = np.linalg.lstsq(weighted, observed[used] * np.sqrt(ivar[used]), rcond=None)[0]
        amplitude = float(coefficient[0])
        error = float(np.sqrt(covariance[0, 0]))
        level = max(float(np.percentile(observed[side], 75)), float(np.percentile(continuum[side], 75)), 0.0)
        contrast = amplitude / max(level, error)
        snr = amplitude / error
        j = min(n_bins - 1, int(np.searchsorted(edges, wave[i], side='right') - 1))
        contrasts[j].append(contrast)
        snrs[j].append(snr)
        counts[j] += int(contrast >= settings.line_min_contrast and snr >= settings.line_min_snr)
    return dict(
        line_information_pass=bool(np.all(counts >= settings.line_min_per_bin)),
        quality_bin_counts=counts,
        quality_bin_contrast=np.array([np.median(v) if v else np.nan for v in contrasts]),
        quality_bin_snr=np.array([np.median(v) if v else np.nan for v in snrs]),
    )
