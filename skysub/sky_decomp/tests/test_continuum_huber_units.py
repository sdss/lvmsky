"""What `continuum_fit_weights` calls an outlier.

The Huber factor is multiplied INTO ivar, so the outlier test and the
weighting should be in the same units. The historical form tests the RAW
residual against a channel-wide MAD while weighting by inverse variance;
those agree only when ivar is constant across the channel, which it is not
once the fit is photon-weighted.
"""

import numpy as np
import pytest

from skysub.sky_decomp.lsf_surface_iterative import (
    LSF_CHANNELS, _channel_mask, _robust_mad, continuum_fit_weights)


N = 4000


def _grid():
    wave = np.linspace(3700.0, 9600.0, N)
    skyline = np.zeros(N, dtype=bool)
    return wave, skyline


def test_flat_ivar_reduces_to_the_raw_residual_form():
    """With constant ivar the noise-unit score is a pure rescaling of the raw
    residual, and the per-channel renormalisation removes even that -- so an
    UNWEIGHTED fit is bit-unaffected by the change."""
    rng = np.random.default_rng(0)
    wave, skyline = _grid()
    resid = rng.normal(size=N)
    ivar = np.full(N, 7.0)                       # constant, but not 1
    w_new, _ = continuum_fit_weights(wave, resid, ivar, skyline)
    w_old = _legacy_weights(wave, resid, ivar, skyline)
    np.testing.assert_allclose(w_new, w_old, rtol=1e-10, atol=1e-12)


def _legacy_weights(wave, residual, ivar, skyline, k=3.0):
    """The pre-2026-09-20 form, kept HERE rather than in production.

    It scored outliers on the RAW residual against a channel-wide MAD while
    the factor it produced was multiplied into ivar. Reproduced locally so the
    tests can still show the direction of the change without carrying a flag
    through the fitter.
    """
    valid = np.isfinite(residual) & np.isfinite(ivar) & (ivar > 0.0) & ~skyline
    g = _robust_mad(residual[valid])
    g = g if np.isfinite(g) and g > 0 else 1.0
    mult = np.ones_like(residual)
    for _ch, lo, hi in LSF_CHANNELS:
        use = valid & _channel_mask(wave, lo, hi)
        if not use.any():
            continue
        sig = _robust_mad(residual[use])
        sig = sig if np.isfinite(sig) and sig > 0 else g
        a = np.abs(residual[use]); thr = k * sig
        h = np.ones_like(a); out = a > thr; h[out] = thr / a[out]
        mult[use] *= h
    w = ivar * mult
    for _ch, lo, hi in LSF_CHANNELS:
        use = np.isfinite(w) & (w > 0.0) & _channel_mask(wave, lo, hi)
        if use.any():
            w[use] /= float(np.mean(w[use]))
    w[~np.isfinite(w) | (w < 0.0)] = 0.0
    return w


def test_varying_ivar_separates_it_from_the_raw_residual_form():
    """A bright pixel with a large residual but correspondingly large noise is
    an outlier under the old test and not under the current one."""
    rng = np.random.default_rng(1)
    wave, skyline = _grid()
    sigma = np.linspace(1.0, 10.0, N)        # noise varies 10x, as photon weighting does
    resid = rng.normal(size=N) * sigma
    ivar = 1.0 / sigma ** 2
    w_new, noise_new = continuum_fit_weights(wave, resid, ivar, skyline)
    w_old = _legacy_weights(wave, resid, ivar, skyline)
    assert not np.allclose(w_new, w_old)
    # Residuals are drawn AT their own sigma, so in noise units the MAD is ~1
    # and almost nothing should be flagged.
    for ch, sd in noise_new.items():
        assert 0.5 < sd < 2.0, f"channel {ch} noise-unit MAD {sd} is not ~1"


def test_it_stops_penalising_the_noisy_end():
    """The pixels the old form down-weighted are the HIGH-NOISE ones, which is
    backwards: they are not discrepant, merely uncertain."""
    rng = np.random.default_rng(2)
    wave, skyline = _grid()
    sigma = np.linspace(1.0, 10.0, N)
    resid = rng.normal(size=N) * sigma
    ivar = 1.0 / sigma ** 2
    w_new, _ = continuum_fit_weights(wave, resid, ivar, skyline)
    w_old = _legacy_weights(wave, resid, ivar, skyline)
    noisy = sigma > np.percentile(sigma, 90)
    quiet = sigma < np.percentile(sigma, 10)
    ratio_noisy = np.median(w_new[noisy]) / np.median(w_old[noisy])
    ratio_quiet = np.median(w_new[quiet]) / np.median(w_old[quiet])
    assert ratio_noisy > ratio_quiet, (
        f"expected weight restored at the noisy end relative to the quiet "
        f"end, got {ratio_noisy:.4f} vs {ratio_quiet:.4f}")


def test_a_genuine_outlier_is_still_down_weighted():
    """The point is to change WHICH pixels are outliers, not to stop clipping."""
    rng = np.random.default_rng(3)
    wave, skyline = _grid()
    sigma = np.full(N, 2.0)
    resid = rng.normal(size=N) * sigma
    resid[N // 2] = 60.0 * sigma[N // 2]          # 60 sigma
    ivar = 1.0 / sigma ** 2
    w_new, _ = continuum_fit_weights(wave, resid, ivar, np.zeros(N, bool))
    assert w_new[N // 2] < 0.1 * np.median(w_new)
