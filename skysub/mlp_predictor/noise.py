"""Photon-noise model for the per-pixel flux-space loss.

The spectra carry no per-pixel error array, so until now the flux-space term in
``trainer.compressed_loss`` weighted every wavelength equally.  That is wrong by
a large factor: the instrument throughput varies by 4.3x across the b channel
alone, so a fixed flux error is a very different number of photons at 3600 A
than at 5000 A, and the unweighted loss was fitting the poorly-measured blue
end with the same authority as the well-measured middle.

The model
---------
``mean-sens-{b,r,z}-v1.1.csv`` in ``$LVMCORE_DIR/sensitivity`` give the mean
standard-star sensitivity per channel: flux per count, i.e. the inverse of the
total throughput (atmosphere x telescope x instrument).  Then

    counts(lambda)  = flux(lambda) / sens(lambda)
    sigma_counts    = sqrt(counts)                      (Poisson)
    sigma_flux      = sigma_counts * sens = sqrt(flux(lambda) * sens(lambda))

so the inverse-variance pixel weight is ``1 / (flux * sens)``.  ``flux`` is the
TOTAL observed flux in the fibre -- the whole photon budget sets the noise, not
the component being fitted -- and a global constant on ``sens`` (exposure time,
fibre count in the median stack, gain) rescales every weight together and drops
out once the weights are normalised.

The one thing that does NOT drop out: the per-arm normalisation
--------------------------------------------------------------
``avgsens-1-1.py`` normalises each channel's curve to unit mean inside a
Gaussian band at 4500 / 6500 / 8500 A BEFORE taking the median over exposures,
so the three published curves are shapes on three independent scales.  Joining
them takes two numbers and they matter: naive concatenation puts factor-2 and
factor-6 steps at the arm joins straight into the noise model.

``SENS_ARM_SCALE`` comes from the total-throughput curve (including atmosphere)
at the three band centres, ~0.375 / 0.50 / 0.545, because a published
throughput is a direct measurement.  Three other routes were tried; they
disagree, and the disagreement is the honest uncertainty on this model:

* Matching the arms in their OVERLAPS (5775-5800, 7520-7570 A) gives r = 0.520,
  z = 0.158.  REJECTED: the overlaps sit at the extreme detector edges where
  the two arms have opposite slopes and disagree in SHAPE -- the b/r ratio
  swings 0.44-0.63 across 51 pixels -- and on a balanced set of line-free
  windows this leaves a 3.22x arm-to-arm residual against 1.36x for the values
  adopted here.
* A low envelope of measured ``sigma^2/flux`` over 20 A windows: r = 0.539,
  z = 0.442.
* The fit RESIDUALS: r = 0.679, z = 4.474, and unusable -- in the NIR they are
  dominated by model misfit on the bright OH lines rather than photon noise
  (log-scatter 88-182% against 23% for the direct estimator).

No noise-based route can settle it, because they all rest on a
second-difference estimator that assumes WHITE noise while the DRP resamples
onto the common 0.5 A grid -- and the LSF is 3.13 / 2.62 / 3.50 pixels FWHM in
b / r / z, so the pixel-to-pixel correlation, and hence the estimator's bias,
is arm-dependent in precisely the comparison being made.  So treat the
arm-to-arm weighting as good to ~1.5x, and the WITHIN-arm shape -- which
carries the 4.3x dynamic range and comes from hundreds of standard-star
exposures -- as good.  Exact band-centre throughputs would go straight into
SENS_ARM_SCALE and nothing else would change.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

# (name, blue edge, red edge) of each published curve, in Angstrom.
SENS_ARMS = (("b", 3600.0, 5800.0), ("r", 5775.0, 7570.0), ("z", 7520.0, 9800.0))

# Per-arm normalisation relative to b; undoes the independent unit-mean
# normalisation each published curve carries.  See the module docstring.
SENS_ARM_SCALE = {"b": 1.0, "r": 0.74, "z": 0.69}

SENS_FILENAME = "mean-sens-{arm}-v1.1.csv"


def sensitivity_dir(sens_dir=None):
    """Directory holding the mean sensitivity CSVs."""
    if sens_dir is not None:
        return Path(sens_dir).expanduser()
    root = os.environ.get("LVMCORE_DIR")
    if not root:
        raise RuntimeError(
            "LVMCORE_DIR is not set and sens_dir was not given; cannot locate "
            "mean-sens-{b,r,z}-v1.1.csv")
    return Path(root).expanduser() / "sensitivity"


def load_relative_sensitivity(wave, sens_dir=None, arm_scale=None, verbose=False):
    """Joined relative sensitivity (flux per count, arbitrary global scale).

    Each published curve is interpolated onto ``wave``, scaled by its
    ``SENS_ARM_SCALE`` factor, and CROSSFADED linearly across the overlap so the
    join is continuous.  A plain average over the overlap leaves a visible step
    (0.10 and 0.056, against a typical adjacent-pixel change of 0.00015)
    because the two arms cross there with opposite slopes; the crossfade also
    puts the weight on whichever arm is further from its own edge, which is the
    one worth trusting.

    Raises if any pixel ends up uncovered: a silent NaN here becomes a NaN
    pixel weight and a NaN loss.
    """
    wave = np.asarray(wave, dtype=np.float64)
    scale = dict(SENS_ARM_SCALE if arm_scale is None else arm_scale)
    directory = sensitivity_dir(sens_dir)
    num = np.zeros_like(wave)
    den = np.zeros_like(wave)
    for arm, _lo, _hi in SENS_ARMS:
        table = np.loadtxt(directory / SENS_FILENAME.format(arm=arm), delimiter=",")
        w_a, s_a = table[:, 0], table[:, 1]
        inside = (wave >= w_a[0]) & (wave <= w_a[-1])
        if not np.any(inside):
            continue
        lam = wave[inside]
        vals = np.interp(lam, w_a, s_a) * float(scale[arm])
        weight = np.ones(lam.size, dtype=np.float64)
        for other, olo, ohi in SENS_ARMS:
            if other == arm:
                continue
            a, b = max(w_a[0], olo), min(w_a[-1], ohi)
            if b <= a:
                continue
            frac = np.clip((lam - a) / (b - a), 0.0, 1.0)
            # Ramp up away from this arm's own edge: if the overlap starts at
            # this arm's blue end the arm is untrustworthy at `a`, else at `b`.
            ramp = frac if a <= w_a[0] + 1e-9 else 1.0 - frac
            weight = np.where((lam >= a) & (lam <= b), ramp * weight, weight)
        num[inside] += vals * weight
        den[inside] += weight
    bad = ~(den > 0) | ~np.isfinite(num)
    if np.any(bad):
        raise RuntimeError(
            f"relative sensitivity is undefined on {int(bad.sum())} of "
            f"{wave.size} pixels (first at {wave[np.flatnonzero(bad)[0]]:.1f} A); "
            f"the published curves do not cover this wavelength grid")
    sens = num / den
    if verbose:
        print(f"  [noise] relative sensitivity {sens.min():.4g}-{sens.max():.4g} "
              f"over {wave[0]:.0f}-{wave[-1]:.0f} A; arm scales "
              + ", ".join(f"{k}={v:g}" for k, v in scale.items()))
    return sens


def photon_pixel_weight(flux, sens, floor_frac=0.05):
    """Inverse-variance pixel weights ``1/(flux*sens)``, row-normalised to mean 1.

    ``flux`` is (n_row, n_pix) TOTAL observed flux; ``sens`` is (n_pix,).  The
    variance is floored at ``floor_frac`` times the row's median variance before
    inverting, which does two jobs: it stops pixels whose observed flux went
    non-positive on noise from taking infinite weight, and it stands in for the
    read-noise term this model does not have (between the OH bands in the NIR
    the sky is faint enough for read noise to matter).

    Normalising each ROW to mean weight 1 is deliberate.  It keeps the relative
    weighting across wavelength, which is the point, while leaving the relative
    weighting BETWEEN rows exactly as it was: row-to-row noise differences are
    real, but they also track sky brightness and exposure depth, and letting
    them in here would silently re-weight the training set on top of the
    existing per-row weights.
    """
    flux = np.asarray(flux, dtype=np.float64)
    sens = np.asarray(sens, dtype=np.float64)
    if flux.ndim != 2 or sens.shape != (flux.shape[1],):
        raise ValueError(
            f"shape mismatch: flux {flux.shape} against sens {sens.shape}")
    var = flux * sens[None, :]
    var = np.where(np.isfinite(var) & (var > 0.0), var, np.nan)
    med = np.nanmedian(var, axis=1, keepdims=True)
    med = np.where(np.isfinite(med) & (med > 0.0), med, 1.0)
    var = np.where(np.isfinite(var), var, med)
    var = np.maximum(var, float(floor_frac) * med)
    w = 1.0 / var
    w /= np.mean(w, axis=1, keepdims=True)
    return w.astype(np.float32)


__all__ = [
    "SENS_ARMS",
    "SENS_ARM_SCALE",
    "load_relative_sensitivity",
    "photon_pixel_weight",
    "sensitivity_dir",
]
