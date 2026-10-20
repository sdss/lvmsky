"""Correct a predicted science sky by the sky arms' own decomposition residual.

Why this works
--------------
The decomposition leaves a residual in every arm that is mostly NOT noise. In
the blue it is a fixed pattern of solar absorption lines (Balmer series, Ca H&K,
the G band) that the moon and zodi basis does not reproduce: on bright-continuum
rows the science arm's fractional residual correlates at 0.98 with one fixed
template, with an RMS of ~10% of the continuum. The science prediction inherits
that residual, because the ML reproduces the science decomposition's
coefficients and not the science spectrum. Both sky arms carry the same pattern
at the same time, so their residual is a direct, simultaneous measurement of
what the prediction is missing.

The estimator
-------------
For each sky arm ``a`` (near, far):

    e_a(lambda) = [obs_a - model_a](lambda) * k_a(lambda)
    k_a         = M_sci(lambda) / M_a(lambda)       (scale='model', the default)

``M_sci`` is the PREDICTED science sky model and ``M_a`` the arm's own model,
both total (every family), taken pixel by pixel and clipped to [0.2, 5]. The
residual at a pixel scales with whatever light makes up that pixel: sunlight in
the blue continuum, airglow in the diffuse, OH in a line core. The total-model
ratio follows each of them, so it needs no choice of family. Measured on 457
uncrowded held-out rows (full-band correction, median chi2):

    k = 1            0.287      moon+zodi ratio          0.276
    model ratio      0.269      on OH-line pixels: model 0.848, moon+zodi 0.934

In the blue the model ratio and the moon+zodi ratio are equivalent (0.148 both).
``scale='solar'`` keeps the moon+zodi ratio, smoothed over ``ratio_smoothing``
Angstrom. The arms are then combined with inverse-variance weights, using each
arm's photon-noise variance times ``k_a**2``.

A narrow feature in ONE sky field (nebular emission, a star, a cosmic) is not
part of the shared pattern, but the weighted mean would carry about half of it
to the science.  Where the two arms' scaled residuals disagree by more than
``consistency_clip`` (5) sigma, the arm with the smaller residual is used
instead.  Measured on 2000 rows (solarblue corpus, a FWHM 2.5 A line injected
into the near arm at 0.3/1/3x the continuum): the transferred fraction falls
from 48% to 7/1.4/0.3%, for blue chi2 0.129 -> 0.138 (1.8% of pixels clipped).
Smoothing does much worse for the same purpose: an 8 A running median still
transfers 15-25% and raises the chi2 to 0.26.

When ``max_wavelength`` cuts the band short, the correction is applied in full
up to it and then tapered to zero with a cosine over the next ``taper_width``
Angstrom, so that no step is left in the corrected spectrum. The ramp sits
OUTSIDE the requested range on purpose: a ramp ending at the boundary halves
the correction inside it (4900-5000 A: chi2 0.079 -> 0.251), while the
correction just redward is as good as anywhere (5000-5100 A: 0.354 -> 0.078
uncorrected -> corrected), so the ramp costs nothing there.

Two safeguards, both measured as necessary on held-out rows:

* An arm whose decomposition failed (blue fractional residual RMS above
  ``gate``; the median is 0.064, failures sit near 1.0) is dropped. Without
  this gate one bad far-arm fit raises the mean chi2 by x40.
* On faint rows the arms' noise competes with the pattern. A 3-pixel running
  median helps there, but costs on bright rows (x0.007 -> x0.06). So
  ``smoothing='auto'`` smooths only rows whose solar-continuum S/N is below
  ``auto_snr``.

Measured on the held-out split of the palacecorr corpus (uncrowded fields, test
rows, blue < 5000 A, median single-fibre photon chi2): uncorrected ML 1.140;
near arm only and unsmoothed 0.236; near arm with a 10 A running median 0.915;
this estimator 0.146 (with moon+zodi scaling; the model ratio gives the same
in the blue). The decomposition's own floor is 0.966.

The full band
-------------
Applied redward of 5000 A as it stands (total-model scaling, OH included, no
high-pass), the same estimator removes most of what the prediction misses in
R and Z as well: the moonlight's telluric bands, the solar pattern, and the OH
residual the line scaling cannot reach.  Held-out gaia-stars-solarblue rows
(2852, OH line scaling applied after it, median single-fibre photon chi2):

                      green   red    z      red OH px   z OH px
    5000 A (before)   0.139   0.244  0.366  0.493       0.234
    full band         0.077   0.111  0.125  0.064       0.075

91% of rows improve.  Two alternatives lose: a 50 A high-pass on the arm
residual (z 0.285, the low frequencies ARE shared with the prediction), and
leaving the OH lines out of the transfer (no gain on the continuum, OH pixels
x4 worse).  The OH band scales fitted afterwards move by a median of 0.005,
so the two corrections are complementary.

Protecting science lines
------------------------
Emission common to BOTH sky fields -- diffuse Galactic H-alpha, [N II],
[S II], [S III], [O III], [O II] -- passes the consistency clip, and the
transfer carries it to the science at 100-115% (the full band and the 5000 A
band alike).  ``line_masks`` (from :func:`line_protection_masks`) suppresses
it inside a window around every line of ``TRANSFER_LINES`` --
+/-max(500 km/s, 3 FWHM) redward of 5800 A, +/-max(300 km/s, 2 FWHM)
blueward -- without simply dropping the correction there (chi2 inside the
windows would go back from 0.093 to 0.236):

* blueward of ``line_split_wavelength`` (5800 A) the correction is
  interpolated across each window between robust anchors, which keeps the
  slowly varying solar pattern.  In-window chi2, scored on the same +/-500
  km/s pixels: 0.240 dropped, 0.174 interpolated at 500 km/s, 0.114 at 300
  km/s, against 0.076 unmasked (green 0.149 / 0.080 / 0.070 / 0.064).  200
  and 150 km/s reach 0.099 and 0.093 but let a line at +/-150 km/s onto the
  anchors (p90 23% and 67% transferred), 300 km/s keeps that at 0.2%;
* redward, inside a window only the OH-dominated part is kept, weighted by the
  science OH fraction ``sci_oh / sci_model``, and every line core
  (+/-(150 km/s + 1.5 FWHM)) is zeroed (red 0.204, z 0.245, against 0.235 /
  0.309 with the window zeroed).  Interpolating there is worse than zeroing:
  the windows hold OH and [O I] airglow.

Injected into both arms at rest or +/-150 km/s, none of 18 lines (including
[O I] 6300 on the airglow line and [S III] 9069/9531 on OH) is transferred,
maximum 0.00 over 400 rows.  Beyond the cores the red windows pass lines on OH
(at +260 km/s up to 25% of H-alpha and [O I] 6300).  Per-row velocity
detection in the arm residuals is not usable for this -- OH residual under
H-alpha exceeds 5 sigma on every row -- so the cores are centred at rest and,
for a pointing inside the Magellanic Clouds, also at their systemic velocity
(:func:`protection_velocities`).  ``line_mode='windows'`` zeroes every red
window instead, velocity-proof to the window edge.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Union

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d

__all__ = ["SkyArm", "arm_photon_variance", "sky_arm_residual_correction", "TRANSFER_LINES",
           "line_protection_masks", "protect_science_lines", "protection_velocities"]

_C_KM_S = 299792.458

# Every nebular line LVM measures, at rest (air): the decomposition's science
# mask, the red lines of sky_line_scaling.NEBULAR_LINES_RED, and the blue and
# green lines neither carries.
TRANSFER_LINES = (
    ("H12", 3750.15), ("H11", 3770.63), ("H10", 3797.90), ("H9", 3835.39),
    ("[OII]3726", 3726.03), ("[OII]3729", 3728.82), ("[NeIII]3869", 3868.76), ("H8+HeI", 3888.65),
    ("[NeIII]3967", 3967.47), ("Heps", 3970.07), ("[SII]4069", 4068.60), ("[SII]4076", 4076.35),
    ("Hdelta", 4101.74), ("Hgamma", 4340.47), ("[OIII]4363", 4363.21), ("HeI4471", 4471.48),
    ("HeII4686", 4685.68), ("[ArIV]4711", 4711.37), ("[ArIV]4740", 4740.16), ("Hbeta", 4861.33),
    ("HeI4922", 4921.93), ("[OIII]4959", 4958.91), ("[OIII]5007", 5006.84), ("HeI5016", 5015.68),
    ("[NI]5198", 5197.90), ("[NI]5200", 5200.26), ("[ClIII]5518", 5517.71), ("[ClIII]5538", 5537.88),
    ("[NII]5755", 5754.64), ("HeI5876", 5875.62), ("[OI]6300", 6300.30), ("[SIII]6312", 6312.06),
    ("[OI]6364", 6363.78), ("[NII]6548", 6548.05), ("Halpha", 6562.80), ("[NII]6583", 6583.45),
    ("HeI6678", 6678.15), ("[SII]6716", 6716.44), ("[SII]6731", 6730.82), ("HeI7065", 7065.19),
    ("[ArIII]7136", 7135.79), ("[FeII]7155", 7155.16), ("HeI7281", 7281.35), ("[OII]7320", 7319.99),
    ("[OII]7330", 7330.20), ("[NiII]7378", 7377.83), ("[ArIII]7751", 7751.06), ("P16", 8502.48),
    ("P15", 8545.38), ("P14", 8598.39), ("[FeII]8617", 8616.95), ("P13", 8665.02),
    ("P12", 8750.47), ("P11", 8862.78), ("P10", 9014.91), ("[SIII]9069", 9068.60),
    ("P9", 9229.01), ("[SIII]9531", 9530.60), ("P8", 9545.97),
)
# (name, RA, Dec, radius in deg, systemic velocity in km/s); footprints as
# mlp_predictor.data.LMC_EXCLUSION / SMC_EXCLUSION.
MAGELLANIC_CLOUDS = (("LMC", 81.0, -69.7, 8.0, 262.0), ("SMC", 14.0, -73.0, 4.0, 146.0))
LINE_WINDOW_KM_S = 500.0
LINE_WINDOW_FWHM = 3.0
# Blueward of LINE_SPLIT_WAVELENGTH the windows are interpolated across, so a
# narrower window keeps more of the solar pattern; 300 km/s is the narrowest
# that still keeps a line at +/-150 km/s off the interpolation anchors.
LINE_WINDOW_BLUE_KM_S = 300.0
LINE_WINDOW_BLUE_FWHM = 2.0
LINE_SPLIT_WAVELENGTH = 5800.0
LINE_CORE_KM_S = 150.0
LINE_CORE_FWHM = 1.5


@dataclass
class SkyArm:
    """One sky arm's inputs; every array is (n_wave,) or (n_rows, n_wave).

    observed  the arm's observed spectrum.
    model     the arm's reconstruction from its OWN decomposition coefficients.
    solar     the solar-continuum part of that model: moon + zodi.
    variance  the photon-noise variance of ``observed`` (see
              :func:`arm_photon_variance`).
    """
    observed: np.ndarray
    model: np.ndarray
    solar: np.ndarray
    variance: np.ndarray


def arm_photon_variance(flux, wave, n_fibres=None, exptime=900.0, sens=None):
    """Photon variance of a sky-arm median stack, with the shared floor.

    ``n_fibres`` is META's ``fibers_sky_{near,far}_used``. Without it the
    variance is the single-fibre one, which over-states the arm's noise by the
    ~50 stacked fibres and so under-weights the arm.  ``sens`` is the
    absolute sensitivity on ``wave``; pass it when calling per row to avoid
    re-reading the sensitivity tables.
    """
    from mlp_predictor.noise import (floor_variance, load_absolute_sensitivity,
                                     photon_variance_absolute)
    wave = np.asarray(wave, dtype=np.float64)
    flux2 = np.atleast_2d(np.asarray(flux, dtype=np.float64))
    sens = np.asarray(load_absolute_sensitivity(wave) if sens is None else sens,
                      dtype=np.float64)
    nf = None if n_fibres is None else np.atleast_1d(np.asarray(n_fibres, dtype=np.float64))
    var = floor_variance(photon_variance_absolute(
        np.clip(np.nan_to_num(flux2), 0.0, None), sens, exptime=float(exptime),
        dwave=float(np.median(np.diff(wave))), n_fibres=nf))
    return var if np.ndim(flux) > 1 else var[0]


def _width_px(wave, width_a):
    return max(3, int(round(float(width_a) / float(np.median(np.diff(wave))))) | 1)


def _running_median(x, wave, width_a):
    return median_filter(x, size=(1, _width_px(wave, width_a)), mode="nearest")


def _running_mean(x, wave, width_a):
    return uniform_filter1d(x, _width_px(wave, width_a), axis=-1, mode="nearest")


def _taper(wave, max_wavelength, taper_width):
    """1 up to ``max``, cosine down to 0 at ``max + width``, 0 beyond.

    All ones when ``max_wavelength`` covers the whole grid.
    """
    lmax = float(max_wavelength)
    if lmax >= wave.max():
        return np.ones_like(wave)
    w = max(float(taper_width), 0.0)
    if w == 0.0:
        return (wave < lmax).astype(np.float64)
    x = np.clip((lmax + w - wave) / w, 0.0, 1.0)    # 1 at the boundary, 0 a width out
    return np.sin(0.5 * np.pi * x) ** 2


def protection_velocities(ra_deg, dec_deg):
    """Systemic velocities (km/s) of the Magellanic Clouds whose footprint holds the pointing."""
    ra, dec = np.radians(float(ra_deg)), np.radians(float(dec_deg))
    out = []
    for _, ra0, dec0, radius, v in MAGELLANIC_CLOUDS:
        r0, d0 = np.radians(ra0), np.radians(dec0)
        cos = np.sin(dec) * np.sin(d0) + np.cos(dec) * np.cos(d0) * np.cos(ra - r0)
        if np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))) <= radius:
            out.append(v)
    return tuple(out)


def _line_fwhm(wave, fwhm, lam):
    if fwhm is None:
        return 0.0
    fwhm = np.broadcast_to(np.asarray(fwhm, dtype=np.float64), wave.shape)
    use = (np.abs(wave - lam) <= 25.0) & np.isfinite(fwhm) & (fwhm > 0)
    return float(np.median(fwhm[use])) if use.any() else 0.0


def line_protection_masks(wave, fwhm=None, velocities=(), lines=TRANSFER_LINES, *,
                          window_km_s=LINE_WINDOW_KM_S, window_fwhm=LINE_WINDOW_FWHM,
                          blue_window_km_s=LINE_WINDOW_BLUE_KM_S, blue_window_fwhm=LINE_WINDOW_BLUE_FWHM,
                          split_wavelength=LINE_SPLIT_WAVELENGTH,
                          core_km_s=LINE_CORE_KM_S, core_fwhm=LINE_CORE_FWHM):
    """``(windows, cores)`` boolean masks over ``wave`` for :func:`protect_science_lines`.

    Windows are +/-max(window_km_s, window_fwhm x FWHM) around every line,
    +/-max(blue_window_km_s, blue_window_fwhm x FWHM) for a line blueward of
    ``split_wavelength``, and cores +/-(core_km_s + core_fwhm x FWHM); all at
    rest and shifted to each of ``velocities``.  ``fwhm`` is the detector FWHM
    in Angstrom, per pixel or scalar.
    """
    wave = np.asarray(wave, dtype=np.float64)
    windows = np.zeros(wave.shape, dtype=bool)
    cores = np.zeros(wave.shape, dtype=bool)
    for _, lam in lines:
        f = _line_fwhm(wave, fwhm, lam)
        kms, kf = (blue_window_km_s, blue_window_fwhm) if lam < split_wavelength else (window_km_s, window_fwhm)
        w = max(lam * kms / _C_KM_S, kf * f)
        c = lam * core_km_s / _C_KM_S + core_fwhm * f
        for v in {0.0, *map(float, velocities)}:
            centre = lam * (1.0 + v / _C_KM_S)
            windows |= np.abs(wave - centre) <= w
            cores |= np.abs(wave - centre) <= c
    return windows, cores & windows


def _interpolate_across(y, mask, skyline, edge_px=8, reach_px=20):
    """A line across each masked run between the medians of clean pixels beside it.

    Anchors use up to ``edge_px`` pixels within ``reach_px`` of the run that are
    neither masked nor dominated by a sky line; a run with none on either side
    is set to zero.
    """
    y = np.array(y, dtype=np.float64)
    n = y.size
    clean = ~mask & ~skyline & np.isfinite(y)
    edges = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    for a, b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        lo = max(0, a - reach_px)
        left = np.flatnonzero(clean[lo:a])[-edge_px:] + lo
        right = np.flatnonzero(clean[b:min(n, b + reach_px)])[:edge_px] + b
        if not left.size and not right.size:
            y[a:b] = 0.0
            continue
        lv = np.median(y[left]) if left.size else np.median(y[right])
        rv = np.median(y[right]) if right.size else lv
        lx = left.mean() if left.size else a - 1.0
        rx = right.mean() if right.size else float(b)
        y[a:b] = lv + (rv - lv) * (np.arange(a, b) - lx) / max(rx - lx, 1.0)
    return y


def protect_science_lines(correction, wave, windows, cores, sci_oh_fraction, *, mode="cores",
                          split_wavelength=LINE_SPLIT_WAVELENGTH, skyline_fraction=0.3):
    """The correction with what the sky arms' emission lines could carry removed.

    Blueward of ``split_wavelength`` each window is interpolated across.
    Redward, ``mode='cores'`` keeps the correction inside a window only in
    proportion to ``sci_oh_fraction`` and zeroes the cores; ``mode='windows'``
    zeroes the whole window.  Rows of a 2-D ``correction`` are treated
    independently, ``sci_oh_fraction`` having the same shape.
    """
    if mode not in ("cores", "windows"):
        raise ValueError(f"mode must be 'cores' or 'windows'; got {mode!r}")
    wave = np.asarray(wave, dtype=np.float64)
    corr = np.atleast_2d(np.asarray(correction, dtype=np.float64))
    frac = np.broadcast_to(np.atleast_2d(np.asarray(sci_oh_fraction, dtype=np.float64)), corr.shape)
    blue = wave < float(split_wavelength)
    out = np.empty_like(corr)
    for i in range(corr.shape[0]):
        f = np.clip(np.nan_to_num(frac[i]), 0.0, 1.0)
        interp = _interpolate_across(corr[i], windows, f > float(skyline_fraction))
        if mode == "cores":
            red = np.where(cores, 0.0, np.where(windows, corr[i] * f, corr[i]))
        else:
            red = np.where(windows, 0.0, corr[i])
        out[i] = np.where(blue, interp, red)
    return out[0] if np.ndim(correction) == 1 else out


def sky_arm_residual_correction(
    wave,
    sci_model,
    sci_solar,
    arms: Sequence[SkyArm],
    *,
    scale: Optional[str] = "model",
    max_wavelength: float = 5000.0,
    taper_width: float = 100.0,
    smoothing: Union[str, float, None] = "auto",
    auto_snr: float = 100.0,
    auto_width: float = 1.5,
    gate: float = 0.25,
    gate_max_wavelength: float = 5000.0,
    consistency_clip: Optional[float] = 5.0,
    ratio_smoothing: float = 100.0,
    line_masks=None,
    sci_oh=None,
    line_mode: str = "cores",
    line_split_wavelength: float = LINE_SPLIT_WAVELENGTH,
    return_info: bool = False,
):
    """Correction to SUBTRACT from (science observed - predicted science sky).

    Equivalently, ADD it to the predicted sky. The returned array has the shape
    of ``sci_model``.

    Parameters
    ----------
    wave : (n_wave,) wavelength grid in Angstrom, shared by every array.
    sci_model : the PREDICTED science sky, total (every family). It sets
        ``k_a`` when ``scale='model'``.
    sci_solar : the PREDICTED science moon + zodi. It sets the row's S/N for
        ``smoothing='auto'``, and ``k_a`` when ``scale='solar'``. Using the
        decomposition's science instead of the prediction changes the result by
        under 1%, so the prediction's errors do not matter here.
    arms : the sky arms to combine, usually ``(near, far)``. One arm works too.
    scale : ``'model'`` (total-model ratio per pixel), ``'solar'`` (smoothed
        moon + zodi ratio), or None (each arm's residual as it stands).
    max_wavelength : the correction is applied in full up to here. The full
        band also helps (median chi2 1.27 -> 0.27 on held-out rows); pass
        ``np.inf`` for it.
    taper_width : width in Angstrom of the cosine taper that starts at
        ``max_wavelength`` and reaches zero ``taper_width`` redward of it.
        Unused when ``max_wavelength`` covers the grid.
    smoothing : ``'auto'`` (3-pixel median below ``auto_snr``), a width in
        Angstrom (running median on every row), or 0 / None (unsmoothed).
    consistency_clip : where both arms are used and their scaled residuals
        differ by more than this many sigma (their combined photon noise),
        take the arm whose residual is SMALLER instead of the weighted mean.
        The solar-line pattern is common to both arms; a narrow feature in
        one sky field (emission, a star, a cosmic) is not, and would otherwise
        be transferred to the science with that arm's weight (~half).  None
        disables it.
    gate : drop an arm whose blue fractional residual RMS,
        ``|obs - model| / solar`` below ``gate_max_wavelength``, exceeds this.
    line_masks : ``(windows, cores)`` from :func:`line_protection_masks`, or
        None (default) for no science-line protection.  Applied last, by
        :func:`protect_science_lines`.
    sci_oh : the PREDICTED science OH component, shaped like ``sci_model``;
        required with ``line_masks``.
    line_mode, line_split_wavelength : passed to :func:`protect_science_lines`.
    return_info : also return a dict with per-arm usage, row S/N, smoothed rows.
    """
    wave = np.asarray(wave, dtype=np.float64)
    one_row = np.ndim(sci_model) == 1
    m_sci = np.atleast_2d(np.asarray(sci_model, dtype=np.float64))
    n_rows, n_wave = m_sci.shape
    s_sci = np.broadcast_to(np.atleast_2d(np.asarray(sci_solar, dtype=np.float64)),
                            (n_rows, n_wave))
    if wave.shape != (n_wave,):
        raise ValueError(f"wave {wave.shape} does not match sci_model {m_sci.shape}")
    if not arms:
        raise ValueError("need at least one sky arm")
    if scale not in ("model", "solar", None):
        raise ValueError(f"scale must be 'model', 'solar' or None; got {scale!r}")
    gate_band = wave < float(gate_max_wavelength)
    if scale == "solar":
        s_sci_smooth = _running_mean(np.nan_to_num(s_sci), wave, ratio_smoothing)

    num = np.zeros((n_rows, n_wave))
    den = np.zeros((n_rows, n_wave))
    used, frac_rms = [], []
    per_arm = []                                  # (estimate, variance) per arm, NaN where unusable
    for arm in arms:
        obs, mod, sol, var = (np.broadcast_to(np.atleast_2d(np.asarray(x, dtype=np.float64)),
                                              (n_rows, n_wave))
                              for x in (arm.observed, arm.model, arm.solar, arm.variance))
        resid = obs - mod
        with np.errstate(divide="ignore", invalid="ignore"):
            fr = np.sqrt(np.nanmedian(np.where(gate_band, (resid / sol) ** 2, np.nan), axis=1))
            if scale == "model":
                k = m_sci / mod
                k = np.where(np.isfinite(k) & (m_sci > 0) & (mod > 0),
                             np.clip(k, 0.2, 5.0), 1.0)
            elif scale == "solar":
                sol_smooth = _running_mean(np.nan_to_num(sol), wave, ratio_smoothing)
                k = np.where(sol_smooth > 0, s_sci_smooth / sol_smooth, np.nan)
            else:
                k = np.ones_like(resid)
            ok_row = np.isfinite(fr) & (fr < float(gate))
            ok = ok_row[:, None] & np.isfinite(resid) & np.isfinite(k) & np.isfinite(var) & (var > 0)
            w = np.where(ok, 1.0 / np.where(ok, var * k ** 2, 1.0), 0.0)
        num += w * np.where(ok, resid * k, 0.0)
        den += w
        per_arm.append((np.where(ok, resid * k, np.nan), np.where(ok, var * k ** 2, np.nan)))
        used.append(ok_row)
        frac_rms.append(fr)

    corr = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
    clipped = np.zeros((n_rows, n_wave), dtype=bool)
    if consistency_clip is not None and len(per_arm) == 2:
        (e1, v1), (e2, v2) = per_arm
        with np.errstate(invalid="ignore"):
            clipped = np.abs(e1 - e2) > float(consistency_clip) * np.sqrt(v1 + v2)
        corr = np.where(clipped, np.where(np.abs(e1) < np.abs(e2), e1, e2), corr)
    with np.errstate(invalid="ignore"):
        snr = np.nanmedian(np.where(gate_band, s_sci * np.sqrt(den), np.nan), axis=1)
    smoothed = np.zeros(n_rows, dtype=bool)
    if isinstance(smoothing, str):
        if smoothing != "auto":
            raise ValueError(f"smoothing must be 'auto', a width in A, or None; got {smoothing!r}")
        smoothed = ~(snr >= float(auto_snr))       # NaN S/N (no usable arm) counts as faint
        if smoothed.any():
            corr[smoothed] = _running_median(corr[smoothed], wave, auto_width)
    elif smoothing:
        smoothed[:] = True
        corr = _running_median(corr, wave, float(smoothing))
    corr = corr * _taper(wave, max_wavelength, taper_width)[None, :]
    if line_masks is not None:
        if sci_oh is None:
            raise ValueError("line_masks needs sci_oh, the predicted science OH component")
        with np.errstate(divide="ignore", invalid="ignore"):
            oh_frac = np.broadcast_to(np.atleast_2d(np.asarray(sci_oh, dtype=np.float64)), m_sci.shape) / m_sci
        corr = protect_science_lines(corr, wave, *line_masks, np.where(m_sci > 0, oh_frac, 0.0),
                                     mode=line_mode, split_wavelength=line_split_wavelength)

    out = corr[0] if one_row else corr
    if not return_info:
        return out
    info = dict(arm_used=np.stack(used, axis=0), arm_frac_rms=np.stack(frac_rms, axis=0),
                snr=snr, smoothed=smoothed, clipped_fraction=clipped.mean(axis=1))
    if one_row:
        info = {k: (v[:, 0] if v.ndim == 2 else v[0]) for k, v in info.items()}
    return out, info
