#!/usr/bin/env python3

"""
Run sky spectral decomposition on a median-stacked LVM frame.

Usage:
    python decompose_parallel.py <data_file> [palace_dir] [options]

Example:
    python decompose_parallel.py lvmsframe_median_stack.fits ../ --n-workers 8
"""

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
# clarabel (Rust QP solver) and any other Rust/Rayon library ignore OMP_NUM_THREADS.
os.environ["RAYON_NUM_THREADS"] = "1"
os.environ["POLARS_MAX_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["TBB_NUM_THREADS"] = "1"

import argparse
import queue as queue_mod
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import multiprocessing as mp

import numpy as np
from astropy.io import fits
from tqdm import tqdm

# ``python /path/to/skysub/decompose_parallel.py`` puts only ``skysub/`` on
# sys.path.  Add the repository/package root so direct-script and ``-m``
# execution use the same fully-qualified package imports.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from skysub.sky_decomp.result_io import results_to_fits
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT as DEFAULT_MOON_ZODI_DATA_ROOT,
    validate_decomposition_data_root,
)

try:
    from threadpoolctl import threadpool_limits
except ImportError:  # threadpoolctl is optional; env vars are the fallback.
    threadpool_limits = None


def _clamp_native_threads(n=1):
    """Force every loaded thread pool (BLAS/OpenMP/Rayon/TBB/etc.) to `n` threads."""
    # Redundant with the env vars but catches lazy imports and fork-inherited pools.
    for var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "RAYON_NUM_THREADS",
        "POLARS_MAX_THREADS",
        "NUMBA_NUM_THREADS",
        "TBB_NUM_THREADS",
    ):
        os.environ[var] = str(n)
    if threadpool_limits is not None:
        threadpool_limits(limits=n)


_clamp_native_threads(1)


_WORKER_DECOMPOSER = None
_WORKER_FACTOR = 1.0
_WORKER_HDU = None
_WORKER_FLUX = {}
_WORKER_LSF = {}
_WORKER_META = None
_WORKER_PROGRESS_QUEUE = None
_WORKER_FIT_MODEL = "baseline"
_WORKER_EXPOSURE_SECONDS = 900.0
_WORKER_SCIENCE_LINE_MASK = None
_WORKER_SCIENCE_LINE_FWHM = None
_WORKER_SCIENCE_LINE_CENTRE = True  # overwritten by init_worker

FIT_MODEL_SUFFIXES = {
    "baseline": "",
    "lsf-surface-iterative": "_lsf_surface_iterative",
    "lsf-surface-iterative-split-zodi": "_lsf_surface_iterative_split_zodi",
    "moon-zodi-lsf-surface-iterative": "_moon_zodi_lsf_surface_iterative",
}
MOON_ZODI_FIT_MODEL = "moon-zodi-lsf-surface-iterative"
SPLIT_ZODI_FIT_MODEL = "lsf-surface-iterative-split-zodi"

# Defaults for the SkyDecompLSFSurfaceIterative(split_zodi=True) knobs; match the
# settings validated on the p40_p70 every10 identifiability notebook.
SPLIT_ZODI_N_KNOTS_DEFAULT = 1
SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT = 1.0e-1
SPLIT_ZODI_MOON_ALBEDO_PHASE_DEG = 30.0
SPLIT_ZODI_COLOR_EXPONENT = 0.26
# --- split-zodi identifiability defaults (validated 2026-09-03) --------------
# Without these the split is degenerate in a way fit quality cannot see: on 200
# lunation-stratified sky spectra the deployed configuration put the moon and
# zodi colours in the WRONG ORDER on 164 of 168 moon-up spectra (fitted moon
# log-log slope +0.21 where physics says -3.7, zodi -4.04 where physics says
# -0.3), left the moon family holding 45% of the continuum with the moon 37 deg
# BELOW the horizon, and had the fitted zodi tracking lunar illumination at
# rho = +0.95 while retaining only rho = 0.09 of its Leinert B500 dependence.
#
# The four knobs below fix those, and they are not interchangeable:
#   * the RATIO BOUNDS fix the colour ordering (reversals 164 -> 2 of 168) and
#     stop a spline zeroing out mid-band, which is how the old fit reached a
#     lower rms -- it used the moon family as piecewise scratch space;
#   * the FRACTION bracket fixes dark time (moon share 0.45 -> 0.02);
#   * the absolute LEINERT ANCHOR fixes the amplitude geometry
#     (rho(zodi, B500) 0.09 -> 0.90, rho(zodi, FLI) 0.94 -> 0.12).
# Dropping the anchor keeps the reversals fixed but leaves the zodi amplitude
# tracking the moon; dropping the bounds re-opens the reversals.
#
# Cost: median rms x1.01 over the sample, concentrated entirely at bright moon
# (x1.49 median for FLI > 0.8, ~1% of the continuum, blue-weighted).  Part of
# that is the baseline overfitting via the spline hole described above.
SPLIT_ZODI_MOON_RATIO_BOUND = 0.7
SPLIT_ZODI_ZODI_RATIO_BOUND = 0.7
SPLIT_ZODI_AMP_PRIOR_TOL = 3.0
SPLIT_ZODI_ZODI_AMP_BOUND = 2.0

# Diffuse species-ratio bracket, ON by default since 2026-09-10.  The three
# diffuse species are individually unidentifiable in the LVM band: on the
# canonical-PALACE basis the free fit spreads log10(FeO/HO2) over 10.2 dex and
# zeroes HO2 on 23% of rows, and the THREE ARMS OF ONE EXPOSURE -- same sky, a
# few degrees apart -- disagree by 0.633 dex at the median.  Airglow does not
# vary by x4 in a species ratio over 5 deg, so that spread is fitting noise,
# and in dark time the amplitude it moves leaks straight into the zodi.
#
# Validated on 500 rows stratified over moon state x |ecl_beta|, two full
# runs: FeO/HO2 5-95 span 10.25 -> 0.40 dex, HO2 near-zero rows 23.2% -> 0.2%,
# for +2.46% of sci blue chi2 (+1.60% near, +0.13% far) and no measurable
# full-band cost.  It moves the SPLIT, not the sum: per-row log10(prior/free)
# on sci is zodi p90 0.029, diffuse p90 0.117, zodi+diffuse SUM p90 0.012.
# Bright-moon rows do not move at all -- their zodi is pinned at the Leinert
# ceiling -- so the bracket acts only where the diffuse dominates.
#
# NOMINAL: FLUX shares (HO2, FeO, O2Ac), the geometric median of the fitted
# corpus, NOT PALACE's own reference shares -- centring on PALACE costs 25.7%
# of the blue chi2 against 0.67% for the corpus median, because PALACE is
# calibrated for Cerro Paranal and LVM observes from LCO.  RE-MEASURE THIS on
# the corpus being fitted if the basis or the corpus changes; the value below
# was measured on 60 every10 rows of gaia-stars-mask on the _canonhyb_v1
# basis.
#
# KNOWN FAILURE MODE: with c >= 0 a ratio bound makes the block all-positive
# or all-zero, so a row whose fit wants HO2 = 0 loses its entire diffuse
# block (measured 0.2-0.4% of rows).  `diffuse_zeroed_keep_mask` catches
# those downstream; count them after a run.
SPLIT_ZODI_DIFFUSE_RATIO_BOUND_DEX = 0.2
SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL = (0.0396, 0.7026, 0.2578)

# Moon-gated upper bound on the DIFFUSE BLOCK relative to OH.
# The three diffuse species are mesospheric chemiluminescence and cannot
# depend on the moon, yet on gaia-stars-mask-cont their amplitude relative to
# OH rises with moon_frac_po -- rho +0.716 for the block (FeO alone +0.693,
# and HO2/O2Ac inherit it through the species-ratio bracket) -- while airmass
# and van Rhijn give +0.015.  The templates are absorbing scattered
# moonlight: ~74% of the fitted FeO on full-moon rows, 3-5% of the fitted
# continuum.  OH itself is moon-independent (rho +0.088), so it is a clean
# normaliser.
#
# SCOPE IS THE BLOCK, not FeO alone.  An FeO-only cap was implemented and
# validated first; it worked (chi2-free, 99% of the released flux to the
# moon) but on ~40% of gated rows the +/-0.2 dex species-ratio bracket became
# the binding constraint instead -- log10(FeO/HO2) lower-edge occupancy rose
# 25.8% -> 40.5% -- so FeO could not fall further and the Noll anchor barely
# moved.  Capping the block lets the ratio bracket distribute the reduction
# rather than block it.  The block ratio is also better behaved: dark-time
# robust sigma 0.216 dex against 0.306 for FeO, and Theil-Sen slope against
# OH +1.16 against +1.40.
#
# CENTRE is the dark-time median of log10(A_diffuse/A_OH) over diffuse-live
# rows, where there is no moon to leak.  Re-measure it per corpus.
#
# ONE-SIDED AND GATED, both forced by measurement: the dark-time scatter is
# real -- clipping it cost 15-32% of the blue chi2 in the FeO-only sizing
# test -- so a two-sided or ungated bound is unaffordable.
#
# RELAX = 0, i.e. a flat width above the gate.  The ramp was there to avoid a
# discontinuity at the gate, but it is unnecessary: the measured excess is
# only -0.04 to +0.04 dex just above the gate, well inside the allowed 0.25,
# so the bound is naturally inactive there and turns on where the excess
# exceeds the legitimate scatter.  Median excess by moon_frac_po:
# +0.04 dex at 0.61-0.82, +0.22 at 0.82-0.91, +0.37 at 0.91-0.95,
# +0.52 at 0.95-1.00.
#
# BOUND TIGHTENED 0.25 -> 0.15 dex on 2026-09-11.  At 0.25 the bound sat at a
# ratio of 0.3991, ABOVE the observed 0.244 and 0.374 in the two lower gated
# bins, so the real moon dependence over moon_frac_po 0.6-0.9 went untouched
# and the residual correlation stayed at +0.644.  0.15 puts the bound at
# 0.317, which bites in those bins too.  Affordable because the cap measured
# chi2-FREE at 0.25: blue -0.21% on binding rows, 0.00% elsewhere, zero dark
# rows touched, collapse rate unchanged at 2.40%.
#
# KNOWN RISK: 0.15 dex is ~0.7x the dark-time robust sigma of
# log10(A_diffuse/A_OH) (0.216 dex), which is our only estimate of the
# intrinsic spread.  On gated rows the intrinsic and moon-driven parts cannot
# be separated, so this may clip legitimate variation.  Watch the blue chi2 on
# gated-but-not-binding rows and the diffuse-collapse rate.
SPLIT_ZODI_DIFFUSE_OH_CENTRE_LOG10 = -0.6489
SPLIT_ZODI_DIFFUSE_OH_BOUND_DEX = 0.15
SPLIT_ZODI_DIFFUSE_OH_GATE_FRAC = 0.6
SPLIT_ZODI_DIFFUSE_OH_RELAX_DEX = 0.0
# Absolute recentring of the Leinert anchor.  The anchor brackets the fitted
# zodi total to [Z_pred/kappa_z, kappa_z * Z_pred], and Z_pred comes from
# _physics_only_model, whose learned scale factors are deliberately zeroed --
# so nothing has ever calibrated its ABSOLUTE normalisation.  Measured on the
# new-oh-2 ML test split (1113 rows) by comparing each fitted zodi total with
# that prediction, the anchor turned out to be SATURATED: 67.7% of all rows and
# 93.1% of moon-up rows sat exactly on the ceiling, every percentile p10-p99 of
# log10(Z_fit/Z_pred) equal to +0.3010 = log10(2.0) to four decimals.  The fit
# was not measuring zodi there, it was reporting kappa_z.
#
# The gap is multiplicative, not a pedestal: log10(Z_fit/Z_pred) has slope
# +0.08 (rho = +0.07) against log10(Z_pred), and a two-parameter c*Z + p fit
# does no better than pure c (23.2% vs 23.5% median error) with an unphysical
# NEGATIVE p.  So Zodi_bs is not absorbing a non-zodiacal continuum; the
# normalisation is simply low.
#
# The size must be measured where there is no moonlight to leak, and with the
# censoring undone -- the interior rows are interior BECAUSE |r| < log10(2), so
# their median is biased toward 1.  A censored-Gaussian MLE gives:
#     moon down  1.61x (sigma 0.274 dex)   <- leakage-free, this is the number
#     moon up    8.32x                     <- not calibration: a calibration
#                                             offset cannot depend on the moon
#                                             (per-quartile: 3.1/9.4/3.8/3.8x)
# The moon carrier is independently 1.35x low (uncensored, IQR 0.120 dex), so
# one shared physical_to_fit_flux_scale error of ~1.4x explains both families.
# That is why the correction is applied to the zodi TOTAL only: the moon
# FRACTION is a ratio through the same conversion, so a shared factor cancels
# and SPLIT_ZODI_AMP_PRIOR_TOL keeps policing moon-into-zodi leakage unchanged.
#
# Checked on the same 200 lunation-stratified sky spectra the bounds were
# adopted on.  Every guardrail holds -- reversals 2/167 (was 2/164), dark-time
# moon share 0.0234 (unchanged to four digits), median rms 0.987x -- and the
# improvement is concentrated where the anchor was worst: rows pinned to a
# bound fall 77% -> 64% overall and 86% -> 71% with the moon up.
#
# Choosing the WIDTH matters as much as the centre, because the bracket is
# [c/kappa_z, c*kappa_z] * Z_raw and Z_true ~ 1.6 * Z_raw, so kappa_z sets the
# bracket in physical units.  Swept at c = 1.6:
#   kappa_z  physical bracket   pinned all/up/dark   rms    rms(FLI>0.8)
#     1.25   [0.80, 1.25]         88% / 90% / 75%   1.0000     1.0000
#     2.00   [0.50, 2.00]         64% / 71% / 21%   0.9867     0.9805
#     3.00   [0.33, 3.00]         44% / 51% /  0%   0.9783     0.9558
# kappa_z = 2.0 is kept: it is ~1.1 sigma of the measured 0.274 dex dark-time
# spread, so it BOUNDS the zodi without dictating it.  1.25 reproduces today's
# fits almost exactly (rms 1.0000 in every lunation bin) because everything
# still sits on a barely-moved ceiling, and it makes 75% of dark-time targets
# synthetic; 3.0 frees dark time completely but gives bright-moon zodi 3x
# headroom against a QP that already demands 8.3x.
#
# Two things NOT to conclude from the surrounding diagnostics.  (a) Tightening
# SPLIT_ZODI_AMP_PRIOR_TOL does not substitute for this: at kappa_f 1.5 and 1.2
# reversals rose to 6 and 10 of 167, and at c = 1.0 tightening it changed
# nothing at all (zodi_tot x1.000) because the anchor already pins the rows it
# would act on.  (b) rho(zodi, B500) and rho(zodi, FLI) are NOT trustworthy
# while the anchor binds: on a bound, zodi_tot == kappa_z * c * Z_pred exactly,
# so those correlations partly measure the constraint.  On rows interior in
# both c = 1.0 and c = 1.6 the fitted zodi is identical to machine precision
# (1.000x, IQR 0.0000) and every config gives the same rho(B500) ~ 0.74,
# rho(FLI) ~ 0.08 and partial rho(FLI | B500) ~ -0.12.  Judge this constraint
# by pinning fraction, reversals and rms, not by those correlations.
SPLIT_ZODI_ZODI_PRIOR_CALIBRATION = 1.6
# Moon_bs interior-knot count.  Deliberately NOT SkyDecomp.__init__'s default
# (25, with n_zodi_spline_knots 3): the deployed corpus and every measurement
# behind the SPLIT_ZODI_* bounds above use 11 moon / 1 zodi interior knots.
# The ratio bounds are per ADJACENT KNOT PAIR, so the same beta is looser the
# more knots there are -- changing these without re-validating the bounds
# changes how much colour freedom each family actually has.
MOON_N_KNOTS_DEFAULT = 11

# --- Science emission-line mask -------------------------------------------
# Nebular lines from the SCIENCE field are not sky, and NONE of them exists in
# the 388-component basis, so the QP is forced to absorb them into whatever it
# has.  Measured on every10 row 773 (expnum 39622, galactic b = +0.48, an
# inner-plane H II region): Halpha equivalent width 62.4 A in the science fibre
# against 10.0 and 2.8 in the two sky arms, and the FITTED OH component runs
# 6.7x its sideband level inside the Halpha/[NII] window and 2.0x inside
# [SII].  So the absorber is OH -- the only narrow-line family with lines
# there -- not the moon or zodi splines, which is why the science-continuum
# colour gate in mlp_predictor.data does not see these rows and why nebular
# equivalent width shows no corpus-wide correlation with moon or zodi
# distortion (rho -0.107 and +0.044).
#
# Zeroing IVAR is the whole mechanism: SkyDecompBase._fit_design selects on
# `np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)`, so masked pixels leave
# the QP, chi2, the dof count and R2 together, with no interpolation and no
# change to the wavelength grid.
#
# COST, measured against the EXACT deployed basis -- design matrix rebuilt
# from the stored per-row LSF surface and validated to 7e-16 against the
# stored COMP_OH.  A design matrix built from the INITIAL LSF instead is wrong
# by a factor of 12 on OH and must not be used for this; it also mislocates OH
# line peaks, which made [OIII]4959 look coincident with OH_314 when it is not.
#
# At the deployed widths (1.5 x FWHM, so +/-2.04 A at Halpha) the mask covers
# 85/12401 pixels = 0.69% and removes 0.96% of the OH model, 0.95% of the
# moon, 0.70% of the zodi, 0.76% of the diffuse block and 0.00% of ATOM.
# Three OH components lose more than 20% of their own support -- OH_316 64%,
# OH_144 45%, OH_341 32% -- carrying 0.238% of the OH flux between them
# (ranks 81-217 of 357), and each keeps support elsewhere: OH_341 at
# 4034-5475 A, OH_316 at 3964-8175, and OH_144 its 6555.5 A line, which the
# window leaves outside while masking its 6561.5 A one.
#
# Almost all of that cost is Halpha and [NII]6583.  [OII]3726/3729,
# [OIII]4959 and [OIII]5007 cost NOTHING -- no component loses even 10%.
# WIDTH IS THE SENSITIVE KNOB, and it is asymmetric: widening buys very little
# line flux (the core is already covered) and costs OH support quickly --
# 2.0 x FWHM takes OH from 0.96% to 1.09% and pulls OH_285 past 20%, and a
# +/-8 A window would take OH_144 to 94% (its whole support is two lines, at
# 6555.5 and 6561.5 A) and add OH_150.  Re-measure before widening.
#
# Both [NII] and both [SII] and [OIII]4959 are included even though only the
# brighter partner was asked for: they are fixed-ratio partners of lines
# already masked (1/3 of [NII]6583 and of [OIII]5007), they sit inside or
# beside the same windows, and leaving one of a doublet unmasked leaves the
# contaminant in the fit at a third of its strength.
SCIENCE_EMISSION_LINES = (
    ("[OII]3726", 3726.03),
    ("[OII]3729", 3728.82),
    ("Hbeta", 4861.33),
    ("[OIII]4959", 4958.91),
    ("[OIII]5007", 5006.84),
    ("[NII]6548", 6548.05),
    ("Halpha", 6562.80),
    ("[NII]6583", 6583.45),
    ("[SII]6716", 6716.44),
    ("[SII]6731", 6730.82),
)
# Half-width = max(MIN, FWHM_MULTIPLE * FWHM(lambda)) + lambda * v/c, scaled
# off the detector FWHM at each line rather than one assumed resolution.
#
# 1.5 x FWHM (3.53 sigma) masks the bright core, which is where essentially all
# the contaminating flux is, and deliberately leaves the wings in: doubling the
# window to 2.0 x FWHM buys almost no extra line flux while taking OH from
# 0.96% to 1.09% and pulling a fourth component past 20% support loss.
#
# A displaced line is not a wing effect -- the core itself moves -- so instead
# of widening the windows the mask SLIDES them, using a Halpha velocity
# measured per row (SCIENCE_LINE_MASK_CENTRE_ON_HALPHA).  The extra widening
# term SCIENCE_LINE_MASK_VELOCITY_KM_S is therefore 0 by default.
#
# THE CENTROID MUST BE MEASURED AGAINST A SKY REFERENCE.  There are OH lines
# directly under Halpha -- OH_144 at 6555.5/6561.5 A, OH_157 at 6559.5/6562.5,
# OH_285 at 6557.5/6561.5, OH_316 at 6563.0 -- and they are bright enough to
# drag the centroid: the two SKY fibres, which are essentially pure OH in this
# window, centroid at -3.08 and -4.20 A.  Measured on 2697 science fibres with
# Halpha EW excess > 10 A:
#
#   raw science flux      median -0.71 A (-32 km/s), p1 -4.33, p99 +4.12,
#                         87.9% within +/-2.04 A
#   minus a scaled sky arm  median +0.17 A ( +8 km/s), p1 -2.07, p99 +5.83,
#                         94.4% within +/-2.04 A
#
# So the true nebular velocities sit much closer to rest than the raw numbers
# suggest, and most of the apparent -32 km/s was OH.
#
# WHAT CENTRING ACTUALLY BUYS, validated on those 2697 rows: the pipeline
# measures a velocity on 92.5% of them (median +8.7 km/s over the whole
# every10 sample, matching the +8 km/s above; the other 73% of all rows have
# no detectable nebular Halpha and correctly fall back to rest).  Core
# coverage goes 89.0% -> 92.5%: of the 297 rows the rest-frame mask missed,
# centring covers 102, while 8 rows (0.30%) that WERE covered are now missed
# through a mis-measured velocity.  It is a tail fix worth roughly +3.5pp, not
# a bulk one -- do not expect it to change aggregate metrics.
#
# Masked pixel count stays 83-87 against the rest-frame 85, so dof moves by at
# most a couple of pixels and reduced_chi2 stays comparable across rows.
SCIENCE_LINE_MASK_ENABLED = True
SCIENCE_LINE_MASK_FWHM_MULTIPLE = 1.5
SCIENCE_LINE_MASK_VELOCITY_KM_S = 0.0
SCIENCE_LINE_MASK_MIN_HALF_WIDTH_A = 2.0
# Per-row centring: measure the nebular velocity from Halpha in the SCIENCE
# fibre and slide every window to the observed wavelength instead of widening
# it.  One velocity per row, applied to all ten lines and all three arms --
# the shift is a property of the emitting gas, not of the line, and keeping it
# common across arms means the three fits still exclude exactly the same
# pixels, which is what makes their coefficients comparable.
SCIENCE_LINE_MASK_CENTRE_ON_HALPHA = True
SCIENCE_LINE_MASK_MAX_SHIFT_KM_S = 300.0
SCIENCE_LINE_MASK_CENTRE_MIN_SNR = 5.0
_C_KM_S = 299792.458


def measure_halpha_velocity(
    wave,
    flux_sci,
    flux_sky=None,
    max_shift_km_s=SCIENCE_LINE_MASK_MAX_SHIFT_KM_S,
    min_snr=SCIENCE_LINE_MASK_CENTRE_MIN_SNR,
    search_half_width_a=12.0,
):
    """Nebular Halpha velocity of one row, in km/s, or 0.0 if not measurable.

    ``flux_sky`` is a sky-fibre spectrum used to cancel the OH lines that sit
    directly under Halpha; pass ``None`` only if none is available, and expect
    a blueward bias of a few km/s to tens of km/s if you do (the sky fibres
    centroid at -3.1 and -4.2 A in this window because of OH alone).  The sky
    is scaled by the ratio of positive flux in the window, capped at 1, so it
    can only remove the shared airglow, never add a negative pedestal deeper
    than the science spectrum itself.

    Returns 0.0 -- i.e. leave the windows at rest -- whenever the line is not
    convincingly detected, the centroid lands outside the search window, or
    the implied shift exceeds ``max_shift_km_s``.  Failing closed matters:
    a spurious shift moves the mask off a line that WAS being masked.
    """
    wave = np.asarray(wave, dtype=np.float64)
    lam0 = 6562.80
    lo, hi = lam0 - float(search_half_width_a), lam0 + float(search_half_width_a)
    i0, i1 = np.searchsorted(wave, [lo, hi])
    if int(i1) - int(i0) < 5:
        return 0.0
    sl = slice(int(i0), int(i1))
    w = wave[sl]
    y = np.asarray(flux_sci, dtype=np.float64)[sl]
    if flux_sky is not None:
        sky = np.asarray(flux_sky, dtype=np.float64)[sl]
        num = float(np.nansum(np.clip(y, 0.0, None)))
        den = float(np.nansum(np.clip(sky, 0.0, None)))
        if den > 0.0:
            y = y - min(num / den, 1.0) * sky
    # Local continuum from sidebands either side of the Halpha/[NII] complex.
    cw = (np.abs(wave - lam0) > 25.0) & (np.abs(wave - lam0) <= 70.0)
    base = np.asarray(flux_sci, dtype=np.float64)[cw]
    cont = float(np.nanmedian(base)) if np.any(np.isfinite(base)) else 0.0
    net = y - cont
    good = np.isfinite(net)
    if not np.any(good):
        return 0.0
    pos = np.where(good & (net > 0.0), net, 0.0)
    total = float(pos.sum())
    if total <= 0.0:
        return 0.0
    # Detection test against the scatter of the negative excursions, which is
    # what the window looks like when there is no line.
    neg = net[good & (net < 0.0)]
    noise = float(np.std(neg)) if neg.size >= 3 else 0.0
    peak = float(np.nanmax(pos))
    if noise > 0.0 and peak < float(min_snr) * noise:
        return 0.0
    centroid = float((pos * w).sum() / total)
    if not np.isfinite(centroid) or centroid <= lo or centroid >= hi:
        return 0.0
    v = (centroid - lam0) / lam0 * _C_KM_S
    if not np.isfinite(v) or abs(v) > float(max_shift_km_s):
        return 0.0
    return v


def science_line_mask(
    wave,
    lsf_fwhm=None,
    lines=SCIENCE_EMISSION_LINES,
    fwhm_multiple=SCIENCE_LINE_MASK_FWHM_MULTIPLE,
    velocity_km_s=SCIENCE_LINE_MASK_VELOCITY_KM_S,
    min_half_width_a=SCIENCE_LINE_MASK_MIN_HALF_WIDTH_A,
    centre_velocity_km_s=0.0,
):
    """Boolean mask over ``wave``, True where a science emission line sits.

    Each window is centred on ``lambda * (1 + centre_velocity_km_s/c)`` --
    normally the per-row Halpha velocity from ``measure_halpha_velocity``,
    applied to every line because the shift belongs to the emitting gas -- and
    has half-width
    ``max(min_half_width_a, fwhm_multiple * FWHM(lambda)) + lambda * v/c``.

    ``lsf_fwhm`` is the DETECTOR LSF FWHM in Angstrom -- the LSF_* arrays in
    the input FITS, which is what the pipeline actually carries (median 1.57 A
    on gaia1over100, 1.36 A at Halpha).  Do NOT pass ``--lsf-sigma`` here: that
    is a scalar Gaussian SIGMA defaulting to 0.5 A, a different quantity by a
    factor of 2.35.  A per-pixel array, a scalar, or ``None`` are all accepted;
    ``None`` (or an array with no usable pixel near a line) falls back to
    ``min_half_width_a``.

    Returns ``(mask, widths)``, ``widths`` being the per-line half-width in
    Angstrom so the caller can report what it actually masked.
    """
    wave = np.asarray(wave, dtype=np.float64)
    fwhm = None
    if lsf_fwhm is not None:
        arr = np.asarray(lsf_fwhm, dtype=np.float64)
        if arr.ndim == 0:
            arr = np.full(wave.shape, float(arr))
        if arr.shape == wave.shape:
            fwhm = arr
    mask = np.zeros(wave.shape, dtype=bool)
    widths = {}
    shift = 1.0 + float(centre_velocity_km_s) / _C_KM_S
    for name, lam in lines:
        centre = float(lam) * shift
        half = float(min_half_width_a)
        if fwhm is not None:
            near = np.abs(wave - centre) <= 25.0
            usable = near & np.isfinite(fwhm) & (fwhm > 0)
            if np.any(usable):
                half = max(half,
                           float(fwhm_multiple) * float(np.median(fwhm[usable])))
        half += float(lam) * float(velocity_km_s) / _C_KM_S
        widths[name] = half
        mask |= (wave >= centre - half) & (wave <= centre + half)
    return mask, widths



def init_worker(
    wave,
    lsf_sigma,
    base_dir,
    factor,
    data_file,
    progress_queue=None,
    fit_model="baseline",
    n_refinement_cycles=5,
    worker_counter=None,
    pin_cpu=False,
    diagnose_threads=False,
    palace_suffix=None,
    palace_oh_suffix=None,
    palace_diffuse_suffix=None,
    exposure_seconds=900.0,
    moon_zodi_data_root=None,
    n_spline_knots=MOON_N_KNOTS_DEFAULT,
    n_zodi_spline_knots=SPLIT_ZODI_N_KNOTS_DEFAULT,
    zodi_smooth_lambda=SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT,
    mask_science_lines=SCIENCE_LINE_MASK_ENABLED,
    centre_on_halpha=SCIENCE_LINE_MASK_CENTRE_ON_HALPHA,
    diffuse_ratio_bound_dex=SPLIT_ZODI_DIFFUSE_RATIO_BOUND_DEX,
    diffuse_ratio_nominal=SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL,
    diffuse_oh_centre_log10=SPLIT_ZODI_DIFFUSE_OH_CENTRE_LOG10,
    diffuse_oh_bound_dex=SPLIT_ZODI_DIFFUSE_OH_BOUND_DEX,
):
    """Initialise one SkyDecomp instance per worker process."""
    global \
        _WORKER_DECOMPOSER, \
        _WORKER_FACTOR, \
        _WORKER_HDU, \
        _WORKER_FLUX, \
        _WORKER_LSF, \
        _WORKER_META, \
        _WORKER_PROGRESS_QUEUE, \
        _WORKER_FIT_MODEL, \
        _WORKER_EXPOSURE_SECONDS, \
        _WORKER_SCIENCE_LINE_MASK, \
        _WORKER_SCIENCE_LINE_FWHM, \
        _WORKER_SCIENCE_LINE_CENTRE

    _clamp_native_threads(1)

    worker_rank = 0
    if worker_counter is not None:
        with worker_counter.get_lock():
            worker_rank = int(worker_counter.value)
            worker_counter.value = worker_rank + 1
    if pin_cpu and hasattr(os, "sched_setaffinity"):
        try:
            available = sorted(os.sched_getaffinity(0))
            if available:
                target = available[worker_rank % len(available)]
                os.sched_setaffinity(0, {target})
        except OSError as exc:
            print(f"[worker pid={os.getpid()}] pin_cpu failed: {exc}", flush=True)

    _WORKER_FACTOR = float(factor)
    _WORKER_FIT_MODEL = fit_model
    _WORKER_EXPOSURE_SECONDS = float(exposure_seconds)
    # Keep worker-local memmapped access to flux tables to avoid large IPC payloads.
    _WORKER_HDU = fits.open(data_file, memmap=True)
    _WORKER_PROGRESS_QUEUE = progress_queue
    _WORKER_FLUX = {
        "sci": np.asarray(_WORKER_HDU["FLUX_SCI"].data),
        "sky1": np.asarray(_WORKER_HDU["FLUX_SKY_NEAR"].data),
        "sky2": np.asarray(_WORKER_HDU["FLUX_SKY_FAR"].data),
    }
    # split-zodi needs the same LSF + META as the model-based mode: its
    # amplitude priors are geometry predictions, one per spectrum.
    if fit_model in (MOON_ZODI_FIT_MODEL, SPLIT_ZODI_FIT_MODEL):
        _WORKER_LSF = {
            "sci": np.asarray(_WORKER_HDU["LSF_SCI"].data),
            "sky1": np.asarray(_WORKER_HDU["LSF_SKY_NEAR"].data),
            "sky2": np.asarray(_WORKER_HDU["LSF_SKY_FAR"].data),
        }
        _WORKER_META = _WORKER_HDU["META"].data
    else:
        _WORKER_LSF = {}
        _WORKER_META = None
    if fit_model == "baseline":
        from skysub.sky_decomp.fit import SkyDecomp

        _WORKER_DECOMPOSER = SkyDecomp(
            wave,
            lsf_sigma=lsf_sigma,
            base_dir=base_dir,
            palace_suffix=palace_suffix,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            moon_smooth_lambda=0.1,
            moon_interline_boost=10000.0,
            moon_interline_red_min=6000.0,
            moon_interline_exclusion_a=2.5,
            moon_interline_line_flux_threshold=0.01,
        )
    elif fit_model == "lsf-surface-iterative":
        from skysub.sky_decomp.lsf_surface_iterative import (
            LSFSurfaceIterativeConfig,
            SkyDecompLSFSurfaceIterative,
        )

        _WORKER_DECOMPOSER = SkyDecompLSFSurfaceIterative(
            wave,
            lsf_sigma=lsf_sigma,
            base_dir=base_dir,
            palace_suffix=palace_suffix,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            moon_smooth_lambda=0.1,
            moon_interline_boost=0.0,
            n_spline_knots=int(n_spline_knots),
            config=LSFSurfaceIterativeConfig(
                n_refinement_cycles=n_refinement_cycles,
            ),
        )
    elif fit_model == SPLIT_ZODI_FIT_MODEL:
        from skysub.sky_decomp.lsf_surface_iterative import (
            LSFSurfaceIterativeConfig,
            SkyDecompLSFSurfaceIterative,
        )

        _WORKER_DECOMPOSER = SkyDecompLSFSurfaceIterative(
            wave,
            lsf_sigma=lsf_sigma,
            base_dir=base_dir,
            palace_suffix=palace_suffix,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            moon_smooth_lambda=0.1,
            moon_interline_boost=0.0,
            n_spline_knots=int(n_spline_knots),
            split_zodi=True,
            n_zodi_spline_knots=int(n_zodi_spline_knots),
            zodi_smooth_lambda=float(zodi_smooth_lambda),
            moon_albedo_fiducial_phase_deg=SPLIT_ZODI_MOON_ALBEDO_PHASE_DEG,
            zodi_color_exponent=SPLIT_ZODI_COLOR_EXPONENT,
            moon_ratio_bound=SPLIT_ZODI_MOON_RATIO_BOUND,
            zodi_ratio_bound=SPLIT_ZODI_ZODI_RATIO_BOUND,
            amp_prior_tol=SPLIT_ZODI_AMP_PRIOR_TOL,
            zodi_amp_bound=SPLIT_ZODI_ZODI_AMP_BOUND,
            diffuse_ratio_bound_dex=float(diffuse_ratio_bound_dex),
            diffuse_ratio_nominal=diffuse_ratio_nominal,
            diffuse_oh_centre_log10=diffuse_oh_centre_log10,
            diffuse_oh_bound_dex=float(diffuse_oh_bound_dex),
            diffuse_oh_gate_frac=SPLIT_ZODI_DIFFUSE_OH_GATE_FRAC,
            diffuse_oh_relax_dex=SPLIT_ZODI_DIFFUSE_OH_RELAX_DEX,
            diffuse_oh_scope="block",
            config=LSFSurfaceIterativeConfig(
                n_refinement_cycles=n_refinement_cycles,
            ),
        )
    elif fit_model == MOON_ZODI_FIT_MODEL:
        from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
        from skysub.sky_decomp.moon_zodi_lsf_surface_iterative import (
            SkyDecompMoonZodiLSFSurfaceIterative,
        )
        _WORKER_DECOMPOSER = SkyDecompMoonZodiLSFSurfaceIterative(
            wave,
            lsf_sigma=lsf_sigma,
            data_root=(
                DEFAULT_MOON_ZODI_DATA_ROOT
                if moon_zodi_data_root is None
                else moon_zodi_data_root
            ),
            palace_suffix=palace_suffix,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            moon_smooth_lambda=0.1,
            moon_interline_boost=0.0,
            physical_to_fit_flux_scale=float(factor),
            config=LSFSurfaceIterativeConfig(
                n_refinement_cycles=n_refinement_cycles,
            ),
        )
    else:
        raise ValueError(f"Unknown fit model: {fit_model}")

    # Science emission-line mask, built once per worker.  Deliberately built
    # from the init-time LSF rather than per row: the row-to-row LSF variation
    # is small next to the +/-150 km/s velocity term, and a row-dependent mask
    # would make the number of fitted pixels vary from row to row, which
    # reduced_chi2 and the dof count would then carry.
    _WORKER_SCIENCE_LINE_MASK = None
    _WORKER_SCIENCE_LINE_FWHM = None
    _WORKER_SCIENCE_LINE_CENTRE = bool(centre_on_halpha)
    if mask_science_lines:
        # Use the DETECTOR LSF FWHM from the input FITS, median-combined over
        # rows and arms, not the scalar `lsf_sigma` argument -- that is a
        # 0.5 A Gaussian sigma by default while the real FWHM is ~1.57 A, so
        # feeding it here would size every window off the wrong quantity.
        # _WORKER_LSF is empty for fit models that do not need it; then the
        # instrumental term is simply dropped.
        _fwhm_ref = None
        if _WORKER_LSF:
            _stack = [np.asarray(v, dtype=np.float64) for v in _WORKER_LSF.values()]
            _stack = [a if a.ndim == 1 else np.nanmedian(a, axis=0) for a in _stack]
            _stack = [a for a in _stack if a.shape == np.shape(wave)]
            if _stack:
                _fwhm_ref = np.nanmedian(np.vstack(_stack), axis=0)
        # Use the `wave` ARGUMENT, not _WORKER_DECOMPOSER.wave: every
        # decomposer is constructed from this same array, and reading it off
        # the object couples worker init to the concrete decomposer class
        # (test_palace_suffix's FakeDecomposer has no `.wave`).
        _mask, _widths = science_line_mask(wave, _fwhm_ref)
        _WORKER_SCIENCE_LINE_MASK = _mask
        # Kept so fit_chunk_worker can rebuild the mask per row at the
        # measured Halpha velocity; the static mask above stays as the
        # fallback for rows where the line is not measurable.
        _WORKER_SCIENCE_LINE_FWHM = _fwhm_ref
        if worker_rank == 0:
            print(f"[science-line mask] {int(_mask.sum())}/{_mask.size} pixels "
                  f"({100.0 * _mask.mean():.2f}%) excluded via IVAR=0 in "
                  f"{len(SCIENCE_EMISSION_LINES)} windows: "
                  + ", ".join(f"{n} +/-{w:.1f}A" for n, w in _widths.items())
                  + ("; centred per row on the measured Halpha velocity"
                     if _WORKER_SCIENCE_LINE_CENTRE else "; fixed at rest"),
                  flush=True)

    # After all heavy imports, clamp once more and (optionally) report per-worker state.
    _clamp_native_threads(1)
    if diagnose_threads and worker_rank == 0:
        _report_thread_diagnostics()


def _report_thread_diagnostics():
    lines = [f"[worker pid={os.getpid()}] thread diagnostics:"]
    lines.append(
        f"  affinity_cores={sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else 'n/a'}"
    )
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "RAYON_NUM_THREADS",
        "NUMBA_NUM_THREADS",
        "TBB_NUM_THREADS",
    ):
        lines.append(f"  {var}={os.environ.get(var, 'unset')}")
    if threadpool_limits is not None:
        try:
            from threadpoolctl import threadpool_info

            for entry in threadpool_info():
                lines.append(
                    f"  loaded_pool: {entry.get('user_api', '?'):8s} "
                    f"{entry.get('prefix', '?'):18s} threads={entry.get('num_threads', '?')}"
                )
        except Exception as exc:
            lines.append(f"  threadpool_info failed: {exc}")
    else:
        lines.append("  (threadpoolctl not installed; cannot enumerate loaded pools)")
    print("\n".join(lines), flush=True)


def _text_value(value):
    return value.decode().strip() if isinstance(value, bytes) else str(value).strip()


def _moon_zodi_observation(kind, row_index):
    from skysub.sky_decomp.moon_zodi_model import MoonZodiObservation

    role_contract = {
        "sci": ("sci", "sci_ra", "sci_dec"),
        "sky1": ("sky_near", "sky_near_ra", "sky_near_dec"),
        "sky2": ("sky_far", "sky_far_ra", "sky_far_dec"),
    }
    role, ra_column, dec_column = role_contract[kind]
    row = _WORKER_META[row_index]
    names = set(_WORKER_META.dtype.names or ())
    exposure = None
    for column in ("exposure_seconds", "exptime"):
        if column in names:
            candidate = float(row[column])
            if np.isfinite(candidate) and candidate > 0.0:
                exposure = candidate
                break
    if exposure is None:
        exposure = _WORKER_EXPOSURE_SECONDS
        exposure_source = "assumed_900s"
    else:
        exposure_source = "metadata"
    return MoonZodiObservation(
        expnum=int(row["expnum"]),
        date_obs=_text_value(row["date_obs"]),
        role=role,
        target_ra_deg=float(row[ra_column]),
        target_dec_deg=float(row[dec_column]),
        exposure_seconds=float(exposure),
        exposure_seconds_source=exposure_source,
    )


_LSF_REPAIR_COUNT = {"rows": 0, "pixels": 0, "reported": False}


def _sanitised_lsf_row(kind, row_index):
    """Return this row's detector LSF FWHM with unusable pixels interpolated.

    The LSF is produced upstream and is normally clean: new-oh-3 has zero bad
    pixels in 17260 rows x 3 arms.  The gaia1over100 selection exposed a
    different set of fibres, 9 of which carry a single LSF pixel of exactly 0.0
    at a spectrograph arm join -- 8 at 5800.0 A (b/r) and 1 at 7570.0 A (r/z),
    10 bad pixels in 179 million.  ``MoonZodiPhysicalModel.predict`` rightly
    rejects a non-positive FWHM, which killed the whole worker chunk.

    Repairing by interpolation is safe HERE specifically because the caller
    wants scalar band integrals (a moon fraction and a zodi total) out of the
    prediction, so a one-pixel correction at an arm edge cannot move them
    measurably.  Do not reuse this to paper over a genuinely broken LSF: if a
    row has no usable pixels at all it is returned as None so the caller can
    skip the prior rather than fit a fabricated one.
    """
    lsf = np.asarray(_WORKER_LSF[kind][row_index], dtype=np.float64)
    good = np.isfinite(lsf) & (lsf > 0.0)
    if good.all():
        return lsf
    if not good.any():
        return None
    idx = np.arange(lsf.size)
    repaired = lsf.copy()
    repaired[~good] = np.interp(idx[~good], idx[good], lsf[good])
    _LSF_REPAIR_COUNT["rows"] += 1
    _LSF_REPAIR_COUNT["pixels"] += int((~good).sum())
    if not _LSF_REPAIR_COUNT["reported"]:
        _LSF_REPAIR_COUNT["reported"] = True
        _bad_w = _WORKER_DECOMPOSER.wave[~good]
        print(
            f"  [lsf-repair] {kind} row {row_index}: interpolated "
            f"{int((~good).sum())} non-positive/non-finite LSF pixel(s) at "
            f"{', '.join(f'{x:.1f}' for x in _bad_w[:4])} A"
            f"{' ...' if _bad_w.size > 4 else ''}.  Further repairs in this "
            f"worker are silent; the total is reported at the end.",
            flush=True,
        )
    return repaired


def _install_split_zodi_amplitude_prior(kind, row_index):
    """Install this spectrum's geometry amplitude prior before fitting.

    The split-zodi priors are per-spectrum: the moon-share bracket and the
    absolute Leinert zodi bracket both come from a geometry prediction for this
    exposure and this telescope.  Only scalars are installed, so the design
    matrix is untouched and the basis stays identical for every row -- which is
    what lets the ML side reconstruct from coefficients with a single
    decomposer.

    Geometry that cannot be modelled (target below the horizon) clears the
    prior for that row instead of failing it: the fit then falls back to the
    shape bounds alone, which is exactly the pre-prior behaviour.

    The zodi total is recentred by SPLIT_ZODI_ZODI_PRIOR_CALIBRATION before it
    is installed; see that constant for the measurement.  The moon fraction is
    passed through untouched -- it is calibration-free by construction, so
    scaling it here would corrupt the one constraint that is not.
    """
    from skysub.sky_decomp.moon_zodi_model import (
        MoonZodiInvalidObservationError,
        geometry_amplitude_prior,
    )

    if _WORKER_META is None or not _WORKER_LSF:
        return
    _lsf = _sanitised_lsf_row(kind, row_index)
    if _lsf is None:
        # No usable LSF anywhere in this row: fall back to the shape bounds
        # alone, exactly as for geometry that cannot be modelled.
        _WORKER_DECOMPOSER.set_amplitude_prior(None, None)
        return
    try:
        fraction, zodi_total, _target_airmass = geometry_amplitude_prior(
            _WORKER_DECOMPOSER.wave,
            _lsf,
            _moon_zodi_observation(kind, row_index),
            physical_to_fit_flux_scale=float(_WORKER_FACTOR),
        )
    except MoonZodiInvalidObservationError:
        _WORKER_DECOMPOSER.set_amplitude_prior(None, None)
        return
    _WORKER_DECOMPOSER.set_amplitude_prior(
        fraction, zodi_total * SPLIT_ZODI_ZODI_PRIOR_CALIBRATION
    )


def _science_line_mask_for_row(row_index):
    """Science-line mask for one row, centred on its measured Halpha velocity.

    Falls back to the worker's static rest-frame mask whenever centring is
    disabled, the science flux is unavailable, or the velocity cannot be
    measured -- ``measure_halpha_velocity`` returns 0.0 in that case, which
    reproduces the static mask exactly.  The window WIDTH never changes, so
    the number of masked pixels is constant to within a pixel or two and the
    dof count stays comparable from row to row.
    """
    if _WORKER_SCIENCE_LINE_MASK is None:
        return None
    if not _WORKER_SCIENCE_LINE_CENTRE:
        return _WORKER_SCIENCE_LINE_MASK
    sci = _WORKER_FLUX.get("sci")
    if sci is None:
        return _WORKER_SCIENCE_LINE_MASK
    sky = _WORKER_FLUX.get("sky1")
    if sky is None:
        sky = _WORKER_FLUX.get("sky2")
    velocity = measure_halpha_velocity(
        _WORKER_DECOMPOSER.wave,
        np.asarray(sci[row_index], dtype=np.float64),
        None if sky is None else np.asarray(sky[row_index], dtype=np.float64),
    )
    if velocity == 0.0:
        return _WORKER_SCIENCE_LINE_MASK
    mask, _ = science_line_mask(
        _WORKER_DECOMPOSER.wave,
        _WORKER_SCIENCE_LINE_FWHM,
        centre_velocity_km_s=velocity,
    )
    return mask


def fit_chunk_worker(args):
    """Fit one chunk of spectra using the worker-local SkyDecomp instance."""
    global _WORKER_DECOMPOSER, _WORKER_FACTOR, _WORKER_PROGRESS_QUEUE, _WORKER_FIT_MODEL, \
        _WORKER_SCIENCE_LINE_MASK, _WORKER_SCIENCE_LINE_FWHM, \
        _WORKER_SCIENCE_LINE_CENTRE
    if _WORKER_DECOMPOSER is None:
        raise RuntimeError("Worker SkyDecomp has not been initialised.")
    kind, idx0, idx1 = args
    flux_chunk = np.asarray(_WORKER_FLUX[kind][idx0:idx1], dtype=np.float64)
    out = []
    for j in range(flux_chunk.shape[0]):
        idx = idx0 + j
        flux_row = flux_chunk[j] * _WORKER_FACTOR
        ivar_row = np.ones_like(flux_row)
        # Science emission lines: excluded from every arm, not just the science
        # one.  The sky fibres also sit on Galactic diffuse ionised gas -- on
        # row 773 the near arm carries Halpha EW 10.0 A -- and masking only the
        # science arm would make the three fits differ in which pixels they
        # used, which is exactly the asymmetry the ML transfer is trying to
        # measure.  The velocity is likewise measured once from the SCIENCE
        # fibre and shared by all three arms, so the masked pixel set stays
        # identical across them.
        if _WORKER_SCIENCE_LINE_MASK is not None:
            ivar_row[_science_line_mask_for_row(idx)] = 0.0
        if _WORKER_FIT_MODEL == "baseline":
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                verbose=False,
                n_lsf_refits=3,
            )
        elif _WORKER_FIT_MODEL == "lsf-surface-iterative":
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                verbose=False,
            )
        elif _WORKER_FIT_MODEL == SPLIT_ZODI_FIT_MODEL:
            _install_split_zodi_amplitude_prior(kind, idx)
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                verbose=False,
            )
        elif _WORKER_FIT_MODEL == "moon-zodi-lsf-surface-iterative":
            # Preserve invalid source pixels; zero IVAR excludes them without
            # interpolating, imputing, cropping, or changing the native grid.
            ivar_row = np.isfinite(flux_row).astype(np.float64)
            if _WORKER_SCIENCE_LINE_MASK is not None:
                ivar_row[_science_line_mask_for_row(idx)] = 0.0
            # Same one-pixel arm-join LSF holes as the split-zodi prior path.
            # Here the LSF drives the CONVOLUTION, not just a scalar prior, so
            # a row with no usable LSF at all cannot be fitted in this mode --
            # _sanitised_lsf_row returns None and the row is failed explicitly
            # rather than fitted against a fabricated LSF.
            lsf_row = _sanitised_lsf_row(kind, idx)
            if lsf_row is None:
                raise ValueError(
                    f"row {idx} ({kind}) has no finite positive LSF pixel; "
                    f"cannot fit with {MOON_ZODI_FIT_MODEL}")
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                observation=_moon_zodi_observation(kind, idx),
                detector_lsf_fwhm=lsf_row,
                verbose=False,
            )
        else:
            raise RuntimeError(f"Worker has unsupported fit model: {_WORKER_FIT_MODEL}")
        out.append((idx, result))
        if _WORKER_PROGRESS_QUEUE is not None:
            _WORKER_PROGRESS_QUEUE.put(1)
    return kind, out


def resolve_base_dir(path_arg):
    """Accept either the project base dir or the palace dir and return SkyDecomp base_dir."""
    path = Path(path_arg).expanduser().resolve()

    if (path / "palace" / "PMD").exists():
        return path
    if path.name == "palace" and (path / "PMD").exists():
        return path.parent

    raise FileNotFoundError(
        "Could not resolve a valid SkyDecomp base directory from "
        f"{path}. Expected either a base dir containing palace/PMD or the palace directory itself."
    )


def resolve_runtime_data_roots(
    fit_model,
    palace_dir=None,
    moon_zodi_data_root=None,
):
    """Resolve only the data contract used by the selected fit model."""
    if fit_model == MOON_ZODI_FIT_MODEL:
        candidate = (
            moon_zodi_data_root
            if moon_zodi_data_root is not None
            else palace_dir
            if palace_dir is not None
            else DEFAULT_MOON_ZODI_DATA_ROOT
        )
        data_root = Path(candidate).expanduser().resolve()
        validate_decomposition_data_root(str(data_root))
        return data_root, data_root

    if fit_model == SPLIT_ZODI_FIT_MODEL and (
        moon_zodi_data_root is not None or palace_dir is None
    ):
        candidate = (
            moon_zodi_data_root
            if moon_zodi_data_root is not None
            else DEFAULT_MOON_ZODI_DATA_ROOT
        )
        data_root = Path(candidate).expanduser().resolve()
        validate_decomposition_data_root(str(data_root))
        return data_root, data_root

    if fit_model == SPLIT_ZODI_FIT_MODEL and palace_dir is not None:
        candidate = Path(palace_dir).expanduser().resolve()
        if (candidate / "bundle_manifest.json").is_file():
            validate_decomposition_data_root(str(candidate))
            return candidate, candidate

    if palace_dir is None:
        raise ValueError(
            "palace_dir is required for baseline and non-split "
            "lsf-surface-iterative fits"
        )
    return resolve_base_dir(palace_dir), None


def _iter_chunk_tasks(n_rows, chunk_size):
    for kind in ("sci", "sky1", "sky2"):
        for i0 in range(0, n_rows, chunk_size):
            i1 = min(i0 + chunk_size, n_rows)
            # Sending only index ranges keeps per-task IPC tiny.
            yield (kind, i0, i1)


def run(
    data_file,
    palace_dir,
    n_workers,
    lsf_sigma,
    factor,
    output_dir,
    chunk_size,
    max_in_flight,
    fit_model="baseline",
    n_refinement_cycles=5,
    limit=None,
    pin_workers=False,
    diagnose_threads=False,
    palace_suffix=None,
    palace_oh_suffix=None,
    palace_diffuse_suffix=None,
    exposure_seconds=900.0,
    moon_zodi_data_root=None,
    n_spline_knots=MOON_N_KNOTS_DEFAULT,
    n_zodi_spline_knots=SPLIT_ZODI_N_KNOTS_DEFAULT,
    zodi_smooth_lambda=SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT,
    diffuse_ratio_bound_dex=SPLIT_ZODI_DIFFUSE_RATIO_BOUND_DEX,
    diffuse_ratio_nominal=SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL,
    diffuse_oh_centre_log10=SPLIT_ZODI_DIFFUSE_OH_CENTRE_LOG10,
    diffuse_oh_bound_dex=SPLIT_ZODI_DIFFUSE_OH_BOUND_DEX,
    mask_science_lines=SCIENCE_LINE_MASK_ENABLED,
    centre_on_halpha=SCIENCE_LINE_MASK_CENTRE_ON_HALPHA,
):
    base_dir, resolved_moon_zodi_data_root = resolve_runtime_data_roots(
        fit_model,
        palace_dir=palace_dir,
        moon_zodi_data_root=moon_zodi_data_root,
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {data_file} ...")
    wave = fits.getdata(data_file, "WAVE").astype(np.float64)
    with fits.open(data_file, memmap=True) as hdul:
        n_rows_total = int(hdul["FLUX_SCI"].data.shape[0])

    n_rows = n_rows_total if limit is None else min(n_rows_total, int(limit))

    if limit is not None and n_rows < n_rows_total:
        print(f"  {n_rows} spectra (limited from {n_rows_total}), {len(wave)} wavelength pixels")
    else:
        print(f"  {n_rows} spectra, {len(wave)} wavelength pixels")
    print(f"  n_workers={n_workers}, lsf_sigma={lsf_sigma}, factor={factor}")
    print(f"  chunk_size={chunk_size}, max_in_flight={max_in_flight}")
    print(f"  fit_model={fit_model}, n_refinement_cycles={n_refinement_cycles}")
    print(f"  palace_suffix={palace_suffix!r}")
    print(f"  palace_oh_suffix={palace_oh_suffix!r}")
    print(f"  palace_diffuse_suffix={palace_diffuse_suffix!r}")
    print(f"  exposure_seconds_fallback={exposure_seconds}")
    if resolved_moon_zodi_data_root is not None:
        print(f"  bundled_data_root={str(resolved_moon_zodi_data_root)!r}")
    if fit_model in ("lsf-surface-iterative", SPLIT_ZODI_FIT_MODEL):
        print(f"  n_spline_knots={n_spline_knots} "
              f"(Moon_bs basis = {int(n_spline_knots) + 4})")
    if fit_model == SPLIT_ZODI_FIT_MODEL:
        print(f"  n_zodi_spline_knots={n_zodi_spline_knots}, "
              f"zodi_smooth_lambda={zodi_smooth_lambda}")
    print(f"  pin_workers={pin_workers}, diagnose_threads={diagnose_threads}")
    print(f"  base_dir={base_dir}")

    n_tasks = int(np.ceil(n_rows / chunk_size)) * 3
    results = {kind: [None] * n_rows for kind in ("sci", "sky1", "sky2")}
    completed = 0

    t0 = time.perf_counter()

    # spawn everywhere: `fork` inherits parent BLAS pools and undermines thread limits.
    mp_context = mp.get_context("spawn")
    progress_queue = mp_context.Queue()
    worker_counter = mp_context.Value("i", 0)

    def _drain_progress_queue():
        increment = 0
        while True:
            try:
                increment += progress_queue.get_nowait()
            except queue_mod.Empty:
                break
        if increment:
            pbar.update(increment)

    with ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=mp_context,
        initializer=init_worker,
        initargs=(
            wave,
            lsf_sigma,
            str(base_dir),
            float(factor),
            str(data_file),
            progress_queue,
            fit_model,
            n_refinement_cycles,
            worker_counter,
            bool(pin_workers),
            bool(diagnose_threads),
            palace_suffix,
            palace_oh_suffix,
            palace_diffuse_suffix,
            exposure_seconds,
            (
                None
                if resolved_moon_zodi_data_root is None
                else str(resolved_moon_zodi_data_root)
            ),
            int(n_spline_knots),
            int(n_zodi_spline_knots),
            float(zodi_smooth_lambda),
            bool(mask_science_lines),
            bool(centre_on_halpha),
            float(diffuse_ratio_bound_dex),
            diffuse_ratio_nominal,
            diffuse_oh_centre_log10,
            float(diffuse_oh_bound_dex),
        ),
    ) as executor:
        pbar = tqdm(
            total=3 * n_rows,
            desc="decomp",
            unit=" decomp",
            mininterval=0.2,
            position=0,
            leave=True,
        )
        pbar.set_postfix(chunks=f"0/{n_tasks}")
        pbar.refresh()

        task_iter = iter(_iter_chunk_tasks(n_rows, chunk_size))
        pending = set()

        def _submit_until_full():
            while len(pending) < max_in_flight:
                try:
                    task = next(task_iter)
                except StopIteration:
                    return
                pending.add(executor.submit(fit_chunk_worker, task))

        _submit_until_full()
        while pending:
            done, pending = wait(pending, timeout=0.2, return_when=FIRST_COMPLETED)
            _drain_progress_queue()
            for future in done:
                kind, chunk_results = future.result()
                for idx, result in chunk_results:
                    results[kind][idx] = result
                completed += 1
                pbar.set_postfix(chunks=f"{completed}/{n_tasks}")
            _submit_until_full()
        _drain_progress_queue()
        pbar.close()

    progress_queue.close()
    progress_queue.join_thread()

    elapsed = time.perf_counter() - t0
    print(f"Fitting done in {elapsed:.1f}s ({elapsed / n_rows:.2f}s per spectrum)")

    stem = Path(data_file).stem
    suffix = FIT_MODEL_SUFFIXES[fit_model]
    for kind in ("sci", "sky1", "sky2"):
        results_to_fits(results[kind], output_dir / f"{stem}_decomp_{kind}{suffix}.fits")


def _copy_hdu_with_name(hdu, extname):
    """Return a copy of an HDU with a new extension name."""
    header = hdu.header.copy()
    header["EXTNAME"] = extname
    return type(hdu)(data=hdu.data, header=header, name=extname)


def _infer_decomp_label(path, index):
    """Infer a stable label for a decomposition file from its filename."""
    name = Path(path).name.lower()
    for label in ("sky1", "sky2", "sci"):
        if label in name:
            return label.upper()
    return f"DEC{index}"


def extract_meta_and_coef_products(
    input_fits_path,
    decomp_fits_path_1,
    decomp_fits_path_2,
    decomp_fits_path_3,
    meta_output_path=None,
    sky1_output_path=None,
    sky2_output_path=None,
    sci_output_path=None,
):
    """Write compact FITS products containing only selected extensions.

    The first output contains only the META extension from `input_fits_path`.
    Each decomposition input gets its own output FITS containing META and COEF.
    Extended LSF-surface products are copied when present. Default output
    paths are written in the current working directory.
    """
    input_path = Path(input_fits_path)
    cwd = Path.cwd()
    if meta_output_path is None:
        meta_output_path = str(cwd / f"{input_path.stem}_meta_only{input_path.suffix}")

    decomp_files = [decomp_fits_path_1, decomp_fits_path_2, decomp_fits_path_3]
    decomp_outputs = [sky1_output_path, sky2_output_path, sci_output_path]

    with fits.open(input_fits_path) as hdul_in:
        if "META" not in hdul_in:
            raise KeyError(f"Missing META extension in {input_fits_path}")
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                _copy_hdu_with_name(hdul_in["META"], "META"),
            ]
        ).writeto(meta_output_path, overwrite=True)

    resolved_outputs = []
    for index, decomp_path in enumerate(decomp_files, start=1):
        label = _infer_decomp_label(decomp_path, index)
        out_path = decomp_outputs[index - 1]
        if out_path is None:
            stem_lower = Path(decomp_path).stem.lower()
            if "moon_zodi_lsf_surface_iterative" in stem_lower:
                variant = "_moon_zodi_lsf_surface_iterative"
            elif "lsf_surface_iterative" in stem_lower:
                variant = "_lsf_surface_iterative"
            else:
                variant = ""
            out_path = str(
                cwd / (f"{input_path.stem}_{label.lower()}_meta_coef{variant}{input_path.suffix}")
            )
        with fits.open(decomp_path) as hdul_dec:
            for extname in ("META", "COEF"):
                if extname not in hdul_dec:
                    raise KeyError(f"Missing {extname} extension in {decomp_path}")
            compact_hdus = [
                fits.PrimaryHDU(),
                _copy_hdu_with_name(hdul_dec["META"], "META"),
                _copy_hdu_with_name(hdul_dec["COEF"], "COEF"),
            ]
            if "COEF_ERR" in hdul_dec:
                compact_hdus.append(_copy_hdu_with_name(hdul_dec["COEF_ERR"], "COEF_ERR"))
            for cov_name in ("COEF_COV_MOON", "COEF_COV_ZODI"):
                if cov_name in hdul_dec:
                    compact_hdus.append(_copy_hdu_with_name(hdul_dec[cov_name], cov_name))
            lsf_extensions = ("LSF_COEF", "LSF_KNOTS", "LSF_META")
            present = [name in hdul_dec for name in lsf_extensions]
            if any(present) and not all(present):
                raise KeyError(f"Incomplete LSF-surface extensions in {decomp_path}")
            if all(present):
                compact_hdus.extend(
                    _copy_hdu_with_name(hdul_dec[name], name) for name in lsf_extensions
                )
            moon_zodi_extensions = ("MZ_MODEL", "MZ_ASSETS", "MZ_KNOTS", "MZ_META")
            moon_zodi_present = [name in hdul_dec for name in moon_zodi_extensions]
            if any(moon_zodi_present) and not all(moon_zodi_present):
                raise KeyError(f"Incomplete Moon/Zodi extensions in {decomp_path}")
            if all(moon_zodi_present):
                compact_hdus.extend(
                    _copy_hdu_with_name(hdul_dec[name], name)
                    for name in moon_zodi_extensions
                )
            fits.HDUList(compact_hdus).writeto(out_path, overwrite=True)
        print(f"Wrote {label} META/COEF file -> {out_path}")
        resolved_outputs.append(out_path)

    print(f"Wrote META-only file -> {meta_output_path}")
    return (meta_output_path, *resolved_outputs)

def thin_fits_every_n(input_path, output_path, n, row_hdu_name="META"):
    """Write a new FITS with every n-th row-like element kept.

    The function preserves HDU structure and headers. It identifies the row
    count from `row_hdu_name` (default: META), then slices any table HDU with
    that row count and any image HDU whose first axis matches that row count.

    Tables that reference the row axis indirectly through a `spectrum_index`
    column (such as `LSF_META`, which has one row per channel per spectrum)
    are filtered to rows whose `spectrum_index` is in the kept set, and the
    `spectrum_index` values are remapped to the new 0-based row positions so
    downstream loaders can still index a thinned coefficient cube directly.
    """
    if n < 1:
        raise ValueError("n must be >= 1")

    with fits.open(input_path) as hdul:
        if row_hdu_name not in hdul:
            raise KeyError(f"HDU '{row_hdu_name}' not found in {input_path}")

        n_rows = len(hdul[row_hdu_name].data)
        indices = np.arange(n_rows, dtype=int)[::n]
        remap = {int(old): int(new) for new, old in enumerate(indices)}
        keep = slice(None, None, n)
        global_hdus = {"MZ_MODEL", "MZ_ASSETS", "MZ_KNOTS"}

        out_hdus = []
        for hdu in hdul:
            header = hdu.header.copy()

            if isinstance(hdu, fits.PrimaryHDU):
                data = hdu.data
                if data is not None and getattr(data, "ndim", 0) >= 1 and data.shape[0] == n_rows:
                    data = data[keep, ...]
                out_hdus.append(fits.PrimaryHDU(data=data, header=header))

            elif isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                data = hdu.data
                if hdu.name == "LSF_META" and data is not None:
                    selected = np.isin(np.asarray(data["spectrum_index"], dtype=int), indices)
                    data = data[selected].copy()
                    for row in data:
                        row["spectrum_index"] = remap[int(row["spectrum_index"])]
                elif hdu.name == "MZ_META" and data is not None:
                    data = data[keep].copy()
                    data["spectrum_index"] = np.arange(len(data), dtype=int)
                elif hdu.name not in global_hdus and data is not None and len(data) == n_rows:
                    data = data[keep]
                elif (
                    data is not None
                    and "spectrum_index" in data.dtype.names
                    and len(data) % n_rows == 0
                ):
                    # Multi-row-per-spectrum table (e.g. LSF_META has one row
                    # per channel per spectrum). Filter by spectrum_index and
                    # remap to the thinned cube's 0-based positions.
                    si = np.asarray(data["spectrum_index"], dtype=np.int64)
                    mask = np.isin(si, indices)
                    data = data[mask].copy()
                    remapped = np.fromiter(
                        (remap[int(v)] for v in data["spectrum_index"]),
                        dtype=np.int64,
                        count=len(data),
                    )
                    data["spectrum_index"] = remapped
                out_hdus.append(type(hdu)(data=data, header=header, name=hdu.name))

            elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                data = hdu.data
                if (
                    hdu.name not in global_hdus
                    and data is not None
                    and getattr(data, "ndim", 0) >= 1
                    and data.shape[0] == n_rows
                ):
                    data = data[keep, ...]
                out_hdus.append(type(hdu)(data=data, header=header, name=hdu.name))

            else:
                out_hdus.append(hdu.copy())

        fits.HDUList(out_hdus).writeto(output_path, overwrite=True)


def main():
    parser = argparse.ArgumentParser(description="LVM sky spectral decomposition")
    parser.add_argument("data_file", help="Input FITS file (median stacked LVM frame)")
    parser.add_argument(
        "palace_dir",
        nargs="?",
        default=None,
        help=(
            "Legacy project/PALACE root, or a complete bundled data root containing "
            "bundle_manifest.json. Required by baseline and non-split "
            "lsf-surface-iterative; optional for split-zodi and Moon/Zodi modes."
        ),
    )
    parser.add_argument(
        "--n-workers", type=int, default=4, help="Number of parallel worker processes (default: 4)"
    )
    parser.add_argument(
        "--lsf-sigma", type=float, default=0.5, help="LSF Gaussian sigma in Å (default: 0.5)"
    )
    parser.add_argument(
        "--factor", type=float, default=1e14, help="Flux scaling factor (default: 1e14)"
    )
    parser.add_argument(
        "--chunk-size", type=int, default=64, help="Rows per worker task chunk (default: 64)"
    )
    parser.add_argument(
        "--max-in-flight",
        type=int,
        default=None,
        help="Max submitted chunks waiting/running at once (default: n-workers)",
    )
    parser.add_argument(
        "--output-dir", default=".", help="Output directory for result FITS files (default: .)"
    )
    parser.add_argument(
        "--fit-model",
        choices=tuple(FIT_MODEL_SUFFIXES),
        default="baseline",
        help="Fit implementation (default: baseline)",
    )
    parser.add_argument(
        "--n-refinement-cycles",
        type=int,
        default=5,
        help="Continuum/LSF/line cycles for lsf-surface-iterative (default: 5)",
    )
    parser.add_argument(
        "--palace-suffix",
        default=None,
        help=(
            "Optional suffix for versioned pmd_popmodel_OH and pmd_refcont files. "
            "The suffix is appended exactly; for example, '_adam_v1' selects "
            "pmd_popmodel_OH_adam_v1.dat and "
            "pmd_refcont_adam_v1.dat (legacy default: canonical unsuffixed "
            "files; bundled Moon/Zodi modes use their manifest defaults)."
        ),
    )
    parser.add_argument(
        "--palace-oh-suffix",
        default=None,
        help="Optional exact suffix overriding --palace-suffix for pmd_popmodel_OH only.",
    )
    parser.add_argument(
        "--palace-diffuse-suffix",
        default=None,
        help="Optional exact suffix overriding --palace-suffix for pmd_refcont only.",
    )
    parser.add_argument(
        "--exposure-seconds",
        type=float,
        default=900.0,
        help=(
            "Exposure duration used only when META has no exposure duration "
            "(v1 provenance contract requires the default 900 s assumption)."
        ),
    )
    parser.add_argument(
        "--moon-zodi-data-root",
        type=Path,
        default=None,
        help=(
            "Complete data root containing moon_zodi/ and palace/PMD for the "
            "Moon/Zodi or split-zodi method (default for both: packaged "
            "skysub/sky_decomp/data)."
        ),
    )
    parser.add_argument(
        "--n-spline-knots",
        type=int,
        default=MOON_N_KNOTS_DEFAULT,
        help=(
            f"Interior B-spline knots for Moon_bs when --fit-model is "
            f"lsf-surface-iterative or {SPLIT_ZODI_FIT_MODEL} "
            f"(default: {MOON_N_KNOTS_DEFAULT}; Moon_bs basis = knots + 4)."
        ),
    )
    parser.add_argument(
        "--n-zodi-spline-knots",
        type=int,
        default=SPLIT_ZODI_N_KNOTS_DEFAULT,
        help=(
            f"Interior B-spline knots for Zodi_bs when "
            f"--fit-model={SPLIT_ZODI_FIT_MODEL} (default: {SPLIT_ZODI_N_KNOTS_DEFAULT})."
        ),
    )
    parser.add_argument(
        "--zodi-smooth-lambda",
        type=float,
        default=SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT,
        help=(
            f"Curvature penalty on the Zodi_bs spline when "
            f"--fit-model={SPLIT_ZODI_FIT_MODEL} (default: {SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT})."
        ),
    )
    parser.add_argument(
        "--diffuse-ratio-bound-dex",
        type=float,
        default=SPLIT_ZODI_DIFFUSE_RATIO_BOUND_DEX,
        help=(
            "Half-width in dex of the bracket on the diffuse species ratios "
            "FeO/HO2 and O2Ac/HO2, about --diffuse-ratio-nominal. 0 disables. "
            "The three species are individually unidentifiable in the LVM band "
            "-- the three arms of one exposure disagree by 0.54 dex at the "
            "median -- so this removes freedom the data cannot measure. "
            "Requires --diffuse-ratio-nominal."
        ),
    )
    parser.add_argument(
        "--diffuse-oh-bound-dex",
        type=float,
        default=SPLIT_ZODI_DIFFUSE_OH_BOUND_DEX,
        help=(
            "Half-width in dex of the MOON-GATED upper bound on A_FeO / A_OH, "
            "above --diffuse-oh-centre-log10. 0 disables. FeO is mesospheric and "
            "cannot depend on the moon, but its amplitude relative to OH rises "
            "a factor 4.5 with moon_frac_po, so the template is absorbing "
            "scattered moonlight. One-sided and gated because the dark-time "
            "scatter of the ratio is 0.306 dex and real."
        ),
    )
    parser.add_argument(
        "--diffuse-oh-centre-log10",
        type=float,
        default=SPLIT_ZODI_DIFFUSE_OH_CENTRE_LOG10,
        help=(
            "log10(A_FeO / A_OH) centre for the bound above: the DARK-TIME "
            "median over FeO-live rows, where there is no moon to leak. "
            "Re-measure per corpus."
        ),
    )
    parser.add_argument(
        "--diffuse-ratio-nominal",
        type=str,
        default=",".join(str(v) for v in SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL),
        help=(
            "Comma-separated HO2,FeO,O2Ac FLUX shares giving the centre of the "
            "ratio bracket, e.g. '0.0396,0.7026,0.2578'. MEASURE THIS ON THE "
            "CORPUS BEING FITTED: centring on PALACE's own shares costs 25.7%% "
            "of the blue chi2 against 0.67%% for the corpus median, because "
            "PALACE is calibrated for Paranal and LVM observes from LCO."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N input rows (default: process all rows)",
    )
    parser.add_argument(
        "--pin-workers",
        action="store_true",
        help=(
            "Pin each worker to a single CPU core via sched_setaffinity (Linux). "
            "Nuclear option: hard-caps observed CPU usage per worker to 100%% even if a "
            "library ignores thread-limit env vars."
        ),
    )
    parser.add_argument(
        "--diagnose-threads",
        action="store_true",
        help="Print per-library thread pool counts from worker 0 after all imports finish.",
    )
    parser.add_argument(
        "--no-science-line-mask",
        action="store_true",
        help=(
            "Fit the science emission-line windows instead of excluding them. "
            "By default [OII]3727, Hbeta, [OIII]4959/5007, Halpha, [NII]6548/6583 "
            "and [SII]6716/6731 are masked via IVAR=0 in every arm, because none "
            "of them is in the basis and the fitted OH component absorbs them "
            "(measured 6.7x the sideband OH level inside Halpha/[NII] on an "
            "inner-plane H II region field).  Use this to A/B the mask."
        ),
    )
    parser.add_argument(
        "--no-halpha-centring",
        action="store_true",
        help=(
            "Keep the science-line windows at rest wavelength instead of "
            "sliding them to the per-row Halpha velocity.  Centring lifts core "
            "coverage on strong-Halpha rows from 89.0%% to 92.5%%; use this to "
            "A/B it.  Ignored when --no-science-line-mask is given."
        ),
    )
    parser.add_argument(
        "--only-thin",
        action="store_true",
        help=(
            "Skip decomposition and extract-compact steps; regenerate only the "
            "every10-thinned FITS from already-existing decomp files."
        ),
    )
    args = parser.parse_args()

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")
    if args.max_in_flight is None:
        args.max_in_flight = args.n_workers
    if args.max_in_flight < 1:
        raise ValueError("--max-in-flight must be >= 1")
    if args.n_refinement_cycles < 1:
        raise ValueError("--n-refinement-cycles must be >= 1")
    if args.limit is not None and args.limit < 1:
        raise ValueError("--limit must be >= 1")
    if not np.isfinite(args.exposure_seconds) or args.exposure_seconds <= 0.0:
        raise ValueError("--exposure-seconds must be positive and finite")
    if args.fit_model == MOON_ZODI_FIT_MODEL and args.exposure_seconds != 900.0:
        raise ValueError(
            "Moon/Zodi v1 records missing META exposure time as 'assumed_900s'; "
            "--exposure-seconds must therefore remain 900"
        )
    if (
        not args.only_thin
        and args.fit_model not in (MOON_ZODI_FIT_MODEL, SPLIT_ZODI_FIT_MODEL)
        and args.palace_dir is None
    ):
        parser.error(
            "palace_dir is required for baseline and non-split "
            "lsf-surface-iterative fits"
        )

    suffix = FIT_MODEL_SUFFIXES[args.fit_model]
    stem = Path(args.data_file).stem
    output_dir = Path(args.output_dir)

    if not args.only_thin:
        run(
            data_file=args.data_file,
            palace_dir=args.palace_dir,
            n_workers=args.n_workers,
            lsf_sigma=args.lsf_sigma,
            factor=args.factor,
            output_dir=args.output_dir,
            chunk_size=args.chunk_size,
            max_in_flight=args.max_in_flight,
            fit_model=args.fit_model,
            n_refinement_cycles=args.n_refinement_cycles,
            limit=args.limit,
            pin_workers=args.pin_workers,
            diagnose_threads=args.diagnose_threads,
            palace_suffix=args.palace_suffix,
            palace_oh_suffix=args.palace_oh_suffix,
            palace_diffuse_suffix=args.palace_diffuse_suffix,
            exposure_seconds=args.exposure_seconds,
            moon_zodi_data_root=args.moon_zodi_data_root,
            n_spline_knots=args.n_spline_knots,
            n_zodi_spline_knots=args.n_zodi_spline_knots,
            zodi_smooth_lambda=args.zodi_smooth_lambda,
            diffuse_ratio_bound_dex=args.diffuse_ratio_bound_dex,
            diffuse_oh_bound_dex=args.diffuse_oh_bound_dex,
            diffuse_oh_centre_log10=args.diffuse_oh_centre_log10,
            diffuse_ratio_nominal=(
                None if not args.diffuse_ratio_nominal
                else tuple(float(v) for v in args.diffuse_ratio_nominal.split(","))
            ),
            mask_science_lines=not args.no_science_line_mask,
            centre_on_halpha=not args.no_halpha_centring,
        )

        extract_meta_and_coef_products(
            input_fits_path=args.data_file,
            decomp_fits_path_1=output_dir / f"{stem}_decomp_sky1{suffix}.fits",
            decomp_fits_path_2=output_dir / f"{stem}_decomp_sky2{suffix}.fits",
            decomp_fits_path_3=output_dir / f"{stem}_decomp_sci{suffix}.fits",
            meta_output_path=output_dir / f"{stem}_meta_only.fits",
            sky1_output_path=output_dir / f"{stem}_sky1_meta_coef{suffix}.fits",
            sky2_output_path=output_dir / f"{stem}_sky2_meta_coef{suffix}.fits",
            sci_output_path=output_dir / f"{stem}_sci_meta_coef{suffix}.fits",
        )
    else:
        required = [
            Path(args.data_file),
            output_dir / f"{stem}_decomp_sci{suffix}.fits",
            output_dir / f"{stem}_decomp_sky1{suffix}.fits",
            output_dir / f"{stem}_decomp_sky2{suffix}.fits",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "--only-thin requires these files to already exist: " + ", ".join(missing)
            )
        print(
            "--only-thin: skipping decomposition and extract-compact; "
            "regenerating thinned files only"
        )

    for kind in ("sci", "sky1", "sky2"):
        thin_fits_every_n(
            output_dir / f"{stem}_decomp_{kind}{suffix}.fits",
            output_dir / f"{stem}_every10_decomp_{kind}{suffix}.fits",
            10,
        )
    thin_fits_every_n(
        args.data_file,
        output_dir / f"{stem}_every10.fits",
        10,
    )

if __name__ == "__main__":
    main()
