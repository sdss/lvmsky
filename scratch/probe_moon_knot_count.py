"""Moon_bs / Zodi_bs knot-count probe.

Runs SkyDecompLSFSurfaceIterative(split_zodi=True) on a small moon-regime
sample from moon_zodi_spline2/*_every10.fits, sweeping one component's
knot count while holding the other at the current deployed default.

Reports for each config:
  - flux-space pRMSE per row (median across sample)
  - median peak |coef| in the swept block
  - median adjacent-knot correlation across the sample
  - peak-SNR (peak |coef| / min sigma)
  - runtime

Regimes:
  bright  moon_alt > 30 deg AND moon_fli > 0.85    -- unambiguous moon
  faint   moon_alt < -5 deg                        -- moon below horizon;
                                                      Moon_bs should be near
                                                      zero; zodi should be
                                                      the dominant diffuse.
Select target with --sweep {moon,zodi} and regime with --regime {bright,faint}.
Defaults: --sweep moon --regime bright.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

REPO = Path("/Users/droryn/prog/lvm/lvmsky")
sys.path.insert(0, str(REPO / "skysub"))

from sky_decomp.lsf_surface_iterative import (  # noqa: E402
    LSFSurfaceIterativeConfig,
    SkyDecompLSFSurfaceIterative,
)

DATA_DIR = REPO / "skysub" / "moon_zodi_spline2"
STEM = "lvmsframe_median_stack_1.2.1_p40_p70_every10"
FITS_EVERY10 = DATA_DIR / f"{STEM}.fits"

N_ROWS_TO_TEST = 20
N_REFINEMENT_CYCLES = 3

# Deployed defaults (used for the non-swept component).
DEFAULT_MOON_KNOTS = 25
DEFAULT_ZODI_KNOTS = 3

# Sweep grids.
MOON_KNOTS_SWEEP = [11, 15, 25]
ZODI_KNOTS_SWEEP = [1, 3, 5]

# Prefix used to identify basis columns in res.design_names.
BLOCK_PREFIX = {"moon": "Moon_bs", "zodi": "Zodi_bs"}

PHYSICAL_TO_FIT_FLUX_SCALE = 1e14


def _load_wave_and_meta():
    with fits.open(FITS_EVERY10) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        meta = Table(hdul["META"].data)
    return wave, meta


def _pick_rows_by_regime(meta, n, regime):
    moon_alt = np.asarray(meta["moon_alt"], dtype=np.float64)
    moon_fli = np.asarray(meta["moon_fli"], dtype=np.float64) \
        if "moon_fli" in meta.colnames else None
    if moon_fli is None:
        moon_phase = np.asarray(meta["moon_phase"], dtype=np.float64)
        moon_fli = 0.5 * (1.0 - np.cos(np.deg2rad(moon_phase)))
    finite = np.isfinite(moon_alt) & np.isfinite(moon_fli)
    if regime == "bright":
        mask = finite & (moon_alt > 30.0) & (moon_fli > 0.85)
    elif regime == "faint":
        mask = finite & (moon_alt < -5.0)
    else:
        raise ValueError(f"unknown regime {regime!r}")
    candidates = np.flatnonzero(mask)
    print(f"regime={regime!r} candidate rows: {candidates.size}")
    if candidates.size == 0:
        raise RuntimeError(f"No rows found for regime={regime!r}")
    rng = np.random.default_rng(42)
    n = int(min(n, candidates.size))
    return np.sort(rng.choice(candidates, size=n, replace=False))


def _load_row_arrays(rows):
    with fits.open(FITS_EVERY10) as hdul:
        flux_sci = np.asarray(hdul["FLUX_SCI"].data[rows], dtype=np.float64)
        lsf_sci = np.asarray(hdul["LSF_SCI"].data[rows], dtype=np.float64)
        if "FLUX_SIGMA_TOTAL" in hdul:
            sigma = np.asarray(hdul["FLUX_SIGMA_TOTAL"].data[rows],
                                dtype=np.float64)
        else:
            sigma = np.abs(flux_sci) * 0.05 + 1.0  # crude fallback
    return flux_sci, lsf_sci, sigma


def _build_decomp(wave, lsf, n_moon_knots, n_zodi_knots):
    lsf_safe = np.where(np.isfinite(lsf) & (lsf > 0), lsf, np.nanmedian(lsf))
    return SkyDecompLSFSurfaceIterative(
        wave=wave,
        lsf_sigma=lsf_safe,
        n_spline_knots=n_moon_knots,
        base_dir=str(REPO / "skysub" / "sky_decomp" / "data"),
        palace_oh_suffix="_joint_v2_updated",
        palace_diffuse_suffix="_joint_native_adam_invsky_p2_10000iter",
        moon_smooth_lambda=1.0e-3,
        split_zodi=True,
        n_zodi_spline_knots=n_zodi_knots,
        zodi_smooth_lambda=1.0e-1,
        moon_albedo_fiducial_phase_deg=30.0,
        zodi_color_exponent=0.26,
        config=LSFSurfaceIterativeConfig(
            n_refinement_cycles=N_REFINEMENT_CYCLES),
    )


def _fit_row(decomp, obs_native, sigma_native):
    obs = obs_native * PHYSICAL_TO_FIT_FLUX_SCALE
    sigma = np.asarray(sigma_native, dtype=np.float64)
    sigma_pos = sigma[np.isfinite(sigma) & (sigma > 0)]
    sigma_floor = 0.05 * float(np.median(sigma_pos)) if sigma_pos.size else 0.0
    sigma_eff = np.where(np.isfinite(sigma),
                          np.maximum(sigma, sigma_floor), np.nan)
    ivar = np.where((sigma_eff > 0) & np.isfinite(sigma_eff),
                     1.0 / (sigma_eff ** 2), 0.0)
    ivar = np.where(np.isfinite(obs), ivar, 0.0)
    result = decomp.fit(obs, ivar, verbose=False)
    return obs, result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", choices=("moon", "zodi"), default="moon")
    ap.add_argument("--regime", choices=("bright", "faint"), default="bright")
    ap.add_argument("--n-rows", type=int, default=N_ROWS_TO_TEST)
    args = ap.parse_args()

    sweep = args.sweep
    prefix = BLOCK_PREFIX[sweep]
    knot_grid = MOON_KNOTS_SWEEP if sweep == "moon" else ZODI_KNOTS_SWEEP
    fixed_knots_label = ("Zodi_bs (n_zodi_spline_knots="
                         f"{DEFAULT_ZODI_KNOTS})" if sweep == "moon"
                         else f"Moon_bs (n_spline_knots={DEFAULT_MOON_KNOTS})")

    wave, meta = _load_wave_and_meta()
    print(f"WAVE_A: {wave.size} pixels, {wave[0]:.1f} -> {wave[-1]:.1f} A")
    print(f"sweeping {prefix}, holding {fixed_knots_label} at deployed default")
    rows = _pick_rows_by_regime(meta, args.n_rows, args.regime)
    print(f"picked {rows.size} {args.regime}-moon rows: {rows.tolist()}")
    flux_sci, lsf_sci, sigma_sci = _load_row_arrays(rows)

    summary_rows = []
    for K in knot_grid:
        n_basis = K + 4  # cubic B-spline: n_interior_knots + degree + 1
        n_moon_knots = K if sweep == "moon" else DEFAULT_MOON_KNOTS
        n_zodi_knots = K if sweep == "zodi" else DEFAULT_ZODI_KNOTS
        print(f"\n=== {sweep}: {K} interior knots  ({prefix} count = "
              f"{n_basis}) ===")
        t0 = time.perf_counter()

        block_coef_matrix = []       # (n_row, n_basis)
        block_sigma_matrix = []      # (n_row, n_basis)
        flux_pRMSE = np.full(rows.size, np.nan)

        for i, r in enumerate(rows):
            decomp = _build_decomp(wave, lsf_sci[i], n_moon_knots,
                                     n_zodi_knots)
            try:
                obs, res = _fit_row(decomp, flux_sci[i], sigma_sci[i])
            except Exception as e:
                print(f"  row {int(r)}: fit failed ({type(e).__name__}: {e})")
                block_coef_matrix.append(np.full(n_basis, np.nan))
                block_sigma_matrix.append(np.full(n_basis, np.nan))
                continue
            residual = obs - res.bestfit_lsf
            flux_pRMSE[i] = float(
                np.sqrt(np.nanmean(residual ** 2))) / PHYSICAL_TO_FIT_FLUX_SCALE

            _names = res.design_names
            block_idx_local = np.array(
                [j for j, n in enumerate(_names) if n.startswith(prefix)],
                dtype=int)
            _c = np.asarray(res.coef[block_idx_local], dtype=np.float64)
            _s = np.asarray(res.coef_err[block_idx_local], dtype=np.float64)
            _s = np.where(np.isfinite(_s), _s, np.inf)
            block_coef_matrix.append(_c)
            block_sigma_matrix.append(_s)

            if (i + 1) % 5 == 0:
                print(f"  ... {i + 1}/{rows.size} fits done in "
                      f"{time.perf_counter() - t0:.1f} s")

        block_coef_matrix = np.vstack(block_coef_matrix)
        block_sigma_matrix = np.vstack(block_sigma_matrix)
        elapsed = time.perf_counter() - t0

        finite_prmse = flux_pRMSE[np.isfinite(flux_pRMSE)]
        med_prmse = float(np.nanmedian(finite_prmse)) if finite_prmse.size \
            else float("nan")
        p95_prmse = float(np.nanpercentile(finite_prmse, 95)) \
            if finite_prmse.size else float("nan")

        _c_abs = np.abs(block_coef_matrix)
        _peak = np.nanmax(_c_abs, axis=1)
        med_peak_abs_coef = float(np.nanmedian(_peak[np.isfinite(_peak)]))
        _min_sigma = np.min(np.where(block_sigma_matrix > 0,
                                       block_sigma_matrix, np.inf), axis=1)
        _min_sigma_finite = np.where(np.isfinite(_min_sigma), _min_sigma,
                                       np.nan)
        peak_snr = _peak / _min_sigma_finite
        med_peak_snr = float(np.nanmedian(peak_snr[np.isfinite(peak_snr)]))

        adj_corrs = np.full(block_coef_matrix.shape[1] - 1, np.nan)
        for k in range(block_coef_matrix.shape[1] - 1):
            a = block_coef_matrix[:, k]
            b = block_coef_matrix[:, k + 1]
            m = np.isfinite(a) & np.isfinite(b)
            if m.sum() >= 5 and a[m].std() > 0 and b[m].std() > 0:
                adj_corrs[k] = float(np.corrcoef(a[m], b[m])[0, 1])
        adj_corrs_valid = adj_corrs[np.isfinite(adj_corrs)]
        med_adj_abs_corr = float(np.nanmedian(np.abs(adj_corrs_valid))) \
            if adj_corrs_valid.size else float("nan")
        n_high_corr = int(np.sum(np.abs(adj_corrs_valid) > 0.90))

        summary_rows.append(dict(
            n_knots=K,
            n_basis=n_basis,
            n_rows_fit=int(finite_prmse.size),
            elapsed_s=float(elapsed),
            med_flux_pRMSE=med_prmse,
            p95_flux_pRMSE=p95_prmse,
            med_peak_SNR=med_peak_snr,
            med_peak_abs_coef=med_peak_abs_coef,
            med_adj_abs_corr=med_adj_abs_corr,
            n_pairs_corr_gt_090=n_high_corr,
            n_pairs=int(adj_corrs_valid.size),
        ))
        print(f"  {K:>2d} knots ({n_basis:>2d} basis): "
              f"med_pRMSE={med_prmse:.3g}, "
              f"peak_|coef|={med_peak_abs_coef:.3g}, "
              f"peak_SNR={med_peak_snr:.2f}, "
              f"med_adj_|corr|={med_adj_abs_corr:.3f}, "
              f"{n_high_corr}/{adj_corrs_valid.size} pairs > 0.90")

    print(f"\n=== SUMMARY  sweep={sweep}  regime={args.regime} ===")
    print("n_knots  n_basis  n_ok  time_s   med_pRMSE  p95_pRMSE   "
          "peak|coef|   peak_SNR  adj_|corr|  n>0.90")
    for s in summary_rows:
        print(f"{s['n_knots']:>7d}  {s['n_basis']:>7d}  "
              f"{s['n_rows_fit']:>4d}  {s['elapsed_s']:>6.1f}  "
              f"{s['med_flux_pRMSE']:>10.4g}  {s['p95_flux_pRMSE']:>10.4g}  "
              f"{s['med_peak_abs_coef']:>10.4g}  "
              f"{s['med_peak_SNR']:>8.2f}  {s['med_adj_abs_corr']:>10.3f}  "
              f"{s['n_pairs_corr_gt_090']:>3d}/{s['n_pairs']:>3d}")


if __name__ == "__main__":
    main()
