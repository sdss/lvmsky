"""Fit the physics-to-data scale factors k_moon, k_zodi from an existing corpus.

The `SkyDecompLSFSurfaceIterative` QP already accepts scalar amplitude priors
(:mod:`sky_decomp.moon_zodi_priors` for the physics side). What is missing is
the corpus-dependent multiplier that converts the KS moon proxy and the
Leinert V-band lookup to units matched to ``sum_i w_i c_i`` on the fitted
splines. This script produces that JSON.

Usage
-----

    python skysub/calibrate_moon_zodi_priors.py \\
        skysub/moon_zodi_spline/lvmsframe_median_stack_1.2.1_p40_p70.fits \\
        skysub/moon_zodi_spline \\
        skysub/sky_decomp/data \\
        --suffix _lsf_surface_iterative_split_zodi \\
        --arm sci \\
        --output skysub/moon_zodi_spline/moon_zodi_priors_calibration.json

The input FITS supplies the META columns needed to build a
:class:`MoonZodiObservation` per row; the ``corpus_dir`` locates the
existing decomposition products with the given ``suffix`` and per-arm
``kind`` (``sci`` / ``sky1`` / ``sky2``). The ``data_root`` locates the
Leinert lookup and ephemeris.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


SKYSUB = Path(__file__).resolve().parent
sys.path.insert(0, str(SKYSUB))
from sky_decomp.fit import SkyDecomp
from sky_decomp.moon_zodi_priors import (
    compute_geometry_proxies,
    fit_calibration,
)
from sky_decomp.moon_zodi_model import MoonZodiObservation


ARM_ROLES = {
    "sci": ("sci", "sci_ra", "sci_dec"),
    "sky1": ("sky_near", "sky_near_ra", "sky_near_dec"),
    "sky2": ("sky_far", "sky_far_ra", "sky_far_dec"),
}


def _text(v):
    return v.decode().strip() if isinstance(v, bytes) else str(v).strip()


def _make_observations(meta, arm: str) -> list[MoonZodiObservation]:
    role, ra_col, dec_col = ARM_ROLES[arm]
    obs = []
    for row in meta:
        exposure = None
        for column in ("exposure_seconds", "exptime"):
            if column in meta.dtype.names:
                candidate = float(row[column])
                if np.isfinite(candidate) and candidate > 0.0:
                    exposure = candidate
                    break
        if exposure is None:
            exposure = 900.0
            exposure_source = "assumed_900s"
        else:
            exposure_source = "metadata"
        obs.append(MoonZodiObservation(
            expnum=int(row["expnum"]),
            date_obs=_text(row["date_obs"]),
            role=role,
            target_ra_deg=float(row[ra_col]),
            target_dec_deg=float(row[dec_col]),
            exposure_seconds=float(exposure),
            exposure_seconds_source=exposure_source,
        ))
    return obs


def _data_amps(corpus_path: Path, wave: np.ndarray, data_root: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return per-row ``sum_i w_i c_i`` for moon and zodi from an existing FITS.

    We recompute the moon and zodi amplitude weights on the shared wave grid
    using a lightweight `SkyDecomp` build to avoid smuggling the whole
    LSF-iterative machinery through the calibration step; the amplitude
    weights only depend on the basis + wavelength grid, both fully defined
    by that class.
    """
    model = SkyDecomp(
        wave=wave, lsf_sigma=1.0, base_dir=str(data_root),
        palace_oh_suffix='_joint_v2_updated',
        palace_diffuse_suffix='_joint_native_adam_invsky_p2_10000iter',
        split_zodi=True, n_spline_knots=25, n_zodi_spline_knots=3,
    )
    with fits.open(corpus_path, memmap=False) as h:
        coef = h["COEF"].data
        names = coef.dtype.names
        moon_names = [n for n in names if n.startswith("Moon_bs")]
        zodi_names = [n for n in names if n.startswith("Zodi_bs")]
        moon_coef = np.stack([coef[n] for n in moon_names], axis=1).astype(np.float64)
        zodi_coef = np.stack([coef[n] for n in zodi_names], axis=1).astype(np.float64)

    w_moon = np.asarray(model.moon_amp_weights, dtype=np.float64)
    w_zodi = np.asarray(model.zodi_amp_weights, dtype=np.float64)
    if moon_coef.shape[1] != w_moon.size:
        raise ValueError(
            f"Moon coefficient count {moon_coef.shape[1]} does not match "
            f"basis size {w_moon.size}; use --n-spline-knots consistent with "
            "the corpus."
        )
    if zodi_coef.shape[1] != w_zodi.size:
        raise ValueError(
            f"Zodi coefficient count {zodi_coef.shape[1]} does not match "
            f"basis size {w_zodi.size}; use --n-zodi-spline-knots consistent."
        )
    moon_amp = moon_coef @ w_moon
    zodi_amp = zodi_coef @ w_zodi
    return moon_amp, zodi_amp


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_fits", type=Path, help="Spectra FITS with META HDU (source of ra/dec/mjd).")
    p.add_argument("corpus_dir", type=Path, help="Directory holding the *_decomp_{sci,sky1,sky2}{suffix}.fits.")
    p.add_argument("data_root", type=Path, help="sky_decomp/data directory (Leinert, ephemeris, PMD).")
    p.add_argument("--suffix", default="_lsf_surface_iterative_split_zodi",
                   help="Suffix that identifies the decomp FITS files under corpus_dir.")
    p.add_argument("--arm", default="sci", choices=list(ARM_ROLES),
                   help="Which arm's decomposition to use for the data-side amplitude side of the fit.")
    p.add_argument("--output", type=Path, required=True, help="Where to write the calibration JSON.")
    p.add_argument("--moon-smooth-gate-k", type=float, default=100.0,
                   help="Gate strength for moon smooth_lambda (0 disables gating).")
    p.add_argument("--zodi-smooth-gate-k", type=float, default=0.0,
                   help="Gate strength for zodi smooth_lambda (0 disables; zodi is already heavily smoothed).")
    p.add_argument("--moon-amp-prior-lambda", type=float, default=1.0e-3)
    p.add_argument("--zodi-amp-prior-lambda", type=float, default=1.0e-3)
    p.add_argument("--moon-smooth-lambda-base", type=float, default=0.1)
    p.add_argument("--zodi-smooth-lambda-base", type=float, default=0.1)
    p.add_argument(
        "--zodi-calib-max-moon-proxy",
        type=float,
        default=1.0,
        help=(
            "Exclude rows with KS moon proxy above this value from the zodi "
            "calibration.  On bright-moon rows the pre-prior fits absorb some "
            "moon signal into Zodi_bs, inflating the ratio "
            "sum_i w_i c_zodi_i / Leinert_b500.  Restricting the zodi fit to "
            "moon-quiet rows (default: KS < 1) recovers the true physics ratio."
        ),
    )
    args = p.parse_args()

    input_fits: Path = args.input_fits.resolve()
    corpus_dir: Path = args.corpus_dir.resolve()
    data_root: Path = args.data_root.resolve()

    if not input_fits.is_file():
        raise SystemExit(f"input_fits does not exist: {input_fits}")
    if not corpus_dir.is_dir():
        raise SystemExit(f"corpus_dir does not exist: {corpus_dir}")
    if not (data_root / "moon_zodi").is_dir():
        raise SystemExit(
            f"data_root does not contain moon_zodi/: {data_root}"
        )

    corpus_stem = input_fits.stem  # e.g. lvmsframe_median_stack_..._every10
    decomp_path = corpus_dir / f"{corpus_stem}_decomp_{args.arm}{args.suffix}.fits"
    if not decomp_path.is_file():
        raise SystemExit(f"Decomp product missing: {decomp_path}")

    with fits.open(input_fits, memmap=False) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        meta = hdul["META"].data
    print(f"loaded meta: {len(meta)} rows; wave range [{wave.min():.1f}, {wave.max():.1f}]")

    observations = _make_observations(meta, args.arm)
    print(f"computing geometry proxies for {len(observations)} rows ...")
    proxies = compute_geometry_proxies(observations, data_dir=data_root)
    print(f"  valid rows: {int(proxies['valid'].sum())}/{len(observations)}")

    print(f"computing data-side integrated amplitudes from {decomp_path.name} ...")
    moon_data_amp, zodi_data_amp = _data_amps(decomp_path, wave, data_root)

    # Zodi calibration mask: exclude rows where the KS moon proxy indicates
    # a bright moon, so pre-prior moon->zodi cross-talk does not bias k_zodi.
    zodi_valid = proxies["valid"].copy()
    if args.zodi_calib_max_moon_proxy > 0.0:
        moon_low = np.isfinite(proxies["moon_amp_proxy"]) & (
            proxies["moon_amp_proxy"] < args.zodi_calib_max_moon_proxy
        )
        zodi_valid = zodi_valid & moon_low
        print(
            f"  zodi calibration restricted to {int(zodi_valid.sum())} moon-quiet "
            f"rows (KS proxy < {args.zodi_calib_max_moon_proxy})"
        )

    cal = fit_calibration(
        moon_proxy=proxies["moon_amp_proxy"],
        zodi_proxy=proxies["zodi_b500"],
        moon_data_amp=moon_data_amp,
        zodi_data_amp=zodi_data_amp,
        valid=proxies["valid"],
        moon_smooth_gate_k=args.moon_smooth_gate_k,
        zodi_smooth_gate_k=args.zodi_smooth_gate_k,
        moon_amp_prior_lambda=args.moon_amp_prior_lambda,
        zodi_amp_prior_lambda=args.zodi_amp_prior_lambda,
        moon_smooth_lambda_base=args.moon_smooth_lambda_base,
        zodi_smooth_lambda_base=args.zodi_smooth_lambda_base,
    )

    # Recompute k_zodi and zodi_ref_amp restricted to moon-quiet rows.
    if args.zodi_calib_max_moon_proxy > 0.0 and int(zodi_valid.sum()) > 32:
        from sky_decomp.moon_zodi_priors import robust_linear_slope
        import dataclasses
        zodi_proxy = np.asarray(proxies["zodi_b500"], dtype=np.float64)
        zodi_data = np.asarray(zodi_data_amp, dtype=np.float64)
        m = zodi_valid & np.isfinite(zodi_proxy) & np.isfinite(zodi_data)
        k_zodi_clean = robust_linear_slope(zodi_proxy[m], zodi_data[m])
        pos = zodi_proxy[m]
        pos = pos[pos > 0.0]
        zodi_ref_clean = (
            float(np.percentile(pos, 50.0)) if pos.size else cal.zodi_ref_amp
        )
        print(
            f"  moon-quiet zodi calibration: k_zodi {cal.k_zodi:.3f} -> "
            f"{k_zodi_clean:.3f}, zodi_ref_amp {cal.zodi_ref_amp:.3f} -> "
            f"{zodi_ref_clean:.3f} (using {int(m.sum())} rows)"
        )
        cal = dataclasses.replace(cal, k_zodi=float(k_zodi_clean), zodi_ref_amp=zodi_ref_clean)

    print(f"calibration:")
    for k, v in cal.as_dict().items():
        print(f"  {k}: {v}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cal.to_json(args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
