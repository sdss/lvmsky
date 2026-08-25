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
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import multiprocessing as mp

import numpy as np
from astropy.io import fits
from tqdm import tqdm

try:
    from sky_decomp.result_io import results_to_fits
except ModuleNotFoundError:
    from skysub.sky_decomp.result_io import results_to_fits

try:
    from sky_decomp.moon_zodi_model import (
        DEFAULT_DATA_ROOT as DEFAULT_MOON_ZODI_DATA_ROOT,
        validate_decomposition_data_root,
    )
except ModuleNotFoundError:
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
SPLIT_ZODI_N_KNOTS_DEFAULT = 3
SPLIT_ZODI_SMOOTH_LAMBDA_DEFAULT = 1.0e-1
SPLIT_ZODI_MOON_ALBEDO_PHASE_DEG = 30.0
SPLIT_ZODI_COLOR_EXPONENT = 0.26
# Moon_bs interior-knot count; matches SkyDecomp.__init__ default in fit.py.
MOON_N_KNOTS_DEFAULT = 25


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
    moon_zodi_priors_bundle=None,
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
        _WORKER_PRIORS_BUNDLE

    _WORKER_PRIORS_BUNDLE = moon_zodi_priors_bundle

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
    if fit_model == MOON_ZODI_FIT_MODEL:
        _WORKER_LSF = {
            "sci": np.asarray(_WORKER_HDU["LSF_SCI"].data),
            "sky1": np.asarray(_WORKER_HDU["LSF_SKY_NEAR"].data),
            "sky2": np.asarray(_WORKER_HDU["LSF_SKY_FAR"].data),
        }
        _WORKER_META = _WORKER_HDU["META"].data
    elif fit_model == SPLIT_ZODI_FIT_MODEL and moon_zodi_priors_bundle is not None:
        _WORKER_LSF = {}
        _WORKER_META = _WORKER_HDU["META"].data
    else:
        _WORKER_LSF = {}
        _WORKER_META = None
    if fit_model == "baseline":
        from sky_decomp.fit import SkyDecomp

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
        from sky_decomp.lsf_surface_iterative import (
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
        from sky_decomp.lsf_surface_iterative import (
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
            config=LSFSurfaceIterativeConfig(
                n_refinement_cycles=n_refinement_cycles,
            ),
        )
    elif fit_model == MOON_ZODI_FIT_MODEL:
        from sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
        from sky_decomp.moon_zodi_lsf_surface_iterative import (
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
    from sky_decomp.moon_zodi_model import MoonZodiObservation

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


def _split_zodi_prior_kwargs_for_row(kind: str, row_idx: int) -> dict:
    """Look up per-row split_zodi prior overrides from ``_WORKER_PRIORS_BUNDLE``.

    Returns an empty dict when no priors bundle is loaded, so the fit call
    reduces to the legacy signature. When a bundle is loaded and this row's
    geometry has not been computed yet, computes it lazily inside the worker
    (parallelized across workers) and caches the result. Applies the
    calibrated slopes to produce the actual `fit` kwargs.
    """
    bundle = _WORKER_PRIORS_BUNDLE
    if bundle is None:
        return {}
    cache = bundle["arm_cache"].get(kind)
    if cache is None:
        return {}

    if not cache["computed"][row_idx]:
        # Lazy per-row geometry compute.
        obs = _moon_zodi_observation(kind, row_idx)
        from sky_decomp.moon_zodi_priors import (
            compute_geometry_proxies,
            resolve_row_targets,
        )
        proxies = compute_geometry_proxies([obs], data_dir=bundle["data_root"])
        moon_prox = float(proxies["moon_amp_proxy"][0])
        zodi_prox = float(proxies["zodi_b500"][0])
        targets = resolve_row_targets(
            moon_proxy=moon_prox,
            zodi_proxy=zodi_prox,
            calibration=bundle["calibration"],
        )
        cache["moon_amp_prior"][row_idx] = (
            targets.moon_amp_prior if targets.moon_amp_prior is not None else np.nan
        )
        cache["zodi_amp_prior"][row_idx] = (
            targets.zodi_amp_prior if targets.zodi_amp_prior is not None else np.nan
        )
        cache["moon_smooth_lambda_override"][row_idx] = (
            targets.moon_smooth_lambda_override
            if targets.moon_smooth_lambda_override is not None
            else np.nan
        )
        cache["zodi_smooth_lambda_override"][row_idx] = (
            targets.zodi_smooth_lambda_override
            if targets.zodi_smooth_lambda_override is not None
            else np.nan
        )
        cache["computed"][row_idx] = True

    kw = {}
    ma = cache["moon_amp_prior"][row_idx]
    if np.isfinite(ma):
        kw["moon_amp_prior"] = float(ma)
        kw["moon_amp_prior_lambda"] = float(bundle["moon_amp_prior_lambda"])
    za = cache["zodi_amp_prior"][row_idx]
    if np.isfinite(za):
        kw["zodi_amp_prior"] = float(za)
        kw["zodi_amp_prior_lambda"] = float(bundle["zodi_amp_prior_lambda"])
    msl = cache["moon_smooth_lambda_override"][row_idx]
    if np.isfinite(msl):
        kw["moon_smooth_lambda_override"] = float(msl)
    zsl = cache["zodi_smooth_lambda_override"][row_idx]
    if np.isfinite(zsl):
        kw["zodi_smooth_lambda_override"] = float(zsl)
    return kw


def fit_chunk_worker(args):
    """Fit one chunk of spectra using the worker-local SkyDecomp instance."""
    global _WORKER_DECOMPOSER, _WORKER_FACTOR, _WORKER_PROGRESS_QUEUE, _WORKER_FIT_MODEL
    if _WORKER_DECOMPOSER is None:
        raise RuntimeError("Worker SkyDecomp has not been initialised.")
    kind, idx0, idx1 = args
    flux_chunk = np.asarray(_WORKER_FLUX[kind][idx0:idx1], dtype=np.float64)
    out = []
    for j in range(flux_chunk.shape[0]):
        idx = idx0 + j
        flux_row = flux_chunk[j] * _WORKER_FACTOR
        ivar_row = np.ones_like(flux_row)
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
            prior_kwargs = _split_zodi_prior_kwargs_for_row(kind, idx)
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                verbose=False,
                **prior_kwargs,
            )
        elif _WORKER_FIT_MODEL == "moon-zodi-lsf-surface-iterative":
            # Preserve invalid source pixels; zero IVAR excludes them without
            # interpolating, imputing, cropping, or changing the native grid.
            ivar_row = np.isfinite(flux_row).astype(np.float64)
            lsf_row = np.asarray(_WORKER_LSF[kind][idx], dtype=np.float64)
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

    if palace_dir is None:
        raise ValueError(
            "palace_dir is required for baseline and lsf-surface-iterative fits"
        )
    return resolve_base_dir(palace_dir), None


def _iter_chunk_tasks(n_rows, chunk_size):
    for kind in ("sci", "sky1", "sky2"):
        for i0 in range(0, n_rows, chunk_size):
            i1 = min(i0 + chunk_size, n_rows)
            # Sending only index ranges keeps per-task IPC tiny.
            yield (kind, i0, i1)


def _load_moon_zodi_priors_bundle(
    *,
    calibration_path: Path,
    data_file: str | Path,
    n_rows: int,
    moon_zodi_data_root: Path | None,
) -> dict:
    """Package a moon/zodi priors bundle for worker consumption.

    Loads the calibration JSON on the main process; the actual per-row
    geometry proxies (KS moon proxy + Leinert B500) are computed lazily
    inside each worker on the rows it processes, which parallelises the
    astropy SkyCoord + ephemeris calls that dominate the precompute cost.
    """
    from sky_decomp.moon_zodi_priors import MoonZodiPriorCalibration
    calibration = MoonZodiPriorCalibration.from_json(calibration_path)
    print(f"  moon_zodi_priors: loaded {calibration_path}")
    print(f"    k_moon={calibration.k_moon:.4g}, k_zodi={calibration.k_zodi:.4g}, "
          f"moon_ref_amp={calibration.moon_ref_amp:.3g}, zodi_ref_amp={calibration.zodi_ref_amp:.3g}, "
          f"moon_smooth_gate_k={calibration.moon_smooth_gate_k}, "
          f"zodi_smooth_gate_k={calibration.zodi_smooth_gate_k}")
    if moon_zodi_data_root is None:
        raise ValueError(
            "moon_zodi_priors requires a resolvable moon_zodi data root; "
            "pass --moon-zodi-data-root or accept the default palace layout."
        )
    arm_cache = {}
    for kind in ("sci", "sky1", "sky2"):
        arm_cache[kind] = {
            "computed": np.zeros(n_rows, dtype=bool),
            "moon_amp_prior": np.full(n_rows, np.nan, dtype=np.float64),
            "zodi_amp_prior": np.full(n_rows, np.nan, dtype=np.float64),
            "moon_smooth_lambda_override": np.full(n_rows, np.nan, dtype=np.float64),
            "zodi_smooth_lambda_override": np.full(n_rows, np.nan, dtype=np.float64),
        }
    return {
        "calibration": calibration,
        "data_root": str(moon_zodi_data_root),
        "moon_amp_prior_lambda": float(calibration.moon_amp_prior_lambda),
        "zodi_amp_prior_lambda": float(calibration.zodi_amp_prior_lambda),
        "arm_cache": arm_cache,
    }


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
    moon_zodi_priors_path=None,
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
    if fit_model == MOON_ZODI_FIT_MODEL:
        print(f"  moon_zodi_data_root={str(resolved_moon_zodi_data_root)!r}")
    if fit_model in ("lsf-surface-iterative", SPLIT_ZODI_FIT_MODEL):
        print(f"  n_spline_knots={n_spline_knots} "
              f"(Moon_bs basis = {int(n_spline_knots) + 4})")
    if fit_model == SPLIT_ZODI_FIT_MODEL:
        print(f"  n_zodi_spline_knots={n_zodi_spline_knots}, "
              f"zodi_smooth_lambda={zodi_smooth_lambda}")
    print(f"  pin_workers={pin_workers}, diagnose_threads={diagnose_threads}")
    print(f"  base_dir={base_dir}")

    moon_zodi_priors_bundle = None
    if moon_zodi_priors_path is not None:
        if fit_model != SPLIT_ZODI_FIT_MODEL:
            raise ValueError(
                f"--moon-zodi-priors is only supported when "
                f"--fit-model={SPLIT_ZODI_FIT_MODEL}."
            )
        priors_data_root = resolved_moon_zodi_data_root
        if priors_data_root is None and moon_zodi_data_root is not None:
            candidate = Path(moon_zodi_data_root).expanduser().resolve()
            if (candidate / "moon_zodi").is_dir():
                priors_data_root = candidate
            elif candidate.is_dir() and (candidate.name == "moon_zodi"):
                priors_data_root = candidate
        if priors_data_root is None:
            candidate = base_dir / "moon_zodi"
            if not candidate.is_dir():
                raise ValueError(
                    f"Could not locate moon_zodi assets for --moon-zodi-priors; "
                    f"tried {candidate}. Pass --moon-zodi-data-root explicitly."
                )
            priors_data_root = candidate
        moon_zodi_priors_bundle = _load_moon_zodi_priors_bundle(
            calibration_path=moon_zodi_priors_path,
            data_file=data_file,
            n_rows=n_rows,
            moon_zodi_data_root=priors_data_root,
        )

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
            moon_zodi_priors_bundle,
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
            "Legacy project/PALACE root. Required by baseline and "
            "lsf-surface-iterative; optional bundle root for the Moon/Zodi mode."
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
            "pmd_refcont_adam_v1.dat (default: canonical unsuffixed files)."
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
            "Moon/Zodi method (default: packaged skysub/sky_decomp/data)."
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
        "--moon-zodi-priors",
        default=None,
        help=(
            "Path to a calibration JSON produced by "
            "`skysub/calibrate_moon_zodi_priors.py`. When set (and "
            f"--fit-model={SPLIT_ZODI_FIT_MODEL}), each row's fit is nudged toward "
            "geometry-derived moon and zodi amplitude targets and the moon "
            "spline curvature penalty is gated by the KS moon-brightness "
            "proxy so it goes flat when the moon should be near zero."
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
        and args.fit_model != MOON_ZODI_FIT_MODEL
        and args.palace_dir is None
    ):
        parser.error(
            "palace_dir is required for baseline and lsf-surface-iterative fits"
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
            moon_zodi_priors_path=args.moon_zodi_priors,
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
