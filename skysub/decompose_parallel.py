#!/usr/bin/env python3

"""
Run sky spectral decomposition on a median-stacked LVM frame.

Usage:
    python decompose_parallel.py <data_file> <palace_dir> [options]

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
from astropy.table import Table
from tqdm import tqdm

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
_WORKER_PROGRESS_QUEUE = None
_WORKER_FIT_MODEL = "baseline"


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
):
    """Initialise one SkyDecomp instance per worker process."""
    global \
        _WORKER_DECOMPOSER, \
        _WORKER_FACTOR, \
        _WORKER_HDU, \
        _WORKER_FLUX, \
        _WORKER_PROGRESS_QUEUE, \
        _WORKER_FIT_MODEL

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
    # Keep worker-local memmapped access to flux tables to avoid large IPC payloads.
    _WORKER_HDU = fits.open(data_file, memmap=True)
    _WORKER_PROGRESS_QUEUE = progress_queue
    _WORKER_FLUX = {
        "sci": np.asarray(_WORKER_HDU["FLUX_SCI"].data),
        "sky1": np.asarray(_WORKER_HDU["FLUX_SKY_NEAR"].data),
        "sky2": np.asarray(_WORKER_HDU["FLUX_SKY_FAR"].data),
    }
    if fit_model == "baseline":
        from sky_decomp.fit import SkyDecomp

        _WORKER_DECOMPOSER = SkyDecomp(
            wave,
            lsf_sigma=lsf_sigma,
            base_dir=base_dir,
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
            moon_smooth_lambda=0.1,
            moon_interline_boost=0.0,
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
        else:
            result = _WORKER_DECOMPOSER.fit(
                flux_row,
                ivar_row,
                verbose=False,
            )
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


def results_to_fits(results, filename):
    """Write a list of SkyDecompResult objects to a FITS file."""
    if not results:
        raise ValueError("results must contain at least one decomposition result")
    rows = {
        "t_o2": [r.t_o2 for r in results],
        "t_o2_err": [r.t_o2_err for r in results],
        "reduced_chi2": [r.reduced_chi2 for r in results],
        "r2": [r.r2 for r in results],
        "rms_resid": [r.rms_resid for r in results],
        "resid_level": [r.resid_level for r in results],
        "fit_status": [r.fit_status for r in results],
        "fit_summary": [r.fit_summary for r in results],
        "fit_elapsed_sec": [r.fit_elapsed_sec for r in results],
        "peak_memory_mb": [r.peak_memory_mb for r in results],
        "o2_fit_status": [r.o2_fit_status for r in results],
        "o2_fit_summary": [r.o2_fit_summary for r in results],
        "o2_fit_elapsed_sec": [r.o2_fit_elapsed_sec for r in results],
        "o2_valid_frac": [r.o2_valid_frac for r in results],
    }
    t = Table(rows)

    def stack(attr):
        return np.vstack([getattr(r, attr) for r in results])

    coef_arr = stack("coef")
    design_names = results[0].design_names
    if coef_arr.shape[1] != len(design_names):
        raise ValueError(
            f"Coefficient count ({coef_arr.shape[1]}) does not match "
            f"number of design names ({len(design_names)})."
        )
    coef_table = Table({name: coef_arr[:, i] for i, name in enumerate(design_names)})
    coef_hdu = fits.BinTableHDU(coef_table, name="COEF")

    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(t, name="META"),
            coef_hdu,
            fits.ImageHDU(stack("bestfit"), name="BESTFIT"),
            fits.ImageHDU(stack("bestfit_lsf"), name="BESTFIT_LSF"),
            fits.ImageHDU(stack("resid"), name="RESID"),
        ]
    )

    comp_keys = list(results[0].components.keys())
    for key in comp_keys:
        arr = np.vstack([r.components[key] for r in results])
        hdul.append(fits.ImageHDU(arr, name=f"COMP_{key.upper()}"))

    iterative_result_types = []
    try:
        from sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeResult

        iterative_result_types.append(LSFSurfaceIterativeResult)
    except ModuleNotFoundError:
        pass
    try:
        from skysub.sky_decomp.lsf_surface_iterative import (
            LSFSurfaceIterativeResult,
        )

        iterative_result_types.append(LSFSurfaceIterativeResult)
    except ModuleNotFoundError:
        pass
    iterative_result_types = tuple(set(iterative_result_types))
    is_iterative = [
        bool(iterative_result_types) and isinstance(result, iterative_result_types)
        for result in results
    ]
    if any(is_iterative) and not all(is_iterative):
        raise ValueError("Cannot mix baseline and LSF-surface iterative results")
    if all(is_iterative):
        states = [result.lsf_state for result in results]
        if any(state is None for state in states):
            raise ValueError("Every iterative result must contain a compact LSF state")
        try:
            from sky_decomp._lsf_surface_fits import build_lsf_hdus
        except ModuleNotFoundError:
            from skysub.sky_decomp._lsf_surface_fits import build_lsf_hdus

        hdul.extend(build_lsf_hdus(states))

    hdul.writeto(filename, overwrite=True)
    print(
        f"Wrote {len(results)} results, {coef_arr.shape[1]} coefs, "
        f"{len(comp_keys)} components → {filename}"
    )


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
):
    base_dir = resolve_base_dir(palace_dir)
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
        ),
    ) as executor:
        pbar = tqdm(
            total=3 * n_rows,
            desc="Arm-fits completed",
            unit=" arm-fits",
            mininterval=0.2,
            position=0,
            leave=True,
        )
        pbar.set_postfix(chunks_done=f"0/{n_tasks}")
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
                pbar.set_postfix(chunks_done=f"{completed}/{n_tasks}")
            _submit_until_full()
        _drain_progress_queue()
        pbar.close()

    progress_queue.close()
    progress_queue.join_thread()

    elapsed = time.perf_counter() - t0
    print(f"Fitting done in {elapsed:.1f}s ({elapsed / n_rows:.2f}s per spectrum)")

    stem = Path(data_file).stem
    suffix = "" if fit_model == "baseline" else "_lsf_surface_iterative"
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
            variant = (
                "_lsf_surface_iterative"
                if "lsf_surface_iterative" in Path(decomp_path).stem.lower()
                else ""
            )
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
            lsf_extensions = ("LSF_COEF", "LSF_KNOTS", "LSF_META")
            present = [name in hdul_dec for name in lsf_extensions]
            if any(present) and not all(present):
                raise KeyError(f"Incomplete LSF-surface extensions in {decomp_path}")
            if all(present):
                compact_hdus.extend(
                    _copy_hdu_with_name(hdul_dec[name], name) for name in lsf_extensions
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
    """
    if n < 1:
        raise ValueError("n must be >= 1")

    with fits.open(input_path) as hdul:
        if row_hdu_name not in hdul:
            raise KeyError(f"HDU '{row_hdu_name}' not found in {input_path}")

        n_rows = len(hdul[row_hdu_name].data)
        keep = slice(None, None, n)

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
                if data is not None and len(data) == n_rows:
                    data = data[keep]
                out_hdus.append(type(hdu)(data=data, header=header, name=hdu.name))

            elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                data = hdu.data
                if data is not None and getattr(data, "ndim", 0) >= 1 and data.shape[0] == n_rows:
                    data = data[keep, ...]
                out_hdus.append(type(hdu)(data=data, header=header, name=hdu.name))

            else:
                out_hdus.append(hdu.copy())

        fits.HDUList(out_hdus).writeto(output_path, overwrite=True)


def main():
    parser = argparse.ArgumentParser(description="LVM sky spectral decomposition")
    parser.add_argument("data_file", help="Input FITS file (median stacked LVM frame)")
    parser.add_argument(
        "palace_dir", help="Path to the project base directory or directly to the palace directory"
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
        choices=("baseline", "lsf-surface-iterative"),
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
    )

    suffix = "" if args.fit_model == "baseline" else "_lsf_surface_iterative"
    extract_meta_and_coef_products(
        input_fits_path=args.data_file,
        decomp_fits_path_1=Path(args.output_dir)
        / f"{Path(args.data_file).stem}_decomp_sky1{suffix}.fits",
        decomp_fits_path_2=Path(args.output_dir)
        / f"{Path(args.data_file).stem}_decomp_sky2{suffix}.fits",
        decomp_fits_path_3=Path(args.output_dir)
        / f"{Path(args.data_file).stem}_decomp_sci{suffix}.fits",
    )

    thin_fits_every_n(f"{Path(args.data_file).stem}_decomp_sci{suffix}.fits", f"{Path(args.data_file).stem}_every10_decomp_sci{suffix}.fits", 10)
    thin_fits_every_n(f"{Path(args.data_file).stem}_decomp_sky1{suffix}.fits", f"{Path(args.data_file).stem}_every10_decomp_sky1{suffix}.fits", 10)
    thin_fits_every_n(f"{Path(args.data_file).stem}_decomp_sky2{suffix}.fits", f"{Path(args.data_file).stem}_every10_decomp_sky2{suffix}.fits", 10)

if __name__ == "__main__":
    main()
