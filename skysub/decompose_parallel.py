#!/usr/bin/env python3

"""
Run sky spectral decomposition on a median-stacked LVM frame.

Usage:
    python decompose_parallel.py <data_file> <palace_dir> [--n-workers N] [--lsf-sigma S] [--factor F] [--chunk-size N] [--output-dir DIR]

Example:
    python decompose_parallel.py lvmsframe_median_stack.fits ../ --n-workers 8
"""

import argparse
import platform
import sys
import time
import queue as queue_mod
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
import multiprocessing as mp

import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm


_WORKER_DECOMPOSER = None
_WORKER_FACTOR = 1.0
_WORKER_HDU = None
_WORKER_FLUX = {}
_WORKER_PROGRESS_QUEUE = None


def init_worker(wave, lsf_sigma, base_dir, factor, data_file, progress_queue=None):
    """Initialise one SkyDecomp instance per worker process."""
    global _WORKER_DECOMPOSER, _WORKER_FACTOR, _WORKER_HDU, _WORKER_FLUX, _WORKER_PROGRESS_QUEUE

    # pid = mp.current_process().pid
    # print(f"[worker {pid}] initializing SkyDecomp", file=sys.stderr, flush=True)

    from sky_decomp.fit import SkyDecomp
    _WORKER_FACTOR = float(factor)
    # Keep worker-local memmapped access to flux tables to avoid large IPC payloads.
    _WORKER_HDU = fits.open(data_file, memmap=True)
    _WORKER_PROGRESS_QUEUE = progress_queue
    _WORKER_FLUX = {
        "sci": np.asarray(_WORKER_HDU["FLUX_SCI"].data),
        "sky1": np.asarray(_WORKER_HDU["FLUX_SKY_NEAR"].data),
        "sky2": np.asarray(_WORKER_HDU["FLUX_SKY_FAR"].data),
    }
    _WORKER_DECOMPOSER = SkyDecomp(wave, lsf_sigma=lsf_sigma, base_dir=base_dir,\
                                    moon_smooth_lambda=0.1,\
                                    moon_interline_boost=10000.0,\
                                    moon_interline_red_min=6000.0,\
                                    moon_interline_exclusion_a=2.5,\
                                    moon_interline_line_flux_threshold=0.01)
    # print(f"[worker {pid}] SkyDecomp ready", file=sys.stderr, flush=True)


def fit_chunk_worker(args):
    """Fit one chunk of spectra using the worker-local SkyDecomp instance."""
    global _WORKER_DECOMPOSER, _WORKER_FACTOR, _WORKER_PROGRESS_QUEUE
    if _WORKER_DECOMPOSER is None:
        raise RuntimeError("Worker SkyDecomp has not been initialised.")

    kind, idx0, idx1 = args
    flux_chunk = np.asarray(_WORKER_FLUX[kind][idx0:idx1], dtype=np.float64)
    out = []
    for j in range(flux_chunk.shape[0]):
        idx = idx0 + j
        flux_row = flux_chunk[j] * _WORKER_FACTOR
        ivar_row = np.ones_like(flux_row)
        result = _WORKER_DECOMPOSER.fit(flux_row, ivar_row, verbose=False, n_lsf_refits=3)
        out.append((idx, result))
        if _WORKER_PROGRESS_QUEUE is not None:
            _WORKER_PROGRESS_QUEUE.put(1)
    return kind, out


def worker_ping(_=None):
    """Return a worker pid once initialization has completed."""
    if _WORKER_DECOMPOSER is None:
        raise RuntimeError("Worker SkyDecomp has not been initialised.")
    return mp.current_process().pid


def _get_mp_context():
    if platform.system() == "Darwin":
        return mp.get_context("spawn")
    return mp.get_context()


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
    rows = {
        "t_o2":               [r.t_o2 for r in results],
        "t_o2_err":           [r.t_o2_err for r in results],
        "reduced_chi2":       [r.reduced_chi2 for r in results],
        "r2":                 [r.r2 for r in results],
        "rms_resid":          [r.rms_resid for r in results],
        "resid_level":        [r.resid_level for r in results],
        "fit_status":         [r.fit_status for r in results],
        "fit_summary":        [r.fit_summary for r in results],
        "fit_elapsed_sec":    [r.fit_elapsed_sec for r in results],
        "peak_memory_mb":     [r.peak_memory_mb for r in results],
        "o2_fit_status":      [r.o2_fit_status for r in results],
        "o2_fit_summary":     [r.o2_fit_summary for r in results],
        "o2_fit_elapsed_sec": [r.o2_fit_elapsed_sec for r in results],
        "o2_valid_frac":      [r.o2_valid_frac for r in results],
    }
    t = Table(rows)

    def stack(attr):
        return np.vstack([getattr(r, attr) for r in results])

    coef_arr = stack("coef")
    design_names = results[0].design_names
    if coef_arr.shape[1] != len(design_names):
        raise ValueError(
            f"Coefficient count ({coef_arr.shape[1]}) does not match number of design names ({len(design_names)})."
        )
    coef_table = Table({name: coef_arr[:, i] for i, name in enumerate(design_names)})
    coef_hdu = fits.BinTableHDU(coef_table, name="COEF")

    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(t, name="META"),
        coef_hdu,
        fits.ImageHDU(stack("bestfit"),     name="BESTFIT"),
        fits.ImageHDU(stack("bestfit_lsf"), name="BESTFIT_LSF"),
        fits.ImageHDU(stack("resid"),       name="RESID"),
    ])

    comp_keys = list(results[0].components.keys())
    for key in comp_keys:
        arr = np.vstack([r.components[key] for r in results])
        hdul.append(fits.ImageHDU(arr, name=f"COMP_{key.upper()}"))

    hdul.writeto(filename, overwrite=True)
    print(f"Wrote {len(results)} results, {coef_arr.shape[1]} coefs, {len(comp_keys)} components → {filename}")


def _iter_chunk_tasks(n_rows, chunk_size):
    for kind in ("sci", "sky1", "sky2"):
        for i0 in range(0, n_rows, chunk_size):
            i1 = min(i0 + chunk_size, n_rows)
            # Sending only index ranges keeps per-task IPC tiny.
            yield (kind, i0, i1)


def run(data_file, palace_dir, n_workers, lsf_sigma, factor, output_dir, chunk_size, max_in_flight):
    base_dir = resolve_base_dir(palace_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {data_file} ...")
    wave     = fits.getdata(data_file, "WAVE").astype(np.float64)
    with fits.open(data_file, memmap=True) as hdul:
        n_rows = int(hdul["FLUX_SCI"].data.shape[0])

    print(f"  {n_rows} spectra, {len(wave)} wavelength pixels")
    print(f"  n_workers={n_workers}, lsf_sigma={lsf_sigma}, factor={factor}")
    print(f"  chunk_size={chunk_size}, max_in_flight={max_in_flight}")
    print(f"  base_dir={base_dir}")

    idxs = list(range(n_rows))
    n_tasks = int(np.ceil(n_rows / chunk_size)) * 3

    result_sci  = [None] * len(idxs)
    result_sky1 = [None] * len(idxs)
    result_sky2 = [None] * len(idxs)
    completed = 0

    t0 = time.perf_counter()

    mp_context = _get_mp_context()
    progress_queue = mp_context.Queue()

    def _drain_progress_queue():
        inc = 0
        while True:
            try:
                inc += progress_queue.get_nowait()
            except queue_mod.Empty:
                break
        if inc:
            pbar.update(inc)

    with ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=mp_context,
        initializer=init_worker,
        initargs=(wave, lsf_sigma, str(base_dir), float(factor), str(data_file), progress_queue),
    ) as executor:
        warmup_pbar = tqdm(total=n_workers, desc="Worker pings", mininterval=0.1, position=0, leave=False)
        seen_pids = set()
        ping_futures = [executor.submit(worker_ping) for _ in range(n_workers)]
        for future in as_completed(ping_futures):
            pid = future.result()
            seen_pids.add(pid)
            warmup_pbar.update(1)
            warmup_pbar.set_postfix(unique=f"{len(seen_pids)}/{n_workers}")
        warmup_pbar.close()

        pbar = tqdm(total=3 * n_rows, desc="Fits completed", mininterval=0.2, position=0, leave=True)
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
                    if kind == "sci":
                        result_sci[idx] = result
                    elif kind == "sky1":
                        result_sky1[idx] = result
                    else:
                        result_sky2[idx] = result
                completed += 1
                pbar.set_postfix(chunks_done=f"{completed}/{n_tasks}")
            _submit_until_full()
        _drain_progress_queue()
        pbar.close()

    progress_queue.close()
    progress_queue.join_thread()

    elapsed = time.perf_counter() - t0
    print(f"Fitting done in {elapsed:.1f}s ({elapsed/len(idxs):.2f}s per spectrum)")

    stem = Path(data_file).stem
    results_to_fits(result_sci,  output_dir / f"{stem}_decomp_sci.fits")
    results_to_fits(result_sky1, output_dir / f"{stem}_decomp_sky1.fits")
    results_to_fits(result_sky2, output_dir / f"{stem}_decomp_sky2.fits")


def _copy_hdu_with_name(hdu, extname):
    """Return a copy of an HDU with a new extension name."""
    header = hdu.header.copy()
    header["EXTNAME"] = extname
    return type(hdu)(data=hdu.data, header=header, name=extname)


def _infer_decomp_label(path, index):
    """Infer a stable label for a decomposition file from its filename."""
    from pathlib import Path

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
    Each decomposition input gets its own output FITS containing just META and
    COEF. Default output paths are written in the current working directory.
    """
    from pathlib import Path

    input_path = Path(input_fits_path)
    cwd = Path.cwd()
    if meta_output_path is None:
        meta_output_path = str(cwd / f"{input_path.stem}_meta_only{input_path.suffix}")

    decomp_files = [decomp_fits_path_1, decomp_fits_path_2, decomp_fits_path_3]
    decomp_outputs = [sky1_output_path, sky2_output_path, sci_output_path]

    with fits.open(input_fits_path) as hdul_in:
        if "META" not in hdul_in:
            raise KeyError(f"Missing META extension in {input_fits_path}")
        fits.HDUList([
            fits.PrimaryHDU(),
            _copy_hdu_with_name(hdul_in["META"], "META"),
        ]).writeto(meta_output_path, overwrite=True)

    resolved_outputs = []
    for index, decomp_path in enumerate(decomp_files, start=1):
        label = _infer_decomp_label(decomp_path, index)
        out_path = decomp_outputs[index - 1]
        if out_path is None:
            out_path = str(cwd / f"{input_path.stem}_{label.lower()}_meta_coef{input_path.suffix}")
        with fits.open(decomp_path) as hdul_dec:
            for extname in ("META", "COEF"):
                if extname not in hdul_dec:
                    raise KeyError(f"Missing {extname} extension in {decomp_path}")
            fits.HDUList([
                fits.PrimaryHDU(),
                _copy_hdu_with_name(hdul_dec["META"], "META"),
                _copy_hdu_with_name(hdul_dec["COEF"], "COEF"),
            ]).writeto(out_path, overwrite=True)
        print(f"Wrote {label} META/COEF file -> {out_path}")
        resolved_outputs.append(out_path)

    print(f"Wrote META-only file -> {meta_output_path}")
    return (meta_output_path, *resolved_outputs)

def main():
    parser = argparse.ArgumentParser(description="LVM sky spectral decomposition")
    parser.add_argument("data_file",    help="Input FITS file (median stacked LVM frame)")
    parser.add_argument("palace_dir",   help="Path to the project base directory or directly to the palace directory")
    parser.add_argument("--n-workers",  type=int,   default=4,    help="Number of parallel worker processes (default: 4)")
    parser.add_argument("--lsf-sigma",  type=float, default=0.5,  help="LSF Gaussian sigma in Å (default: 0.5)")
    parser.add_argument("--factor",     type=float, default=1e14, help="Flux scaling factor (default: 1e14)")
    parser.add_argument("--chunk-size", type=int,   default=64,   help="Rows per worker task chunk (default: 64)")
    parser.add_argument("--max-in-flight", type=int, default=None, help="Max submitted chunks waiting/running at once (default: n-workers)")
    parser.add_argument("--output-dir", default=".",              help="Output directory for result FITS files (default: .)")
    args = parser.parse_args()

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")
    if args.max_in_flight is None:
        args.max_in_flight = args.n_workers
    if args.max_in_flight < 1:
        raise ValueError("--max-in-flight must be >= 1")

    run(
        data_file  = args.data_file,
        palace_dir = args.palace_dir,
        n_workers  = args.n_workers,
        lsf_sigma  = args.lsf_sigma,
        factor     = args.factor,
        output_dir = args.output_dir,
        chunk_size = args.chunk_size,
        max_in_flight = args.max_in_flight,
    )

    extract_meta_and_coef_products(
        input_fits_path=args.data_file,
        decomp_fits_path_1=Path(args.output_dir) / f"{Path(args.data_file).stem}_decomp_sky1.fits",
        decomp_fits_path_2=Path(args.output_dir) / f"{Path(args.data_file).stem}_decomp_sky2.fits",
        decomp_fits_path_3=Path(args.output_dir) / f"{Path(args.data_file).stem}_decomp_sci.fits",
    )

if __name__ == "__main__":
    main()