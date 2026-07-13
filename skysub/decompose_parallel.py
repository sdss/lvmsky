"""
Run sky spectral decomposition on a median-stacked LVM frame.

Usage:
    python run_sky_decomp.py <data_file> <palace_dir> [--n-workers N] [--lsf-sigma S] [--factor F] [--output-dir DIR]

Example:
    python run_sky_decomp.py lvmsframe_median_stack.fits ../ --n-workers 8
"""

import argparse
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm


def fit_chunk_worker(args):
    """Worker function: builds its own SkyDecomp and fits a chunk of spectra."""
    chunk_idxs, flx_sci_chunk, flx_sky1_chunk, flx_sky2_chunk, flx_ivar_chunk, wave, lsf_sigma, base_dir = args

    from sky_decomp.fit import SkyDecomp
    local_decomposer = SkyDecomp(wave, lsf_sigma=lsf_sigma, base_dir=base_dir)

    local_sci, local_sky1, local_sky2 = [], [], []
    for i in range(len(chunk_idxs)):
        local_sci.append(local_decomposer.fit(flx_sci_chunk[i],  flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))
        local_sky1.append(local_decomposer.fit(flx_sky1_chunk[i], flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))
        local_sky2.append(local_decomposer.fit(flx_sky2_chunk[i], flx_ivar_chunk[i], verbose=False, n_lsf_refits=3))

    return chunk_idxs, local_sci, local_sky1, local_sky2


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
    coef_hdu = fits.ImageHDU(coef_arr, name="COEF")
    design_names = results[0].design_names
    for i, name in enumerate(design_names):
        coef_hdu.header[f"COEF{i:04d}"] = name

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


def run(data_file, palace_dir, n_workers, lsf_sigma, factor, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {data_file} ...")
    wave     = fits.getdata(data_file, "WAVE").astype(np.float64)
    flx_sky1 = fits.getdata(data_file, "FLUX_SKY_NEAR").astype(np.float64) * factor
    flx_sky2 = fits.getdata(data_file, "FLUX_SKY_FAR").astype(np.float64) * factor
    flx_sci  = fits.getdata(data_file, "FLUX_SCI").astype(np.float64) * factor
    flx_ivar = np.ones_like(flx_sci)

    print(f"  {len(flx_sci)} spectra, {len(wave)} wavelength pixels")
    print(f"  n_workers={n_workers}, lsf_sigma={lsf_sigma}, factor={factor}")

    idxs   = list(range(len(flx_sci)))
    chunks = np.array_split(idxs, n_workers)

    worker_args = [
        (
            list(chunk),
            flx_sci[chunk],
            flx_sky1[chunk],
            flx_sky2[chunk],
            flx_ivar[chunk],
            wave,
            lsf_sigma,
            palace_dir,
        )
        for chunk in chunks
    ]

    result_sci  = [None] * len(idxs)
    result_sky1 = [None] * len(idxs)
    result_sky2 = [None] * len(idxs)

    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(fit_chunk_worker, args): tid
                   for tid, args in enumerate(worker_args)}
        for future in tqdm(as_completed(futures), total=n_workers, desc="Chunks done"):
            chunk_idxs, sci, sky1, sky2 = future.result()
            for i, idx in enumerate(chunk_idxs):
                result_sci[idx]  = sci[i]
                result_sky1[idx] = sky1[i]
                result_sky2[idx] = sky2[i]

    elapsed = time.perf_counter() - t0
    print(f"Fitting done in {elapsed:.1f}s ({elapsed/len(idxs):.2f}s per spectrum)")

    stem = Path(data_file).stem
    results_to_fits(result_sci,  output_dir / f"{stem}_decomp_sci.fits")
    results_to_fits(result_sky1, output_dir / f"{stem}_decomp_sky1.fits")
    results_to_fits(result_sky2, output_dir / f"{stem}_decomp_sky2.fits")


def main():
    parser = argparse.ArgumentParser(description="LVM sky spectral decomposition")
    parser.add_argument("data_file",    help="Input FITS file (median stacked LVM frame)")
    parser.add_argument("palace_dir",   help="Path to PALACE data directory")
    parser.add_argument("--n-workers",  type=int,   default=4,    help="Number of parallel worker processes (default: 4)")
    parser.add_argument("--lsf-sigma",  type=float, default=0.5,  help="LSF Gaussian sigma in Å (default: 0.5)")
    parser.add_argument("--factor",     type=float, default=1e14, help="Flux scaling factor (default: 1e14)")
    parser.add_argument("--output-dir", default=".",              help="Output directory for result FITS files (default: .)")
    args = parser.parse_args()

    run(
        data_file  = args.data_file,
        palace_dir = args.palace_dir,
        n_workers  = args.n_workers,
        lsf_sigma  = args.lsf_sigma,
        factor     = args.factor,
        output_dir = args.output_dir,
    )


if __name__ == "__main__":
    main()