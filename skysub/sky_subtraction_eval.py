#! /usr/bin/env python

import os

# Maximum number of sky-subtraction subprocesses running concurrently.
N_SUBPROCESS = 8

# Number of fibers extracted from each exposure.
N_FIBERS = 10

# Number of wavelength pixels in an LVM spectrum.
N_PIX = 12401

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["NUMEXPR_MAX_THREADS"] = "1"
# clarabel (Rust QP solver) and any other Rust/Rayon library ignore OMP_NUM_THREADS.
os.environ["RAYON_NUM_THREADS"] = "1"
os.environ["POLARS_MAX_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["TBB_NUM_THREADS"] = "1"

from astropy.io import fits
from astropy.table import Table
import numpy as np
import os.path as path
import sys
import getopt
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

LVM_ROOT='/uufs/chpc.utah.edu/common/home/sdss50/'
OUTPUT_DIR=os.getenv("SAS_BASE_DIR")+"/sdsswork/users/u6045701/skyml_block"
DRPALL = os.getenv('LVM_SPECTRO_REDUX')+'/1.3.2/drpall-1.3.2.fits'

t = Table.read(DRPALL)

select = ((t['mjd']>60355) * (t['mjd']<60402) * (t['exptime']==900) * (t['tileid']>11111) * (t['fluxcal']=='MOD'))
t = t[select]
t = t[2:-1]

print(len(t))

def select_fibers(filepath, fibers=None, n=N_FIBERS, quiet_fraction=0.4):
    """Select `n` science fibers without significant continuum (i.e. no star).

    Continuum is measured as the median of the FLUX spectrum over wavelength,
    which ignores emission lines. The valid science fibers are sorted by this
    level and the `n` fibers starting at rank `quiet_fraction` (0 = lowest,
    1 = highest) are chosen, so the selection is deterministic and avoids the
    extreme fibers.

    If `fibers` is given (row indices into FLUX), no selection is done and
    those fibers are extracted instead. Use this to read the same fibers from
    the other products (ML sky, DRP sky) of the same exposure.

    Returns (fibers, flux, ivar), with flux and ivar of shape (n, npix) taken
    from the FLUX and IVAR extensions for the same rows.
    """
    with fits.open(filepath) as hdul:
        flux = hdul['FLUX'].data
        ivar = hdul['IVAR'].data

        if fibers is None:
            slitmap = Table(hdul['SLITMAP'].data)
            ok = np.asarray(slitmap['targettype'] == 'science')
            if 'fibstatus' in slitmap.colnames:
                ok &= np.asarray(slitmap['fibstatus'] == 0)

            # Fibers with no usable pixels give all-NaN rows; skip them rather
            # than call nanmedian on them (it warns, and warnings.catch_warnings
            # is not thread-safe).
            good = (ivar > 0) & np.isfinite(flux)
            ok &= good.any(axis=1)
            level = np.full(len(flux), np.nan)
            rows = np.flatnonzero(ok)
            level[rows] = np.nanmedian(np.where(good[rows], flux[rows], np.nan), axis=1)
            ok &= np.isfinite(level)

            candidates = np.flatnonzero(ok)
            if len(candidates) < n:
                raise ValueError(f"Only {len(candidates)} valid science fibers in {filepath}")

            order = candidates[np.argsort(level[candidates])]
            start = min(int(quiet_fraction * len(order)), len(order) - n)
            fibers = order[start:start + n]

        fibers = np.asarray(fibers)
        return fibers, np.array(flux[fibers], dtype=float), np.array(ivar[fibers], dtype=float)


def process_observation(obs):
    """Extract the same N fibers from the ML sky, DRP sky and CFrame files.

    Any file that is missing or unreadable is filled with NaN. The fibers are
    chosen from the first readable file of ML sky, CFrame, SFrame; if none can
    be read, the fiber indices are -1 and everything is NaN.
    """
    sframe = obs['location']
    expnum = obs['expnum']
    # f = obs['filename']
    files = {
        'mlsky': OUTPUT_DIR+f'/lvmSFrame-skyml-block-new-{expnum}.fits',
        'drpsky': sframe,
        'flux': sframe.replace("SFrame", "CFrame"),
    }
    # Selection order: ML sky first, then the CFrame, then the SFrame.
    order = ['mlsky', 'flux', 'drpsky']

    nan = np.full((N_FIBERS, N_PIX), np.nan)
    out = {name: (nan, nan) for name in files}  # (flux, ivar) per product
    fibers = None

    for name in order:
        try:
            if fibers is None:
                fibers, f, iv = select_fibers(files[name])
            else:
                _, f, iv = select_fibers(files[name], fibers=fibers)
        except (OSError, KeyError, ValueError) as err:
            print(f"{expnum}: cannot read {files[name]}, filling with NaN: {err}")
            continue
        out[name] = (f, iv)

    if fibers is None:
        fibers = np.full(N_FIBERS, -1)

    return fibers, out['mlsky'][0], out['drpsky'][0], out['flux'][0], out['flux'][1]


def main():
    nobs = len(t)
    shape = (nobs * N_FIBERS, N_PIX)
    mlsky = np.zeros(shape)
    drpsky = np.zeros(shape)
    cframe = np.zeros(shape)
    ivar = np.zeros(shape)
    expnum = np.repeat(np.asarray(t['expnum']), N_FIBERS)
    fiber = np.zeros(nobs * N_FIBERS, dtype=int)

    with ThreadPoolExecutor(max_workers=N_SUBPROCESS) as executor:
        # Map each future to its observation so results can arrive in any order
        # and still land in the same N rows of every extension.
        futures = {executor.submit(process_observation, obs): i for i, obs in enumerate(t)}

        for future in as_completed(futures):
            i = futures[future]
            fibers, ml, drp, cf, iv = future.result()
            rows = slice(i * N_FIBERS, (i + 1) * N_FIBERS)
            fiber[rows] = fibers
            mlsky[rows] = ml
            drpsky[rows] = drp
            cframe[rows] = cf
            ivar[rows] = iv

    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU(data=mlsky, name='MLSKY'))
    hdul.append(fits.ImageHDU(data=drpsky, name='DRPSKY'))
    hdul.append(fits.ImageHDU(data=cframe, name='FLUX'))
    hdul.append(fits.ImageHDU(data=ivar, name='IVAR'))
    hdul.append(fits.BinTableHDU(Table({'expnum': expnum, 'fiber': fiber}), name='FIBERS'))
    hdul.writeto('ml_sky_evaluation.fits', overwrite=True)

if __name__ == "__main__":
    main()
