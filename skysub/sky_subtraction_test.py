#! /usr/bin/env python

import os

# Maximum number of sky-subtraction subprocesses running concurrently.
N_SUBPROCESS = 4

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
DRPALL = os.getenv('LVM_SPECTRO_REDUX')+'/1.3.2/drpall-1.3.2.fits'
t = Table.read(DRPALL)

select = ((t['mjd']>60355) * (t['mjd']<60402) * (t['exptime']==900) * (t['tileid']>11111) * (t['fluxcal']=='MOD'))
t = t[select]
t = t[2:-1]

print(len(t))

def process_observation(obs):
    """Run the original loop body for one observation."""
    d = obs['location']
    expnum = obs['expnum']
    f = obs['filename']
    d = d.replace("SFrame", "CFrame")
    cmd = ["python",\
           "-m", "skysub.sky_subtraction_block",\
           os.getenv("SAS_BASE_DIR")+'/'+d,\
           "--ensemble",\
           "./mlp_ensemble_split_zodi_cont-1.3.2.pt",\
           "--gaia-dir",\
           os.getenv("SAS_BASE_DIR")+"/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/gaia/fibers",\
           "--jobs", "32",\
           "--output",\
           os.getenv("SAS_BASE_DIR")+f"/sdsswork/users/u6045701/skyml_block/lvmSFrame-skyml-block-{expnum}.fits"]

    print(cmd)
    subprocess.run(cmd)
    return expnum


def main():
    with ThreadPoolExecutor(max_workers=N_SUBPROCESS) as executor:
        futures = []
        for obs in t[1:10]:
            futures.append(executor.submit(process_observation, obs))

        for line, future in enumerate(as_completed(futures), start=1):
            expnum = future.result()
            print(line, expnum)


if __name__ == "__main__":
    main()
