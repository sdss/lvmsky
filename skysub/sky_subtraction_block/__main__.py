"""Create a compact SFrame with a sky-only global LSF per spectrograph/channel."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import logging


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input');p.add_argument('--output');p.add_argument('--seed',help='Same-exposure compact SkyML preparation; reuse medians, ML and corrections')
    p.add_argument('--ensemble',help='Trusted ML archive, required when creating preparation')
    p.add_argument('--gaia-dir');p.add_argument('--jobs',type=int,default=8)
    p.add_argument('--slit-count',type=int,default=4)
    p.add_argument('--centered',action='store_true',help='Also constrain every block profile to zero first moment (default: unconstrained centroid)')
    p.add_argument('--keep-preparation',action='store_true',help='Also write the in-memory preparation next to the output so later runs can reuse it with --seed (default: not written)')
    a=p.parse_args();logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    from .pipeline import run
    print(run(a.input,output=a.output,seed=a.seed,ensemble=a.ensemble,gaia_dir=a.gaia_dir,jobs=a.jobs,slit_count=a.slit_count,centered=a.centered,keep_preparation=a.keep_preparation))


if __name__=='__main__':main()
