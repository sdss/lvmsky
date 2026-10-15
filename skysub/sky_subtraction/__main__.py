"""Run one CFrame without overwriting a DRP product."""

import os

os.environ.setdefault('MPLBACKEND', 'Agg')
# One native thread per numerical worker; must be set before numpy is imported.
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse  # noqa: E402
import logging  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='lvmCFrame FITS input')
    parser.add_argument('--ensemble', required=True, help='Trusted trained ensemble .pt archive')
    parser.add_argument('--gaia-dir', help='Existing per-fiber Gaia cache; ARI TAP otherwise')
    parser.add_argument('--lsf-calibration', help='Survey ASCII LSF directory (bundled by default)')
    parser.add_argument('--output', help='Output path; default is adjacent lvmSFrame-skyml-*.fits')
    parser.add_argument('--jobs', type=int, default=1, help='Fiber fit worker processes')
    parser.add_argument('--survey-only', action='store_true', help='Use survey R/Z for the validation baseline')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    from .pipeline import run

    output = run(
        args.input,
        args.ensemble,
        gaia_dir=args.gaia_dir,
        lsf_calibration=args.lsf_calibration,
        output=args.output,
        jobs=args.jobs,
        survey_only=args.survey_only,
    )
    logging.getLogger(__name__).info('Published %s', output)


if __name__ == '__main__':
    main()
