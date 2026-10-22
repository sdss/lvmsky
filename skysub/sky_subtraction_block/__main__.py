"""Subtract sky from one CFrame with a sky-only global block LSF."""
from __future__ import annotations

import argparse
import logging

from .options import DEFAULTS, HELP
from .resources import limit_environment

limit_environment()


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


class LegacyParser(argparse.ArgumentParser):
    def _get_option_tuples(self, option_string):
        matches = super()._get_option_tuples(option_string)
        previous = [m for m in matches if m[1] not in ("--ml-model-file", "--reuse-state")]
        # Preserve formerly unambiguous abbreviations such as --r for --reference-prior.
        return previous if len(previous) == 1 else matches


def build_parser():
    description = __doc__ + """

Examples:
  python -m skysub.sky_subtraction_block /data/lvmCFrame-00011017.fits \\
    --ml-model-file /models/sky.pt --jobs 1 \\
    --output /results/lvmSFrame-skyml-block-00011017.fits

  python -m skysub.sky_subtraction_block /data/lvmCFrame-00011017.fits \\
    --ensemble /models/sky.pt --keep-preparation

The original one-exposure interface remains compatible. For explicit single and
batch commands: python -m skysub.sky_subtraction_block.cli --help
Ordinary operation computes intermediate state in memory and writes the final
SFrame; saved-state reuse is a special optional mode.
"""
    parser = LegacyParser(description=description, formatter_class=HelpFormatter,
                          usage="python -m skysub.sky_subtraction_block INPUT [OPTIONS]")
    parser.add_argument("input", help="One observed-flux lvmCFrame FITS path.")
    parser.add_argument("--output", help="One final lvmSFrame-skyml-block-* FITS; default beside INPUT. Existing files are protected.")
    parser.add_argument("--ml-model-file", "--ensemble", dest="ensemble", help=HELP["ensemble"])
    parser.add_argument("--gaia-dir", help=HELP["gaia_dir"])
    parser.add_argument("--reuse-state", "--seed", dest="seed", help=HELP["seed"])
    parser.add_argument("--jobs", type=int, default=8, help=HELP["jobs"])
    for name in ("slit_count", "offset_count", "wave_count"):
        parser.add_argument("--" + name.replace("_", "-"), type=int, default=DEFAULTS[name], help=HELP[name])
    for name in ("offset_roughness", "wave_roughness", "slit_roughness", "reference_prior"):
        parser.add_argument("--" + name.replace("_", "-"), type=float, default=DEFAULTS[name], help=HELP[name])
    parser.add_argument("--centered", action="store_true", help=HELP["centered"])
    parser.add_argument("--keep-preparation", action="store_true", help=HELP["keep_preparation"])
    return parser


def main():
    args = build_parser().parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    from .pipeline import run
    kwargs = vars(args).copy()
    input_path = kwargs.pop("input")
    print(run(input_path, **kwargs))


if __name__ == "__main__":
    main()
