#!/usr/bin/env python
"""Scan a DRP redux tree for lvmSFrame FITS files."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


DRP_VERSION = "1.2.1"
DEFAULT_RELATIVE_ROOT = f"sdsswork/lvm/spectro/redux/{DRP_VERSION}"
DEFAULT_RELATIVE_LIST = f"sdsswork/users/u6058164/lvmsframe_{DRP_VERSION}.dat"


def sas_path(relative_path: str) -> Path:
    sas_base = os.environ.get("SAS_BASE_DIR")
    if not sas_base:
        raise RuntimeError("SAS_BASE_DIR is not set")
    return Path(sas_base) / relative_path


def default_input_root() -> Path:
    return sas_path(DEFAULT_RELATIVE_ROOT)


def default_output_list() -> Path:
    return sas_path(DEFAULT_RELATIVE_LIST)


def scan_sframes(input_root: Path) -> list[Path]:
    if not input_root.exists():
        raise FileNotFoundError(f"input root does not exist: {input_root}")
    if not input_root.is_dir():
        raise NotADirectoryError(f"input root is not a directory: {input_root}")

    paths: list[Path] = []
    for root, dirs, files in os.walk(input_root):
        dirs.sort()
        for filename in sorted(files):
            if filename.startswith("lvmSFrame-") and filename.endswith(".fits"):
                paths.append(Path(root) / filename)
    return sorted(paths, key=lambda path: path.as_posix())


def write_list(paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{path}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a deterministic list of lvmSFrame FITS files."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=None,
        help=f"DRP root to scan. Default: $SAS_BASE_DIR/{DEFAULT_RELATIVE_ROOT}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"Path to write. Default: $SAS_BASE_DIR/{DEFAULT_RELATIVE_LIST}",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_root = args.input_root if args.input_root is not None else default_input_root()
    output_path = args.output if args.output is not None else default_output_list()

    paths = scan_sframes(input_root)
    write_list(paths, output_path)
    print(f"wrote {len(paths)} paths to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
