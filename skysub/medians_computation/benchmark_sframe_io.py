#!/usr/bin/env python
"""Benchmark parallel reads of lvmSFrame FITS files."""

from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import median

import numpy as np
from astropy.io import fits


DRP_VERSION = "1.2.1"
DEFAULT_RELATIVE_LIST = f"sdsswork/users/u6058164/lvmsframe_{DRP_VERSION}.dat"
DEFAULT_WORKER_COUNTS = (1, 2, 4, 8, 16, 32)


def sas_path(relative_path: str) -> Path:
    sas_base = os.environ.get("SAS_BASE_DIR")
    if not sas_base:
        raise RuntimeError("SAS_BASE_DIR is not set")
    return Path(sas_base) / relative_path


def default_file_list() -> Path:
    return sas_path(DEFAULT_RELATIVE_LIST)


def read_file_list(path: Path, limit: int | None = None) -> list[Path]:
    with path.open("r", encoding="utf-8") as handle:
        paths = [Path(line.strip()) for line in handle if line.strip()]
    if limit is not None:
        return paths[:limit]
    return paths


def force_read(path: Path) -> tuple[int, float]:
    """Read the large SFrame arrays and return bytes touched plus a checksum."""
    checksum = 0.0
    bytes_read = 0
    with fits.open(path, memmap=True) as hdul:
        for name in ("WAVE", "FLUX", "SKY", "LSF"):
            data = np.asarray(hdul[name].data)
            bytes_read += int(data.nbytes)
            checksum += float(np.nanmean(data))
        bytes_read += int(hdul["SLITMAP"].data.nbytes)
    return bytes_read, checksum


def run_once(paths: list[Path], workers: int) -> dict[str, float | int]:
    workers = max(1, min(workers, len(paths)))
    start = time.perf_counter()
    bytes_read = 0
    checksum = 0.0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(force_read, path) for path in paths]
        for future in as_completed(futures):
            file_bytes, file_checksum = future.result()
            bytes_read += file_bytes
            checksum += file_checksum
    elapsed = max(time.perf_counter() - start, 1.0e-9)
    return {
        "workers": workers,
        "files": len(paths),
        "seconds": elapsed,
        "mb_read": bytes_read / 1_000_000.0,
        "mb_per_second": bytes_read / 1_000_000.0 / elapsed,
        "checksum": checksum,
    }


def recommend_worker_count(results: list[dict[str, float | int]]) -> int:
    max_rate = max(float(row["mb_per_second"]) for row in results)
    threshold = 0.90 * max_rate
    eligible = [
        int(row["workers"])
        for row in results
        if float(row["mb_per_second"]) >= threshold
    ]
    return min(eligible)


def parse_worker_counts(value: str) -> list[int]:
    counts = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        counts.append(int(part))
    if not counts:
        raise argparse.ArgumentTypeError("worker count list is empty")
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark parallel FITS I/O on a sample of SFrames."
    )
    parser.add_argument(
        "--input-list",
        type=Path,
        default=None,
        help=f"Input .dat list. Default: $SAS_BASE_DIR/{DEFAULT_RELATIVE_LIST}",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=24,
        help="Number of files from the start of the list to benchmark.",
    )
    parser.add_argument(
        "--workers",
        type=parse_worker_counts,
        default=list(DEFAULT_WORKER_COUNTS),
        help="Comma-separated worker counts to test.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Number of times to repeat each worker count; median MB/s is reported.",
    )
    parser.add_argument(
        "--recommendation-file",
        type=Path,
        default=None,
        help="Optional file to write the recommended worker count.",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        help="Optional path to write full benchmark results as JSON.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_list = args.input_list if args.input_list is not None else default_file_list()
    if args.sample_size <= 0:
        raise ValueError("--sample-size must be positive")
    if args.repeats <= 0:
        raise ValueError("--repeats must be positive")

    paths = read_file_list(input_list, limit=args.sample_size)
    if not paths:
        raise RuntimeError(f"no input paths found in {input_list}")

    results: list[dict[str, float | int | list[float]]] = []
    for requested_workers in args.workers:
        repeat_rows = [run_once(paths, requested_workers) for _ in range(args.repeats)]
        mbps_values = [float(row["mb_per_second"]) for row in repeat_rows]
        row = dict(repeat_rows[0])
        row["requested_workers"] = requested_workers
        row["mb_per_second"] = median(mbps_values)
        row["repeat_mb_per_second"] = mbps_values
        results.append(row)
        print(
            f"workers={row['workers']:>2} files={row['files']:>3} "
            f"read={row['mb_read']:.1f} MB rate={row['mb_per_second']:.1f} MB/s"
        )

    recommendation = recommend_worker_count(results)
    print(f"recommended_workers={recommendation}")
    if args.recommendation_file is not None:
        args.recommendation_file.parent.mkdir(parents=True, exist_ok=True)
        args.recommendation_file.write_text(f"{recommendation}\n", encoding="utf-8")
        print(f"wrote recommendation to {args.recommendation_file}")

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        payload = {"recommendation": recommendation, "results": results}
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"wrote JSON results to {args.json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
