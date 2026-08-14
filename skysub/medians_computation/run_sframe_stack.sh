#!/usr/bin/env bash
set -euo pipefail

DRP_VERSION="1.2.1"

if [[ -z "${SAS_BASE_DIR:-}" ]]; then
  echo "SAS_BASE_DIR is not set" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_ROOT="${SAS_BASE_DIR}/sdsswork/lvm/spectro/redux/${DRP_VERSION}"
FILE_LIST="${SAS_BASE_DIR}/sdsswork/users/u6058164/lvmsframe_${DRP_VERSION}.dat"
OUTPUT=""
SAMPLE_SIZE=24
WORKER_COUNTS="1,2,4,8,16,32"
WORKERS=""
LIMIT=""
OVERWRITE=0
FAINT_FIBER_PERCENTILE="70"
SCI_FAINT_FIBER_PERCENTILE=""
SKY_FAINT_FIBER_PERCENTILE=""
EVERY_NTH="1"
SKIP_SCAN=0
SKIP_BENCHMARK=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --input-root PATH       Input DRP root.
  --file-list PATH        SFrame .dat list to write/read.
  --output PATH           Output FITS path. Default is derived automatically.
  --faint-fiber-percentile P
                          Fallback faint-fiber percentile applied to any telescope
                          group without a per-group override. Default: 70.
  --sci-faint-fiber-percentile P
                          Overrides --faint-fiber-percentile for the science telescope.
  --sky-faint-fiber-percentile P
                          Overrides --faint-fiber-percentile for both sky telescopes.
  --every-nth N           Process every Nth input entry. Default: 1.
  --sample-size N         Benchmark sample size. Default: ${SAMPLE_SIZE}.
  --worker-counts LIST    Benchmark worker counts. Default: ${WORKER_COUNTS}.
  --workers N             Skip recommendation and build with N workers.
  --limit N               Build only the first N rows.
  --overwrite             Replace an existing output FITS.
  --skip-scan             Reuse an existing file list.
  --skip-benchmark        Skip benchmark; requires --workers.
  -h, --help              Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-root)
      INPUT_ROOT="$2"
      shift 2
      ;;
    --file-list)
      FILE_LIST="$2"
      shift 2
      ;;
    --output)
      OUTPUT="$2"
      shift 2
      ;;
    --sample-size)
      SAMPLE_SIZE="$2"
      shift 2
      ;;
    --worker-counts)
      WORKER_COUNTS="$2"
      shift 2
      ;;
    --workers)
      WORKERS="$2"
      shift 2
      ;;
    --faint-fiber-percentile)
      FAINT_FIBER_PERCENTILE="$2"
      shift 2
      ;;
    --sci-faint-fiber-percentile)
      SCI_FAINT_FIBER_PERCENTILE="$2"
      shift 2
      ;;
    --sky-faint-fiber-percentile)
      SKY_FAINT_FIBER_PERCENTILE="$2"
      shift 2
      ;;
    --every-nth)
      EVERY_NTH="$2"
      shift 2
      ;;
    --limit)
      LIMIT="$2"
      shift 2
      ;;
    --overwrite)
      OVERWRITE=1
      shift
      ;;
    --skip-scan)
      SKIP_SCAN=1
      shift
      ;;
    --skip-benchmark)
      SKIP_BENCHMARK=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "${SKIP_SCAN}" -eq 0 ]]; then
  python "${SCRIPT_DIR}/scan_sframes.py" \
    --input-root "${INPUT_ROOT}" \
    --output "${FILE_LIST}"
fi

if [[ -z "${WORKERS}" ]]; then
  if [[ "${SKIP_BENCHMARK}" -eq 0 ]]; then
    BENCHMARK_OUTPUT="$(python "${SCRIPT_DIR}/benchmark_sframe_io.py" \
      --input-list "${FILE_LIST}" \
      --sample-size "${SAMPLE_SIZE}" \
      --workers "${WORKER_COUNTS}")"
    printf '%s\n' "${BENCHMARK_OUTPUT}"
    WORKERS="$(printf '%s\n' "${BENCHMARK_OUTPUT}" | awk -F= '/^recommended_workers=/{print $2}' | tail -n 1)"
  else
    echo "--skip-benchmark requires --workers N; no recommendation file is used" >&2
    exit 2
  fi

  if [[ -z "${WORKERS}" ]]; then
    echo "could not parse recommended_workers from benchmark output" >&2
    exit 2
  fi
fi

BUILD_ARGS=(
  --input-list "${FILE_LIST}"
  --workers "${WORKERS}"
  --faint-fiber-percentile "${FAINT_FIBER_PERCENTILE}"
  --every-nth "${EVERY_NTH}"
)

if [[ -n "${SCI_FAINT_FIBER_PERCENTILE}" ]]; then
  BUILD_ARGS+=(--sci-faint-fiber-percentile "${SCI_FAINT_FIBER_PERCENTILE}")
fi

if [[ -n "${SKY_FAINT_FIBER_PERCENTILE}" ]]; then
  BUILD_ARGS+=(--sky-faint-fiber-percentile "${SKY_FAINT_FIBER_PERCENTILE}")
fi

if [[ -n "${OUTPUT}" ]]; then
  BUILD_ARGS+=(--output "${OUTPUT}")
fi

if [[ -n "${LIMIT}" ]]; then
  BUILD_ARGS+=(--limit "${LIMIT}")
fi

if [[ "${OVERWRITE}" -eq 1 ]]; then
  BUILD_ARGS+=(--overwrite)
fi

python "${SCRIPT_DIR}/build_sframe_stack.py" "${BUILD_ARGS[@]}"
