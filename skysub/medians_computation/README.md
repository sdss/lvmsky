# LVM medians CLI

`lvm-medians` creates reproducible summary FITS products from LVM SFrame
exposures. It discovers input files, downloads and caches optional Gaia DR3
matches, can combine Gaia caches into one large table, and builds median spectra
from selected fibers in each exposure.

## What each command does

| Command | Purpose |
| --- | --- |
| `scan` | Find SFrames and write their sorted paths to a plain-text SFrame list. It does not process spectra. |
| `fetch-gaia` | Query TAP and write separate raw-source and derived-fiber cache files for each exposure. It does not combine those files. |
| `combine-gaia` | Combine all available per-exposure Gaia caches selected by the SFrame list into one FITS or Parquet table. It makes no network requests. |
| `build-medians` | Select fibers and build median spectra (or faint-fiber rows) from the SFrames. |
| `status` | Show the latest run state, processing counts, and Gaia cache coverage. |

The **SFrame list** is a plain-text file with one SFrame path per line.
`scan` creates it as `WORK_DIR/sframes.txt` by default. Every later command reads
the same list, so it processes the same ordered set of exposures without
scanning the data directory again.

## Run from this directory

Python 3.11 or newer is required.

From `/home/ik52/sdss5/lvmsky/skysub/medians_computation`, run the CLI without
installing the package:

```bash
PYTHONPATH="$PWD/src" python -m lvm_medians --help
```

The active Python environment must contain the dependencies. To create an
isolated environment and install the command instead:

```bash
cd /home/ik52/sdss5/lvmsky/skysub/medians_computation
python -m venv .venv
. .venv/bin/activate
python -m pip install -e .
lvm-medians --help
```

Use `python -m pip install -e '.[parquet]'` when Parquet output is needed.
FITS output is available in the base installation. On another server, use the
same commands after changing the repository and SFrame paths.
In the examples below, replace `lvm-medians` with
`PYTHONPATH="$PWD/src" python -m lvm_medians` to run without installation.

## Quick smoke test on three SFrames

The following commands were tested with the data path shown here. The first
smoke test is entirely local and does not contact a TAP service.

```bash
cd /home/ik52/sdss5/lvmsky/skysub/medians_computation

# Find files, sort them, and keep only three paths in the SFrame list.
lvm-medians --work-dir /tmp/lvm-medians-smoke scan \
  --sframes-root /data/sas/sdsswork/lvm/spectro/redux/1.2.1dev_telluric_correction/0011XX/11111 \
  --limit 3

# Inspect the exact input list.
cat /tmp/lvm-medians-smoke/sframes.txt

# Build median spectra from those three SFrames.
lvm-medians --work-dir /tmp/lvm-medians-smoke build-medians \
  --workers 1 \
  --output /tmp/lvm-medians-smoke/smoke-median.fits \
  --overwrite

lvm-medians --work-dir /tmp/lvm-medians-smoke status
```

A successful build reports `complete: true`, `ok: 3`, and `error: 0`.
`scan --limit 3` still traverses the selected directory so paths can be sorted,
but only three are written to the SFrame list. Add `--every-nth 10` to sample every
tenth sorted path before applying the limit.

### Add Gaia to the same smoke test

```bash
# Create or resume separate cache files for all three exposures.
lvm-medians --work-dir /tmp/lvm-medians-smoke fetch-gaia \
  --query-workers 5

# Combine derived per-fiber caches into one large FITS table.
lvm-medians --work-dir /tmp/lvm-medians-smoke combine-gaia \
  --table fibers \
  --output /tmp/lvm-medians-smoke/gaia-fibers-all.fits \
  --overwrite

# Build Gaia-filtered median spectra from the same cache.
lvm-medians --work-dir /tmp/lvm-medians-smoke build-medians \
  --workers 1 \
  --gaia-ratio-threshold 0.1 \
  --output /tmp/lvm-medians-smoke/smoke-median-gaia.fits \
  --overwrite
```

To combine the raw Gaia matches instead, pass `--table sources`. To write
Parquet, install the `parquet` extra and use a `.parquet` output name:

```bash
lvm-medians --work-dir /tmp/lvm-medians-smoke combine-gaia \
  --table fibers \
  --output /tmp/lvm-medians-smoke/gaia-fibers-all.parquet \
  --overwrite
```

`fetch-gaia --retry-failed` reads `gaia-failures.jsonl` and processes only
previously failed exposures. A normal rerun also reuses valid caches, but checks
every SFrame-list entry. `combine-gaia` can then rebuild a single table entirely
offline.

## Command reference

Global options must precede the command:

```text
lvm-medians [GLOBAL OPTIONS] COMMAND [COMMAND OPTIONS]
```

| Global option | Meaning |
| --- | --- |
| `--work-dir PATH` | State/cache directory; default `./lvm-medians-work`. |
| `--log-file PATH` | Log file; default `WORK_DIR/logs/lvm-medians.log`. |
| `--no-progress` | Hide progress bars while retaining logs and status updates. |
| `--version` | Print the package version and exit. |
| `-h`, `--help` | Show help and exit. |

### `scan`

```text
lvm-medians scan (--sframes-root PATH | --input-list PATH) [OPTIONS]
```

| Option | Meaning |
| --- | --- |
| `--sframes-root PATH` | Recursively find `lvmSFrame-*.fits*` files below this directory. |
| `--input-list PATH` | Validate and normalize paths from an existing text SFrame list. |
| `--sframe-list PATH` | Output SFrame list; default `WORK_DIR/sframes.txt`. |
| `--every-nth INTEGER` | Keep every Nth path after sorting; default `1`. |
| `--limit INTEGER` | Keep at most this many paths in the output SFrame list. |
| `-h`, `--help` | Show command help. |

Relative paths in `--input-list` are resolved relative to that list file.
Sampling is applied after sorting, so the same input produces the same SFrame list.

### `fetch-gaia`

```text
lvm-medians fetch-gaia [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--sframe-list PATH` | `WORK_DIR/sframes.txt` | Plain-text list of SFrame paths. |
| `--cache-dir PATH` | `WORK_DIR/gaia` | Raw and derived Gaia cache root. |
| `--tap-service ALIAS_OR_URL` | `ari` | Built-in alias or a complete VO TAP base URL. |
| `--query-workers INTEGER` | `5` | Concurrent TAP requests. |
| `--retries INTEGER` | `3` | Retries after each network/query error. |
| `--retry-failed` | off | Process only exposures recorded in `gaia-failures.jsonl`. |
| `--timeout FLOAT` | `120` | HTTP timeout in seconds. |
| `--maxrec INTEGER` | `1000000` | TAP `MAXREC` limit. |
| `--token-env NAME` | unset | Read an optional bearer token from this environment variable. |
| `--passband PATH` | bundled Gaia DR3 G | Override the passband VOTable. |
| `--every-nth INTEGER` | `1` | Process every Nth SFrame path in the list. |
| `--limit INTEGER` | all | Process at most this many selected entries. |
| `-h`, `--help` |  | Show command help. |

Each successful exposure creates:

- `gaia/sources/lvmGAIA-sources-EXPNUM.fits`: one row per Gaia source matched
  to an exposure fiber;
- `gaia/fibers/lvmGAIA-fibers-EXPNUM.fits`: one row per exposure fiber with
  derived Gaia and synthetic LVM G-band measurements.

Built-in TAP services:

| Alias | URL |
| --- | --- |
| `ari` | `https://gaia.ari.uni-heidelberg.de/tap` |
| `esa` | `https://gea.esac.esa.int/tap-server/tap` |
| `aip` | `https://gaia.aip.de/tap` |

### `combine-gaia`

```text
lvm-medians combine-gaia --output PATH [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--output PATH` | required | Combined `.fits`, `.fit`, or `.parquet` table. |
| `--table fibers\|sources` | `fibers` | Combine derived fiber rows or raw source matches. |
| `--sframe-list PATH` | `WORK_DIR/sframes.txt` | Defines the exposures and deterministic row order. |
| `--cache-dir PATH` | `WORK_DIR/gaia` | Per-exposure Gaia cache root. |
| `--every-nth INTEGER` | `1` | Process every Nth SFrame path in the list. |
| `--limit INTEGER` | all | Process at most this many selected entries. |
| `--overwrite`, `--no-overwrite` | off | Permit or prevent replacement of an existing table. |
| `-h`, `--help` |  | Show command help. |

The command writes a partial table with `COMPLETE = false` and exits non-zero
when some selected cache files are unavailable or corrupt. Details are written
to the run log. FITS includes checksums and primary-header provenance; Parquet
stores equivalent metadata in the table schema.

### `build-medians`

```text
lvm-medians build-medians --output PATH [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--output PATH` | required | Destination FITS product. |
| `--sframe-list PATH` | `WORK_DIR/sframes.txt` | Plain-text list of SFrame paths. |
| `--mode median\|faint-fibers` | `median` | Output aggregation mode. |
| `--workers INTEGER` | `4` | Local worker processes. |
| `--every-nth INTEGER` | `1` | Process every Nth SFrame path in the list. |
| `--limit INTEGER` | all | Process at most this many selected entries. |
| `--sci-percentile FLOAT` | `70` | Faintest science-fiber percentage in median mode. |
| `--sky-percentile FLOAT` | `70` | Faintest SkyE/SkyW fiber percentage in median mode. |
| `--gaia-sigma 1\|3` | unset | Keep Gaia flux below 1 or 3 times the formal LVM G-band uncertainty. |
| `--gaia-ratio-threshold FLOAT` | unset | Keep Gaia/LVM `(FLUX + SKY)` G-band ratio at or below this value; range `(0, 1)`. |
| `--gaia-fibers-dir PATH` | `WORK_DIR/gaia/fibers` | Per-exposure derived Gaia cache directory. |
| `--fibers-per-telescope INTEGER` | `3` | Fibers retained per Sci/SkyE/SkyW exposure in `faint-fibers` mode. |
| `--temp-dir PATH` | output directory | Directory for temporary arrays. |
| `--keep-temp`, `--no-keep-temp` | off | Retain or remove temporary arrays. |
| `--overwrite`, `--no-overwrite` | off | Permit or prevent replacement of an existing product. |
| `-h`, `--help` |  | Show command help. |

`--gaia-ratio-threshold 0.1` keeps a fiber when its complete Gaia G-band
catalog flux divided by its synthetic LVM `(FLUX + SKY)` G-band flux is at
most `0.1` (10%). Fibers with no Gaia source are also kept; incomplete Gaia
matches are rejected. This option and `--gaia-sigma` are mutually exclusive.
Gaia selections apply only to `median` mode.

### `status`

```text
lvm-medians status [--sframe-list PATH] [--cache-dir PATH]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--sframe-list PATH` | `WORK_DIR/sframes.txt` | SFrame list used to calculate coverage. |
| `--cache-dir PATH` | `WORK_DIR/gaia` | Gaia cache root. |
| `-h`, `--help` |  | Show command help. |

`SFrames in list` is the number selected by `scan`; `ready` is the number with
a per-exposure Gaia cache; `awaiting download` is the number still needing a
successful `fetch-gaia` request; and `failures` is the number currently listed
in `gaia/gaia-failures.jsonl`.

## Runtime files and recovery

```text
lvm-medians-work/
├── sframes.txt
├── run-status.json
├── logs/lvm-medians.log        # detailed run log
└── gaia/
    ├── gaia-failures.jsonl     # unresolved network/query failures
    ├── sources/                # raw per-exposure TAP results
    └── fibers/                 # derived per-exposure fiber measurements
```

Combined Gaia tables are written exactly to the `combine-gaia --output` path.
Progress bars and `run-status.json` show completed, cached, combined, skipped,
and failed counts. Detailed exceptions go to the log.
`gaia-failures.jsonl` records unresolved TAP exposures; run
`fetch-gaia --retry-failed` and then rerun `combine-gaia`.

Median products contain `WAVE`, seven flux/LSF image extensions, `META`, and
`INPUT_STATUS`. Faint-fiber products contain `WAVE`, `FLUX`, combined `IVAR`,
`META`, and `INPUT_STATUS`. FITS checksums and atomic replacement protect output
files. Partial products remain available for inspection and are marked with
`COMPLETE = false`.
