"""Typer command-line interface for portable LVM median-stack production."""

from __future__ import annotations

import json
import logging
import os
import tempfile
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from enum import Enum
from pathlib import Path
from typing import Any

import typer
from tqdm import tqdm

from . import __version__
from .gaia import TAP_SERVICES, cache_status, combine_gaia_tables, failed_sframes, fetch_gaia
from .stack import build_stack, read_manifest, scan_sframes, write_manifest


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as fh:
        json.dump(value, fh, indent=2, sort_keys=True)
        fh.write("\n")
        temporary = Path(fh.name)
    os.replace(temporary, path)


class RunStatus:
    def __init__(self, path: Path, command: str, log_file: Path):
        self.path = path
        self.data: dict[str, Any] = {
            "command": command,
            "state": "running",
            "started": _now(),
            "updated": _now(),
            "log_file": str(log_file),
            "counts": {},
        }
        self.last_write = 0.0
        self.write(force=True)

    def write(
        self,
        counts: dict[str, Any] | None = None,
        *,
        state: str | None = None,
        error: str | None = None,
        force: bool = False,
    ) -> None:
        if counts is not None:
            self.data["counts"] = counts
        if state is not None:
            self.data["state"] = state
        if error is not None:
            self.data["error"] = error
        now = time.monotonic()
        if not force and state is None and now - self.last_write < 1.0:
            return
        self.data["updated"] = _now()
        _atomic_json(self.path, self.data)
        self.last_write = now


def _logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("lvm_medians")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    return logger


def _progress(
    total: int, description: str, disabled: bool, status: RunStatus
) -> tuple[tqdm[Any], Callable[[dict[str, int]], None]]:
    bar = tqdm(total=total, desc=description, unit="exposure", disable=disabled, dynamic_ncols=True)
    last = 0

    def update(counts: dict[str, int]) -> None:
        nonlocal last
        completed = counts.get("completed", 0)
        bar.update(completed - last)
        last = completed
        shown = {key: value for key, value in counts.items() if key not in {"total", "completed"}}
        bar.set_postfix(shown, refresh=False)
        status.write(counts)

    return bar, update


SERVICE_HELP = ", ".join(f"{name}={url}" for name, url in TAP_SERVICES.items())

app = typer.Typer(
    name="lvm-medians",
    help="Build reproducible LVM SFrame products with resumable Gaia TAP caching.",
    epilog=f"Known TAP services: {SERVICE_HELP}. A custom TAP URL is also accepted.",
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode="rich",
    pretty_exceptions_enable=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


class BuildMode(str, Enum):
    median = "median"
    faint_fibers = "faint-fibers"


class GaiaSigma(str, Enum):
    one = "1"
    three = "3"


class GaiaTable(str, Enum):
    fibers = "fibers"
    sources = "sources"


def _version(value: bool) -> None:
    if value:
        typer.echo(f"lvm-medians {__version__}")
        raise typer.Exit()


@app.callback()
def options(
    ctx: typer.Context,
    work_dir: Path = typer.Option(
        Path("lvm-medians-work"), help="State, cache, and default SFrame-list directory."
    ),
    log_file: Path | None = typer.Option(
        None, help="Log file (default: WORK_DIR/logs/lvm-medians.log)."
    ),
    no_progress: bool = typer.Option(
        False, "--no-progress", help="Disable progress bars; logs and status remain enabled."
    ),
    version: bool = typer.Option(
        False, "--version", callback=_version, is_eager=True, help="Show the version and exit."
    ),
) -> None:
    """Set options shared by all commands."""
    ctx.obj = {"work_dir": work_dir, "log_file": log_file, "no_progress": no_progress}


@contextmanager
def _tracked(ctx: typer.Context, command: str) -> Iterator[tuple[Path, logging.Logger, RunStatus]]:
    work_dir = ctx.obj["work_dir"].expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    log_file = (ctx.obj["log_file"] or work_dir / "logs/lvm-medians.log").expanduser().resolve()
    logger = _logging(log_file)
    status = RunStatus(work_dir / "run-status.json", command, log_file)
    logger.info("starting command=%s", command)
    try:
        yield work_dir, logger, status
    except typer.Exit:
        raise
    except KeyboardInterrupt:
        status.write(state="interrupted", error="interrupted by user", force=True)
        logger.warning("command interrupted by user")
        typer.echo(f"interrupted\nlog: {log_file}", err=True)
        raise typer.Exit(130) from None
    except Exception as exc:
        message = f"{type(exc).__name__}: {exc}"
        status.write(state="failed", error=message, force=True)
        logger.exception("command failed")
        typer.echo(f"error: {message}\nlog: {log_file}", err=True)
        raise typer.Exit(1) from None


def _complete(
    command: str, result: dict[str, Any], logger: logging.Logger, status: RunStatus
) -> None:
    state = (
        "complete"
        if not result.get("failed") and not result.get("error")
        else "complete_with_errors"
    )
    status.write(result, state=state, force=True)
    logger.info("completed command=%s result=%s", command, result)
    typer.echo(json.dumps(result, indent=2, sort_keys=True))
    if state != "complete":
        raise typer.Exit(1)


@app.command()
def scan(
    ctx: typer.Context,
    sframes_root: Path | None = typer.Option(
        None, help="Directory tree containing lvmSFrame FITS files."
    ),
    input_list: Path | None = typer.Option(
        None, help="Validate and normalize an existing path list."
    ),
    sframe_list: Path | None = typer.Option(
        None, help="SFrame list (default: WORK_DIR/sframes.txt)."
    ),
    every_nth: int = typer.Option(1, min=1, help="Keep every Nth discovered path."),
    limit: int | None = typer.Option(None, min=1, help="Keep at most this many paths."),
) -> None:
    """Find SFrames and write their paths to a deterministic list."""
    if (sframes_root is None) == (input_list is None):
        raise typer.BadParameter(
            "provide exactly one source", param_hint="--sframes-root / --input-list"
        )
    with _tracked(ctx, "scan") as (work_dir, logger, status):
        destination = (sframe_list or work_dir / "sframes.txt").expanduser().resolve()
        with tqdm(
            desc="Scanning", unit="SFrame", disable=ctx.obj["no_progress"], dynamic_ncols=True
        ) as bar:
            if sframes_root is not None:
                last = 0

                def found(count: int) -> None:
                    nonlocal last
                    bar.update(count - last)
                    last = count
                    status.write({"discovered": count})

                paths = scan_sframes(sframes_root.expanduser().resolve(), found)
            else:
                paths = [path for _, path in read_manifest(input_list.expanduser().resolve())]
                bar.update(len(paths))
                status.write({"discovered": len(paths)}, force=True)
        discovered = len(paths)
        paths = paths[::every_nth]
        if limit is not None:
            paths = paths[:limit]
        status.write({"discovered": discovered, "selected": len(paths)}, force=True)
        if not paths:
            raise RuntimeError("no lvmSFrame files found")
        write_manifest(paths, destination)
        _complete("scan", {"files": len(paths), "sframe_list": str(destination)}, logger, status)


@app.command("fetch-gaia")
def fetch_gaia_command(
    ctx: typer.Context,
    sframe_list: Path | None = typer.Option(
        None, help="SFrame list (default: WORK_DIR/sframes.txt)."
    ),
    cache_dir: Path | None = typer.Option(None, help="Gaia cache root (default: WORK_DIR/gaia)."),
    tap_service: str = typer.Option(
        "ari", help=f"TAP alias or URL. Built-in aliases: {SERVICE_HELP}."
    ),
    query_workers: int = typer.Option(5, min=1, help="Parallel TAP requests."),
    retries: int = typer.Option(3, min=0, help="Retries after a network or query error."),
    retry_failed: bool = typer.Option(
        False, "--retry-failed", help="Process only exposures in gaia-failures.jsonl."
    ),
    timeout: float = typer.Option(120.0, min=0.1, help="HTTP timeout in seconds."),
    maxrec: int = typer.Option(1_000_000, min=1, help="TAP MAXREC value."),
    token_env: str | None = typer.Option(
        None, help="Environment variable containing an optional bearer token."
    ),
    passband: Path | None = typer.Option(None, help="Override the bundled Gaia DR3 G VOTable."),
    every_nth: int = typer.Option(1, min=1, help="Process every Nth SFrame-list entry."),
    limit: int | None = typer.Option(None, min=1, help="Process at most this many entries."),
) -> None:
    """Download missing Gaia data and build per-fiber cache files."""
    with _tracked(ctx, "fetch-gaia") as (work_dir, logger, status):
        sframe_list = (sframe_list or work_dir / "sframes.txt").expanduser().resolve()
        cache_dir = (cache_dir or work_dir / "gaia").expanduser().resolve()
        selected = (
            failed_sframes(sframe_list, cache_dir, every_nth, limit)
            if retry_failed
            else read_manifest(sframe_list, every_nth, limit)
        )
        total = len(selected)
        bar, callback = _progress(total, "Gaia", ctx.obj["no_progress"], status)
        try:
            token = os.environ.get(token_env) if token_env else None
            if token_env and token is None:
                raise ValueError(f"token environment variable is not set: {token_env}")
            result = fetch_gaia(
                sframe_list,
                cache_dir,
                service=tap_service,
                workers=query_workers,
                retries=retries,
                timeout=timeout,
                maxrec=maxrec,
                token=token,
                passband=passband,
                every_nth=every_nth,
                limit=limit,
                retry_failed=retry_failed,
                progress=callback,
            )
        finally:
            bar.close()
        _complete("fetch-gaia", result, logger, status)


@app.command("combine-gaia")
def combine_gaia_command(
    ctx: typer.Context,
    output: Path = typer.Option(..., help="Combined .fits or .parquet table."),
    table: GaiaTable = typer.Option(
        GaiaTable.fibers,
        "--table",
        help="Combine derived per-fiber rows or raw Gaia source matches.",
    ),
    sframe_list: Path | None = typer.Option(
        None, help="SFrame list (default: WORK_DIR/sframes.txt)."
    ),
    cache_dir: Path | None = typer.Option(None, help="Gaia cache root (default: WORK_DIR/gaia)."),
    every_nth: int = typer.Option(1, min=1, help="Process every Nth SFrame-list entry."),
    limit: int | None = typer.Option(None, min=1, help="Process at most this many entries."),
    overwrite: bool = typer.Option(False, help="Replace an existing combined table."),
) -> None:
    """Combine per-exposure Gaia caches into one FITS or Parquet table."""
    with _tracked(ctx, "combine-gaia") as (work_dir, logger, status):
        sframe_list = (sframe_list or work_dir / "sframes.txt").expanduser().resolve()
        cache_dir = (cache_dir or work_dir / "gaia").expanduser().resolve()
        total = len(read_manifest(sframe_list, every_nth, limit))
        bar, callback = _progress(total, "Combining Gaia", ctx.obj["no_progress"], status)
        try:
            result = combine_gaia_tables(
                sframe_list,
                cache_dir,
                output.expanduser().resolve(),
                table_kind=table.value,
                every_nth=every_nth,
                limit=limit,
                overwrite=overwrite,
                progress=callback,
            )
        finally:
            bar.close()
        _complete("combine-gaia", result, logger, status)


@app.command("build-medians")
def build_medians(
    ctx: typer.Context,
    output: Path = typer.Option(..., help="Destination FITS file."),
    sframe_list: Path | None = typer.Option(
        None, help="SFrame list (default: WORK_DIR/sframes.txt)."
    ),
    mode: BuildMode = typer.Option(BuildMode.median, help="Output aggregation mode."),
    workers: int = typer.Option(4, min=1, help="Local worker-process count."),
    every_nth: int = typer.Option(1, min=1, help="Process every Nth SFrame-list entry."),
    limit: int | None = typer.Option(None, min=1, help="Process at most this many entries."),
    sci_percentile: float = typer.Option(70.0, help="Faintest science-fiber percentile."),
    sky_percentile: float = typer.Option(70.0, help="Faintest sky-fiber percentile."),
    gaia_sigma: GaiaSigma | None = typer.Option(None, help="Gaia flux cut: 1 or 3 sigma."),
    gaia_ratio_threshold: float | None = typer.Option(
        None, help="Maximum Gaia/LVM (FLUX + SKY) G-band flux ratio, for example 0.1."
    ),
    gaia_fibers_dir: Path | None = typer.Option(
        None, help="Per-fiber Gaia cache (default: WORK_DIR/gaia/fibers)."
    ),
    fibers_per_telescope: int = typer.Option(3, min=1, help="Fibers per telescope in faint mode."),
    temp_dir: Path | None = typer.Option(None, help="Directory for temporary array files."),
    keep_temp: bool = typer.Option(False, help="Keep temporary arrays after the build."),
    overwrite: bool = typer.Option(False, help="Replace an existing output file."),
) -> None:
    """Build median spectra or faint-fiber rows into a FITS product."""
    if gaia_sigma is not None and gaia_ratio_threshold is not None:
        raise typer.BadParameter(
            "Gaia sigma and ratio selections are mutually exclusive",
            param_hint="--gaia-sigma / --gaia-ratio-threshold",
        )
    if gaia_ratio_threshold is not None and not 0 < gaia_ratio_threshold < 1:
        raise typer.BadParameter(
            "must be greater than 0 and less than 1", param_hint="--gaia-ratio-threshold"
        )
    with _tracked(ctx, "build-medians") as (work_dir, logger, status):
        sframe_list = (sframe_list or work_dir / "sframes.txt").expanduser().resolve()
        total = len(read_manifest(sframe_list, every_nth, limit))
        bar, callback = _progress(total, "Medians", ctx.obj["no_progress"], status)
        try:
            result = build_stack(
                sframe_list,
                output.expanduser().resolve(),
                mode=mode.value,
                workers=workers,
                every_nth=every_nth,
                limit=limit,
                sci_percentile=sci_percentile,
                sky_percentile=sky_percentile,
                gaia_dir=(gaia_fibers_dir or work_dir / "gaia/fibers").expanduser().resolve(),
                gaia_sigma=int(gaia_sigma.value) if gaia_sigma is not None else None,
                gaia_ratio_threshold=gaia_ratio_threshold,
                fibers_per_telescope=fibers_per_telescope,
                temp_dir=temp_dir.expanduser().resolve() if temp_dir else None,
                keep_temp=keep_temp,
                overwrite=overwrite,
                progress=callback,
            )
        finally:
            bar.close()
        _complete("build-medians", result, logger, status)


def _print_status(work_dir: Path, sframe_list: Path, cache_dir: Path) -> None:
    run_path = work_dir / "run-status.json"
    if run_path.exists():
        run = json.loads(run_path.read_text(encoding="utf-8"))
        counts = dict(run.get("counts", {}))
        if "manifest" in counts:  # Read status files created by versions before 0.2.
            counts["sframe_list"] = counts.pop("manifest")
        typer.echo(f"Last run: {run.get('command', 'unknown')}")
        typer.echo(f"State: {run.get('state', 'unknown')}")
        typer.echo(f"Updated: {run.get('updated', 'unknown')}")
        if counts:
            typer.echo("Result:")
            for key, value in counts.items():
                label = {
                    "files": "SFrames selected",
                    "sframe_list": "SFrame list",
                }.get(key, key.replace("_", " "))
                typer.echo(f"  {label}: {value}")
        if run.get("error"):
            typer.echo(f"Error: {run['error']}")
        if run.get("log_file"):
            typer.echo(f"Log: {run['log_file']}")
    else:
        typer.echo("Last run: none")
    if sframe_list.exists():
        values = cache_status(sframe_list, cache_dir)
        typer.echo("Gaia cache:")
        typer.echo(f"  SFrames in list: {values['sframes']}")
        typer.echo(f"  ready: {values['ready']}")
        typer.echo(f"  skipped (not flux calibrated): {values['skipped']}")
        typer.echo(f"  awaiting download: {values['awaiting_download']}")
        typer.echo(f"  failures: {values['failures']}")
    else:
        typer.echo(f"SFrame list has not been created: {sframe_list}")


@app.command()
def status(
    ctx: typer.Context,
    sframe_list: Path | None = typer.Option(
        None, help="SFrame list (default: WORK_DIR/sframes.txt)."
    ),
    cache_dir: Path | None = typer.Option(None, help="Gaia cache root (default: WORK_DIR/gaia)."),
) -> None:
    """Show the last run and Gaia cache coverage."""
    work_dir = ctx.obj["work_dir"].expanduser().resolve()
    _print_status(
        work_dir,
        (sframe_list or work_dir / "sframes.txt").expanduser().resolve(),
        (cache_dir or work_dir / "gaia").expanduser().resolve(),
    )


if __name__ == "__main__":
    app()
