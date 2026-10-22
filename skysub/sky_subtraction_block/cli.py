"""Modern single-exposure and batch CLI; the original argparse entry point remains supported."""
from __future__ import annotations

from enum import Enum
import logging
import os
from pathlib import Path
import subprocess
from typing import Annotated

import typer

from .options import DEFAULTS, HELP, validate_model
from .resources import apply_affinity, check_runtime_dependencies, cpu_budget, limit_environment

limit_environment()

ROOT_HELP = """Subtract sky from calibrated CFrame exposures with the sky-only block LSF method.

Examples:

    python -m skysub.sky_subtraction_block.cli single \\
        /data/lvmCFrame-00011017.fits --ml-model-file \\
        /models/sky.pt

    python -m \\
        skysub.sky_subtraction_block.cli batch --expnum 11017,11018 \\
        --ml-model-file /models/sky.pt --cpus 2

    python -m skysub.sky_subtraction_block.cli batch --input-list \\
        cframes.txt --ml-model-file /models/sky.pt --flat-output \\
        --output-dir /results

Single is the explicit compatible name for the original one-exposure operation:
python -m skysub.sky_subtraction_block INPUT [OPTIONS].
Batch writes a separate final product for each selected exposure. The trained model option is --ml-model-file; --ensemble remains a compatible alias.
Each command has its own --help with examples first and detailed options below.
"""
SINGLE_HELP = """Subtract sky from one observed-flux CFrame. Same scientific runner and defaults as the original command.

Examples:

    python -m skysub.sky_subtraction_block.cli single \\
        /data/lvmCFrame-00011017.fits --ml-model-file \\
        /models/sky.pt --nproc 1 --cpus 1 --output \\
        /results/lvmSFrame-skyml-block-00011017.fits

    python -m skysub.sky_subtraction_block.cli single \\
        /data/lvmCFrame-00011017.fits --ensemble /models/sky.pt \\
        --keep-preparation --output \\
        /results/reference/lvmSFrame-skyml-block-00011017.fits

    python -m skysub.sky_subtraction_block.cli single \\
        /data/lvmCFrame-00011017.fits --reuse-state \\
        /results/reference/lvmSFrame-skyml-preparation-00011017.fits \\
        --offset-roughness 0.01 --output \\
        /results/trial/lvmSFrame-skyml-block-00011017.fits

Ordinary runs calculate intermediate state in memory and write only the final product.
Saved-state reuse is an optional controlled-rerun mode. Existing --ensemble and --seed spellings remain supported.
"""
BATCH_HELP = """Subtract sky from selected CFrames, writing one final compact SFrame per exposure. Default: one exposure, one computational thread.

Examples:

    python -m \\
        skysub.sky_subtraction_block.cli batch --expnum 11017,11018 \\
        --ml-model-file /models/sky.pt --cpus 2 \\
        --output-dir /results --preserve-drp-tree

    python -m \\
        skysub.sky_subtraction_block.cli batch --mjd 60325:60330 \\
        --ml-model-file /models/sky.pt --max-exposures 3 --dry-run

    python -m skysub.sky_subtraction_block.cli batch --input-list \\
        cframes.txt --ml-model-file /models/sky.pt --flat-output \\
        --output-dir /results --nproc 2 --cpus 4

    python -m skysub.sky_subtraction_block.cli batch \\
        --redux-dir /data/redux/1.3.2 --input-list exposures.txt \\
        --ml-model-file /models/sky.pt --cpus 4

Commas/repeated values mean OR within a selector; different selectors mean AND. Ranges are inclusive.
Input-list files contain either CFrame paths or exposure IDs; do not mix them.
Set --cpus for the total CPU budget and --nproc for workers per exposure.
Concurrent exposures are computed automatically as floor(cpus / nproc).
Serial stages remain serial. Missing/identity-mismatched Gaia cache queries TAP; malformed FITS/schema is an error.
"""


class Existing(str, Enum):
    error = "error"
    skip = "skip"
    overwrite = "overwrite"


class Missing(str, Enum):
    error = "error"
    skip = "skip"


class OnError(str, Enum):
    continue_ = "continue"
    stop = "stop"


class Level(str, Enum):
    debug = "DEBUG"
    info = "INFO"
    warning = "WARNING"
    error = "ERROR"
    critical = "CRITICAL"


app = typer.Typer(name="skyml-block", help=ROOT_HELP, rich_markup_mode="markdown",
                  no_args_is_help=True, context_settings={"help_option_names": ["-h", "--help"]})

ModelFile = Annotated[Path | None, typer.Option("--ml-model-file", "--ensemble",
                      help=HELP["ensemble"], rich_help_panel="Model and Gaia")]
Gaia = Annotated[Path | None, typer.Option("--gaia-dir", help=HELP["gaia_dir"], rich_help_panel="Model and Gaia")]
Seed = Annotated[Path | None, typer.Option("--reuse-state", "--seed", help=HELP["seed"], rich_help_panel="Optional saved-state reuse")]
Keep = Annotated[bool, typer.Option("--keep-preparation", help=HELP["keep_preparation"], rich_help_panel="Optional saved-state reuse")]
Jobs = Annotated[int, typer.Option("--nproc", "--jobs", min=1, help=HELP["jobs"], rich_help_panel="CPU budget")]
CPUs = Annotated[int | None, typer.Option("--cpus", min=1, help=(
    "TOTAL allowed logical-CPU count for the whole command and descendants, including coordination. "
    "Default: batch uses nproc CPUs; single inherits allowed CPUs. Must cover nproc. "
    "Batch automatically runs floor(cpus / nproc) exposures concurrently."), rich_help_panel="CPU budget")]
CPUList = Annotated[str | None, typer.Option("--cpu-list", help=(
    "Optional exact CPU IDs, e.g. --cpus 4 --cpu-list 8-11 or 0,2,4,6. "
    "Supply --cpus; list must contain exactly that many distinct scheduler-allowed IDs."),
    rich_help_panel="CPU budget")]
Slit = Annotated[int, typer.Option(min=4, help=HELP["slit_count"], rich_help_panel="R/Z block LSF")]
Offset = Annotated[int, typer.Option(min=5, help=HELP["offset_count"], rich_help_panel="R/Z block LSF")]
Wave = Annotated[int, typer.Option(min=4, help=HELP["wave_count"], rich_help_panel="R/Z block LSF")]
OffsetPenalty = Annotated[float, typer.Option(min=0, help=HELP["offset_roughness"], rich_help_panel="R/Z regularization")]
WavePenalty = Annotated[float, typer.Option(min=0, help=HELP["wave_roughness"], rich_help_panel="R/Z regularization")]
SlitPenalty = Annotated[float, typer.Option(min=0, help=HELP["slit_roughness"], rich_help_panel="R/Z regularization")]
Prior = Annotated[float, typer.Option(min=0, help=HELP["reference_prior"], rich_help_panel="R/Z regularization")]
Centered = Annotated[bool, typer.Option("--centered", help=HELP["centered"], rich_help_panel="R/Z block LSF")]
LogLevel = Annotated[Level, typer.Option("--log-level", case_sensitive=False, help=(
    "DEBUG: detailed diagnostics; INFO: normal stages; WARNING: attention; "
    "ERROR: execution errors; CRITICAL: fatal application conditions."),
    rich_help_panel="Run control and logging")]


def _model(values):
    return {name: values[name] for name in ("ensemble", "gaia_dir", "jobs", *DEFAULTS)}


def _logging(level):
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")


def _version(value: bool):
    if value:
        revision = subprocess.run(["git", "-C", str(Path(__file__).resolve().parents[2]),
                                   "rev-parse", "--short", "HEAD"], text=True, capture_output=True)
        typer.echo("skyml-block CLI 1; revision " + (revision.stdout.strip() or "unavailable"))
        raise typer.Exit()


@app.callback()
def root(version: Annotated[bool, typer.Option("--version", callback=_version, is_eager=True,
                                              help="Print CLI version and code revision, then exit.")] = False):
    pass


@app.command(help=SINGLE_HELP)
def single(
    input: Annotated[Path, typer.Argument(exists=True, dir_okay=False, help="One local observed-flux lvmCFrame FITS.")],
    output: Annotated[Path | None, typer.Option(help=(
        "One final FITS file, retaining the lvmSFrame-skyml-block- prefix. Default: beside INPUT, "
        "replace lvmCFrame- with lvmSFrame-skyml-block-; remove trailing .gz. Existing files are protected by default."),
        rich_help_panel="Input and final output")] = None,
    ensemble: ModelFile = None, gaia_dir: Gaia = None, seed: Seed = None,
    on_existing: Annotated[Existing, typer.Option(help=(
        "error: destination conflict; skip: verify/reuse a matching product; "
        "overwrite: recompute and atomically replace the final FITS only on success."))] = Existing.error,
    jobs: Jobs = 8, cpus: CPUs = None, cpu_list: CPUList = None,
    slit_count: Slit = DEFAULTS["slit_count"], offset_count: Offset = DEFAULTS["offset_count"],
    wave_count: Wave = DEFAULTS["wave_count"],
    offset_roughness: OffsetPenalty = DEFAULTS["offset_roughness"],
    wave_roughness: WavePenalty = DEFAULTS["wave_roughness"],
    slit_roughness: SlitPenalty = DEFAULTS["slit_roughness"],
    reference_prior: Prior = DEFAULTS["reference_prior"],
    centered: Centered = False, keep_preparation: Keep = False,
    log_level: LogLevel = Level.info,
):
    model = _model(locals())
    _logging(log_level.value)
    try:
        validate_model(model)
        budget = cpu_budget(jobs=jobs, cpus=cpus, cpu_list=cpu_list, strict=False)
        apply_affinity(budget["cpu_ids"])
        from .planning import exposure_from_file, prepare_exposure
        entry = exposure_from_file(input)
        entry.output = output.expanduser().resolve() if output else entry.input_path.with_name(
            entry.input_path.name.replace("lvmCFrame-", "lvmSFrame-skyml-block-").removesuffix(".gz"))
        prepare_exposure(entry, model, seed=seed, on_existing=on_existing.value)
    except (ValueError, OSError, KeyError) as error:
        raise typer.BadParameter(str(error)) from error
    if entry.status == "reused":
        typer.echo(f"Reused: {entry.output}")
        return
    logging.info("Input: %s; %s: %s; output: %s; CPU budget: %s",
                 entry.input_path, entry.state_mode, entry.seed or "in-memory", entry.output, budget)
    try:
        from .resources import limit_numerics, runtime_resources
        limit_numerics()
        from .pipeline import run
        metadata = dict(fingerprint=entry.fingerprint, state_mode=entry.state_mode,
                        budget=budget, resources=runtime_resources())
        result = run(entry.input_path, output=entry.output, seed=entry.seed,
                     execution_metadata=metadata, overwrite=entry.overwrite, **model)
        typer.echo(str(result))
    except KeyboardInterrupt:
        raise typer.Exit(130) from None
    except Exception as error:
        logging.exception("Exposure %s: %s", entry.expnum, error)
        raise typer.Exit(1) from error


@app.command(help=BATCH_HELP)
def batch(
    inputs: Annotated[list[Path] | None, typer.Argument(help="Explicit CFrame paths; alternative to catalogue selection.")] = None,
    input_list: Annotated[Path | None, typer.Option(help=(
        "File containing CFrame paths (one per line, relative to this file) or exposure IDs "
        "(one per line, whitespace/comma-separated also accepted; requires --redux-dir or SAS environment). "
        "Do not mix IDs and paths. Blank lines/# comments allowed."),
                         rich_help_panel="Input and SDSS paths")] = None,
    sas_base_dir: Annotated[Path | None, typer.Option(help=(
        "SAS mirror root containing sdsswork/, e.g. /data/sas. Default: SAS_BASE_DIR environment. "
        "An explicit value also overrides SAS_BASE_DIR for scientific assets. Use --redux-dir for a partial data mirror."),
        rich_help_panel="Input and SDSS paths")] = None,
    cframe_root: Annotated[Path | None, typer.Option("--redux-dir", "--cframe-root", help=(
        "Input reduction directory containing drpall-VERSION.fits and TILEGROUP/TILE/MJD. "
        "The catalogue determines the DRP version; the directory name is unrestricted. "
        "Overrides SAS/environment data paths. --cframe-root is a compatible alias."),
        rich_help_panel="Input and SDSS paths")] = None,
    drpall: Annotated[Path | None, typer.Option(help=(
        "Explicit complete drpall-VERSION.fits catalogue. Default: discover exactly one "
        "drpall-*.fits in the input directory, excluding *-incomplete.fits. "
        "A relocated catalogue does not relocate an explicit --redux-dir."),
        rich_help_panel="Input and SDSS paths")] = None,
    expnum: Annotated[list[str] | None, typer.Option(help=(
        "Exposure numbers: 11017; 11017,11018; 11017:11025; or 11017,11020:11025. "
        "Repeatable; ranges inclusive. Default: no exposure-number cut."),
        rich_help_panel="Catalogue selection")] = None,
    mjd: Annotated[list[str] | None, typer.Option(help=(
        "Observing nights: 60325,60326 or 60325:60330. Repeatable, inclusive. Default: no night cut."),
        rich_help_panel="Catalogue selection")] = None,
    tileid: Annotated[list[str] | None, typer.Option(help=(
        "Tiles: 1029016,1029017 or open range 11112: (>=11112). Repeatable. Default: no tile cut."),
        rich_help_panel="Catalogue selection")] = None,
    exptime: Annotated[float | None, typer.Option(min=0, help=(
        "Match recorded exposure duration in seconds, e.g. 900; does not rescale flux. Default: no duration cut."),
        rich_help_panel="Catalogue selection")] = None,
    catalog_fluxcal: Annotated[str | None, typer.Option(help=(
        "Match drpall.fluxcal, e.g. MOD. Selects rows only, performs no flux calibration. "
        "Default: all recorded calibration modes."),
        rich_help_panel="Catalogue selection")] = None,
    all_exposures: Annotated[bool, typer.Option("--all-exposures", help=(
        "Select every exposure row in the chosen catalogue. Off by default. Use instead of filters."),
        rich_help_panel="Catalogue selection")] = False,
    max_exposures: Annotated[int | None, typer.Option(min=1, help=(
        "First N selected entries after filtering/deduplication and version/MJD/expnum/tile/path sorting. "
        "Not a quality ranking. Later skips are not replaced. Default: the whole selected set."),
        rich_help_panel="Catalogue selection")] = None,
    output_dir: Annotated[Path, typer.Option(help="Generated product root; relative to the launch directory.",
                         rich_help_panel="Final output")] = Path("skyml-block-output"),
    preserve_drp_tree: Annotated[bool, typer.Option("--preserve-drp-tree/--flat-output", help=(
        "Preserve: OUTPUT_DIR/TILEGROUP/TILE/MJD/lvmSFrame-skyml-block-EXPNUM.fits. "
        "Flat: write final FITS directly in OUTPUT_DIR. EXPNUM has eight digits."),
        rich_help_panel="Final output")] = True,
    on_existing: Annotated[Existing, typer.Option(help=(
        "error: destination conflict; skip: verify complete FITS, input/scientific fingerprint and checksums. "
        "Mismatches remain conflicts in skip mode. overwrite: recompute and atomically replace the final "
        "FITS only on success; retain the old output on failure. No automatic filename versions."),
        rich_help_panel="Final output")] = Existing.error,
    missing_input: Annotated[Missing, typer.Option(help="error: stop planning; skip: record unavailable local CFrame and continue.",
                         rich_help_panel="Catalogue selection")] = Missing.error,
    ensemble: ModelFile = None, gaia_dir: Gaia = None,
    reuse_state_dir: Annotated[Path | None, typer.Option(help=(
        "Advanced: existing snapshots from earlier --keep-preparation runs. Preserve-tree mode expects "
        "DIR/TILEGROUP/TILE/MJD/lvmSFrame-skyml-preparation-00011017.fits; flat mode "
        "expects that basename directly under DIR. Each executing exposure needs matching state. "
        "Unset: ordinary in-memory processing or compatible adjacent-state lookup."),
        rich_help_panel="Optional saved-state reuse")] = None,
    keep_preparation: Keep = False,
    jobs: Jobs = 1,
    cpus: CPUs = None, cpu_list: CPUList = None,
    slit_count: Slit = DEFAULTS["slit_count"], offset_count: Offset = DEFAULTS["offset_count"],
    wave_count: Wave = DEFAULTS["wave_count"],
    offset_roughness: OffsetPenalty = DEFAULTS["offset_roughness"],
    wave_roughness: WavePenalty = DEFAULTS["wave_roughness"],
    slit_roughness: SlitPenalty = DEFAULTS["slit_roughness"],
    reference_prior: Prior = DEFAULTS["reference_prior"], centered: Centered = False,
    run_dir: Annotated[Path | None, typer.Option(help="Plan/status/results/logs directory; default OUTPUT_DIR/.batch/<UTC run ID>.",
                         rich_help_panel="Run control and logging")] = None,
    on_error: Annotated[OnError, typer.Option(help="continue: schedule independent remaining work; stop: finish active work and stop submissions.",
                         rich_help_panel="Run control and logging")] = OnError.continue_,
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Print resolved selection, paths, state mode and CPU budget; no fitting, TAP queries or output creation.",
                         rich_help_panel="Run control and logging")] = False,
    log_level: LogLevel = Level.info,
    no_progress: Annotated[bool, typer.Option("--no-progress", help="Hide terminal progress; keep logs/status. Noninteractive output is plain.",
                         rich_help_panel="Run control and logging")] = False,
):
    model = _model(locals())
    _logging(log_level.value)
    if sas_base_dir is not None:
        os.environ["SAS_BASE_DIR"] = str(sas_base_dir.expanduser().resolve())
    try:
        validate_model(model)
        budget = cpu_budget(jobs=jobs, cpus=cpus, cpu_list=cpu_list)
        apply_affinity(budget["cpu_ids"])
        from .planning import make_plan, select_exposures
        entries, context = select_exposures(
            inputs=inputs, input_list=input_list, drpall=drpall, sas_base_dir=sas_base_dir,
            cframe_root=cframe_root, expnum=expnum, mjd=mjd,
            tileid=tileid, exptime=exptime, catalog_fluxcal=catalog_fluxcal,
            all_exposures=all_exposures, max_exposures=max_exposures, missing_input=missing_input.value)
        make_plan(entries, model, output_dir=output_dir, flat_output=not preserve_drp_tree,
                  reuse_state_dir=reuse_state_dir, on_existing=on_existing.value,
                  missing_input=missing_input.value)
        if dry_run or log_level == Level.debug:
            typer.echo("Resolved paths: " + str(context))
            typer.echo("CPU allocation: " + str(budget))
            for entry in entries:
                typer.echo(f"{entry.expnum:08d} [{entry.status}] {entry.input_path} -> {entry.output}; "
                           f"{entry.state_mode}: {entry.seed or 'in-memory'}")
        else:
            typer.echo(f"Selected {len(entries)} exposures; up to {budget['processes']} concurrent, "
                       f"{jobs} worker(s) per exposure")
            typer.echo(f"Output: {output_dir.expanduser().resolve()}")
        if dry_run:
            typer.echo("Dry run complete; output directories/files were not created.")
            return
        if any(entry.status == "queued" for entry in entries):
            check_runtime_dependencies()
        from .batch import run_batch
        result = run_batch(entries, model, budget, output_dir=output_dir, run_dir=run_dir,
                           context=context, on_error=on_error.value, log_level=log_level.value,
                           no_progress=no_progress)
    except KeyboardInterrupt:
        raise typer.Exit(130) from None
    except (ValueError, OSError, KeyError) as error:
        raise typer.BadParameter(str(error)) from error
    counts = result["counts"]
    summary = ", ".join(f"{label}={counts.get(status, 0)}" for status, label in (
        ("complete", "OK"), ("error", "failed"), ("reused", "reused"), ("skipped_input", "skipped"),
        ("not_started", "not started")))
    typer.echo(f"Batch {result['state']}: {summary}")
    if result["state"] != "complete":
        raise typer.Exit(1)


def main():
    app()


if __name__ == "__main__":
    main()
