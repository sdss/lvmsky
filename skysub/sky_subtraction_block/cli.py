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
from .resources import apply_affinity, cpu_budget, limit_environment

limit_environment()

ROOT_HELP = """Subtract sky from calibrated CFrame exposures with the sky-only block LSF method.

Examples:

    python -m skysub.sky_subtraction_block.cli single \\
        /data/lvmCFrame-00011017.fits --ml-model-file \\
        /models/sky.pt

    SAS_BASE_DIR=/data/sas python -m \\
        skysub.sky_subtraction_block.cli batch --expnum 11017,11018 \\
        --ml-model-file /models/sky.pt --processes 2 --cpus 2

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
        /models/sky.pt --jobs 1 --cpus 1 --output \\
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

    SAS_BASE_DIR=/data/sas python -m \\
        skysub.sky_subtraction_block.cli batch --expnum 11017,11018 \\
        --ml-model-file /models/sky.pt --processes 2 --cpus 2 \\
        --output-dir /results --preserve-drp-tree

    SAS_BASE_DIR=/data/sas python -m \\
        skysub.sky_subtraction_block.cli batch --mjd 60325:60330 \\
        --ml-model-file /models/sky.pt --max-exposures 3 --dry-run

    python -m skysub.sky_subtraction_block.cli batch --input-list \\
        cframes.txt --ml-model-file /models/sky.pt --flat-output \\
        --output-dir /results --processes 2 --jobs 2 --cpus 4

    SAS_BASE_DIR=/data/sas python -m \\
        skysub.sky_subtraction_block.cli batch --all-exposures \\
        --ml-model-file /models/sky.pt --dry-run

Commas/repeated values mean OR within a selector; different selectors mean AND. Ranges are inclusive.
Total CPUs never increase internal workers automatically: processes=10, cpus=20 still uses jobs=1.
Request jobs=2 explicitly to permit two CPUs per exposure; require processes * jobs <= cpus.
Serial stages remain serial. Missing/identity-mismatched Gaia cache queries TAP; malformed FITS/schema is an error.
"""


class Existing(str, Enum):
    error = "error"
    skip = "skip"


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
Jobs = Annotated[int, typer.Option(min=1, help=HELP["jobs"], rich_help_panel="CPU budget")]
CPUs = Annotated[int | None, typer.Option("--cpus", min=1, help=(
    "TOTAL allowed logical-CPU count for the whole command and descendants, including coordination. "
    "Default: batch processes * jobs; single inherits allowed CPUs. Must cover processes * jobs. "
    "Extra CPUs never increase jobs automatically."), rich_help_panel="CPU budget")]
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
        "replace lvmCFrame- with lvmSFrame-skyml-block-; remove trailing .gz. Existing files are protected."),
        rich_help_panel="Input and final output")] = None,
    ensemble: ModelFile = None, gaia_dir: Gaia = None, seed: Seed = None,
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
        prepare_exposure(entry, model, seed=seed)
    except (ValueError, OSError, KeyError) as error:
        raise typer.BadParameter(str(error)) from error
    logging.info("Input: %s; %s: %s; output: %s; CPU budget: %s",
                 entry.input_path, entry.state_mode, entry.seed or "in-memory", entry.output, budget)
    try:
        from .resources import limit_numerics, runtime_resources
        limit_numerics()
        from .pipeline import run
        metadata = dict(fingerprint=entry.fingerprint, state_mode=entry.state_mode,
                        budget=budget, resources=runtime_resources())
        result = run(entry.input_path, output=entry.output, seed=entry.seed,
                     execution_metadata=metadata, **model)
        typer.echo(str(result))
    except KeyboardInterrupt:
        raise typer.Exit(130) from None
    except Exception as error:
        logging.exception("Exposure %s: %s", entry.expnum, error)
        raise typer.Exit(1) from error


@app.command(help=BATCH_HELP)
def batch(
    inputs: Annotated[list[Path] | None, typer.Argument(help="Explicit CFrame paths; alternative to catalogue selection.")] = None,
    input_list: Annotated[Path | None, typer.Option(help="One CFrame path per line; relative to the list file. Blank lines/# comments allowed.",
                         rich_help_panel="Input and SDSS paths")] = None,
    sas_base_dir: Annotated[Path | None, typer.Option(help=(
        "SAS mirror root containing sdsswork/, e.g. /data/sas. Default: SAS_BASE_DIR environment. "
        "An explicit value also overrides SAS_BASE_DIR for scientific assets. Use --cframe-root for a partial data mirror."),
        rich_help_panel="Input and SDSS paths")] = None,
    drp_version: Annotated[str, typer.Option(help="Chosen catalogue/data version; current supported native grid is DRP 1.3.2.",
                         rich_help_panel="Input and SDSS paths")] = "1.3.2",
    cframe_root: Annotated[Path | None, typer.Option(help=(
        "Advanced: version directory containing TILEGROUP/TILE/MJD, e.g. /scratch/lvm/1.3.2. "
        "Precedence: this option; explicit --sas-base-dir plus standard suffix; "
        "LVM_SPECTRO_REDUX/version; SAS_BASE_DIR/sdsswork/lvm/spectro/redux/version."),
        rich_help_panel="Input and SDSS paths")] = None,
    drpall: Annotated[Path | None, typer.Option(help=(
        "LVM summary FITS. Default in catalogue mode: CFRAME_ROOT/drpall-VERSION.fits. "
        "A relocated catalogue does not relocate the CFrame root."),
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
        "Preserve: OUTPUT_DIR/VERSION/TILEGROUP/TILE/MJD/lvmSFrame-skyml-block-EXPNUM.fits. "
        "Flat: write final FITS directly in OUTPUT_DIR. EXPNUM has eight digits."),
        rich_help_panel="Final output")] = True,
    on_existing: Annotated[Existing, typer.Option(help=(
        "error: destination conflict; skip: verify complete FITS, input/scientific fingerprint and checksums. "
        "Mismatches remain conflicts. No automatic filename versions."),
        rich_help_panel="Final output")] = Existing.error,
    missing_input: Annotated[Missing, typer.Option(help="error: stop planning; skip: record unavailable local CFrame and continue.",
                         rich_help_panel="Catalogue selection")] = Missing.error,
    ensemble: ModelFile = None, gaia_dir: Gaia = None,
    reuse_state_dir: Annotated[Path | None, typer.Option(help=(
        "Advanced: existing snapshots from earlier --keep-preparation runs. Preserve-tree mode expects "
        "DIR/VERSION/TILEGROUP/TILE/MJD/lvmSFrame-skyml-preparation-00011017.fits; flat mode "
        "expects that basename directly under DIR. Each executing exposure needs matching state. "
        "Unset: ordinary in-memory processing or compatible adjacent-state lookup."),
        rich_help_panel="Optional saved-state reuse")] = None,
    keep_preparation: Keep = False,
    jobs: Jobs = 1,
    processes: Annotated[int, typer.Option(min=1, help="Maximum exposures active simultaneously, separate from internal --jobs.",
                         rich_help_panel="CPU budget")] = 1,
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
        budget = cpu_budget(processes, jobs, cpus, cpu_list)
        apply_affinity(budget["cpu_ids"])
        from .planning import make_plan, select_exposures
        entries, context = select_exposures(
            inputs=inputs, input_list=input_list, drpall=drpall, sas_base_dir=sas_base_dir,
            cframe_root=cframe_root, drp_version=drp_version, expnum=expnum, mjd=mjd,
            tileid=tileid, exptime=exptime, catalog_fluxcal=catalog_fluxcal,
            all_exposures=all_exposures, max_exposures=max_exposures, missing_input=missing_input.value)
        make_plan(entries, model, output_dir=output_dir, flat_output=not preserve_drp_tree,
                  reuse_state_dir=reuse_state_dir, on_existing=on_existing.value,
                  missing_input=missing_input.value)
        typer.echo("Resolved paths: " + str(context))
        typer.echo("CPU allocation: " + str(budget))
        for entry in entries:
            typer.echo(f"{entry.expnum:08d} [{entry.status}] {entry.input_path} -> {entry.output}; "
                       f"{entry.state_mode}: {entry.seed or 'in-memory'}")
        if dry_run:
            typer.echo("Dry run complete; output directories/files were not created.")
            return
        from .batch import run_batch
        result = run_batch(entries, model, budget, output_dir=output_dir, run_dir=run_dir,
                           context=context, on_error=on_error.value, log_level=log_level.value,
                           no_progress=no_progress)
    except KeyboardInterrupt:
        raise typer.Exit(130) from None
    except (ValueError, OSError, KeyError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(f"Batch {result['state']}: {result['counts']}; records: {result['run_dir']}")
    if result["state"] != "complete":
        raise typer.Exit(1)


def main():
    app()


if __name__ == "__main__":
    main()
