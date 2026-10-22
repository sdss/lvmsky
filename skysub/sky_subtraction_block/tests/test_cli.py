"""CLI selection, compatibility, reuse and resource contracts on small FITS fixtures."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from astropy.io import fits
from astropy.table import Table
import numpy as np
import pytest
from typer.testing import CliRunner

from skysub.sky_subtraction_block import cli, planning, resources
from skysub.sky_subtraction_block.__main__ import build_parser
from skysub.sky_subtraction_block.options import DEFAULTS

RUNNER = CliRunner()


@pytest.fixture
def cframe(tmp_path):
    def make(expnum=11017, mjd=60325, tileid=1029016, root=None):
        path = (root or tmp_path) / f"lvmCFrame-{expnum:08d}.fits"
        path.parent.mkdir(parents=True, exist_ok=True)
        header = fits.Header(dict(EXPOSURE=expnum, MJD=mjd, TILE_ID=tileid, DRPVER="1.3.2"))
        fits.HDUList([fits.PrimaryHDU(header=header),
                      fits.ImageHDU(np.arange(4.0), name="WAVE"),
                      fits.BinTableHDU(Table(dict(fiberid=[1, 2])), name="SLITMAP")]).writeto(path)
        return path
    return make


@pytest.fixture
def model(tmp_path):
    path = tmp_path / "trained-model.pt"
    path.write_bytes(b"model identity only; numerical models are not loaded in planning tests")
    return dict(ensemble=path, gaia_dir=None, jobs=1, **DEFAULTS)


@pytest.fixture
def catalogue(tmp_path, cframe):
    root = tmp_path / "redux/1.3.2"
    rows = []
    for exp, night, tile, duration, calibration in (
        (11019, 60326, 1029016, 900.0, "MOD"),
        (11018, 60325, 1028962, 600.0, "STD"),
        (11017, 60325, 1029016, 900.0, "MOD"),
    ):
        group = f"{tile // 1000:04d}XX"
        path = cframe(exp, night, tile, root / group / str(tile) / str(night))
        rows.append(dict(expnum=exp, mjd=night, tileid=tile, tilegrp=group, exptime=duration,
                         fluxcal=calibration, drpver="1.3.2",
                         location=f"sdsswork/lvm/spectro/redux/1.3.2/{group}/{tile}/{night}/lvmSFrame-{exp:08d}.fits"))
    root.mkdir(parents=True, exist_ok=True)
    cat = root / "drpall-1.3.2.fits"
    Table(rows=rows).write(cat)
    return root, cat


def model_arguments(model):
    return ["--ml-model-file", str(model["ensemble"])]


def test_legacy_defaults_and_aliases():
    old = vars(build_parser().parse_args(["lvmCFrame-00011017.fits", "--ensemble", "model.pt", "--seed", "state.fits"]))
    new = vars(build_parser().parse_args(["lvmCFrame-00011017.fits", "--ml-model-file", "model.pt", "--reuse-state", "state.fits"]))
    assert old == new
    assert old["jobs"] == 8 and old["seed"] == "state.fits"
    assert {name: old[name] for name in DEFAULTS} == DEFAULTS


@pytest.mark.parametrize("command,panel", [([], "Options"), (["single"], "Arguments"), (["batch"], "Arguments")])
def test_help_examples_before_options(command, panel):
    result = RUNNER.invoke(cli.app, [*command, "--help"])
    assert result.exit_code == 0, result.output
    assert result.output.index("Usage:") < result.output.index("Examples:") < result.output.index(panel)
    if command:
        assert "--ml-model-file" in result.output and "--ensemble" in result.output
    else:
        assert "single" in result.output and "batch" in result.output


def test_help_does_not_import_the_scientific_runtime():
    code = (
        "from typer.testing import CliRunner; from skysub.sky_subtraction_block.cli import app; "
        "import sys; r=CliRunner().invoke(app,['batch','--help']); "
        "assert r.exit_code==0,r.output; "
        "assert 'torch' not in sys.modules; "
        "assert 'skysub.sky_subtraction_block.pipeline' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.mark.parametrize("bad", ["1:0", ":", "1,,2", "one", "1:2:3", "-1"])
def test_bad_selectors(bad):
    with pytest.raises(ValueError):
        planning.parse_values([bad])


def test_inclusive_union_and_open_ranges():
    ranges = planning.parse_values(["11017,11020:11025", "11112:"])
    assert planning.includes(11017, ranges)
    assert planning.includes(11025, ranges)
    assert not planning.includes(11026, ranges)
    assert planning.includes(10**12, ranges)


def test_catalogue_and_intersection_sort_maximum(catalogue):
    root, cat = catalogue
    selected, context = planning.select_exposures(cframe_root=root, drpall=cat,
        expnum=["11017:11019"], mjd=["60325:60326"], catalog_fluxcal="MOD", max_exposures=1)
    assert [v.expnum for v in selected] == [11017]
    assert context["cframe_root"] == str(root)
    selected, _ = planning.select_exposures(cframe_root=root, drpall=cat, all_exposures=True)
    assert [v.expnum for v in selected] == [11017, 11018, 11019]


def test_selection_requires_explicit_intent(catalogue):
    root, cat = catalogue
    with pytest.raises(ValueError, match="Select"):
        planning.select_exposures(cframe_root=root, drpall=cat)
    with pytest.raises(ValueError, match="instead"):
        planning.select_exposures(cframe_root=root, all_exposures=True, mjd=["60325"])
    with pytest.raises(ValueError, match="absent"):
        planning.select_exposures(cframe_root=root, expnum=["11017:11020"])


def test_manifest_relative_paths_and_deduplication(tmp_path, cframe):
    path = cframe()
    manifest = tmp_path / "cframes.txt"
    manifest.write_text(f"# inputs\n{path.name}\n\n{path.name}\n")
    selected, _ = planning.select_exposures(input_list=manifest)
    assert len(selected) == 1 and selected[0].input_path == path


def test_explicit_and_catalogue_sources_are_exclusive(cframe):
    with pytest.raises(ValueError, match="not both"):
        planning.select_exposures(inputs=[cframe()], expnum=["11017"])


def test_sas_path_precedence(tmp_path, monkeypatch):
    monkeypatch.setenv("SAS_BASE_DIR", str(tmp_path / "sas"))
    monkeypatch.delenv("LVM_SPECTRO_REDUX", raising=False)
    root, cat = planning.paths()
    assert root == tmp_path / "sas/sdsswork/lvm/spectro/redux/1.3.2"
    monkeypatch.setenv("LVM_SPECTRO_REDUX", str(tmp_path / "redux"))
    assert planning.paths()[0] == tmp_path / "redux/1.3.2"
    assert planning.paths(sas_base_dir=tmp_path / "explicit")[0] == tmp_path / "explicit/sdsswork/lvm/spectro/redux/1.3.2"
    assert planning.paths(cframe_root=tmp_path / "direct")[0] == tmp_path / "direct"
    assert planning.paths(cframe_root=tmp_path / "direct", drpall=tmp_path / "moved.fits")[1] == tmp_path / "moved.fits"


def test_output_tree_flat_and_missing_policy(catalogue, model, tmp_path):
    root, _ = catalogue
    entries, _ = planning.select_exposures(cframe_root=root, expnum=["11017,11018"])
    planning.make_plan(entries, model, output_dir=tmp_path / "out")
    assert entries[0].output == tmp_path / "out/1.3.2/1029XX/1029016/60325/lvmSFrame-skyml-block-00011017.fits"
    assert entries[1].output == tmp_path / "out/1.3.2/1028XX/1028962/60325/lvmSFrame-skyml-block-00011018.fits"
    assert not (tmp_path / "out").exists()
    entries[1].input_path.unlink()
    planning.make_plan(entries, model, output_dir=tmp_path / "flat", flat_output=True, missing_input="skip")
    assert entries[0].output.parent == tmp_path / "flat"
    assert entries[1].status == "skipped_input"


def test_output_collision_rejected_before_computation(cframe, model, tmp_path):
    one = planning.exposure_from_file(cframe())
    two = deepcopy(one)
    two.version = "other-version"
    with pytest.raises(ValueError, match="same output"):
        planning.make_plan([one, two], model, output_dir=tmp_path / "flat", flat_output=True)


def test_dry_run_uses_new_and_compatible_model_flags(cframe, model, tmp_path, monkeypatch):
    path = cframe()
    monkeypatch.setattr(cli, "apply_affinity", lambda _: None)
    destination = tmp_path / "dry-output"
    for flag in ("--ml-model-file", "--ensemble"):
        result = RUNNER.invoke(cli.app, ["batch", str(path), flag, str(model["ensemble"]),
            "--flat-output", "--output-dir", str(destination), "--processes", "1", "--cpus", "1", "--dry-run"])
        assert result.exit_code == 0, result.output
        assert str(destination / "lvmSFrame-skyml-block-00011017.fits") in result.output
        assert "compute state in memory" in result.output
    assert not destination.exists()


@pytest.mark.parametrize("options", [dict(offset_count=6), dict(offset_roughness=float("nan")), dict(reference_prior=-1)])
def test_scientific_settings_validated_early(cframe, model, tmp_path, options):
    model.update(options)
    item = planning.exposure_from_file(cframe())
    with pytest.raises(ValueError):
        planning.make_plan([item], model, output_dir=tmp_path / "output")


def make_product(item):
    shape = (2, 4)
    hdus = [fits.PrimaryHDU(header=fits.Header(dict(SKYMETH="skyml-block")))]
    for name in ("FLUX", "IVAR", "MASK", "LSF", "SKY"):
        hdus.append(fits.ImageHDU(np.zeros(shape), name=name))
    hdus.extend([fits.ImageHDU(np.arange(4.0), name="WAVE"),
                 fits.BinTableHDU(Table(dict(fiberid=[1, 2])), name="SLITMAP"),
                 fits.ImageHDU(np.ones((1, 1)), name="BLOCK_LSF_COEF"),
                 fits.BinTableHDU(Table(dict(adopted=[True])), name="BLOCK_LSF_META")])
    from skysub.sky_subtraction.product import SHARED_ORDER
    present = {h.name for h in hdus}
    for name in (*SHARED_ORDER, "BLOCK_SLIT_KNOTS", "BLOCK_WAVE_KNOTS", "SKY_FIBER_COEF"):
        if name not in present and name not in {"GAIA_SOURCES", "FIT_LSF_HISTORY"}:
            hdus.append(fits.ImageHDU(np.zeros(1), name=name))
    payload = json.dumps(dict(cli_run=dict(fingerprint=item.fingerprint))).encode()
    hdus.append(fits.ImageHDU(np.frombuffer(payload, dtype="u1"), name="SKYML_META"))
    item.output.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(item.output, checksum=True)


def test_verified_existing_product_and_changed_settings(cframe, model, tmp_path):
    item = planning.exposure_from_file(cframe())
    planning.make_plan([item], model, output_dir=tmp_path / "output", flat_output=True)
    make_product(item)
    with pytest.raises(FileExistsError):
        planning.make_plan([item], model, output_dir=tmp_path / "output", flat_output=True)
    fresh = planning.exposure_from_file(item.input_path)
    planning.make_plan([fresh], model, output_dir=tmp_path / "output", flat_output=True, on_existing="skip")
    assert fresh.status == "reused"
    changed = dict(model, reference_prior=0.1)
    with pytest.raises(ValueError, match="fingerprint"):
        planning.make_plan([planning.exposure_from_file(item.input_path)], changed,
                           output_dir=tmp_path / "output", flat_output=True, on_existing="skip")


def test_cpu_budget_and_no_implicit_jobs(monkeypatch):
    monkeypatch.setattr(resources.os, "sched_getaffinity", lambda _: set(range(20)))
    allocation = resources.cpu_budget(processes=10, cpus=20)
    assert allocation["jobs"] == 1
    assert all(len(v) == 1 for v in allocation["subsets"])
    parallel = resources.cpu_budget(processes=10, jobs=2, cpus=20)
    assert all(len(v) == 2 for v in parallel["subsets"])
    with pytest.raises(ValueError, match="Need"):
        resources.cpu_budget(processes=10, jobs=2, cpus=10)


def test_cpu_ids_ranges_duplicates_and_available_set(monkeypatch):
    monkeypatch.setattr(resources.os, "sched_getaffinity", lambda _: {8, 9, 10, 11})
    assert resources.cpu_budget(2, 2, 4, "8-11")["cpu_ids"] == [8, 9, 10, 11]
    for bad in ("8,8", "0-3", "8-1000000000"):
        with pytest.raises(ValueError):
            resources.cpu_budget(1, 4, 4, bad)


def test_single_aliases_forward_same_scientific_options(cframe, model, tmp_path, monkeypatch):
    path = cframe()
    calls = []
    fake = SimpleNamespace(run=lambda input, **kwargs: calls.append((input, kwargs)) or kwargs["output"])
    monkeypatch.setitem(sys.modules, "skysub.sky_subtraction_block.pipeline", fake)
    monkeypatch.setattr(resources, "limit_numerics", lambda: None)
    monkeypatch.setattr(resources, "runtime_resources", lambda: {})
    monkeypatch.setattr(cli, "apply_affinity", lambda _: None)
    for flag in ("--ml-model-file", "--ensemble"):
        result = RUNNER.invoke(cli.app, ["single", str(path), flag, str(model["ensemble"]),
                                        "--output", str(tmp_path / "lvmSFrame-skyml-block-00011017.fits")])
        assert result.exit_code == 0, result.output
    assert calls[0] == calls[1]
    assert calls[0][1]["jobs"] == 8


def test_explicit_and_adjacent_saved_state(cframe, model, tmp_path):
    source = cframe()
    root = tmp_path / "states"
    root.mkdir()
    saved = root / "lvmSFrame-skyml-preparation-00011017.fits"
    payload = json.dumps(dict(input_sha256=planning.file_sha256(source))).encode()
    with fits.open(source) as original:
        hdus = [fits.PrimaryHDU(), original["WAVE"].copy(), original["SLITMAP"].copy()]
    for name in ("MED_COEF", "MED_LABELS", "ML_COEF", "MED_LSF_COEF", "RENDER_META", "LSF_FIBERS"):
        hdus.append(fits.ImageHDU(np.zeros(1), name=name))
    hdus.append(fits.ImageHDU(np.frombuffer(payload, dtype="u1"), name="SKYML_META"))
    fits.HDUList(hdus).writeto(saved)
    item = planning.exposure_from_file(source)
    item.output = root / "lvmSFrame-skyml-block-00011017.fits"
    planning.prepare_exposure(item, dict(model, ensemble=None))
    assert item.seed == saved and item.state_mode == "reuse adjacent saved state"
    explicit = planning.exposure_from_file(source)
    explicit.output = tmp_path / "other/lvmSFrame-skyml-block-00011017.fits"
    planning.prepare_exposure(explicit, dict(model, ensemble=None), seed=saved)
    assert explicit.state_mode == "reuse explicitly selected state"
    with source.open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="different CFrame"):
        planning.prepare_exposure(explicit, model, seed=saved)


def test_catalogue_identity_disagreement(catalogue, tmp_path, model):
    root, _ = catalogue
    entries, _ = planning.select_exposures(cframe_root=root, expnum=["11017"])
    with fits.open(entries[0].input_path, mode="update") as h:
        h[0].header["MJD"] = 60326
    with pytest.raises(ValueError, match="identity disagree"):
        planning.make_plan(entries, model, output_dir=tmp_path / "output")


@pytest.mark.parametrize("policy,expected", [("continue", ["error", "complete", "complete"]),
                                           ("stop", ["error", "not_started", "not_started"])])
def test_batch_error_policies_and_records(cframe, model, tmp_path, monkeypatch, policy, expected):
    from skysub.sky_subtraction_block import batch as scheduler

    entries = [planning.exposure_from_file(cframe(v)) for v in (11017, 11018, 11019)]
    for item in entries:
        item.output = tmp_path / f"product-{item.expnum}.fits"
    real_popen = subprocess.Popen
    probe = """
import json,sys,time
from pathlib import Path
request=json.loads(Path(sys.argv[1]).read_text())
entry=request['exposure']
time.sleep(.1)
code=1 if entry['expnum']==11017 else 0
if not code:
    Path(entry['output']).touch()
Path(request['result_path']).write_text(json.dumps(dict(
    state='error' if code else 'complete', error='probe error' if code else '',
    elapsed_seconds=.1)))
raise SystemExit(code)
"""
    def launch(command, **kwargs):
        return real_popen([sys.executable, "-c", probe, command[-1]], **kwargs)
    monkeypatch.setattr(scheduler.subprocess, "Popen", launch)
    monkeypatch.setattr(scheduler, "apply_affinity", lambda _: None)
    budget = dict(processes=1, jobs=1, cpus=1, cpu_ids=[0], subsets=[[0]])
    run_dir = tmp_path / "records"
    result = scheduler.run_batch(entries, model, budget, output_dir=tmp_path,
                                  run_dir=run_dir, on_error=policy, no_progress=True)
    assert [v.status for v in entries] == expected
    assert result["state"] == "complete_with_errors"
    assert json.loads((run_dir / "status.json").read_text())["counts"]["error"] == 1
    assert (run_dir / "plan.ecsv").is_file() and (run_dir / "results.ecsv").is_file()


def test_batch_interrupt_terminates_owned_process(cframe, model, tmp_path, monkeypatch):
    from skysub.sky_subtraction_block import batch as scheduler

    entry = planning.exposure_from_file(cframe())
    entry.output = tmp_path / "output.fits"
    real_popen = subprocess.Popen
    children = []
    def launch(command, **kwargs):
        process = real_popen([sys.executable, "-c", "import time;time.sleep(60)"], **kwargs)
        children.append(process)
        return process
    original_sleep = scheduler.time.sleep
    def interrupt(_):
        monkeypatch.setattr(scheduler.time, "sleep", original_sleep)
        raise KeyboardInterrupt
    monkeypatch.setattr(scheduler.subprocess, "Popen", launch)
    monkeypatch.setattr(scheduler, "apply_affinity", lambda _: None)
    monkeypatch.setattr(scheduler.time, "sleep", interrupt)
    budget = dict(processes=1, jobs=1, cpus=1, cpu_ids=[0], subsets=[[0]])
    run_dir = tmp_path / "interrupted"
    with pytest.raises(KeyboardInterrupt):
        scheduler.run_batch([entry], model, budget, output_dir=tmp_path,
                            run_dir=run_dir, no_progress=True)
    assert all(p.poll() is not None for p in children)
    assert json.loads((run_dir / "status.json").read_text())["state"] == "interrupted"

def test_legacy_long_option_abbreviations_remain_compatible():
    args = build_parser().parse_args(["lvmCFrame-00011017.fits", "--r", "0.25"])
    assert args.reference_prior == 0.25


def test_snapshot_saved_by_computation_can_be_resumed(cframe, model, tmp_path):
    source = cframe()
    item = planning.exposure_from_file(source)
    planning.make_plan([item], model, output_dir=tmp_path / "out", flat_output=True)
    # Create a same-input optional snapshot after the original compute plan.
    item.output.parent.mkdir()
    snapshot = item.output.parent / "lvmSFrame-skyml-preparation-00011017.fits"
    with fits.open(source) as original:
        hdus = [fits.PrimaryHDU(), original["WAVE"].copy(), original["SLITMAP"].copy()]
    for name in ("MED_COEF", "MED_LABELS", "ML_COEF", "MED_LSF_COEF", "RENDER_META", "LSF_FIBERS"):
        hdus.append(fits.ImageHDU(np.zeros(1), name=name))
    payload = json.dumps(dict(input_sha256=item.input_sha256)).encode()
    hdus.append(fits.ImageHDU(np.frombuffer(payload, dtype="u1"), name="SKYML_META"))
    fits.HDUList(hdus).writeto(snapshot)
    make_product(item)
    with fits.open(item.output, mode="update") as product:
        meta = dict(cli_run=dict(fingerprint=item.fingerprint, state_mode="compute state in memory"),
                    preparation_sha256=planning.file_sha256(snapshot))
        product["SKYML_META"].data = np.frombuffer(json.dumps(meta).encode(), dtype="u1")
        for hdu in product:
            hdu.add_checksum()
    fresh = planning.exposure_from_file(source)
    planning.make_plan([fresh], model, output_dir=item.output.parent, flat_output=True, on_existing="skip")
    assert fresh.status == "reused"


def test_checkout_launcher_help_outside_repository(tmp_path):
    import os
    root = Path(cli.__file__).resolve().parents[2]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(str(Path(v).resolve()) for v in env.get("PYTHONPATH", "").split(os.pathsep) if v)
    result = subprocess.run([sys.executable, str(root / "skyml-block"), "--help"],
                            cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "single" in result.stdout and "batch" in result.stdout

def test_runtime_resource_check_observes_excess_threads(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(
        get_num_threads=lambda: 1, get_num_interop_threads=lambda: 1))
    monkeypatch.setitem(sys.modules, "threadpoolctl", SimpleNamespace(
        threadpool_info=lambda: [dict(num_threads=2)]))
    with pytest.raises(RuntimeError, match="exceeds one thread"):
        resources.runtime_resources()


@pytest.mark.parametrize("manifest_mode", [False, True])
def test_explicit_missing_input_skip(cframe, model, tmp_path, manifest_mode):
    source = cframe()
    absent = tmp_path / "lvmCFrame-00011018.fits"
    arguments = [str(source), str(absent)]
    if manifest_mode:
        manifest = tmp_path / "inputs.txt"
        manifest.write_text(source.name + "\n" + absent.name + "\n")
        arguments = ["--input-list", str(manifest)]
    result = RUNNER.invoke(cli.app, ["batch", *arguments, *model_arguments(model),
        "--missing-input", "skip", "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert "skipped_input" in result.output and "queued" in result.output
    with pytest.raises(FileNotFoundError):
        planning.select_exposures(inputs=[source, absent])


def test_default_output_avoids_launcher(cframe, model, monkeypatch):
    from skysub.sky_subtraction_block import batch as scheduler
    captured = {}
    def run(entries, model, budget, **kwargs):
        captured.update(kwargs)
        return dict(state="complete", counts={}, run_dir="records")
    monkeypatch.setattr(scheduler, "run_batch", run)
    monkeypatch.setattr(cli, "apply_affinity", lambda _: None)
    result = RUNNER.invoke(cli.app, ["batch", str(cframe()), *model_arguments(model)])
    assert result.exit_code == 0, result.output
    assert captured["output_dir"] == Path("skyml-block-output")


@pytest.mark.parametrize("extension", ["ML_COEF", "LSF_COEF", "MED_LABELS", "BLOCK_WAVE_KNOTS"])
def test_reuse_requires_reconstruction_extensions(cframe, model, tmp_path, extension):
    item = planning.exposure_from_file(cframe())
    planning.make_plan([item], model, output_dir=tmp_path / "out")
    make_product(item)
    with fits.open(item.output) as product:
        reduced = fits.HDUList([h.copy() for h in product if h.name != extension])
    reduced.writeto(item.output, overwrite=True, checksum=True)
    with pytest.raises(ValueError, match="complete block product"):
        planning.verify_product(item)


def test_all_explicit_inputs_skipped_are_recorded(tmp_path, model, monkeypatch):
    from skysub.sky_subtraction_block import batch as scheduler
    path = tmp_path / 'lvmCFrame-00011018.fits'
    entries, _ = planning.select_exposures(inputs=[path], missing_input='skip')
    planning.make_plan(entries, model, output_dir=tmp_path / 'out', missing_input='skip')
    monkeypatch.setattr(scheduler, 'apply_affinity', lambda _: None)
    budget = dict(processes=1, jobs=1, cpus=1, cpu_ids=[0], subsets=[[0]])
    result = scheduler.run_batch(entries, model, budget, output_dir=tmp_path / 'out',
                                 no_progress=True)
    assert result['state'] == 'complete'
    assert result['counts'] == {'skipped_input': 1}
    records = Table.read(Path(result['run_dir']) / 'results.ecsv')
    assert str(records['input'][0]) == str(path)
    assert entries[0].output is None and entries[0].mjd == -1
