"""Batch and Sky Subtraction share the selected physical-fiber CSV reference."""
from copy import deepcopy
import json
from types import SimpleNamespace
import shutil
import sys

import numpy as np
import pytest
from astropy.io import fits

from skysub import decompose_parallel as dp
from skysub.sky_decomp import lsf_reference as ref
from skysub.sky_subtraction import pipeline


def test_shared_mean_matches_previous_sky_subtraction_and_batch(monkeypatch):
    states = {i: ref.load_survey_ascii_lsf(ref.BUNDLED_CALIBRATION, i, allow_warning=True)
              for i in (84, 85, 500)}
    original = deepcopy(states)
    expected = deepcopy(states[84])
    for ch in "BRZ":
        values = [s.coefficients[ch] for s in states.values() if s.metrics[ch]["quality"] == "measured"]
        expected.coefficients[ch] = np.mean(values, axis=0)
    assert pipeline.mean_reference is ref.mean_reference
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE_FILE", ref.BUNDLED_CALIBRATION / "manifest.json")
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE", states)
    monkeypatch.setattr(dp, "_WORKER_META", np.array([(0,)], dtype=[("exposure", "i8")]))
    roles = np.zeros((1, 1944), dtype=np.uint8)
    roles[0, np.array([84, 85, 500]) - 1] = 1
    monkeypatch.setattr(dp, "_WORKER_FIBER_ROLE", roles)
    actual = dp._blue_reference_lsf("sci", 0)
    for ch in "BRZ":
        np.testing.assert_array_equal(actual.coefficients[ch], expected.coefficients[ch])
        np.testing.assert_array_equal(actual.coefficients[ch], ref.mean_reference(reversed(list(states.values()))).coefficients[ch])
    assert actual.metrics["B"]["n_fibers"] == 3
    assert actual.metrics["Z"]["n_fibers"] == 1
    for fiberid in states:
        for ch in "BRZ":
            np.testing.assert_array_equal(states[fiberid].coefficients[ch], original[fiberid].coefficients[ch])


def test_loader_grid_and_current_bundle_checksums(tmp_path, monkeypatch):
    directory = tmp_path / "lsf"
    shutil.copytree(ref.BUNDLED_CALIBRATION, directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    monkeypatch.setattr(ref, "_wave_fingerprint", lambda wave: tuple(manifest["wave_fingerprint"]))
    loaded = ref.load_fiber_references([500, 84, 500], np.arange(3), directory)
    assert set(loaded) == {84, 500}
    monkeypatch.setattr(ref, "_wave_fingerprint", lambda wave: (3, 0., 2., "wrong"))
    with pytest.raises(ValueError, match="Target grid"):
        ref.load_fiber_references([500], np.arange(3), directory)
    # A previous successful load must not hide a subsequently damaged CSV.
    with (directory / "survey_lsf_B.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="checksum"):
        ref.survey_calibration_manifest(directory)


def test_csv_membership_is_required_and_counts_are_checked():
    meta = fits.BinTableHDU(np.array([(2,)], dtype=[("fibers_sci_used", "i8")]), name="META")
    stack = fits.HDUList([fits.PrimaryHDU(), meta])
    with pytest.raises(ValueError, match="requires FIBER_ROLE"):
        dp._validate_csv_membership(stack)
    roles = np.zeros((1, 1944), dtype=np.uint8)
    roles[0, :4] = [1, 1, 2, 3]
    stack.append(fits.ImageHDU(roles, name="FIBER_ROLE"))
    dp._validate_csv_membership(stack)
    stack["FIBER_ROLE"].data[0, 1] = 0
    with pytest.raises(ValueError, match="counts differ"):
        dp._validate_csv_membership(stack)


@pytest.mark.parametrize("mutate", [False, True])
def test_worker_holds_explicit_blue_reference(monkeypatch, mutate):
    state = ref.load_survey_ascii_lsf(ref.BUNDLED_CALIBRATION, 500)
    fitted = deepcopy(state)
    if mutate:
        fitted.coefficients["B"][0, 0] += 0.01
    result = SimpleNamespace(fit_status="Solved", lsf_state=fitted)
    model = SimpleNamespace()
    monkeypatch.setattr(dp, "_science_line_mask_for_row", lambda idx: np.zeros(3, bool))
    monkeypatch.setattr(dp, "_telluric_decomposer", lambda *args: model)
    monkeypatch.setattr(dp, "_blue_reference_lsf", lambda *args: pytest.fail("explicit state must win"))
    monkeypatch.setattr(dp, "_install_split_zodi_amplitude_prior", lambda *args: None)
    monkeypatch.setattr(dp, "_fit_split_zodi_with_reversal_retry", lambda *args, **kw: (result, False, np.nan, False))
    monkeypatch.setattr(dp, "_reliability_flags", lambda *args, **kw: {})
    if mutate:
        with pytest.raises(RuntimeError, match="Blue LSF changed"):
            dp._fit_worker_row("sci", 0, np.ones(3), np.ones(3), initial_lsf_state=state)
    else:
        dp._fit_worker_row("sci", 0, np.ones(3), np.ones(3), initial_lsf_state=state)
        assert model.blue_reference_lsf is state
        assert fitted.metrics["B"]["reference_n_fibers"] == 1


def test_csv_provenance_and_fingerprint_change_with_manifest(tmp_path, monkeypatch):
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}")
    first = dp._fit_model_primary_meta("reference", manifest)
    manifest.write_text('{"changed": true}')
    second = dp._fit_model_primary_meta("reference", manifest)
    assert first["LSFSRC"] == "csv"
    assert first["LSFREFSH"] != second["LSFREFSH"]
    assert first["LSFMEAN"] == "selected_fiber_equal_mean"
    monkeypatch.setattr(dp, "file_sha256", lambda path: str(path))
    a = dp._compact_run_provenance(tmp_path / "stack.fits", np.arange(3.), tmp_path,
                                   {"blue_lsf_reference_sha256": first["LSFREFSH"][0]})
    b = dp._compact_run_provenance(tmp_path / "stack.fits", np.arange(3.), tmp_path,
                                   {"blue_lsf_reference_sha256": second["LSFREFSH"][0]})
    assert a["run_fingerprint"] != b["run_fingerprint"]


@pytest.mark.parametrize("compact", [False, True])
def test_empty_reference_is_a_failed_input_row(monkeypatch, compact):
    failed = SimpleNamespace(fit_status="failed_input")
    model = SimpleNamespace(failed_input_result=lambda reason: failed)
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE_FILE", "manifest.json")
    monkeypatch.setattr(dp, "_science_line_mask_for_row", lambda idx: np.zeros(3, bool))
    monkeypatch.setattr(dp, "_telluric_decomposer", lambda *args: model)
    def missing(*args):
        raise ValueError("No usable Z survey reference")
    monkeypatch.setattr(dp, "_blue_reference_lsf", missing)
    monkeypatch.setattr(dp, "_reliability_flags", lambda *args, **kw: {})
    result, flags = dp._fit_worker_row("sci", 0, np.ones(3), np.ones(3))
    assert result is failed
    monkeypatch.setattr(dp, "_WORKER_DECOMPOSER", object())
    monkeypatch.setattr(dp, "_WORKER_FLUX", {"sci": np.ones((1, 3))})
    monkeypatch.setattr(dp, "_WORKER_COMPACT_CACHE_DIR", "cache" if compact else None)
    monkeypatch.setattr(dp, "_WORKER_PROGRESS_QUEUE", None)
    monkeypatch.setattr(dp, "_compact_cache_is_current", lambda *args: False)
    monkeypatch.setattr(dp, "_fit_ivar_row", lambda *args: np.ones(3))
    saved = []
    monkeypatch.setattr(dp, "_save_compact_cache", lambda *args, **kw: saved.append(kw))
    kind, rows = dp.fit_chunk_worker(("sci", 0, 1))
    assert kind == "sci" and len(rows) == 1
    assert (saved[0]["result"] if compact else rows[0][1]) is failed


def test_worker_preloads_supported_ids_without_scanning_unprocessed_rows(monkeypatch, tmp_path):
    for name, value in list(vars(dp).items()):
        if name.startswith("_WORKER_"):
            monkeypatch.setattr(dp, name, value)
    state = ref.load_survey_ascii_lsf(ref.BUNDLED_CALIBRATION, 500)
    roles = np.zeros((2, 1944), dtype=np.uint8)
    roles[0, 499] = 1
    roles[1, 1943] = 1  # Unsupported in the test bundle; may lie beyond --limit.
    hdus = [fits.PrimaryHDU(), fits.ImageHDU(roles, name="FIBER_ROLE"),
            fits.BinTableHDU(np.array([(0,), (1,)], dtype=[("exposure", "i8")]), name="META")]
    for name in ("FLUX_SCI", "LSF_SCI"):
        hdus.append(fits.ImageHDU(np.ones((2, 3)), name=name))
    calls = []
    def preload(ids, wave, directory):
        calls.append(ids)
        assert ids is None
        return {500: state}
    monkeypatch.setattr(dp, "load_fiber_references", preload)
    class Ready(Exception):
        pass
    def stop_after_reference_loading():
        raise Ready
    monkeypatch.setitem(sys.modules, "lvmdrp.core.fluxcal", SimpleNamespace(
        TelluricCalculator=stop_after_reference_loading))
    with pytest.raises(Ready):
        dp.init_worker(np.arange(3.), .5, tmp_path, 1., fits.HDUList(hdus),
                       kinds=("sci",), lsf_reference_file=tmp_path / "manifest.json")
    assert calls == [None]
    np.testing.assert_array_equal(dp._blue_reference_lsf("sci", 0).coefficients["B"], state.coefficients["B"])
    with pytest.raises(KeyError):
        dp._blue_reference_lsf("sci", 1)
