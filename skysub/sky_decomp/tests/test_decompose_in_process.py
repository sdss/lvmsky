"""The in-process single-exposure path must set a worker up exactly like `run`."""
import inspect

from skysub import decompose_parallel as dp


def _initargs_by_name(**overrides):
    kw = dict(
        wave="WAVE", lsf_sigma=0.5, base_dir="BASE", factor=1e14, data_file="F",
        progress_queue="Q", n_refinement_cycles=5,
        worker_counter="C", pin_workers=False, diagnose_threads=False,
        exposure_seconds=900.0, compact_cache_root="CACHE", run_fingerprint="FP")
    kw.update(overrides)
    names = [p for p in inspect.signature(dp.init_worker).parameters]
    return dict(zip(names, dp._worker_initargs(**kw))), kw


def test_initargs_land_on_the_right_init_worker_parameters():
    by_name, kw = _initargs_by_name()
    assert len(by_name) == len(inspect.signature(dp.init_worker).parameters)
    # Every value arrives under the parameter of the same meaning.
    assert by_name["pin_cpu"] is False and by_name["worker_counter"] == "C"
    assert by_name["n_refinement_cycles"] == 5 and by_name["exposure_seconds"] == 900.0
    assert by_name["compact_cache_dir"] == "CACHE" and by_name["run_fingerprint"] == "FP"
    assert by_name["kinds"] == ("sci", "sky1", "sky2") and by_name["lsf_reference_file"] is None


def test_parser_defaults_are_the_deployed_configuration():
    args = dp.build_arg_parser().parse_args(["stack.fits"])
    kw = dp._run_kwargs_from_args(args)
    assert kw["fit_model"] == dp.FIT_MODEL
    assert kw["n_refinement_cycles"] == dp.N_REFINEMENT_CYCLES
    assert kw["blue_lsf"] == "solar"


def test_individual_fiber_cli_and_input_pixel_weights(monkeypatch):
    import numpy as np

    args = dp.build_arg_parser().parse_args(
        ["stack.fits", "--compact-only", "--kinds", "sci"])
    assert dp._run_kwargs_from_args(args)["kinds"] == ("sci",)
    monkeypatch.setattr(dp, "_WORKER_INPUT_IVAR", {"sci": np.array([[4.0, 9.0, np.nan]])})
    monkeypatch.setattr(dp, "_WORKER_INPUT_MASK", {"sci": np.array([[0, 1, 0]])})
    monkeypatch.setattr(dp, "_WORKER_FACTOR", 2.0)
    got = dp._fit_ivar_row("sci", 0, np.array([1.0, 1.0, 1.0]))
    assert np.array_equal(got, [1.0, 0.0, 0.0])


def test_mixed_individual_fibers_use_their_physical_pointings(monkeypatch):
    import numpy as np

    dtype = [(name, "U16") for name in
             ("telescope", "sky_near_label")]
    meta = np.array([("Sci", "SkyE"), ("SkyE", "SkyE"),
                     ("SkyW", "SkyE")], dtype=dtype)
    monkeypatch.setattr(dp, "_WORKER_META", meta)
    assert dp._row_role("sci", 0) == ("sci", "sci_ra", "sci_dec", "sci")
    assert dp._row_role("sci", 1) == ("sky_near", "skye_ra", "skye_dec", "skye")
    assert dp._row_role("sci", 2) == ("sky_far", "skyw_ra", "skyw_dec", "skyw")
    monkeypatch.setattr(dp, "_WORKER_SCIENCE_LINE_MASK", np.array([True, False]))
    assert not dp._science_line_mask_for_row(1).any()


def test_stacked_sky_arms_point_at_their_telescope(monkeypatch):
    import numpy as np

    meta = np.array([("SkyW", "SkyE")], dtype=[("sky_near_label", "U8"), ("sky_far_label", "U8")])
    monkeypatch.setattr(dp, "_WORKER_META", meta)
    # The pointing names the telescope, so the airmass column is skyw_/skye_airmass.
    assert dp._row_role("sky1", 0) == ("sky_near", "sky_near_ra", "sky_near_dec", "skyw")
    assert dp._row_role("sky2", 0) == ("sky_far", "sky_far_ra", "sky_far_dec", "skye")


def test_blue_lsf_options_and_reference_file(tmp_path):
    import pytest

    kw = dp._run_kwargs_from_args(dp.build_arg_parser().parse_args(["stack.fits"]))
    assert kw["blue_lsf"] == "solar"
    assert kw["lsf_calibration"] == dp.DEFAULT_LSF_CALIBRATION_DIR
    kw = dp._run_kwargs_from_args(dp.build_arg_parser().parse_args(
        ["stack.fits", "--blue-lsf", "reference", "--lsf-calibration", str(tmp_path)]))
    assert kw["blue_lsf"] == "reference" and kw["lsf_calibration"] == tmp_path
    assert dp._survey_lsf_reference_file("solar", tmp_path) is None
    with pytest.raises(FileNotFoundError):
        dp._survey_lsf_reference_file("reference", tmp_path)
    survey = tmp_path / dp.SURVEY_LSF_REFERENCE
    survey.parent.mkdir()
    survey.write_bytes(b"")
    assert dp._survey_lsf_reference_file("reference", tmp_path) == survey.resolve()
    by_name, _ = _initargs_by_name(lsf_reference_file=survey)
    assert by_name["lsf_reference_file"] == str(survey)


def test_blue_reference_follows_the_stacked_telescope_or_the_fiber(monkeypatch):
    import numpy as np

    means = {"Sci": "SCI", "SkyE": "SKYE", "SkyW": "SKYW"}
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE", dict(means))
    stacked = np.array([("SkyE", "SkyW"), ("skyw", "skye")],
                       dtype=[("sky_near_label", "U8"), ("sky_far_label", "U8")])
    monkeypatch.setattr(dp, "_WORKER_META", stacked)
    assert dp._blue_reference_lsf("sci", 0) == "SCI"
    assert dp._blue_reference_lsf("sky1", 0) == "SKYE"
    assert dp._blue_reference_lsf("sky2", 0) == "SKYW"
    assert dp._blue_reference_lsf("sky1", 1) == "SKYW"
    individual = np.array([(937, "SkyW", "SkyE", "SkyW")], dtype=[
        ("fiberid", "i8"), ("telescope", "U8"), ("sky_near_label", "U8"), ("sky_far_label", "U8")])
    monkeypatch.setattr(dp, "_WORKER_META", individual)
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE", {**means, 937: "FIBER937"})
    assert dp._blue_reference_lsf("sci", 0) == "FIBER937"


def test_blue_reference_averages_the_fiber_role_fibers(monkeypatch):
    import numpy as np
    import pytest

    class FakeReferences:
        fiberids = np.array([1, 2, 3, 5])

        def mean(self, fiberids, label=""):
            return ("MEAN", tuple(int(f) for f in fiberids), label)

    means = {"Sci": "SCI", "SkyE": "SKYE", "SkyW": "SKYW"}
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE", dict(means))
    monkeypatch.setattr(dp, "_WORKER_BLUE_REFERENCE_SET", FakeReferences())
    meta = np.array([("SkyE", "SkyW")] * 2,
                    dtype=[("sky_near_label", "U8"), ("sky_far_label", "U8")])
    monkeypatch.setattr(dp, "_WORKER_META", meta)
    # Column j is fiberid j + 1; 1 = sci, 2 = near, 3 = far.
    roles = np.array([[1, 1, 2, 0, 3, 0], [1, 0, 0, 2, 0, 0]], dtype=np.uint8)
    monkeypatch.setattr(dp, "_WORKER_FIBER_ROLE", roles)
    assert dp._blue_reference_lsf("sci", 0) == ("MEAN", (1, 2), "sci row 0")
    assert dp._blue_reference_lsf("sky1", 0) == ("MEAN", (3,), "sky1 row 0")
    assert dp._blue_reference_lsf("sky2", 0) == ("MEAN", (5,), "sky2 row 0")
    # Row 1: the near-arm fiber (4) has no reference, the far arm has no fibers.
    with pytest.warns(RuntimeWarning):
        assert dp._blue_reference_lsf("sky1", 1) == "SKYE"
    with pytest.warns(RuntimeWarning):
        assert dp._blue_reference_lsf("sky2", 1) == "SKYW"
    monkeypatch.setattr(dp, "_WORKER_FIBER_ROLE", None)
    assert dp._blue_reference_lsf("sky1", 0) == "SKYE"


def test_an_in_memory_stack_is_passed_through_unchanged():
    from astropy.io import fits
    hdul = fits.HDUList([fits.PrimaryHDU()])
    by_name, _ = _initargs_by_name(data_file=hdul)
    assert by_name["data_file"] is hdul          # not stringified into a path
    by_name, _ = _initargs_by_name(data_file="stack.fits")
    assert by_name["data_file"] == "stack.fits"


def test_solarblue_cli_routing_and_provenance():
    args = dp.build_arg_parser().parse_args(["stack.fits", "--fit-model", dp.FIT_MODEL])
    kw = dp._run_kwargs_from_args(args)
    assert kw["fit_model"] == dp.FIT_MODEL
    assert dp.PALACE_OH_SUFFIX == dp.SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    assert dp.FIT_MODEL_SUFFIX.endswith("_solarblue")
    meta = dp._fit_model_primary_meta()
    assert meta["BLUELSF"] == "solar_continuum_joint"
    assert meta["OHRIDGE"] == dp.SKYFAR_LINEAR_RIDGE_LAMBDA


def test_solarblue_cache_fingerprint_tracks_implementation(tmp_path, monkeypatch):
    import numpy as np
    monkeypatch.setattr(dp, "file_sha256", lambda p: str(p))
    before = dp._compact_run_provenance(tmp_path / "stack.fits", np.arange(10.),
        tmp_path, {})
    assert "solar_lsf.py" in before["source_sha256"]
    monkeypatch.setattr(dp, "file_sha256", lambda p: "changed" if str(p).endswith("solar_lsf.py") else str(p))
    after = dp._compact_run_provenance(tmp_path / "stack.fits", np.arange(10.),
        tmp_path, {})
    assert before["run_fingerprint"] != after["run_fingerprint"]
