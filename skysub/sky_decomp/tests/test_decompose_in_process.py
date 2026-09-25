"""The in-process single-exposure path must set a worker up exactly like `run`."""
import inspect

from skysub import decompose_parallel as dp


def _initargs_by_name(**overrides):
    kw = dict(
        wave="WAVE", lsf_sigma=0.5, base_dir="BASE", factor=1e14, data_file="F",
        progress_queue="Q", fit_model="M", n_refinement_cycles=5,
        worker_counter="C", pin_workers=False, diagnose_threads=False,
        palace_suffix="ps", palace_oh_suffix="pos", palace_diffuse_suffix="pds",
        exposure_seconds=900.0, moon_zodi_data_root="ROOT", n_spline_knots=11,
        n_zodi_spline_knots=1, zodi_smooth_lambda=0.1, moon_smooth_lambda=0.2,
        mask_science_lines=True, centre_on_halpha=True, fit_pixel_weights=True,
        fit_pixel_weight_clip=None, reversal_retry_bound=0.85,
        diffuse_ratio_bound_dex=0.2, diffuse_ratio_nominal=(1, 2, 3),
        diffuse_oh_centre_log10=-0.6, diffuse_oh_bound_dex=0.15,
        compact_cache_root="CACHE", run_fingerprint="FP")
    kw.update(overrides)
    names = [p for p in inspect.signature(dp.init_worker).parameters]
    return dict(zip(names, dp._worker_initargs(**kw))), kw


def test_initargs_land_on_the_right_init_worker_parameters():
    by_name, kw = _initargs_by_name()
    assert len(by_name) == len(inspect.signature(dp.init_worker).parameters)
    # Every value arrives under the parameter of the same meaning.
    assert by_name["pin_cpu"] is False and by_name["worker_counter"] == "C"
    assert by_name["moon_zodi_data_root"] == "ROOT"
    assert by_name["moon_smooth_lambda"] == 0.2 and by_name["zodi_smooth_lambda"] == 0.1
    assert by_name["reversal_retry_bound"] == 0.85
    assert by_name["diffuse_ratio_nominal"] == (1, 2, 3)
    assert by_name["compact_cache_dir"] == "CACHE" and by_name["run_fingerprint"] == "FP"


def test_parser_defaults_are_the_deployed_configuration():
    args = dp.build_arg_parser().parse_args(["stack.fits"])
    kw = dp._run_kwargs_from_args(args)
    assert kw["fit_model"] == dp.PALACECORR_VNF_SPLIT_ZODI_FIT_MODEL
    assert kw["mask_science_lines"] and kw["centre_on_halpha"] and kw["fit_pixel_weights"]
    assert kw["reversal_retry_bound"] == dp.SPLIT_ZODI_REVERSAL_RETRY_BOUND
    assert kw["n_spline_knots"] == dp.MOON_N_KNOTS_DEFAULT
    assert kw["diffuse_ratio_nominal"] == tuple(float(v) for v in dp.SPLIT_ZODI_DIFFUSE_RATIO_NOMINAL)


def test_an_in_memory_stack_is_passed_through_unchanged():
    from astropy.io import fits
    hdul = fits.HDUList([fits.PrimaryHDU()])
    by_name, _ = _initargs_by_name(data_file=hdul)
    assert by_name["data_file"] is hdul          # not stringified into a path
    by_name, _ = _initargs_by_name(data_file="stack.fits")
    assert by_name["data_file"] == "stack.fits"
