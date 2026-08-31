from __future__ import annotations

from dataclasses import asdict, fields
from hashlib import sha256
import json
from numbers import Real
import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from skysub.decompose_parallel import results_to_fits
from skysub.sky_decomp.fit import SkyDecompResult
from skysub.sky_decomp.lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    LSFSurfaceIterativeResult,
    SkyDecompLSFSurfaceIterative,
    evaluate_lsf_surface,
    load_lsf_surface_state,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
BASE_DIR = REPO_ROOT / "skysub"
DATA_DIR = Path(__file__).with_name("data")
MANIFEST_PATH = DATA_DIR / "lsf_surface_iterative_row837_n5_golden.json"
FIXTURE_PATH = DATA_DIR / "lsf_surface_iterative_row837_n5_golden.npz"
COMPONENT_ORDER = (
    "oh",
    "moon",
    "ho2",
    "feo",
    "o2ac",
    "atom",
    "orc",
    "o2",
    "diffuse",
)
LSF_META_COLUMNS = (
    "spectrum_index",
    "channel",
    "available",
    "lower",
    "upper",
    "degree",
    "n_basis",
    "n_knots",
    "tap_offsets",
    "status",
    "reason",
    "schema_version",
    "requested_cycles",
    "completed_cycles",
    "fit_status",
    "failure_reason",
    "final_continuum_status",
    "final_line_status",
    "knot_strategy",
    "legacy_kernel_representation",
    "wave_n",
    "wave_min",
    "wave_max",
    "wave_sha256",
    "line_weight",
    "skyline_cumulative_fraction",
    "skyline_half_width_angstrom",
    "huber_transition_sigma",
    "config_n_basis",
    "config_degree",
    "roughness_fraction",
    "fallback_prior_fraction",
    "information_prior_max_boost",
    "background_degree",
    "blue_fit_lower",
    "n_pixels",
    "n_valid_pixels",
    "chi2_red",
    "rms_resid",
    "center_pix_min",
    "center_pix_median",
    "center_pix_max",
    "sigma_pix_min",
    "sigma_pix_median",
    "sigma_pix_max",
    "runtime_sec",
    "fit_lower",
    "fit_upper",
    "knots_fixed",
    "model",
)


def _golden_input_unavailable(message):
    if os.environ.get("SKYSUB_REQUIRE_LSF_GOLDEN") == "1":
        pytest.fail(message)
    pytest.skip(message)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stable_json(value):
    if isinstance(value, dict):
        return {
            str(key): _stable_json(item)
            for key, item in value.items()
            if "runtime" not in str(key).lower() and "elapsed" not in str(key).lower()
        }
    if isinstance(value, (list, tuple)):
        return [_stable_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return _stable_json(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        if np.isnan(value):
            return "nan"
        return "inf" if value > 0.0 else "-inf"
    return value


def _assert_json_close(actual, expected, *, rtol, atol, path="root"):
    actual = _stable_json(actual)
    expected = _stable_json(expected)
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        assert set(actual) == set(expected), path
        for key in expected:
            _assert_json_close(
                actual[key],
                expected[key],
                rtol=rtol,
                atol=atol,
                path=f"{path}.{key}",
            )
        return
    if isinstance(expected, list):
        assert isinstance(actual, list), path
        assert len(actual) == len(expected), path
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected, strict=True)):
            _assert_json_close(
                actual_item,
                expected_item,
                rtol=rtol,
                atol=atol,
                path=f"{path}[{index}]",
            )
        return
    if (
        isinstance(actual, Real)
        and not isinstance(actual, bool)
        and isinstance(expected, Real)
        and not isinstance(expected, bool)
    ):
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=rtol,
            atol=atol,
            equal_nan=True,
            err_msg=path,
        )
        return
    assert actual == expected, path


def _actual_arrays(wave, flux, ivar, model, result):
    state = result.lsf_state
    arrays = {
        "wave": wave,
        "flux": flux,
        "ivar": ivar,
        "coef": result.coef,
        "design_names": np.asarray(result.design_names),
        "bestfit": result.bestfit,
        "resid": result.resid,
        "bestfit_lsf": result.bestfit_lsf,
        "skyline_mask": model.skyline_mask,
        "continuum_weights": model.continuum_weights,
        "final_continuum": model.final_continuum,
        "final_line_model": model.final_line_model,
        "lsf_surface": evaluate_lsf_surface(state, wave),
        "tap_offsets": state.tap_offsets,
        "moon_knots": result.moon_knots,
        "moon_boosted_pixels": result.moon_boosted_pixels,
    }
    arrays.update({f"component_{name}": values for name, values in result.components.items()})
    for channel in ("B", "R", "Z"):
        arrays[f"lsf_coefficient_{channel}"] = state.coefficients[channel]
        arrays[f"lsf_knots_{channel}"] = state.knot_vectors[channel]
        arrays[f"lsf_kernel_{channel}"] = result.lsf_kernels[channel]
    return arrays


@pytest.fixture(scope="module")
def golden_fit():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 1
    assert _file_sha256(FIXTURE_PATH) == manifest["fixture"]["sha256"]

    with np.load(FIXTURE_PATH, allow_pickle=False) as stored:
        assert set(stored.files) == set(manifest["array_keys"])
        golden = {name: stored[name] for name in stored.files}

    for name, specification in manifest["fixture"]["arrays"].items():
        assert list(golden[name].shape) == specification["shape"]
        assert str(golden[name].dtype) == specification["dtype"]

    for relative_path, expected_hash in manifest["source_hashes"].items():
        input_path = BASE_DIR / relative_path
        if not input_path.is_file():
            _golden_input_unavailable(f"golden input is unavailable: {input_path}")
        if _file_sha256(input_path) != expected_hash:
            _golden_input_unavailable(
                f"golden input provenance does not match: {input_path}"
            )

    input_specification = manifest["input"]
    data_path = BASE_DIR / input_specification["data_path"]
    with fits.open(data_path, memmap=True) as hdul:
        assert [hdu.name for hdu in hdul] == input_specification["input_hdus"]
        wave = np.asarray(
            hdul[input_specification["wave_extension"]].data,
            dtype=np.float64,
        ).copy()
        flux = np.asarray(
            hdul[input_specification["flux_hdu"]].data[input_specification["row"]],
            dtype=np.float64,
        ).copy()
    flux *= input_specification["flux_scale"]
    ivar = np.ones_like(wave, dtype=np.float64)

    assert wave.shape == (input_specification["wave_n"],)
    assert np.array_equal(wave, golden["wave"])
    assert np.array_equal(flux, golden["flux"])
    assert np.array_equal(ivar, golden["ivar"])

    config = LSFSurfaceIterativeConfig()
    assert asdict(config) == manifest["model"]["config"]
    assert config.n_refinement_cycles == 5
    model = SkyDecompLSFSurfaceIterative(
        wave,
        base_dir=BASE_DIR / manifest["model"]["base_dir"],
        config=config,
        **manifest["model"]["kwargs"],
    )
    result = model.fit(flux, ivar)
    return manifest, golden, wave, flux, ivar, model, result


def test_default_n5_fit_matches_real_data_golden(golden_fit):
    manifest, golden, wave, flux, ivar, model, result = golden_fit
    rtol = manifest["comparison"]["rtol"]
    atol = manifest["comparison"]["atol"]
    actual = _actual_arrays(wave, flux, ivar, model, result)

    assert set(actual) == set(manifest["array_keys"])
    for name in manifest["array_keys"]:
        assert actual[name].shape == golden[name].shape, name
        assert actual[name].dtype == golden[name].dtype, name
        if name in {"wave", "flux", "ivar", "design_names", "skyline_mask", "tap_offsets"}:
            assert np.array_equal(actual[name], golden[name]), name
        elif name == "coef":
            # The constrained decomposition is rank-degenerate: solver-equivalent
            # coefficient vectors can differ substantially while every persisted
            # physical component and the full model remain unchanged. Those
            # scientifically identifiable arrays are compared below.
            assert np.all(np.isfinite(actual[name]))
            assert np.all(actual[name] >= -1.0e-8)
        else:
            np.testing.assert_allclose(
                actual[name],
                golden[name],
                rtol=rtol,
                atol=atol,
                equal_nan=True,
                err_msg=name,
            )

    baseline_fields = [field.name for field in fields(SkyDecompResult)]
    iterative_fields = [field.name for field in fields(LSFSurfaceIterativeResult)]
    assert baseline_fields == manifest["public_contract"]["baseline_result_fields"]
    assert iterative_fields == manifest["public_contract"]["iterative_result_fields"]
    assert result.design_names == manifest["public_contract"]["design_names"]

    expected_result = manifest["result"]
    assert result.fit_status == expected_result["fit_status"]
    assert list(result.components) == list(COMPONENT_ORDER)
    assert sorted(result.components) == expected_result["component_names"]
    assert list(result.lsf_kernels) == expected_result["lsf_kernel_names"]
    for name in (
        "reduced_chi2",
        "r2",
        "rms_resid",
        "resid_level",
        "t_o2",
        "t_o2_err",
        "o2_valid_frac",
    ):
        np.testing.assert_allclose(
            getattr(result, name),
            expected_result[name],
            rtol=rtol,
            atol=atol,
            equal_nan=True,
            err_msg=name,
        )
    assert result.o2_fit_status == expected_result["o2_fit_status"]
    _assert_json_close(
        result.lsf_metrics,
        expected_result["lsf_metrics"],
        rtol=rtol,
        atol=atol,
        path="result.lsf_metrics",
    )

    state = result.lsf_state
    expected_state = manifest["lsf_state"]
    state_scalars = {
        name: getattr(state, name)
        for name in (
            "degrees",
            "channel_bounds",
            "config",
            "metrics",
            "requested_cycles",
            "completed_cycles",
            "wave_n",
            "wave_min",
            "wave_max",
            "wave_sha256",
            "fit_status",
            "failure_reason",
            "final_continuum_status",
            "final_line_status",
            "knot_strategy",
            "legacy_kernel_representation",
            "schema_version",
        )
    }
    _assert_json_close(
        state_scalars,
        expected_state,
        rtol=rtol,
        atol=atol,
        path="lsf_state",
    )
    assert state.requested_cycles == 5
    assert state.completed_cycles == 5
    assert not hasattr(state, "wave")
    assert not hasattr(state, "kernel_surface")
    _assert_json_close(
        model.channel_noise,
        manifest["runtime_state"]["channel_noise"],
        rtol=rtol,
        atol=atol,
        path="model.channel_noise",
    )
    _assert_json_close(
        model.stage_metrics,
        manifest["runtime_state"]["stage_metrics"],
        rtol=rtol,
        atol=atol,
        path="model.stage_metrics",
    )


def test_default_n5_state_has_semantic_fits_round_trip(golden_fit, tmp_path):
    manifest, _, wave, _, _, _, result = golden_fit
    rtol = manifest["comparison"]["rtol"]
    atol = manifest["comparison"]["atol"]
    output_path = tmp_path / "lsf_surface_iterative_row837_n5.fits"
    results_to_fits([result], output_path)

    expected_hdus = [
        "PRIMARY",
        "META",
        "COEF",
        "COEF_ERR",
        "BESTFIT",
        "BESTFIT_LSF",
        "FLUX_SIGMA_TOTAL",
        "RESID",
        "VECTOR_O2",
        "COEF_COV_MOON",
        *[f"COMP_{name.upper()}" for name in COMPONENT_ORDER],
        "LSF_COEF",
        "LSF_KNOTS",
        "LSF_META",
    ]
    with fits.open(output_path, memmap=False) as hdul:
        hdu_names = [hdu.name for hdu in hdul]
        assert hdu_names == expected_hdus
        np.testing.assert_allclose(
            np.asarray(hdul["COEF_COV_MOON"].data, dtype=float),
            np.asarray(result.coef_cov_moon, dtype=float)[None, ...],
            rtol=rtol,
            atol=atol,
            equal_nan=True,
        )
        assert "COEF_COV_ZODI" not in hdu_names
        assert tuple(hdul["LSF_META"].columns.names) == LSF_META_COLUMNS
        assert "LSF_WAVE" not in hdu_names
        assert not any("SURFACE" in name for name in hdu_names)

    original = result.lsf_state
    restored = load_lsf_surface_state(output_path)
    for name in (
        "degrees",
        "channel_bounds",
        "config",
        "requested_cycles",
        "completed_cycles",
        "wave_n",
        "wave_min",
        "wave_max",
        "wave_sha256",
        "fit_status",
        "failure_reason",
        "final_continuum_status",
        "final_line_status",
        "knot_strategy",
        "legacy_kernel_representation",
        "schema_version",
    ):
        _assert_json_close(
            getattr(restored, name),
            getattr(original, name),
            rtol=rtol,
            atol=atol,
            path=f"restored.{name}",
        )
    for channel, restored_metrics in restored.metrics.items():
        persisted_original_metrics = {
            name: original.metrics[channel][name] for name in restored_metrics
        }
        _assert_json_close(
            restored_metrics,
            persisted_original_metrics,
            rtol=rtol,
            atol=atol,
            path=f"restored.metrics.{channel}",
        )
    assert np.array_equal(restored.tap_offsets, original.tap_offsets)
    for channel in ("B", "R", "Z"):
        np.testing.assert_allclose(
            restored.coefficients[channel],
            original.coefficients[channel],
            rtol=rtol,
            atol=atol,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            restored.knot_vectors[channel],
            original.knot_vectors[channel],
            rtol=rtol,
            atol=atol,
            equal_nan=True,
        )
    np.testing.assert_allclose(
        evaluate_lsf_surface(restored, wave),
        evaluate_lsf_surface(original, wave),
        rtol=rtol,
        atol=atol,
        equal_nan=True,
    )
    assert not hasattr(restored, "wave")
    assert not hasattr(restored, "kernel_surface")
