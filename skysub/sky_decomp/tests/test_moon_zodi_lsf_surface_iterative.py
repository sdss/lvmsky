"""Integration tests for the shared Moon/Zodi spline and FITS contract."""

from __future__ import annotations

import copy
from dataclasses import replace
from pathlib import Path

from astropy.io import fits
import numpy as np
import pytest

from skysub import decompose_parallel
from skysub.sky_decomp.result_io import load_moon_zodi_state
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.moon_zodi_lsf_surface_iterative import (
    MoonZodiLSFSurfaceIterativeResult,
    SkyDecompMoonZodiLSFSurfaceIterative,
)
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT, MoonZodiObservation


REFERENCE = Path(__file__).parent / "data" / "moon_zodi_predictor_reference_v1.npz"
BASE_DIR = Path(__file__).resolve().parents[2]
EXPECTED_COMPONENTS = ("moon", "zodi", "diffuse", "oh", "atom", "orc", "o2")
EXPECTED_HDU_ORDER = [
    "PRIMARY",
    "META",
    "COEF",
    "COEF_ERR",
    "BESTFIT",
    "BESTFIT_LSF",
    "FLUX_SIGMA_TOTAL",
    "RESID",
    "VECTOR_O2",
    "COMP_MOON",
    "COMP_ZODI",
    "COMP_DIFFUSE",
    "COMP_OH",
    "COMP_ATOM",
    "COMP_ORC",
    "COMP_O2",
    "LSF_COEF",
    "LSF_KNOTS",
    "LSF_META",
    "MZ_MODEL",
    "MZ_ASSETS",
    "MZ_KNOTS",
    "MZ_META",
]


@pytest.fixture(scope="module")
def fitted_case():
    with np.load(REFERENCE, allow_pickle=False) as reference:
        index = 1
        wave = np.asarray(reference["wave"], dtype=np.float64)
        moon = np.asarray(reference["moon"][index], dtype=np.float64)
        zodi = np.asarray(reference["zodi"][index], dtype=np.float64)
        lsf = np.asarray(reference["lsf"][index], dtype=np.float64)
        observation = MoonZodiObservation(
            int(reference["expnum"][index]),
            str(reference["date_obs"][index]),
            "sky_far",
            float(reference["ra_deg"][index]),
            float(reference["dec_deg"][index]),
            900.0,
            "assumed_900s",
        )
    decomposer = SkyDecompMoonZodiLSFSurfaceIterative(
        wave,
        lsf_sigma=0.5,
        moon_smooth_lambda=0.1,
        physical_to_fit_flux_scale=1.0e14,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=1),
    )
    assert decomposer.data_root == DEFAULT_DATA_ROOT.resolve()
    assert decomposer.base_dir == DEFAULT_DATA_ROOT.resolve()
    assert decomposer.pmd_dir == DEFAULT_DATA_ROOT.resolve() / "palace" / "PMD"
    assert decomposer.palace_oh_suffix == "_joint_v2_updated"
    assert (
        decomposer.palace_diffuse_suffix
        == "_joint_native_adam_invsky_p2_10000iter"
    )
    # A finite native-grid spectrum with a nonzero independent diffuse block.
    flux = moon + zodi + 0.02 * np.sum(decomposer.matrix_diffuse, axis=0)
    result = decomposer.fit(
        np.asarray(flux, dtype=np.float64),
        np.ones(wave.size, dtype=np.float64),
        observation=observation,
        detector_lsf_fwhm=lsf,
    )
    return wave, lsf, observation, decomposer, result


def test_shared_correction_and_component_closure(fitted_case):
    wave, _, _, decomposer, result = fitted_case
    assert isinstance(result, MoonZodiLSFSurfaceIterativeResult)
    assert tuple(result.components) == EXPECTED_COMPONENTS
    correction_indices = [
        index for index, name in enumerate(result.design_names) if name.startswith("MoonZodi_bs")
    ]
    assert correction_indices
    assert [result.design_names[index] for index in correction_indices] == [
        f"MoonZodi_bs{index:03d}" for index in range(len(correction_indices))
    ]
    assert np.all(result.coef[correction_indices] >= 0.0)
    assert decomposer.moon_smooth_lambda == 0.1
    assert result.lsf_state.completed_cycles == 1
    assert tuple(result.components) == EXPECTED_COMPONENTS
    total = sum(result.components.values(), np.zeros_like(wave))
    tolerance = 1.0e-10 * max(1.0, float(np.max(np.abs(result.bestfit_lsf))))
    assert np.max(np.abs(total - result.bestfit_lsf)) <= tolerance
    assert np.max(np.abs(result.components["moon"] + result.components["zodi"])) > 0.0
    o2_index = result.design_names.index("O2_b01")
    np.testing.assert_allclose(
        result.components["o2"],
        result.coef[o2_index] * result.vector_o2,
        rtol=0.0,
        atol=tolerance,
    )


def test_diffuse_block_is_not_multiplied_by_moon_zodi_spline(fitted_case):
    _, _, _, decomposer, result = fitted_case
    slices = decomposer._component_slices(
        decomposer._matrix_bundle(
            decomposer.matrix_oh,
            decomposer.matrix_moon,
            decomposer.matrix_diffuse,
            decomposer.matrix_atom,
            decomposer.matrix_orc,
            decomposer.matrix_o2,
        )
    )
    expected = decomposer.matrix_diffuse.T @ result.coef[slices["diffuse"]]
    np.testing.assert_allclose(result.components["diffuse"], expected, rtol=0.0, atol=1.0e-12)


def test_default_contract_requests_exactly_five_cycles(fitted_case):
    wave, _, _, _, _ = fitted_case
    decomposer = SkyDecompMoonZodiLSFSurfaceIterative(
        wave,
        physical_to_fit_flux_scale=1.0e14,
    )
    assert decomposer.config.n_refinement_cycles == 5


def test_fit_rejects_precision_grid_and_ivar_violations(fitted_case):
    wave, lsf, observation, decomposer, _ = fitted_case
    flux = np.ones_like(wave)
    with pytest.raises(ValueError, match="float64"):
        decomposer.fit(
            flux.astype(np.float32),
            np.ones_like(wave),
            observation=observation,
            detector_lsf_fwhm=lsf,
        )
    invalid = flux.copy()
    invalid[-1] = np.nan
    with pytest.raises(ValueError, match="zero ivar"):
        decomposer.fit(
            invalid,
            np.ones_like(wave),
            observation=observation,
            detector_lsf_fwhm=lsf,
        )


def test_fits_schema_roundtrip_thinning_and_named_o2(fitted_case, tmp_path):
    wave, _, _, _, result = fitted_case
    output = tmp_path / "moon_zodi.fits"
    decompose_parallel.results_to_fits([copy.deepcopy(result) for _ in range(3)], output)
    with fits.open(output) as hdul:
        assert [hdu.name for hdu in hdul] == EXPECTED_HDU_ORDER
        assert len(hdul["MZ_MODEL"].data) == 12
        assert len(hdul["MZ_ASSETS"].data) == 4
        assert hdul["MZ_META"].data.shape[0] == 3
        coef = hdul["COEF"].data
        assert "O2_b01" in coef.names
        np.testing.assert_allclose(
            hdul["COMP_O2"].data,
            np.asarray(coef["O2_b01"])[:, None] * hdul["VECTOR_O2"].data,
            rtol=0.0,
            atol=1.0e-10,
        )
        full = sum(
            (np.asarray(hdul[f"COMP_{name.upper()}"].data) for name in EXPECTED_COMPONENTS),
            np.zeros((3, wave.size), dtype=np.float64),
        )
        np.testing.assert_allclose(full, hdul["BESTFIT_LSF"].data, rtol=0.0, atol=1.0e-10)
    loaded = load_moon_zodi_state(output, spectrum_index=1)
    assert loaded == result.moon_zodi_state

    thinned = tmp_path / "moon_zodi_thinned.fits"
    decompose_parallel.thin_fits_every_n(output, thinned, 2)
    with fits.open(output) as original, fits.open(thinned) as reduced:
        assert len(reduced["META"].data) == 2
        assert len(reduced["MZ_META"].data) == 2
        np.testing.assert_array_equal(reduced["MZ_META"].data["spectrum_index"], [0, 1])
        assert len(reduced["LSF_META"].data) == 2 * 3
        np.testing.assert_array_equal(np.unique(reduced["LSF_META"].data["spectrum_index"]), [0, 1])
        for name in ("MZ_MODEL", "MZ_ASSETS", "MZ_KNOTS"):
            np.testing.assert_array_equal(reduced[name].data, original[name].data)


def test_writer_rejects_mixed_names_and_model_hashes(fitted_case, tmp_path):
    _, _, _, _, result = fitted_case
    with pytest.raises(ValueError, match="same concrete type"):
        decompose_parallel.results_to_fits([result, object()], tmp_path / "mixed.fits")
    changed_names = copy.deepcopy(result)
    changed_names.design_names = list(changed_names.design_names)
    changed_names.design_names[0] = "changed"
    with pytest.raises(ValueError, match="ordered design_names"):
        decompose_parallel.results_to_fits([result, changed_names], tmp_path / "names.fits")
    changed_model = copy.deepcopy(result)
    changed_model.moon_zodi_state = replace(
        changed_model.moon_zodi_state,
        checkpoint_sha256="0" * 64,
    )
    with pytest.raises(ValueError, match="incompatible model"):
        decompose_parallel.results_to_fits([result, changed_model], tmp_path / "model.fits")


def test_batch_role_coordinate_and_lsf_contract(monkeypatch):
    dtype = [
        ("expnum", "i8"),
        ("date_obs", "U30"),
        ("sci_ra", "f8"),
        ("sci_dec", "f8"),
        ("sky_near_ra", "f8"),
        ("sky_near_dec", "f8"),
        ("sky_far_ra", "f8"),
        ("sky_far_dec", "f8"),
    ]
    meta = np.zeros(1, dtype=dtype)
    meta[0] = (42, "2024-01-02T03:04:05", 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    monkeypatch.setattr(decompose_parallel, "_WORKER_META", meta)
    monkeypatch.setattr(decompose_parallel, "_WORKER_EXPOSURE_SECONDS", 900.0)
    expected = {
        "sci": ("sci", 1.0, 2.0),
        "sky1": ("sky_near", 3.0, 4.0),
        "sky2": ("sky_far", 5.0, 6.0),
    }
    for kind, values in expected.items():
        observation = decompose_parallel._moon_zodi_observation(kind, 0)
        assert (observation.role, observation.target_ra_deg, observation.target_dec_deg) == values
        assert observation.exposure_seconds == 900.0
        assert observation.exposure_seconds_source == "assumed_900s"
    assert decompose_parallel.FIT_MODEL_SUFFIXES == {
        "baseline": "",
        "lsf-surface-iterative": "_lsf_surface_iterative",
        "moon-zodi-lsf-surface-iterative": "_moon_zodi_lsf_surface_iterative",
    }


def test_runtime_data_roots_are_selected_by_fit_model(monkeypatch, tmp_path):
    validated = []

    def reject_legacy_path(_path):
        raise AssertionError("Moon/Zodi mode must not resolve the legacy PALACE path")

    monkeypatch.setattr(
        decompose_parallel,
        "validate_decomposition_data_root",
        validated.append,
    )
    monkeypatch.setattr(decompose_parallel, "resolve_base_dir", reject_legacy_path)

    named_root = tmp_path / "named_bundle"
    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        decompose_parallel.MOON_ZODI_FIT_MODEL,
        palace_dir=tmp_path / "invalid_legacy_path",
        moon_zodi_data_root=named_root,
    )
    assert base_dir == named_root.resolve()
    assert data_root == named_root.resolve()
    assert validated == [str(named_root.resolve())]

    positional_root = tmp_path / "positional_bundle"
    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        decompose_parallel.MOON_ZODI_FIT_MODEL,
        palace_dir=positional_root,
    )
    assert base_dir == positional_root.resolve()
    assert data_root == positional_root.resolve()

    default_root = tmp_path / "default_bundle"
    monkeypatch.setattr(
        decompose_parallel,
        "DEFAULT_MOON_ZODI_DATA_ROOT",
        default_root,
    )
    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        decompose_parallel.MOON_ZODI_FIT_MODEL,
    )
    assert base_dir == default_root.resolve()
    assert data_root == default_root.resolve()


def test_runtime_data_roots_preserve_legacy_resolution(monkeypatch, tmp_path):
    resolved = tmp_path / "legacy_root"
    calls = []

    def resolve(path):
        calls.append(path)
        return resolved

    monkeypatch.setattr(decompose_parallel, "resolve_base_dir", resolve)
    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        "lsf-surface-iterative",
        palace_dir="legacy-palace",
        moon_zodi_data_root=tmp_path / "ignored_moon_zodi_bundle",
    )
    assert base_dir == resolved
    assert data_root is None
    assert calls == ["legacy-palace"]
