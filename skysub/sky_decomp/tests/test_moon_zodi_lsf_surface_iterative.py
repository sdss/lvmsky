"""Integration tests for the shared Moon/Zodi spline and FITS contract."""

from __future__ import annotations

import copy
from dataclasses import replace
from pathlib import Path

from astropy.io import fits
import numpy as np
import pytest

from skysub import decompose_parallel
from skysub.sky_decomp import moon_zodi_model
from skysub.sky_decomp.result_io import (
    INVALID_OBSERVATION_FIT_STATUS,
    load_lsf_surface_state,
    load_moon_zodi_state,
)
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.moon_zodi_lsf_surface_iterative import (
    MoonZodiLSFSurfaceIterativeResult,
    SkyDecompMoonZodiLSFSurfaceIterative,
)
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_OH_SUFFIX,
    SKYFAR_LINEAR_RIDGE_LAMBDA,
    SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX,
    MoonZodiObservation,
)


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
    "COEF_COV_MOON",
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


def test_science_mask_lsf_reference_preserves_arm_median_weighting():
    wave = np.array([6562.8, 7000.0])
    lsf = {
        "sci": np.array([[1.0, 10.0], [100.0, 10.0]]),
        "sky1": np.array([[2.0, 20.0], [2.0, 20.0]]),
        "sky2": np.array([[3.0, 30.0], [3.0, 30.0]]),
    }

    reference = decompose_parallel._science_mask_lsf_reference(wave, lsf)

    assert reference[0] == 3.0
    assert np.isnan(reference[1])


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
    # Defaults combine the telluric-aware OH table with the hybrid continuum.
    assert (
        decomposer.palace_oh_suffix
        == DEFAULT_PALACE_OH_SUFFIX
        == "_telluric_upper_parity_lsf_adam_25000_v1"
    )
    # The hybrid table uses canonical PALACE
    # fcHO2/fcFeO with the native-LVM fcO2Ac.  Asserted against the module
    # constant AND the literal so a silent default change still trips here.
    assert (
        decomposer.palace_diffuse_suffix
        == moon_zodi_model.DEFAULT_PALACE_DIFFUSE_SUFFIX
        == "_canonhyb_v1"
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


@pytest.fixture(scope="module")
def invalid_case(fitted_case):
    wave, lsf, observation, decomposer, _ = fitted_case
    invalid_observation = replace(
        observation,
        target_ra_deg=(observation.target_ra_deg + 180.0) % 360.0,
        target_dec_deg=-observation.target_dec_deg,
    )
    result = decomposer.fit(
        np.zeros_like(wave),
        np.ones_like(wave),
        observation=invalid_observation,
        detector_lsf_fwhm=lsf,
    )
    return invalid_observation, result


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


def test_invalid_observation_result_is_same_schema_and_all_nan(
    fitted_case,
    invalid_case,
):
    _, _, _, decomposer, _ = fitted_case
    observation, result = invalid_case
    assert isinstance(result, MoonZodiLSFSurfaceIterativeResult)
    assert result.fit_status == INVALID_OBSERVATION_FIT_STATUS
    assert "reason=target_below_horizon" in result.fit_summary
    assert tuple(result.components) == EXPECTED_COMPONENTS
    for value in (
        result.coef,
        result.coef_err,
        result.bestfit,
        result.bestfit_lsf,
        result.bestfit_lsf_sigma,
        result.resid,
        result.vector_o2,
        *result.components.values(),
    ):
        assert np.all(np.isnan(value))
    assert all(
        np.all(np.isnan(coefficient))
        for coefficient in result.lsf_state.coefficients.values()
    )
    assert result.moon_zodi_state.observation == observation
    assert result.moon_zodi_state.geometry.target_altitude_deg <= 0.0
    assert np.isnan(result.moon_zodi_state.geometry.target_airmass)
    assert "invalid_observation" in result.moon_zodi_state.flags
    assert decomposer._prediction_state is result.moon_zodi_state
    assert decomposer.fit_status == INVALID_OBSERVATION_FIT_STATUS
    assert np.all(np.isnan(decomposer.bestfit_lsf))


def test_valid_fit_recovers_after_invalid_placeholder(fitted_case, invalid_case):
    wave, lsf, observation, decomposer, _ = fitted_case
    _ = invalid_case
    recovered = decomposer.fit(
        np.zeros_like(wave),
        np.ones_like(wave),
        observation=observation,
        detector_lsf_fwhm=lsf,
    )
    assert recovered.fit_status != INVALID_OBSERVATION_FIT_STATUS
    assert np.all(np.isfinite(recovered.bestfit_lsf))


def test_near_sun_leinert_rejection_returns_nan_result(fitted_case, monkeypatch):
    wave, lsf, observation, decomposer, _ = fitted_case

    def reject_near_sun(*_args):
        raise moon_zodi_model._LeinertDomainError(
            "zodi_invalid_near_sun_cell",
            "Invalid near-Sun Leinert cell: lon=6.770, lat=0.000",
        )

    monkeypatch.setattr(moon_zodi_model, "_interpolate_leinert", reject_near_sun)
    result = decomposer.fit(
        np.zeros_like(wave),
        np.ones_like(wave),
        observation=observation,
        detector_lsf_fwhm=lsf,
    )
    assert result.fit_status == INVALID_OBSERVATION_FIT_STATUS
    assert "reason=zodi_invalid_near_sun_cell" in result.fit_summary
    assert np.all(np.isnan(result.bestfit_lsf))
    assert np.isnan(result.moon_zodi_state.geometry.zodi_b500)


def test_decomposition_zeroes_only_moon_below_horizon(fitted_case):
    _, original_lsf, original_observation, decomposer, _ = fitted_case
    with np.load(REFERENCE, allow_pickle=False) as reference:
        index = 0
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
    decomposer._install_prediction(observation, lsf)
    assert decomposer._prediction_state.geometry.moon_altitude_deg < 0.0
    assert np.all(decomposer.physical_moon_prediction == 0.0)
    assert np.max(decomposer.physical_zodi_prediction) > 0.0
    assert "moon_zeroed_for_decomposition" in decomposer._prediction_state.flags

    # Leave the shared module-scoped fixture in its original valid state.
    decomposer._install_prediction(original_observation, original_lsf)


def test_success_and_invalid_rows_share_one_fits_contract(
    fitted_case,
    invalid_case,
    tmp_path,
):
    wave, _, _, _, valid = fitted_case
    _, invalid = invalid_case
    output = tmp_path / "mixed_validity.fits"
    decompose_parallel.results_to_fits([valid, invalid], output)
    with fits.open(output) as hdul:
        assert [hdu.name for hdu in hdul] == EXPECTED_HDU_ORDER
        assert len(hdul["META"].data) == 2
        assert len(hdul["MZ_META"].data) == 2
        assert hdul["META"].data[1]["fit_status"].strip() == (
            INVALID_OBSERVATION_FIT_STATUS
        )
        for name in (
            "BESTFIT",
            "BESTFIT_LSF",
            "FLUX_SIGMA_TOTAL",
            "RESID",
            "VECTOR_O2",
            *(f"COMP_{component.upper()}" for component in EXPECTED_COMPONENTS),
        ):
            assert np.all(np.isnan(hdul[name].data[1]))
        assert np.all(np.isnan(np.asarray(list(hdul["COEF"].data[1]))))
        assert np.all(np.isnan(hdul["LSF_COEF"].data[1]))
        assert hdul["MZ_META"].data[1]["target_altitude_deg"] <= 0.0
        assert np.isnan(hdul["MZ_META"].data[1]["target_airmass"])
        assert hdul["BESTFIT_LSF"].data.shape == (2, wave.size)
    loaded_mz = load_moon_zodi_state(output, spectrum_index=1)
    loaded_lsf = load_lsf_surface_state(output, spectrum_index=1)
    assert "invalid_observation" in loaded_mz.flags
    assert loaded_lsf.fit_status == INVALID_OBSERVATION_FIT_STATUS
    assert all(
        np.all(np.isnan(coefficient))
        for coefficient in loaded_lsf.coefficients.values()
    )


def test_writer_rejects_non_nan_invalid_placeholder(invalid_case, tmp_path):
    _, invalid = invalid_case
    changed = copy.deepcopy(invalid)
    changed.components["zodi"][0] = 0.0
    with pytest.raises(ValueError, match="only NaN in component zodi"):
        decompose_parallel.results_to_fits(
            [changed],
            tmp_path / "invalid_placeholder.fits",
        )


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
        "lsf-surface-iterative-split-zodi": "_lsf_surface_iterative_split_zodi",
        "lsf-spline2d-split-zodi": "_lsf_spline2d_split_zodi",
        "moon-zodi-lsf-surface-iterative": "_moon_zodi_lsf_surface_iterative",
        "adam25k-telluric-lsf-spline2d": "_adam25k_telluric_lsf_spline2d",
        "palace-aijc-vnf-line-amplitude-pca30": (
            "_palace_aijc_vnf_line_amplitude_pca30"
        ),
        "adam25k-telluric-split-zodi-lsf-spline2d": (
            "_adam25k_telluric_split_zodi_lsf_spline2d"
        ),
        "palace-aijc-vnf-split-zodi-lsf-spline2d": (
            "_palace_aijc_vnf_split_zodi_lsf_spline2d"
        ),
        "palacecorr-aijc-vnf-split-zodi-lsf-spline2d": (
            "_palacecorr_aijc_vnf_split_zodi_lsf_spline2d"
        ),
        "palace-aijc-vnf-pca30-split-zodi-lsf-spline2d": (
            "_palace_aijc_vnf_pca30_split_zodi_lsf_spline2d"
        ),
        "adam25k-telluric-niv-continuum": "_adam25k_telluric_niv_continuum",
        "palace-aijc-vnf-pca30-niv-continuum": (
            "_palace_aijc_vnf_pca30_niv_continuum"
        ),
    }
    assert (
        decompose_parallel.ADAM25K_NIV_CONTINUUM_FIT_MODEL
        == decompose_parallel.LEGACY_ADAM25K_SPLIT_ZODI_FIT_MODEL
    )
    assert (
        decompose_parallel.PALACE_VNF_PCA30_NIV_CONTINUUM_FIT_MODEL
        == decompose_parallel.LEGACY_PALACE_VNF_PCA30_SPLIT_ZODI_FIT_MODEL
    )


def test_palacecorr_suffix_and_primary_provenance_are_explicit():
    fit_model = decompose_parallel.PALACECORR_VNF_SPLIT_ZODI_FIT_MODEL
    suffix = decompose_parallel._resolved_palace_oh_suffix(fit_model, None)
    assert suffix == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    assert decompose_parallel._fit_model_primary_meta(fit_model, suffix) == {
        "DECOMPM": fit_model,
        "OHFILE": f"pmd_popmodel_OH{suffix}.dat",
        "OHRIDGE": SKYFAR_LINEAR_RIDGE_LAMBDA,
    }
    with pytest.raises(ValueError, match="palacecorr requires"):
        decompose_parallel._resolved_palace_oh_suffix(fit_model, "_other")


@pytest.mark.parametrize(
    ("fit_model", "class_path"),
    [
        (
            decompose_parallel.ADAM25K_TELLURIC_FIT_MODEL,
            "skysub.sky_decomp.telluric_corrected_lines."
            "SkyDecompAdam25kTelluricLSFSpline2D",
        ),
        (
            decompose_parallel.PALACE_VNF_PCA30_FIT_MODEL,
            "skysub.sky_decomp.residual_pca."
            "SkyDecompPalaceAijcVNFLineAmplitudePCA",
        ),
        (
            decompose_parallel.ADAM25K_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.telluric_corrected_lines."
            "SkyDecompAdam25kTelluricSplitZodiLSFSpline2D",
        ),
        (
            decompose_parallel.PALACE_VNF_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.residual_pca."
            "SkyDecompPalaceAijcVNFSplitZodiLSFSpline2D",
        ),
        (
            decompose_parallel.PALACECORR_VNF_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.residual_pca."
            "SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D",
        ),
        (
            decompose_parallel.PALACE_VNF_PCA30_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.residual_pca."
            "SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30",
        ),
        (
            decompose_parallel.LEGACY_ADAM25K_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.telluric_corrected_lines."
            "SkyDecompAdam25kTelluricSplitZodiLSFSpline2D",
        ),
        (
            decompose_parallel.LEGACY_PALACE_VNF_PCA30_SPLIT_ZODI_FIT_MODEL,
            "skysub.sky_decomp.residual_pca."
            "SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30",
        ),
    ],
)
@pytest.mark.parametrize(("pwv_med", "expected_pwv"), [(4.2, 4.2), (-999.9, 15.0)])
def test_telluric_cli_models_use_role_lsf_pwv_and_airmass(
    monkeypatch,
    tmp_path,
    fit_model,
    class_path,
    pwv_med,
    expected_pwv,
):
    from astropy.table import Table
    import lvmdrp.core.fluxcal

    wave = np.array([5000.0, 5001.0, 5002.0])
    flux = np.array([[1.0, 2.0, 3.0]])
    meta = Table(
        {
            "pwv_med": [pwv_med],
            "sci_airmass": [1.3],
            "skye_airmass": [1.4],
            "skyw_airmass": [1.5],
            "sky_near_label": ["SkyW"],
            "sky_far_label": ["SkyE"],
        }
    )
    input_path = tmp_path / "stack.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            *(fits.ImageHDU(flux, name=name) for name in (
                "FLUX_SCI",
                "FLUX_SKY_NEAR",
                "FLUX_SKY_FAR",
            )),
            *(fits.ImageHDU(np.full_like(flux, value), name=name) for name, value in (
                ("LSF_SCI", 1.1),
                ("LSF_SKY_NEAR", 1.2),
                ("LSF_SKY_FAR", 1.3),
            )),
            fits.BinTableHDU(meta, name="META"),
        ]
    ).writeto(input_path)

    transmissions = []

    class DummyTelluricCalculator:
        def __init__(self, path=None):
            assert path is None

        def match_to_data(
            self,
            target_wave,
            lsf,
            pwv,
            *,
            airmass,
            lsf_in_wavelength,
        ):
            transmissions.append((lsf.copy(), pwv, airmass, lsf_in_wavelength))
            return np.full_like(target_wave, 0.9)

    constructor_calls = []

    class DummyDecomposer:
        def __init__(self, model_wave, **kwargs):
            np.testing.assert_array_equal(model_wave, wave)
            constructor_calls.append(kwargs)

        def fit(self, row_flux, row_ivar, *, verbose):
            assert not verbose
            np.testing.assert_array_equal(row_flux, flux[0] * 2.0)
            np.testing.assert_array_equal(row_ivar, np.ones(wave.size))
            return len(constructor_calls)

    monkeypatch.setattr(lvmdrp.core.fluxcal, "TelluricCalculator", DummyTelluricCalculator)
    monkeypatch.setattr(class_path, DummyDecomposer)
    monkeypatch.setattr(
        decompose_parallel,
        "_install_split_zodi_amplitude_prior",
        lambda decomposer, kind, row_index: None,
    )
    fallback_warnings = []
    monkeypatch.setattr(
        decompose_parallel.warnings,
        "warn",
        lambda message, category, stacklevel: fallback_warnings.append(
            (str(message), category, stacklevel)
        ),
    )
    decompose_parallel.init_worker(
        wave,
        0.5,
        DEFAULT_DATA_ROOT,
        2.0,
        input_path,
        fit_model=fit_model,
        # 3-pixel synthetic grid, and this test is about the telluric kwargs
        # rather than the photon weights (on by default since 2026-09-18).
        fit_pixel_weights=False,
        n_zodi_spline_knots=3,
        zodi_smooth_lambda=0.25,
        diffuse_ratio_bound_dex=0.12,
        diffuse_ratio_nominal=(0.1, 0.6, 0.3),
        diffuse_oh_centre_log10=-0.5,
        diffuse_oh_bound_dex=0.08,
    )
    try:
        for kind in ("sci", "sky1", "sky2"):
            returned_kind, rows = decompose_parallel.fit_chunk_worker((kind, 0, 1))
            assert returned_kind == kind
            # fit_chunk_worker rows are (index, result, reliability_columns)
            assert [(i, r) for i, r, _f in rows] == [(0, len(constructor_calls))]
    finally:
        decompose_parallel._WORKER_HDU.close()

    assert [call["source_airmass"] for call in constructor_calls] == [1.3, 1.5, 1.4]
    assert all(call["pwv_mm"] == expected_pwv for call in constructor_calls)
    assert all(np.array_equal(call["drp_transmission"], np.full(3, 0.9)) for call in constructor_calls)
    assert [float(values[0][0]) for values in transmissions] == [1.1, 1.2, 1.3]
    assert all(values[1:] == (expected_pwv, 1.3, True) for values in transmissions)
    if pwv_med > 0.0:
        assert fallback_warnings == []
    else:
        assert len(fallback_warnings) == 1
        assert "using the LVM DRP default PWV=15.0 mm" in fallback_warnings[0][0]
        assert fallback_warnings[0][1:] == (RuntimeWarning, 2)
    if fit_model in (
        decompose_parallel.PALACE_VNF_PCA30_FIT_MODEL,
        *decompose_parallel.PALACE_VNF_PCA30_SPLIT_ZODI_FIT_MODELS,
    ):
        assert all(call["n_line_amplitude_pca_components"] == 30 for call in constructor_calls)
    if fit_model == decompose_parallel.PALACECORR_VNF_SPLIT_ZODI_FIT_MODEL:
        assert all(
            call["palace_oh_suffix"] == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
            for call in constructor_calls
        )
    if fit_model in decompose_parallel.SPLIT_ZODI_TELLURIC_FIT_MODELS:
        expected = {
            "n_zodi_spline_knots": 3,
            "zodi_smooth_lambda": 0.25,
            "diffuse_ratio_bound_dex": 0.12,
            "diffuse_ratio_nominal": (0.1, 0.6, 0.3),
            "diffuse_oh_centre_log10": -0.5,
            "diffuse_oh_bound_dex": 0.08,
        }
        assert all(
            all(call[key] == value for key, value in expected.items())
            for call in constructor_calls
        )


def test_invalid_airmass_marks_only_that_row_failed(monkeypatch):
    from skysub.sky_decomp import telluric_corrected_lines

    dtype = [
        ("pwv_med", "f8"),
        ("sci_airmass", "f8"),
        ("skye_airmass", "f8"),
        ("skyw_airmass", "f8"),
        ("sky_near_label", "U8"),
        ("sky_far_label", "U8"),
    ]
    meta = np.array([(4.2, -999.9, 1.4, 1.5, "SkyW", "SkyE")], dtype=dtype)
    reasons = []
    constructor_calls = []
    transmission_airmasses = []
    sentinel = object()

    class TemplateDecomposer:
        def __init__(self, model_wave, **kwargs):
            constructor_calls.append((model_wave.copy(), kwargs))

        def failed_input_result(self, reason):
            reasons.append(reason)
            return sentinel

    def calculate_transmission(wave, lsf, pwv, airmass, calculator):
        transmission_airmasses.append(airmass)
        return np.ones_like(wave)

    monkeypatch.setattr(
        telluric_corrected_lines,
        "calculate_drp_transmission",
        calculate_transmission,
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", TemplateDecomposer)
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_FIT_MODEL",
        decompose_parallel.PALACE_VNF_SPLIT_ZODI_FIT_MODEL,
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_META", meta)
    monkeypatch.setattr(decompose_parallel, "_WORKER_WAVE", np.arange(3.0))
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_LSF",
        {"sci": np.ones((1, 3))},
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_TELLURIC_CALCULATOR", object())
    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER_KWARGS", {})
    monkeypatch.setattr(decompose_parallel, "_WORKER_SCIENCE_LINE_MASK", None)

    result, _flags = decompose_parallel._fit_worker_row(
        "sci", 0, np.ones(3), np.ones(3)
    )

    assert result is sentinel
    assert reasons == [
        "invalid_airmass: row=0, role=sci, "
        "sci_airmass=-999.9, source_airmass=-999.9"
    ]
    assert transmission_airmasses == [1.0]
    assert constructor_calls[0][1]["source_airmass"] == 1.0

    meta["sci_airmass"] = 1.3
    meta["sky_near_label"] = "unknown"
    with pytest.raises(ValueError, match="Unknown sky_near_label"):
        decompose_parallel._fit_worker_row(
            "sky1", 0, np.ones(3), np.ones(3)
        )


def test_lsf_row_defect_classifies_what_can_and_cannot_be_repaired():
    clean = np.ones(6)
    good, reason = decompose_parallel._lsf_row_defect(clean)
    assert reason is None and good.all()

    # One isolated interior bad pixel is repairable, so no reason is given.
    repairable = np.array([1.0, 1.0, 0.0, 1.0, 1.0, 1.0])
    good, reason = decompose_parallel._lsf_row_defect(repairable)
    assert reason is None and not good.all()

    for row, expected in (
        (np.zeros(6), "no finite positive pixel"),
        (np.full(6, np.nan), "no finite positive pixel"),
        (np.array([0.0, 1.0, 1.0, 1.0, 1.0, 1.0]), "a bad pixel at the array edge"),
        (np.array([1.0, 1.0, 1.0, 1.0, 1.0, -1.0]), "a bad pixel at the array edge"),
        (np.array([1.0, 1.0, 0.0, 0.0, 1.0, 1.0]), "adjacent bad pixels"),
    ):
        _good, reason = decompose_parallel._lsf_row_defect(row)
        assert reason == expected, row


def test_unusable_lsf_marks_only_that_row_failed(monkeypatch):
    """A row with no usable detector LSF must fail alone, with NaN coefficients.

    `_sanitised_lsf_row` refuses to fabricate an LSF curve for such a row, and
    the telluric decomposition cannot be constructed without one.  That used
    to raise out of `_fit_worker_row` and kill the whole chunk; it must instead
    produce a same-schema failed_input row, exactly as an invalid airmass does.
    """
    from skysub.sky_decomp import telluric_corrected_lines

    dtype = [
        ("pwv_med", "f8"),
        ("sci_airmass", "f8"),
        ("skye_airmass", "f8"),
        ("skyw_airmass", "f8"),
        ("sky_near_label", "U8"),
        ("sky_far_label", "U8"),
    ]
    meta = np.array(
        [(4.2, 1.3, 1.4, 1.5, "SkyW", "SkyE")] * 2,
        dtype=dtype,
    )
    reasons = []
    fitted_rows = []
    transmission_lsf = []
    sentinel = object()
    fitted = object()

    class TemplateDecomposer:
        def __init__(self, model_wave, **kwargs):
            self.kwargs = kwargs

        def failed_input_result(self, reason):
            reasons.append(reason)
            return sentinel

        def fit(self, flux, ivar, **kwargs):
            fitted_rows.append(np.asarray(flux).copy())
            return fitted

    def calculate_transmission(wave, lsf, pwv, airmass, calculator):
        transmission_lsf.append(np.asarray(lsf).copy())
        return np.ones_like(wave)

    monkeypatch.setattr(
        telluric_corrected_lines,
        "calculate_drp_transmission",
        calculate_transmission,
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", TemplateDecomposer)
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_FIT_MODEL",
        decompose_parallel.PALACE_VNF_SPLIT_ZODI_FIT_MODEL,
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_META", meta)
    monkeypatch.setattr(decompose_parallel, "_WORKER_WAVE", np.arange(3.0))
    # Row 0 has no usable LSF pixel at all; row 1 is clean.
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_LSF",
        {"sci": np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])},
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_TELLURIC_CALCULATOR", object())
    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER_KWARGS", {})
    monkeypatch.setattr(decompose_parallel, "_WORKER_SCIENCE_LINE_MASK", None)
    monkeypatch.setattr(decompose_parallel, "_LSF_UNUSABLE_REPORTED", False)
    # The geometry prior is exercised elsewhere; stub it so the clean row below
    # reaches `fit` without needing the full META geometry columns.
    monkeypatch.setattr(
        decompose_parallel,
        "_install_split_zodi_amplitude_prior",
        lambda decomposer, kind, idx: None,
    )

    with pytest.warns(RuntimeWarning, match="unusable_lsf"):
        result, _flags = decompose_parallel._fit_worker_row(
            "sci", 0, np.ones(3), np.ones(3)
        )

    assert result is sentinel
    assert reasons == [
        "unusable_lsf: row=0, role=sci, finite_positive_pixels=0/3, "
        "no finite positive pixel"
    ]
    # The schema-only construction still needs a positive LSF to reach the
    # result shape; it must be a placeholder, never the row's own values.
    assert transmission_lsf and np.all(transmission_lsf[0] > 0.0)
    assert not fitted_rows

    # The next row is untouched by the failure and is fitted normally.
    assert decompose_parallel._fit_worker_row(
        "sci", 1, np.full(3, 2.0), np.ones(3)
    )[0] is fitted
    assert len(reasons) == 1


def test_unusable_lsf_marks_only_that_row_failed_for_moon_zodi(monkeypatch):
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
    reasons = []
    sentinel = object()

    class TemplateDecomposer:
        def failed_input_result(self, reason):
            reasons.append(reason)
            return sentinel

        def fit(self, *args, **kwargs):
            raise AssertionError("an unusable-LSF row must never be fitted")

    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", TemplateDecomposer())
    monkeypatch.setattr(
        decompose_parallel, "_WORKER_FIT_MODEL", decompose_parallel.MOON_ZODI_FIT_MODEL
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_META", meta)
    monkeypatch.setattr(decompose_parallel, "_WORKER_WAVE", np.arange(4.0))
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_LSF",
        {"sci": np.array([[1.0, 0.0, 0.0, 1.0]])},
    )
    monkeypatch.setattr(decompose_parallel, "_WORKER_SCIENCE_LINE_MASK", None)
    monkeypatch.setattr(decompose_parallel, "_LSF_UNUSABLE_REPORTED", False)

    with pytest.warns(RuntimeWarning, match="unusable_lsf"):
        result, _flags = decompose_parallel._fit_worker_row(
            "sci", 0, np.ones(4), np.ones(4)
        )

    assert result is sentinel
    assert reasons == [
        "unusable_lsf: row=0, role=sci, finite_positive_pixels=2/4, "
        "adjacent bad pixels"
    ]


def test_failed_input_result_is_all_nan_for_an_unusable_lsf_row():
    """The recorded row must carry NaN coefficients, not a fabricated fit."""
    from skysub.sky_decomp.lsf_surface_iterative import (
        LSFSurfaceIterativeConfig,
        SkyDecompLSFSurfaceIterative,
    )
    from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT
    from skysub.sky_decomp.result_io import FAILED_INPUT_FIT_STATUS

    wave = np.linspace(3600.0, 9800.0, 400)
    model = SkyDecompLSFSurfaceIterative(
        wave,
        base_dir=DEFAULT_DATA_ROOT,
        lsf_sigma=0.5,
        moon_smooth_lambda=0.1,
        moon_interline_boost=0.0,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=1),
    )
    reason = "unusable_lsf: row=5166, role=sci, finite_positive_pixels=0/400, "\
             "no finite positive pixel"
    result = model.failed_input_result(reason)

    assert result.fit_status == FAILED_INPUT_FIT_STATUS
    assert reason in result.fit_summary
    assert result.coef.shape == (len(result.design_names),)
    for values in (result.coef, result.coef_err, result.bestfit, result.resid):
        assert np.all(np.isnan(values))
    for name, component in result.components.items():
        assert np.all(np.isnan(component)), name
    assert np.isnan(result.r2) and np.isnan(result.reduced_chi2)


def test_batch_preserves_placeholder_and_propagates_unexpected_errors(monkeypatch):
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
    meta = np.zeros(2, dtype=dtype)
    meta[0] = (42, "2024-01-02T03:04:05", 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    meta[1] = (43, "2024-01-02T03:19:05", 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)
    sentinels = {42: object(), 43: object()}

    class PlaceholderResult:
        def fit(self, *_args, **kwargs):
            return sentinels[kwargs["observation"].expnum]

    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", PlaceholderResult())
    monkeypatch.setattr(decompose_parallel, "_WORKER_FACTOR", 1.0)
    monkeypatch.setattr(decompose_parallel, "_WORKER_FIT_MODEL", decompose_parallel.MOON_ZODI_FIT_MODEL)
    monkeypatch.setattr(decompose_parallel, "_WORKER_FLUX", {"sky2": np.ones((2, 2))})
    monkeypatch.setattr(decompose_parallel, "_WORKER_LSF", {"sky2": np.ones((2, 2))})
    monkeypatch.setattr(decompose_parallel, "_WORKER_META", meta)
    monkeypatch.setattr(decompose_parallel, "_WORKER_PROGRESS_QUEUE", None)
    monkeypatch.setattr(decompose_parallel, "_WORKER_SCIENCE_LINE_MASK", None)
    kind, rows = decompose_parallel.fit_chunk_worker(("sky2", 0, 2))
    assert kind == "sky2"
    assert [(i, r) for i, r, _f in rows] == [
        (0, sentinels[42]), (1, sentinels[43])]

    class UnexpectedFailure(PlaceholderResult):
        def fit(self, *_args, **_kwargs):
            raise RuntimeError("unexpected implementation error")

    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", UnexpectedFailure())
    with pytest.raises(RuntimeError, match="unexpected implementation error"):
        decompose_parallel.fit_chunk_worker(("sky2", 0, 1))


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

    validated.clear()
    for fit_model in decompose_parallel.TELLURIC_FIT_MODELS:
        base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(fit_model)
        assert base_dir == default_root.resolve()
        assert data_root == default_root.resolve()
    assert validated == [str(default_root.resolve())] * len(
        decompose_parallel.TELLURIC_FIT_MODELS
    )


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
