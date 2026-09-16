import hashlib
import json

import numpy as np
import pytest
from astropy.table import Table

from skysub.sky_decomp.residual_pca import (
    LINE_AMPLITUDE_PCA_ASSET,
    RESIDUAL_PCA_ASSET,
    VN_LINE_AMPLITUDE_PCA_ASSET,
    VNF_COEFFICIENT_LINE_AMPLITUDE_PCA_ASSET,
    VNF_COEFFICIENT_PCA_PREP_ASSET,
    VNF_LINE_ADJOINT_PCA_ASSET,
    SPLIT_ZODI_VNF_LINE_AMPLITUDE_PCA_ASSET,
    SkyDecompPalaceAijcVNFLineAmplitudePCA,
    SkyDecompPalaceAijcVNFNivContinuumLineAmplitudePCA,
    SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30,
    SkyDecompTelluricCorrectedLinesLineAmplitudePCA,
    SkyDecompTelluricCorrectedLinesPalaceAijcVN,
    SkyDecompTelluricCorrectedLinesResidualPCA,
    SkyDecompTelluricCorrectedLinesVNLineAmplitudePCA,
    _load_basis,
    _load_line_amplitude_basis,
    _load_oh_coefficient_basis,
)
from skysub.sky_decomp.fit import CAP_WAVE, decode_hitran_id, vac_to_air
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT, NATIVE_GRID_SHA256


def test_split_zodi_pca30_uses_the_neutral_production_asset_and_legacy_alias():
    assert SPLIT_ZODI_VNF_LINE_AMPLITUDE_PCA_ASSET.endswith(
        "palace_aijc_vnf_split_zodi_line_amplitude_pca30_v1.npz"
    )
    assert (DEFAULT_DATA_ROOT / SPLIT_ZODI_VNF_LINE_AMPLITUDE_PCA_ASSET).is_file()
    assert (
        SkyDecompPalaceAijcVNFNivContinuumLineAmplitudePCA
        is SkyDecompPalaceAijcVNFSplitZodiLineAmplitudePCA30
    )


def test_split_zodi_pca30_manifest_preserves_the_frozen_basis_identity():
    path = DEFAULT_DATA_ROOT / SPLIT_ZODI_VNF_LINE_AMPLITUDE_PCA_ASSET
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"].item()))
    manifest = json.loads((DEFAULT_DATA_ROOT / "bundle_manifest.json").read_text())
    contract = manifest["split_zodi_vnf_line_amplitude_pca_contract"]

    assert contract["basis_id"] == metadata["basis_id"]
    assert contract["production_basis_id"] == (
        "palace-aijc-vnf-split-zodi-line-amplitude-pca30-v1"
    )


def test_frozen_basis_contract():
    path = DEFAULT_DATA_ROOT / RESIDUAL_PCA_ASSET
    wave, mean, components, metadata = _load_basis(str(path))

    assert wave.shape == mean.shape == (12_401,)
    assert components.shape == (30, 12_401)
    assert metadata["wave_sha256"] == NATIVE_GRID_SHA256
    assert metadata["training_spectra"] == 1000
    assert metadata["available_component_counts"] == [10, 20, 30]
    assert metadata["excluded_expnums"] == [45851]


def test_load_basis_preserves_the_exact_float64_grid(tmp_path):
    wave = np.arange(35.0, dtype=np.float64)
    components = np.linalg.qr(
        np.random.default_rng(7).normal(size=(35, 30))
    )[0].T
    source_oh = tmp_path / "palace/PMD/source.dat"
    source_oh.parent.mkdir(parents=True)
    source_oh.write_text("test\n", encoding="utf-8")
    metadata = {
        "wave_sha256": hashlib.sha256(
            np.ascontiguousarray(wave.astype("<f8", copy=False)).tobytes()
        ).hexdigest(),
        "source_oh_asset": "palace/PMD/source.dat",
        "source_oh_asset_sha256": hashlib.sha256(source_oh.read_bytes()).hexdigest(),
    }
    path = tmp_path / "residual_pca/basis.npz"
    path.parent.mkdir()
    np.savez_compressed(
        path,
        wave=wave,
        residual_mean=np.zeros_like(wave),
        components=components,
        metadata_json=np.asarray(json.dumps(metadata)),
    )

    loaded_wave, mean, loaded_components, loaded_metadata = _load_basis(str(path))

    np.testing.assert_array_equal(loaded_wave, wave)
    np.testing.assert_array_equal(mean, np.zeros_like(wave))
    np.testing.assert_allclose(loaded_components, components)
    assert loaded_metadata == metadata


def test_frozen_line_amplitude_basis_contract():
    path = DEFAULT_DATA_ROOT / LINE_AMPLITUDE_PCA_ASSET
    wave, names, line_wave, line_group, mean, components, metadata = (
        _load_line_amplitude_basis(str(path))
    )

    assert wave.shape == (12_401,)
    assert names.shape == line_wave.shape == line_group.shape == mean.shape == (11_552,)
    assert components.shape == (30, 11_552)
    assert metadata["wave_sha256"] == NATIVE_GRID_SHA256
    assert metadata["training_spectra"] == 1000
    assert metadata["available_component_counts"] == [10, 20, 30]
    assert metadata["amplitude_constraints"] == (
        "none; all 11552 individual amplitudes signed"
    )
    assert metadata["physical_grouping_used"] is False
    assert metadata["observable_line_transitions"] == 11_538
    assert metadata["excluded_expnums"] == [45851]


def test_frozen_vn_line_amplitude_basis_contract():
    path = DEFAULT_DATA_ROOT / VN_LINE_AMPLITUDE_PCA_ASSET
    wave, names, line_wave, line_group, mean, components, metadata = (
        _load_line_amplitude_basis(str(path))
    )

    assert wave.shape == (12_401,)
    assert names.shape == line_wave.shape == line_group.shape == mean.shape == (11_552,)
    assert components.shape == (30, 11_552)
    assert metadata["wave_sha256"] == NATIVE_GRID_SHA256
    assert metadata["source_oh_group_keys"] == ["v_upper", "N_upper"]
    assert metadata["source_oh_groups"] == 188
    assert metadata["source_total_emission_groups"] == 196
    assert metadata["input_spectra"] == 1000
    assert metadata["training_spectra"] == 999
    assert metadata["preexcluded_expnums"] == [45851]
    assert metadata["rms_outlier_expnums"] == [45858]
    assert metadata["observable_line_transitions"] == 11_538
    assert metadata["physical_grouping_used"] is False


def test_load_line_amplitude_basis_validates_names_and_sources(tmp_path):
    wave = np.arange(35.0, dtype=np.float64)
    names = np.asarray([f"line_{index}" for index in range(35)])
    components = np.linalg.qr(
        np.random.default_rng(9).normal(size=(35, 30))
    )[0].T
    source = tmp_path / "palace/PMD/source.dat"
    source.parent.mkdir(parents=True)
    source.write_text("test\n", encoding="utf-8")
    metadata = {
        "wave_sha256": hashlib.sha256(
            np.ascontiguousarray(wave.astype("<f8", copy=False)).tobytes()
        ).hexdigest(),
        "line_names_sha256": hashlib.sha256(
            "\n".join(names.tolist()).encode()
        ).hexdigest(),
        "source_assets_sha256": {
            "palace/PMD/source.dat": hashlib.sha256(source.read_bytes()).hexdigest()
        },
    }
    path = tmp_path / "residual_pca/basis.npz"
    path.parent.mkdir()
    np.savez_compressed(
        path,
        wave=wave,
        line_names=names,
        line_wave=np.arange(names.size, dtype=np.float64),
        line_group=np.arange(names.size, dtype=np.int64),
        amplitude_mean=np.zeros(names.size),
        components=components,
        metadata_json=np.asarray(json.dumps(metadata)),
    )

    loaded = _load_line_amplitude_basis(str(path))

    np.testing.assert_array_equal(loaded[0], wave)
    np.testing.assert_array_equal(loaded[1], names)
    np.testing.assert_array_equal(loaded[2], np.arange(names.size, dtype=np.float64))
    np.testing.assert_array_equal(loaded[3], np.arange(names.size, dtype=np.int64))
    np.testing.assert_array_equal(loaded[4], np.zeros(names.size))
    np.testing.assert_allclose(loaded[5], components)
    assert loaded[6] == metadata


def test_load_oh_coefficient_basis_preserves_orthonormal_components(tmp_path):
    wave = np.arange(40.0, dtype=np.float64)
    names = np.asarray([f"OH_{index:03d}" for index in range(357)])
    components = np.linalg.qr(
        np.random.default_rng(11).normal(size=(357, 3))
    )[0].T
    source = tmp_path / "palace/PMD/source.dat"
    source.parent.mkdir(parents=True)
    source.write_text("test\n", encoding="utf-8")
    metadata = {
        "wave_sha256": hashlib.sha256(
            np.ascontiguousarray(wave.astype("<f8", copy=False)).tobytes()
        ).hexdigest(),
        "source_assets_sha256": {
            "palace/PMD/source.dat": hashlib.sha256(source.read_bytes()).hexdigest()
        },
    }
    path = tmp_path / "residual_pca/basis.npz"
    path.parent.mkdir()
    np.savez_compressed(
        path,
        wave=wave,
        coefficient_names=names,
        coefficient_mean=np.arange(357, dtype=np.float64),
        components=components,
        metadata_json=np.asarray(json.dumps(metadata)),
    )

    loaded = _load_oh_coefficient_basis(str(path))

    np.testing.assert_array_equal(loaded[0], wave)
    np.testing.assert_array_equal(loaded[1], names)
    np.testing.assert_array_equal(loaded[2], np.arange(357, dtype=np.float64))
    np.testing.assert_allclose(loaded[3], components)
    assert loaded[4] == metadata


def test_frozen_vnf_coefficient_pca_prep_contract():
    wave, names, mean, components, metadata = _load_oh_coefficient_basis(
        str(DEFAULT_DATA_ROOT / VNF_COEFFICIENT_PCA_PREP_ASSET)
    )

    assert wave.shape == (12_401,)
    assert names.shape == mean.shape == (357,)
    assert components.shape == (18, 357)
    assert metadata["source_oh_group_keys"] == ["v_upper", "N_upper", "F_upper"]
    assert metadata["input_spectra"] == 1000
    assert metadata["training_spectra"] == 994
    assert metadata["selected_components"] == 18
    assert metadata["rms_outlier_expnums"] == [4114, 6374, 45858, 13197, 23762, 39028]


def test_frozen_vnf_coefficient_line_amplitude_pca_contract():
    wave, names, line_wave, line_group, mean, components, metadata = (
        _load_line_amplitude_basis(
            str(DEFAULT_DATA_ROOT / VNF_COEFFICIENT_LINE_AMPLITUDE_PCA_ASSET)
        )
    )

    assert wave.shape == (12_401,)
    assert names.shape == line_wave.shape == line_group.shape == mean.shape == (11_552,)
    assert components.shape == (200, 11_552)
    assert metadata["training_spectra"] == 994
    assert metadata["selected_components"] == 50
    assert metadata["variance_target_components"] == 993
    assert metadata["production_default_reaches_variance_target"] is False
    assert metadata["observable_line_transitions"] == 11_538


def test_frozen_vnf_line_adjoint_pca_contract():
    wave, names, line_wave, line_group, mean, components, metadata = (
        _load_line_amplitude_basis(
            str(DEFAULT_DATA_ROOT / VNF_LINE_ADJOINT_PCA_ASSET)
        )
    )

    assert wave.shape == (12_401,)
    assert names.shape == line_wave.shape == line_group.shape == mean.shape == (11_552,)
    assert components.shape == (500, 11_552)
    assert metadata["training_spectra"] == 994
    assert metadata["input_spectra"] == 1000
    assert metadata["selected_components"] == 500
    assert metadata["source_model"] == "SkyDecompTelluricCorrectedLinesPalaceAijc"
    assert metadata["full_dictionary_amplitude_fit"] is False
    assert metadata["feature_transform"] == "diag(A.T @ A)^-1 @ A.T @ residual"
    assert metadata["observable_line_transitions"] == 11_538


def test_palace_aijc_is_the_only_oh_strength_used_by_the_pca_method():
    group = Table(
        {"Aij": [100.0, 200.0], "Aijc": [2.0, 3.0], "gi": [5.0, 7.0]}
    )
    np.testing.assert_array_equal(
        SkyDecompTelluricCorrectedLinesResidualPCA._oh_amplitude(group),
        np.array([10.0, 21.0]),
    )
    assert SkyDecompPalaceAijcVNFLineAmplitudePCA.oh_group_keys == (
        "v_upper",
        "N_upper",
        "F_upper",
    )


def test_vn_grouping_has_188_oh_groups_on_the_frozen_native_grid():
    with np.load(DEFAULT_DATA_ROOT / RESIDUAL_PCA_ASSET, allow_pickle=False) as data:
        wave = np.asarray(data["wave"], dtype=np.float64)
    table = Table.read(
        DEFAULT_DATA_ROOT
        / "palace/PMD/pmd_popmodel_OH_telluric_upper_parity_lsf_adam_25000_v1.dat",
        format="ascii.basic",
        guess=False,
        comment="#",
        fast_reader=False,
    )
    table["wave"] = vac_to_air(np.asarray(table["lam"], dtype=float) * 1.0e4)
    use = (table["wave"] >= wave[0] - CAP_WAVE) & (
        table["wave"] <= wave[-1] + CAP_WAVE
    )
    table = decode_hitran_id(table[use])

    assert SkyDecompTelluricCorrectedLinesPalaceAijcVN.oh_group_keys == (
        "v_upper",
        "N_upper",
    )
    assert len(table.group_by(("v_upper", "N_upper")).groups) == 188
    assert len(table.group_by(("v_upper", "N_upper", "F_upper")).groups) == 357


def test_only_frozen_component_counts_are_accepted():
    with pytest.raises(ValueError, match="must be 10, 20, or 30"):
        SkyDecompTelluricCorrectedLinesResidualPCA(
            np.arange(4.0), n_residual_pca_components=3
        )
    with pytest.raises(ValueError, match="must be 10, 20, or 30"):
        SkyDecompTelluricCorrectedLinesLineAmplitudePCA(
            np.arange(4.0), n_line_amplitude_pca_components=3
        )
    with pytest.raises(ValueError, match="must be 10, 20, or 30"):
        SkyDecompTelluricCorrectedLinesVNLineAmplitudePCA(
            np.arange(4.0), n_line_amplitude_pca_components=3
        )
