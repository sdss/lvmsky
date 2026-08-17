"""Contract tests for the frozen production Moon/Zodi predictor."""

from __future__ import annotations

import ast
from dataclasses import replace
from pathlib import Path
import shutil

import numpy as np
import pytest

from skysub.sky_decomp import moon_zodi_model
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_DIR,
    DEFAULT_DATA_ROOT,
    MODEL_ID,
    NATIVE_GRID_SHA256,
    MoonZodiInvalidObservationError,
    MoonZodiObservation,
    MoonZodiPhysicalModel,
    file_sha256,
    validate_decomposition_data_root,
    wave_sha256,
)


REFERENCE = Path(__file__).parent / "data" / "moon_zodi_predictor_reference_v1.npz"
EXPECTED_ASSETS = {
    "jpl_de432s_short_planetary_ephemeris.bsp": "363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b",
    "meftah_solar_hrs_disk_integrated_v1_1_vacuum_velocity_step_2kms.npz": "cab38bff4cffca72855cd7f305712c2743df4af8ba66d60a1575cdeba5ad8fba",
    "eso_skycalc_rolo_moon_albedo.dat": "86b9f9860fabb283de6659aabee8959186dc03e9d081aaa7c2761c2869ff16cc",
    "eso_skycalc_leinert_zodiacal_light.dat": "93d97a0635f2b49a96b2a7f840303f36d9254d97a8aafdeb24ddbec01a4bdbbe",
}
GEOMETRY_NAMES = (
    "target_altitude_deg",
    "target_airmass",
    "sun_altitude_deg",
    "moon_altitude_deg",
    "moon_airmass",
    "moon_separation_deg",
    "signed_phase_deg",
    "moon_distance_km",
    "sun_moon_distance_km",
    "solar_elongation_deg",
    "ecliptic_lon_relative_deg",
    "ecliptic_latitude_deg",
    "moon_velocity_kms",
    "zodi_velocity_kms",
    "zodi_b500",
)


def _normalized_max(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(
        np.max(np.abs(actual - expected))
        / max(float(np.max(np.abs(expected))), 1.0e-30)
    )


def test_bundle_checksums_and_manifest_contract():
    bundle = validate_decomposition_data_root(str(DEFAULT_DATA_ROOT))
    assert bundle["bundle_id"] == "sky_decomp_moon_zodi_lsf_surface_iterative_v1"
    model = MoonZodiPhysicalModel()
    assert model.manifest["model_id"] == MODEL_ID
    assert model.manifest["scientific_status"] == (
        "diagnostic_only_decomposition_nonregression_failed"
    )
    assert len(model.parameter_values) == 12
    assert dict((name, digest) for name, digest, _ in model.asset_records) == EXPECTED_ASSETS
    for name, expected in EXPECTED_ASSETS.items():
        assert file_sha256(DEFAULT_DATA_DIR / name) == expected


@pytest.mark.parametrize("failure", ["missing", "changed"])
def test_bundle_loader_reports_exact_asset_failure(tmp_path, failure):
    copied = tmp_path / failure / "moon_zodi"
    shutil.copytree(DEFAULT_DATA_DIR, copied)
    target = copied / "eso_skycalc_leinert_zodiacal_light.dat"
    if failure == "missing":
        target.unlink()
        match = "asset is missing.*eso_skycalc_leinert_zodiacal_light.dat"
    else:
        target.write_bytes(target.read_bytes() + b"\n# changed\n")
        match = "Checksum mismatch.*eso_skycalc_leinert_zodiacal_light.dat"
    with pytest.raises((FileNotFoundError, ValueError), match=match):
        MoonZodiPhysicalModel(data_dir=copied)


@pytest.mark.parametrize("failure", ["missing", "changed"])
def test_complete_data_root_reports_exact_palace_failure(tmp_path, failure):
    copied = tmp_path / failure / "data"
    shutil.copytree(DEFAULT_DATA_ROOT, copied)
    target = copied / "palace" / "PMD" / "pmd_popmodel_O2.dat"
    if failure == "missing":
        target.unlink()
        match = "PALACE asset is missing.*pmd_popmodel_O2.dat"
    else:
        target.write_bytes(target.read_bytes() + b"\n# changed\n")
        match = "Checksum mismatch.*pmd_popmodel_O2.dat"
    with pytest.raises((FileNotFoundError, ValueError), match=match):
        validate_decomposition_data_root.cache_clear()
        validate_decomposition_data_root(str(copied))


def test_predictor_matches_three_frozen_jax_reference_cases():
    with np.load(REFERENCE, allow_pickle=False) as reference:
        assert tuple(str(value) for value in reference["geometry_names"]) == GEOMETRY_NAMES
        wave = np.asarray(reference["wave"], dtype=np.float64)
        model = MoonZodiPhysicalModel()
        assert wave.shape == (12_401,)
        assert wave.dtype == np.float64
        assert wave_sha256(wave) == NATIVE_GRID_SHA256
        for index in range(3):
            observation = MoonZodiObservation(
                expnum=int(reference["expnum"][index]),
                date_obs=str(reference["date_obs"][index]),
                role="sky_far",
                target_ra_deg=float(reference["ra_deg"][index]),
                target_dec_deg=float(reference["dec_deg"][index]),
                exposure_seconds=900.0,
                exposure_seconds_source="assumed_900s",
            )
            prediction = model.predict(
                wave,
                np.asarray(reference["lsf"][index], dtype=np.float64),
                observation,
                physical_to_fit_flux_scale=1.0e14,
            )
            assert prediction.moon.dtype == np.float64
            assert prediction.zodi.dtype == np.float64
            assert _normalized_max(prediction.moon, reference["moon"][index]) <= 1.0e-9
            assert _normalized_max(prediction.zodi, reference["zodi"][index]) <= 1.0e-9
            actual_geometry = np.asarray(
                [getattr(prediction.state.geometry, name) for name in GEOMETRY_NAMES],
                dtype=np.float64,
            )
            np.testing.assert_allclose(
                actual_geometry,
                reference["geometry"][index],
                rtol=1.0e-12,
                atol=1.0e-10,
            )
        assert "moon_below_horizon" in model.predict(
            wave,
            np.asarray(reference["lsf"][0], dtype=np.float64),
            MoonZodiObservation(
                int(reference["expnum"][0]),
                str(reference["date_obs"][0]),
                "sky_far",
                float(reference["ra_deg"][0]),
                float(reference["dec_deg"][0]),
                900.0,
                "assumed_900s",
            ),
            physical_to_fit_flux_scale=1.0e14,
        ).state.flags


def test_target_below_horizon_raises_typed_rejection_with_partial_geometry():
    with np.load(REFERENCE, allow_pickle=False) as reference:
        index = 1
        wave = np.asarray(reference["wave"], dtype=np.float64)
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
    invalid = replace(
        observation,
        target_ra_deg=(observation.target_ra_deg + 180.0) % 360.0,
        target_dec_deg=-observation.target_dec_deg,
    )
    with pytest.raises(MoonZodiInvalidObservationError) as caught:
        MoonZodiPhysicalModel().predict(
            wave,
            lsf,
            invalid,
            physical_to_fit_flux_scale=1.0e14,
        )
    error = caught.value
    assert error.reason == "target_below_horizon"
    assert error.observation == invalid
    assert error.geometry.target_altitude_deg <= 0.0
    assert np.isnan(error.geometry.target_airmass)
    assert error.geometry.midpoint_utc


def test_native_grid_and_precision_are_fail_fast():
    with np.load(REFERENCE, allow_pickle=False) as reference:
        wave = np.asarray(reference["wave"], dtype=np.float64)
    model = MoonZodiPhysicalModel()
    model.validate_wave(wave)
    with pytest.raises(ValueError, match="float64"):
        model.validate_wave(wave.astype(np.float32))
    with pytest.raises(ValueError, match="12401"):
        model.validate_wave(wave[:-1])
    changed = wave.copy()
    changed[100] = np.nextafter(changed[100], np.inf)
    with pytest.raises(ValueError, match="fingerprint changed"):
        model.validate_wave(changed)
    with pytest.raises(ValueError, match="fingerprint changed"):
        model.validate_wave(wave[::-1].copy())


def test_production_module_has_no_jax_or_experiment_imports():
    tree = ast.parse(Path(moon_zodi_model.__file__).read_text(encoding="utf-8"))
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert not any(name == "jax" or name.startswith("jax.") for name in imported)
    assert not any("experiment" in name for name in imported)
