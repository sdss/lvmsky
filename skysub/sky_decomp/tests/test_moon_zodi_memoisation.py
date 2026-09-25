"""Opt-in memoisation of MoonZodiPhysicalModel.predict.

The cache builder turns memoisation on so the fitted and physics-only
evaluations of an arm share one LSF projection operator and one ephemeris.  It
must give the same answer as the plain path, and it must be off by default so a
process-wide cache cannot outlive a change to module state.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from skysub.sky_decomp import moon_zodi_model
from skysub.sky_decomp.moon_zodi_model import (
    MoonZodiObservation,
    MoonZodiPhysicalModel,
    geometry_amplitude_prior,
)

REFERENCE = Path(__file__).parent / "data" / "moon_zodi_predictor_reference_v1.npz"


def _case(index):
    with np.load(REFERENCE, allow_pickle=False) as reference:
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
    return wave, lsf, observation


@pytest.fixture(autouse=True)
def _memoisation_off_afterwards():
    yield
    moon_zodi_model.set_memoisation(False)


def test_memoisation_is_off_by_default():
    assert moon_zodi_model._MEMOISE is False


def test_memoised_predict_matches_plain_predict():
    wave, lsf, observation = _case(1)
    model = MoonZodiPhysicalModel()
    plain = model.predict(wave, lsf, observation, physical_to_fit_flux_scale=1e14)
    plain_prior = geometry_amplitude_prior(
        wave, lsf, observation, physical_to_fit_flux_scale=1e14)

    moon_zodi_model.set_memoisation(True)
    for _ in range(2):   # second pass is served from the caches
        memo = model.predict(wave, lsf, observation, physical_to_fit_flux_scale=1e14)
        memo_prior = geometry_amplitude_prior(
            wave, lsf, observation, physical_to_fit_flux_scale=1e14)
        np.testing.assert_allclose(memo.moon, plain.moon, rtol=1e-13, atol=0.0)
        np.testing.assert_allclose(memo.zodi, plain.zodi, rtol=1e-13, atol=0.0)
        np.testing.assert_allclose(memo_prior, plain_prior, rtol=1e-13, atol=0.0)
        assert memo.state.geometry == plain.state.geometry
    assert len(moon_zodi_model._PROJECTION_CACHE) == 1


def test_set_memoisation_empties_the_caches(monkeypatch):
    wave, lsf, observation = _case(1)
    moon_zodi_model.set_memoisation(True)
    MoonZodiPhysicalModel().predict(wave, lsf, observation, physical_to_fit_flux_scale=1e14)
    assert moon_zodi_model._PROJECTION_CACHE
    moon_zodi_model.set_memoisation(False)
    assert not moon_zodi_model._PROJECTION_CACHE
    assert moon_zodi_model._cached_midpoint_geometry.cache_info().currsize == 0

    # With memoisation off, module state changes take effect immediately.
    def reject(*_args):
        raise moon_zodi_model._LeinertDomainError("zodi_invalid_near_sun_cell", "x")

    monkeypatch.setattr(moon_zodi_model, "_interpolate_leinert", reject)
    with pytest.raises(moon_zodi_model.MoonZodiInvalidObservationError):
        MoonZodiPhysicalModel().predict(
            wave, lsf, observation, physical_to_fit_flux_scale=1e14)
