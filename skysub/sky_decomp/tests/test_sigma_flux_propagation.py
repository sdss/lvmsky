"""Unit tests for the LSF-aware sigma propagator (2026-08-16)."""
import numpy as np
import pytest

from skysub.sky_decomp.fit import (
    SkyDecomp,
)


@pytest.fixture(scope="module")
def _model_setup():
    """Small wave grid + a SkyDecomp model whose matrices we can inspect."""
    wave = np.linspace(4000.0, 9500.0, 400)
    model = SkyDecomp(wave, lsf_sigma=1.5, n_spline_knots=25)
    mats = model._matrix_bundle(
        model.matrix_oh,
        model.matrix_moon,
        model.matrix_diffuse,
        model.matrix_atom,
        model.matrix_orc,
        model.matrix_o2,
    )
    return model, mats, wave


def test_zero_coef_err_gives_zero_sigma(_model_setup):
    """coef_err = 0 -> component sigma = 0 everywhere."""
    model, mats, _ = _model_setup
    n_par = sum(m.shape[0] for m in mats.values())
    err = np.zeros(n_par)
    sigmas = model._components_sigma_from_coef_err(err, mats)
    for k in ("oh", "moon", "diffuse", "atom", "orc", "o2",
              "ho2", "feo", "o2ac"):
        np.testing.assert_allclose(sigmas[k], 0.0, atol=1e-15)


def test_nan_coef_err_treated_as_zero(_model_setup):
    """NaN entries in coef_err contribute nothing to the propagated variance."""
    model, mats, _ = _model_setup
    n_par = sum(m.shape[0] for m in mats.values())
    err_all_nan = np.full(n_par, np.nan)
    sigmas = model._components_sigma_from_coef_err(err_all_nan, mats)
    for k in ("oh", "moon", "diffuse", "atom", "orc", "o2"):
        np.testing.assert_allclose(sigmas[k], 0.0, atol=1e-15)


def test_single_coef_matches_column(_model_setup):
    """Setting one coefficient's sigma to 1 gives sigma = |M[j, :]|."""
    model, mats, _ = _model_setup
    sl = model._component_slices(mats)
    n_par = sum(m.shape[0] for m in mats.values())
    err = np.zeros(n_par)
    # Pick a mid oh coefficient and set sigma to 2.5.
    oh_lo = sl["oh"].start
    j = oh_lo + mats["oh"].shape[0] // 2
    err[j] = 2.5
    sigmas = model._components_sigma_from_coef_err(err, mats)
    expected_oh = 2.5 * np.abs(mats["oh"][j - oh_lo, :])
    np.testing.assert_allclose(sigmas["oh"], expected_oh, atol=1e-14)
    # Other component sigmas are zero (float-precision tolerance).
    for k in ("moon", "atom", "orc", "o2", "ho2", "feo", "o2ac", "diffuse"):
        assert float(np.max(np.abs(sigmas[k]))) < 1e-12


def test_length_mismatch_raises(_model_setup):
    model, mats, _ = _model_setup
    with pytest.raises(ValueError, match="coef_err length mismatch"):
        model._components_sigma_from_coef_err(np.array([1.0, 2.0]), mats)


def test_diffuse_sigma_is_quadrature_sum(_model_setup):
    """Component-level ``diffuse`` variance = ho2² + feo² + o2ac²."""
    model, mats, _ = _model_setup
    sl = model._component_slices(mats)
    n_par = sum(m.shape[0] for m in mats.values())
    err = np.zeros(n_par)
    diff_lo = sl["diffuse"].start
    err[diff_lo + 0] = 0.5
    err[diff_lo + 1] = 0.2
    err[diff_lo + 2] = 0.7
    sigmas = model._components_sigma_from_coef_err(err, mats)
    expected = np.sqrt(
        sigmas["ho2"] ** 2 + sigmas["feo"] ** 2 + sigmas["o2ac"] ** 2
    )
    np.testing.assert_allclose(sigmas["diffuse"], expected, rtol=1e-12)
