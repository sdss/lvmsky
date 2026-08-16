import numpy as np
import pytest

from skysub.sky_decomp.fit import SkyDecomp


def test_per_column_chi2_matches_hand_formula():
    """Verify per-column chi2 matches sum(w_j * ivar * resid**2)/sum(w_j)."""
    rng = np.random.default_rng(42)
    n_pix, n_col = 100, 5
    a_full = np.abs(rng.normal(size=(n_pix, n_col)))
    # design_matrix has shape (n_col, n_pix) per _fit_design convention.
    design_matrix = a_full.T
    resid = rng.normal(size=n_pix)
    ivar = np.abs(rng.normal(size=n_pix)) + 0.1
    good = np.ones(n_pix, dtype=bool)
    fallback = 3.14

    out = SkyDecomp._per_column_chi2_from_residuals(
        design_matrix, resid, ivar, good, fallback=fallback,
    )
    expected = np.zeros(n_col)
    for j in range(n_col):
        col_max = np.abs(a_full[good, j]).max()
        w_j = np.abs(a_full[good, j]) / col_max
        num = np.sum(w_j * ivar[good] * resid[good] ** 2)
        den = np.sum(w_j)
        expected[j] = num / den

    np.testing.assert_allclose(out, expected, rtol=1e-12)


def test_per_column_chi2_falls_back_for_zero_support():
    """Columns with all-zero design fall back to the aggregate reduced_chi2."""
    a_full = np.array(
        [
            [1.0, 0.0, 2.0],
            [0.5, 0.0, 1.0],
            [0.2, 0.0, 3.0],
        ],
        dtype=float,
    )
    design_matrix = a_full.T
    resid = np.array([0.1, 0.2, -0.1], dtype=float)
    ivar = np.array([1.0, 2.0, 4.0], dtype=float)
    good = np.array([True, True, True])
    fallback = 2.71

    out = SkyDecomp._per_column_chi2_from_residuals(
        design_matrix, resid, ivar, good, fallback=fallback,
    )
    assert out[1] == pytest.approx(fallback)
    assert out[0] > 0.0
    assert out[2] > 0.0


def test_per_column_chi2_respects_good_mask():
    """Pixels outside `good` are excluded from the pseudo-chi2 sum."""
    rng = np.random.default_rng(7)
    n_pix, n_col = 20, 3
    a_full = np.abs(rng.normal(size=(n_pix, n_col)))
    design_matrix = a_full.T
    resid = rng.normal(size=n_pix)
    ivar = np.abs(rng.normal(size=n_pix)) + 0.1
    good = np.zeros(n_pix, dtype=bool)
    good[:5] = True

    out = SkyDecomp._per_column_chi2_from_residuals(
        design_matrix, resid, ivar, good, fallback=1.0,
    )
    for j in range(n_col):
        col_max = np.abs(a_full[good, j]).max()
        w_j = np.abs(a_full[good, j]) / col_max
        num = np.sum(w_j * ivar[good] * resid[good] ** 2)
        den = np.sum(w_j)
        assert out[j] == pytest.approx(num / den)


def test_per_column_chi2_all_bad_returns_fallback():
    """An all-False good mask returns the fallback for every column."""
    design_matrix = np.random.default_rng(0).normal(size=(4, 10))
    resid = np.zeros(10)
    ivar = np.ones(10)
    good = np.zeros(10, dtype=bool)
    out = SkyDecomp._per_column_chi2_from_residuals(
        design_matrix, resid, ivar, good, fallback=5.0,
    )
    np.testing.assert_array_equal(out, np.full(4, 5.0))


def test_coef_err_active_set_uses_per_column_chi2():
    """Per-column inflation multiplies the diagonal variance element-wise."""
    n_par = 3
    coef = np.array([1.0, 2.0, 0.0], dtype=float)  # last coefficient inactive
    # One solve group: identity-ish Hessian in scaled space.
    p_solve = np.eye(3, dtype=float) * 4.0  # sigma_scaled ~ sqrt(1/4)
    col_scale = np.ones(3, dtype=float)
    data_scale = 1.0
    source_p = [p_solve] * n_par
    source_col_scale = [col_scale] * n_par
    source_data_scale = [data_scale] * n_par
    source_local_index = [0, 1, 2]

    reduced_chi2 = 1.0  # no aggregate inflation
    per_col = np.array([1.0, 9.0, 1.0], dtype=float)  # inflate col 1 by 9

    err_scalar = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2,
    )
    err_percol = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2, per_column_chi2=per_col,
    )
    # Column 0 unchanged, column 1 sigma inflated by sqrt(9) = 3x,
    # column 2 stays NaN (inactive).
    assert err_percol[0] == pytest.approx(err_scalar[0])
    assert err_percol[1] == pytest.approx(err_scalar[1] * 3.0)
    assert np.isnan(err_percol[2])
    assert np.isnan(err_scalar[2])


def test_coef_err_active_set_per_column_never_shrinks():
    """per_column_chi2 < 1 is clipped up to 1 (never below the Cramer-Rao floor)."""
    coef = np.array([1.0, 1.0], dtype=float)
    p_solve = np.eye(2, dtype=float)
    source_p = [p_solve, p_solve]
    source_col_scale = [np.ones(2), np.ones(2)]
    source_data_scale = [1.0, 1.0]
    source_local_index = [0, 1]

    per_col = np.array([0.01, 100.0], dtype=float)
    err = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2=1.0, per_column_chi2=per_col,
    )
    baseline = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2=1.0,
    )
    # Column 0's inflation was clipped to 1 -> matches baseline.
    assert err[0] == pytest.approx(baseline[0])
    # Column 1's inflation is 100 -> sigma multiplied by 10.
    assert err[1] == pytest.approx(baseline[1] * 10.0)


def test_coef_err_active_set_scalar_fallback_preserved():
    """per_column_chi2=None reproduces the pre-2026-08-16 scalar inflation."""
    coef = np.array([1.0, 1.0], dtype=float)
    p_solve = np.eye(2, dtype=float)
    source_p = [p_solve, p_solve]
    source_col_scale = [np.ones(2), np.ones(2)]
    source_data_scale = [1.0, 1.0]
    source_local_index = [0, 1]

    err_scalar = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2=4.0,
    )
    err_none = SkyDecomp._coef_err_active_set(
        coef, source_p, source_col_scale, source_data_scale,
        source_local_index, reduced_chi2=4.0, per_column_chi2=None,
    )
    np.testing.assert_allclose(err_scalar, err_none)
    # sqrt(reduced_chi2)=2, applied to sigma from unit Hessian -> sigma=2.
    np.testing.assert_allclose(err_scalar, [2.0, 2.0])
