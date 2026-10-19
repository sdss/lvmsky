"""The sparse Gram product in ``_fit_design`` must pose the same QP as the dense one."""

import clarabel
import numpy as np
import pytest

import skysub.sky_decomp.fit as fit_module
from skysub.sky_decomp.fit import SkyDecomp

WAVE = np.linspace(4000.0, 9500.0, 400)
SOLVER = clarabel.DefaultSolver


@pytest.fixture(scope="module")
def spectrum():
    model = SkyDecomp(WAVE, lsf_sigma=1.5, n_spline_knots=25)
    truth = 1e-16 * (model.matrix_moon.sum(0) + 0.3 * model.matrix_oh.sum(0) + 0.2 * model.matrix_diffuse.sum(0))
    sigma = 0.01 * np.median(truth)
    flux = truth + np.random.default_rng(0).normal(0.0, sigma, WAVE.size)
    return flux, np.full_like(WAVE, sigma**-2)


def _fit(monkeypatch, density, flux, ivar):
    """Fit with every design forced through one Gram path; record each QP's P and q."""
    problems = []

    def recording(p, q, *args):
        problems.append((p.toarray(), np.array(q)))
        return SOLVER(p, q, *args)

    monkeypatch.setattr(fit_module, "_SPARSE_GRAM_DENSITY", density)
    monkeypatch.setattr(fit_module.clarabel, "DefaultSolver", recording)
    result = SkyDecomp(WAVE, lsf_sigma=1.5, n_spline_knots=25).fit(flux, ivar, verbose=False)
    return result, problems


def test_sparse_and_dense_gram_give_the_same_fit(monkeypatch, spectrum):
    dense, dense_problems = _fit(monkeypatch, 0.0, *spectrum)
    sparse, sparse_problems = _fit(monkeypatch, 1.1, *spectrum)

    # Identical inputs: the first QP agrees to rounding.  Later QPs start from a
    # solution that is itself only converged to the solver tolerance.
    (p_dense, q_dense), (p_sparse, q_sparse) = dense_problems[0], sparse_problems[0]
    np.testing.assert_allclose(p_sparse, p_dense, rtol=0.0, atol=1e-12 * np.abs(p_dense).max())
    np.testing.assert_allclose(q_sparse, q_dense, rtol=0.0, atol=1e-12 * np.abs(q_dense).max())
    assert len(sparse_problems) == len(dense_problems)

    # Clarabel stops at a 1e-8 relative gap, so the fitted values agree to that level.
    assert sparse.fit_status == dense.fit_status == "Solved"
    np.testing.assert_allclose(sparse.coef, dense.coef, rtol=0.0, atol=1e-7 * np.abs(dense.coef).max())
    finite = np.isfinite(dense.coef_err)
    np.testing.assert_array_equal(np.isfinite(sparse.coef_err), finite)
    np.testing.assert_allclose(
        sparse.coef_err[finite], dense.coef_err[finite], rtol=0.0, atol=1e-7 * np.abs(dense.coef_err[finite]).max()
    )
