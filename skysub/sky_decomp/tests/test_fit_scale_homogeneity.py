"""Brightness-homogeneity of the constrained least-squares core.

The decomposition is run over a corpus spanning more than an order of
magnitude in sky brightness, and the moon/zodi smoothing strengths are
calibrated once for all of it. That calibration is only meaningful if the
solver is positively homogeneous in the data: scaling one row's flux by a
constant (at fixed ivar) must scale every fitted coefficient by the same
constant, so that `moon_smooth_lambda` means the same thing on a bright row
as on a faint one.

`_fit_design` rescales the problem by `data_scale` (target RMS -> 1) and
`col_scale` (columns -> unit norm) before handing it to Clarabel. Because
`col_scale` is computed *after* the `data_scale` division it scales as
1/data_scale, so a D2 penalty built as `D / col_scale` carries an extra
`data_scale**2` -- and `data_scale` is additionally clamped at 1.0, which
breaks the proportionality outright on faint rows. Both effects are
measured in the xfailing tests below.

Scope: only the constraints with a zero right-hand side are homogeneous by
construction (nonnegativity, the adjacent-knot ratio bounds, the moon-share
bracket, the diffuse species bracket). The FeO/OH block cap and the absolute
Leinert zodi bracket are deliberately *not* -- they are physical priors in
native flux units -- so the tests here leave those two priors uninstalled.
"""

import numpy as np
import pytest

from skysub.sky_decomp.fit import SkyDecomp, _build_d2_operator


N_OTHER = 4
N_MOON = 8
N_ZODI = 4
N_PAR = N_OTHER + N_MOON + N_ZODI
MOON_SLICE = slice(N_OTHER, N_OTHER + N_MOON)
ZODI_SLICE = slice(N_OTHER + N_MOON, N_PAR)


def _make_solver(*, moon_lambda=0.0, zodi_lambda=0.0, ratio_bound=0.0):
    """A bare SkyDecomp carrying only what `_fit_design` reads.

    The real constructor needs template FITS files on disk; `_fit_design`
    needs a dozen scalars, so build the object directly and set those.
    """
    self = SkyDecomp.__new__(SkyDecomp)
    self.moon_interline_boost = 0.0
    self.moon_smooth_lambda = float(moon_lambda)
    self.zodi_smooth_lambda = float(zodi_lambda)
    self._d2_moon = _build_d2_operator(N_MOON)
    self._d2_zodi = _build_d2_operator(N_ZODI)
    self.moon_ratio_bound = float(ratio_bound)
    self.zodi_ratio_bound = float(ratio_bound)
    self._moon_relaxed_basis = None
    # Homogeneous priors: off unless a test turns them on.
    self._amp_prior_moon_fraction = None
    self.amp_prior_tol = 0.0
    self.amp_prior_floor = 0.02
    self.diffuse_ratio_bound_dex = 0.0
    self.diffuse_ratio_nominal = None
    # Inhomogeneous priors: always off here (see module docstring).
    self._diffuse_oh_amp = None
    self.diffuse_oh_centre_log10 = None
    self.diffuse_oh_bound_dex = 0.0
    self._amp_prior_zodi_total = None
    self.zodi_amp_bound = 0.0
    return self


def _make_problem(seed=0, n_pix=200, flux_level=1.0):
    """A well-conditioned all-positive design with a strictly interior truth.

    Columns are narrow enough (cond(A) ~ 2.4) that Clarabel reproduces the
    solution to ~1e-8; a broad, collinear basis blurs the tolerances to
    ~1e-4 and would hide the effect being measured.
    """
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, 1.0, n_pix)
    cols = []
    for j in range(N_PAR):
        centre = (j + 0.5) / N_PAR
        cols.append(np.exp(-0.5 * ((x - centre) * N_PAR / 0.5) ** 2))
    design_matrix = np.asarray(cols, dtype=float)          # (n_par, n_pix)
    truth = rng.uniform(0.5, 2.0, size=N_PAR)
    flux = flux_level * (design_matrix.T @ truth)
    flux = flux + flux_level * 0.01 * rng.normal(size=n_pix)
    ivar = np.full(n_pix, 1.0)
    return design_matrix, flux, ivar


def _at_scale(flux, ivar, target_rms):
    """`flux` rescaled so its ivar-weighted RMS -- the raw data_scale -- is
    exactly `target_rms`, which is what the legacy clamp compares against."""
    rms = float(np.sqrt(np.mean((flux * np.sqrt(ivar)) ** 2)))
    return flux * (target_rms / rms)


def _fit(solver, design_matrix, flux, ivar):
    return solver._fit_design(
        design_matrix, flux, ivar,
        moon_slice=MOON_SLICE, zodi_slice=ZODI_SLICE, diffuse_slice=None,
    )["coef"]


# --- control: with no curvature penalty the solve must already be homogeneous

@pytest.mark.parametrize("alpha", [2.0, 0.25])
def test_homogeneous_without_curvature_penalty(alpha):
    """lambda = 0 leaves a pure rescaled NNLS, which is exactly homogeneous."""
    solver = _make_solver(moon_lambda=0.0, zodi_lambda=0.0)
    dm, flux, ivar = _make_problem()
    c1 = _fit(solver, dm, flux, ivar)
    c2 = _fit(solver, dm, alpha * flux, ivar)
    np.testing.assert_allclose(c2, alpha * c1, rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("alpha", [2.0, 0.25])
def test_homogeneous_with_zero_rhs_constraints_active(alpha):
    """The adjacent-knot ratio bounds are homogeneous even when they bind.

    Looser than the unconstrained case: with an active set the two solves
    no longer follow bit-identical interior-point paths, and ~1e-6 is
    Clarabel's own reproducibility floor on this problem.
    """
    solver = _make_solver(moon_lambda=0.0, zodi_lambda=0.0, ratio_bound=0.7)
    dm, flux, ivar = _make_problem(seed=3)
    c1 = _fit(solver, dm, flux, ivar)
    c2 = _fit(solver, dm, alpha * flux, ivar)
    np.testing.assert_allclose(c2, alpha * c1, rtol=1e-5, atol=1e-9)


# --- the property the data_scale**2 normalisation exists to guarantee

@pytest.mark.parametrize("alpha", [2.0, 0.25, 10.0])
def test_homogeneous_with_curvature_penalty(alpha):
    """Scaling the data must scale the coefficients, penalty or no penalty.

    This is the acceptance criterion for the data_scale**2 normalisation:
    with it, `moon_smooth_lambda` means the same thing on every row.
    """
    solver = _make_solver(moon_lambda=1e-1, zodi_lambda=1.0)
    dm, flux, ivar = _make_problem()
    c1 = _fit(solver, dm, flux, ivar)
    c2 = _fit(solver, dm, alpha * flux, ivar)
    np.testing.assert_allclose(c2, alpha * c1, rtol=1e-5, atol=1e-9)


def test_homogeneous_across_the_data_scale_clamp():
    """Two rows either side of the legacy clamp must be proportional.

    The bright row's raw data_scale is well above 1 and the faint row's well
    below, so under the legacy path exactly one of the two is pinned. This is
    the failure mode the clamp adds on top of the data_scale**2 coupling, and
    it is why normalising by the CLAMPED scale would not have been enough --
    hence the epsilon guard rather than a floor at 1.
    """
    solver = _make_solver(moon_lambda=1e-1, zodi_lambda=1.0)
    dm, flux, ivar = _make_problem()
    c_bright = _fit(solver, dm, _at_scale(flux, ivar, 10.0), ivar)
    c_faint = _fit(solver, dm, _at_scale(flux, ivar, 0.1), ivar)
    np.testing.assert_allclose(c_faint, 0.01 * c_bright, rtol=1e-5, atol=1e-12)


def test_faint_rows_are_normalised_not_floored():
    """A row far below the legacy clamp is still scaled to unit RMS.

    Guards the epsilon: a floor anywhere near 1.0 would silently restore the
    clamp for faint rows, which is most of the corpus.
    """
    solver = _make_solver(moon_lambda=1e-1, zodi_lambda=1.0)
    dm, flux, ivar = _make_problem()
    c_a = _fit(solver, dm, _at_scale(flux, ivar, 1e-3), ivar)
    c_b = _fit(solver, dm, _at_scale(flux, ivar, 1e-5), ivar)
    np.testing.assert_allclose(c_b, 0.01 * c_a, rtol=1e-4, atol=1e-16)


def test_degenerate_all_zero_target_does_not_blow_up():
    """The epsilon guard's only job: an all-zero row must not divide by zero."""
    solver = _make_solver(moon_lambda=1e-1, zodi_lambda=1.0)
    dm, flux, ivar = _make_problem()
    coef = _fit(solver, dm, np.zeros_like(flux), ivar)
    assert np.all(np.isfinite(coef))
    # Interior-point, so the answer approaches the c = 0 vertex rather than
    # landing on it; a normal fit of this problem gives coefficients ~1.
    np.testing.assert_allclose(coef, 0.0, atol=1e-4)
