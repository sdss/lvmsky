from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.interpolate import BSpline

from skysub.sky_decomp import solar_lsf
from skysub.sky_decomp.fit import LSF_CHANNELS, read_static_matrix, vac_to_air
from skysub.sky_decomp.lsf_spline2d import (
    _integrated_components,
    _line_design,
    _OFFSET_KNOTS,
    evaluate_lsf_diagnostics,
    mspline_basis,
    native_pixel_edges,
)
from skysub.sky_decomp.lsf_surface_iterative import (
    LSFChannelSplineConfig,
    _channel_mask,
    _configured_knot_vector,
    _fit_bspline_channel,
    evaluate_bspline_basis,
)
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT

GOLDEN = Path(__file__).parent / "data/lsf_surface_iterative_row837_n5_golden.npz"
SOLAR = DEFAULT_DATA_ROOT / "Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt"
CENTERS = 0.5 * (_OFFSET_KNOTS[:-4] + _OFFSET_KNOTS[4:])


def _uniform_knots(wave, n_basis=4, degree=3):
    spline = LSFChannelSplineConfig(n_basis=n_basis, degree=degree, knot_strategy="uniform")
    return _configured_knot_vector(wave, spline)


def _synthetic_solar(path, scale=None):
    """Write a Meftah-format file (vacuum nm) covering the grid, with Fraunhofer-like lines."""
    nm = np.arange(399.80, 403.60, 0.0021)
    air = vac_to_air(nm * 10.0)
    flux = 1.0 + 0.1 * np.sin(air / 7.0)
    for centre, depth, width in ((4005.3, 0.8, 0.15), (4012.9, 0.5, 0.3), (4021.1, 0.9, 0.1)):
        flux -= depth * np.exp(-0.5 * ((air - centre) / width) ** 2)
    if scale is not None:
        flux = flux * scale(air)
    np.savetxt(path, np.column_stack([nm, flux]), header="; synthetic", comments="")
    return air, flux


def _theta(sigmas):
    """Unimodal, slightly asymmetric kernel per wavelength basis; every column sums to one."""
    kernel = np.exp(-0.5 * ((CENTERS[:, None] - 0.05) / np.asarray(sigmas)[None, :]) ** 2)
    kernel *= 1.0 + 0.3 * np.tanh(CENTERS)[:, None]
    return kernel / kernel.sum(axis=0)


@pytest.mark.parametrize("strength", [None, 0., 0.1])
def test_constructor_defaults_and_regularization_override(monkeypatch, strength):
    captured = {}
    def initialize(self, wave, *args, **kwargs):
        captured.update(kwargs)
        self.wave = wave
        self.config = kwargs["config"]
    monkeypatch.setattr(solar_lsf.SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D,
        "__init__", initialize)
    options = {} if strength is None else {"offset_roughness_fraction": strength}
    solar_lsf.SkyDecompPalaceCorrSolarBlueLSF(np.linspace(3600., 9800., 20), **options)
    config, splines = captured["config"], captured["spline_config"]
    assert config.n_basis == 4 and config.n_refinement_cycles == 30
    assert [splines.for_channel(c).n_basis for c in ("B", "R", "Z")] == [4, 4, 4]
    assert config.roughness_fraction == config.fallback_prior_fraction == 0.
    assert captured["offset_half_width_angstrom"] == 2.5
    assert captured["offset_roughness_fraction"] == (0.01 if strength is None else strength)


def test_profile_light_widths_and_asymmetry():
    wave = np.array([4000., 4200.])
    knots = _uniform_knots(wave)
    coefficient = np.zeros((11, 4))
    coefficient[5] = 1.0
    state = SimpleNamespace(coefficients={"B": coefficient}, knot_vectors={"B": knots}, degrees={"B": 3})
    metric = evaluate_lsf_diagnostics(state, wave)
    np.testing.assert_allclose(metric["centroid_angstrom"], 0., atol=1e-14)
    np.testing.assert_allclose(metric["peak_angstrom"], 0., atol=1e-14)
    np.testing.assert_allclose(metric["asymmetry"], 0., atol=1e-12)
    assert np.all(metric["w90_angstrom"] > metric["w50_angstrom"])
    # Independent dense trapezoid CDF checks the equal-tail definitions.
    delta = np.linspace(-3., 3., 24001)
    profile = mspline_basis(delta)[:, 5]
    cdf = np.r_[0., np.cumsum((profile[:-1] + profile[1:]) * np.diff(delta) / 2)]
    q = np.interp([.05, .25, .75, .95], cdf, delta)
    np.testing.assert_allclose(metric["w50_angstrom"], q[2]-q[1], atol=5e-5)
    np.testing.assert_allclose(metric["w90_angstrom"], q[3]-q[0], atol=5e-5)


def test_offset_only_penalty_reduces_coefficient_curvature():
    rng = np.random.default_rng(27)
    wave = np.linspace(4000., 4100., 500)
    design = rng.normal(size=(500, 44))
    truth = np.zeros((11, 4))
    truth[4:7] = np.array([.1, .8, .1])[:, None]
    target = design @ truth.ravel()
    solutions = []
    for strength in (0., 1.):
        _, coefficient, _, metric = _fit_bspline_channel(wave, np.ones(500), target,
            np.ones(500), truth[:, 0], n_basis=4, degree=3, knot_vector=_uniform_knots(wave),
            kernel_design=design, roughness_fraction=0., fallback_prior_fraction=0.,
            offset_roughness_fraction=strength, numerical_ridge=0.)
        assert metric["status"] in {"Solved", "AlmostSolved"}
        np.testing.assert_allclose(coefficient.sum(axis=0), 1., atol=1e-12)
        solutions.append(coefficient)
    roughness = [np.sum(np.diff(c, n=2, axis=0)**2) for c in solutions]
    assert roughness[1] < roughness[0] * .9


def _direct_g(wave, knots, air, weight, pixels, count=11, half_width=3.0):
    """Independent integral: Gauss-Legendre on every piece of the piecewise-cubic M-splines."""
    edges = native_pixel_edges(wave)
    nodes, gl_weight = np.polynomial.legendre.leggauss(3)
    out = np.zeros((len(pixels), count, knots.size - 4))
    for row, pixel in enumerate(pixels):
        lo, hi = edges[pixel], edges[pixel + 1]
        for index in np.flatnonzero(np.abs(air - wave[pixel]) < 4.0):
            cuts = np.unique(np.r_[lo, hi, np.clip(air[index] + np.linspace(-half_width, half_width, count+4), lo, hi)])
            half, mid = 0.5 * np.diff(cuts), 0.5 * (cuts[1:] + cuts[:-1])
            x = (mid[:, None] + half[:, None] * nodes[None, :]).ravel()
            w = (half[:, None] * gl_weight[None, :]).ravel()
            kernel = w @ mspline_basis(x - air[index], n_basis=count, half_width=half_width) / (hi - lo)
            basis = evaluate_bspline_basis(
                np.clip(air[index : index + 1], knots[3], knots[-4]), knots, 3
            )[0]
            out[row] += weight[index] * kernel[:, None] * basis[None, :]
    return out


@pytest.fixture(scope="module")
def synthetic_grid(tmp_path_factory):
    path = tmp_path_factory.mktemp("solar") / "solar.txt"
    air, flux = _synthetic_solar(path)
    wave = np.linspace(4000.0, 4030.0, 61) + 0.01 * np.sin(np.arange(61))
    return path, wave, air, flux


@pytest.mark.parametrize("count, half_width", [(11, 3.0), (15, 2.5)])
def test_solar_tensor_matches_direct_double_integral(synthetic_grid, count, half_width):
    path, wave, air, flux = synthetic_grid
    knots = _uniform_knots(wave)
    tensor = solar_lsf._solar_tensor(wave, knots, 3, str(path), chunk=700, offset_basis_count=count, offset_half_width_angstrom=half_width)
    weight = flux / np.nanmedian(flux) * np.diff(native_pixel_edges(air))
    pixels = [3, 30, 57]
    direct = _direct_g(wave, knots, air, weight, pixels, count, half_width)
    error = np.max(np.abs(tensor[pixels] - direct)) / np.max(np.abs(direct))
    print(f"G vs direct integral, max relative error {error:.2e}")
    assert error < 1.0e-10
    again = solar_lsf.solar_lsf_tensor(wave, knots, 3, path, count, half_width)
    assert solar_lsf.solar_lsf_tensor(wave, knots, 3, path, count, half_width) is again
    np.testing.assert_allclose(again, tensor, rtol=1e-12, atol=1e-12)


def test_center_offsets_gives_the_same_constrained_solution():
    rng = np.random.default_rng(3)
    wave = np.linspace(4000.0, 4100.0, 400)
    knots = _uniform_knots(wave, 2, 1)
    basis = evaluate_bspline_basis(wave, knots, 1)
    truth = _theta([0.4, 0.55])
    # Common term plus an offset-dependent part.  A much larger common term only degrades
    # the uncentered solve (that is why the centering exists), so the equivalence is checked here.
    design = 1.0 + rng.normal(0.0, 3.0, (wave.size, 11, 1)) * (1.0 + basis[:, None, :])
    target = np.einsum("pkb,kb->p", design, truth) + rng.normal(0.0, 0.05, wave.size)
    ivar = np.full(wave.size, 400.0)
    fallback = truth.mean(axis=1)
    common = dict(
        n_basis=2, degree=1, knot_vector=knots, roughness_fraction=0.0,
        offset_roughness_fraction=0.0, fallback_prior_fraction=1.0e-8, background_degree=0,
    )
    flat = design.reshape(wave.size, -1)
    plain = _fit_bspline_channel(wave, np.ones(wave.size), target, ivar, fallback, kernel_design=flat, **common)
    centered = _fit_bspline_channel(
        wave, np.ones(wave.size), target, ivar, fallback, kernel_design=flat, center_offsets=True, **common
    )
    assert plain[3]["status"] == centered[3]["status"] == "Solved"
    error = np.max(np.abs(plain[1] - centered[1]))
    print(f"center_offsets coefficient difference {error:.2e}")
    np.testing.assert_allclose(centered[1], plain[1], rtol=1e-6, atol=2e-6)
    with pytest.raises(ValueError, match="center_offsets"):
        _fit_bspline_channel(wave, np.ones(wave.size), target, ivar, fallback, center_offsets=True, **common)


def test_envelope_at_output_pixel_approximates_source_envelope(synthetic_grid, tmp_path):
    path, wave, air, flux = synthetic_grid
    knots = _uniform_knots(wave)
    # Smooth envelope: lambda^-4 times a cubic B-spline bump of 100 A knot spacing
    # (production Moon/Zodi knots are several times wider).
    bump = BSpline.basis_element(np.linspace(3815.0, 4215.0, 5), extrapolate=False)
    envelope = lambda x: (x / 4015.0) ** -4 * np.nan_to_num(bump(x))
    shifted = tmp_path / "enveloped.txt"
    _, enveloped = _synthetic_solar(shifted, scale=envelope)
    plain = solar_lsf._solar_tensor(wave, knots, 3, str(path))
    exact = solar_lsf._solar_tensor(wave, knots, 3, str(shifted))
    # Each file is normalised by its own median; restore the common unit.
    exact *= np.nanmedian(enveloped) / np.nanmedian(flux)
    theta = _theta([0.4, 0.5, 0.6, 0.7])
    approx = envelope(wave)[:, None] * np.einsum("pkb,kb->p", plain, theta)[:, None]
    exact = np.einsum("pkb,kb->p", exact, theta)[:, None]
    bright = np.abs(exact[:, 0]) > 0.5 * np.max(np.abs(exact))
    error = np.max(np.abs(approx - exact)[bright] / np.abs(exact)[bright])
    print(f"envelope at output pixel vs source wavelength, max relative error {error:.2e} ({bright.sum()} bright pixels)")
    assert bright.sum() > 10 and error < 1.0e-3


@pytest.fixture(scope="module")
def real_blue():
    if not SOLAR.is_file() or not GOLDEN.is_file():
        pytest.skip("solar source or golden grid not available")
    wave = np.asarray(np.load(GOLDEN)["wave"], dtype=float)
    blue = _channel_mask(wave, None, LSF_CHANNELS[0][2])
    knots = _uniform_knots(wave[blue])
    return wave, blue, knots, solar_lsf.solar_lsf_tensor(wave, knots, 3, SOLAR)


def _recovered_fwhm(real_blue, continuum_level, seed, noise_fraction, regularisation):
    wave, blue, knots, tensor = real_blue
    rng = np.random.default_rng(seed)
    truth = _theta([0.40, 0.46, 0.52, 0.58])
    line = _integrated_components(wave, np.array([5577.34]), blue)
    basis = evaluate_bspline_basis(np.array([5577.34]), knots, 3)
    line_design = _line_design(line, np.array([5.0e3]), basis)[blue].toarray()
    extra = continuum_level * tensor.reshape(tensor.shape[0], -1)
    clean = (line_design + extra) @ truth.reshape(-1)
    sigma = 2.0 if noise_fraction is None else noise_fraction * np.median(clean)
    flux = clean + rng.normal(0.0, sigma, clean.size)
    _, theta, _, metric = _fit_bspline_channel(
        wave[blue], np.ones(blue.sum()), flux, np.full(blue.sum(), sigma**-2), _theta([0.5] * 4)[:, 0],
        n_basis=4, degree=3, knot_vector=knots, kernel_design=line_design + extra,
        center_offsets=True, roughness_fraction=regularisation, offset_roughness_fraction=regularisation,
        fallback_prior_fraction=regularisation, background_degree=3,
    )
    probe = np.arange(3700.0, 5701.0, 100.0)
    states = [SimpleNamespace(coefficients={"B": value}, knot_vectors={"B": knots}, degrees={"B": 3}) for value in (truth, theta)]
    fwhm_true, fwhm = (evaluate_lsf_diagnostics(state, probe)["fwhm_angstrom"] for state in states)
    return metric, theta, fwhm / fwhm_true - 1.0


def test_synthetic_recovery_bright_solar_continuum(real_blue):
    # Weak priors isolate the joint design; the production priors (1e-4) pull a
    # kernel 20% narrower than the fallback wider by several percent at the blue edge.
    metric, theta, ratio = _recovered_fwhm(real_blue, 100.0, 1, 0.003, 1.0e-6)
    print(f"bright continuum: FWHM recovery max |error| {np.max(np.abs(ratio)):.3%}, status {metric['status']}")
    assert metric["status"] in {"Solved", "AlmostSolved"}
    np.testing.assert_allclose(theta.sum(axis=0), 1.0, atol=1e-12)
    assert np.max(np.abs(ratio)) < 0.02


def test_synthetic_recovery_dark_sky_is_finite_and_valid(real_blue):
    metric, theta, ratio = _recovered_fwhm(real_blue, 1.0e-9, 2, None, 1.0e-4)
    print(f"dark continuum: status {metric['status']} reason {metric['reason']!r}, FWHM error range {np.nanmin(ratio):.3%}..{np.nanmax(ratio):.3%}")
    assert np.all(np.isfinite(theta)) and np.all(theta >= 0.0)
    np.testing.assert_allclose(theta.sum(axis=0), 1.0, atol=1e-12)
    assert np.all(np.isfinite(ratio))


class _FakeTelluric:
    def match_to_data(self, wave, lsf, pwv, *, airmass, lsf_in_wavelength):
        return np.full_like(wave, 0.85, dtype=float)


def test_integration_smoke_fits_one_spectrum(tmp_path, monkeypatch):
    if not SOLAR.is_file() or not GOLDEN.is_file():
        pytest.skip("solar source or golden spectrum not available")
    from skysub import decompose_parallel as dp
    from skysub.sky_decomp.result_io import load_lsf_surface_state

    data = np.load(GOLDEN)
    wave = np.asarray(data["wave"], dtype=float)
    model = solar_lsf.SkyDecompPalaceCorrSolarBlueLSF(
        wave, telluric_calculator=_FakeTelluric(), pwv_mm=1.0, source_airmass=1.2,
        drp_transmission=np.full(wave.size, 0.9), lsf_sigma=0.5, base_dir=DEFAULT_DATA_ROOT,
        moon_smooth_lambda=0.1, moon_interline_boost=0.0,
    )
    result = model.fit(np.asarray(data["flux"], dtype=float), np.asarray(data["ivar"], dtype=float))
    state = result.lsf_state
    metric = state.metrics["B"]
    print(f"smoke: status {result.fit_status}, B metric {metric.get('model')} info fraction {metric.get('extra_information_fraction')}")
    assert all(c.shape == (11, 4) for c in state.coefficients.values())
    assert metric["model"] == "solar_continuum_joint_mspline"
    assert 0.0 <= metric["extra_information_fraction"] <= 1.0
    assert "blue_lsf=solar_continuum_joint_n4" in result.fit_summary
    np.testing.assert_allclose(state.coefficients["B"].sum(axis=0), 1.0, atol=2.0e-12)
    assert state.config["offset_roughness_fraction"] == 0.01
    assert state.config["roughness_fraction"] == state.config["fallback_prior_fraction"] == 0.
    assert model.converged and 2 <= state.completed_cycles < 30
    assert len(model.convergence_history) == state.completed_cycles
    assert not model.convergence_history[0]["converged"]
    assert model.convergence_history[-1]["converged"]
    assert model.convergence_history[-1]["centroid_change_angstrom"] < 0.015
    assert all(not m["converged"] for m in model.convergence_history[:-1])
    # Exercise the same cache and compact FITS writer used by server workers.
    monkeypatch.setattr(dp, "_WORKER_COMPACT_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(dp, "_WORKER_RUN_FINGERPRINT", "solar-defaults")
    dp._save_compact_cache("sky1", 0, result=result)
    output = tmp_path / "solar.fits"
    dp._write_compact_fits(tmp_path / "cache", "sky1", 1, "solar-defaults",
        dp.PALACECORR_SOLARBLUE_FIT_MODEL, output)
    saved = load_lsf_surface_state(output, 0)
    assert saved.config["offset_roughness_fraction"] == 0.01
    assert saved.config["offset_half_width_angstrom"] == 2.5
    assert saved.completed_cycles == state.completed_cycles
    assert all(m["converged"] for m in saved.metrics.values())
    assert all(m["convergence_centroid_atol_angstrom"] == 0.015 for m in saved.metrics.values())
    for channel in ("B", "R", "Z"):
        np.testing.assert_array_equal(saved.coefficients[channel], state.coefficients[channel])
    assert "convergence=converged" in result.fit_summary
    np.testing.assert_allclose(sum(result.components[k] for k in
        ("moon", "zodi", "diffuse", "oh", "atom", "orc", "o2")), result.bestfit_lsf, atol=1e-12)


@pytest.mark.parametrize("value", [0.0, -0.01, np.nan])
def test_invalid_convergence_tolerance_is_rejected(value):
    with pytest.raises(ValueError, match="convergence tolerances"):
        solar_lsf.SkyDecompPalaceCorrSolarBlueLSF(np.linspace(3600, 9800, 20), lsf_convergence_rtol=value)
