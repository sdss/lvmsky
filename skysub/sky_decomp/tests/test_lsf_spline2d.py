import copy
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from scipy.integrate import quad

from skysub.sky_decomp.lsf_spline2d import (
    LSF_SPLINE2D_REPRESENTATION,
    SkyDecompLSFSpline2D,
    _integrated_components,
    evaluate_lsf_density,
    mspline_basis,
    native_pixel_edges,
)
from skysub.sky_decomp.lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    evaluate_lsf_surface,
)
from skysub.sky_decomp.result_io import load_lsf_surface_state, results_to_fits


@pytest.mark.parametrize("strength", [-0.001, np.nan, np.inf])
def test_invalid_offset_regularization_is_rejected(strength):
    with pytest.raises(ValueError, match="offset_roughness_fraction"):
        SkyDecompLSFSpline2D(np.array([4000., 4001.]), offset_roughness_fraction=strength)


@pytest.mark.parametrize("count, half_width", [(11, 3.0), (15, 2.5)])
def test_mspline_contract_and_exact_native_bins(count, half_width):
    endpoints = np.array([-half_width, half_width])
    options = dict(n_basis=count, half_width=half_width)
    assert np.array_equal(mspline_basis(endpoints, **options), np.zeros((2, count)))
    assert np.array_equal(mspline_basis(endpoints, derivative=1, **options), np.zeros((2, count)))
    for index in range(count):
        integral = quad(
            lambda value: mspline_basis(np.array([value]), **options)[0, index],
            -half_width,
            half_width,
            points=np.linspace(-half_width, half_width, count+4),
        )[0]
        assert integral == pytest.approx(1.0, abs=2.0e-12)

    wave = np.array([4996.0, 4996.8, 4997.7, 4998.5, 4999.4, 5000.2, 5001.1, 5002.0, 5003.0, 5004.1])
    design = _integrated_components(wave, np.array([5000.0]), np.ones(wave.size, bool), count, half_width)
    widths = np.diff(native_pixel_edges(wave))
    masses = widths @ design.toarray()
    np.testing.assert_allclose(masses, np.ones(count), rtol=0.0, atol=3.0e-14)


@pytest.fixture(scope="module")
def real_result():
    root = Path(__file__).parents[1]
    data = np.load(root / "tests/data/lsf_surface_iterative_row837_n5_golden.npz")
    wave = np.asarray(data["wave"], dtype=np.float64)
    flux = np.asarray(data["flux"], dtype=np.float64)
    ivar = np.asarray(data["ivar"], dtype=np.float64)
    model = SkyDecompLSFSpline2D(
        wave,
        lsf_sigma=0.5,
        base_dir=root / "data",
        n_spline_knots=25,
        n_zodi_spline_knots=3,
        zodi_smooth_lambda=0.1,
        moon_smooth_lambda=0.1,
        config=LSFSurfaceIterativeConfig(
            n_refinement_cycles=1,
            roughness_fraction=1.0e-4,
        ),
    )
    return wave, model, model.fit(flux, ivar)


def test_real_fit_state_result_and_fits(real_result, tmp_path):
    wave, _, result = real_result
    state = result.lsf_state
    assert result.fit_status == "Solved"
    assert state.completed_cycles == 1
    assert state.legacy_kernel_representation == LSF_SPLINE2D_REPRESENTATION
    assert state.coefficients["B"].shape == (11, 1)
    assert state.coefficients["R"].shape == state.coefficients["Z"].shape == (11, 6)
    for coefficient in state.coefficients.values():
        np.testing.assert_allclose(coefficient.sum(axis=0), 1.0, atol=2.0e-12)
    density = evaluate_lsf_density(state, np.array([5600.0, 6500.0, 8500.0]), np.array([-3.0, 3.0]))
    assert np.array_equal(density, np.zeros((3, 2)))
    with pytest.raises(ValueError, match="discrete 11-tap"):
        evaluate_lsf_surface(state, wave)
    reconstructed = sum(
        (result.components[key] for key in ("oh", "moon", "zodi", "diffuse", "atom", "orc", "o2")),
        np.zeros_like(wave),
    )
    np.testing.assert_allclose(reconstructed, result.bestfit_lsf, rtol=0.0, atol=3.0e-12)

    output = tmp_path / "result.fits"
    results_to_fits([result], output)
    restored = load_lsf_surface_state(output)
    assert restored.legacy_kernel_representation == LSF_SPLINE2D_REPRESENTATION
    with fits.open(output) as hdul:
        assert hdul[0].header["DECOMPM"] == "lsf-spline2d-split-zodi"
        assert hdul["LSF_COEF"].header["BASIS"] == "M-spline"
        assert "ZODICORR" not in hdul[0].header

    # Run-level cards (decompose_parallel's ZODICORR) reach the primary header
    # of the full product, not just the compact one.
    tagged = tmp_path / "tagged.fits"
    results_to_fits([result], tagged,
                    primary_meta={"ZODICORR": ("lvm-ecl-2026-09", "anchor correction")})
    with fits.open(tagged) as hdul:
        assert hdul[0].header["ZODICORR"] == "lvm-ecl-2026-09"
        assert hdul[0].header["DECOMPM"] == "lsf-spline2d-split-zodi"


def test_failed_input_keeps_mspline_contract(real_result):
    _, model, solved = real_result
    failed = model.failed_input_result("invalid_airmass: row=833")

    assert tuple(failed.components) == tuple(solved.components)
    assert failed.lsf_state.legacy_kernel_representation == LSF_SPLINE2D_REPRESENTATION
    assert failed.lsf_state.coefficients["B"].shape == (11, 1)
    assert failed.lsf_state.coefficients["R"].shape == (11, 6)
    assert failed.lsf_state.coefficients["Z"].shape == (11, 6)
    assert all(
        np.all(np.isnan(values))
        for values in failed.lsf_state.coefficients.values()
    )


def test_vn_line_amplitude_pca_fits_header_keeps_compact_lsf(real_result, tmp_path):
    _, _, result = real_result
    candidate = copy.deepcopy(result)
    candidate.design_names = [
        name.replace("OH_", "OHVN_", 1) if name.startswith("OH_") else name
        for name in candidate.design_names
    ]
    pca_names = ["LineAmplitudePCA_mean"] + [
        f"LineAmplitudePCA_{index:02d}" for index in range(1, 31)
    ]
    candidate.design_names.extend(pca_names)
    candidate.coef = np.concatenate((candidate.coef, np.zeros(len(pca_names))))
    candidate.coef_err = np.concatenate(
        (candidate.coef_err, np.zeros(len(pca_names)))
    )
    candidate.components["line_amplitude_pca"] = np.zeros_like(candidate.bestfit_lsf)

    output = tmp_path / "vn-pca.fits"
    results_to_fits([candidate], output)

    with fits.open(output) as hdul:
        assert hdul[0].header["DECOMPM"] == (
            "telluric-corrected-lines-palace-aijc-vn-line-amplitude-pca"
        )
        assert hdul[0].header["LAPCAK"] == 30
        assert hdul[0].header["OHGROUP"] == "v_upper,N_upper"
        assert all(name in hdul for name in ("LSF_COEF", "LSF_KNOTS", "LSF_META"))


def test_vnf_pca_fits_header_keeps_compact_lsf(real_result, tmp_path):
    _, _, result = real_result
    candidate = copy.deepcopy(result)
    old_oh = sum(name.startswith("OH_") for name in candidate.design_names)
    keep = np.arange(old_oh, len(candidate.design_names))
    oh_names = ["OHVNFPCAPrep_mean", "OHVNFPCAPrep_001"]
    candidate.design_names = oh_names + [candidate.design_names[index] for index in keep]
    candidate.coef = np.concatenate((np.ones(2), candidate.coef[keep]))
    candidate.coef_err = np.concatenate((np.ones(2), candidate.coef_err[keep]))
    pca_names = ["LineAmplitudePCA_mean", "LineAmplitudePCA_001"]
    candidate.design_names.extend(pca_names)
    candidate.coef = np.concatenate((candidate.coef, np.zeros(2)))
    candidate.coef_err = np.concatenate((candidate.coef_err, np.zeros(2)))
    candidate.components["line_amplitude_pca"] = np.zeros_like(candidate.bestfit_lsf)

    output = tmp_path / "vnf-pca.fits"
    results_to_fits([candidate], output)

    with fits.open(output) as hdul:
        assert hdul[0].header["DECOMPM"] == (
            "telluric-corrected-lines-palace-aijc-vnf-pca-line-pca"
        )
        assert hdul[0].header["OHPCAK"] == 1
        assert hdul[0].header["LAPCAK"] == 1
        assert hdul[0].header["OHGROUP"] == "v_upper,N_upper,F_upper"
        assert all(name in hdul for name in ("LSF_COEF", "LSF_KNOTS", "LSF_META"))


def test_vnf_line_adjoint_pca_fits_header_keeps_compact_lsf(real_result, tmp_path):
    _, _, result = real_result
    candidate = copy.deepcopy(result)
    pca_names = ["LineAdjointPCA_mean", "LineAdjointPCA_001"]
    candidate.design_names.extend(pca_names)
    candidate.coef = np.concatenate((candidate.coef, np.zeros(2)))
    candidate.coef_err = np.concatenate((candidate.coef_err, np.zeros(2)))
    candidate.components["line_adjoint_pca"] = np.zeros_like(candidate.bestfit_lsf)

    output = tmp_path / "vnf-adjoint-pca.fits"
    results_to_fits([candidate], output)

    with fits.open(output) as hdul:
        assert hdul[0].header["DECOMPM"] == (
            "telluric-corrected-lines-palace-aijc-vnf-line-adjoint-pca"
        )
        assert "OHPCAK" not in hdul[0].header
        assert hdul[0].header["LADPCAK"] == 1
        assert hdul[0].header["OHGROUP"] == "v_upper,N_upper,F_upper"
        assert all(name in hdul for name in ("LSF_COEF", "LSF_KNOTS", "LSF_META"))


def test_flexible_offset_support_roundtrip(real_result, tmp_path):
    from skysub.sky_decomp.result_io import build_lsf_hdus
    from skysub.sky_decomp.lsf_surface_iterative import LSFChannelSplineConfig, _configured_knot_vector
    wave, _, result = real_result
    state = copy.deepcopy(result.lsf_state)
    state.tap_offsets = np.arange(-7, 8)
    state.config["offset_half_width_angstrom"] = 2.5
    state.config["offset_roughness_fraction"] = 0.01
    for channel, lo, hi in (("B", 3600, 5787), ("R", 5787, 7454), ("Z", 7454, 9801)):
        state.coefficients[channel] = np.zeros((15, 8))
        state.coefficients[channel][7] = 1.0
        state.knot_vectors[channel] = _configured_knot_vector(wave[(wave>=lo)&(wave<hi)],
            LSFChannelSplineConfig(n_basis=8, degree=3, knot_strategy="uniform"))
        state.degrees[channel] = 3
    path = tmp_path / "flexible.fits"
    fits.HDUList([fits.PrimaryHDU(), *build_lsf_hdus([state])]).writeto(path)
    restored = load_lsf_surface_state(path)
    assert restored.config["offset_half_width_angstrom"] == 2.5
    assert restored.config["offset_roughness_fraction"] == 0.01
    assert restored.coefficients["B"].shape == (15, 8)
    probes = np.array([4000., 6500., 8500.])
    delta = np.linspace(-3., 3., 601)
    before = evaluate_lsf_density(state, probes, delta)
    np.testing.assert_array_equal(evaluate_lsf_density(restored, probes, delta), before)
    assert np.all(before[:, np.abs(delta)>=2.5] == 0)
    incompatible = copy.deepcopy(state)
    incompatible.config["offset_half_width_angstrom"] = 3.0
    with pytest.raises(ValueError, match="same offset support"):
        build_lsf_hdus([state, incompatible])
    incompatible.config["offset_half_width_angstrom"] = 2.5
    incompatible.config["offset_roughness_fraction"] = 0.1
    with pytest.raises(ValueError, match="same offset regularization"):
        build_lsf_hdus([state, incompatible])
