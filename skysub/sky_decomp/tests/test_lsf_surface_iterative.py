import numpy as np
import pytest
from astropy.io import fits
from dataclasses import fields

from skysub.sky_decomp.lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    LSFSurfaceIterativeResult,
    apply_lsf_channelwise,
    apply_lsf_surface,
    build_bspline_basis,
    build_skyline_mask,
    continuum_fit_weights,
    evaluate_lsf_kernel,
    evaluate_lsf_surface,
    fit_bspline_channel,
    fit_lsf_surface,
    kernel_moments,
    load_lsf_surface_state,
)
from skysub.sky_decomp.fit import SkyDecompResult
from skysub.decompose_parallel import extract_meta_and_coef_products, results_to_fits


def _smooth_test_surface(wave, kernel_size=11, n_basis=6):
    basis, _ = build_bspline_basis(
        wave,
        n_basis,
        3,
        np.ones_like(wave),
    )
    offsets = np.arange(-(kernel_size // 2), kernel_size // 2 + 1)
    coefficient = np.empty((kernel_size, n_basis))
    for column in range(n_basis):
        center = -0.15 + 0.30 * column / (n_basis - 1)
        sigma = 0.85 + 0.35 * column / (n_basis - 1)
        kernel = np.exp(-0.5 * ((offsets - center) / sigma) ** 2)
        coefficient[:, column] = kernel / np.sum(kernel)
    return basis @ coefficient.T


def _baseline_result_kwargs(wave):
    return {
        "t_o2": 191.5,
        "t_o2_err": 1.0,
        "reduced_chi2": 1.0,
        "r2": 0.9,
        "rms_resid": 0.1,
        "resid_level": -0.3,
        "fit_status": "Solved",
        "fit_summary": "test",
        "fit_elapsed_sec": 1.0,
        "peak_memory_mb": 1.0,
        "o2_fit_status": "Solved",
        "o2_fit_summary": "test",
        "o2_fit_elapsed_sec": 0.1,
        "o2_valid_frac": 1.0,
        "coef": np.array([1.0]),
        "coef_err": np.array([np.nan]),
        "design_names": ["component"],
        "bestfit": np.zeros_like(wave),
        "bestfit_lsf": np.zeros_like(wave),
        "resid": np.zeros_like(wave),
        "components": {"oh": np.zeros_like(wave)},
        "lsf_kernels": {},
        "lsf_metrics": {},
        "moon_knots": np.array([], dtype=float),
        "moon_boosted_pixels": np.array([], dtype=float),
        "vector_o2": np.zeros_like(wave),
        "o2_prefit_amp": 0.0,
        "bestfit_lsf_sigma": np.zeros_like(wave),
    }


def _apply_lsf_surface_index_reference(
    wave,
    values,
    kernel_surface,
    channel_bounds,
    tap_offsets,
):
    source = np.atleast_2d(np.asarray(values, dtype=float))
    output = np.zeros_like(source)
    for channel in ("B", "R", "Z"):
        lower, upper = channel_bounds[channel]
        mask = np.ones(wave.size, dtype=bool)
        if lower is not None:
            mask &= wave >= lower
        if upper is not None:
            mask &= wave < upper
        indices = np.flatnonzero(mask)
        for tap, offset in enumerate(tap_offsets):
            if offset >= 0:
                target_index = indices[offset:]
                source_index = indices[: indices.size - offset]
            else:
                target_index = indices[:offset]
                source_index = indices[-offset:]
            if target_index.size:
                output[:, target_index] += (
                    source[:, source_index]
                    * kernel_surface[target_index, tap][None, :]
                )
    return output


def test_fit_recovers_a_smooth_normalized_kernel_surface():
    wave = np.linspace(5800.0, 7453.0, 1200)
    source = np.zeros_like(wave)
    source[np.arange(30, wave.size - 30, 19)] = np.linspace(0.5, 1.5, 60)
    true_surface = _smooth_test_surface(wave)
    target = apply_lsf_surface(wave, source, true_surface)

    fitted_surface, _, _, metrics = fit_bspline_channel(
        wave,
        source,
        target,
        np.ones_like(wave),
        np.median(true_surface, axis=0),
        roughness_fraction=0.0,
        fallback_prior_fraction=0.0,
    )
    fitted = apply_lsf_surface(wave, source, fitted_surface)

    assert metrics["status"] in {"Solved", "AlmostSolved"}
    assert np.min(fitted_surface) >= 0.0
    assert np.allclose(np.sum(fitted_surface, axis=1), 1.0, atol=1.0e-12)
    assert np.sqrt(np.mean((fitted - target) ** 2)) < 1.0e-4


def test_constraints_keep_both_wings_single_peaked():
    wave = np.linspace(7454.0, 9800.0, 900)
    source = np.zeros_like(wave)
    source[25:-25:23] = 1.0
    true_surface = _smooth_test_surface(wave)
    target = apply_lsf_surface(wave, source, true_surface)
    fitted_surface, _, _, _ = fit_bspline_channel(
        wave,
        source,
        target,
        np.ones_like(wave),
        np.median(true_surface, axis=0),
    )

    center = fitted_surface.shape[1] // 2
    assert np.all(np.diff(fitted_surface[:, : center + 1], axis=1) >= -1.0e-10)
    assert np.all(np.diff(fitted_surface[:, center:], axis=1) <= 1.0e-10)
    centroid, sigma = kernel_moments(fitted_surface)
    assert np.all(np.isfinite(centroid))
    assert np.all(sigma > 0.0)


def test_surface_application_does_not_mix_adjacent_channels():
    wave = np.arange(5783.0, 5790.0)
    values = np.zeros(wave.size)
    values[1] = 1.0
    values[3] = 3.0
    surface = np.zeros((wave.size, 3))
    surface[:, 2] = 1.0

    output = apply_lsf_surface(wave, values, surface)

    # The first impulse shifts within B.  The second would cross B -> R and is
    # therefore absent from the result.
    assert np.array_equal(
        output,
        np.array([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
    )


def test_surface_application_matches_index_reference_for_matrix_input():
    rng = np.random.default_rng(9274)
    wave = np.linspace(5000.0, 5012.0, 25)
    values = rng.normal(size=(7, wave.size))
    surface = rng.random(size=(wave.size, 5))
    surface /= np.sum(surface, axis=1, keepdims=True)
    channel_bounds = {
        "B": (None, 5004.0),
        "R": (5004.0, 5008.0),
        "Z": (5008.0, None),
    }
    tap_offsets = np.arange(-2, 3)

    expected = _apply_lsf_surface_index_reference(
        wave,
        values,
        surface,
        channel_bounds,
        tap_offsets,
    )
    actual = apply_lsf_surface(
        wave,
        values,
        surface,
        channel_bounds=channel_bounds,
        tap_offsets=tap_offsets,
    )

    # Sparse mat-mul reorders summations vs. the per-tap index reference; agreement is exact to float64 rounding.
    assert np.allclose(actual, expected, atol=1e-14, rtol=1e-14)


def test_blue_fit_window_still_builds_a_full_channel_surface():
    wave = np.arange(5400.0, 7601.0)
    source = np.zeros_like(wave)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.1) ** 2)
    kernel /= np.sum(kernel)

    state = fit_lsf_surface(
        wave,
        np.zeros_like(wave),
        np.ones_like(wave),
        source,
        np.zeros_like(wave),
        {"B": kernel, "R": kernel, "Z": kernel},
        LSFSurfaceIterativeConfig(),
    )
    reconstructed = evaluate_lsf_surface(state, wave)

    assert np.all(np.isfinite(reconstructed))
    assert np.allclose(np.sum(reconstructed, axis=1), 1.0, atol=1.0e-12)
    assert not hasattr(state, "wave")
    assert not hasattr(state, "kernel_surface")


def test_public_config_defaults_to_five_configurable_cycles():
    assert LSFSurfaceIterativeConfig().n_refinement_cycles == 5
    assert LSFSurfaceIterativeConfig(n_refinement_cycles=2).n_refinement_cycles == 2
    with pytest.raises(ValueError, match="positive"):
        LSFSurfaceIterativeConfig(n_refinement_cycles=0)
    with pytest.raises(TypeError, match="integer"):
        LSFSurfaceIterativeConfig(n_refinement_cycles=True)


def test_scalar_evaluation_and_state_application_use_compact_state():
    wave = np.arange(5400.0, 8801.0)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.1) ** 2)
    kernel /= np.sum(kernel)
    source = np.zeros_like(wave)
    source[100:-100:31] = 1.0
    state = fit_lsf_surface(
        wave,
        apply_lsf_surface(wave, source, np.tile(kernel, (wave.size, 1))),
        np.ones_like(wave),
        source,
        np.zeros_like(wave),
        {channel: kernel for channel in ("B", "R", "Z")},
        LSFSurfaceIterativeConfig(),
    )

    scalar_kernel = evaluate_lsf_kernel(state, 6500.0)
    assert scalar_kernel.shape == (11,)
    assert np.isclose(np.sum(scalar_kernel), 1.0)
    assert np.allclose(
        apply_lsf_channelwise(state, wave, source),
        apply_lsf_surface(wave, source, evaluate_lsf_surface(state, wave)),
    )
    with pytest.raises(ValueError, match="fitted pixel grid"):
        apply_lsf_channelwise(state, wave + 0.01, source)


def test_skyline_mask_keeps_strong_group_lines_and_expands_in_wavelength():
    wave = np.arange(5000.0, 5011.0)
    sticks = np.zeros((2, wave.size))
    sticks[0, 2] = 9.0
    sticks[0, 8] = 1.0
    sticks[1, 6] = 2.0

    mask = build_skyline_mask(
        wave,
        (sticks,),
        cumulative_fraction=0.8,
        half_width_angstrom=1.1,
    )

    assert np.array_equal(np.flatnonzero(mask), np.array([1, 2, 3, 5, 6, 7]))


def test_continuum_weights_suppress_skyline_pixels_and_huber_outliers():
    wave = np.arange(5000.0, 5010.0)
    residual = np.zeros_like(wave)
    residual[7] = 10.0
    skyline_mask = np.zeros(wave.size, dtype=bool)
    skyline_mask[4] = True

    weights, channel_noise = continuum_fit_weights(
        wave,
        residual,
        np.ones_like(wave),
        skyline_mask,
        line_weight=5.0e-4,
        huber_transition_sigma=3.0,
    )

    assert np.isclose(weights[4] / weights[0], 5.0e-4)
    assert weights[7] < weights[0]
    assert channel_noise["B"] == 1.0


def test_later_lsf_cycles_keep_first_cycle_knots_fixed():
    wave = np.arange(5400.0, 8801.0)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.1) ** 2)
    kernel /= np.sum(kernel)
    source_first = np.zeros_like(wave)
    source_first[100:-100:31] = 1.0
    source_second = np.zeros_like(wave)
    source_second[110:-100:47] = 2.0
    fallback = {channel: kernel for channel in ("B", "R", "Z")}

    first = fit_lsf_surface(
        wave,
        apply_lsf_surface(wave, source_first, np.tile(kernel, (wave.size, 1))),
        np.ones_like(wave),
        source_first,
        np.zeros_like(wave),
        fallback,
        LSFSurfaceIterativeConfig(),
    )
    second = fit_lsf_surface(
        wave,
        apply_lsf_surface(wave, source_second, np.tile(kernel, (wave.size, 1))),
        np.ones_like(wave),
        source_second,
        np.zeros_like(wave),
        fallback,
        LSFSurfaceIterativeConfig(),
        previous_state=first,
    )

    for channel in ("B", "R", "Z"):
        assert np.array_equal(second.knot_vectors[channel], first.knot_vectors[channel])
        assert second.metrics[channel]["knots_fixed"] is True


def test_later_channel_failure_preserves_previous_surface():
    wave = np.arange(5400.0, 8801.0)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.1) ** 2)
    kernel /= np.sum(kernel)
    source = np.zeros_like(wave)
    source[100:-100:31] = 1.0
    fallback = {channel: kernel for channel in ("B", "R", "Z")}
    first = fit_lsf_surface(
        wave,
        apply_lsf_surface(wave, source, np.tile(kernel, (wave.size, 1))),
        np.ones_like(wave),
        source,
        np.zeros_like(wave),
        fallback,
        LSFSurfaceIterativeConfig(),
    )

    second = fit_lsf_surface(
        wave,
        np.zeros_like(wave),
        np.ones_like(wave),
        np.zeros_like(wave),
        np.zeros_like(wave),
        fallback,
        LSFSurfaceIterativeConfig(),
        previous_state=first,
    )

    assert np.array_equal(
        evaluate_lsf_surface(second, wave),
        evaluate_lsf_surface(first, wave),
    )
    for channel in ("B", "R", "Z"):
        assert second.metrics[channel]["status"] == "fallback"
        assert second.metrics[channel]["reason"].startswith("previous_surface:")


def test_extended_result_only_adds_lsf_state_to_baseline_contract():
    baseline_fields = [field.name for field in fields(SkyDecompResult)]
    extended_fields = [field.name for field in fields(LSFSurfaceIterativeResult)]

    assert extended_fields == baseline_fields + ["lsf_state"]


def test_lsf_fits_round_trip_reconstructs_surface(tmp_path):
    wave = np.arange(5400.0, 8801.0)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.1) ** 2)
    kernel /= np.sum(kernel)
    source = np.zeros_like(wave)
    source[100:-100:31] = 1.0
    surface = np.tile(kernel, (wave.size, 1))
    state = fit_lsf_surface(
        wave,
        apply_lsf_surface(wave, source, surface),
        np.ones_like(wave),
        source,
        np.zeros_like(wave),
        {channel: kernel for channel in ("B", "R", "Z")},
        LSFSurfaceIterativeConfig(),
    )
    state.completed_cycles = 4
    state.fit_status = "Solved"
    result = LSFSurfaceIterativeResult(
        **_baseline_result_kwargs(wave),
        lsf_state=state,
    )
    output = tmp_path / "extended.fits"

    results_to_fits([result], output)
    restored = load_lsf_surface_state(output)

    with fits.open(output) as hdul:
        assert all(name in hdul for name in ("LSF_COEF", "LSF_KNOTS", "LSF_META"))
        assert "LSF_WAVE" not in hdul
    assert restored.completed_cycles == 4
    assert np.allclose(
        evaluate_lsf_surface(restored, wave),
        evaluate_lsf_surface(state, wave),
        atol=1.0e-10,
    )

    input_path = tmp_path / "input.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU.from_columns([], nrows=1, name="META"),
        ]
    ).writeto(input_path)
    compact_paths = [tmp_path / f"compact_{index}.fits" for index in range(3)]
    extract_meta_and_coef_products(
        input_path,
        output,
        output,
        output,
        meta_output_path=tmp_path / "meta.fits",
        sky1_output_path=compact_paths[0],
        sky2_output_path=compact_paths[1],
        sci_output_path=compact_paths[2],
    )
    with fits.open(compact_paths[0]) as hdul:
        assert [hdu.name for hdu in hdul] == [
            "PRIMARY",
            "META",
            "COEF",
            "COEF_ERR",
            "LSF_COEF",
            "LSF_KNOTS",
            "LSF_META",
        ]


def test_baseline_writer_contract_is_unchanged_and_mixed_results_are_rejected(
    tmp_path,
):
    wave = np.arange(12.0)
    baseline = SkyDecompResult(**_baseline_result_kwargs(wave))
    baseline_path = tmp_path / "baseline.fits"
    results_to_fits([baseline], baseline_path)

    with fits.open(baseline_path) as hdul:
        assert [hdu.name for hdu in hdul] == [
            "PRIMARY",
            "META",
            "COEF",
            "COEF_ERR",
            "BESTFIT",
            "BESTFIT_LSF",
            "FLUX_SIGMA_TOTAL",
            "RESID",
            "VECTOR_O2",
            "COMP_OH",
        ]

    iterative = LSFSurfaceIterativeResult(
        **_baseline_result_kwargs(wave),
        lsf_state=None,
    )
    with pytest.raises(ValueError, match="Cannot mix"):
        results_to_fits([baseline, iterative], tmp_path / "mixed.fits")
