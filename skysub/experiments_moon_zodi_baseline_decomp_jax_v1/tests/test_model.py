from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np

from skysub.experiments_moon_zodi_baseline_decomp_jax_v1 import model


def test_line_subtracted_target_closes_to_observed() -> None:
    observed = np.asarray([3.0, 5.0, 7.0])
    lines = np.asarray([0.5, 4.0, 0.0])
    target = model._line_subtracted_target(observed, lines)
    assert np.array_equal(target + lines, observed)


def test_weighted_loss_suppresses_a_skyline_residual() -> None:
    residual = jnp.asarray([[0.1, 10.0]])
    ordinary = model._weighted_loss(residual, jnp.asarray([[1.0, 1.0]]))
    protected = model._weighted_loss(residual, jnp.asarray([[1.0, 5.0e-4]]))
    assert float(protected[0]) < 0.02 * float(ordinary[0])


def test_weighted_loss_is_jax_differentiable() -> None:
    weights = jnp.asarray([[1.0, 5.0e-4]])

    def objective(value: jax.Array) -> jax.Array:
        return jnp.sum(model._weighted_loss(value[None, :], weights))

    gradient = jax.grad(objective)(jnp.asarray([0.1, 2.0]))
    assert np.all(np.isfinite(np.asarray(gradient)))


def test_zero_weight_removes_a_nonfinite_residual() -> None:
    loss = model._weighted_loss(
        jnp.asarray([[0.1, jnp.nan]]),
        jnp.asarray([[1.0, 0.0]]),
    )
    assert np.all(np.isfinite(np.asarray(loss)))


def test_zero_weight_sanitizes_observed_before_autodiff() -> None:
    observed = model._masked_observed(
        jnp.asarray([[1.0, jnp.nan]]),
        jnp.asarray([[1.0, 0.0]]),
    )
    assert np.array_equal(np.asarray(observed), [[1.0, 0.0]])


def test_full_lvm_bands_are_contiguous() -> None:
    assert model.BANDS[0][1] == 3600.0
    assert model.BANDS[-1][2] == 9800.0
    assert model.BANDS[0][2] == model.BANDS[1][1]
    assert model.BANDS[1][2] == model.BANDS[2][1]


def test_continuum_display_limits_ignore_skyline_outlier() -> None:
    lower, upper = model._continuum_display_limits(
        np.asarray([1.0, 2.0, 3.0, 1000.0]),
        np.asarray([False, False, False, True]),
    )
    assert lower == 0.0
    assert 3.0 < upper < 10.0


def test_residual_rms_ignores_skyline_outlier() -> None:
    rms = model._line_free_residual_rms(
        np.asarray([3.0, 4.0, 1000.0]),
        np.asarray([False, False, True]),
    )
    assert np.isclose(rms, np.sqrt(12.5))


def test_display_nodes_require_five_eligible_native_pixels() -> None:
    wave = np.arange(3600.0, 3620.0, 0.5)
    values = np.arange(wave.size, dtype=float)
    valid = np.ones(wave.size, dtype=bool)
    valid[:36] = False
    centers, nodes = model._display_nodes(wave, values, valid)
    assert centers[0] == 3610.0
    assert np.isnan(nodes[0])
