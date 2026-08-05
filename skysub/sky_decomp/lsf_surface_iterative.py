"""Iterative continuum, line, and wavelength-dependent LSF decomposition.

The blue arm uses one shape-constrained kernel because it has relatively little
isolated line information.  Each red-arm kernel tap is a smooth B-spline of
wavelength.  The quadratic program keeps every evaluated kernel nonnegative,
normalized, and single-peaked at the central tap.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import time
import tracemalloc
from pathlib import Path

import clarabel
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import BSpline

from .fit import LSF_CHANNELS, LSF_KERNEL_SIZE, SkyDecomp, SkyDecompResult


LSF_STATE_SCHEMA_VERSION = 1
LSF_TAP_OFFSETS = np.arange(-(LSF_KERNEL_SIZE // 2), LSF_KERNEL_SIZE // 2 + 1)


@dataclass(frozen=True, slots=True)
class LSFSurfaceIterativeConfig:
    """Configuration for the complete iterative decomposition method."""

    line_weight: float = 5.0e-4
    skyline_cumulative_fraction: float = 0.99
    skyline_half_width_angstrom: float = 2.0
    huber_transition_sigma: float = 3.0
    n_refinement_cycles: int = 4
    n_basis: int = 6
    degree: int = 3
    roughness_fraction: float = 0.01
    fallback_prior_fraction: float = 1.0e-4
    information_prior_max_boost: float = 50.0
    background_degree: int = 3
    blue_fit_lower: float = 5500.0

    def __post_init__(self) -> None:
        if isinstance(self.n_refinement_cycles, (bool, np.bool_)) or not isinstance(
            self.n_refinement_cycles,
            (int, np.integer),
        ):
            raise TypeError("n_refinement_cycles must be an integer")
        if not 0.0 <= self.line_weight <= 1.0:
            raise ValueError("line_weight must lie in [0, 1]")
        if not 0.0 < self.skyline_cumulative_fraction <= 1.0:
            raise ValueError(
                "skyline_cumulative_fraction must lie in (0, 1]"
            )
        if self.skyline_half_width_angstrom <= 0.0:
            raise ValueError("skyline_half_width_angstrom must be positive")
        if self.huber_transition_sigma <= 0.0:
            raise ValueError("huber_transition_sigma must be positive")
        if self.n_refinement_cycles < 1:
            raise ValueError("n_refinement_cycles must be positive")
        if isinstance(self.n_basis, (bool, np.bool_)) or not isinstance(
            self.n_basis,
            (int, np.integer),
        ):
            raise TypeError("n_basis must be an integer")
        if isinstance(self.degree, (bool, np.bool_)) or not isinstance(
            self.degree,
            (int, np.integer),
        ):
            raise TypeError("degree must be an integer")
        if self.degree < 0:
            raise ValueError("degree must be nonnegative")
        if self.n_basis < self.degree + 1:
            raise ValueError("n_basis must be at least degree + 1")
        if self.roughness_fraction < 0.0:
            raise ValueError("roughness_fraction must be nonnegative")
        if self.fallback_prior_fraction < 0.0:
            raise ValueError("fallback_prior_fraction must be nonnegative")
        if self.information_prior_max_boost < 1.0:
            raise ValueError("information_prior_max_boost must be at least one")
        if isinstance(self.background_degree, (bool, np.bool_)) or not isinstance(
            self.background_degree,
            (int, np.integer),
        ):
            raise TypeError("background_degree must be an integer")
        if self.background_degree < 0:
            raise ValueError("background_degree must be nonnegative")
        if not np.isfinite(self.blue_fit_lower):
            raise ValueError("blue_fit_lower must be finite")


@dataclass(slots=True)
class LSFSurfaceState:
    """Compact, self-describing representation of one fitted LSF."""

    coefficients: dict[str, np.ndarray]
    knot_vectors: dict[str, np.ndarray]
    degrees: dict[str, int]
    channel_bounds: dict[str, tuple[float | None, float | None]]
    tap_offsets: np.ndarray
    config: dict[str, float | int]
    metrics: dict[str, dict[str, object]]
    requested_cycles: int
    completed_cycles: int
    wave_n: int
    wave_min: float
    wave_max: float
    wave_sha256: str
    fit_status: str = ""
    failure_reason: str = ""
    final_continuum_status: str = "not_run"
    final_line_status: str = "not_run"
    knot_strategy: str = "first_cycle_adaptive_then_fixed"
    legacy_kernel_representation: str = "native_grid_channel_median"
    schema_version: int = LSF_STATE_SCHEMA_VERSION


@dataclass(slots=True)
class LSFSurfaceIterativeResult(SkyDecompResult):
    """Baseline-compatible decomposition result with one additive LSF field."""

    lsf_state: LSFSurfaceState


def _wave_fingerprint(wave: np.ndarray) -> tuple[int, float, float, str]:
    value = np.asarray(wave, dtype=np.float64)
    canonical = np.ascontiguousarray(value.astype("<f8", copy=False))
    return (
        int(value.size),
        float(value[0]),
        float(value[-1]),
        hashlib.sha256(canonical.tobytes()).hexdigest(),
    )


def _channel_mask(
    wave: np.ndarray,
    lower: float | None,
    upper: float | None,
) -> np.ndarray:
    mask = np.ones(wave.size, dtype=bool)
    if lower is not None:
        mask &= wave >= lower
    if upper is not None:
        mask &= wave < upper
    return mask


def build_skyline_mask(
    wave: np.ndarray,
    line_stick_matrices: tuple[np.ndarray, ...],
    *,
    cumulative_fraction: float = 0.99,
    half_width_angstrom: float = 2.0,
) -> np.ndarray:
    """Mask the strongest reference lines in every physical line group."""
    wave = np.asarray(wave, dtype=float)
    if wave.ndim != 1 or wave.size == 0:
        raise ValueError("wave must be a nonempty one-dimensional array")
    if np.any(np.diff(wave) <= 0.0):
        raise ValueError("wave must be strictly increasing")
    if not 0.0 < cumulative_fraction <= 1.0:
        raise ValueError("cumulative_fraction must lie in (0, 1]")
    if half_width_angstrom <= 0.0:
        raise ValueError("half_width_angstrom must be positive")

    centers = np.zeros(wave.size, dtype=bool)
    for matrix in line_stick_matrices:
        rows = np.atleast_2d(np.asarray(matrix, dtype=float))
        if rows.shape[1] != wave.size:
            raise ValueError("line stick matrices must match wave")
        for row in rows:
            strength = np.clip(np.nan_to_num(np.abs(row)), 0.0, np.inf)
            positive = np.flatnonzero(strength > 0.0)
            if positive.size == 0:
                continue
            order = positive[np.argsort(strength[positive])[::-1]]
            cumulative = np.cumsum(strength[order])
            stop = int(
                np.searchsorted(
                    cumulative,
                    cumulative_fraction * cumulative[-1],
                    side="left",
                )
            )
            centers[order[: stop + 1]] = True

    center_wave = wave[centers]
    if center_wave.size == 0:
        return centers
    lower = np.searchsorted(
        wave,
        center_wave - half_width_angstrom,
        side="left",
    )
    upper = np.searchsorted(
        wave,
        center_wave + half_width_angstrom,
        side="right",
    )
    difference = np.zeros(wave.size + 1, dtype=int)
    np.add.at(difference, lower, 1)
    np.add.at(difference, upper, -1)
    return np.cumsum(difference[:-1]) > 0


def _robust_mad(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    median = float(np.median(finite))
    return float(1.4826 * np.median(np.abs(finite - median)))


def continuum_fit_weights(
    wave: np.ndarray,
    residual: np.ndarray,
    ivar: np.ndarray,
    skyline_mask: np.ndarray,
    *,
    line_weight: float = 5.0e-4,
    huber_transition_sigma: float = 3.0,
) -> tuple[np.ndarray, dict[str, float]]:
    """Return skyline-protected weights for a continuum-dominated solve."""
    wave = np.asarray(wave, dtype=float)
    residual = np.asarray(residual, dtype=float)
    ivar = np.asarray(ivar, dtype=float)
    skyline_mask = np.asarray(skyline_mask, dtype=bool)
    if not (
        wave.shape == residual.shape == ivar.shape == skyline_mask.shape
    ):
        raise ValueError("wave, residual, ivar, and skyline_mask must match")
    if not 0.0 <= line_weight <= 1.0:
        raise ValueError("line_weight must lie in [0, 1]")
    if huber_transition_sigma <= 0.0:
        raise ValueError("huber_transition_sigma must be positive")

    valid = (
        np.isfinite(residual)
        & np.isfinite(ivar)
        & (ivar > 0.0)
        & ~skyline_mask
    )
    global_sigma = _robust_mad(residual[valid])
    if not np.isfinite(global_sigma) or global_sigma <= 0.0:
        global_sigma = 1.0

    multiplier = np.ones_like(residual, dtype=float)
    multiplier[skyline_mask] = line_weight
    channel_noise: dict[str, float] = {}
    for channel, lower, upper in LSF_CHANNELS:
        channel_use = valid & _channel_mask(wave, lower, upper)
        sigma = _robust_mad(residual[channel_use])
        if not np.isfinite(sigma) or sigma <= 0.0:
            sigma = global_sigma
        channel_noise[channel] = float(sigma)

        absolute = np.abs(residual[channel_use])
        threshold = huber_transition_sigma * sigma
        huber = np.ones_like(absolute)
        outside = absolute > threshold
        huber[outside] = threshold / absolute[outside]
        multiplier[channel_use] *= huber

    weights = ivar * multiplier
    for _, lower, upper in LSF_CHANNELS:
        channel_use = (
            np.isfinite(weights)
            & (weights > 0.0)
            & _channel_mask(wave, lower, upper)
        )
        if np.any(channel_use):
            weights[channel_use] /= float(np.mean(weights[channel_use]))
    weights[~np.isfinite(weights) | (weights < 0.0)] = 0.0
    return weights, channel_noise


def _shifted_source(source: np.ndarray, kernel_size: int) -> np.ndarray:
    offsets = np.arange(-(kernel_size // 2), kernel_size // 2 + 1)
    shifted = np.zeros((source.size, kernel_size), dtype=float)
    for column, offset in enumerate(offsets):
        if offset == 0:
            shifted[:, column] = source
        elif offset > 0:
            shifted[offset:, column] = source[:-offset]
        else:
            shifted[:offset, column] = source[-offset:]
    return shifted


def _information_adaptive_knots(
    wave: np.ndarray,
    information: np.ndarray,
    n_basis: int,
    degree: int,
) -> np.ndarray:
    n_interior = n_basis - degree - 1
    if n_interior <= 0:
        return np.empty(0, dtype=float)

    weight = np.clip(np.nan_to_num(information), 0.0, np.inf)
    if np.sum(weight) == 0.0:
        return np.linspace(wave[0], wave[-1], n_interior + 2)[1:-1]

    cumulative = np.cumsum(weight) / np.sum(weight)
    quantiles = np.arange(1, n_interior + 1) / (n_interior + 1)
    knots = np.interp(quantiles, cumulative, wave)
    spacing = (wave[-1] - wave[0]) * min(0.08, 0.5 / (n_interior + 1))
    lower = wave[0] + spacing * np.arange(1, n_interior + 1)
    upper = wave[-1] - spacing * np.arange(n_interior, 0, -1)
    knots = np.clip(knots, lower, upper)
    for index in range(1, n_interior):
        knots[index] = max(knots[index], knots[index - 1] + spacing)
    for index in range(n_interior - 2, -1, -1):
        knots[index] = min(knots[index], knots[index + 1] - spacing)
    return knots


def build_bspline_basis(
    wave: np.ndarray,
    n_basis: int,
    degree: int,
    information: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a clamped partition-of-unity basis on one channel."""
    interior = _information_adaptive_knots(
        wave,
        information,
        n_basis,
        degree,
    )
    knots = np.concatenate(
        [
            np.repeat(wave[0], degree + 1),
            interior,
            np.repeat(wave[-1], degree + 1),
        ]
    )
    basis = BSpline.design_matrix(wave, knots, degree).toarray()
    basis = np.maximum(basis, 0.0)
    basis /= np.sum(basis, axis=1, keepdims=True)
    return basis, knots


def evaluate_bspline_basis(
    wave: np.ndarray,
    knots: np.ndarray,
    degree: int,
) -> np.ndarray:
    """Evaluate a persisted clamped B-spline basis without moving its knots."""
    wave = np.asarray(wave, dtype=float)
    knots = np.asarray(knots, dtype=float)
    if wave.ndim != 1 or wave.size == 0:
        raise ValueError("wave must be a nonempty one-dimensional array")
    if knots.ndim != 1:
        raise ValueError("knots must be one-dimensional")
    n_basis = knots.size - int(degree) - 1
    if n_basis < 1:
        raise ValueError("knots and degree do not define a B-spline basis")
    basis = BSpline.design_matrix(
        wave,
        knots,
        int(degree),
        extrapolate=False,
    ).toarray()
    basis = np.maximum(basis, 0.0)
    row_sum = np.sum(basis, axis=1, keepdims=True)
    if np.any(row_sum <= 0.0):
        raise ValueError("wave lies outside the persisted B-spline domain")
    return basis / row_sum


def _project_unimodal(kernel: np.ndarray) -> np.ndarray:
    value = np.clip(np.asarray(kernel, dtype=float), 0.0, np.inf)
    center = value.size // 2
    value[: center + 1] = np.maximum.accumulate(value[: center + 1])
    value[center:] = np.maximum.accumulate(value[center:][::-1])[::-1]
    value /= np.sum(value)
    return value


def _sum_constraints(kernel_size: int, n_basis: int) -> sp.csc_matrix:
    rows = np.tile(np.arange(n_basis), kernel_size)
    columns = np.arange(kernel_size * n_basis)
    return sp.coo_matrix(
        (np.ones(columns.size), (rows, columns)),
        shape=(n_basis, columns.size),
    ).tocsc()


def _monotonicity_constraints(
    kernel_size: int,
    n_basis: int,
) -> sp.csc_matrix:
    center = kernel_size // 2
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    row = 0
    for basis_index in range(n_basis):
        for tap in range(center):
            rows.extend([row, row])
            columns.extend(
                [
                    tap * n_basis + basis_index,
                    (tap + 1) * n_basis + basis_index,
                ]
            )
            values.extend([1.0, -1.0])
            row += 1
        for tap in range(center, kernel_size - 1):
            rows.extend([row, row])
            columns.extend(
                [
                    (tap + 1) * n_basis + basis_index,
                    tap * n_basis + basis_index,
                ]
            )
            values.extend([1.0, -1.0])
            row += 1
    return sp.coo_matrix(
        (values, (rows, columns)),
        shape=((kernel_size - 1) * n_basis, kernel_size * n_basis),
    ).tocsc()


def kernel_moments(kernel_surface: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return kernel centroid and standard deviation in pixel units."""
    offsets = np.arange(
        -(kernel_surface.shape[1] // 2),
        kernel_surface.shape[1] // 2 + 1,
        dtype=float,
    )
    center = kernel_surface @ offsets
    variance = np.sum(
        kernel_surface * (offsets[None, :] - center[:, None]) ** 2,
        axis=1,
    )
    return center, np.sqrt(np.maximum(variance, 0.0))


def _fallback_solution(
    basis: np.ndarray,
    fallback_kernel: np.ndarray,
    reason: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | str]]:
    coefficient = np.repeat(
        fallback_kernel[:, None],
        basis.shape[1],
        axis=1,
    )
    surface = basis @ coefficient.T
    center, sigma = kernel_moments(surface)
    metrics: dict[str, float | str] = {
        "status": "fallback",
        "reason": reason,
        "n_basis": float(basis.shape[1]),
        "chi2_red": np.nan,
        "rms_resid": np.nan,
        "center_pix_min": float(np.min(center)),
        "center_pix_median": float(np.median(center)),
        "center_pix_max": float(np.max(center)),
        "center_pix": float(np.median(center)),
        "sigma_pix_min": float(np.min(sigma)),
        "sigma_pix_median": float(np.median(sigma)),
        "sigma_pix_max": float(np.max(sigma)),
        "sigma_pix": float(np.median(sigma)),
        "runtime_sec": 0.0,
    }
    return surface, coefficient, metrics


def fit_bspline_channel(
    wave: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    ivar: np.ndarray,
    fallback_kernel: np.ndarray,
    *,
    n_basis: int = 6,
    degree: int = 3,
    roughness_fraction: float = 0.01,
    fallback_prior_fraction: float = 1.0e-4,
    information_prior_max_boost: float = 50.0,
    background_degree: int = 3,
    free_amplitude: bool = False,
    knot_vector: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float | str]]:
    """Fit one normalized, central-peak B-spline kernel surface."""
    wave = np.asarray(wave, dtype=float)
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    ivar = np.asarray(ivar, dtype=float)
    fallback = _project_unimodal(fallback_kernel)

    shifted = _shifted_source(source, fallback.size)
    information = np.clip(np.nan_to_num(ivar), 0.0, np.inf) * np.sum(
        shifted**2,
        axis=1,
    )
    if knot_vector is None:
        basis, knots = build_bspline_basis(wave, n_basis, degree, information)
    else:
        knots = np.asarray(knot_vector, dtype=float).copy()
        persisted_n_basis = knots.size - degree - 1
        if persisted_n_basis != n_basis:
            raise ValueError(
                "persisted knot vector does not match n_basis and degree"
            )
        basis = evaluate_bspline_basis(wave, knots, degree)
    kernel_design = (shifted[:, :, None] * basis[:, None, :]).reshape(
        wave.size,
        fallback.size * n_basis,
    )

    midpoint = 0.5 * (wave[0] + wave[-1])
    half_width = max(0.5 * (wave[-1] - wave[0]), 1.0)
    scaled_wave = (wave - midpoint) / half_width
    target_scale = max(float(np.sqrt(np.nanmean(target**2))), 1.0e-30)
    background = (
        np.polynomial.legendre.legvander(scaled_wave, background_degree)
        * target_scale
    )
    design = np.column_stack([kernel_design, background])
    valid = (
        np.isfinite(target)
        & np.isfinite(ivar)
        & (ivar > 0.0)
        & np.all(np.isfinite(design), axis=1)
    )

    n_kernel = fallback.size * n_basis
    n_nuisance = background.shape[1]
    n_parameter = n_kernel + n_nuisance
    if (
        np.count_nonzero(valid) <= n_parameter
        or not np.any(np.abs(kernel_design[valid]) > 0.0)
    ):
        surface, coefficient, metrics = _fallback_solution(
            basis,
            fallback,
            "insufficient_active_pixels",
        )
        return surface, coefficient, knots, metrics

    root_weight = np.sqrt(ivar[valid])
    weighted_design = design[valid] * root_weight[:, None]
    weighted_target = target[valid] * root_weight
    data_scale = max(float(np.sqrt(np.mean(weighted_target**2))), 1.0e-30)
    weighted_design /= data_scale
    weighted_target /= data_scale

    hessian = weighted_design.T @ weighted_design
    hessian_scale = max(float(np.trace(hessian)) / n_parameter, 1.0e-12)
    if n_basis >= 3 and roughness_fraction > 0.0:
        d2 = np.diff(np.eye(n_basis), n=2, axis=0)
        curvature = sp.kron(sp.eye(fallback.size), d2, format="csc")
        hessian[:n_kernel, :n_kernel] += (
            roughness_fraction
            * hessian_scale
            * (curvature.T @ curvature).toarray()
        )

    basis_information = np.divide(
        basis.T @ information,
        np.sum(basis, axis=0),
        out=np.zeros(n_basis),
        where=np.sum(basis, axis=0) > 0.0,
    )
    relative_information = basis_information / max(
        float(np.max(basis_information)),
        1.0e-30,
    )
    prior_boost = 1.0 + (information_prior_max_boost - 1.0) * (
        1.0 - np.clip(relative_information, 0.0, 1.0)
    ) ** 2
    prior_weight = np.tile(prior_boost, fallback.size)
    prior_scale = fallback_prior_fraction * hessian_scale
    hessian[:n_kernel, :n_kernel] += prior_scale * np.diag(prior_weight)
    hessian += 1.0e-10 * np.eye(n_parameter)

    fallback_vector = np.repeat(fallback[:, None], n_basis, axis=1).reshape(-1)
    linear = -(weighted_design.T @ weighted_target)
    linear[:n_kernel] -= prior_scale * prior_weight * fallback_vector

    monotonicity_kernel = _monotonicity_constraints(fallback.size, n_basis)
    nonnegative = sp.hstack(
        [-sp.eye(n_kernel), sp.csc_matrix((n_kernel, n_nuisance))],
        format="csc",
    )
    monotonicity = sp.hstack(
        [
            monotonicity_kernel,
            sp.csc_matrix((monotonicity_kernel.shape[0], n_nuisance)),
        ],
        format="csc",
    )
    constraints = []
    rhs = []
    cones = []
    if not free_amplitude:
        equality_kernel = _sum_constraints(fallback.size, n_basis)
        equality = sp.hstack(
            [equality_kernel, sp.csc_matrix((n_basis, n_nuisance))],
            format="csc",
        )
        constraints.append(equality)
        rhs.append(np.ones(n_basis))
        cones.append(clarabel.ZeroConeT(n_basis))
    constraints.extend([nonnegative, monotonicity])
    rhs.extend([np.zeros(n_kernel), np.zeros(monotonicity.shape[0])])
    cones.extend(
        [
            clarabel.NonnegativeConeT(n_kernel),
            clarabel.NonnegativeConeT(monotonicity.shape[0]),
        ]
    )

    settings = clarabel.DefaultSettings()
    settings.verbose = False
    started = time.perf_counter()
    solver = clarabel.DefaultSolver(
        sp.triu(sp.csc_matrix((hessian + hessian.T) * 0.5)).tocsc(),
        np.asarray(linear, dtype=np.float64),
        sp.vstack(constraints, format="csc"),
        np.concatenate(rhs),
        cones,
        settings,
    )
    solution = solver.solve()
    runtime = time.perf_counter() - started
    status = str(solution.status)
    vector = np.asarray(solution.x, dtype=float)
    if status not in {"Solved", "AlmostSolved"} or not np.all(np.isfinite(vector)):
        surface, coefficient, metrics = _fallback_solution(
            basis,
            fallback,
            "solver_failed",
        )
        metrics["runtime_sec"] = runtime
        return surface, coefficient, knots, metrics

    raw_coefficient = vector[:n_kernel].reshape(fallback.size, n_basis)
    raw_coefficient[np.abs(raw_coefficient) < 1.0e-12] = 0.0
    coefficient = np.maximum(raw_coefficient, 0.0)
    amplitude_scale = np.sum(coefficient, axis=0)
    if np.any(amplitude_scale <= 0.0):
        surface, coefficient, metrics = _fallback_solution(
            basis,
            fallback,
            "empty_kernel",
        )
        metrics["runtime_sec"] = runtime
        return surface, coefficient, knots, metrics
    coefficient /= amplitude_scale[None, :]
    surface = basis @ coefficient.T
    surface /= np.sum(surface, axis=1, keepdims=True)

    model = design @ vector
    residual = target[valid] - model[valid]
    chi2 = float(np.sum(residual**2 * ivar[valid]))
    center, sigma = kernel_moments(surface)
    metrics = {
        "status": status,
        "reason": "",
        "n_pixels": float(wave.size),
        "n_valid_pixels": float(np.count_nonzero(valid)),
        "n_basis": float(n_basis),
        "degree": float(degree),
        "chi2_red": chi2 / max(np.count_nonzero(valid) - n_parameter, 1),
        "rms_resid": float(np.sqrt(np.mean(residual**2))),
        "center_pix_min": float(np.min(center)),
        "center_pix_median": float(np.median(center)),
        "center_pix_max": float(np.max(center)),
        "center_pix": float(np.median(center)),
        "sigma_pix_min": float(np.min(sigma)),
        "sigma_pix_median": float(np.median(sigma)),
        "sigma_pix_max": float(np.max(sigma)),
        "sigma_pix": float(np.median(sigma)),
        "fitted_amplitude_scale": float(np.median(amplitude_scale)),
        "runtime_sec": runtime,
    }
    return surface, coefficient, knots, metrics


def fit_lsf_surface(
    wave: np.ndarray,
    flux: np.ndarray,
    ivar: np.ndarray,
    source: np.ndarray,
    fixed_background: np.ndarray,
    fallback_kernels: dict[str, np.ndarray],
    config: LSFSurfaceIterativeConfig,
    *,
    previous_state: LSFSurfaceState | None = None,
) -> LSFSurfaceState:
    """Fit the hybrid LSF, reusing the first fit's knots on later cycles."""
    wave = np.asarray(wave, dtype=float)
    target = np.asarray(flux, dtype=float) - np.asarray(
        fixed_background,
        dtype=float,
    )
    surface = np.zeros((wave.size, LSF_KERNEL_SIZE), dtype=float)
    previous_surface = (
        None
        if previous_state is None
        else evaluate_lsf_surface(previous_state, wave)
    )
    coefficients: dict[str, np.ndarray] = {}
    knots: dict[str, np.ndarray] = {}
    degrees: dict[str, int] = {}
    bounds: dict[str, tuple[float | None, float | None]] = {}
    metrics: dict[str, dict[str, object]] = {}

    for channel, lower, upper in LSF_CHANNELS:
        mask = _channel_mask(wave, lower, upper)
        if not np.any(mask):
            continue
        fit_mask = mask.copy()
        n_basis = config.n_basis
        degree = config.degree
        free_amplitude = False
        if channel == "B":
            fit_mask &= wave >= config.blue_fit_lower
            n_basis = 1
            degree = 0
            free_amplitude = True

        previous_knots = None
        if previous_state is not None:
            previous_knots = previous_state.knot_vectors.get(channel)

        local_surface, coefficient, knot_vector, channel_metrics = fit_bspline_channel(
            wave[fit_mask],
            source[fit_mask],
            target[fit_mask],
            ivar[fit_mask],
            fallback_kernels[channel],
            n_basis=n_basis,
            degree=degree,
            roughness_fraction=(
                0.0 if channel == "B" else config.roughness_fraction
            ),
            fallback_prior_fraction=config.fallback_prior_fraction,
            information_prior_max_boost=(
                1.0 if channel == "B" else config.information_prior_max_boost
            ),
            background_degree=config.background_degree,
            free_amplitude=free_amplitude,
            knot_vector=previous_knots,
        )

        restored_previous = (
            previous_state is not None
            and channel_metrics["status"] == "fallback"
        )
        if restored_previous:
            local_surface = previous_surface[mask].copy()
            coefficient = previous_state.coefficients[channel].copy()
            knot_vector = previous_state.knot_vectors[channel].copy()
            channel_metrics["reason"] = (
                f"previous_surface:{channel_metrics['reason']}"
            )
        if channel == "B":
            surface[mask] = local_surface if restored_previous else local_surface[0]
            channel_metrics["model"] = "constant_blue_kernel"
        else:
            surface[mask] = local_surface
            channel_metrics["model"] = "cubic_bspline_kernel"
        channel_metrics["fit_lower"] = float(np.min(wave[fit_mask]))
        channel_metrics["fit_upper"] = float(np.max(wave[fit_mask]))
        channel_metrics["knots_fixed"] = previous_knots is not None
        coefficients[channel] = coefficient
        knots[channel] = knot_vector
        degrees[channel] = degree
        bounds[channel] = (lower, upper)
        metrics[channel] = channel_metrics

    row_sum = np.sum(surface, axis=1, keepdims=True)
    if np.any(~np.isfinite(row_sum)) or np.any(row_sum <= 0.0):
        raise RuntimeError("LSF fit did not define a kernel at every wavelength")
    surface /= row_sum
    wave_n, wave_min, wave_max, wave_sha256 = _wave_fingerprint(wave)
    return LSFSurfaceState(
        coefficients=coefficients,
        knot_vectors=knots,
        degrees=degrees,
        channel_bounds=bounds,
        tap_offsets=LSF_TAP_OFFSETS.copy(),
        config=asdict(config),
        metrics=metrics,
        requested_cycles=config.n_refinement_cycles,
        completed_cycles=0,
        wave_n=wave_n,
        wave_min=wave_min,
        wave_max=wave_max,
        wave_sha256=wave_sha256,
    )


def apply_lsf_surface(
    wave: np.ndarray,
    values: np.ndarray,
    kernel_surface: np.ndarray,
    *,
    channel_bounds: dict[str, tuple[float | None, float | None]] | None = None,
    tap_offsets: np.ndarray | None = None,
) -> np.ndarray:
    """Apply a target-wavelength-dependent kernel without crossing channels."""
    wave = np.asarray(wave, dtype=float)
    source = np.asarray(values, dtype=float)
    kernel_surface = np.asarray(kernel_surface, dtype=float)
    if wave.ndim != 1 or wave.size == 0 or np.any(np.diff(wave) <= 0.0):
        raise ValueError("wave must be nonempty and strictly increasing")
    one_dimensional = source.ndim == 1
    if one_dimensional:
        source = source[None, :]
    if source.ndim != 2 or source.shape[1] != wave.size:
        raise ValueError("values must end with the wavelength dimension")
    if (
        kernel_surface.ndim != 2
        or kernel_surface.shape[0] != wave.size
        or kernel_surface.shape[1] % 2 != 1
    ):
        raise ValueError("kernel_surface must have shape (wave, odd kernel size)")
    if np.any(~np.isfinite(kernel_surface)) or np.any(kernel_surface < -1.0e-12):
        raise ValueError("kernel_surface must be finite and nonnegative")
    if not np.allclose(np.sum(kernel_surface, axis=1), 1.0, atol=1.0e-10):
        raise ValueError("each LSF kernel must sum to one")
    if tap_offsets is None:
        offsets = np.arange(
            -(kernel_surface.shape[1] // 2),
            kernel_surface.shape[1] // 2 + 1,
        )
    else:
        offsets = np.asarray(tap_offsets)
        if (
            offsets.ndim != 1
            or offsets.size != kernel_surface.shape[1]
            or not np.all(np.isfinite(offsets))
            or not np.all(offsets == offsets.astype(int))
        ):
            raise ValueError("tap_offsets must be finite integer pixel offsets")
        offsets = offsets.astype(int)
    bounds = (
        {channel: (lower, upper) for channel, lower, upper in LSF_CHANNELS}
        if channel_bounds is None
        else channel_bounds
    )
    output = np.zeros_like(source)
    for channel in ("B", "R", "Z"):
        lower, upper = bounds[channel]
        indices = np.flatnonzero(_channel_mask(wave, lower, upper))
        if indices.size == 0:
            continue
        for tap, offset in enumerate(offsets):
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
    return output[0] if one_dimensional else output


def evaluate_lsf_surface(
    state: LSFSurfaceState,
    wave: np.ndarray,
) -> np.ndarray:
    """Reconstruct a fitted LSF surface from persisted knots and coefficients."""
    if state.schema_version != LSF_STATE_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported LSF state schema version: {state.schema_version}"
        )
    evaluation_wave = np.asarray(wave, dtype=float)
    if (
        evaluation_wave.ndim != 1
        or evaluation_wave.size == 0
        or np.any(np.diff(evaluation_wave) <= 0.0)
    ):
        raise ValueError("wave must be nonempty and strictly increasing")

    surface = np.zeros((evaluation_wave.size, state.tap_offsets.size), dtype=float)
    covered = np.zeros(evaluation_wave.size, dtype=bool)
    for channel in ("B", "R", "Z"):
        if channel not in state.channel_bounds:
            continue
        lower, upper = state.channel_bounds[channel]
        mask = _channel_mask(evaluation_wave, lower, upper)
        if not np.any(mask):
            continue
        if channel not in state.coefficients:
            raise ValueError(f"LSF state does not contain channel {channel}")
        coefficient = np.asarray(state.coefficients[channel], dtype=float)
        if (
            coefficient.ndim != 2
            or coefficient.shape[0] != state.tap_offsets.size
            or np.any(~np.isfinite(coefficient))
            or np.any(coefficient < -1.0e-12)
        ):
            raise ValueError(f"Invalid LSF coefficients for channel {channel}")
        degree = int(state.degrees[channel])
        if degree == 0 and coefficient.shape[1] == 1:
            surface[mask] = coefficient[:, 0]
            covered[mask] = True
            continue
        basis = evaluate_bspline_basis(
            evaluation_wave[mask],
            state.knot_vectors[channel],
            degree,
        )
        if basis.shape[1] != coefficient.shape[1]:
            raise ValueError(
                f"LSF coefficient and knot dimensions disagree for channel {channel}"
            )
        surface[mask] = basis @ coefficient.T
        covered[mask] = True

    if not np.all(covered):
        raise ValueError("persisted LSF state does not cover the requested wave")
    if np.any(surface < -1.0e-10) or np.any(~np.isfinite(surface)):
        raise ValueError("persisted LSF state produces invalid kernel weights")
    surface = np.maximum(surface, 0.0)
    row_sum = np.sum(surface, axis=1, keepdims=True)
    if np.any(row_sum <= 0.0) or np.any(~np.isfinite(row_sum)):
        raise ValueError("persisted LSF state does not cover the requested wave")
    if not np.allclose(row_sum, 1.0, atol=1.0e-8, rtol=0.0):
        raise ValueError("persisted LSF kernels are not normalized")
    return surface / row_sum


def evaluate_lsf_kernel(
    state: LSFSurfaceState,
    wavelength: float,
) -> np.ndarray:
    """Evaluate the persisted 11-tap LSF at one wavelength."""
    if not np.isfinite(wavelength):
        raise ValueError("wavelength must be finite")
    return evaluate_lsf_surface(state, np.array([float(wavelength)]))[0]


def apply_lsf_channelwise(
    state: LSFSurfaceState,
    wave: np.ndarray,
    values: np.ndarray,
    *,
    allow_different_grid: bool = False,
) -> np.ndarray:
    """Reconstruct and apply a persisted LSF state to one or more spectra."""
    evaluation_wave = np.asarray(wave, dtype=float)
    if not allow_different_grid:
        wave_n, wave_min, wave_max, wave_sha256 = _wave_fingerprint(evaluation_wave)
        if (
            wave_n != state.wave_n
            or wave_min != state.wave_min
            or wave_max != state.wave_max
            or wave_sha256 != state.wave_sha256
        ):
            raise ValueError(
                "wave does not match the fitted pixel grid; pass "
                "allow_different_grid=True only when pixel-offset semantics are intended"
            )
    source = np.asarray(values, dtype=float)
    if source.ndim < 1 or source.shape[-1] != evaluation_wave.size:
        raise ValueError("values must end with the wavelength dimension")
    original_shape = source.shape
    flattened = source.reshape(-1, evaluation_wave.size)
    applied = apply_lsf_surface(
        evaluation_wave,
        flattened,
        evaluate_lsf_surface(state, evaluation_wave),
        channel_bounds=state.channel_bounds,
        tap_offsets=state.tap_offsets,
    )
    return applied.reshape(original_shape)


def _text_value(value: object) -> str:
    return value.decode().strip() if isinstance(value, bytes) else str(value).strip()


def load_lsf_surface_state(
    filename: str | Path,
    spectrum_index: int = 0,
) -> LSFSurfaceState:
    """Load one compact LSF state from an extended output FITS."""
    from astropy.io import fits

    with fits.open(filename, memmap=False) as hdul:
        for name in ("LSF_COEF", "LSF_KNOTS", "LSF_META"):
            if name not in hdul:
                raise KeyError(f"Missing {name} extension in {filename}")
        coefficient_cube = np.asarray(hdul["LSF_COEF"].data, dtype=float)
        knot_cube = np.asarray(hdul["LSF_KNOTS"].data, dtype=float)
        meta = hdul["LSF_META"].data

    if coefficient_cube.ndim != 4 or knot_cube.ndim != 3:
        raise ValueError("Stored LSF coefficient or knot arrays have invalid shapes")
    if spectrum_index < 0 or spectrum_index >= coefficient_cube.shape[0]:
        raise IndexError("spectrum_index is outside the stored LSF array")
    rows = meta[meta["spectrum_index"] == spectrum_index]
    if len(rows) != len(LSF_CHANNELS):
        raise ValueError("LSF_META does not contain exactly one row per channel")
    if not np.all(np.asarray(rows["available"], dtype=bool)):
        raise ValueError("No fitted LSF state is available for this spectrum")

    coefficients: dict[str, np.ndarray] = {}
    knots: dict[str, np.ndarray] = {}
    degrees: dict[str, int] = {}
    bounds: dict[str, tuple[float | None, float | None]] = {}
    metrics: dict[str, dict[str, object]] = {}
    config: dict[str, float | int] | None = None
    first_row = rows[0]
    tap_offsets = np.asarray(first_row["tap_offsets"], dtype=int)
    if tap_offsets.ndim != 1 or tap_offsets.size != coefficient_cube.shape[2]:
        raise ValueError("Stored tap offsets do not match LSF_COEF")
    schema_version = int(first_row["schema_version"])
    if schema_version != LSF_STATE_SCHEMA_VERSION:
        raise ValueError(f"Unsupported LSF state schema version: {schema_version}")
    channel_to_index = {name: index for index, (name, _, _) in enumerate(LSF_CHANNELS)}
    for row in rows:
        channel = _text_value(row["channel"])
        if channel not in channel_to_index or channel in coefficients:
            raise ValueError(f"Invalid or duplicate LSF channel row: {channel}")
        if not np.array_equal(np.asarray(row["tap_offsets"], dtype=int), tap_offsets):
            raise ValueError("Stored LSF tap offsets differ between channels")
        channel_index = channel_to_index[channel]
        n_basis = int(row["n_basis"])
        n_knots = int(row["n_knots"])
        if n_basis < 1 or n_knots < 2:
            raise ValueError(f"Invalid stored LSF dimensions for channel {channel}")
        coefficients[channel] = coefficient_cube[
            spectrum_index, channel_index, :, :n_basis
        ].copy()
        knots[channel] = knot_cube[
            spectrum_index, channel_index, :n_knots
        ].copy()
        degrees[channel] = int(row["degree"])
        lower_value = float(row["lower"])
        upper_value = float(row["upper"])
        bounds[channel] = (
            None if np.isnan(lower_value) else lower_value,
            None if np.isnan(upper_value) else upper_value,
        )
        metrics[channel] = {
            "status": _text_value(row["status"]),
            "reason": _text_value(row["reason"]),
            "n_pixels": int(row["n_pixels"]),
            "n_valid_pixels": int(row["n_valid_pixels"]),
            "chi2_red": float(row["chi2_red"]),
            "rms_resid": float(row["rms_resid"]),
            "center_pix_min": float(row["center_pix_min"]),
            "center_pix_median": float(row["center_pix_median"]),
            "center_pix_max": float(row["center_pix_max"]),
            "sigma_pix_min": float(row["sigma_pix_min"]),
            "sigma_pix_median": float(row["sigma_pix_median"]),
            "sigma_pix_max": float(row["sigma_pix_max"]),
            "runtime_sec": float(row["runtime_sec"]),
            "fit_lower": float(row["fit_lower"]),
            "fit_upper": float(row["fit_upper"]),
            "knots_fixed": bool(row["knots_fixed"]),
            "model": _text_value(row["model"]),
        }
        if config is None:
            config = {
                "line_weight": float(row["line_weight"]),
                "skyline_cumulative_fraction": float(
                    row["skyline_cumulative_fraction"]
                ),
                "skyline_half_width_angstrom": float(
                    row["skyline_half_width_angstrom"]
                ),
                "huber_transition_sigma": float(row["huber_transition_sigma"]),
                "n_refinement_cycles": int(row["requested_cycles"]),
                "n_basis": int(row["config_n_basis"]),
                "degree": int(row["config_degree"]),
                "roughness_fraction": float(row["roughness_fraction"]),
                "fallback_prior_fraction": float(row["fallback_prior_fraction"]),
                "information_prior_max_boost": float(row["information_prior_max_boost"]),
                "background_degree": int(row["background_degree"]),
                "blue_fit_lower": float(row["blue_fit_lower"]),
            }

    return LSFSurfaceState(
        coefficients=coefficients,
        knot_vectors=knots,
        degrees=degrees,
        channel_bounds=bounds,
        tap_offsets=tap_offsets,
        config={} if config is None else config,
        metrics=metrics,
        requested_cycles=int(first_row["requested_cycles"]),
        completed_cycles=int(first_row["completed_cycles"]),
        wave_n=int(first_row["wave_n"]),
        wave_min=float(first_row["wave_min"]),
        wave_max=float(first_row["wave_max"]),
        wave_sha256=_text_value(first_row["wave_sha256"]),
        fit_status=_text_value(first_row["fit_status"]),
        failure_reason=_text_value(first_row["failure_reason"]),
        final_continuum_status=_text_value(first_row["final_continuum_status"]),
        final_line_status=_text_value(first_row["final_line_status"]),
        knot_strategy=_text_value(first_row["knot_strategy"]),
        legacy_kernel_representation=_text_value(
            first_row["legacy_kernel_representation"]
        ),
        schema_version=schema_version,
    )


def _nominal_lsf_state(
    wave: np.ndarray,
    fallback_kernels: dict[str, np.ndarray],
    config: LSFSurfaceIterativeConfig,
    reason: str,
) -> LSFSurfaceState:
    coefficients: dict[str, np.ndarray] = {}
    knots: dict[str, np.ndarray] = {}
    degrees: dict[str, int] = {}
    bounds: dict[str, tuple[float | None, float | None]] = {}
    metrics: dict[str, dict[str, object]] = {}
    for channel, lower, upper in LSF_CHANNELS:
        mask = _channel_mask(wave, lower, upper)
        if not np.any(mask):
            continue
        kernel = _project_unimodal(fallback_kernels[channel])
        coefficients[channel] = kernel[:, None]
        knots[channel] = np.array([wave[mask][0], wave[mask][-1]], dtype=float)
        degrees[channel] = 0
        bounds[channel] = (lower, upper)
        center, sigma = kernel_moments(kernel[None, :])
        metrics[channel] = {
            "status": "fallback",
            "reason": reason,
            "n_pixels": int(np.sum(mask)),
            "n_valid_pixels": 0,
            "chi2_red": np.nan,
            "rms_resid": np.nan,
            "center_pix_min": float(center[0]),
            "center_pix_median": float(center[0]),
            "center_pix_max": float(center[0]),
            "sigma_pix_min": float(sigma[0]),
            "sigma_pix_median": float(sigma[0]),
            "sigma_pix_max": float(sigma[0]),
            "runtime_sec": 0.0,
            "fit_lower": float(wave[mask][0]),
            "fit_upper": float(wave[mask][-1]),
            "knots_fixed": True,
            "model": "constant_nominal_fallback",
        }
    wave_n, wave_min, wave_max, wave_sha256 = _wave_fingerprint(wave)
    return LSFSurfaceState(
        coefficients=coefficients,
        knot_vectors=knots,
        degrees=degrees,
        channel_bounds=bounds,
        tap_offsets=LSF_TAP_OFFSETS.copy(),
        config=asdict(config),
        metrics=metrics,
        requested_cycles=config.n_refinement_cycles,
        completed_cycles=0,
        wave_n=wave_n,
        wave_min=wave_min,
        wave_max=wave_max,
        wave_sha256=wave_sha256,
        fit_status=f"failed:{reason}",
        failure_reason=reason,
        knot_strategy="nominal_constant_fallback",
    )


class SkyDecompLSFSurfaceIterative(SkyDecomp):
    """Iteratively refine continuum, smooth LSF, and emission-line amplitudes."""

    _CONTINUUM_KEYS = ("moon", "diffuse")
    _LINE_KEYS = ("oh", "atom", "orc", "o2")
    _SOLVED = {"Solved", "AlmostSolved"}

    def __init__(
        self,
        *args,
        config: LSFSurfaceIterativeConfig | None = None,
        **kwargs,
    ) -> None:
        if float(kwargs.get("moon_interline_boost", 0.0)) != 0.0:
            raise ValueError(
                "SkyDecompLSFSurfaceIterative requires moon_interline_boost=0"
            )
        kwargs["moon_interline_boost"] = 0.0
        self.config = config or LSFSurfaceIterativeConfig()
        self.lsf_surface_state: LSFSurfaceState | None = None
        self._lsf_surface: np.ndarray | None = None
        self.stage_metrics: tuple[dict[str, object], ...] = ()
        self.skyline_mask = np.array([], dtype=bool)
        self.continuum_weights = np.array([], dtype=float)
        self.channel_noise: dict[str, float] = {}
        self.final_continuum = np.array([], dtype=float)
        self.final_line_model = np.array([], dtype=float)
        super().__init__(*args, **kwargs)

    @property
    def lsf_surface(self) -> np.ndarray | None:
        return self._lsf_surface

    @staticmethod
    def _stack_matrices(
        matrices: dict[str, np.ndarray],
        keys: tuple[str, ...],
    ) -> np.ndarray:
        return np.vstack([matrices[key] for key in keys])

    def _line_stick_matrices(self) -> tuple[np.ndarray, ...]:
        return (
            self.matrix_oh_stick,
            self.matrix_atom_stick,
            self.matrix_orc_stick,
            self.matrix_o2_stick,
        )

    def _set_lsf_state(self, state: LSFSurfaceState | None) -> None:
        self.lsf_surface_state = state
        if state is None:
            self._lsf_surface = None
            self.lsf_metrics = {}
            self.lsf_kernels = {}
            return
        self._lsf_surface = evaluate_lsf_surface(state, self.wave)
        self.lsf_metrics = {
            channel: {
                **metric,
                "legacy_kernel_representation": state.legacy_kernel_representation,
            }
            for channel, metric in state.metrics.items()
        }
        self.lsf_kernels = {
            channel: np.median(
                self._lsf_surface[_channel_mask(self.wave, lower, upper)],
                axis=0,
            )
            for channel, lower, upper in LSF_CHANNELS
            if channel in state.channel_bounds
        }

    def _fit_lsf_channels(
        self,
        flux: np.ndarray,
        ivar: np.ndarray,
        source: np.ndarray,
        fixed_background: np.ndarray,
    ) -> LSFSurfaceState:
        previous = self.lsf_surface_state
        previous_surface = self._lsf_surface
        fallback = {}
        for channel, lower, upper in LSF_CHANNELS:
            mask = _channel_mask(self.wave, lower, upper)
            if previous is None:
                fallback[channel] = self._default_channel_kernel(mask)
            else:
                fallback[channel] = np.median(
                    previous_surface[mask],
                    axis=0,
                )
        state = fit_lsf_surface(
            self.wave,
            flux,
            ivar,
            source,
            fixed_background,
            fallback,
            self.config,
            previous_state=previous,
        )
        self._set_lsf_state(state)
        return state

    def _convolve_matrix_channelwise(self, matrix: np.ndarray) -> np.ndarray:
        if self._lsf_surface is None:
            return super()._convolve_matrix_channelwise(matrix)
        return apply_lsf_surface(
            self.wave,
            matrix,
            self._lsf_surface,
        )

    def _combine_coefficients(
        self,
        matrices: dict[str, np.ndarray],
        continuum_coefficient: np.ndarray,
        line_coefficient: np.ndarray,
    ) -> np.ndarray:
        slices = self._component_slices(matrices)
        coefficient = np.zeros(slices["o2"].stop, dtype=float)
        offset = 0
        for key in self._CONTINUUM_KEYS:
            size = matrices[key].shape[0]
            coefficient[slices[key]] = continuum_coefficient[offset : offset + size]
            offset += size
        offset = 0
        for key in self._LINE_KEYS:
            size = matrices[key].shape[0]
            coefficient[slices[key]] = line_coefficient[offset : offset + size]
            offset += size
        return coefficient

    def _line_source(self, line_coefficient: np.ndarray) -> np.ndarray:
        return np.vstack(self._line_stick_matrices()).T @ line_coefficient

    @staticmethod
    def _relative_change(new: np.ndarray, old: np.ndarray) -> float:
        denominator = max(float(np.linalg.norm(old)), 1.0e-30)
        return float(np.linalg.norm(new - old) / denominator)

    @staticmethod
    def _metric(
        stage: str,
        status: str,
        flux: np.ndarray,
        model: np.ndarray,
        ivar: np.ndarray,
        skyline_mask: np.ndarray,
        **extra: object,
    ) -> dict[str, object]:
        residual = flux - model
        valid = np.isfinite(residual) & np.isfinite(ivar) & (ivar > 0.0)
        line_free = valid & ~skyline_mask
        skyline = valid & skyline_mask

        def rms(mask: np.ndarray) -> float:
            return float(np.sqrt(np.mean(residual[mask] ** 2))) if np.any(mask) else np.nan

        metric: dict[str, object] = {
            "stage": stage,
            "solver_status": status,
            "chi2_per_valid_pixel": float(
                np.sum(residual[valid] ** 2 * ivar[valid])
                / max(int(np.sum(valid)), 1)
            ),
            "rms": rms(valid),
            "line_free_rms": rms(line_free),
            "line_free_mad": _robust_mad(residual[line_free]),
            "skyline_rms": rms(skyline),
            "skyline_mad": _robust_mad(residual[skyline]),
        }
        metric.update(extra)
        return metric

    def fit(
        self,
        flux: np.ndarray,
        ivar: np.ndarray,
        *,
        verbose: bool = False,
    ) -> LSFSurfaceIterativeResult:
        """Run the approved continuum -> LSF -> lines refinement sequence."""
        flux = np.asarray(flux, dtype=float)
        ivar = np.asarray(ivar, dtype=float)
        if flux.shape != self.wave.shape or ivar.shape != self.wave.shape:
            raise ValueError("flux and ivar must match wave")

        trace_started = False
        if not tracemalloc.is_tracing():
            tracemalloc.start()
            trace_started = True
        tracemalloc.reset_peak()
        started = time.perf_counter()
        config = self.config
        self._set_lsf_state(None)
        self._prefit_o2(flux, ivar)
        matrices = self._matrix_bundle(
            self.matrix_oh,
            self.matrix_moon,
            self.matrix_diffuse,
            self.matrix_atom,
            self.matrix_orc,
            self.matrix_o2,
        )
        full_design = self._stack_matrices(
            matrices,
            ("oh", "moon", "diffuse", "atom", "orc", "o2"),
        )
        slices = self._component_slices(matrices)
        seed = self._fit_design(
            full_design,
            flux,
            ivar,
            moon_slice=slices["moon"],
            diffuse_slice=slices["diffuse"],
        )
        seed_coefficient = np.asarray(seed["coef"], dtype=float)
        continuum_coefficient = np.concatenate(
            [seed_coefficient[slices[key]] for key in self._CONTINUUM_KEYS]
        )
        line_coefficient = np.concatenate(
            [seed_coefficient[slices[key]] for key in self._LINE_KEYS]
        )
        continuum_design = self._stack_matrices(matrices, self._CONTINUUM_KEYS)
        line_design = self._stack_matrices(matrices, self._LINE_KEYS)
        continuum = continuum_design.T @ continuum_coefficient
        line_model = line_design.T @ line_coefficient

        skyline_mask = build_skyline_mask(
            self.wave,
            self._line_stick_matrices(),
            cumulative_fraction=config.skyline_cumulative_fraction,
            half_width_angstrom=config.skyline_half_width_angstrom,
        )
        stage_metrics: list[dict[str, object]] = [
            self._metric(
                "00_nominal_joint_seed",
                str(seed["status"]),
                flux,
                continuum + line_model,
                ivar,
                skyline_mask,
            )
        ]
        continuum_weights = ivar.copy()
        channel_noise = {channel: np.nan for channel, _, _ in LSF_CHANNELS}
        failure: str | None = None
        final_solver_status = str(seed["status"])
        seed_solved = final_solver_status in self._SOLVED
        final_continuum_status = "not_run"
        final_line_status = "not_run"

        if not seed_solved:
            failure = f"nominal_seed:{final_solver_status}"

        for cycle in range(1, config.n_refinement_cycles + 1):
            if failure is not None:
                break
            old_continuum = continuum.copy()
            old_line_model = line_model.copy()
            previous_state = self.lsf_surface_state
            previous_surface = (
                None if self._lsf_surface is None else self._lsf_surface.copy()
            )
            previous_matrices = matrices
            previous_continuum_coefficient = continuum_coefficient.copy()
            previous_line_coefficient = line_coefficient.copy()

            continuum_weights, channel_noise = continuum_fit_weights(
                self.wave,
                flux - continuum - line_model,
                ivar,
                skyline_mask,
                line_weight=config.line_weight,
                huber_transition_sigma=config.huber_transition_sigma,
            )
            continuum_design = self._stack_matrices(matrices, self._CONTINUUM_KEYS)
            n_moon = matrices["moon"].shape[0]
            continuum_fit = self._fit_design(
                continuum_design,
                flux - line_model,
                continuum_weights,
                moon_slice=slice(0, n_moon),
                diffuse_slice=slice(n_moon, continuum_design.shape[0]),
            )
            continuum_status = str(continuum_fit["status"])
            if continuum_status not in self._SOLVED:
                failure = f"cycle_{cycle}_continuum:{continuum_status}"
                stage_metrics.append(
                    self._metric(
                        f"{cycle:02d}_continuum_failed",
                        continuum_status,
                        flux,
                        continuum + line_model,
                        ivar,
                        skyline_mask,
                    )
                )
                break
            candidate_continuum_coefficient = np.asarray(
                continuum_fit["coef"], dtype=float
            )
            candidate_continuum = continuum_design.T @ candidate_continuum_coefficient
            stage_metrics.append(
                self._metric(
                    f"{cycle:02d}_continuum",
                    continuum_status,
                    flux,
                    candidate_continuum + line_model,
                    ivar,
                    skyline_mask,
                    continuum_relative_change=self._relative_change(
                        candidate_continuum, old_continuum
                    ),
                )
            )

            state = self._fit_lsf_channels(
                flux,
                ivar,
                self._line_source(line_coefficient),
                candidate_continuum,
            )
            matrices = self._assemble_refined_matrices()
            line_design = self._stack_matrices(matrices, self._LINE_KEYS)
            line_fit = self._fit_design(
                line_design,
                flux - candidate_continuum,
                ivar,
            )
            line_status = str(line_fit["status"])
            if line_status not in self._SOLVED:
                self._set_lsf_state(previous_state)
                matrices = previous_matrices
                continuum_coefficient = previous_continuum_coefficient
                line_coefficient = previous_line_coefficient
                continuum = old_continuum
                line_model = old_line_model
                failure = f"cycle_{cycle}_lines:{line_status}"
                stage_metrics.append(
                    self._metric(
                        f"{cycle:02d}_lines_failed",
                        line_status,
                        flux,
                        continuum + line_model,
                        ivar,
                        skyline_mask,
                    )
                )
                break

            continuum_coefficient = candidate_continuum_coefficient
            continuum = candidate_continuum
            line_coefficient = np.asarray(line_fit["coef"], dtype=float)
            line_model = line_design.T @ line_coefficient
            state.completed_cycles = cycle
            final_solver_status = line_status
            lsf_change = (
                np.nan
                if previous_state is None
                else self._relative_change(
                    self._lsf_surface,
                    previous_surface,
                )
            )
            stage_metrics.append(
                self._metric(
                    f"{cycle:02d}_lsf_lines",
                    line_status,
                    flux,
                    continuum + line_model,
                    ivar,
                    skyline_mask,
                    line_relative_change=self._relative_change(
                        line_model, old_line_model
                    ),
                    lsf_relative_change=lsf_change,
                    lsf_status={
                        channel: metric["status"]
                        for channel, metric in state.metrics.items()
                    },
                )
            )

        if seed_solved:
            saved_continuum = continuum.copy()
            saved_line_model = line_model.copy()
            saved_continuum_coefficient = continuum_coefficient.copy()
            continuum_weights, channel_noise = continuum_fit_weights(
                self.wave,
                flux - continuum - line_model,
                ivar,
                skyline_mask,
                line_weight=config.line_weight,
                huber_transition_sigma=config.huber_transition_sigma,
            )
            continuum_design = self._stack_matrices(matrices, self._CONTINUUM_KEYS)
            n_moon = matrices["moon"].shape[0]
            continuum_fit = self._fit_design(
                continuum_design,
                flux - line_model,
                continuum_weights,
                moon_slice=slice(0, n_moon),
                diffuse_slice=slice(n_moon, continuum_design.shape[0]),
            )
            continuum_status = str(continuum_fit["status"])
            final_continuum_status = continuum_status
            if continuum_status not in self._SOLVED:
                final_failure = f"final_continuum:{continuum_status}"
                failure = (
                    final_failure
                    if failure is None
                    else f"{failure};{final_failure}"
                )
            else:
                candidate_continuum_coefficient = np.asarray(
                    continuum_fit["coef"], dtype=float
                )
                candidate_continuum = (
                    continuum_design.T @ candidate_continuum_coefficient
                )
                stage_metrics.append(
                    self._metric(
                        "final_continuum",
                        continuum_status,
                        flux,
                        candidate_continuum + line_model,
                        ivar,
                        skyline_mask,
                        continuum_relative_change=self._relative_change(
                            candidate_continuum, continuum
                        ),
                    )
                )
                line_design = self._stack_matrices(matrices, self._LINE_KEYS)
                line_fit = self._fit_design(
                    line_design,
                    flux - candidate_continuum,
                    ivar,
                )
                line_status = str(line_fit["status"])
                final_line_status = line_status
                if line_status not in self._SOLVED:
                    continuum = saved_continuum
                    line_model = saved_line_model
                    continuum_coefficient = saved_continuum_coefficient
                    final_failure = f"final_lines:{line_status}"
                    failure = (
                        final_failure
                        if failure is None
                        else f"{failure};{final_failure}"
                    )
                else:
                    continuum_coefficient = candidate_continuum_coefficient
                    continuum = candidate_continuum
                    line_coefficient = np.asarray(line_fit["coef"], dtype=float)
                    line_model = line_design.T @ line_coefficient
                    final_solver_status = line_status
                    stage_metrics.append(
                        self._metric(
                            "final_lines",
                            line_status,
                            flux,
                            continuum + line_model,
                            ivar,
                            skyline_mask,
                            line_relative_change=self._relative_change(
                                line_model, saved_line_model
                            ),
                        )
                    )

        coefficient = self._combine_coefficients(
            matrices,
            continuum_coefficient,
            line_coefficient,
        )
        components = self._components_from_coef(coefficient, matrices)
        continuum = components["moon"] + components["diffuse"]
        line_model = (
            components["oh"]
            + components["atom"]
            + components["orc"]
            + components["o2"]
        )
        bestfit_lsf = continuum + line_model
        component_total = (
            components["moon"]
            + components["diffuse"]
            + components["oh"]
            + components["atom"]
            + components["orc"]
            + components["o2"]
        )
        if not np.allclose(bestfit_lsf, component_total, atol=1.0e-10, rtol=1.0e-10):
            raise RuntimeError("Final component spectra do not close to bestfit_lsf")
        residual = flux - bestfit_lsf
        valid = np.isfinite(residual) & np.isfinite(ivar) & (ivar > 0.0)
        chi2 = float(np.sum(residual[valid] ** 2 * ivar[valid]))
        reduced_chi2 = chi2 / max(int(np.sum(valid)) - coefficient.size, 1)
        rms_resid = (
            float(np.sqrt(np.mean(residual[valid] ** 2))) if np.any(valid) else np.nan
        )
        y_mean = float(np.average(flux[valid], weights=ivar[valid])) if np.any(valid) else np.nan
        total_sum = (
            float(np.sum((flux[valid] - y_mean) ** 2 * ivar[valid]))
            if np.any(valid)
            else np.nan
        )
        r2 = (
            1.0 - chi2 / total_sum
            if np.isfinite(total_sum) and total_sum > 0.0
            else np.nan
        )
        elapsed = time.perf_counter() - started
        peak_memory_mb = tracemalloc.get_traced_memory()[1] / 1024**2
        if trace_started:
            tracemalloc.stop()

        self.coef = coefficient
        self.bestfit = np.asarray(seed["bestfit"], dtype=float)
        self.bestfit_lsf = bestfit_lsf
        self.fit_status = (
            final_solver_status if failure is None else f"failed:{failure}"
        )
        self.r2 = r2
        self.rms_resid = rms_resid
        self.peak_memory_mb = peak_memory_mb
        self.stage_metrics = tuple(stage_metrics)
        self.skyline_mask = skyline_mask
        self.continuum_weights = continuum_weights
        self.channel_noise = channel_noise
        self.final_continuum = continuum
        self.final_line_model = line_model
        if self.lsf_surface_state is None:
            fallback = {
                channel: self._default_channel_kernel(
                    _channel_mask(self.wave, lower, upper)
                )
                for channel, lower, upper in LSF_CHANNELS
            }
            self._set_lsf_state(
                _nominal_lsf_state(
                    self.wave,
                    fallback,
                    config,
                    failure or "no_completed_lsf_cycle",
                )
            )
        completed_cycles = self.lsf_surface_state.completed_cycles
        self.lsf_surface_state.fit_status = self.fit_status
        self.lsf_surface_state.failure_reason = failure or ""
        self.lsf_surface_state.final_continuum_status = final_continuum_status
        self.lsf_surface_state.final_line_status = final_line_status
        self.fit_summary = (
            f"status={self.fit_status} | npar={coefficient.size} | "
            f"ngood={int(np.sum(valid))} | chi2_red={reduced_chi2:.4g} | "
            f"R2={r2:.5f} | line_weight={config.line_weight:.3g} | "
            f"skyline_fraction={float(np.mean(skyline_mask)):.3f} | "
            f"refinement_cycles={completed_cycles}/{config.n_refinement_cycles} | "
            f"fixed_lsf_knots=True | dt={elapsed:.2f}s"
        )

        if verbose:
            for metric in stage_metrics:
                print(
                    f"{metric['stage']}: chi2/pix={metric['chi2_per_valid_pixel']:.5g}, "
                    f"line-free RMS={metric['line_free_rms']:.5g}, "
                    f"skyline RMS={metric['skyline_rms']:.5g}"
                )

        return LSFSurfaceIterativeResult(
            coef=self.coef,
            bestfit=np.asarray(seed["bestfit"], dtype=float),
            resid=np.asarray(seed["resid"], dtype=float),
            resid_level=float(seed["resid_level"]),
            fit_status=self.fit_status,
            fit_summary=self.fit_summary,
            reduced_chi2=reduced_chi2,
            fit_elapsed_sec=elapsed,
            components=components,
            design_names=self.design_names,
            t_o2=self.t_o2,
            t_o2_err=self.t_o2_err,
            r2=r2,
            rms_resid=rms_resid,
            peak_memory_mb=peak_memory_mb,
            o2_fit_status=self.o2_fit_status,
            o2_fit_summary=self.o2_fit_summary,
            o2_fit_elapsed_sec=self.o2_fit_elapsed_sec,
            o2_valid_frac=self.o2_valid_frac,
            lsf_kernels=self.lsf_kernels,
            lsf_metrics=self.lsf_metrics,
            bestfit_lsf=bestfit_lsf,
            moon_knots=self.moon_knots_used.copy(),
            moon_boosted_pixels=self.moon_boosted_pixels_used.copy(),
            lsf_state=self.lsf_surface_state,
        )


__all__ = [
    "LSFSurfaceIterativeConfig",
    "LSFSurfaceIterativeResult",
    "LSFSurfaceState",
    "SkyDecompLSFSurfaceIterative",
    "apply_lsf_channelwise",
    "apply_lsf_surface",
    "build_bspline_basis",
    "build_skyline_mask",
    "continuum_fit_weights",
    "evaluate_bspline_basis",
    "evaluate_lsf_kernel",
    "evaluate_lsf_surface",
    "fit_bspline_channel",
    "fit_lsf_surface",
    "kernel_moments",
    "load_lsf_surface_state",
]
