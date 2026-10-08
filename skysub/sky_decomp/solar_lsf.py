"""Blue-channel LSF fitted jointly from sky emission lines and solar Fraunhofer lines.

The Moon and Zodi continua carry the solar absorption spectrum, so with the
continuum coefficients fixed the blue-channel model is linear in the LSF
coefficients ``theta[k, b]``:

    y_B = line_design @ theta + C(lambda) * (G @ theta) + diffuse + background,
    G[i, k, b] = (1 / width_i) Int_bin_i dl Int S_sun(l') B_b(l') M_k(l - l') dl',
    C(lambda) = sum_j c_j E_j(lambda).

``G`` integrates every original Meftah sample exactly (no regridding or
smoothing of the solar source) and does not depend on the fitted row.  ``E_j``
are the smooth solar-free envelopes of the Moon and Zodi templates.
"""

from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
from pathlib import Path
import time
from typing import Any

import numpy as np
import scipy.sparse as sp

from .fit import LSF_CHANNELS, _solar_templates, read_static_matrix, vac_to_air
from .lsf_spline2d import (
    _build_integrated_components,
    native_pixel_edges,
    evaluate_lsf_diagnostics,
)
from .lsf_surface_iterative import (
    LSFChannelSplineConfig,
    LSFSurfaceIterativeConfig,
    _channel_mask,
    _configured_knot_vector,
    _resolved_spline_config,
    evaluate_bspline_basis,
)
from .residual_pca import SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D

SOLAR_BLUE_MODEL = "solar_continuum_joint_mspline"
_BLUE = LSF_CHANNELS[0]


def _solar_tensor(
    wave: np.ndarray,
    knots: np.ndarray,
    degree: int,
    solar_path: str,
    chunk: int = 20000,
    offset_basis_count: int = 11, offset_half_width_angstrom: float = 3.0,
) -> np.ndarray:
    """Exact midpoint-source-cell quadrature of ``G`` over every Meftah sample."""
    raw = read_static_matrix(solar_path, ";")
    source_wave = vac_to_air(raw[:, 0] * 10.0)
    # Same unit as the Moon/Zodi templates: median over the whole file.
    flux = raw[:, 1] / np.nanmedian(raw[:, 1])
    weight = flux * np.diff(native_pixel_edges(source_wave))
    output = _channel_mask(wave, *_BLUE[1:])
    rows = np.flatnonzero(output)
    edges = native_pixel_edges(wave)
    # Compact support of the configured M-splines makes all other samples exactly zero.
    active = np.flatnonzero(
        (source_wave >= edges[rows[0]] - offset_half_width_angstrom) & (source_wave <= edges[rows[-1] + 1] + offset_half_width_angstrom)
        & np.isfinite(weight)
    )
    n_basis = knots.size - degree - 1
    tensor = np.zeros((rows.size, offset_basis_count, n_basis))
    for start in range(0, active.size, chunk):
        selected = active[start : start + chunk]
        basis = evaluate_bspline_basis(
            np.clip(source_wave[selected], knots[degree], knots[-degree - 1]), knots, degree
        )
        columns = weight[selected, None] * basis
        operator = _build_integrated_components(wave, source_wave[selected], output, offset_basis_count, offset_half_width_angstrom).tocoo()
        keep = np.isin(operator.row, rows)
        row, col, data = operator.row[keep], operator.col[keep], operator.data[keep]
        row = np.searchsorted(rows, row)
        for offset in range(offset_basis_count):
            pick = col % offset_basis_count == offset
            part = sp.csr_matrix(
                (data[pick], (row[pick], col[pick] // offset_basis_count)),
                shape=(rows.size, selected.size),
            )
            tensor[:, offset] += part @ columns
    return tensor


@lru_cache(maxsize=4)
def _solar_tensor_cached(wave_bytes, knots_bytes, degree, solar_path, size, mtime, offset_basis_count, offset_half_width_angstrom) -> np.ndarray:
    tensor = _solar_tensor(
        np.frombuffer(wave_bytes, dtype=float),
        np.frombuffer(knots_bytes, dtype=float),
        degree,
        solar_path, offset_basis_count=offset_basis_count, offset_half_width_angstrom=offset_half_width_angstrom,
    )
    tensor.setflags(write=False)
    return tensor


def solar_lsf_tensor(
    wave: np.ndarray,
    knot_vector: np.ndarray,
    degree: int,
    solar_path: str | Path,
    offset_basis_count: int = 11, offset_half_width_angstrom: float = 3.0,
) -> np.ndarray:
    """Return read-only ``G`` with axes (blue output pixel, M-spline, B-spline).

    Memoised on the exact grid, knots, degree and solar file identity; the
    result is shared and must not be modified.
    """
    path = Path(solar_path).resolve()
    stat = path.stat()
    return _solar_tensor_cached(
        np.ascontiguousarray(wave, dtype=float).tobytes(),
        np.ascontiguousarray(knot_vector, dtype=float).tobytes(),
        int(degree),
        str(path),
        stat.st_size,
        stat.st_mtime_ns, offset_basis_count, offset_half_width_angstrom,
    )


class SkyDecompPalaceCorrSolarBlueLSF(SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D):
    """PALACE-corrected sky model with solar absorption constraining the blue LSF.

    Defaults: 11 offset by 4 wavelength coefficients in every channel, support
    +/-2.5 Angstrom, and a 0.01 second-difference penalty along offset only.
    Stop when W50/W90 and data chi2 change by <1%, and the barycenter by
    <0.015 Angstrom at every wavelength; at most 30 refinement cycles.
    """

    _close_each_refinement = True
    _lsf_numerical_ridge = 0.0

    def __init__(self, wave: np.ndarray, *args: Any, solar_blue_n_basis: int | None = None,
                 lsf_convergence_rtol: float = 0.01, chi2_convergence_rtol: float = 0.01,
                 **kwargs: Any) -> None:
        self.lsf_convergence_rtol = float(lsf_convergence_rtol)
        self.chi2_convergence_rtol = float(chi2_convergence_rtol)
        if not all(np.isfinite(x) and x > 0 for x in (self.lsf_convergence_rtol, self.chi2_convergence_rtol)):
            raise ValueError("convergence tolerances must be finite and positive")
        config = kwargs.get("config") or LSFSurfaceIterativeConfig(
            n_refinement_cycles=30, n_basis=4, roughness_fraction=0.0, fallback_prior_fraction=0.0)
        if config.blue_fit_lower == LSFSurfaceIterativeConfig().blue_fit_lower:
            # The continuum carries blue information everywhere, not only [OI] 5577.
            config = replace(config, blue_fit_lower=float(np.asarray(wave)[0]))
        kwargs["config"] = config
        kwargs.setdefault("offset_half_width_angstrom", 2.5)
        kwargs.setdefault("offset_roughness_fraction", 0.01)
        blue = LSFChannelSplineConfig(n_basis=config.n_basis if solar_blue_n_basis is None else solar_blue_n_basis,
            degree=3, knot_strategy="uniform")
        kwargs["spline_config"] = replace(
            _resolved_spline_config(config, kwargs.get("spline_config")), b=blue
        )
        self._solar_blue_stash: tuple[np.ndarray, np.ndarray] | None = None
        self._solar_pair: tuple[np.ndarray, np.ndarray] | None = None
        # The nominal seed solve keeps the DRP-Gaussian Moon/Zodi templates.
        super().__init__(wave, *args, **kwargs)
        self._blue = _channel_mask(self.wave, *_BLUE[1:])
        self._blue_fit = self._blue & (self.wave >= self.config.blue_fit_lower)

    def _run_iterations(self, flux, ivar):
        self.convergence_history = []
        self.converged = False
        self._previous_lsf_diagnostics = None
        self._convergence_started = time.perf_counter()
        return super()._run_iterations(flux, ivar)

    def _refinement_converged(self, run, flux, ivar, previous_surface):
        # Close the continuum/line solve on the NEW LSF before measuring chi2.
        fit, coefficient, error, continuum, weights, noise = self._fit_continuum_stage(
            run, flux, ivar, self._convergence_skyline_mask)
        run.continuum_status = str(fit["status"])
        if run.continuum_status not in self._SOLVED:
            run.failure = f"cycle_{self.lsf_surface_state.completed_cycles}_closure_continuum:{run.continuum_status}"
            return True
        design = self._stack_matrices(run.matrices, self._LINE_KEYS)
        line_fit = self._fit_design(design, flux - continuum, ivar)
        run.line_status = str(line_fit["status"])
        if run.line_status not in self._SOLVED:
            run.failure = f"cycle_{self.lsf_surface_state.completed_cycles}_closure_lines:{run.line_status}"
            return True
        run.continuum_coefficient, run.continuum_coef_err = coefficient, error
        run.coef_cov_moon, run.coef_cov_zodi = fit.get("coef_cov_moon"), fit.get("coef_cov_zodi")
        run.line_coefficient = np.asarray(line_fit["coef"], dtype=float)
        run.line_coef_err = np.asarray(line_fit.get("coef_err", np.full_like(run.line_coefficient, np.nan)))
        run.continuum, run.line_model = continuum, design.T @ run.line_coefficient
        run.weights, run.channel_noise, run.solver_status = weights, noise, run.line_status
        valid = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)
        chi2 = float(np.sum((flux[valid] - continuum[valid] - run.line_model[valid]) ** 2 * ivar[valid]))
        current = evaluate_lsf_diagnostics(self.lsf_surface_state, self.wave)
        previous = self._previous_lsf_diagnostics
        changes, centroid_changes = {}, {}
        for channel, lower, upper in LSF_CHANNELS:
            mask = _channel_mask(self.wave, lower, upper)
            changes[channel] = {key: (np.nan if previous is None else float(np.max(
                np.abs(current[key][mask] - previous[key][mask]) / previous[key][mask])))
                for key in ("w50_angstrom", "w90_angstrom")}
            centroid_changes[channel] = (np.nan if previous is None else float(np.max(
                np.abs(current["centroid_angstrom"][mask] - previous["centroid_angstrom"][mask]))))
        lsf_change = float(np.max([value for channel in changes.values() for value in channel.values()]))
        centroid_change = float(np.max(list(centroid_changes.values())))
        self._previous_lsf_diagnostics = current
        previous_chi2 = self.convergence_history[-1]["chi2"] if self.convergence_history else np.nan
        chi2_change = abs(chi2 - previous_chi2) / max(abs(previous_chi2), 1e-30)
        solved = all(m["status"] in self._SOLVED for m in self.lsf_surface_state.metrics.values())
        self.converged = bool(solved and np.isfinite(lsf_change) and np.isfinite(chi2_change)
            and np.isfinite(centroid_change) and centroid_change < 0.015
            and lsf_change < self.lsf_convergence_rtol and chi2_change < self.chi2_convergence_rtol)
        metric = dict(cycle=self.lsf_surface_state.completed_cycles, chi2=chi2,
            lsf_relative_change=lsf_change, chi2_relative_change=float(chi2_change),
            channel_width_relative_change=changes, centroid_change_angstrom=centroid_change,
            channel_centroid_change_angstrom=centroid_changes, converged=self.converged,
            elapsed_seconds=time.perf_counter() - self._convergence_started)
        self.convergence_history.append(metric)
        run.metrics.append(dict(stage=f"{metric['cycle']:02d}_convergence", solver_status=run.solver_status, **metric))
        return self.converged

    def _solar_tensor(self, knots: np.ndarray, degree: int) -> np.ndarray:
        return solar_lsf_tensor(self.wave, knots, degree, self.solar_path,
            self.offset_basis_count, self.offset_half_width_angstrom)

    def _envelopes(self) -> tuple[np.ndarray, np.ndarray]:
        """Solar-free Moon and Zodi envelope rows on blue pixels, from the current row's templates."""
        if self._solar_pair is None:
            self._solar_pair = _solar_templates(self._require_path(self.solar_path), self.wave, self.lsf_sigma)
        solar_hr, solar_rb = (values[self._blue] for values in self._solar_pair)
        # Zodi templates use the median-normalised convolved solar; undo that unit change.
        zodi_scale = solar_rb * float(np.nanmedian(self._solar_pair[0]))
        moon = np.divide(self.matrix_moon_hr[:, self._blue], solar_hr, out=np.zeros((self.matrix_moon_hr.shape[0], solar_hr.size)), where=solar_hr > 0.0)
        zodi = np.divide(self.matrix_zodi_hr[:, self._blue], zodi_scale, out=np.zeros((self.matrix_zodi_hr.shape[0], solar_hr.size)), where=zodi_scale > 0.0)
        return moon, zodi

    def _fit_continuum_stage(self, run, flux, ivar, skyline_mask):
        self._convergence_skyline_mask = skyline_mask
        result = super()._fit_continuum_stage(run, flux, ivar, skyline_mask)
        coefficient = result[1]
        self._solar_blue_stash = None
        if coefficient is not None:
            n_moon = run.matrices["moon"].shape[0]
            n_zodi = run.matrices["zodi"].shape[0]
            # Exact Moon + Zodi part that the LSF stage subtracted as fixed background.
            restored = (
                run.matrices["moon"].T @ coefficient[:n_moon]
                + run.matrices["zodi"].T @ coefficient[n_moon : n_moon + n_zodi]
            )
            self._solar_blue_stash = (np.asarray(coefficient, dtype=float).copy(), restored)
        return result

    def _extra_kernel_design(self, channel, fit_mask, knot_vector, degree):
        if channel != "B" or self._solar_blue_stash is None:
            return None
        coefficient, restored = self._solar_blue_stash
        moon, zodi = self._envelopes()
        n_moon = moon.shape[0]
        continuum = coefficient[:n_moon] @ moon + coefficient[n_moon : n_moon + zodi.shape[0]] @ zodi
        tensor = self._solar_tensor(knot_vector, degree)
        # Blue pixels are a prefix of the grid, so fit pixels index the tensor directly.
        keep = fit_mask[self._blue]
        design = continuum[keep, None] * tensor[keep].reshape(keep.sum(), -1)
        return design, restored[fit_mask]

    def _fit_lsf_channels(self, flux, ivar, source, fixed_background):
        state = super()._fit_lsf_channels(flux, ivar, source, fixed_background)
        metric = state.metrics["B"]
        if "extra_information_fraction" in metric:
            metric["model"] = SOLAR_BLUE_MODEL
        return state

    def _assemble_refined_matrices(self) -> dict[str, np.ndarray]:
        matrices = super()._assemble_refined_matrices()
        state = self.lsf_surface_state
        moon, zodi = self._envelopes()
        tensor = self._solar_tensor(state.knot_vectors["B"], state.degrees["B"])
        # Moon and Zodi in the blue use the same exact solar renderer as the LSF fit.
        profile = tensor.reshape(tensor.shape[0], -1) @ state.coefficients["B"].reshape(-1)
        for key, envelope in (("moon", moon), ("zodi", zodi)):
            matrices[key] = matrices[key].copy()
            matrices[key][:, self._blue] = envelope * profile
        return matrices

    def _finalize_result(self, *args: Any, **kwargs: Any):
        result = super()._finalize_result(*args, **kwargs)
        result.fit_summary += f" | blue_lsf=solar_continuum_joint_n{self.spline_config.b.n_basis}"
        result.fit_summary += f" | convergence={'converged' if self.converged else 'not_converged'}"
        for channel in result.lsf_state.metrics:
            result.lsf_state.metrics[channel].update(converged=self.converged,
                convergence_lsf_rtol=self.lsf_convergence_rtol, convergence_chi2_rtol=self.chi2_convergence_rtol,
                convergence_centroid_atol_angstrom=0.015,
                convergence_lsf_metric="maximum_pointwise_relative_W50_W90")
        self.fit_summary = result.fit_summary
        return result


__all__ = ["SOLAR_BLUE_MODEL", "SkyDecompPalaceCorrSolarBlueLSF", "solar_lsf_tensor"]
