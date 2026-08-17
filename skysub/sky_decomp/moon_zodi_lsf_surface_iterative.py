"""Moon/Zodi-aware iterative LVM sky decomposition."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import time
import tracemalloc
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

from .result_io import load_moon_zodi_state
from .lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    LSFSurfaceIterativeResult,
    SkyDecompLSFSurfaceIterative,
    _FitRun,
)
from .moon_zodi_model import (
    DEFAULT_DATA_DIR,
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_DIFFUSE_SUFFIX,
    DEFAULT_PALACE_OH_SUFFIX,
    MoonZodiObservation,
    MoonZodiPhysicalModel,
    MoonZodiState,
    validate_decomposition_data_root,
)


@dataclass(slots=True)
class MoonZodiLSFSurfaceIterativeResult(LSFSurfaceIterativeResult):
    """Iterative result with one compact Moon/Zodi provenance state."""

    moon_zodi_state: MoonZodiState


class SkyDecompMoonZodiLSFSurfaceIterative(SkyDecompLSFSurfaceIterative):
    """Fit one shared B-spline correction to a physical Moon plus Zodi model."""

    def __init__(
        self,
        wave: np.ndarray,
        *,
        physical_to_fit_flux_scale: float,
        config: LSFSurfaceIterativeConfig | None = None,
        data_root: str | Path = DEFAULT_DATA_ROOT,
        moon_zodi_data_dir: str | Path | None = None,
        **kwargs,
    ) -> None:
        if (
            not np.isfinite(physical_to_fit_flux_scale)
            or physical_to_fit_flux_scale <= 0.0
        ):
            raise ValueError("physical_to_fit_flux_scale must be positive and finite")
        wave_value = np.asarray(wave)
        if wave_value.dtype != np.float64:
            raise ValueError("Moon/Zodi wavelength grid must be float64")
        data_root_value = Path(data_root).expanduser().resolve()
        validate_decomposition_data_root(str(data_root_value))
        physical_data_dir = (
            Path(moon_zodi_data_dir).expanduser().resolve()
            if moon_zodi_data_dir is not None
            else data_root_value / DEFAULT_DATA_DIR.name
        )
        if kwargs.get("base_dir") is None:
            kwargs["base_dir"] = data_root_value
        if kwargs.get("palace_oh_suffix") is None:
            kwargs["palace_oh_suffix"] = DEFAULT_PALACE_OH_SUFFIX
        if kwargs.get("palace_diffuse_suffix") is None:
            kwargs["palace_diffuse_suffix"] = DEFAULT_PALACE_DIFFUSE_SUFFIX
        self.physical_to_fit_flux_scale = float(physical_to_fit_flux_scale)
        self.data_root = data_root_value
        self.physical_model = MoonZodiPhysicalModel(data_dir=physical_data_dir)
        self.physical_model.validate_wave(wave_value)
        self._moon_zodi_matrix_pairs: dict[
            int,
            tuple[np.ndarray, np.ndarray],
        ] = {}
        self._prediction_state: MoonZodiState | None = None
        self.physical_moon_prediction = np.array([], dtype=np.float64)
        self.physical_zodi_prediction = np.array([], dtype=np.float64)
        self.moon_correction_basis = np.empty((0, 0), dtype=np.float64)
        self.moon_correction_knots = np.array([], dtype=np.float64)
        super().__init__(wave_value, config=config, **kwargs)

    def _build_moon(self) -> tuple[np.ndarray, list[str]]:
        """Create the fixed cubic correction basis; physical carriers arrive at fit time."""
        w0, w1 = float(self.wave[0]), float(self.wave[-1])
        interior = self._uniform_moon_knots(w0, w1, self.n_spline_knots)
        full_knots = np.r_[(w0,) * 4, interior, (w1,) * 4]
        basis = np.asarray(
            BSpline.design_matrix(self.wave, full_knots, 3).toarray(),
            dtype=np.float64,
        )
        if (
            basis.shape[0] != self.wave.size
            or np.any(~np.isfinite(basis))
            or np.any(basis < 0.0)
            or not np.allclose(np.sum(basis, axis=1), 1.0, rtol=0.0, atol=1.0e-14)
        ):
            raise RuntimeError("Moon/Zodi correction B-spline basis is invalid")
        self.moon_knots_used = interior.copy()
        self.moon_correction_knots = full_knots.astype(np.float64, copy=True)
        self.moon_correction_basis = basis
        self.vector_moon = np.zeros_like(self.wave, dtype=np.float64)
        matrix = np.zeros((basis.shape[1], self.wave.size), dtype=np.float64)
        self.matrix_moon_hr = matrix.copy()
        names = [f"MoonZodi_bs{index:03d}" for index in range(matrix.shape[0])]
        return matrix, names

    def _install_prediction(
        self,
        observation: MoonZodiObservation,
        detector_lsf_fwhm: np.ndarray,
    ) -> None:
        prediction = self.physical_model.predict(
            self.wave,
            detector_lsf_fwhm,
            observation,
            physical_to_fit_flux_scale=self.physical_to_fit_flux_scale,
        )
        moon = np.asarray(prediction.moon, dtype=np.float64)
        zodi = np.asarray(prediction.zodi, dtype=np.float64)
        basis = self.moon_correction_basis
        matrix_moon = (basis * moon[:, None]).T
        matrix_zodi = (basis * zodi[:, None]).T
        matrix_total = matrix_moon + matrix_zodi
        self.physical_moon_prediction = moon.copy()
        self.physical_zodi_prediction = zodi.copy()
        self.vector_moon = moon + zodi
        self.matrix_moon = matrix_total
        # This method's base detector projection already includes the measured
        # LSF. The iterative 11-tap surface is an additional residual correction.
        self.matrix_moon_hr = matrix_total.copy()
        self.design_matrix = self._assemble_design_matrix()
        self._moon_zodi_matrix_pairs = {
            id(self.matrix_moon): (matrix_moon, matrix_zodi),
        }
        self._prediction_state = replace(
            prediction.state,
            correction_degree=3,
            correction_knots=tuple(float(value) for value in self.moon_correction_knots),
        )

    def _assemble_refined_matrices(self) -> dict[str, np.ndarray]:
        moon = self._convolve_matrix_channelwise(
            (self.moon_correction_basis * self.physical_moon_prediction[:, None]).T
        )
        zodi = self._convolve_matrix_channelwise(
            (self.moon_correction_basis * self.physical_zodi_prediction[:, None]).T
        )
        total = moon + zodi
        matrices = self._matrix_bundle(
            self._convolve_matrix_channelwise(self.matrix_oh_stick),
            total,
            self.matrix_diffuse,
            self._convolve_matrix_channelwise(self.matrix_atom_stick),
            self._convolve_matrix_channelwise(self.matrix_orc_stick),
            self._convolve_matrix_channelwise(self.matrix_o2_stick),
        )
        self._moon_zodi_matrix_pairs[id(total)] = (moon, zodi)
        return matrices

    def _components_from_coef(
        self,
        coef: np.ndarray,
        mats: dict[str, np.ndarray],
    ) -> dict[str, np.ndarray]:
        slices = self._component_slices(mats)
        try:
            matrix_moon, matrix_zodi = self._moon_zodi_matrix_pairs[id(mats["moon"])]
        except KeyError as error:
            raise RuntimeError("Moon/Zodi component matrices do not match committed state") from error
        correction = coef[slices["moon"]]
        moon = matrix_moon.T @ correction
        zodi = matrix_zodi.T @ correction
        moon_zodi_total = mats["moon"].T @ correction
        tolerance = 1.0e-10 * max(1.0, float(np.max(np.abs(moon_zodi_total))))
        if float(np.max(np.abs(moon + zodi - moon_zodi_total))) > tolerance:
            raise RuntimeError("Moon plus Zodi does not close to the shared design block")
        components = {
            "moon": moon,
            "zodi": zodi,
            "diffuse": mats["diffuse"].T @ coef[slices["diffuse"]],
            "oh": mats["oh"].T @ coef[slices["oh"]],
            "atom": mats["atom"].T @ coef[slices["atom"]],
            "orc": mats["orc"].T @ coef[slices["orc"]],
            "o2": mats["o2"].T @ coef[slices["o2"]],
        }
        return components

    def _continuum_from_components(
        self,
        components: dict[str, np.ndarray],
    ) -> np.ndarray:
        return components["moon"] + components["zodi"] + components["diffuse"]

    def _metric(
        self,
        stage: str,
        status: str,
        flux: np.ndarray,
        model: np.ndarray,
        ivar: np.ndarray,
        skyline_mask: np.ndarray,
        **extra: object,
    ) -> dict[str, object]:
        """Attach per-cycle kernel moments without changing baseline diagnostics."""
        metric = SkyDecompLSFSurfaceIterative._metric(
            stage,
            status,
            flux,
            model,
            ivar,
            skyline_mask,
            **extra,
        )
        state = self.lsf_surface_state
        if stage.endswith("_lsf_lines") and state is not None:
            names = (
                "center_pix_min",
                "center_pix_median",
                "center_pix_max",
                "sigma_pix_min",
                "sigma_pix_median",
                "sigma_pix_max",
            )
            metric["lsf_kernel_moments"] = {
                channel: {name: values.get(name, np.nan) for name in names}
                for channel, values in state.metrics.items()
            }
        return metric

    def _finalize_result(
        self,
        run: _FitRun,
        seed: dict[str, object],
        flux: np.ndarray,
        ivar: np.ndarray,
        skyline_mask: np.ndarray,
        started: float,
        trace_started: bool,
        verbose: bool,
    ) -> MoonZodiLSFSurfaceIterativeResult:
        base = super()._finalize_result(
            run,
            seed,
            flux,
            ivar,
            skyline_mask,
            started,
            trace_started,
            verbose,
        )
        if self._prediction_state is None:
            raise RuntimeError("Moon/Zodi prediction state was not initialized")
        expected_keys = ("moon", "zodi", "diffuse", "oh", "atom", "orc", "o2")
        if tuple(base.components) != expected_keys:
            raise RuntimeError(f"Unexpected Moon/Zodi component order: {tuple(base.components)}")
        closure = sum((base.components[name] for name in expected_keys), np.zeros_like(self.wave))
        tolerance = 1.0e-10 * max(1.0, float(np.max(np.abs(base.bestfit_lsf))))
        if float(np.max(np.abs(closure - base.bestfit_lsf))) > tolerance:
            raise RuntimeError("Moon/Zodi refined model failed full component closure")
        base_values = {
            field.name: getattr(base, field.name)
            for field in fields(LSFSurfaceIterativeResult)
        }
        return MoonZodiLSFSurfaceIterativeResult(
            **base_values,
            moon_zodi_state=self._prediction_state,
        )

    def fit(
        self,
        flux: np.ndarray,
        ivar: np.ndarray,
        *,
        observation: MoonZodiObservation,
        detector_lsf_fwhm: np.ndarray,
        verbose: bool = False,
    ) -> MoonZodiLSFSurfaceIterativeResult:
        """Predict Moon/Zodi and run the transactional five-stage LSF workflow."""
        flux_value = np.asarray(flux)
        ivar_value = np.asarray(ivar)
        lsf_value = np.asarray(detector_lsf_fwhm)
        if flux_value.dtype != np.float64 or ivar_value.dtype != np.float64:
            raise ValueError("Moon/Zodi flux and ivar must be float64")
        if flux_value.shape != self.wave.shape or ivar_value.shape != self.wave.shape:
            raise ValueError("flux and ivar must match the native wavelength grid")
        if lsf_value.dtype != np.float64 or lsf_value.shape != self.wave.shape:
            raise ValueError("detector_lsf_fwhm must be float64 and match wave")
        if not np.all(np.isfinite(lsf_value) & (lsf_value > 0.0)):
            raise ValueError("detector_lsf_fwhm must be finite and positive")
        if not np.all(np.isfinite(ivar_value) & (ivar_value >= 0.0)):
            raise ValueError("Moon/Zodi ivar must be finite and nonnegative")
        invalid_flux = ~np.isfinite(flux_value)
        if np.any(invalid_flux & (ivar_value > 0.0)):
            raise ValueError("non-finite Moon/Zodi flux samples must have zero ivar")

        trace_started = False
        if not tracemalloc.is_tracing():
            tracemalloc.start()
            trace_started = True
        tracemalloc.reset_peak()
        started = time.perf_counter()
        self._set_lsf_state(None)
        self._install_prediction(observation, lsf_value)
        self._prefit_o2(flux_value, ivar_value)
        seed, skyline_mask, run = self._run_iterations(flux_value, ivar_value)
        return self._finalize_result(
            run,
            seed,
            flux_value,
            ivar_value,
            skyline_mask,
            started,
            trace_started,
            verbose,
        )


__all__ = [
    "MoonZodiLSFSurfaceIterativeResult",
    "SkyDecompMoonZodiLSFSurfaceIterative",
    "load_moon_zodi_state",
]
