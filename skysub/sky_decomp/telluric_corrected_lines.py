"""Telluric-corrected line strengths on the continuous 2-D LSF model."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .fit import grp2vector, sticks2vector
from .lsf_spline2d import SkyDecompLSFSpline2D
from .moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_OH_SUFFIX,
    validate_decomposition_data_root,
)


TELLURIC_CORRECTED_LINES_FIT_MODEL = "telluric-corrected-lines-lsf-spline2d"


def _positive_scalar(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def calculate_drp_transmission(
    wave: np.ndarray,
    all_fiber_lsf: np.ndarray,
    pwv_mm: float,
    reduction_airmass: float,
    telluric_calculator: Any,
) -> np.ndarray:
    """Reproduce the native-pixel transmission used by the DRP."""
    wave = np.asarray(wave, dtype=float)
    lsf = np.asarray(all_fiber_lsf)
    pwv_mm = _positive_scalar(pwv_mm, "pwv_mm")
    reduction_airmass = _positive_scalar(reduction_airmass, "reduction_airmass")
    if wave.ndim != 1 or np.any(~np.isfinite(wave)) or np.any(np.diff(wave) <= 0.0):
        raise ValueError("wave must be finite and strictly increasing")
    if lsf.ndim != 2 or lsf.shape[1] != wave.size:
        raise ValueError("all_fiber_lsf must have shape (n_fibers, n_wave)")
    lsf_median = np.nanmedian(lsf, axis=0)
    if np.any(~np.isfinite(lsf_median)) or np.any(lsf_median <= 0.0):
        raise ValueError("the all-fiber DRP median LSF must be finite and positive")
    transmission = np.asarray(
        telluric_calculator.match_to_data(
            wave,
            lsf_median,
            pwv_mm,
            airmass=reduction_airmass,
            lsf_in_wavelength=True,
        ),
        dtype=float,
    )
    if transmission.shape != wave.shape:
        raise ValueError("DRP transmission does not match the native wavelength grid")
    if np.any(~np.isfinite(transmission)) or np.any(transmission <= 0.0):
        raise ValueError("DRP transmission must be finite and strictly positive")
    if np.any(transmission > 1.00001):
        raise ValueError("DRP transmission exceeds unity")
    return transmission


def restore_drp_input(
    flux_corrected: np.ndarray,
    ivar_corrected: np.ndarray,
    wave: np.ndarray,
    all_fiber_lsf: np.ndarray,
    pwv_mm: float,
    sci_airmass: float,
    telluric_calculator: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Undo the DRP pixel-space telluric division without changing the grid."""
    flux = np.asarray(flux_corrected, dtype=float)
    ivar = np.asarray(ivar_corrected, dtype=float)
    wave = np.asarray(wave, dtype=float)
    if wave.ndim != 1 or flux.shape != wave.shape or ivar.shape != wave.shape:
        raise ValueError("flux, ivar, and wave must be conformable one-dimensional arrays")
    if np.any(~np.isfinite(flux)) or np.any(~np.isfinite(ivar)) or np.any(ivar <= 0.0):
        raise ValueError("flux and positive ivar must be finite on every native pixel")
    transmission = calculate_drp_transmission(
        wave,
        all_fiber_lsf,
        pwv_mm,
        sci_airmass,
        telluric_calculator,
    )
    return flux * transmission, ivar / transmission**2, transmission


class SkyDecompTelluricLinesLSFSpline2D(SkyDecompLSFSpline2D):
    """Attenuate every configured sky line at its line-centre transmission."""

    def __init__(
        self,
        *args: Any,
        telluric_calculator: Any,
        pwv_mm: float,
        line_airmass: float,
        **kwargs: Any,
    ) -> None:
        self.telluric_calculator = telluric_calculator
        self.pwv_mm = _positive_scalar(pwv_mm, "pwv_mm")
        self.line_airmass = _positive_scalar(line_airmass, "line_airmass")
        self._transmission_hr: np.ndarray | None = None
        base_dir = Path(kwargs.get("base_dir") or DEFAULT_DATA_ROOT).resolve()
        if (base_dir / "bundle_manifest.json").is_file():
            validate_decomposition_data_root(str(base_dir))
        kwargs["base_dir"] = base_dir
        super().__init__(*args, **kwargs)

        self._intrinsic_groups = {
            "oh": [(wave.copy(), amp.copy()) for wave, amp in self._oh_line_groups],
            "atom": [(wave.copy(), amp.copy()) for wave, amp in self._atom_line_groups],
            "orc": [(wave.copy(), amp.copy()) for wave, amp in self._orc_line_groups],
        }
        self._transmission_hr = np.asarray(
            self.telluric_calculator.calc_transmission(
                self.pwv_mm,
                airmass=self.line_airmass,
            ),
            dtype=float,
        )
        if self._transmission_hr.shape != np.asarray(
            self.telluric_calculator.wave_air
        ).shape:
            raise ValueError("high-resolution transmission has an unexpected shape")
        if np.any(~np.isfinite(self._transmission_hr)) or np.any(
            self._transmission_hr < 0.0
        ):
            raise ValueError("high-resolution transmission must be finite and nonnegative")

        self._line_transmission(self._line_wave)
        for family in ("oh", "atom", "orc"):
            self._rebuild_family(family)
        self.design_matrix = self._assemble_design_matrix()

    def _line_transmission(self, line_wave_air: np.ndarray) -> np.ndarray:
        line_wave_air = np.asarray(line_wave_air, dtype=float)
        if self._transmission_hr is None:
            return np.ones_like(line_wave_air)
        transmission = np.interp(
            line_wave_air,
            np.asarray(self.telluric_calculator.wave_air, dtype=float),
            self._transmission_hr,
            left=np.nan,
            right=np.nan,
        )
        if np.any(~np.isfinite(transmission)) or np.any(transmission <= 0.0):
            raise ValueError(
                "the transmission model does not cover every retained line with positive transmission"
            )
        return transmission

    def _line_weights(self) -> np.ndarray:
        return super()._line_weights() * self._line_transmission(self._line_wave)

    def _rebuild_family(self, family: str) -> None:
        groups = self._intrinsic_groups[family]
        matrix = np.zeros((len(groups), self.wave.size), dtype=float)
        matrix_stick = np.zeros_like(matrix)
        for index, (line_wave, intrinsic_amplitude) in enumerate(groups):
            amplitude = intrinsic_amplitude * self._line_transmission(line_wave)
            matrix[index] = grp2vector(line_wave, amplitude, self.wave, self.lsf_sigma)
            matrix_stick[index] = sticks2vector(line_wave, amplitude, self.wave)
        setattr(self, f"matrix_{family}", matrix)
        setattr(self, f"matrix_{family}_stick", matrix_stick)

    def _prefit_o2(self, flux: np.ndarray, ivar: np.ndarray) -> None:
        intrinsic_aij = self.aij_o2.copy()
        self.aij_o2 = intrinsic_aij * self._line_transmission(self.lam_o2)
        try:
            super()._prefit_o2(flux, ivar)
        finally:
            self.aij_o2 = intrinsic_aij


class SkyDecompTelluricCorrectedLinesLSFSpline2D(
    SkyDecompTelluricLinesLSFSpline2D
):
    """Fit a DRP-corrected spectrum with line-centre telluric attenuation."""

    def __init__(
        self,
        *args: Any,
        telluric_calculator: Any,
        pwv_mm: float,
        source_airmass: float,
        drp_transmission: np.ndarray,
        **kwargs: Any,
    ) -> None:
        wave = np.asarray(args[0] if args else kwargs.get("wave"), dtype=float)
        transmission = np.asarray(drp_transmission, dtype=float)
        if transmission.shape != wave.shape:
            raise ValueError("drp_transmission must match the native wavelength grid")
        if np.any(~np.isfinite(transmission)) or np.any(transmission <= 0.0):
            raise ValueError("drp_transmission must be finite and strictly positive")
        self.drp_transmission = transmission.copy()
        self.source_airmass = _positive_scalar(source_airmass, "source_airmass")
        super().__init__(
            *args,
            telluric_calculator=telluric_calculator,
            pwv_mm=pwv_mm,
            line_airmass=self.source_airmass,
            **kwargs,
        )
        inverse = (1.0 / self.drp_transmission)[:, None]
        self._line_components = {
            channel: matrix.multiply(inverse).tocsc()
            for channel, matrix in self._line_components.items()
        }

    def _rebuild_family(self, family: str) -> None:
        super()._rebuild_family(family)
        setattr(
            self,
            f"matrix_{family}",
            getattr(self, f"matrix_{family}") / self.drp_transmission,
        )

    def _prefit_o2(self, flux: np.ndarray, ivar: np.ndarray) -> None:
        moon = self.vector_moon.copy()
        self.vector_moon = moon * self.drp_transmission
        try:
            super()._prefit_o2(
                np.asarray(flux) * self.drp_transmission,
                np.asarray(ivar) / self.drp_transmission**2,
            )
        finally:
            self.vector_moon = moon
        self.vector_o2 /= self.drp_transmission
        self.matrix_o2 /= self.drp_transmission
        self.o2_prefit_bestfit /= self.drp_transmission


class SkyDecompAdam25kTelluricLSFSpline2D(
    SkyDecompTelluricCorrectedLinesLSFSpline2D
):
    """Fit Adam-25k OH ratios with tellurics and a continuous 2-D LSF."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        requested_suffix = kwargs.get("palace_oh_suffix")
        if requested_suffix not in (None, DEFAULT_PALACE_OH_SUFFIX):
            raise ValueError("Adam-25k requires its bundled OH source table")
        kwargs["palace_oh_suffix"] = DEFAULT_PALACE_OH_SUFFIX
        super().__init__(*args, **kwargs)


__all__ = [
    "SkyDecompAdam25kTelluricLSFSpline2D",
    "SkyDecompTelluricCorrectedLinesLSFSpline2D",
    "SkyDecompTelluricLinesLSFSpline2D",
    "TELLURIC_CORRECTED_LINES_FIT_MODEL",
    "calculate_drp_transmission",
    "restore_drp_input",
]
