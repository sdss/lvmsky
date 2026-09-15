"""Frozen continuum configuration shared by the Niv-integrated methods."""

from __future__ import annotations

from typing import Any


NIV_CONTINUUM_KWARGS = {
    "n_spline_knots": 11,
    "moon_smooth_lambda": 1.0e-1,
    "split_zodi": True,
    "n_zodi_spline_knots": 1,
    "zodi_smooth_lambda": 1.0e-1,
    "moon_albedo_fiducial_phase_deg": 30.0,
    "zodi_color_exponent": 0.26,
    "moon_ratio_bound": 0.7,
    "zodi_ratio_bound": 0.7,
    "amp_prior_tol": 3.0,
    "amp_prior_floor": 0.02,
    "zodi_amp_bound": 2.0,
    "diffuse_ratio_bound_dex": 0.2,
    "diffuse_ratio_nominal": (0.0396, 0.7026, 0.2578),
    "diffuse_oh_centre_log10": -0.6489,
    "diffuse_oh_bound_dex": 0.15,
    "diffuse_oh_gate_frac": 0.6,
    "diffuse_oh_relax_dex": 0.0,
    "diffuse_oh_scope": "block",
    "moon_interline_boost": 0.0,
}
NIV_ZODI_PRIOR_CALIBRATION = 1.6


def apply_niv_continuum_contract(kwargs: dict[str, Any]) -> dict[str, Any]:
    """Install the validated continuum settings and reject silent drift."""
    kwargs = dict(kwargs)
    for key, expected in NIV_CONTINUUM_KWARGS.items():
        if key in kwargs:
            actual = kwargs[key]
            matches = (
                tuple(actual) == expected
                if isinstance(expected, tuple)
                else actual == expected
            )
            if not matches:
                raise ValueError(
                    f"Niv continuum method requires {key}={expected!r}, "
                    f"got {actual!r}"
                )
        kwargs[key] = expected
    return kwargs


__all__ = [
    "NIV_CONTINUUM_KWARGS",
    "NIV_ZODI_PRIOR_CALIBRATION",
    "apply_niv_continuum_contract",
]
