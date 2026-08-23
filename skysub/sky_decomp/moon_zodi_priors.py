"""Row-conditional amplitude and smoothness priors for the moon and zodi splines.

The core `sky_decomp.fit.SkyDecomp` QP already accepts scalar amplitude priors
via `moon_amp_prior`, `zodi_amp_prior` and their `_lambda` companions, evaluated
as ``lambda * (sum_i w_i c_i - T)^2``. What was missing was a way to derive
per-row targets `T` from the observation geometry, and a way to gate the
D²-smoothness penalty on `Moon_bs` so it grows when the geometry says the moon
should be near zero. This module supplies both.

Two physics-side proxies are exposed:

* :func:`kelsall_scaled_zodi_b500` returns the Leinert 1998 V-band lookup at the
  row's ecliptic longitude relative to the sun and its ecliptic latitude; this
  becomes the zodi amplitude target once multiplied by the corpus calibration.
* :func:`krisciunas_schaefer_moon_proxy` returns the Krisciunas & Schaefer 1991
  moon-scattered surface brightness in arbitrary units, zero when the moon is
  below the horizon; this becomes the moon amplitude target after calibration.

A calibration JSON glues the physics-side proxies to the data-side integrated
amplitudes ``sum_i w_i c_i`` recorded in an existing corpus. It is produced
once by ``skysub/calibrate_moon_zodi_priors.py`` and consumed at decomposition
time by ``decompose_parallel.py``.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, fields
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .moon_zodi_model import (
    MoonZodiInvalidObservationError,
    MoonZodiObservation,
    ZODIACAL_LIGHT_ASSET,
    _load_leinert,
    _interpolate_leinert,
    compute_midpoint_geometry,
)


@dataclass(frozen=True, slots=True)
class MoonZodiPriorTargets:
    """Per-row prior targets, ready to be passed to `SkyDecompLSFSurfaceIterative.fit`."""

    moon_amp_prior: float | None
    zodi_amp_prior: float | None
    moon_smooth_lambda_override: float | None
    zodi_smooth_lambda_override: float | None


@dataclass(frozen=True, slots=True)
class MoonZodiPriorCalibration:
    """Physics-to-data linear scale factors + smoothing-gate reference amplitudes.

    Attributes
    ----------
    k_moon, k_zodi
        Slopes of the robust linear fit ``target_data ~ k * proxy_physics``
        (no intercept). One-sample multipliers that convert the KS moon
        surface-brightness proxy and Leinert V-band lookup to units matched
        to the corpus's ``sum_i w_i c_i``.
    moon_ref_amp, zodi_ref_amp
        Reference amplitudes (typically the 20th percentile of the
        physics-side proxies on the calibration set) used by
        :func:`smooth_lambda_scale` to lift the smoothing penalty when a
        row's predicted amplitude is small.
    moon_smooth_gate_k, zodi_smooth_gate_k
        Strengths of the smoothing gates; see :func:`smooth_lambda_scale`.
    moon_amp_prior_lambda, zodi_amp_prior_lambda
        Default penalty strengths to pass with the targets; downstream may
        override.
    n_calibration_rows
        Number of rows used to derive `k_moon` / `k_zodi`.
    """

    k_moon: float
    k_zodi: float
    moon_ref_amp: float
    zodi_ref_amp: float
    moon_smooth_gate_k: float
    zodi_smooth_gate_k: float
    moon_amp_prior_lambda: float
    zodi_amp_prior_lambda: float
    n_calibration_rows: int
    moon_smooth_lambda_base: float
    zodi_smooth_lambda_base: float

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, path: str | Path) -> "MoonZodiPriorCalibration":
        payload = json.loads(Path(path).read_text())
        return cls(**payload)


def krisciunas_schaefer_moon_proxy(
    moon_altitude_deg: float,
    moon_separation_deg: float,
    signed_phase_deg: float,
    target_airmass: float,
    *,
    k_ext: float = 0.15,
) -> float:
    """Krisciunas & Schaefer 1991 moon-scattered brightness proxy (arb. units).

    Returns 0.0 when the moon is below the horizon or when any input is
    non-finite. Formula lifted from KS1991 eqs. 15-18 with a single generic
    V-band extinction coefficient `k_ext`; the absolute scale is calibrated
    against the corpus by :func:`fit_calibration`.

    Parameters
    ----------
    moon_altitude_deg
        Moon altitude above the observatory horizon.
    moon_separation_deg
        Angular separation between target and moon.
    signed_phase_deg
        Lunar phase angle (signed for KS's convention; only |phase| is used).
    target_airmass
        Airmass of the target line of sight.
    """
    if not (
        np.isfinite(moon_altitude_deg)
        and np.isfinite(moon_separation_deg)
        and np.isfinite(signed_phase_deg)
        and np.isfinite(target_airmass)
    ):
        return 0.0
    if moon_altitude_deg <= 0.0:
        return 0.0

    alpha = float(np.abs(signed_phase_deg))
    Istar = 10.0 ** (-0.4 * (3.84 + 0.026 * alpha + 4.0e-9 * (alpha ** 4)))

    rho = float(np.abs(moon_separation_deg))
    cos_rho = float(np.cos(np.deg2rad(rho)))
    f_rho = 10.0 ** 5.36 * (1.06 + cos_rho ** 2) + 10.0 ** (6.15 - rho / 40.0)

    zenith_moon = 90.0 - float(moon_altitude_deg)
    denom = np.cos(np.deg2rad(zenith_moon)) + 0.15 * (93.885 - zenith_moon) ** (-1.253)
    if denom <= 0.0:
        return 0.0
    X_moon = 1.0 / denom
    X_target = float(target_airmass)
    if not np.isfinite(X_target) or X_target <= 0.0:
        return 0.0

    return float(
        f_rho * Istar * 10.0 ** (-0.4 * k_ext * X_moon)
        * (1.0 - 10.0 ** (-0.4 * k_ext * X_target))
    )


def kelsall_scaled_zodi_b500(
    ecl_lon_relative_deg: float,
    ecl_latitude_deg: float,
    *,
    leinert_lookup=None,
) -> float:
    """Leinert 1998 V-band zodi surface brightness (S10_sun/deg^2) at (Δλ_ecl, β)."""
    if not (np.isfinite(ecl_lon_relative_deg) and np.isfinite(ecl_latitude_deg)):
        return 0.0
    try:
        return float(_interpolate_leinert(
            float(ecl_lon_relative_deg),
            float(ecl_latitude_deg),
            leinert_lookup,
        ))
    except Exception:
        return 0.0


def load_leinert_lookup(data_dir: str | Path):
    """Return the memoised Leinert lookup used by :func:`kelsall_scaled_zodi_b500`.

    ``data_dir`` may point at either the top-level ``sky_decomp/data`` folder
    or the ``moon_zodi/`` subfolder that actually hosts the asset; both are
    resolved to the same file.
    """
    data_dir = Path(data_dir)
    if data_dir.name != "moon_zodi" and (data_dir / "moon_zodi").is_dir():
        data_dir = data_dir / "moon_zodi"
    return _load_leinert(str(data_dir / ZODIACAL_LIGHT_ASSET))


def _resolve_moon_zodi_root(data_dir: str | Path) -> Path:
    """Return the ``moon_zodi/`` subfolder, whether ``data_dir`` is it or its parent."""
    data_dir = Path(data_dir)
    if data_dir.name == "moon_zodi":
        return data_dir
    if (data_dir / "moon_zodi").is_dir():
        return data_dir / "moon_zodi"
    return data_dir


def compute_geometry_proxies(
    observations: Iterable[MoonZodiObservation],
    *,
    data_dir: str | Path,
) -> dict[str, np.ndarray]:
    """Compute per-row moon and zodi physics-side proxies for an iterable of observations.

    Returns
    -------
    dict with keys ``'moon_amp_proxy'``, ``'zodi_b500'``, ``'valid'``, each a
    1D numpy array with the same length as ``observations``. Rows where the
    geometry call raises ``MoonZodiInvalidObservationError`` (target below
    horizon, near-sun zodi cell, etc.) are marked ``valid = False`` and given
    NaN proxies.
    """
    observations = list(observations)
    n_rows = len(observations)
    moon_arr = np.full(n_rows, np.nan, dtype=np.float64)
    zodi_arr = np.full(n_rows, np.nan, dtype=np.float64)
    valid = np.zeros(n_rows, dtype=bool)
    root = _resolve_moon_zodi_root(data_dir)
    lein = load_leinert_lookup(root)

    for i, obs in enumerate(observations):
        try:
            geom = compute_midpoint_geometry(obs, data_dir=root)
        except MoonZodiInvalidObservationError:
            continue
        moon_arr[i] = krisciunas_schaefer_moon_proxy(
            geom.moon_altitude_deg,
            geom.moon_separation_deg,
            geom.signed_phase_deg,
            geom.target_airmass,
        )
        zodi_arr[i] = kelsall_scaled_zodi_b500(
            geom.ecliptic_lon_relative_deg,
            geom.ecliptic_latitude_deg,
            leinert_lookup=lein,
        )
        valid[i] = np.isfinite(moon_arr[i]) and np.isfinite(zodi_arr[i])

    return {"moon_amp_proxy": moon_arr, "zodi_b500": zodi_arr, "valid": valid}


def robust_linear_slope(
    x: np.ndarray,
    y: np.ndarray,
    *,
    x_min_percentile: float = 20.0,
    trim_ratio_percentile: tuple[float, float] = (5.0, 95.0),
) -> float:
    """Return a robust no-intercept slope ``k`` in ``y = k * x``.

    Estimated as the trimmed median of ``y / x`` on rows where ``x`` is above
    the ``x_min_percentile``-th percentile (to keep the ratio well-defined
    when many rows have near-zero physics-side proxies). ``trim_ratio_percentile``
    trims outlier ratios inside the retained subset before the median.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size != y.size:
        raise ValueError('x and y must be the same length')
    finite = np.isfinite(x) & np.isfinite(y) & (x > 0.0)
    if not finite.any():
        return 0.0
    xf = x[finite]
    yf = y[finite]
    threshold = float(np.percentile(xf, x_min_percentile))
    keep = xf >= threshold
    if not keep.any():
        return float(np.median(yf / xf))
    ratios = yf[keep] / xf[keep]
    lo, hi = np.percentile(ratios, trim_ratio_percentile)
    trimmed = ratios[(ratios >= lo) & (ratios <= hi)]
    return float(np.median(trimmed if trimmed.size else ratios))


def fit_calibration(
    *,
    moon_proxy: np.ndarray,
    zodi_proxy: np.ndarray,
    moon_data_amp: np.ndarray,
    zodi_data_amp: np.ndarray,
    valid: np.ndarray | None = None,
    moon_smooth_gate_k: float = 100.0,
    zodi_smooth_gate_k: float = 0.0,
    moon_amp_prior_lambda: float = 1.0e-3,
    zodi_amp_prior_lambda: float = 1.0e-3,
    moon_smooth_lambda_base: float = 0.1,
    zodi_smooth_lambda_base: float = 0.1,
    moon_ref_percentile: float = 50.0,
    zodi_ref_percentile: float = 50.0,
) -> MoonZodiPriorCalibration:
    """Fit robust physics->data slopes and derive smoothing-gate reference amplitudes.

    ``moon_data_amp[i] = sum_j w_j c_{moon, j}`` from an existing decomposition;
    ``moon_proxy[i]`` is the KS output on the same row.  Same convention for zodi.

    The reference amplitude for the smoothness gate is taken as the given
    percentile of the STRICTLY-POSITIVE physics proxies, so a corpus with a
    lot of moon-below-horizon rows (KS proxy = 0) still yields a positive
    reference — the gate acts smoothly on rows below the typical "moon-up"
    amplitude and is inactive on rows near or above it.
    """
    valid = (
        np.ones_like(moon_proxy, dtype=bool)
        if valid is None
        else np.asarray(valid, dtype=bool)
    )
    m_ok = valid & np.isfinite(moon_proxy) & np.isfinite(moon_data_amp)
    z_ok = valid & np.isfinite(zodi_proxy) & np.isfinite(zodi_data_amp)
    k_moon = robust_linear_slope(moon_proxy[m_ok], moon_data_amp[m_ok])
    k_zodi = robust_linear_slope(zodi_proxy[z_ok], zodi_data_amp[z_ok])

    def _pos_pct(arr, pct):
        pos = arr[np.isfinite(arr) & (arr > 0.0)]
        if pos.size == 0:
            return 0.0
        return float(np.percentile(pos, pct))

    moon_ref_amp = _pos_pct(moon_proxy[m_ok], moon_ref_percentile)
    zodi_ref_amp = _pos_pct(zodi_proxy[z_ok], zodi_ref_percentile)

    return MoonZodiPriorCalibration(
        k_moon=float(k_moon),
        k_zodi=float(k_zodi),
        moon_ref_amp=moon_ref_amp,
        zodi_ref_amp=zodi_ref_amp,
        moon_smooth_gate_k=float(moon_smooth_gate_k),
        zodi_smooth_gate_k=float(zodi_smooth_gate_k),
        moon_amp_prior_lambda=float(moon_amp_prior_lambda),
        zodi_amp_prior_lambda=float(zodi_amp_prior_lambda),
        n_calibration_rows=int(m_ok.sum()),
        moon_smooth_lambda_base=float(moon_smooth_lambda_base),
        zodi_smooth_lambda_base=float(zodi_smooth_lambda_base),
    )


def smooth_lambda_scale(row_proxy: float, ref_proxy: float, gate_k: float) -> float:
    """Return the smoothing-gate scale factor for one row.

    ``factor = 1 + gate_k * relu(1 - proxy_row / ref_proxy)^2``. When the
    row's proxy is at or above the reference, the factor is 1.0. When it is
    zero (moon below horizon, for example), the factor is ``1 + gate_k``.
    """
    if gate_k <= 0.0 or not np.isfinite(row_proxy) or ref_proxy <= 0.0:
        return 1.0
    r = float(row_proxy) / float(ref_proxy)
    deficit = max(0.0, 1.0 - r)
    return 1.0 + float(gate_k) * deficit * deficit


def resolve_row_targets(
    *,
    moon_proxy: float,
    zodi_proxy: float,
    calibration: MoonZodiPriorCalibration,
) -> MoonZodiPriorTargets:
    """Convert row physics proxies + calibration into per-row prior targets."""
    if calibration.k_moon > 0.0 and np.isfinite(moon_proxy):
        moon_amp_prior = float(calibration.k_moon * moon_proxy)
    else:
        moon_amp_prior = None
    if calibration.k_zodi > 0.0 and np.isfinite(zodi_proxy):
        zodi_amp_prior = float(calibration.k_zodi * zodi_proxy)
    else:
        zodi_amp_prior = None

    moon_lambda_scale = smooth_lambda_scale(
        moon_proxy, calibration.moon_ref_amp, calibration.moon_smooth_gate_k
    )
    zodi_lambda_scale = smooth_lambda_scale(
        zodi_proxy, calibration.zodi_ref_amp, calibration.zodi_smooth_gate_k
    )

    return MoonZodiPriorTargets(
        moon_amp_prior=moon_amp_prior,
        zodi_amp_prior=zodi_amp_prior,
        moon_smooth_lambda_override=(
            calibration.moon_smooth_lambda_base * moon_lambda_scale
            if moon_lambda_scale != 1.0 else None
        ),
        zodi_smooth_lambda_override=(
            calibration.zodi_smooth_lambda_base * zodi_lambda_scale
            if zodi_lambda_scale != 1.0 else None
        ),
    )
