"""High-level inference API: predict SCI sky spectrum from minimal per-row inputs.

Wraps

- context construction from raw pointings/time (using astropy),
- ensemble prediction via
  :func:`mlp_predictor.trainer.predict_sci_coefficients_default`,
- per-row confidence via the ensemble spread across seeds,

so callers can go from ``(mjd, sci_ra, sci_dec, sky_near_ra, sky_near_dec,
sky_far_ra, sky_far_dec, coef_near, coef_far)`` straight to a predicted
science-arm coefficient vector plus a per-row confidence estimate, without
constructing FITS files or having to know which of the 39 context features
enter which specialised branch of the model.
"""

from __future__ import annotations

from typing import Mapping, Any

import numpy as np
from astropy.table import Table

from . import data, config as _cfg_module
from .trainer import predict_sci_coefficients_default


# Default context feature list.  Matches the deployed training pipeline:
# 28 base features that expand to 39 after the two augment stages (ecliptic
# geometry and physics-prior moon-scatter proxies).
_DEFAULT_CONTEXT_COLUMNS = tuple(_cfg_module.DataConfig().context_columns)


def _as_1d_f64(arr, name):
    a = np.asarray(arr, dtype=np.float64).ravel()
    if a.size == 0:
        raise ValueError(f"{name} is empty.")
    return a


def _compute_moon_phase_deg(obstime_mjd: np.ndarray) -> np.ndarray:
    """Return per-row moon phase in degrees (0 = new moon, 180 = full moon).

    Matches the convention of the ``MOON_PHASE`` META column consumed by the
    training pipeline: the encoded angle folds through ``moon_fli =
    (1 - cos(phase))/2``, so ``phase = 180 - elongation_from_sun``.
    """
    from astropy.coordinates import get_body, get_sun
    from astropy.time import Time
    import astropy.units as u

    lco = data._lco_earth_location()
    time = Time(np.asarray(obstime_mjd, dtype=np.float64), format="mjd", scale="utc")
    sun = get_sun(time)
    moon = get_body("moon", time, location=lco)
    elong_deg = moon.separation(sun).to_value(u.deg)
    return (180.0 - elong_deg).astype(np.float64)


def build_triplet_from_pointings(
    *,
    obstime_mjd,
    sci_ra, sci_dec,
    sky_near_ra, sky_near_dec,
    sky_far_ra,  sky_far_dec,
    moon_phase_deg=None,
    sky_near_label=None,
    sky_far_label=None,
    context_columns=None,
) -> dict:
    """Build a triplet dict with 39-dim per-arm ctx from minimal inputs.

    All ``ra`` / ``dec`` are in degrees (ICRS).  ``obstime_mjd`` is UT MJD.
    All arrays are broadcast to a common length (per-row); scalars are
    accepted for single-row predictions.

    Parameters
    ----------
    obstime_mjd:
        UT MJD of the observation.  Everything time-related derives from this
        (cyclic time features, solar-activity indices from the F10.7/Kp cache,
        moon/sun ephemerides).
    sci_ra, sci_dec:
        Pointing of the science fiber (degrees).
    sky_near_ra, sky_near_dec, sky_far_ra, sky_far_dec:
        Pointings of the two sky arms (degrees).  Which arm is "near" vs "far"
        is determined by the caller — the training pipeline defines "near" as
        the sky arm at smaller angular separation from the science pointing.
    moon_phase_deg:
        Optional per-row moon phase in degrees using the convention
        ``0 = new moon`` / ``180 = full moon``.  If omitted, computed from
        astropy sun/moon ephemerides.
    sky_near_label, sky_far_label:
        Optional per-row string labels 'SKYE' or 'SKYW' selecting the
        east/west assignment for the ``ew`` context feature.  Defaults to
        ('SKYE', 'SKYW') — this only affects the sign of a single ctx dim
        and its impact on predictions is small.
    context_columns:
        Override the default 28-feature base list.  Must be a subset of what
        ``mlp_predictor.data`` knows how to compute; see
        :attr:`DataConfig.context_columns` for the deployed default.

    Returns
    -------
    dict
        Keys: ``ctx_names``, ``ctx_sci`` / ``ctx_near`` / ``ctx_far`` (each
        of shape ``(n_rows, 39)``), and ``obstime_mjd``.  Ready to hand to
        :func:`predict_sky_from_minimal_inputs` or, together with the two
        sky-arm ``coef`` vectors, to
        :func:`~mlp_predictor.trainer.predict_sci_coefficients_default`
        directly.
    """
    obstime_mjd = _as_1d_f64(obstime_mjd, "obstime_mjd")
    n_rows = obstime_mjd.size
    _bcast = lambda a, name: np.broadcast_to(_as_1d_f64(a, name), (n_rows,)).astype(np.float64).copy()

    sci_ra_arr  = _bcast(sci_ra,  "sci_ra")
    sci_dec_arr = _bcast(sci_dec, "sci_dec")
    n1_ra  = _bcast(sky_near_ra,  "sky_near_ra")
    n1_dec = _bcast(sky_near_dec, "sky_near_dec")
    f1_ra  = _bcast(sky_far_ra,   "sky_far_ra")
    f1_dec = _bcast(sky_far_dec,  "sky_far_dec")

    if moon_phase_deg is None:
        moon_phase_arr = _compute_moon_phase_deg(obstime_mjd)
    else:
        moon_phase_arr = _bcast(moon_phase_deg, "moon_phase_deg")

    if sky_near_label is None:
        near_lbl = np.array(["SKYE"] * n_rows)
    else:
        near_lbl = np.asarray(sky_near_label).astype(str)
        if near_lbl.size == 1:
            near_lbl = np.repeat(near_lbl, n_rows)
    if sky_far_label is None:
        far_lbl = np.array(["SKYW"] * n_rows)
    else:
        far_lbl = np.asarray(sky_far_label).astype(str)
        if far_lbl.size == 1:
            far_lbl = np.repeat(far_lbl, n_rows)

    meta = Table({
        "OBSTIME":         obstime_mjd,
        "SCI_RA":          sci_ra_arr,
        "SCI_DEC":         sci_dec_arr,
        "SKY_NEAR_RA":     n1_ra,
        "SKY_NEAR_DEC":    n1_dec,
        "SKY_FAR_RA":      f1_ra,
        "SKY_FAR_DEC":     f1_dec,
        "MOON_PHASE":      moon_phase_arr,
        "SKY_NEAR_LABEL":  near_lbl,
        "SKY_FAR_LABEL":   far_lbl,
    })

    ctx_cols = list(context_columns) if context_columns is not None else list(_DEFAULT_CONTEXT_COLUMNS)
    ctx_sci,  ctx_names   = data._build_context_matrix(meta, ctx_cols, "sci")
    ctx_near, ctx_names_a = data._build_context_matrix(meta, ctx_cols, "sky1")
    ctx_far,  ctx_names_b = data._build_context_matrix(meta, ctx_cols, "sky2")
    if ctx_names != ctx_names_a or ctx_names != ctx_names_b:
        raise RuntimeError(
            "Context feature ordering diverged between arms — this is a bug in "
            "_build_context_matrix and should never happen with matching context_columns.")

    triplet: dict = {
        "ctx_names":   list(ctx_names),
        "ctx_sci":     np.asarray(ctx_sci,  dtype=np.float32),
        "ctx_near":    np.asarray(ctx_near, dtype=np.float32),
        "ctx_far":     np.asarray(ctx_far,  dtype=np.float32),
        "obstime_mjd": obstime_mjd,
        # Per-arm RA/Dec needed by the ecliptic augment.
        "sci_ra":  sci_ra_arr,  "sci_dec":  sci_dec_arr,
        "near_ra": n1_ra,       "near_dec": n1_dec,
        "far_ra":  f1_ra,       "far_dec":  f1_dec,
    }
    data._augment_triplet_with_ecliptic(triplet, force=True)
    data._augment_triplet_with_physics_priors(triplet, force=True)
    return triplet


def predict_sky_from_minimal_inputs(
    ensemble: Mapping[str, Any],
    *,
    obstime_mjd,
    sci_ra, sci_dec,
    sky_near_ra, sky_near_dec,
    sky_far_ra,  sky_far_dec,
    coef_near, coef_far,
    moon_phase_deg=None,
    sky_near_label=None, sky_far_label=None,
    return_per_seed: bool = False,
) -> dict:
    """Predict the science-arm decomposition coefficients + a confidence score.

    Parameters
    ----------
    ensemble:
        The dict returned by :func:`~mlp_predictor.serialization.load_ensemble`.
        Must have ``is_ensemble=True``.
    obstime_mjd, sci_ra, sci_dec, sky_near_ra, sky_near_dec, sky_far_ra, sky_far_dec:
        The minimal per-row observation metadata; see
        :func:`build_triplet_from_pointings`.
    coef_near, coef_far:
        The two sky-arm decomposition coefficient vectors, shape ``(n_rows,
        n_coef)`` where ``n_coef`` matches ``len(ensemble['coef_names'])``.
        These come from running the QP decomposition on the two sky-arm
        spectra separately.
    moon_phase_deg, sky_near_label, sky_far_label:
        Optional overrides passed straight to :func:`build_triplet_from_pointings`.
    return_per_seed:
        If True, additionally include the ``(n_seeds, n_rows, n_coef)`` array
        of per-seed predictions under key ``'per_seed'``.

    Returns
    -------
    dict
        Keys:

        - ``coef``:        ``(n_rows, n_coef)`` predicted SCI coefficient
          vector (ensemble mean, with the Jensen post-training lift applied).
        - ``coef_std``:    ``(n_rows, n_coef)`` per-coefficient std across
          the 10 ensemble members (epistemic uncertainty, physical units).
        - ``confidence``:  ``(n_rows,)`` scalar confidence in ``(0, 1]``,
          computed as ``1 / (1 + median_k(coef_std_k / |coef_k|))``.
          Higher = more agreement across seeds = more confident.
        - ``coef_names``:  the list of coefficient names matching the column
          order of ``coef`` and ``coef_std``.
        - ``triplet``:     the built context dict (useful for debugging).
        - ``per_seed`` (optional): raw per-seed predictions.

    Notes
    -----
    * The "confidence" summary is a simple monotone rescaling of the median
      relative std across coefficient dimensions.  For detailed downstream
      logic use ``coef_std`` directly.
    * All Jensen post-training bias corrections are applied inside
      :func:`predict_sci_coefficients_default`, so ``coef`` is already the
      final calibrated prediction.
    """
    if not ensemble.get("is_ensemble", False):
        raise ValueError("`ensemble` must be an ensemble artifact (is_ensemble=True).")

    triplet = build_triplet_from_pointings(
        obstime_mjd=obstime_mjd,
        sci_ra=sci_ra, sci_dec=sci_dec,
        sky_near_ra=sky_near_ra, sky_near_dec=sky_near_dec,
        sky_far_ra=sky_far_ra, sky_far_dec=sky_far_dec,
        moon_phase_deg=moon_phase_deg,
        sky_near_label=sky_near_label, sky_far_label=sky_far_label,
    )

    coef_near_arr = np.asarray(coef_near, dtype=np.float32)
    coef_far_arr  = np.asarray(coef_far,  dtype=np.float32)
    n_coef = len(ensemble["coef_names"])
    for name, arr in (("coef_near", coef_near_arr), ("coef_far", coef_far_arr)):
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.shape[1] != n_coef:
            raise ValueError(
                f"{name} has {arr.shape[1]} columns but the ensemble expects "
                f"{n_coef} (coef_names length).")
    if coef_near_arr.ndim == 1:
        coef_near_arr = coef_near_arr.reshape(1, -1)
    if coef_far_arr.ndim == 1:
        coef_far_arr = coef_far_arr.reshape(1, -1)

    ctx_near = triplet["ctx_near"]
    ctx_far  = triplet["ctx_far"]
    ctx_sci  = triplet["ctx_sci"]

    per_seed = []
    for member in ensemble["members"]:
        pred = predict_sci_coefficients_default(
            member,
            coef_near_phys=coef_near_arr, coef_far_phys=coef_far_arr,
            ctx_near_phys=ctx_near, ctx_far_phys=ctx_far, ctx_sci_phys=ctx_sci,
        )
        per_seed.append(np.asarray(pred, dtype=np.float64))
    per_seed_arr = np.stack(per_seed, axis=0)  # (n_seeds, n_rows, n_coef)

    coef_mean = per_seed_arr.mean(axis=0).astype(np.float32)
    coef_std  = per_seed_arr.std(axis=0, ddof=1).astype(np.float32)

    # Row-level scalar confidence: 1 / (1 + median_k(std_k / |mean_k|)).
    _denom = np.maximum(np.abs(coef_mean), 1e-12)
    _rel_std = coef_std / _denom
    _row_med = np.median(_rel_std, axis=1)
    confidence = (1.0 / (1.0 + _row_med)).astype(np.float32)

    result = {
        "coef":         coef_mean,
        "coef_std":     coef_std,
        "confidence":   confidence,
        "coef_names":   list(ensemble["coef_names"]),
        "triplet":      triplet,
    }
    if return_per_seed:
        result["per_seed"] = per_seed_arr.astype(np.float32)
    return result


__all__ = [
    "build_triplet_from_pointings",
    "predict_sky_from_minimal_inputs",
]
