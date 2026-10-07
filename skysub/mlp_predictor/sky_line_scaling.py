"""Rescale the predicted sky emission lines on the spectrum being sky-subtracted.

Why this works
--------------
The prediction transfers the sky lines from the sky arms to the science
pointing.  The line SHAPES transfer well, but their BRIGHTNESS does not: the OH
emission differs by a few percent between the science and sky directions, and
across the IFU.  OH brightness is most of the flux-error tail.  A sky line,
though, is 2-3 A wide and stands out from whatever science signal lies under
it.  So its brightness can be measured on the very spectrum it will be
subtracted from, while the broad science continuum, which the sky model cannot
separate from the target, is left alone.

The model
---------
The sky lines are split into a few templates built from the PREDICTED
coefficients (:func:`line_templates`):

* one per OH vibrational band v' (each OH coefficient is the population of one
  upper level v', N', F', so a band is a group of coefficients);
* optionally one rotational-temperature tilt per band: the same lines weighted
  by (E_upper - <E_upper>), so a temperature difference between science and
  sky is a linear term;
* one per atomic line family (ATOM_Og [O I] 5577, ATOM_Or [O I] 6300/6364,
  ATOM_Na Na D, ...), and the O2 b band.

On the science spectrum, the residual ``r = observed - sky`` is modelled as

    HP(r) = sum_b  d_b * HP(T_b)  +  noise  +  science features

where ``HP`` is the same LINEAR high-pass (x minus its running mean over
``highpass_A``) applied to data and templates alike.  The science continuum is
broad, so HP removes it; any line leakage into the running mean happens
identically to the templates, so the model stays exact.  ``d_b`` is the
fractional brightness change of template b (``s_b = 1 + d_b``).

The fit is weighted least squares with two safeguards:

* robust (Huber) iteration on standardised residuals, so that narrow science
  features (stellar absorption lines, unmasked emission lines) are
  down-weighted instead of pulling the scales;
* a Gaussian prior on every ``d_b`` (``prior_sigma``), so that a template with
  little signal on this row stays at the prediction.

The correction ``sum_b d_b * T_b`` (unfiltered templates) is ADDED to the
predicted sky, just like :func:`~mlp_predictor.sky_arm_correction.sky_arm_residual_correction`.

Na D and K I are NOT rescaled by default (``exclude``): stars have the same
lines in absorption at the same wavelengths, so starlight under the sky line
biases its scale directly.  With a G-star spectrum injected at 5x the sky
continuum, rescaling Na D raised the chi2 in the Na D window from 0.16 to 5.3.
Left at the prediction it stays at 0.26 either way.

Measured (600 held-out rows of the palacecorr corpus, after the full-band
sky-arm correction, median single-fibre photon chi2, science lines masked):

    red >= 6000 A     0.459 -> 0.284      decomposition's own fit 0.725
    OH-line pixels    1.187 -> 0.203                              0.614
    full band         0.364 -> 0.260                              0.931

* The fitted scales match the science decomposition's own band brightness:
  the median error is 0.1-0.6%, against the 1.5-2.7% the prediction misses.
* A 5x G star moves the OH bands by <= 1.5%.
* The result barely depends on the settings: high-pass 10-50 A, prior
  0.03-0.3 and Huber on/off all agree to 0.003.
* The tilt helps on OH pixels (0.223 -> 0.203).  A single scale for all OH is
  clearly worse (0.293).
"""
from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np
from scipy.ndimage import uniform_filter1d

__all__ = ["oh_group_labels", "line_templates", "sky_line_scaling_correction"]

_E_SCALE = 1000.0          # cm^-1: the tilt regressor's unit of upper-level energy


def oh_group_labels(model):
    """``(v_upper, E_upper)`` per OH coefficient, in the basis order of ``model``.

    Rebuilds the grouping ``sky_decomp.fit._oh_line_catalog`` uses -- same
    table, same wavelength cut, same keys -- and reads the vibrational level and
    the upper-level energy (cm^-1) of each group.  Checked against the model's
    own line groups, so a change in the catalog cannot silently mislabel bands.
    """
    from sky_decomp.fit import CAP_WAVE, decode_hitran_id, read_static_table, vac_to_air
    oh = read_static_table(str(model._pmd_path("pmd_popmodel_OH.dat")))
    oh["wave"] = vac_to_air(np.asarray(oh["lam"], float) * 1e4)
    lo, hi = float(model.wave.min()), float(model.wave.max())
    oh = decode_hitran_id(oh[(oh["wave"] >= lo - CAP_WAVE) & (oh["wave"] <= hi + CAP_WAVE)])
    groups = oh.group_by(list(model.oh_group_keys)).groups
    v = np.asarray([int(g["v_upper"][0]) for g in groups])
    e = np.asarray([float(np.mean(g["Ei"])) for g in groups])
    ref = getattr(model, "_oh_line_groups", None)
    if ref is not None:
        if len(ref) != len(groups):
            raise ValueError(f"OH catalog has {len(groups)} groups, the model {len(ref)}")
        for (wl, _), g in zip(ref, groups):
            if not np.allclose(np.sort(wl), np.sort(np.asarray(g["wave"], float))):
                raise ValueError("OH catalog grouping does not match the model's line groups")
    return v, e


def line_templates(model, mats, coef, *, tilt=True, min_band_fraction=0.01,
                   atoms=True, o2=True, exclude=("ATOM_Na", "ATOM_K")):
    """Per-family sky-line templates from ONE row's predicted coefficients.

    ``model`` is the row's decomposer and ``mats`` its assembled matrices (as
    used for the reconstruction), so the templates carry the same LSF and
    telluric transmission as the predicted sky.  Returns ``{name: (n_wave,)}``
    in the units of the reconstruction.  Bands holding less than
    ``min_band_fraction`` of the predicted OH flux are merged into their
    neighbour, so that every template has lines to fit.  Families named in
    ``exclude`` get no template, so they keep their predicted brightness (see
    the module notes on Na D).
    """
    coef = np.asarray(coef, dtype=np.float64).ravel()
    sl = model._component_slices(mats)
    m_oh = np.asarray(mats["oh"], dtype=np.float64)
    c_oh = coef[sl["oh"]]
    v, e = oh_group_labels(model)
    if v.size != m_oh.shape[0]:
        raise ValueError(f"{v.size} OH labels for {m_oh.shape[0]} OH basis rows")
    flux_g = np.clip(c_oh, 0.0, None) * m_oh.sum(axis=1)       # flux carried by each group
    total = float(flux_g.sum())
    bands = sorted(set(v.tolist()))
    # Merge faint bands into the nearest brighter one.
    keep = [b for b in bands if total > 0 and flux_g[v == b].sum() >= min_band_fraction * total]
    if not keep:
        keep = bands[:1]
    assign = {b: min(keep, key=lambda k: (abs(k - b), k)) for b in bands}
    out = {}
    for k in keep:
        idx = np.flatnonzero(np.isin(v, [b for b in bands if assign[b] == k]))
        out[f"OH_v{k}"] = m_oh[idx].T @ c_oh[idx]
        if tilt:
            w = flux_g[idx]
            e_mean = float(np.sum(w * e[idx]) / w.sum()) if w.sum() > 0 else float(np.mean(e[idx]))
            out[f"OH_v{k}_tilt"] = m_oh[idx].T @ (c_oh[idx] * (e[idx] - e_mean) / _E_SCALE)
    if atoms:
        m_at = np.asarray(mats["atom"], dtype=np.float64)
        c_at = coef[sl["atom"]]
        for j, name in enumerate(model.atom_names):
            if c_at[j] > 0 and name not in exclude:
                out[name] = m_at[j] * c_at[j]
    if o2 and "o2" in sl:
        c_o2 = coef[sl["o2"]]
        if np.any(c_o2 > 0):
            out["O2"] = np.asarray(mats["o2"], dtype=np.float64).T @ c_o2
    return out


def _highpass(x, mask, width_px):
    """x minus its masked running mean: linear in x for a fixed mask."""
    m = mask.astype(np.float64)
    num = uniform_filter1d(np.where(mask, x, 0.0), width_px, mode="nearest")
    den = uniform_filter1d(m, width_px, mode="nearest")
    return np.where(mask, x - num / np.where(den > 0, den, 1.0), 0.0)


def sky_line_scaling_correction(
    wave,
    sci_observed,
    sky_model,
    templates: Mapping[str, np.ndarray],
    *,
    variance=None,
    mask=None,
    highpass_A: float = 25.0,
    prior_sigma: float = 0.1,
    huber_k: float = 2.0,
    n_iter: int = 8,
    line_fraction: float = 0.02,
    return_info: bool = False,
):
    """Correction to ADD to the predicted sky, from rescaling its sky lines.

    Parameters
    ----------
    wave : (n_wave,) wavelength grid in Angstrom.
    sci_observed : the science spectrum being sky-subtracted (one fibre, or here
        the science median), (n_wave,).
    sky_model : the sky that will be subtracted from it before this step: the
        prediction, plus the sky-arm correction when that is applied.
    templates : ``{name: (n_wave,)}`` from :func:`line_templates`, in the units
        of ``sky_model``.
    variance : photon variance of ``sci_observed`` (only its SHAPE matters: the
        noise scale is measured from the residual).  None = uniform.
    mask : True where pixels must not be used (the science-line mask).
    highpass_A : width of the running mean the high-pass removes.  Wide enough
        to keep a line and its wings (~10x the LSF FWHM), narrow enough that the
        science continuum under the lines is gone.
    prior_sigma : 1-sigma prior on each fractional brightness change.  The tilt
        terms use the same width per 1000 cm^-1 of upper-level energy.
    huber_k : Huber threshold, in robust standard deviations.
    line_fraction : pixels where the high-passed templates reach this fraction
        of their row's 99th percentile define the robust noise scale.
    return_info : also return ``{name: scale}``, its 1-sigma, and fit stats.
    """
    wave = np.asarray(wave, dtype=np.float64)
    y_obs = np.asarray(sci_observed, dtype=np.float64)
    sky = np.asarray(sky_model, dtype=np.float64)
    names = [k for k, t in templates.items() if np.any(np.asarray(t) != 0)]
    if not names:
        zero = np.zeros_like(sky)
        return (zero, dict(scales={}, sigma={}, n_pix=0)) if return_info else zero
    T = np.vstack([np.asarray(templates[k], dtype=np.float64) for k in names])
    if T.shape[1] != wave.size or y_obs.shape != wave.shape or sky.shape != wave.shape:
        raise ValueError("wave, sci_observed, sky_model and templates must share one grid")
    good = np.isfinite(y_obs) & np.isfinite(sky) & np.all(np.isfinite(T), axis=0)
    if mask is not None:
        good &= ~np.asarray(mask, dtype=bool)
    w0 = np.ones_like(wave)
    if variance is not None:
        var = np.asarray(variance, dtype=np.float64)
        good &= np.isfinite(var) & (var > 0)
        w0 = np.where(good, 1.0 / np.where(good, var, 1.0), 0.0)
        w0 /= np.median(w0[good]) if good.any() else 1.0
    px = max(3, int(round(highpass_A / float(np.median(np.diff(wave))))) | 1)
    y = _highpass(y_obs - sky, good, px)
    X = np.vstack([_highpass(t, good, px) for t in T])
    strength = np.sum(np.abs(X), axis=0)
    line_pix = good & (strength > line_fraction * np.percentile(strength[good], 99)) if good.any() else good
    n_line = int(line_pix.sum())
    if n_line < 3 * len(names):
        zero = np.zeros_like(sky)
        return (zero, dict(scales={k: 1.0 for k in names}, sigma={k: np.inf for k in names},
                           n_pix=n_line)) if return_info else zero

    # Prior precision in data units: the residual's robust scale sets the noise
    # unit, so the prior's pull is the same whatever the flux units.
    hub = np.ones_like(wave)
    d = np.zeros(len(names))
    for _ in range(int(n_iter)):
        e = y - d @ X
        z = np.sqrt(w0) * e
        s = 1.4826 * np.median(np.abs(z[line_pix])) or 1.0
        u = np.abs(z) / s
        hub = huber_k / np.maximum(u, huber_k)          # 1 inside the threshold
        w = np.where(good, w0 * hub, 0.0) / s ** 2
        A = (X * w) @ X.T + np.eye(len(names)) / prior_sigma ** 2
        b = (X * w) @ y
        d_new = np.linalg.solve(A, b)
        if np.max(np.abs(d_new - d)) < 1e-5:
            d = d_new
            break
        d = d_new
    corr = d @ T
    if not return_info:
        return corr
    cov = np.linalg.inv(A)
    info = dict(scales={k: 1.0 + float(v) for k, v in zip(names, d)},
                sigma={k: float(np.sqrt(cov[i, i])) for i, k in enumerate(names)},
                n_pix=n_line, noise_scale=float(s),
                downweighted=float(np.mean(hub[line_pix] < 1.0)))
    return corr, info
