"""Night-long field of the zenith [O I] 6300 measured by the two sky telescopes (inputs of ``oi6300`` with
``sky_field=True``).

Every exposure of a night gives two samples of the emitting layer, at the 285 km pierce points of SkyNear and SkyFar.
A Gaussian process is fitted to all samples of the night and read at the Science and SkyNear pierce points of each
exposure.  Nothing measured on the Science spectrum enters.  The field needs the whole night, later exposures
included, so it is an input for post-processing only.

    v = log10 z_6300 = m + f(x, y + u t, t) + exposure term + telescope term + noise
    k = sT^2 exp(-dt^2 / 2 tauT^2)                        the whole field rising and falling through the night
      + sS^2 exp(-dr^2 / 2 l^2 - dt^2 / 2 tauS^2)         horizontal structure of the layer
      + sS2^2 exp(-dr^2 / 2 l2^2 - dt^2 / 2 tau2^2)       small-scale structure
      + sX^2 [same exposure] + sE^2 s_i s_j + s0^2        factor common to both telescopes of an exposure;
                                                          SkyE/SkyW offset of the night (s = +-1); night level

x, y   magnetic meridional / zonal pierce-point offsets from the observatory [km] (corpus ``<arm>_ipp_dx_mag_*``).
t      hours since 0 UT of the evening date.
u      zonal drift of the frame [km/h]: 0 (fixed to the Earth) or ``V_ROT`` (fixed in local time).

``fit`` learns the kernel parameters of the two frames from the corpus: they maximise the likelihood of the
sky-telescope samples of ``N_FIT_NIGHTS`` random nights.  They are the same for every night, so nothing is fitted on
a night with few exposures.  The likelihood constrains some terms weakly and they move between fits on
different nights: read the parameters as a fitted kernel, not as physics.
``oi6300_build`` fits them when it builds a model and stores them on it (``field_hyper``); ``night_field`` must be
given the parameters of the model it feeds.
Method from ``experiments_external_data_v1/WP8_lvm_internal/skyfield.py``.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

from .oi6300 import SKY_FIELD as COLUMNS
from .oi6300 import log_zenith

V_ROT = 2 * np.pi * (6371.0 + 285.0) * np.cos(np.radians(29.0)) / 24.0  # [km/h] at the latitude of LCO
S0, CLIP, MIN_SAMPLES, COUNT_WINDOW = 0.3, 4.0, 3, 1.0  # night level [dex]; clip [robust sigma]; samples; [h]
FRAMES = {'earth': 0.0, 'corot': V_ROT}  # frame -> u
PARAMS = ('sT', 'tauT', 'sS', 'l', 'tauS', 'sS2', 'l2', 'tau2', 'sE', 'sX', 'sn')  # [dex, h, dex, km, h, ...]
START = ((0.27, 1.4, 0.25, 590.0, 0.85, 0.02, 40.0, 2.0, 0.003, 0.054, 0.024),  # starting points of the fit
         (0.3, 1.0, 0.3, 1000.0, 2.0, 0.01, 15.0, 0.5, 0.01, 0.05, 0.02))
N_FIT_NIGHTS, MIN_FIT_EXPOSURES = 80, 8


def _spatial(p, a, b, u):
    """The two horizontal-structure terms between point sets a, b = (x, y, t, ...), and the offsets dx, dy."""
    dx, dy = a[0][:, None] - b[0][None], (a[1] + u * a[2])[:, None] - (b[1] + u * b[2])[None]
    dt2 = (a[2][:, None] - b[2][None]) ** 2
    terms = [p[s] ** 2 * np.exp(-0.5 * (dx**2 + dy**2) / p[ln] ** 2 - 0.5 * dt2 / p[tau] ** 2)
             for s, ln, tau in (('sS', 'l', 'tauS'), ('sS2', 'l2', 'tau2'))]
    return terms, dx, dy, dt2


def _kernel(p, a, b, u):
    """Covariance between a and b = (x, y, t, s, exposure index)."""
    terms, _, _, dt2 = _spatial(p, a, b, u)
    return (p['sT'] ** 2 * np.exp(-0.5 * dt2 / p['tauT'] ** 2) + terms[0] + terms[1]
            + p['sX'] ** 2 * (a[4][:, None] == b[4][None]) + p['sE'] ** 2 * a[3][:, None] * b[3][None] + S0**2)


def _clip(p, d, v, u):
    """Drop samples whose leave-one-out residual is beyond ``CLIP`` robust sigma (failed decompositions)."""
    for _ in range(2):
        if v.size < 5:
            break
        ki = np.linalg.inv(_kernel(p, d, d, u) + p['sn'] ** 2 * np.eye(v.size))
        z = (ki @ (v - np.median(v))) / np.sqrt(np.diag(ki))
        keep = np.abs(z) < CLIP * max(1.4826 * np.median(np.abs(z)), 1e-9)
        if keep.all():
            break
        d, v = [q[keep] for q in d], v[keep]
    return d, v


def _night(p, u, d, v, sci, near, t):
    """Field at Science minus field at SkyNear, and the SkyNear gradient times the Science-SkyNear offset [dex]."""
    al = cho_solve(cho_factor(_kernel(p, d, d, u) + p['sn'] ** 2 * np.eye(v.size), lower=True), v - np.median(v))
    z, none = np.zeros_like(t), np.full(t.size, -1)  # no telescope term, no exposure term at the read-out points
    k_sci, k_near = (_kernel(p, (*q, t, z, none), d, u) for q in (sci, near))
    terms, dx, dy, _ = _spatial(p, (*near, t), d, u)
    gx = -sum((dx / p[ln] ** 2 * k) @ al for k, ln in zip(terms, ('l', 'l2'), strict=True))
    gy = -sum((dy / p[ln] ** 2 * k) @ al for k, ln in zip(terms, ('l', 'l2'), strict=True))
    return (k_sci - k_near) @ al, gx * (sci[0] - near[0]) + gy * (sci[1] - near[1])


def _nights(tab, usable):
    """Per night: row indices, t [h], the samples (x, y, t, s, exposure index) and values, and the pierce points."""
    g = lambda c: np.asarray(tab[c], float)  # noqa: E731
    night, t_all = np.asarray(tab['night_mjd']), (g('obstime_mjd') - g('night_mjd')) * 24.0
    skye = np.asarray(tab['near_is_skye'], bool)
    xy = {a: (g(f'{a}_ipp_dx_mag_mer'), g(f'{a}_ipp_dx_mag_zon')) for a in ('sci', 'near', 'far')}
    val = {a: log_zenith(tab, a, '6300') for a in ('near', 'far')}
    for n in np.unique(night):
        r = np.flatnonzero(night == n)
        t = t_all[r]
        ok = {a: np.asarray(usable)[r] & np.isfinite(val[a][r] + xy[a][0][r] + xy[a][1][r]) for a in ('near', 'far')}
        cat = lambda f, ok=ok: np.concatenate([f(a)[ok[a]] for a in ('near', 'far')])  # noqa: E731
        d = [cat(lambda a, r=r: xy[a][0][r]), cat(lambda a, r=r: xy[a][1][r]), cat(lambda a, t=t: t),
             cat(lambda a, r=r: np.where(skye[r] == (a == 'near'), 1.0, -1.0)), cat(lambda a, r=r: np.arange(r.size))]
        yield r, t, d, cat(lambda a, r=r: val[a][r]), {a: [q[r] for q in xy[a]] for a in ('sci', 'near')}


def _nll(p, d, v, u):
    """Negative log likelihood of the samples of one night (constant dropped)."""
    if v.size < MIN_SAMPLES:
        return 0.0
    cf = cho_factor(_kernel(p, d, d, u) + p['sn'] ** 2 * np.eye(v.size), lower=True)
    vv = v - np.median(v)
    return 0.5 * vv @ cho_solve(cf, vv) + np.log(np.diag(cf[0])).sum()


def fit(tab, usable, seed=0, maxiter=4000):
    """Kernel parameters of both frames, frame -> dict, from ``N_FIT_NIGHTS`` random nights of ``tab`` that hold at
    least ``MIN_FIT_EXPOSURES`` usable exposures (all of them when fewer).  Two passes: all samples, then the samples
    left by the clipping with the first-pass parameters.  About 1 min per frame on 80 nights.

    The entry ``'fit'`` of the result records the nights and samples used and, per frame, the final negative log
    likelihood and whether the optimiser converged within ``maxiter``."""
    data = [(d, v) for _, _, d, v, _ in _nights(tab, usable) if np.unique(d[4]).size >= MIN_FIT_EXPOSURES]
    if not data:
        raise ValueError(f'no night with {MIN_FIT_EXPOSURES} usable exposures: the night field cannot be fitted')
    pick = np.random.default_rng(seed).choice(len(data), min(N_FIT_NIGHTS, len(data)), replace=False)
    data = [data[i] for i in pick]
    hyper = {'fit': {'nights': len(data), 'samples': int(sum(v.size for _, v in data))}}
    opt = dict(maxiter=maxiter, xatol=1e-3, fatol=0.1, adaptive=True)
    for frame, u in FRAMES.items():
        def total(q, sample, u=u):
            p = dict(zip(PARAMS, np.exp(q), strict=True))
            try:  # parameters far from any sensible kernel overflow or leave it singular: a point to avoid
                nll = sum(_nll(p, d, v, u) for d, v in sample)
            except (np.linalg.LinAlgError, ValueError):
                return 1e12
            return nll if np.isfinite(nll) else 1e12
        first = min((minimize(total, np.log(q0), args=(data,), method='Nelder-Mead', options=opt) for q0 in START),
                    key=lambda res: res.fun)
        clipped = [_clip(dict(zip(PARAMS, np.exp(first.x), strict=True)), d, v, u) for d, v in data]
        res = minimize(total, first.x, args=(clipped,), method='Nelder-Mead', options=opt)
        if not np.isfinite(res.fun) or res.fun >= 1e12:
            raise ValueError(f'night-field fit of the {frame} frame failed: likelihood {res.fun}')
        hyper[frame] = {k: float(x) for k, x in zip(PARAMS, np.exp(res.x), strict=True)}
        hyper['fit'][frame] = {'nll': float(res.fun), 'converged': bool(res.success)}
    return hyper


def night_field(tab, usable, hyper):
    """Field inputs for every row of ``tab`` (``COLUMNS``; NaN on a night with fewer than 3 usable samples).

    tab     corpus-style table holding whole nights: ``night_mjd``, ``obstime_mjd``, ``near_is_skye`` and, per arm,
            ``ipp_dx_mag_mer``, ``ipp_dx_mag_zon``, ``coef_6300``, ``vanrhijn_285km``.
    usable  (n_rows,) bool: exposures whose SkyNear and SkyFar medians can be trusted (enough fibres).
    hyper   kernel parameters from ``fit`` (those stored on the model that reads the inputs).
    ``field_n_1h`` is the number of samples within 1 h of the exposure that the field rests on.
    """
    out = {c: np.full(len(tab['night_mjd']), np.nan) for c in COLUMNS}
    for r, t, d0, v0, xy in _nights(tab, usable):
        for frame, u in FRAMES.items():
            d, v = _clip(hyper[frame], d0, v0, u)
            if v.size < MIN_SAMPLES:
                continue
            out[f'field_sci_near_{frame}'][r], out[f'field_grad_{frame}'][r] = _night(
                hyper[frame], u, d, v, xy['sci'], xy['near'], t)
            if frame == 'earth':
                out['field_n_1h'][r] = (np.abs(t[:, None] - d[2][None]) <= COUNT_WINDOW).sum(1)
    return out
