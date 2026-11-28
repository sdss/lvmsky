"""Science [O I] 6300 of single exposures from a saved ``oi6300.OI6300Model``: the preparatory step of SkyML.

``table`` builds, from the time, the three pointings and the decomposition coefficients of the three medians, the
input table that ``oi6300_build.build_table`` builds for a corpus: the observing context as SkyML computes it for any
exposure (``inference.build_triplet_from_pointings``), the ionospheric context (``iono_corpus.exposure_columns``),
the lines, the control-line ratios, and the night field.

The night field needs every exposure of the night, so it is read from the model's corpus table by exposure number
(``OI6300Model.corpus``).  A model with the night field gives no value for an exposure that was not in its corpus;
use a model built with ``--no-sky-field`` for such exposures.
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord

from . import inference, iono_corpus, oi6300
from .oi6300_build import O2_COLUMN

PAIRS = (('sci', 'near'), ('sci', 'far'), ('near', 'far'))


def table(model, obstime_mjd, expnum, radec, coef, coef_names, near_is_skye, download=True):
    """Input table of ``model`` for n exposures.

    obstime_mjd, expnum (n,); radec: arm -> (ra, dec) [deg], each (n,), for the arms ``sci``, ``near``, ``far``;
    coef: arm -> (n, n_coef) decomposition coefficients of the median of that telescope, columns ``coef_names``;
    near_is_skye (n,) bool: the near sky telescope is SkyE.  ``download`` fetches missing IONEX maps and GNSS
    station-days.  The column ``in_corpus`` tells whether the exposure is in the model's corpus table.
    """
    mjd, expnum = np.atleast_1d(np.asarray(obstime_mjd, float)), np.atleast_1d(expnum)
    pos = {a: tuple(np.broadcast_to(np.asarray(x, float), mjd.shape) for x in radec[a]) for a in oi6300.ARMS}
    label = np.where(np.broadcast_to(near_is_skye, mjd.shape), 'SKYE', 'SKYW')
    t = inference.build_triplet_from_pointings(
        obstime_mjd=mjd, sci_ra=pos['sci'][0], sci_dec=pos['sci'][1], sky_near_ra=pos['near'][0],
        sky_near_dec=pos['near'][1], sky_far_ra=pos['far'][0], sky_far_dec=pos['far'][1], sky_near_label=label,
        sky_far_label=np.where(label == 'SKYE', 'SKYW', 'SKYE'), verbose=False)
    names = [str(n) for n in t['ctx_names']]
    ctx = {a: np.asarray(t[f'ctx_{a}'], float) for a in oi6300.ARMS}
    col = lambda a, c: ctx[a][:, names.index(c)]  # noqa: E731
    alt = {a: col(a, 'alt') for a in ctx}
    az = {a: np.degrees(np.arctan2(col(a, 'az_sin'), col(a, 'az_cos'))) % 360.0 for a in ctx}
    tab = iono_corpus.exposure_columns(mjd, alt, az, download)
    tab['expnum'] = expnum
    cn = [str(n) for n in coef_names]
    c = {a: np.asarray(coef[a], float).reshape(mjd.size, len(cn)) for a in oi6300.ARMS}
    for a in oi6300.ARMS:
        tab |= {f'{a}_ctx_{k}': col(a, k) for k in model.ctx_names}
        tab |= {f'{a}_coef_{line}': c[a][:, cn.index(column)] for column, line in iono_corpus.LINES.items()}
    sky = {a: SkyCoord(*pos[a], unit='deg') for a in oi6300.ARMS}
    for a, b in PAIRS:
        tab[f'sep_{a}_{b}'] = sky[a].separation(sky[b]).deg
    v = tab['sci_vanrhijn_95km'], tab['near_vanrhijn_95km']
    for line, columns in (('oh', list(model.oh_groups)), ('o2', [O2_COLUMN])):
        if line in model.sci_lines:
            grab = lambda a: c[a][:, [cn.index(k) for k in columns]]  # noqa: E731, B023
            tab[f'sci_near_log_{line}'] = oi6300.control_ratio(grab('sci'), grab('near'), *v)
    row = model.corpus_row(expnum) if model.corpus is not None else np.full(mjd.size, -1)
    tab['in_corpus'] = row >= 0
    if model.sky_field:
        field = lambda k: model.corpus[k][np.maximum(row, 0)] if tab['in_corpus'].any() else np.nan  # noqa: E731
        tab |= {k: np.where(tab['in_corpus'], field(k), np.nan) for k in oi6300.SKY_FIELD}
    return tab


def predict(model, **exposures):
    """Prediction of ``model`` for the exposures of ``table(model, **exposures)``.

    Returns ``coef`` (ATOM_Or of the Science telescope; NaN where the model gives no value), ``log10_R`` (slant
    brightness), ``out_of_range`` and ``missing`` (inputs outside the training range, inputs not finite),
    ``lines_used`` (name of each guarded Science line -> used or not), ``in_corpus`` and the input ``table``.
    """
    tab = table(model, **exposures)
    log_r, flag = model.predict(tab)
    if model.sky_field:  # without its night an exposure is outside what the model was built for
        log_r = np.where(tab['in_corpus'], log_r, np.nan)
    used = model.lines_used(tab)
    return {'log10_R': log_r, 'coef': 10**log_r / oi6300.K['6300'], 'out_of_range': flag,
            'missing': model.quality(tab)['missing'], 'in_corpus': tab['in_corpus'],
            'lines_used': {n: used[:, i] for i, n in enumerate(model.guarded)}, 'table': tab}
