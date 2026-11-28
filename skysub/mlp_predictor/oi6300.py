"""[O I] 6300 along the Science line of sight, predicted from the two sky telescopes.

Model (the form of the SkyML ionospheric group: a blend of the sky arms plus a learned correction)::

    log10 I_sci = log10 V_sci + a log10 z_near + (1 - a) log10 z_far + delta(x)

I_sci   slant [O I] 6300 brightness of the Science telescope [R]; I = K c, c = decomposition coefficient ATOM_Or.
z_arm   zenith brightness of a sky telescope, K c_arm / V_arm.
V       van Rhijn factor of a thin layer at 285 km (corpus column ``<arm>_vanrhijn_285km``).
a       one scalar fitted on the training rows (minimum mean |residual|); 1.0 on the knox 1.3.2 corpus.
delta   gradient boosting with absolute loss on x: SkyNear/SkyFar lines, observing context and ionosphere of all
        three telescopes, and the geometry between them.  By default no quantity measured on the Science spectrum
        enters x.

Optional Science-spectrum airglow lines (``sci_lines``; never [O I] 6300 itself), as Science/SkyNear zenith log ratios:
``'oh'``, ``'o2'``  mesospheric control lines (``control_ratio``); a nebula does not emit them.  They carry the
                    brightness factor common to all lines that differs between the two telescopes.
``'5577'``, ``'7774'``, ``'8446'``, ``'5199'``  lines a nebula can emit; each is guarded.  Its ratio is predicted
                    from x and the control lines, and the excess over that prediction is expressed in units of its
                    test-night scatter (``line_levels``).  ``predict`` drops a line that is missing or whose excess is
                    above ``guard`` sigma and uses the model trained without it (one model per subset of the
                    guarded lines).

Optional night field of the sky telescopes (``sky_field=True``): the columns ``SKY_FIELD`` of the table, computed
by ``sky_field.night_field`` from every exposure of the night with the kernel parameters learned at build time
(``field_hyper``).  Post-processing only: the field of an exposure uses
the later exposures of its night.

Inputs are kept only if the training rows cover the application domain: >= ``coverage`` of the domain rows inside the
1-99% training range.  Outside that range trees return their edge value, so a prediction stays bounded on any
exposure; ``predict`` also returns the number of inputs outside the range as a quality flag.  A missing input stays
missing (NaN; the trees route it), ``quality`` counts them.  The prediction is NaN when the sky-telescope [O I] 6300
that the base rests on is not valid.
With ``fit(..., members=5)`` the correction is an ensemble of 5 models, each trained without one random group of
nights, and ``predict`` averages them.

``oi6300_build`` builds the input table from a decomposed median stack, applies the data cuts (``clean_rows``),
cross-validates by nights and saves a model with its provenance (``meta``) and its corpus table (``corpus``: the
night field and the prediction of every corpus exposure).  ``oi6300_exposure`` applies a saved model to single
exposures (SkyML).  The studies behind the method are in
``experiments_ionospheric_oi_model`` and ``experiments_external_data_v1`` (see ``EXPERIMENTS_CHANGELOG.md``).
"""

from __future__ import annotations

import pickle
from itertools import combinations

import numpy as np

K = {'6300': 1.298, '7774': 34.5, '8446': 15.46, '5199': 1.43}  # R per unit decomposition coefficient
ARMS = ('sci', 'near', 'far')
# Ionosphere and thermosphere at the 285 km pierce point (corpus columns <arm>_<name>); log10 taken where named.
IONO = ('log10_model6300', 'log10_model7774', 'h6300', 'wind_w', 'wind_mer_mag', 'vtec', 'hmf2', 'log10_nmf2',
        'b_bot', 'log10_o2', 'log10_o', 'tn', 'dlog10vtec_dt', 'dlog10vtec_dlat', 'dlog10vtec_dlon', 'dip_lat',
        'bubble_p', 'roti_max')
GEOM = ('ipp_lat', 'ipp_lon', 'local_time', 'airmass')
CLEAN_MIN_FIBRES, CLEAN_MIN_LOG_EW = 10, 1.2  # a median of >= 10 fibres, log10 EW([O I] 6300) > 1.2 [A]
SCI_CONTROL, SCI_GUARDED = ('oh', 'o2'), ('5577', '7774', '8446', '5199')  # allowed entries of sci_lines
ROTI_FLOOR = 0.005  # [TECU/min] floor of a measured ROTI before the logarithm; a missing ROTI stays missing
# Raise when ``inputs``, ``science_lines``, ``base`` or the night field change meaning: saved models must be rebuilt.
INPUTS_VERSION = 2
HGB = dict(loss='absolute_error', max_iter=400, learning_rate=0.03, max_leaf_nodes=15, min_samples_leaf=30,
           l2_regularization=1.0, early_stopping=False, random_state=0)
SKY_FIELD = ('field_sci_near_earth', 'field_grad_earth', 'field_sci_near_corot', 'field_grad_corot', 'field_n_1h')
HC, OMEGA = 6.62607015e-27 * 2.99792458e10, np.pi * (0.5 * 35.3 / 206265.0) ** 2  # erg cm; fibre solid angle [sr]


def log_zenith(tab, arm, line):
    """log10 zenith brightness [R] of a line; NaN outside 0.1-1e4 R (a line the fit drove to ~0)."""
    z = K[line] * np.asarray(tab[f'{arm}_coef_{line}'], float) / np.asarray(tab[f'{arm}_vanrhijn_285km'], float)
    return np.log10(np.where((z > 0.1) & (z < 1e4), z, np.nan))


def control_ratio(coef_sci, coef_near, v_sci, v_near):
    """Science/SkyNear zenith log10 ratio of a mesospheric line group: the median over lines of
    log10(c_sci / V_sci) - log10(c_near / V_near).  coef_* (n_rows,) or (n_rows, n_lines); V = van Rhijn at 95 km."""
    lz = lambda c, v: np.log10(np.where(c > 0, c, np.nan) / np.asarray(v, float)[:, None])  # noqa: E731
    c_s, c_n = (np.asarray(c, float).reshape(len(v_sci), -1) for c in (coef_sci, coef_near))
    return np.nanmedian(lz(c_s, v_sci) - lz(c_n, v_near), axis=1)


def science_lines(tab, names):
    """Science-spectrum inputs, name -> (n_rows, k); the last column is the Science/SkyNear zenith log ratio.

    Control lines are read from ``tab['sci_near_log_oh']`` / ``tab['sci_near_log_o2']`` (see ``control_ratio``)."""
    f = {}
    for n in names:
        if n in SCI_CONTROL:
            f[n] = np.asarray(tab[f'sci_near_log_{n}'], float)[:, None]
        elif n == '5577':  # 95 km layer
            f[n] = control_ratio(tab['sci_coef_5577'], tab['near_coef_5577'], tab['sci_vanrhijn_95km'],
                                 tab['near_vanrhijn_95km'])[:, None]
        else:
            z = log_zenith(tab, 'sci', n)
            f[n] = np.column_stack([z, z - log_zenith(tab, 'near', n)])
    return f


def _night_groups(tab, rows, n):
    night = np.asarray(tab['night_mjd'])[rows]
    u = np.unique(night)
    g_of = dict(zip(np.random.default_rng(0).permutation(u), np.arange(u.size) % n, strict=True))
    return np.array([g_of[x] for x in night])


def equivalent_width(wave, flux, coef_6300):
    """EW of [O I] 6300 [A]: decomposition line flux / median of the observed spectrum in 6200-6285, 6315-6400 A.

    wave [A] (n_pix,), flux [erg s-1 cm-2 A-1 per fibre] (n_rows, n_pix) on the native grid, coef_6300 (n_rows,).
    """
    side = ((wave > 6200) & (wave < 6285)) | ((wave > 6315) & (wave < 6400))
    line = K['6300'] * np.asarray(coef_6300, float) / (6300.3e-8 / HC * 4 * np.pi * 1e-6 / OMEGA)
    return line / np.nanmedian(np.asarray(flux, float)[:, side], axis=1)


def clean_rows(tab, n_fibres, log_ew):
    """Rows usable for training and validation: every arm a >= 10-fibre median, log10 EW > 1.2, all lines valid.

    n_fibres, log_ew: dicts arm -> (n_rows,).  Few-fibre medians in crowded fields are dominated by stellar
    continuum and their weak-line coefficients fail.
    """
    ok = True
    for a in ARMS:
        ok = ok & (np.asarray(n_fibres[a]) >= CLEAN_MIN_FIBRES) & (np.asarray(log_ew[a]) > CLEAN_MIN_LOG_EW)
        ok = ok & np.all([np.isfinite(log_zenith(tab, a, line)) for line in K], axis=0)
    return ok


def _arm(tab, a):
    g = lambda c: np.asarray(tab[f'{a}_{c}'], float)  # noqa: E731
    q = {c: g(c) for c in IONO + GEOM[:3]}
    q['vtec'], q['roti_max'] = np.log10(q['vtec']), np.log10(np.maximum(q['roti_max'], ROTI_FLOOR))  # NaN stays NaN
    q['airmass'] = 1.0 / np.sin(np.radians(g('alt')))
    return q


def inputs(tab, ctx_names, ni=True, sky_field=False):
    """Input dict x: name -> (n_rows,).  tab maps corpus-style names (``<arm>_<column>``, ``<arm>_ctx_<name>``,
    ``sep_<a>_<b>`` [deg], ``f107``, ``f107a``, ``kp``, ``sun_alt_conjugate``) to arrays.  ni=False drops [N I] 5199;
    sky_field=True adds the ``SKY_FIELD`` columns."""
    lines = [line for line in K if ni or line != '5199']
    f = {f'{a} log z{line}': log_zenith(tab, a, line) for a in ('near', 'far') for line in lines}
    f |= {f'far-near log z{line}': f[f'far log z{line}'] - f[f'near log z{line}'] for line in lines}
    q = {a: _arm(tab, a) for a in ARMS}
    for a in ARMS:
        f |= {f'{a} ctx {c}': np.asarray(tab[f'{a}_ctx_{c}'], float) for c in ctx_names}
        f |= {f'{a} {k}': v for k, v in q[a].items()}
    f |= {'sun_alt_conjugate': np.asarray(tab['sun_alt_conjugate'], float), 'f107': np.asarray(tab['f107'], float),
          'log10 f107a': np.log10(np.asarray(tab['f107a'], float)), 'kp': np.asarray(tab['kp'], float)}
    for a, b in (('sci', 'near'), ('sci', 'far'), ('near', 'far')):
        f[f'sep {a}-{b}'] = np.asarray(tab[f'sep_{a}_{b}'], float)
        f |= {f'{a}-{b} {c}': q[a][c] - q[b][c] for c in ('vtec', 'dip_lat') + GEOM}
    if sky_field:
        f |= {c.replace('_', ' '): np.asarray(tab[c], float) for c in SKY_FIELD}
    return f


class OI6300Model:
    """Science [O I] 6300 from SkyNear + SkyFar + context; see the module docstring."""

    sci_lines, guarded, guard, level, fallback = (), (), None, {}, {}  # defaults for models saved before sci_lines
    sky_field = False  # default for models saved before the night field
    oh_groups, meta = (), {}  # set by oi6300_build: OH groups behind tab['sci_near_log_oh']; provenance and scores
    field_hyper = None  # set by oi6300_build: kernel parameters of the night field (sky_field.fit)
    inputs_version = 1  # models saved before INPUTS_VERSION
    # set by oi6300_build, per corpus exposure: expnum, obstime_mjd, log10_pred (test-night prediction where
    # pred_from_cv, else this model's) and, with the night field, the SKY_FIELD columns
    corpus = None

    def __init__(self, ni=True, coverage=0.9, ctx_names=(), sci_lines=(), guard=3.0, sky_field=False):
        bad = set(sci_lines) - set(SCI_CONTROL + SCI_GUARDED)
        if bad or len(set(sci_lines)) != len(sci_lines):
            raise ValueError(f'sci_lines {list(sci_lines)}: unknown {sorted(bad)} or repeated entries')
        if guard is not None and not guard > 0:
            raise ValueError(f'guard must be positive or None, not {guard}')
        self.inputs_version = INPUTS_VERSION
        self.ni, self.coverage, self.ctx_names = ni, coverage, tuple(ctx_names)
        self.sci_lines, self.guard, self.sky_field = tuple(sci_lines), guard, sky_field
        self.guarded = tuple(n for n in self.sci_lines if n in SCI_GUARDED)

    def _matrix(self, tab):
        f = inputs(tab, self.ctx_names, self.ni, self.sky_field)
        Z = np.column_stack(list(f.values()))
        Z[~np.isfinite(Z)] = np.nan
        return list(f), Z

    def _x(self, Z, S, sub):
        """Model matrix: kept sky inputs, control lines, and the guarded Science lines in ``sub``."""
        return np.column_stack([Z[:, self.keep]] + [S[n] for n in self.sci_lines if n in SCI_CONTROL or n in sub])

    def base(self, tab, alpha=None):
        """Blend of the sky telescopes rescaled to the Science line of sight; a telescope of zero weight is not read."""
        a = self.alpha if alpha is None else alpha
        near, far = (log_zenith(tab, arm, '6300') for arm in ('near', 'far'))
        blend = near if a == 1 else far if a == 0 else a * near + (1 - a) * far
        return np.log10(np.asarray(tab['sci_vanrhijn_285km'], float)) + blend

    def fit(self, tab, rows, domain=None, members=1):
        """Train on boolean/index ``rows``; ``domain`` = rows of the application domain for the coverage cut.

        members > 1 trains an ensemble: the training nights (``tab['night_mjd']``) are split at random into
        ``members`` groups, one model is trained without each group, and ``predict`` averages them.
        """
        from sklearn.ensemble import HistGradientBoostingRegressor

        # night groups: one per ensemble member, and 5 for the level models of the guarded lines
        need = max(members if members > 1 else 0, 5 if self.guarded and self.guard is not None else 0)
        if members < 1 or (need and np.unique(np.asarray(tab['night_mjd'])[rows]).size < need):
            raise ValueError(f'members = {members}: needs members >= 1 and at least {need} training nights')
        y = np.log10(K['6300'] * np.asarray(tab['sci_coef_6300'], float))
        names, Z = self._matrix(tab)
        lo, hi = np.nanpercentile(Z[rows], [1, 99], axis=0)
        keep = np.isfinite(lo) & (hi > lo)
        if domain is not None:
            d = Z[domain]
            cov = ((d >= lo) & (d <= hi)).sum(0) / np.maximum(np.isfinite(d).sum(0), 1)
            keep &= cov >= self.coverage
        self.keep, self.names, self.lo, self.hi = np.flatnonzero(keep), names, lo, hi
        grid = np.linspace(0, 1, 101)
        self.alpha = grid[np.argmin([np.nanmean(np.abs(y[rows] - self.base(tab, a)[rows])) for a in grid])]
        S, t = science_lines(tab, self.sci_lines), (y - self.base(tab))[rows]
        if not np.all(np.isfinite(t)):
            raise ValueError('training rows without a valid [O I] 6300 in Science, SkyNear or SkyFar (see clean_rows)')
        hgb = lambda X, v: HistGradientBoostingRegressor(**HGB).fit(X, v)  # noqa: E731
        # members == 1: one model on every row
        group = _night_groups(tab, rows, members) if members > 1 else np.full(t.size, -1)

        def train(sub):
            X = self._x(Z, S, sub)[rows]
            return [hgb(X[group != k], t[group != k]) for k in range(members)]

        self.models, self.fallback, self.level = train(self.guarded), {}, {}
        if not self.guarded or self.guard is None:
            return self
        self.fallback = {sub: train(sub) for n in range(len(self.guarded)) for sub in combinations(self.guarded, n)}
        X, g5 = self._x(Z, S, ())[rows], _night_groups(tab, rows, 5)  # the same night groups as an ensemble of 5
        for line in self.guarded:  # level model of the line ratio; sigma from nights left out of the fit
            v = S[line][rows, -1]
            ok, r = np.isfinite(v), np.full(v.size, np.nan)
            for k in range(5):
                r[g5 == k] = v[g5 == k] - hgb(X[ok & (g5 != k)], v[ok & (g5 != k)]).predict(X[g5 == k])
            self.level[line] = (hgb(X[ok], v[ok]), 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r))))
        return self

    def line_levels(self, tab):
        """Excess of each guarded Science line over its prediction from the sky telescopes and the control lines,
        in units of its test-night scatter: name -> (n_rows,).  The scatter [dex] is ``self.level[name][1]``."""
        _, Z = self._matrix(tab)
        S = science_lines(tab, self.sci_lines)
        X = self._x(Z, S, ())
        return {n: (S[n][:, -1] - m.predict(X)) / s for n, (m, s) in self.level.items()}

    def lines_used(self, tab):
        """(n_rows, n_guarded) bool, columns in the order of ``self.guarded``: the guarded Science lines that
        ``predict`` uses for each exposure, those measured and not above ``guard`` sigma."""
        if not self.level:
            return np.ones((len(np.asarray(tab['night_mjd'])), len(self.guarded)), bool)
        return np.column_stack([np.isfinite(v) & ~(v > self.guard) for v in self.line_levels(tab).values()])

    def quality(self, tab):
        """Per-row counts over the kept inputs: ``out_of_range`` (outside the 1-99% training range, the flag of
        ``predict``) and ``missing`` (not finite; Science lines and controls included)."""
        _, Z = self._matrix(tab)
        X = self._x(Z, science_lines(tab, self.sci_lines), self.guarded)
        Zk = Z[:, self.keep]
        return {'out_of_range': ((Zk < self.lo[self.keep]) | (Zk > self.hi[self.keep])).sum(1),
                'missing': np.isnan(X).sum(1)}

    def predict(self, tab):
        """log10 slant Science [O I] 6300 [R] and the quality flag (inputs outside the 1-99% training range).

        Guarded Science lines that are missing or above ``guard`` sigma (``lines_used``) are not used for that
        exposure.  NaN where the sky-telescope [O I] 6300 of the base is not valid."""
        _, Z = self._matrix(tab)
        S = science_lines(tab, self.sci_lines)
        Zk = Z[:, self.keep]
        flag = ((Zk < self.lo[self.keep]) | (Zk > self.hi[self.keep])).sum(1)

        def mean(models, sub, sel):
            return np.mean([m.predict(self._x(Z, S, sub)[sel]) for m in models], axis=0)

        if not self.level:
            return self.base(tab) + mean(self.models, self.guarded, slice(None)), flag
        ok = self.lines_used(tab)  # a missing line goes to the model without it
        delta = np.full(len(Z), np.nan)
        for pat in np.unique(ok, axis=0):
            sel = (ok == pat).all(1)
            sub = tuple(n for n, b in zip(self.guarded, pat, strict=True) if b)  # the guarded lines kept
            delta[sel] = mean(self.models if sub == self.guarded else self.fallback[sub], sub, sel)
        return self.base(tab) + delta, flag

    def coefficient(self, tab):
        """Predicted decomposition coefficient ATOM_Or of the Science telescope."""
        return 10 ** self.predict(tab)[0] / K['6300']

    def corpus_row(self, expnum):
        """Row of ``self.corpus`` of each exposure number, -1 for an exposure that was not in the build corpus."""
        e, expnum = self.corpus['expnum'], np.atleast_1d(expnum)
        order = np.argsort(e)
        k = np.minimum(np.searchsorted(e[order], expnum), e.size - 1)
        return np.where(e[order][k] == expnum, order[k], -1)

    def corpus_coefficient(self, expnum):
        """ATOM_Or predicted for exposures of the build corpus (NaN elsewhere), as an input for training another
        model on that corpus: the test-night prediction for the rows this model was trained on."""
        if self.corpus is None or not self.corpus['pred_from_cv'].any():
            raise ValueError('this model holds no test-night predictions of its corpus: rebuild it without --no-cv')
        k = self.corpus_row(expnum)
        return np.where(k >= 0, 10 ** self.corpus['log10_pred'][np.maximum(k, 0)] / K['6300'], np.nan)

    def save(self, path):
        with open(path, 'wb') as fh:
            pickle.dump(self, fh)

    @staticmethod
    def load(path):
        with open(path, 'rb') as fh:
            m = pickle.load(fh)
        if m.inputs_version != INPUTS_VERSION:
            raise ValueError(f'{path} was built for inputs version {m.inputs_version}, this code builds version '
                             f'{INPUTS_VERSION}: rebuild the model (oi6300_build)')
        return m
