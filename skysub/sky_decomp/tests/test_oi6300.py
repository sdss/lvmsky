"""Science [O I] 6300 model: base blend, training on a synthetic corpus, quality flag, save/load."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import oi6300  # noqa: E402


def synthetic(n=600, seed=0):
    rng = np.random.default_rng(seed)
    tab = {'sun_alt_conjugate': rng.uniform(-60, 0, n), 'f107': rng.uniform(70, 200, n),
           'f107a': rng.uniform(70, 200, n), 'kp': rng.uniform(0, 5, n), 'night_mjd': np.arange(n) // 10}
    zen = 10 ** rng.uniform(1.3, 2.5, n)  # zenith 6300 [R], common to the three telescopes
    for a in oi6300.ARMS:
        tab[f'{a}_alt'] = rng.uniform(35, 85, n)
        tab[f'{a}_vanrhijn_285km'] = 1 / np.sqrt(1 - (6371 / 6656) ** 2 * np.cos(np.radians(tab[f'{a}_alt'])) ** 2)
        for c in oi6300.IONO + oi6300.GEOM[:3]:
            tab[f'{a}_{c}'] = rng.uniform(0.5, 2.0, n)
        tab[f'{a}_ctx_moon_alt'] = rng.uniform(-60, 60, n)
        for line in oi6300.K:
            tab[f'{a}_coef_{line}'] = zen * tab[f'{a}_vanrhijn_285km'] / oi6300.K[line]
    tab['sci_coef_6300'] = tab['sci_coef_6300'] * 10 ** (0.02 * (tab['sci_wind_w'] - 1.25))  # learnable correction
    for a, b in (('sci', 'near'), ('sci', 'far'), ('near', 'far')):
        tab[f'sep_{a}_{b}'] = rng.uniform(1, 30, n)
    return tab


def test_base_is_van_rhijn_scaled_skynear():
    tab = synthetic()
    m = oi6300.OI6300Model(ctx_names=('moon_alt',))
    m.alpha = 1.0
    expect = np.log10(oi6300.K['6300'] * tab['near_coef_6300'] / tab['near_vanrhijn_285km'] * tab['sci_vanrhijn_285km'])
    assert m.base(tab) == pytest.approx(expect)


def test_fit_predict_flag_and_roundtrip(tmp_path):
    pytest.importorskip('sklearn')
    tab = synthetic()
    rows = np.arange(400)
    m = oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, rows, domain=np.arange(400, 600))
    y = np.log10(oi6300.K['6300'] * tab['sci_coef_6300'])
    pred, flag = m.predict(tab)
    blend_only = m.base(tab)
    test = np.arange(400, 600)
    assert np.std((y - pred)[test]) < 0.5 * np.std((y - blend_only)[test])  # the correction is learned
    assert flag.shape == pred.shape and flag.min() >= 0
    m.save(tmp_path / 'm.pkl')
    assert oi6300.OI6300Model.load(tmp_path / 'm.pkl').predict(tab)[0] == pytest.approx(pred)
    assert m.coefficient(tab) == pytest.approx(10 ** pred / oi6300.K['6300'])


def test_ensemble_averages_members():
    pytest.importorskip('sklearn')
    tab = synthetic()
    rows = np.arange(400)
    ens = oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, rows, members=3)
    assert len(ens.models) == 3
    Z = ens._matrix(tab)[1][:, ens.keep]
    expect = ens.base(tab) + np.mean([m.predict(Z) for m in ens.models], axis=0)
    assert ens.predict(tab)[0] == pytest.approx(expect)


def with_science_lines(tab, seed=1):
    """Add a telescope factor common to all Science lines, the control lines and [O I] 5577 to the synthetic corpus."""
    rng = np.random.default_rng(seed)
    n = tab['kp'].size
    common = rng.normal(0, 0.03, n)  # dex, Science relative to SkyNear
    for a in oi6300.ARMS:
        tab[f'{a}_vanrhijn_95km'] = np.ones(n)
        tab[f'{a}_coef_5577'] = np.full(n, 100.0)
    for line in list(oi6300.K) + ['5577']:
        own = 0 if line == '6300' else rng.normal(0, 0.003, n)
        tab[f'sci_coef_{line}'] = tab[f'sci_coef_{line}'] * 10 ** (common + own)
    oh = np.tile(common[:, None], (1, 5)) + rng.normal(0, 0.003, (n, 5))
    tab['sci_near_log_oh'] = oi6300.control_ratio(10 ** oh, np.ones((n, 5)), np.ones(n), np.ones(n))
    tab['sci_near_log_o2'] = common + rng.normal(0, 0.003, n)
    return tab


def test_control_ratio_is_median_zenith_log_ratio():
    c_sci, c_near = np.array([[2.0, 20.0, 400.0], [1.0, -1.0, 3.0]]), np.array([[1.0, 10.0, 100.0], [1.0, 1.0, 1.0]])
    got = oi6300.control_ratio(c_sci, c_near, np.array([2.0, 1.0]), np.array([1.0, 1.0]))
    assert got == pytest.approx([np.log10(2.0) - np.log10(2.0), np.log10(3.0) / 2])  # non-positive coefficient skipped


def test_science_lines_help_and_guard_drops_contaminated_line():
    pytest.importorskip('sklearn')
    tab = with_science_lines(synthetic())
    rows, test = np.arange(400), np.arange(400, 600)
    y = np.log10(oi6300.K['6300'] * tab['sci_coef_6300'])
    plain = oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, rows)
    m = oi6300.OI6300Model(ctx_names=('moon_alt',), sci_lines=('oh', 'o2', '5577', '7774')).fit(tab, rows)
    assert m.guarded == ('5577', '7774') and set(m.fallback) == {(), ('5577',), ('7774',)}
    pred = m.predict(tab)[0]
    assert np.std((y - pred)[test]) < 0.6 * np.std((y - plain.predict(tab)[0])[test])  # the common factor is learned
    assert np.mean(m.line_levels(tab)['7774'][test] > m.guard) < 0.1
    dirty = dict(tab, sci_coef_7774=tab['sci_coef_7774'] * 10.0)  # a nebula adds 1 dex to Science 7774
    assert np.all(m.line_levels(dirty)['7774'][test] > m.guard)
    keep = ~(m.line_levels(tab)['5577'] > m.guard)
    Z, S = m._matrix(tab)[1], oi6300.science_lines(tab, m.sci_lines)
    expect = m.base(tab) + m.fallback[('5577',)][0].predict(m._x(Z, S, ('5577',)))
    assert m.predict(dirty)[0][keep] == pytest.approx(expect[keep])  # the model without 7774, unaffected by it
    missing = dict(tab, sci_coef_7774=np.where(np.arange(600) < 500, tab['sci_coef_7774'], np.nan))
    assert m.predict(missing)[0][keep & (np.arange(600) >= 500)] == pytest.approx(
        expect[keep & (np.arange(600) >= 500)])  # a missing line is handled like a dropped one
    assert m.predict(missing)[0][:500] == pytest.approx(pred[:500])
    unguarded = oi6300.OI6300Model(ctx_names=('moon_alt',), sci_lines=m.sci_lines, guard=None).fit(tab, rows)
    assert not unguarded.fallback and unguarded.predict(tab)[0][test] == pytest.approx(m.models[0].predict(
        m._x(Z, S, m.guarded))[test] + m.base(tab)[test])
    with pytest.raises(ValueError):
        oi6300.OI6300Model(sci_lines=('6300',))


def test_build_row_sets_folds_and_cross_validation():
    pytest.importorskip('sklearn')
    from mlp_predictor import oi6300_build as ob

    tab = synthetic()
    i = np.arange(tab['kp'].size)
    aux = {'n_fibres': {a: np.full(i.size, 20) for a in oi6300.ARMS},
           'log_ew': {a: np.full(i.size, 2.0) for a in oi6300.ARMS},
           'valid_coords': i != 0, 'lmc': i < 200, 'smc': np.zeros(i.size, bool), 'gal_b': np.full(i.size, 40.0),
           'halpha_excess': np.where(i % 2, 0.0, 5.0)}  # [R]; even rows have gas
    sets = ob.row_sets(tab, aux)
    assert not sets['clean'][0] and sets['lmcsmc'].sum() == 199  # the row with broken coordinates is dropped
    assert sets['gaspoor'].sum() == 200 and not np.any(sets['gaspoor'] & sets['lmcsmc'])
    assert np.array_equal(sets['train'], sets['lmcsmc'] | sets['gaspoor'])
    fold = ob.night_folds(tab)
    assert set(fold) == set(range(ob.N_FOLDS))
    assert all(np.unique(fold[tab['night_mjd'] == x]).size == 1 for x in np.unique(tab['night_mjd']))  # whole nights
    pred, dropped = ob.cross_validate(tab, sets['train'], sets['domain'], ctx_names=('moon_alt',))
    assert np.all(np.isfinite(pred[sets['train']])) and np.all(np.isnan(pred[~sets['train']])) and not dropped.any()
    seen = []  # the inputs learned from the corpus are rebuilt from the other nights of each run
    ob.cross_validate(tab, sets['train'], sets['domain'], prepare=seen.append, ctx_names=('moon_alt',))
    assert len(seen) == ob.N_FOLDS and all(not np.any(s & (fold == k)) and np.all(s | (fold == k))
                                           for k, s in enumerate(seen))
    s = ob.scores(tab, pred, sets)
    assert s['all']['rows'] == 399 and s['LMC+SMC']['rows'] == 199
    base = oi6300.OI6300Model().base(tab, alpha=1.0)  # SkyNear rescaled, no correction
    y = np.log10(oi6300.K['6300'] * tab['sci_coef_6300'])
    assert s['all']['sigma'] < 0.7 * ob.sigma((y - base)[sets['train']])


def test_fit_rejects_rows_without_a_valid_line():
    pytest.importorskip('sklearn')
    tab = synthetic()
    tab['near_coef_6300'][3] = 0.0
    with pytest.raises(ValueError):
        oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, np.arange(400))


def test_invalid_options_are_rejected_before_the_fit():
    pytest.importorskip('sklearn')
    tab = with_science_lines(synthetic())
    for kw in ({'sci_lines': ('7774', '7774')}, {'guard': 0.0}):
        with pytest.raises(ValueError):
            oi6300.OI6300Model(**kw)
    with pytest.raises(ValueError, match='members'):
        oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, np.arange(400), members=0)
    with pytest.raises(ValueError, match='nights'):  # one night cannot hold the night groups of the guard
        oi6300.OI6300Model(ctx_names=('moon_alt',), sci_lines=('oh', '7774')).fit(tab, np.arange(10))


def test_missing_inputs_stay_missing_and_a_zero_weight_telescope_is_not_read(tmp_path):
    pytest.importorskip('sklearn')
    tab = synthetic()
    tab['near_roti_max'][:5] = [np.nan, 0.0, 0.001, 0.005, 0.2]  # missing, below the floor, at it, ordinary
    got = oi6300._arm(tab, 'near')['roti_max'][:5]
    assert np.isnan(got[0]) and got[1:] == pytest.approx(np.log10([0.005, 0.005, 0.005, 0.2]))
    m = oi6300.OI6300Model(ctx_names=('moon_alt',)).fit(tab, np.arange(5, 400))
    assert np.isnan(m._matrix(tab)[1][0, m.names.index('near roti_max')])
    pred, flag = m.predict(tab)
    no_far = dict(tab, far_coef_6300=np.where(np.arange(600) == 500, np.nan, tab['far_coef_6300']))
    m.alpha = 1.0  # SkyNear alone: a missing SkyFar line must not make the base NaN
    assert np.isfinite(m.predict(no_far)[0][500])
    q = m.quality(no_far)
    assert q['missing'][500] > m.quality(tab)['missing'][500]  # the SkyFar line is counted as missing
    assert np.array_equal(m.quality(tab)['out_of_range'], flag)
    m.alpha = 0.5
    assert np.isnan(m.predict(no_far)[0][500])
    no_near = dict(tab, near_coef_6300=np.where(np.arange(600) == 500, 0.0, tab['near_coef_6300']))
    assert np.isnan(m.predict(no_near)[0][500])  # no valid sky line: no prediction
    m.inputs_version = oi6300.INPUTS_VERSION - 1  # a model built for other inputs is refused
    m.save(tmp_path / 'old.pkl')
    with pytest.raises(ValueError, match='rebuild'):
        oi6300.OI6300Model.load(tmp_path / 'old.pkl')


def night_table(n_nights=40, per_night=12, seed=1):
    """Synthetic corpus whose layer has a meridional gradient that differs from night to night."""
    tab = synthetic(n_nights * per_night, seed)
    rng, n = np.random.default_rng(seed), n_nights * per_night
    tab['night_mjd'] = 60000 + np.arange(n) // per_night
    tab['obstime_mjd'] = tab['night_mjd'] + (2.0 + 0.5 * (np.arange(n) % per_night)) / 24.0
    tab['near_is_skye'] = rng.random(n) < 0.5
    slope = rng.normal(0, 2e-4, n_nights)[np.arange(n) // per_night]  # [dex / km]
    for a in oi6300.ARMS:
        tab[f'{a}_ipp_dx_mag_mer'], tab[f'{a}_ipp_dx_mag_zon'] = rng.uniform(-300, 300, n), rng.uniform(-300, 300, n)
        zen = 50.0 * 10 ** (slope * tab[f'{a}_ipp_dx_mag_mer'])  # [R]; the same layer for the three telescopes
        tab[f'{a}_coef_6300'] = zen * tab[f'{a}_vanrhijn_285km'] / oi6300.K['6300']
    return tab, slope


def test_night_field_recovers_a_gradient_and_improves_the_model():
    pytest.importorskip('sklearn')
    from mlp_predictor import sky_field

    tab, slope = night_table()
    n = tab['kp'].size
    usable = np.ones(n, bool)
    usable[:11] = False  # the first night keeps one exposure: two samples, no field
    hyper = sky_field.fit(tab, usable, maxiter=150)  # a short fit; the parameters are learned from the sky arms
    assert set(hyper) == {*sky_field.FRAMES, 'fit'} and hyper['fit']['nights'] == 39  # not the 1-exposure night
    assert all(np.isfinite(list(hyper[f].values())).all() for f in sky_field.FRAMES)
    assert not any(hyper['fit'][f]['converged'] for f in sky_field.FRAMES)  # 150 steps are recorded as too few
    with pytest.raises(ValueError):
        sky_field.fit(tab, np.zeros(n, bool))  # no usable sample: nothing to fit
    f = sky_field.night_field(tab, usable, hyper)
    assert set(f) == set(oi6300.SKY_FIELD) and np.all(np.isnan(f['field_sci_near_earth'][:12]))
    truth = slope * (tab['sci_ipp_dx_mag_mer'] - tab['near_ipp_dx_mag_mer'])
    assert np.corrcoef(f['field_sci_near_earth'][12:], truth[12:])[0, 1] > 0.9
    assert f['field_n_1h'][12:].min() >= 4 and f['field_n_1h'][12:].max() <= 10  # +-1 h holds up to 5 exposures
    from mlp_predictor import oi6300_build as ob

    aux = {'n_fibres': {a: np.full(n, 20) for a in oi6300.ARMS}, 'valid_coords': np.arange(n) != 20}
    broken = dict(tab, far_coef_6300=np.where(np.arange(n) == 20, 1e3, tab['far_coef_6300']),
                  far_ipp_dx_mag_mer=np.where(np.arange(n) == 20, 5e3, tab['far_ipp_dx_mag_mer']))
    fields = []
    for t in (tab, broken):  # an exposure with a broken pointing gives no sample to its night, whatever it holds
        t = dict(t)
        ob.night_field(t, aux, hyper)
        fields.append(t['field_sci_near_earth'])
    assert fields[0][12:] == pytest.approx(fields[1][12:], abs=1e-12)
    tab |= f
    rows, test = np.arange(n) < 360, np.arange(n) >= 360
    y = np.log10(oi6300.K['6300'] * tab['sci_coef_6300'])
    err = {k: np.abs(y - oi6300.OI6300Model(ctx_names=('moon_alt',), sky_field=k).fit(tab, rows).predict(tab)[0])[test]
           for k in (False, True)}
    assert np.median(err[True]) < 0.5 * np.median(err[False])
