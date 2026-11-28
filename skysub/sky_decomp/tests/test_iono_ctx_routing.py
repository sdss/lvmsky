"""Ionospheric ctx columns: the augment appends vertical-only features; iono_ctx_groups routes them to F-region
branches only; checkpoints saved without the option still load."""

import sys
from pathlib import Path

import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import data  # noqa: E402
from mlp_predictor import ionosphere as io  # noqa: E402
from mlp_predictor.model import DualEncoderGroupHeadMLPCompressed  # noqa: E402

BASE = ['alt', 'az_sin', 'az_cos', 'airmass', 'moon_up_smooth', 'ecl_beta_deg', 'zodi_log10_v']
GROUPS = {'moon': 2, 'zodi': 2, 'continuum': 2, 'mesospheric': 3, 'iono_red': 1, 'iono_rec': 2, 'oi_green': 1}


def _model(ctx_names, iono_ctx_groups=None, seed=0):
    torch.manual_seed(seed)
    return DualEncoderGroupHeadMLPCompressed(
        n_score=sum(GROUPS.values()),
        n_ctx=len(ctx_names),
        group_score_dims=GROUPS,
        ctx_names=ctx_names,
        encoder_dims=(16,),
        ctx_dims=(8,),
        trunk_dims=(16,),
        head_dim=8,
        continuum_branch_dims=(8, 8),
        zodi_ctx_restriction=('zodi_log10_v', 'ecl_beta_deg'),
        continuum_ctx_restriction=('alt', 'airmass'),
        moon_zodi_ctx_restriction=('moon_up_smooth', 'airmass'),
        moon_zodi_coupling_dims=(8, 8),
        isolated_branches={'iono_red': (8,), 'iono_rec': (8,), 'oi_green': (8,)},
        isolated_branch_inputs=('iono_red', 'iono_rec', 'oi_green'),
        iono_ctx_groups=iono_ctx_groups,
    ).eval()


def _inputs(n_ctx, rng):
    t = lambda *s: torch.as_tensor(rng.normal(size=s), dtype=torch.float32)  # noqa: E731
    n_score = sum(GROUPS.values())
    return [t(4, n_score), t(4, n_score), t(4, n_ctx), t(4, n_ctx), t(4, n_ctx)]


def test_iono_columns_reach_only_the_listed_branches():
    names = BASE + list(data.IONOSPHERE_FEATURE_NAMES)
    model = _model(names, iono_ctx_groups=('iono_red', 'iono_rec'))
    rng = np.random.default_rng(1)
    x = _inputs(len(names), rng)
    with torch.no_grad():
        ref = model(*x)
        for k in (2, 3, 4):  # perturb only the iono_* columns of every arm
            x[k][:, len(BASE) :] += 5.0
        new = model(*x)
    changed = {g for g in GROUPS if not torch.allclose(ref[g], new[g])}
    assert changed == {'iono_red', 'iono_rec'}


def test_without_routing_every_path_sees_iono_columns():
    names = BASE + list(data.IONOSPHERE_FEATURE_NAMES)
    model = _model(names)
    x = _inputs(len(names), np.random.default_rng(2))
    with torch.no_grad():
        ref = model(*x)
        for k in (2, 3, 4):
            x[k][:, len(BASE) :] += 5.0
        new = model(*x)
    assert not torch.allclose(ref['mesospheric'], new['mesospheric'])
    assert not torch.allclose(ref['oi_green'], new['oi_green'])


def test_old_checkpoint_loads_strictly():
    old = _model(BASE)
    new = _model(BASE, seed=3)
    new.load_state_dict(old.state_dict(), strict=True)  # shared_ctx_idx is not persistent
    x = _inputs(len(BASE), np.random.default_rng(4))
    with torch.no_grad():
        for g, v in old(*x).items():
            assert torch.equal(v, new(*x)[g])


def test_routing_rejects_missing_iono_columns():
    with pytest.raises(ValueError, match='no iono_'):
        _model(BASE, iono_ctx_groups=('iono_red',))


def test_ionosphere_augment_appends_vertical_features(tmp_path, monkeypatch):
    from test_ionosphere import DAY, LAT, LON, _write_ionex

    grid = np.full((len(LAT), len(LON)), 250)
    (tmp_path / 'UQRG').mkdir()
    _write_ionex(tmp_path / 'UQRG' / 'uqrg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(io, 'SOURCES', {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i'})
    io._day_maps.cache_clear()
    alt = np.array([[80.0], [50.0], [30.0]])
    az = np.radians([[10.0], [120.0], [250.0]])
    ctx = lambda a, b: np.hstack([a, np.sin(b), np.cos(b)]).astype(np.float32)  # noqa: E731
    triplet = {
        'ctx_names': ['alt', 'az_sin', 'az_cos'],
        'obstime_mjd': np.array([DAY + 0.2]),
        'ctx_sci': ctx(alt[0:1], az[0:1]),
        'ctx_near': ctx(alt[1:2], az[1:2]),
        'ctx_far': ctx(alt[2:3], az[2:3]),
    }
    data._augment_triplet_with_ionosphere(triplet)
    data._augment_triplet_with_ionosphere(triplet)  # force=True strips and recomputes, no duplicates
    io._day_maps.cache_clear()
    names = triplet['ctx_names']
    assert names == ['alt', 'az_sin', 'az_cos', *data.IONOSPHERE_FEATURE_NAMES]
    j = names.index('iono_log10_vtec')
    for arm in ('ctx_sci', 'ctx_near', 'ctx_far'):  # uniform 25 TECU map: identical for every elevation
        assert triplet[arm][0, j] == pytest.approx(np.log10(25.0), abs=1e-6)
    assert list(triplet['ionex_source']) == ['UQRG']


def test_augment_reads_bundled_table_first(tmp_path, monkeypatch):
    """A row matching the bundled table (obstime and every arm's alt/az) needs no IONEX; others are computed."""
    from test_ionosphere import DAY, LAT, LON, _write_ionex

    grid = np.full((len(LAT), len(LON)), 250)
    (tmp_path / 'UQRG').mkdir()
    _write_ionex(tmp_path / 'UQRG' / 'uqrg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(io, 'SOURCES', {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i'})
    alt, az = {'sci': 80.0, 'near': 50.0, 'far': 30.0}, {'sci': 10.0, 'near': 120.0, 'far': 250.0}
    table = {'obstime_mjd': np.array([DAY + 0.2]), 'sun_alt_conjugate': np.array([-7.0], np.float32)}
    for r in io.ROLES:
        table |= {
            f'{r}_alt': np.array([alt[r]], np.float32),
            f'{r}_az': np.array([az[r]], np.float32),
            f'{r}_ionex_source': np.array(['TEST']),
        }
        table |= {f'{r}_{f}': np.array([1000.0 if f == 'vtec' else 0.5], np.float32) for f in io._TABLE_FEATURES}
    tables = tmp_path / 'tables'
    tables.mkdir()
    np.savez(tables / 'corpus.npz', **table)
    monkeypatch.setattr(io, 'CONTEXT_TABLE_DIR', tables)
    io._bundled_tables.cache_clear()
    io._day_maps.cache_clear()

    def ctx(r, da=0.0):
        return np.array([[alt[r] + da, np.sin(np.radians(az[r])), np.cos(np.radians(az[r]))]], np.float32)

    def make(da):
        return {
            'ctx_names': ['alt', 'az_sin', 'az_cos'],
            'obstime_mjd': np.array([DAY + 0.2]),
            **{f'ctx_{r}': ctx(r, da) for r in io.ROLES},
        }

    hit, miss = make(0.0), make(0.5)  # a 0.5 deg pointing change is not the tabled exposure
    data._augment_triplet_with_ionosphere(hit)
    data._augment_triplet_with_ionosphere(miss)
    io._bundled_tables.cache_clear()
    io._day_maps.cache_clear()
    j = hit['ctx_names'].index('iono_log10_vtec')
    assert hit['ctx_sci'][0, j] == pytest.approx(3.0) and list(hit['ionex_source']) == ['TEST']
    assert hit['ctx_sci'][0, -1] == pytest.approx(-7.0)
    assert miss['ctx_sci'][0, j] == pytest.approx(np.log10(25.0), abs=1e-6) and list(miss['ionex_source']) == ['UQRG']


def test_inference_explains_missing_ionospheric_inputs(monkeypatch):
    from mlp_predictor import inference

    def missing(triplet, force=True):
        raise FileNotFoundError('no IONEX source covers MJD 61400')

    monkeypatch.setattr(data, '_augment_triplet_with_ionosphere', missing)
    with pytest.raises(RuntimeError, match='IONEX maps appear 1-2 days after the night'):
        inference.build_triplet_from_pointings(
            obstime_mjd=[60605.05],
            sci_ra=[40.0],
            sci_dec=[-30.0],
            sky_near_ra=[41.0],
            sky_near_dec=[-30.0],
            sky_far_ra=[60.0],
            sky_far_dec=[-40.0],
            ionosphere=True,
            verbose=False,
        )


def test_thermosphere_augment_records_provenance(monkeypatch):
    from mlp_predictor import thermosphere

    n = 3
    tri = {'ctx_names': BASE + list(data.IONOSPHERE_FEATURE_NAMES), 'obstime_mjd': np.full(n, 60605.05)}
    for arm in ('ctx_near', 'ctx_far', 'ctx_sci'):
        tri[arm] = np.ones((n, len(tri['ctx_names'])), np.float32)
    monkeypatch.setattr(io, 'context_from_tables', lambda *a, **k: (np.zeros(n, bool), np.zeros(0, int)))
    monkeypatch.setattr(
        thermosphere, 'features', lambda mjd, *a: {k: np.full(np.size(mjd), 1.0) for k in thermosphere.FEATURES}
    )
    monkeypatch.setattr(data, '_solar_activity_lookup', lambda mjd, f: np.full(np.size(mjd), 150.0))
    data._augment_triplet_with_thermosphere(tri, force=True)
    assert list(tri['ctx_names'][-6:]) == list(data.THERMOSPHERE_FEATURE_NAMES)
    assert tri['thermosphere_source'].tolist() == ['computed'] * n
    assert tri['thermosphere_models']['hwm14_commit'] == thermosphere.HWM_COMMIT


def test_sky_line_features_use_only_the_sky_arms():
    n = 2
    names = ['alt', 'az_sin', 'az_cos', 'vanrhijn_285km', 'iono_log10_vtec']
    tri = {'ctx_names': list(names)}
    for arm, vtec in (('ctx_near', 1.0), ('ctx_far', 1.2), ('ctx_sci', 1.1)):
        tri[arm] = np.tile(np.array([45.0, 0.0, 1.0, 2.0, vtec], np.float32), (n, 1))
    coef_names = ['ATOM_Or', 'ATOM_Orc_OI0777']
    near = np.array([[200.0, 2.0], [200.0, 2.0]])  # zenith-equivalent 100 and 1
    far = np.array([[200.0, 8.0], [400.0, 2.0]])
    data._augment_triplet_with_sky_lines(tri, near, far, coef_names)
    i_a, i_q = tri['ctx_names'].index('iono_sky_altimeter'), tri['ctx_names'].index('iono_sky_q')
    assert tri['ctx_near'][:, i_a] == pytest.approx([2.0, 2.0])  # log 100 - 0.5 log 1
    assert tri['ctx_far'][:, i_a] == pytest.approx([2.0 - 0.5 * np.log10(4), np.log10(200)])
    assert tri['ctx_sci'][:, i_a] == pytest.approx(tri['ctx_near'][:, i_a])  # never Science's own lines
    # near/far 7774 ratio beyond 2 x the map's log VTEC difference (-0.2): row 0 |log(1/4) + 0.4|, row 1 |0 + 0.4|
    assert tri['ctx_far'][:, i_q] == pytest.approx([abs(np.log10(0.25) + 0.4), 0.4], abs=1e-6)
    assert np.array_equal(tri['ctx_sci'][:, i_q], tri['ctx_near'][:, i_q])


def test_oi6300_model_feature_is_the_square_root_with_a_sky_fallback():
    tri = {'ctx_names': ['alt']}
    for arm in ('ctx_near', 'ctx_far', 'ctx_sci'):
        tri[arm] = np.full((3, 1), 45.0, np.float32)
    near = np.array([[1.0, 9.0], [1.0, 16.0], [1.0, 25.0]])  # columns: another line, ATOM_Or
    data._augment_triplet_with_oi6300(tri, [4.0, np.nan, -1.0], near, ['OH_x', 'ATOM_Or'])
    assert tri['ctx_names'] == ['alt', *data.OI6300_MODEL_NAMES]
    for arm in ('ctx_near', 'ctx_far', 'ctx_sci'):  # the model value, else SkyNear's own line; flag 1 / 0
        assert tri[arm][:, 1:].tolist() == [[2.0, 1.0], [4.0, 0.0], [5.0, 0.0]]
    data._augment_triplet_with_oi6300(tri, [1.0, 1.0, 1.0], near, ['OH_x', 'ATOM_Or'])  # replaces, does not repeat
    assert tri['ctx_sci'].shape == (3, 3) and tri['ctx_sci'][:, 1:].tolist() == [[1.0, 1.0]] * 3
