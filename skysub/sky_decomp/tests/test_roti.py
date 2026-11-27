"""GNSS ROTI: RINEX parsing, broadcast orbit, elevation, ROTI from phase and the per-exposure lookup."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import roti  # noqa: E402

HEADER = (
    '     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE\n'
    '  1984093.3050 -5068867.9860 -3314488.0340                  APPROX POSITION XYZ\n'
    'G    4 C1C L1C L2X L2W                                      SYS / # / OBS TYPES\n'
    'R    2 C1C L1C                                              SYS / # / OBS TYPES\n'
    '                                                            END OF HEADER\n'
)


def _field(v, lli=' '):
    return f'{v:14.3f}{lli} '


def _rinex(epochs):
    out = [HEADER]
    for k, (l1, l2, lli) in enumerate(epochs):
        out.append(f'> 2024 10 22 00 {k // 2:02d} {30.0 * (k % 2):10.7f}  0  2\n')
        out.append('G05' + _field(2.0e7) + _field(l1, lli) + ' ' * 16 + _field(l2) + '\n')
        out.append('R07' + _field(2.0e7) + _field(1.0e8) + '\n')
    return ''.join(out)


def test_parse_obs_reads_gps_l1c_l2w_and_lli():
    obs = roti.parse_obs(_rinex([(1.0e8, 7.0e7, ' '), (1.0e8 + 5, 7.0e7 + 3, '1')]))
    assert obs['prn'].tolist() == [5, 5]
    assert obs['l1'].tolist() == [1.0e8, 1.0e8 + 5]
    assert obs['l2'].tolist() == [7.0e7, 7.0e7 + 3]  # L2W, not the empty L2X column
    assert obs['lli'].tolist() == [False, True]
    assert obs['mjd'][1] - obs['mjd'][0] == pytest.approx(30 / 86400)
    assert obs['rx'] == pytest.approx([1984093.305, -5068867.986, -3314488.034])


def test_tec_scale_matches_dual_frequency_definition():
    # 1 TECU delays L2 relative to L1 by 40.3e16 (1/f2^2 - 1/f1^2) metres, about 0.105 m
    assert 1 / roti.TECU_PER_M == pytest.approx(0.1050, abs=2e-4)


def _circular_eph(prn=5, a_m=26_560e3):
    keys = (
        'prn',
        'af0',
        'af1',
        'af2',
        'iode',
        'crs',
        'dn',
        'm0',
        'cuc',
        'e',
        'cus',
        'sqrta',
        'toe',
        'cic',
        'omega0',
        'cis',
        'i0',
        'crc',
        'w',
        'omegadot',
        'idot',
        'codes',
        'week',
        'l2p',
        'acc',
        'health',
    )
    eph = {k: np.zeros(1) for k in keys}
    eph.update(
        prn=np.array([prn]), sqrta=np.array([np.sqrt(a_m)]), i0=np.array([np.radians(55.0)]), toe=np.array([86400.0])
    )
    return eph


def test_broadcast_orbit_radius_and_validity_window():
    eph = _circular_eph()
    pos = roti.sat_ecef(eph, 5, np.array([86400.0, 90000.0, 86400.0 + 5 * 3600]))
    assert np.linalg.norm(pos[:2], axis=1) == pytest.approx([26_560e3, 26_560e3], rel=1e-9)
    assert np.isnan(pos[2]).all()  # more than 4 h from toe
    assert np.isnan(roti.sat_ecef(eph, 9, np.array([86400.0]))).all()  # no ephemeris


def test_elevation_of_zenith_and_horizon():
    rx = np.array([6378137.0, 0.0, 0.0])
    el, _ = roti.el_az(rx, np.array([[2.0e7 + 6378137.0, 0.0, 0.0], [6378137.0, 1.0e7, 0.0]]))
    assert el == pytest.approx([90.0, 0.0], abs=1e-6)


def test_roti_is_std_of_rot_and_slips_start_new_arcs(monkeypatch):
    t = 60605 + np.arange(20) * 30 / 86400  # ten minutes of 30 s epochs, two 5-min bins
    rot = np.where(np.arange(19) % 2 == 0, 1.0, -1.0)  # TECU/min alternating: std 1 in each bin
    stec = np.concatenate([[0.0], np.cumsum(rot * 0.5)])
    l1 = stec / roti.TECU_PER_M / (roti.C_LIGHT / roti.F1)  # all of the TEC on L1 phase
    obs = dict(
        rx=np.array([1984093.3, -5068868.0, -3314488.0]),
        mjd=t,
        prn=np.full(20, 5),
        l1=l1,
        l2=np.zeros(20),
        lli=np.zeros(20, bool),
    )
    monkeypatch.setattr(roti, 'sat_ecef', lambda eph, prn, sow: np.array([[2.6e7, -1.0e7, -1.5e7]]))
    monkeypatch.setattr(roti, 'el_az', lambda rx, sat: (np.array([60.0]), np.array([0.0])))
    s = roti.roti_samples(obs, eph=None)
    assert s['roti'].size >= 1
    assert s['roti'] == pytest.approx(np.ones(s['roti'].size), abs=0.06)
    obs['lli'][5:15] = True  # loss of lock kills most of both bins
    assert roti.roti_samples(obs, eph=None)['roti'].size == 0


def test_features_at_radius_and_window():
    smp = {
        'mjd': np.array([0.0, 0.0, 0.02]),
        'ipp_lat': np.array([-29.0, -35.0, -29.0]),
        'ipp_lon': np.array([-70.0, -70.0, -70.0]),
        'roti': np.array([0.5, 2.0, 3.0]),
        'el': np.zeros(3),
    }
    f = roti.features_at(smp, np.array([0.0, 0.5]), np.array([-29.5, -29.5]), np.array([-70.0, -70.0]))
    assert f[0, :3].tolist() == [0.5, 0.5, 1.0]  # 6 deg away and 29 min later are excluded
    assert np.isnan(f[1, 0]) and f[1, 2] == 0


def test_exposure_flag_tristate_and_never_raises(monkeypatch):
    days = {}

    def fake_station_day(sta, day, download=True):
        if sta == 'BAD':
            raise OSError('server down')
        return days.get((sta, day), {**roti.roti_samples_empty(), 'missing': np.array(True)})

    monkeypatch.setattr(roti, 'station_day', fake_station_day)
    ipp = {'Sci': (-29.5, -70.5), 'SkyE': (-26.0, -71.0)}
    f = roti.exposure_flag(60605.05, ipp, stations=('BAD', 'OK'))
    assert f['bubble'] == -1 and f['stations'] == [] and 'BAD' in f['errors'][0]
    days[('OK', 60605)] = {
        'mjd': np.array([60605.05, 60605.05]),
        'prn': np.zeros(2),
        'el': np.full(2, 60.0),
        'ipp_lat': np.array([-29.4, -10.0]),
        'ipp_lon': np.array([-70.5, -70.0]),
        'roti': np.array([0.04, 3.0]),
        'missing': np.array(False),
    }
    f = roti.exposure_flag(60605.05, ipp, stations=('BAD', 'OK'))
    assert f['bubble'] == 0 and f['stations'] == ['OK']  # the 3.0 sample is 2000 km away
    assert f['telescopes']['Sci']['n_samples'] == 1 and f['telescopes']['SkyE']['n_samples'] == 0
    days[('OK', 60605)]['roti'][0] = 1.5
    assert roti.exposure_flag(60605.05, ipp, stations=('OK',))['bubble'] == 1


def test_highpass_removes_trend_and_keeps_a_depletion():
    t = np.arange(0, 3 * 3600, 30.0)
    bubble = -8.0 * np.exp(-0.5 * ((t - 5400) / 600.0) ** 2)  # 8 TECU dip, ~20 min wide
    d = roti._arc_highpass(t, 40.0 + 2e-3 * t + bubble, np.zeros(t.size, int))
    centre = np.argmin(np.abs(t - 5400))
    assert d[centre] == pytest.approx(-8.0 * 0.7, abs=1.5)  # most of the dip survives the 60 min window
    far = np.abs(t - 5400) > 3000
    assert np.nanmax(np.abs(d[far])) < 0.5  # the linear trend is gone away from it
    assert np.isnan(d[:30]).all() and np.isfinite(d[60])  # arc ends carry no estimate
    short = roti._arc_highpass(t[:60], np.ones(60), np.zeros(60, int))
    assert np.isnan(short).all()  # arcs shorter than the window give no estimate


def test_mapping_is_one_at_zenith_and_grows_to_the_horizon():
    assert roti.mapping(90.0) == pytest.approx(1.0)
    assert 1.7 < roti.mapping(30.0) < 1.9


def test_aligned_distance_favours_the_field_direction():
    north = roti._aligned_dist_km(-29.0, -70.0, -29.0 + 600 / 111.6, -70.0)
    east = roti._aligned_dist_km(-29.0, -70.0, -29.0, -70.0 + 600 / (111.6 * np.cos(np.radians(29.0))))
    assert north == pytest.approx(600 / roti.MER_SCALE * (6656 / 6656), rel=0.05)
    assert east == pytest.approx(600, rel=0.05)


def test_recent_404_is_not_cached(tmp_path, monkeypatch):
    import urllib.error

    def not_found(urls, **fmt):
        raise urllib.error.HTTPError('u', 404, 'Not Found', None, None)

    monkeypatch.setattr(roti, 'ROTI_DIR', tmp_path)
    monkeypatch.setattr(roti, '_get_first', not_found)
    today = int(roti.time.time() / 86400.0 + 40587.0)
    with pytest.raises(urllib.error.HTTPError):
        roti.station_day('SANT00CHL', today - 1)  # may still appear on the archive: not remembered as missing
    assert not list(tmp_path.glob('*.npz'))
    assert bool(roti.station_day('SANT00CHL', today - 30)['missing'])  # old enough: cached as missing


def test_roti_at_groups_by_night_and_tolerates_missing_days(monkeypatch):
    def station_day(sta, day, download=False):
        if sta != 'OAFA00ARG' or day != 60000:
            raise FileNotFoundError(day)
        return {'mjd': np.array([60000.75]), 'prn': np.zeros(1), 'el': np.full(1, 60.0),
                'ipp_lat': np.array([-29.0]), 'ipp_lon': np.array([-70.7]), 'roti': np.array([0.8]),
                'dvtec': np.full(1, np.nan)}

    monkeypatch.setattr(roti, 'station_day', station_day)
    f = roti.roti_at([60000.75, 60010.75], [-29.0, -29.0], [-70.7, -70.7], aligned=True)
    assert f[0].tolist() == [0.8, 0.8, 1.0]
    assert np.isnan(f[1, 0]) and f[1, 2] == 0  # no station-day at all that night


def test_bundled_roti_table_is_row_aligned_with_ionosphere():
    tab = roti.bundled_table()
    if tab is None:
        pytest.skip('no bundled tables')
    from mlp_predictor import ionosphere

    assert np.array_equal(tab['obstime_mjd'], ionosphere._bundled_tables()['obstime_mjd'])
    assert np.mean(tab['sci_roti_n'] > 0) > 0.85


def test_roti_augment_columns_and_missing_rows(monkeypatch):
    from mlp_predictor import data

    monkeypatch.setattr(roti, 'bundled_table', lambda: None)
    monkeypatch.setattr(roti, 'roti_at', lambda mjd, lat, lon, **kw: np.array([[0.4, 0.1, 5.0], [np.nan, np.nan, 0.0]]))
    names = ['alt', 'az_sin', 'az_cos']
    ctx = np.array([[60.0, 0.0, 1.0], [70.0, 1.0, 0.0]], np.float32)
    tri = {'ctx_names': names, 'obstime_mjd': np.array([1.0, 2.0]), 'ctx_sci': ctx, 'ctx_near': ctx, 'ctx_far': ctx}
    data._augment_triplet_with_roti(tri)
    assert tri['ctx_names'][-3:] == list(data.ROTI_NAMES)
    np.testing.assert_allclose(tri['ctx_sci'][:, 3:], [[np.log10(0.4), -1.0, 1.0], [np.log10(data.ROTI_QUIET)] * 2 + [0.0]],
                               rtol=1e-6)
    data._augment_triplet_with_roti(tri)  # force=True replaces, does not duplicate
    assert len(tri['ctx_names']) == 6


def test_dvtec_event_takes_the_nearest_high_ray():
    smp = {'mjd': np.array([0.0, 0.0, 0.0, 0.01]), 'el': np.array([60.0, 60.0, 20.0, 60.0]),
           'ipp_lat': np.array([-29.0, -30.0, -29.0, -29.0]), 'ipp_lon': np.array([-70.5, -70.0, -70.0, -70.0]),
           'dvtec': np.array([-1.0, -2.0, -3.0, -4.0])}
    f = roti.dvtec_event(smp, np.array([0.0, 0.5]), np.array([-29.0, -29.0]), np.array([-70.0, -70.0]))
    assert f[0, 0] == -2.0  # 111 km north counts as ~28 km along the field; the low ray and the late one are skipped
    assert np.isnan(f[1, 0])


def test_navigation_download_skips_an_incomplete_file(monkeypatch):
    import gzip
    import urllib.error

    nav = lambda n: gzip.compress(('END OF HEADER\n' + ''.join(f'G{p:02d} 2023 09 01\n' for p in range(1, n + 1))).encode())  # noqa: E731
    files = {'BRDC00WRD_S': nav(9), 'BRDC00IGS_R': nav(31)}  # BRDC00WRD_R is absent that day

    def get(url, timeout=180):
        for name, raw in files.items():
            if name in url:
                return raw
        raise urllib.error.HTTPError(url, 404, 'not found', None, None)

    monkeypatch.setattr(roti, '_get', get)
    assert roti._get_nav(2023, 244) == files['BRDC00IGS_R']
    del files['BRDC00IGS_R']
    assert roti._get_nav(2023, 244) == files['BRDC00WRD_S']  # nothing complete: the most complete file
    files.clear()
    with pytest.raises(urllib.error.HTTPError):
        roti._get_nav(2023, 244)
