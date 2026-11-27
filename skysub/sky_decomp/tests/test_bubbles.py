"""Plasma-bubble context: apex geometry, local time, IBP evaluation, IBI parsing and the reach fit."""

import io as _io
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import bubbles  # noqa: E402


def test_apex_height_geometry():
    assert bubbles.apex_height_km(0.0, 285.0) == pytest.approx(285.0)
    assert bubbles.apex_height_km(-16.0, 285.0) == pytest.approx(832.0, abs=2.0)  # LCO pierce point
    assert bubbles.apex_height_km(-16.0, 285.0) == pytest.approx(bubbles.apex_height_km(16.0, 285.0))


def test_local_time_folds_around_midnight():
    lt = bubbles.local_time_h(np.array([60000.0 + 1 / 24, 60000.5]), np.array([-60.0, 0.0]))
    assert lt == pytest.approx([-3.0, 12.0])


def test_ibp_grouped_matches_single_calls():
    ibpmodel = pytest.importorskip('ibpmodel')
    mjd = np.array([60604.95, 60605.1, 60605.2, 60605.1])
    lon = np.array([-70.0, -70.0, -70.0, -72.4])
    p = bubbles.ibp_probability(mjd, lon, 200.0)
    for i in range(mjd.size):
        doy = (bubbles.ionosphere._MJD0 + bubbles.ionosphere.dt.timedelta(days=mjd[i])).timetuple().tm_yday
        one = ibpmodel.calculateIBPindex(
            day_month=doy,
            longitude=[int(round(lon[i]))],
            local_time=[float(bubbles.local_time_h(mjd[i], lon[i]))],
            f107=[200.0],
        )
        assert p[i] == pytest.approx(float(one['IBP'].iloc[0]))


def test_fetch_keeps_valid_sector_samples(tmp_path, monkeypatch):
    csv = (
        '2024-10-22T01:00:00.000Z,-20.0,-60.0,6884000.0,1\n'
        '2024-10-22T01:00:01.000Z,-20.1,-60.0,6884000.0,-1\n'  # not applicable
        '2024-10-22T01:00:02.000Z,-20.2,-30.0,6884000.0,0\n'  # outside the sector
        '2024-10-22T01:00:03.000Z,-20.3,-60.0,6884000.0,0\n'
    )

    class _Resp(_io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(bubbles, 'IBI_DIR', tmp_path)
    monkeypatch.setattr(bubbles.urllib.request, 'urlopen', lambda req, timeout: _Resp(csv.encode()))
    s = bubbles.fetch_ibi_night(60604, 'B')
    assert s['ibi'].tolist() == [1, 0]
    assert s['h_km'] == pytest.approx([513.0, 513.0])
    assert s['mjd'][0] == pytest.approx(60605 + 1 / 24)
    assert (tmp_path / 'ibi_B_60604.npz').exists()
    with pytest.raises(FileNotFoundError):
        bubbles.fetch_ibi_night(60605, 'B', download=False)


def test_fit_reach_recovers_apex_slope(monkeypatch):
    rng = np.random.default_rng(0)
    n = 400
    monkeypatch.setattr(bubbles, 'ibp_probability', lambda mjd, lon, f: np.full(np.shape(mjd), 0.5))
    top = rng.uniform(400.0, 1400.0, n)  # bubble top apex; flagged on every field line below it
    passes = {
        'night': np.zeros(n),
        'mjd': np.zeros(n),
        'lon': np.full(n, -70.0),
        'lt': np.zeros(n),
        'apex_cov_km': np.full(n, 3000.0),
        'apex_bubble_km': np.where(top > 500, top, np.nan),
        'sat': np.zeros(n),
    }
    calib, (x, y, _) = bubbles.fit_reach(passes, {0: 150.0})
    c = np.array(calib['coef'])
    assert c[2] < 0  # higher apex is reached less often
    p = 1 / (1 + np.exp(-x @ c))
    a = bubbles.APEX_GRID_KM
    emp = np.array([(top >= ai).mean() for ai in a])
    model = np.array([p[np.isclose(x[:, 2], (ai - bubbles.APEX_REF_KM) / bubbles.APEX_SCALE_KM)].mean() for ai in a])
    assert np.max(np.abs(model - emp)) < 0.12
