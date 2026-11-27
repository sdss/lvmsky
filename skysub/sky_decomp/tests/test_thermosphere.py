"""F-layer / thermosphere context: emission chemistry, wind geometry, TEC anchoring and table alignment."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import thermosphere as th  # noqa: E402


def test_kp_to_ap_standard_table():
    assert th.kp_to_ap(np.array([0.0, 1.0, 3.0, 9.0])) == pytest.approx([0, 4, 15, 400])


def test_6300_emission_scales_with_o2_and_is_quenched():
    tn, ne = 1000.0, 1e6
    low = th.ver6300(ne, 1e7, 1e8, tn)
    assert th.ver6300(ne, 1e5, 1e6, tn) / th.ver6300(ne, 1e6, 1e7, tn) == pytest.approx(0.1, rel=0.05)  # thin air
    assert th.ver6300(ne, 1e7, 1e8, tn) / th.ver6300(ne, 1e6, 1e7, tn) < 9.0  # quenching bites low down
    assert low > 0
    assert th.ver6300(2e4, 1e6, 1e7, tn) / th.ver6300(1e4, 1e6, 1e7, tn) == pytest.approx(2, rel=0.01)  # ~ n_e
    assert th.ver6300(2e6, 1e6, 1e7, tn) / th.ver6300(1e6, 1e6, 1e7, tn) < 1.9  # electron quenching at the peak


def test_iri_profile_is_scaled_to_the_measured_tec():
    pytest.importorskip('PyIRI')
    pytest.importorskip('pymsis')
    try:
        th._hwm()
    except ImportError:
        pytest.skip('HWM14 not built')
    mjd = np.array([60605.05, 60605.05])
    f = th.features(
        mjd,
        np.array([-29.0, -29.0]),
        np.array([-70.7, -70.7]),
        np.array([30.0, 60.0]),
        np.full(2, 180.0),
        np.full(2, 170.0),
        np.full(2, 7.0),
    )
    assert f['log10_nmf2'][1] - f['log10_nmf2'][0] == pytest.approx(np.log10(2), abs=1e-6)  # amplitude anchored
    assert f['hmf2'][0] == f['hmf2'][1]  # height from IRI, not from the TEC
    assert f['log10_model6300'][1] - f['log10_model6300'][0] == pytest.approx(np.log10(2), abs=0.02)  # ~ n_e
    assert f['log10_model7774'][1] - f['log10_model7774'][0] == pytest.approx(2 * np.log10(2), abs=1e-6)  # n_e^2
    assert 200 < f['h6300'][0] < f['hmf2'][0]  # 6300 comes from the bottomside


def test_equatorward_wind_lifts_the_southern_layer(monkeypatch):
    class _Hwm:
        @staticmethod
        def hwm14(*args):
            return np.array([50.0, 0.0])  # 50 m/s toward geographic north

    monkeypatch.setattr(th, '_hwm', lambda: _Hwm)
    u, w = th.wind(
        np.array([60605.05]),
        np.array([-29.0]),
        np.array([-70.7]),
        np.array([180.0]),
        np.array([170.0]),
        np.array([7.0]),
    )
    assert u[0] == pytest.approx(50.0, rel=0.01)  # declination ~-2 deg
    assert 15.0 < w[0] < 25.0  # upward: 50 * sin(30.5) cos(30.5) ~ 22 m/s at dip -30.5 deg


def test_bundled_table_is_row_aligned_with_ionosphere():
    tab = th.bundled_table()
    if tab is None:
        pytest.skip('no bundled tables')
    from mlp_predictor import ionosphere

    assert np.array_equal(tab['obstime_mjd'], ionosphere._bundled_tables()['obstime_mjd'])
    assert all(np.isfinite(tab[f'near_{f}']).all() for f in th.FEATURES)
