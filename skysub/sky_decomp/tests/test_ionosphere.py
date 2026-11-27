"""Ionospheric context: IONEX parsing, pierce-point geometry, and the single-count van Rhijn contract."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import data  # noqa: E402
from mlp_predictor import ionosphere as io

LAT = (-25.0, -27.5, -30.0, -32.5)
LON = (-80.0, -75.0, -70.0, -65.0, -60.0)
DAY = 60409  # 2024-04-09


def _label(text, label):
    return f'{text:<60s}{label}'


def _write_ionex(path, tec_maps, epochs_hms, exponent=-1):
    """Minimal IONEX 1.0 file; ``tec_maps[k]`` holds raw integers, shape (len(LAT), len(LON))."""
    out = [
        _label('     1.0            IONOSPHERE MAPS     GPS', 'IONEX VERSION / TYPE'),
        _label('  COSZ', 'MAPPING FUNCTION'),
        _label('   450.0 450.0   0.0', 'HGT1 / HGT2 / DHGT'),
        _label(f'  {LAT[0]:6.1f}{LAT[-1]:6.1f}{LAT[1] - LAT[0]:6.1f}', 'LAT1 / LAT2 / DLAT'),
        _label(f'  {LON[0]:6.1f}{LON[-1]:6.1f}{LON[1] - LON[0]:6.1f}', 'LON1 / LON2 / DLON'),
        _label(f'{exponent:6d}', 'EXPONENT'),
        _label('', 'END OF HEADER'),
    ]
    for kind in ('TEC', 'RMS'):
        for k, (grid, (h, m)) in enumerate(zip(tec_maps, epochs_hms, strict=True), start=1):
            out.append(_label(f'{k:6d}', f'START OF {kind} MAP'))
            out.append(_label(f'  2024     4     9{h:6d}{m:6d}     0', 'EPOCH OF CURRENT MAP'))
            for i, lat in enumerate(LAT):
                out.append(
                    _label(
                        f'  {lat:6.1f}{LON[0]:6.1f}{LON[-1]:6.1f}{LON[1] - LON[0]:6.1f} 450.0', 'LAT/LON1/LON2/DLON/H'
                    )
                )
                vals = np.asarray(grid)[i] if kind == 'TEC' else np.full(len(LON), 30)
                out.append(''.join(f'{int(v):5d}' for v in vals))
            out.append(_label(f'{k:6d}', f'END OF {kind} MAP'))
    path.write_text('\n'.join(out) + '\n')


@pytest.fixture
def uniform_cache(tmp_path, monkeypatch):
    """Cache with one day of spatially uniform 25.0 TECU maps at 00:00 and 24:00 (hour 24, as IONEX writes)."""
    grid = np.full((len(LAT), len(LON)), 250)
    (tmp_path / 'UQRG').mkdir()
    _write_ionex(tmp_path / 'UQRG' / 'uqrg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(io, 'SOURCES', {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i'})
    io._day_maps.cache_clear()
    yield tmp_path
    io._day_maps.cache_clear()


def test_read_ionex_scales_and_masks(tmp_path):
    grid = np.arange(len(LAT) * len(LON)).reshape(len(LAT), len(LON)) * 10
    grid[0, 0] = 9999
    _write_ionex(tmp_path / 'f.24i', [grid, grid], [(0, 0), (24, 0)])
    m = io.read_ionex(tmp_path / 'f.24i')
    assert m['tec'].shape == (2, len(LAT), len(LON))
    assert np.isnan(m['tec'][0, 0, 0])
    assert m['tec'][0, 1, 2] == pytest.approx(grid[1, 2] * 0.1)
    np.testing.assert_allclose(m['epoch'], [DAY, DAY + 1.0])
    np.testing.assert_allclose(m['lat'], LAT)
    np.testing.assert_allclose(m['lon'], LON)


def test_pierce_point_geometry():
    lat, lon, dn, de = io.pierce_point(90.0, 0.0)
    assert (lat, lon, dn, de) == pytest.approx((io.LCO_LAT_DEG, io.LCO_LON_DEG, 0.0, 0.0), abs=1e-9)
    # Due north: latitude moves by the Earth-central angle, longitude stays.
    el = np.radians(40.0)
    r = io.EARTH_RADIUS_KM / (io.EARTH_RADIUS_KM + io.IPP_HEIGHT_KM)
    psi = np.pi / 2 - el - np.arcsin(r * np.cos(el))
    lat, lon, dn, de = io.pierce_point(40.0, 0.0)
    assert lat - io.LCO_LAT_DEG == pytest.approx(np.degrees(psi))
    assert lon == pytest.approx(io.LCO_LON_DEG)
    assert dn == pytest.approx((io.EARTH_RADIUS_KM + io.IPP_HEIGHT_KM) * psi)
    assert de == pytest.approx(0.0, abs=1e-9)


def test_ionex_mapping_function_is_the_van_rhijn_factor():
    """IONEX's single-layer mapping function at h is the same expression as data._van_rhijn_factor."""
    alt = np.array([90.0, 70.0, 45.0, 30.0])
    z = np.radians(90.0 - alt)
    r = io.EARTH_RADIUS_KM / (io.EARTH_RADIUS_KM + io.IPP_HEIGHT_KM)
    mapping = 1.0 / np.sqrt(1.0 - (r * np.sin(z)) ** 2)
    np.testing.assert_allclose(data._van_rhijn_factor(alt, io.IPP_HEIGHT_KM), mapping, rtol=1e-6)


def test_vtec_features_carry_no_van_rhijn_factor(uniform_cache):
    """A uniform vertical-TEC sky must give the same VTEC for every elevation, while V(z) changes.

    data.airglow_geometry_scale already divides airglow coefficients by V(z, 285 km); a slant
    (V-weighted) TEC here would apply the geometry a second time.
    """
    alts = np.array([90.0, 60.0, 40.0, 30.0])
    mjd = np.full(alts.size, DAY + 0.2)
    feats, source = io.telescope_features(mjd, alts, np.full(alts.size, 90.0))
    np.testing.assert_allclose(feats['vtec'], 25.0)
    for k in ('dvtec_dt', 'dlog10vtec_dt', 'dlog10vtec_dlat', 'dlog10vtec_dlon'):
        np.testing.assert_allclose(feats[k], 0.0, atol=1e-12)
    assert set(source) == {'UQRG'}
    v = data._van_rhijn_factor(alts, io.IPP_HEIGHT_KM)
    assert v[-1] / v[0] > 1.5


def test_vtec_time_interpolation_rotates_with_the_sun(tmp_path, monkeypatch):
    """Two maps 24 h apart with the same longitude gradient: the rotated interpolation is exact at mid-day."""
    grid = np.tile(np.arange(len(LON)) * 10 + 200, (len(LAT), 1))
    (tmp_path / 'UQRG').mkdir()
    _write_ionex(tmp_path / 'UQRG' / 'uqrg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(io, 'SOURCES', {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i'})
    io._day_maps.cache_clear()
    # t = 0.01 d: lon shifts by +3.6 deg for map 0 and -356.4 = +3.6 deg (mod 360) for map 1.
    vals, _, _ = io.vtec_at([DAY + 0.01], [-27.5], [-75.0])
    assert vals[0] == pytest.approx(21.0 + 0.1 * 3.6 / 5.0 * 10, rel=1e-9)
    io._day_maps.cache_clear()


def test_truncated_file_falls_back_to_next_source(tmp_path, monkeypatch):
    """UPC serves a truncated UQRG for some days; the next source is used and reported per row."""
    (tmp_path / 'UQRG').mkdir()
    (tmp_path / 'UPRG').mkdir()
    (tmp_path / 'UQRG' / 'uqrg1000.24i').write_text(_label('IGS GPS stations used', 'COMMENT') + '\n')
    grid = np.full((len(LAT), len(LON)), 300)
    _write_ionex(tmp_path / 'UPRG' / 'uprg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(
        io,
        'SOURCES',
        {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i', 'UPRG': 'file:///nonexistent/uprg{doy:03d}0.{yy}i'},
    )
    io._day_maps.cache_clear()
    vals, _, source = io.vtec_at([DAY + 0.3], [-29.0], [-70.7])
    io._day_maps.cache_clear()
    assert vals[0] == pytest.approx(30.0)
    assert list(source) == ['UPRG']


def test_missing_ionex_raises(uniform_cache):
    with pytest.raises(FileNotFoundError, match='No IONEX map'):
        io.vtec_at([DAY + 3.2], [-29.0], [-70.7])


def test_local_vtec_gradient_at_pierce_point(tmp_path, monkeypatch):
    """A map growing linearly in raw TEC along latitude gives the analytic d log10 VTEC / d lat."""
    grid = np.tile((200 + 20 * np.arange(len(LAT)))[:, None], (1, len(LON)))  # +2 TECU per 2.5 deg southward
    (tmp_path / 'UQRG').mkdir()
    _write_ionex(tmp_path / 'UQRG' / 'uqrg1000.24i', [grid, grid], [(0, 0), (24, 0)])
    monkeypatch.setattr(io, 'CACHE_DIR', tmp_path)
    monkeypatch.setattr(io, 'SOURCES', {'UQRG': 'file:///nonexistent/uqrg{doy:03d}0.{yy}i'})
    io._day_maps.cache_clear()
    feats, _ = io.telescope_features([DAY + 0.2], [90.0], [0.0])
    io._day_maps.cache_clear()
    lat = io.LCO_LAT_DEG
    tec = lambda x: 20.0 + 2.0 * (LAT[0] - x) / 2.5  # noqa: E731
    expected = (np.log10(tec(lat + 1.25)) - np.log10(tec(lat - 1.25))) / 2.5
    assert feats['dlog10vtec_dlat'][0] == pytest.approx(expected, rel=1e-9)
    assert feats['dlog10vtec_dlon'][0] == pytest.approx(0.0, abs=1e-12)


def test_conjugate_point_and_sun_altitude():
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, get_sun
    from astropy.time import Time

    g = np.load(io.GEOMAG_GRID)
    assert np.all((g['conjugate_lat_deg'] > 5.0) & (g['conjugate_lat_deg'] < 12.0))  # northern EIA side of LCO
    assert np.all((g['apex_height_km'] > 800.0) & (g['apex_height_km'] < 1200.0))
    mjd = 61000.0  # 2025-11-21, inside the grid years
    year = 2000.0 + (mjd - 51544.5) / 365.25
    lat, lon = (np.interp(year, g['year'], g[k]) for k in ('conjugate_lat_deg', 'conjugate_lon_deg'))
    t = Time(mjd, format='mjd')
    ref = get_sun(t).transform_to(AltAz(obstime=t, location=EarthLocation(lat=lat * u.deg, lon=lon * u.deg)))
    assert io.sun_alt_conjugate([mjd])[0] == pytest.approx(ref.alt.deg, abs=1e-6)


def test_geomag_grid_at_lco():
    dec, inc, dip = io.geomag_at(io.LCO_LAT_DEG, io.LCO_LON_DEG, 60500.0)
    assert -4.0 < dec[0] < 0.0
    assert -33.0 < inc[0] < -28.0
    assert dip[0] == pytest.approx(np.degrees(np.arctan(np.tan(np.radians(inc[0])) / 2)))


def test_exponent_inside_a_map_overrides_the_header(tmp_path):
    grid = np.full((len(LAT), len(LON)), 300)
    _write_ionex(tmp_path / 'f.24i', [grid, grid], [(0, 0), (24, 0)], exponent=-1)
    text = (
        (tmp_path / 'f.24i')
        .read_text()
        .replace(
            _label('     1', 'START OF TEC MAP') + '\n',
            _label('     1', 'START OF TEC MAP') + '\n' + _label('    -2', 'EXPONENT') + '\n',
            1,
        )
    )
    (tmp_path / 'f.24i').write_text(text)
    m = io.read_ionex(tmp_path / 'f.24i')
    assert m['tec'][0, 0, 0] == pytest.approx(3.0)  # first map: exponent -2 from inside the map
    assert m['tec'][1, 0, 0] == pytest.approx(30.0)  # second map: header exponent -1


def test_nan_rms_is_an_error_not_a_feature(uniform_cache, monkeypatch):
    maps = io._day_maps(DAY)

    def fake(day):
        return maps | {'rms': np.full_like(maps['tec'], np.nan)}

    fake.cache_clear = lambda: None
    monkeypatch.setattr(io, '_day_maps', fake)
    with pytest.raises(ValueError, match='RMS samples are not finite'):
        io.vtec_at(np.array([DAY + 0.1]), np.array([-28.0]), np.array([-70.0]))
