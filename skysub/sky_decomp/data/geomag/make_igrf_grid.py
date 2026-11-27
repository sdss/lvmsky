"""Build igrf_285km_lco.npz: IGRF declination and inclination on a grid around LCO at 285 km, and the
magnetic conjugate point (field line traced to 285 km in the northern hemisphere) of the LCO zenith pierce point.

Dev-only generator (needs ``ppigrf``, which is not a runtime dependency of lvmsky):

    python make_igrf_grid.py

The runtime code (``mlp_predictor/ionosphere.py``) only interpolates the stored arrays.
"""

import datetime as dt
from importlib.metadata import version
from pathlib import Path

import numpy as np
import ppigrf

HEIGHT_KM = 285.0
LAT = np.arange(-45.0, -12.99, 0.5)
LON = np.arange(-90.0, -49.99, 0.5)
YEARS = np.arange(2023.0, 2027.01, 1.0)
LCO_LAT, LCO_LON = -29.00889, -70.68992
R_EARTH_KM = 6371.2  # IGRF reference radius


def _unit_b(pos, date):
    """Unit field vector (dr, dtheta, dphi per km of arc length) in geocentric spherical coordinates."""
    r, theta, phi = pos
    br, bt, bp = (np.asarray(b).ravel()[0] for b in ppigrf.igrf_gc(r, theta, phi, date))
    b = np.sqrt(br * br + bt * bt + bp * bp)
    return np.array([br / b, np.degrees(bt / b / r), np.degrees(bp / b / (r * np.sin(np.radians(theta))))])


def conjugate_point(lat, lon, height_km, date, step_km=5.0):
    """Trace the IGRF field line from (lat, lon, height) along +B (upward in the south) back down to the same
    radius; RK4 in geocentric spherical coordinates, latitude treated as geocentric (error ~0.2 deg)."""
    r0 = R_EARTH_KM + height_km
    pos = np.array([r0, 90.0 - lat, lon])
    apex = r0
    for _ in range(100000):
        k1 = _unit_b(pos, date)
        k2 = _unit_b(pos + 0.5 * step_km * k1, date)
        k3 = _unit_b(pos + 0.5 * step_km * k2, date)
        k4 = _unit_b(pos + step_km * k3, date)
        new = pos + step_km / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        apex = max(apex, new[0])
        if apex > r0 + 1.0 and new[0] <= r0:
            f = (pos[0] - r0) / (pos[0] - new[0])
            end = pos + f * (new - pos)
            return 90.0 - end[1], end[2], apex - R_EARTH_KM
        pos = new
    raise RuntimeError('field line did not return to the start radius')


def main():
    glon, glat = np.meshgrid(LON, LAT)
    dec = np.empty((YEARS.size, LAT.size, LON.size), np.float32)
    inc = np.empty_like(dec)
    for k, year in enumerate(YEARS):
        be, bn, bu = ppigrf.igrf(glon, glat, HEIGHT_KM, dt.datetime(int(year), 1, 1))
        be, bn, bu = (np.asarray(b).reshape(glat.shape) for b in (be, bn, bu))
        dec[k] = np.degrees(np.arctan2(be, bn))
        inc[k] = np.degrees(np.arctan2(-bu, np.hypot(be, bn)))
    conj = np.array([conjugate_point(LCO_LAT, LCO_LON, HEIGHT_KM, dt.datetime(int(y), 1, 1)) for y in YEARS])
    out = Path(__file__).with_name('igrf_285km_lco.npz')
    np.savez(out, lat=LAT, lon=LON, year=YEARS, declination_deg=dec, inclination_deg=inc,
             conjugate_lat_deg=conj[:, 0], conjugate_lon_deg=conj[:, 1], apex_height_km=conj[:, 2],
             height_km=HEIGHT_KM, model=f"IGRF via ppigrf {version('ppigrf')}")
    print(out, dec.shape, 'conjugate (lat, lon, apex km) per year:', np.round(conj, 2).tolist())


if __name__ == '__main__':
    main()
