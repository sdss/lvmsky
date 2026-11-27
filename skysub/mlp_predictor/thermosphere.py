"""F-layer and thermosphere context for SkyML: what sets [O I] 6300 besides the electron column.

[O I] 6300 is made by O+ + O2 -> O2+ followed by O2+ + e -> O(1D), and O(1D) is quenched by N2, O2 and electrons,
so its zenith brightness is roughly n_e x [O2] at the layer, which falls steeply as the layer rises.  The IONEX
column (``ionosphere.py``) gives no height.  For each telescope's 285 km pierce point and exposure time:

- PyIRI (IRI-2020 F2 parameters, CCIR) gives the electron-density profile; it is scaled by one factor so that its
  vertical TEC (100-2000 km) equals the measured UQRG VTEC at the same point and time.  hmF2 and the bottomside
  thickness are IRI's (not measured); NmF2 and the profile amplitude are TEC-anchored.
- NRLMSIS 2.1 (pymsis) gives [O2], [N2], [O] and the neutral temperature, driven by the previous day's F10.7, its
  81-day mean and the 3-hour ap from the GFZ table already used by SkyML.
- HWM14 gives the neutral wind at 285 km; its magnetic-meridional part U drives the plasma along the field, which
  moves the layer vertically by w = -U sin(I) cos(I) (I = IGRF inclination; at LCO an equatorward wind lifts it).
  HWM14 is compiled from a pinned upstream commit by ``python -m mlp_predictor.thermosphere build-hwm`` (gfortran).
- The 6300 volume emission rate follows Link & Cogger (1988), with ion temperature = neutral temperature and
  [O+] = n_e:  eta = (A6300 / A1D) * beta * k1 [O2] n_e / (1 + (k2 [N2] + k3 [O2] + k4 n_e) / A1D),
  integrated over 150-700 km into Rayleigh; ``h6300`` is the emission-weighted height.  O I 7774 is radiative
  recombination, ~ integral of n_e^2 dz, given as a relative proxy.

Every quantity is vertical (zenith column or a local value at the layer), never along the slant path, so the van
Rhijn correction SkyML already applies to the airglow coefficients is not counted twice.

Run from ``skysub/``: ``python -m mlp_predictor.thermosphere build-hwm`` once, then
``python -m mlp_predictor.thermosphere table --iono-table PATH --out PATH`` (needs the ``thermosphere`` extra).
"""

import argparse
import datetime as dt
import os
from functools import lru_cache

import numpy as np

from mlp_predictor import ionosphere

ROLES = ('sci', 'near', 'far')
ALT_KM = np.concatenate([np.arange(100.0, 700.1, 5.0), np.arange(740.0, 2000.1, 40.0)])
EMIT = (ALT_KM >= 150.0) & (ALT_KM <= 700.0)
A6300, A1D, BETA = 5.63e-3, 7.45e-3, 1.1  # s^-1, s^-1, O(1D) yield of O2+ dissociative recombination
FEATURES = (
    'hmf2',
    'log10_nmf2',
    'b_bot',
    'log10_tec_scale',
    'log10_o2',
    'log10_n2',
    'log10_o',
    'tn',
    'log10_model6300',
    'h6300',
    'log10_model7774',
    'wind_mer_mag',
    'wind_w',
)
HWM_DIR = ionosphere.CACHE_DIR / 'hwm14'
HWM_COMMIT = '115fd958f4d196d6eb89bf6d997b5930355ada98'  # github.com/rilma/pyHWM14 (MIT)
HWM_FILES = {  # repository path -> SHA256
    'source/hwm14.f90': '7cd8bc52ccead92e29aac34980b8aaad26b238709384cb67c73b75038d7fd3bc',
    'data/hwm123114.bin': '6e445f8337c7efc815b7ff7f9f967d8a9a16b469d93df9c46e54930553beb906',
    'data/dwm07b104i.dat': 'f2b8eff002d55b0f6d49d202c73b7a9cb6685f2c6bb57f1291a2adc7aa07cf4f',
    'data/gd2qd.dat': '6bb1f2384e30b409240ee92c32726ed2a3d73eb550c05e0969d1178032319804',
}
# Standard Kp -> ap conversion (Kp in thirds: 0, 0+, 1-, 1, ...).
_KP = np.arange(0, 28) / 3.0
_AP = np.array(
    [0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400],
    float,
)


def kp_to_ap(kp):
    return np.interp(kp, _KP, _AP)


def ver6300(ne_cm3, o2, n2, tn):
    """6300 volume emission rate [photons cm^-3 s^-1] (Link & Cogger 1988; T_i = T_e = T_n)."""
    r = 300.0 / tn
    k1 = 3.23e-12 * np.exp(3.72 * r - 1.87 * r**2)
    k2 = 2.0e-11 * np.exp(107.8 / tn)
    k3 = 2.9e-11 * np.exp(67.5 / tn)
    k4 = 1.6e-12 * tn**0.91
    return (A6300 / A1D) * BETA * k1 * o2 * ne_cm3 / (1.0 + (k2 * n2 + k3 * o2 + k4 * ne_cm3) / A1D)


def model_versions():
    """Versions of the models behind the features, recorded in SkyML provenance."""
    from importlib.metadata import PackageNotFoundError, version

    def v(name):
        try:
            return version(name)
        except PackageNotFoundError:
            return None

    return {
        'iri': f"PyIRI {v('PyIRI')} (IRI-2020, CCIR)",
        'msis': f"pymsis {v('pymsis')} (NRLMSIS 2.1)",
        'hwm14_commit': HWM_COMMIT,
        'chemistry': 'Link & Cogger 1988',
        'ipp_height_km': ionosphere.IPP_HEIGHT_KM,
    }


def build_hwm():
    """Download the pinned HWM14 source and coefficients, check SHA256, compile with f2py into HWM_DIR / 'src'."""
    import hashlib
    import subprocess
    import sys
    import urllib.request

    src = HWM_DIR / 'src'
    src.mkdir(parents=True, exist_ok=True)
    for path, digest in HWM_FILES.items():
        dest = src / path.rsplit('/', 1)[1]
        if not dest.exists():
            url = f'https://raw.githubusercontent.com/rilma/pyHWM14/{HWM_COMMIT}/{path}'
            req = urllib.request.Request(url, headers={'User-Agent': 'lvmsky/1.0'})
            with urllib.request.urlopen(req, timeout=60) as r:
                tmp = dest.with_name(f'.{dest.name}.{os.getpid()}.part')
                tmp.write_bytes(r.read())
            tmp.replace(dest)
        if hashlib.sha256(dest.read_bytes()).hexdigest() != digest:
            raise ValueError(f'{dest}: SHA256 mismatch')
    # numpy 1.26 f2py builds with numpy.distutils, which needs the standard-library distutils (Python <= 3.11).
    cmd = [sys.executable, '-m', 'numpy.f2py', '-c', 'hwm14.f90', '-m', 'hwm14', '--f90flags=-w']
    subprocess.run(
        cmd, cwd=src, env={**os.environ, 'SETUPTOOLS_USE_DISTUTILS': 'stdlib'}, check=True, capture_output=True
    )


@lru_cache(maxsize=1)
def _hwm():
    import importlib
    import sys

    src = HWM_DIR / 'src'
    if not list(src.glob('hwm14*.so')):
        raise ImportError(f'HWM14 is not built in {src}; run python -m mlp_predictor.thermosphere build-hwm')
    os.environ['HWMPATH'] = str(src)
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    return importlib.import_module('hwm14')


def wind(mjd, lat, lon, f107_prev, f107a, ap, h_km=ionosphere.IPP_HEIGHT_KM):
    """HWM14 magnetic-meridional wind [m/s, toward magnetic north] and the vertical plasma drift it drives [m/s]."""
    hwm = _hwm()
    mer, zon = np.empty(np.size(mjd)), np.empty(np.size(mjd))
    for i, (m, la, lo) in enumerate(zip(mjd, lat, lon, strict=True)):
        sec = (m % 1.0) * 86400.0
        doy = (ionosphere._MJD0 + dt.timedelta(days=float(m))).timetuple().tm_yday
        stl = (sec / 3600.0 + lo / 15.0) % 24.0
        mer[i], zon[i] = hwm.hwm14(doy, sec, h_km, la, lo, stl, f107a[i], f107_prev[i], np.array([-1.0, ap[i]]))
    dec, inc, _ = ionosphere.geomag_at(lat, lon, mjd)
    u = mer * np.cos(np.radians(dec)) + zon * np.sin(np.radians(dec))
    return u, -u * np.sin(np.radians(inc)) * np.cos(np.radians(inc))


def _iri_day(day, ut_h, lat, lon, f107):
    """IRI profiles and F2 parameters at (ut_h[i], lat[i], lon[i]) for one UT day (grid evaluated, diagonal kept)."""
    import PyIRI
    import PyIRI.main_library as ml

    ut, iu = np.unique(ut_h, return_inverse=True)
    pts, ip = np.unique(np.column_stack([lon, lat]), axis=0, return_inverse=True)
    res = ml.IRI_density_1day(
        day.year, day.month, day.day, ut, pts[:, 0], pts[:, 1], ALT_KM, float(f107), PyIRI.coeff_dir, 0
    )
    f2, edp = res[0], res[-1]
    return (edp[iu, :, ip.ravel()], f2['hm'][iu, ip.ravel()], f2['Nm'][iu, ip.ravel()], f2['B_bot'][iu, ip.ravel()])


def features(mjd, lat, lon, vtec, f107_prev, f107a, ap):
    """All FEATURES for points (mjd, lat, lon) with measured VTEC [TECU]; arrays of one length."""
    from pymsis import msis

    mjd, lat, lon, vtec = (np.asarray(x, np.float64) for x in (mjd, lat, lon, vtec))
    n = mjd.size
    edp = np.empty((n, ALT_KM.size))
    hm, nm, bb = np.empty(n), np.empty(n), np.empty(n)
    day = np.floor(mjd).astype(int)
    for d in np.unique(day):
        k = np.flatnonzero(day == d)
        date = ionosphere._MJD0 + dt.timedelta(days=int(d))
        f = float(np.clip(np.nanmedian(f107_prev[k]), 60.0, 300.0))
        edp[k], hm[k], nm[k], bb[k] = _iri_day(date, (mjd[k] - d) * 24.0, lat[k], lon[k], f)
    scale = vtec / _tec(edp)
    ne = edp * scale[:, None]  # m^-3
    when = (np.datetime64('1858-11-17') + (mjd * 86400e3).astype('timedelta64[ms]')).astype('datetime64[ms]')
    alts = ALT_KM[EMIT]
    nz = alts.size
    rep = lambda x: np.repeat(np.asarray(x), nz)  # noqa: E731
    out = msis.calculate(
        rep(when),
        rep(lon),
        rep(lat),
        np.tile(alts, n),
        rep(f107_prev),
        rep(f107a),
        np.repeat(np.column_stack([ap] * 7), nz, axis=0),
    )
    n2, o2, o, tn = (out[:, i].reshape(n, nz) * (1e-6 if i < 10 else 1.0) for i in (1, 2, 3, 10))  # cm^-3, K
    ne_cm3 = ne[:, EMIT] * 1e-6
    v = ver6300(ne_cm3, o2, n2, tn)
    dz = np.gradient(alts) * 1e5  # cm
    i6300 = (v * dz).sum(axis=1) / 1e6  # Rayleigh
    h6300 = (v * alts * dz).sum(axis=1) / (v * dz).sum(axis=1)
    j = np.argmin(np.abs(alts[None, :] - hm[:, None]), axis=1)
    pick = lambda a: a[np.arange(n), j]  # noqa: E731
    return {
        'hmf2': hm,
        'log10_nmf2': np.log10(nm * scale),
        'b_bot': bb,
        'log10_tec_scale': np.log10(scale),
        'log10_o2': np.log10(pick(o2)),
        'log10_n2': np.log10(pick(n2)),
        'log10_o': np.log10(pick(o)),
        'tn': pick(tn),
        'log10_model6300': np.log10(i6300),
        'h6300': h6300,
        'log10_model7774': np.log10((ne_cm3**2 * dz).sum(axis=1) / 1e30),
        **dict(zip(('wind_mer_mag', 'wind_w'), wind(mjd, lat, lon, f107_prev, f107a, ap), strict=True)),
    }


def _tec(edp):
    """Vertical TEC [TECU] of profiles on ALT_KM (trapezoid, 100-2000 km)."""
    return np.trapz(edp, ALT_KM * 1e3, axis=1) / 1e16


TABLE_DIR = ionosphere.CONTEXT_TABLE_DIR / 'thermosphere'


@lru_cache(maxsize=1)
def bundled_table():
    """Bundled thermosphere tables, row-aligned with ``ionosphere._bundled_tables()`` (same sources, same order)."""
    tables = [dict(np.load(p, allow_pickle=False)) for p in sorted(TABLE_DIR.glob('*.npz'))]
    iono = ionosphere._bundled_tables()
    if not tables or iono is None:
        return None
    merged = {k: np.concatenate([t[k] for t in tables]) for k in tables[0] if np.ndim(tables[0][k]) == 1}
    order = np.argsort(merged['obstime_mjd'], kind='stable')
    merged = {k: v[order] for k, v in merged.items()}
    if not np.array_equal(merged['obstime_mjd'], iono['obstime_mjd']):
        raise ValueError('bundled thermosphere and ionosphere tables are not row-aligned; rebuild the former')
    return merged


def build_table(iono_table, out):
    """Per-role thermosphere features for every row of a bundled ionosphere context table."""
    from mlp_predictor import data

    t = np.load(iono_table)
    mjd = t['obstime_mjd'].astype(np.float64)
    f107_prev = data._solar_activity_lookup(mjd - 1.0, 'f107').astype(np.float64)
    f107a = data._solar_activity_lookup(mjd, 'f107_81d').astype(np.float64)
    ap = kp_to_ap(data._solar_activity_lookup(mjd, 'kp').astype(np.float64))
    table = {'obstime_mjd': mjd, 'f107_prev': f107_prev, 'f107a': f107a, 'ap': ap}
    for role in ROLES:
        f = features(mjd, t[f'{role}_ipp_lat'], t[f'{role}_ipp_lon'], t[f'{role}_vtec'], f107_prev, f107a, ap)
        table |= {f'{role}_{k}': np.asarray(v, np.float32) for k, v in f.items()}
        print(f'{role}: done', flush=True)
    np.savez_compressed(out, **table)
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=('build-hwm', 'table'))
    ap.add_argument('--iono-table')
    ap.add_argument('--out')
    args = ap.parse_args()
    if args.step == 'build-hwm':
        build_hwm()
    else:
        build_table(args.iono_table, args.out)


if __name__ == '__main__':
    main()
