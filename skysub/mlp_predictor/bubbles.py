"""Plasma-bubble context for SkyML: chance that an equatorial plasma bubble reaches each telescope's pierce point.

Bubbles rise over the magnetic equator and fill whole magnetic flux tubes, so a bubble whose top reaches apex
height A covers every field line with apex <= A, in both hemispheres.  A pierce point at height h and dip
latitude d sits on the field line with apex (R + h) / cos^2(d) - R (dipole), about 830 km for LCO at 285 km.

Prior: the ESA Swarm IBP climatological model (``ibpmodel``; Stolle et al. 2024, Space Weather 22,
e2023SW003809) gives the probability of a post-sunset bubble at a longitude, local time, day of year and
F10.7, but has no latitude dependence.  We calibrate the latitude reach on Swarm in-situ bubble-index (IBI)
passes over the LCO sector (VirES HAPI, public):

    P(bubble reaches apex A) = sigmoid(c0 + c1 * logit(IBP) + c2 * (A - 830 km) / 100 km + c3 * (F10.7 - 150) / 50)

The feature is a probability per arm; it carries no line-of-sight factor, so van Rhijn stays counted once.

Run from ``skysub/``: ``python -m mlp_predictor.bubbles fetch --oi-table PATH`` downloads and reduces the IBI
passes of the corpus nights; ``fit --oi-table PATH`` fits and writes ``CALIBRATION``.  Fetching and fitting
need the ``geomag`` extra (IGRF for Swarm sample positions); the per-exposure probability needs ``bubbles``.
"""

import argparse
import http.client
import io as _io
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from mlp_predictor import ionosphere

HAPI = (
    'https://vires.services/hapi/data?dataset=SW_OPER_IBI{sat}TMS_2F&parameters=Latitude,Longitude,Radius,Bubble_Index'
    '&start={start}&stop={stop}&format=csv'
)
SATS = ('A', 'B')  # Swarm C flies next to A
SECTOR_LON = (-90.0, -50.0)
SECTOR_LAT = (-45.0, 25.0)
IBI_DIR = ionosphere.CACHE_DIR / 'swarm_ibi'
CALIBRATION = ionosphere.CONTEXT_TABLE_DIR / 'bubble_reach_calibration.json'
APEX_REF_KM, APEX_SCALE_KM = 830.0, 100.0
APEX_GRID_KM = np.arange(500.0, 1601.0, 50.0)


def apex_height_km(dip_lat_deg, h_km):
    """Dipole apex height of the field line through height ``h_km`` at dip latitude ``dip_lat_deg``."""
    r = ionosphere.EARTH_RADIUS_KM
    return (r + h_km) / np.cos(np.radians(dip_lat_deg)) ** 2 - r


def local_time_h(mjd, lon_deg):
    """Solar local time in hours, folded to the IBP model range [-6, 18) around midnight."""
    lt = ((np.asarray(mjd) % 1.0) * 24.0 + np.asarray(lon_deg) / 15.0) % 24.0
    return np.where(lt >= 18.0, lt - 24.0, lt)


def _iso(mjd):
    return (ionosphere._MJD0 + ionosphere.dt.timedelta(days=float(mjd))).strftime('%Y-%m-%dT%H:%M:%SZ')


def fetch_ibi_night(night, sat, download=True):
    """Valid IBI samples of one Swarm satellite over the LCO sector, 18-06 LT of evening MJD ``night``."""
    path = IBI_DIR / f'ibi_{sat}_{int(night)}.npz'
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        url = HAPI.format(sat=sat, start=_iso(night + 22 / 24), stop=_iso(night + 1 + 10 / 24))
        req = urllib.request.Request(url, headers={'User-Agent': 'lvmsky/1.0'})
        with urllib.request.urlopen(req, timeout=300) as r:
            raw = r.read().decode()
        t = np.genfromtxt(_io.StringIO(raw), delimiter=',', usecols=0, dtype='U24')
        v = np.genfromtxt(_io.StringIO(raw), delimiter=',', usecols=(1, 2, 3, 4)).reshape(-1, 4)
        keep = (
            (v[:, 3] >= 0)
            & (v[:, 1] >= SECTOR_LON[0])
            & (v[:, 1] <= SECTOR_LON[1])
            & (v[:, 0] >= SECTOR_LAT[0])
            & (v[:, 0] <= SECTOR_LAT[1])
        )
        mjd = (np.char.rstrip(t[keep], 'Z').astype('datetime64[ms]') - np.datetime64('1858-11-17')).astype(
            float
        ) / 86400e3
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.part.npz')
        np.savez_compressed(
            tmp,
            mjd=mjd,
            lat=v[keep, 0],
            lon=v[keep, 1],
            h_km=v[keep, 2] / 1e3 - ionosphere.EARTH_RADIUS_KM,
            ibi=v[keep, 3].astype(np.int8),
        )
        tmp.replace(path)
    return dict(np.load(path))


def _dip_lat(lat, lon, h_km, mjd):
    """IGRF dip latitude atan(tan I / 2) at the given points (needs the ``geomag`` extra)."""
    import ppigrf

    date = ionosphere._MJD0 + ionosphere.dt.timedelta(days=float(np.median(mjd)))
    be, bn, bu = (np.asarray(b).ravel() for b in ppigrf.igrf(lon, lat, h_km, date))
    inc = np.arctan2(-bu, np.hypot(be, bn))
    return np.degrees(np.arctan(np.tan(inc) / 2.0))


def ibi_passes(nights, sats=SATS):
    """One row per sector pass: time, longitude, local time, apex range covered, highest apex flagged as bubble."""
    rows = []
    for night in nights:
        for sat in sats:
            s = fetch_ibi_night(night, sat)
            if s['mjd'].size == 0:
                continue
            cut = np.flatnonzero(np.diff(s['mjd']) * 86400 > 60) + 1
            for k in np.split(np.arange(s['mjd'].size), cut):
                if k.size < 60:
                    continue
                apex = apex_height_km(_dip_lat(s['lat'][k], s['lon'][k], s['h_km'][k], s['mjd'][k]), s['h_km'][k])
                bub = s['ibi'][k] == 1
                t, lon = float(np.median(s['mjd'][k])), float(np.median(s['lon'][k]))
                rows.append(
                    (
                        night,
                        t,
                        lon,
                        float(local_time_h(t, lon)),
                        apex.max(),
                        apex[bub].max() if bub.any() else np.nan,
                        SATS.index(sat),
                    )
                )
    names = ('night', 'mjd', 'lon', 'lt', 'apex_cov_km', 'apex_bubble_km', 'sat')
    a = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: a[:, i] for i, n in enumerate(names)}


def ibp_probability(mjd, lon_deg, f107):
    """Swarm IBP climatological bubble probability at the given times, longitudes and F10.7 (vectorised)."""
    import ibpmodel

    mjd, lon, f107 = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (mjd, lon_deg, f107)))
    doy = np.array(
        [(ionosphere._MJD0 + ionosphere.dt.timedelta(days=float(m))).timetuple().tm_yday for m in mjd.ravel()]
    )
    lt = local_time_h(mjd, lon).ravel()
    lon_i = np.clip(np.round(lon.ravel()).astype(int), -180, 179)
    f = np.clip(f107.ravel(), 80.0, 200.0)  # the IBP model is valid for 80-200 sfu
    doy = np.minimum(doy, 365)
    out = np.full(mjd.size, np.nan)
    for key in set(zip(doy, lon_i, f, strict=True)):  # ibpmodel evaluates a grid: one call per (day, longitude, F10.7)
        sel = np.flatnonzero((doy == key[0]) & (lon_i == key[1]) & (f == key[2]))
        lts = np.unique(lt[sel])
        df = ibpmodel.calculateIBPindex(
            day_month=int(key[0]), longitude=[int(key[1])], local_time=lts.tolist(), f107=[float(key[2])]
        )
        out[sel] = df['IBP'].to_numpy()[np.searchsorted(lts, lt[sel])]
    return out.reshape(mjd.shape)


def _design(ibp, apex_km, f107):
    p = np.clip(ibp, 1e-3, 1 - 1e-3)
    return np.column_stack(
        [
            np.ones_like(p),
            np.log(p / (1 - p)),
            (apex_km - APEX_REF_KM) / APEX_SCALE_KM,
            (np.asarray(f107) - 150.0) / 50.0,
        ]
    )


def reach_probability(mjd, lon_deg, apex_km, f107, coef=None):
    """Calibrated probability that a bubble reaches field-line apex ``apex_km`` at this time and longitude."""
    coef = np.asarray(coef if coef is not None else json.loads(CALIBRATION.read_text())['coef'])
    x = _design(
        ibp_probability(mjd, lon_deg, f107).ravel(), np.ravel(apex_km), np.ravel(np.broadcast_to(f107, np.shape(mjd)))
    )
    return (1.0 / (1.0 + np.exp(-x @ coef))).reshape(np.shape(mjd))


def fit_reach(passes, f107_by_night):
    """Logistic fit of 'bubble flagged at apex >= A' over passes x apex grid (A within each pass's coverage)."""
    from scipy.optimize import minimize

    f = np.array([f107_by_night[int(n)] for n in passes['night']])
    ibp = ibp_probability(passes['mjd'], passes['lon'], f)
    X, y, g = [], [], []
    for i in range(passes['mjd'].size):
        a = APEX_GRID_KM[APEX_GRID_KM <= passes['apex_cov_km'][i]]
        X.append(_design(np.full(a.size, ibp[i]), a, np.full(a.size, f[i])))
        y.append(np.nan_to_num(passes['apex_bubble_km'][i], nan=-1.0) >= a)
        g.append(np.full(a.size, passes['night'][i]))
    X, y, g = np.vstack(X), np.concatenate(y).astype(float), np.concatenate(g)

    def nll(c):
        z = X @ c
        return np.sum(np.logaddexp(0, z) - y * z) + 1e-3 * c @ c

    c = minimize(nll, np.zeros(X.shape[1]), method='BFGS').x
    return dict(
        coef=c.tolist(),
        n_passes=int(passes['mjd'].size),
        n_rows=int(y.size),
        n_nights=int(np.unique(g).size),
        bubble_passes=int(np.isfinite(passes['apex_bubble_km']).sum()),
        apex_ref_km=APEX_REF_KM,
        apex_scale_km=APEX_SCALE_KM,
        model='sigmoid(c0 + c1 logit(IBP) + c2 (A-830)/100 + c3 (F10.7-150)/50)',
    ), (X, y, g)


def exposure_table(oi_table, iono_table, coef=None):
    """Per exposure and telescope: IBP probability, pierce-point apex height and calibrated reach probability."""
    d, t = np.load(oi_table), np.load(iono_table)
    if not np.array_equal(d['expnum'], t['expnum']):
        raise ValueError('oi_table and ionosphere table rows differ')
    f107 = d['ctx_sci'][:, [str(n) for n in d['ctx_names']].index('f107')]
    coef = np.asarray(coef if coef is not None else json.loads(CALIBRATION.read_text())['coef'])
    out = {'expnum': t['expnum'], 'obstime_mjd': t['obstime_mjd']}
    for tel in ionosphere.TELESCOPES:
        apex = apex_height_km(t[f'{tel}_dip_lat'], float(t['ipp_height_km']))
        ibp = ibp_probability(t['obstime_mjd'], t[f'{tel}_ipp_lon'], f107)
        out.update(
            {
                f'{tel}_ibp': ibp,
                f'{tel}_apex_km': apex,
                f'{tel}_p_reach': 1.0 / (1.0 + np.exp(-_design(ibp, apex, f107) @ coef)),
            }
        )
    return out


def _fetch_retry(job, tries=4):
    for k in range(tries):
        try:
            fetch_ibi_night(*job)
            return True
        except (OSError, http.client.HTTPException) as e:  # noqa: PERF203 - transient server errors
            print(f'night {job[0]} sat {job[1]} try {k + 1}: {e!r}', flush=True)
            time.sleep(5 * (k + 1))
    return False


def _corpus_nights(oi_table):
    d = np.load(oi_table)
    names = [str(n) for n in d['ctx_names']]
    night = np.floor(d['obstime_mjd'] - 0.5).astype(int)  # MJD of the evening: UT 12h-36h is one night
    f = d['ctx_sci'][:, names.index('f107')]
    return {int(n): float(np.median(f[night == n])) for n in np.unique(night)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=('fetch', 'fit'))
    ap.add_argument('--oi-table', required=True)
    args = ap.parse_args()
    f107 = _corpus_nights(args.oi_table)
    if args.step == 'fetch':
        jobs = [(n, sat) for n in sorted(f107) for sat in SATS]
        with ThreadPoolExecutor(4) as pool:
            failed = [j for j, ok in zip(jobs, pool.map(_fetch_retry, jobs), strict=True) if not ok]
        print(f'{len(jobs) - len(failed)}/{len(jobs)} night-satellite files cached; failed: {failed}')
        return
    passes = ibi_passes(sorted(f107))
    np.savez(IBI_DIR / 'passes.npz', **passes)
    calib, _ = fit_reach(passes, f107)
    CALIBRATION.write_text(json.dumps(calib, indent=1) + '\n')
    print(json.dumps(calib, indent=1))


if __name__ == '__main__':
    main()
