"""GNSS ROTI near the LVM pierce points: small-scale ionospheric irregularities (plasma bubbles) seen by GPS.

For each ground receiver and GPS satellite, the geometry-free carrier phase gives slant TEC up to a constant:
STEC = (lambda1 L1 - lambda2 L2) / (40.3e16 (1/f2^2 - 1/f1^2)) [TECU].  ROT is its rate between consecutive 30 s
epochs [TECU/min]; ROTI is the standard deviation of ROT in fixed 5 min bins (Pi et al. 1997).  Each ROTI sample is
placed at the ray's pierce point on the same 285 km shell used for the LVM telescopes.  Satellites below 30 deg
elevation are dropped (multipath); a loss-of-lock flag, a data gap or |ROT| > 20 TECU/min starts a new arc.

ROTI describes structure along the ray, not a column amount: it carries no van Rhijn factor and enters SkyML
only as a local irregularity indicator per telescope, so the single geometry correction is untouched.

Data: IGS daily 30 s RINEX 3 (Hatanaka) observation files and the merged broadcast navigation file from the BKG
GNSS Data Center (public, no login).  Each station-day is reduced to its ROTI samples and cached in ``ROTI_DIR``;
the raw RINEX is not kept.  Run from ``skysub/``: ``python -m mlp_predictor.roti fetch --nights FILE``, then
``python -m mlp_predictor.roti table --iono-table PATH --out PATH`` for the per-exposure table SkyML reads.
"""

import argparse
import gzip
import http.client
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from mlp_predictor import ionosphere

BKG = 'https://igs.bkg.bund.de/root_ftp/IGS'
OBS_URLS = tuple(BKG + '/obs/{y}/{doy:03d}/{sta}_' + k + '_{y}{doy:03d}0000_01D_30S_MO.crx.gz' for k in ('R', 'S'))
NAV_URLS = tuple(
    BKG + '/BRDC/{y}/{doy:03d}/' + f + '_{y}{doy:03d}0000_01D_' + t + '.rnx.gz'
    for f, t in (('BRDC00WRD_S', 'MN'), ('BRDC00WRD_R', 'MN'), ('BRDC00IGS_R', 'MN'), ('BRDC00WRD_R', 'GN'))
)
MIN_NAV_GPS = 24  # GPS satellites of a complete daily file (~31); BRDC00WRD_S has 8-11 on MJD 60181-60213
STATIONS = (
    'OAFA00ARG',
    'SANT00CHL',
    'CORD00ARG',
    'UNSA00ARG',
    'MGUE00ARG',
)  # San Juan, Santiago, Cordoba, Salta, Malargue
ROTI_DIR = ionosphere.CACHE_DIR / 'roti'
F1, F2, C_LIGHT = 1575.42e6, 1227.60e6, 299792458.0
TECU_PER_M = 1.0 / (40.3e16 * (1.0 / F2**2 - 1.0 / F1**2))
ELEV_MASK_DEG, ROT_SLIP, BIN_S, MIN_ROT = 30.0, 20.0, 300.0, 8
DETREND_MIN = 60.0  # high-pass window of the arc TEC: a 200 km bubble at 100 m/s takes ~35 min to pass
MER_SCALE = 4.0  # bubbles are field-aligned: a magnetic-meridian (~north-south at LCO) offset counts 1/4
SAMPLE_KEYS = ('mjd', 'prn', 'el', 'ipp_lat', 'ipp_lon', 'roti', 'dvtec')
MISSING_AFTER_DAYS = 4.0  # a 404 for a younger day is not cached (BKG publishes daily RINEX with a delay)
GM, OMEGA_E = 3.986005e14, 7.2921151467e-5
_GPS0_MJD = 44244.0  # 1980-01-06


def _get(url, timeout=180):
    req = urllib.request.Request(url, headers={'User-Agent': 'lvmsky/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _get_first(urls, **fmt):
    """Bytes of the first URL that exists; HTTPError 404 if none does (BKG naming varies by day)."""
    for i, u in enumerate(urls):
        try:
            return _get(u.format(**fmt))
        except urllib.error.HTTPError as e:
            if e.code != 404 or i == len(urls) - 1:
                raise
    raise AssertionError


def _get_nav(y, doy):
    """Bytes of the first navigation file with >= MIN_NAV_GPS GPS satellites, else of the most complete one."""
    best, n_best, missing = None, -1, None
    for u in NAV_URLS:
        try:
            raw = _get(u.format(y=y, doy=doy))
        except urllib.error.HTTPError as e:
            if e.code != 404:
                raise
            missing = e
            continue
        n = len({ln[:3] for ln in gzip.decompress(raw).decode().splitlines() if ln.startswith('G')})
        if n >= MIN_NAV_GPS:
            return raw
        if n > n_best:
            best, n_best = raw, n
    if best is None:
        raise missing
    return best


def _day(mjd_day):
    d = ionosphere._MJD0 + ionosphere.dt.timedelta(days=int(mjd_day))
    return d.year, d.timetuple().tm_yday


def parse_nav(text):
    """GPS broadcast ephemerides of a RINEX 3 navigation file, as a dict of arrays (one row per record)."""
    body = text[text.index('END OF HEADER') :].splitlines()[1:]
    rows = []
    for i, line in enumerate(body):
        if not line.startswith('G') or i + 7 >= len(body):
            continue
        rec = [line[4 + 19 * k : 4 + 19 * (k + 1)] for k in range(1, 4)]
        for ln in body[i + 1 : i + 8]:
            rec += [ln[4 + 19 * k : 4 + 19 * (k + 1)] for k in range(4)]
        v = [float(x.replace('D', 'E')) if x.strip() else 0.0 for x in rec]
        rows.append([int(line[1:3])] + v)
    a = np.array(rows)
    names = (
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
    return {n: a[:, k] for k, n in enumerate(names)}


def sat_ecef(eph, prn, gps_sow):
    """ECEF positions [m] of satellite ``prn`` at GPS seconds of week ``gps_sow`` (IS-GPS-200 broadcast model)."""
    sel = np.flatnonzero((eph['prn'] == prn) & (eph['health'] == 0))
    if sel.size == 0:
        return np.full((np.size(gps_sow), 3), np.nan)
    t = np.atleast_1d(gps_sow)
    dt_all = t[:, None] - eph['toe'][sel][None, :]
    dt_all = (dt_all + 302400) % 604800 - 302400
    k = sel[np.argmin(np.abs(dt_all), axis=1)]
    e = {n: v[k] for n, v in eph.items()}
    tk = (t - e['toe'] + 302400) % 604800 - 302400
    a = e['sqrta'] ** 2
    m = e['m0'] + (np.sqrt(GM / a**3) + e['dn']) * tk
    ecc = m.copy()
    for _ in range(10):
        ecc = m + e['e'] * np.sin(ecc)
    nu = np.arctan2(np.sqrt(1 - e['e'] ** 2) * np.sin(ecc), np.cos(ecc) - e['e'])
    phi = nu + e['w']
    u = phi + e['cus'] * np.sin(2 * phi) + e['cuc'] * np.cos(2 * phi)
    r = a * (1 - e['e'] * np.cos(ecc)) + e['crs'] * np.sin(2 * phi) + e['crc'] * np.cos(2 * phi)
    inc = e['i0'] + e['idot'] * tk + e['cis'] * np.sin(2 * phi) + e['cic'] * np.cos(2 * phi)
    om = e['omega0'] + (e['omegadot'] - OMEGA_E) * tk - OMEGA_E * e['toe']
    x, y = r * np.cos(u), r * np.sin(u)
    out = np.column_stack(
        [x * np.cos(om) - y * np.cos(inc) * np.sin(om), x * np.sin(om) + y * np.cos(inc) * np.cos(om), y * np.sin(inc)]
    )
    out[np.abs(tk) > 4 * 3600] = np.nan
    return out


def _geodetic(xyz):
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    x, y, z = xyz
    lon, p = np.arctan2(y, x), np.hypot(x, y)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(6):
        n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        lat = np.arctan2(z + e2 * n * np.sin(lat), p)
    return lat, lon


def el_az(rx, sat):
    """Elevation and azimuth [deg] of satellite ECEF positions ``sat`` (N x 3) from receiver ECEF ``rx``."""
    lat, lon = _geodetic(rx)
    d = sat - np.asarray(rx)
    sl, cl, so, co = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    e = -so * d[:, 0] + co * d[:, 1]
    n = -sl * co * d[:, 0] - sl * so * d[:, 1] + cl * d[:, 2]
    u = cl * co * d[:, 0] + cl * so * d[:, 1] + sl * d[:, 2]
    return np.degrees(np.arctan2(u, np.hypot(e, n))), np.degrees(np.arctan2(e, n)) % 360.0


def parse_obs(text):
    """GPS L1/L2 carrier phase [cycles] and loss-of-lock flags of a RINEX 3 observation file."""
    head, body = text.split('END OF HEADER', 1)
    types, rx, system = [], None, None
    for line in head.splitlines():
        if line[60:].startswith('SYS / # / OBS TYPES'):
            system = line[0] if line[0] != ' ' else system
            if system == 'G':
                types += line[7:60].split()
        elif line[60:].startswith('APPROX POSITION XYZ'):
            rx = np.array(line[:42].split(), float)
    i2 = next((types.index(t) for t in ('L2W', 'L2X', 'L2L', 'L2S') if t in types), None)
    if 'L1C' not in types or i2 is None:
        raise ValueError('RINEX file has no GPS L1C/L2 carrier-phase observations')
    i1 = types.index('L1C')
    t_l, prn_l, l1_l, l2_l, lli_l = [], [], [], [], []
    sow = None
    for line in body.splitlines()[1:]:
        if line.startswith('>'):
            f = line[1:].split()
            y, mo, d, h, mi, sec = int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), float(f[5])
            mjd = ionosphere._mjd(y, mo, d) + (h * 3600 + mi * 60 + sec) / 86400.0
            sow = ((mjd - _GPS0_MJD) % 7) * 86400.0
            continue
        if line[:1] != 'G' or sow is None:
            continue
        f1, f2 = line[3 + 16 * i1 : 3 + 16 * i1 + 16], line[3 + 16 * i2 : 3 + 16 * i2 + 16]
        if not f1[:14].strip() or not f2[:14].strip():
            continue
        t_l.append(mjd)
        prn_l.append(int(line[1:3]))
        l1_l.append(float(f1[:14]))
        l2_l.append(float(f2[:14]))
        lli_l.append(any(f[14:15].strip() not in ('', '0') and int(f[14]) & 1 for f in (f1, f2)))
    return dict(
        rx=rx, mjd=np.array(t_l), prn=np.array(prn_l), l1=np.array(l1_l), l2=np.array(l2_l), lli=np.array(lli_l)
    )


def _arc_highpass(t, stec, arc):
    """TEC minus its centred DETREND_MIN running mean within each arc; NaN for arcs shorter than the window."""
    out = np.full(stec.size, np.nan)
    w = int(round(DETREND_MIN * 60.0 / 30.0))
    for a in np.unique(arc):
        k = np.flatnonzero(arc == a)
        if t[k[-1]] - t[k[0]] < DETREND_MIN * 60.0:
            continue
        c = np.concatenate([[0.0], np.cumsum(stec[k])])
        i = np.arange(k.size)
        hw = np.minimum(w // 2, np.minimum(i, k.size - 1 - i))  # symmetric: a linear trend cancels exactly
        d = stec[k] - (c[i + hw + 1] - c[i - hw]) / (2 * hw + 1)
        out[k] = np.where(hw >= w // 4, d, np.nan)  # no estimate within a quarter window of the arc ends
    return out


def mapping(el_deg, h_km=ionosphere.IPP_HEIGHT_KM):
    """Thin-shell slant factor STEC / VTEC at ``h_km`` (the same geometry as the van Rhijn factor of the shell)."""
    r = ionosphere.EARTH_RADIUS_KM
    return 1.0 / np.sqrt(1.0 - (r * np.cos(np.radians(el_deg)) / (r + h_km)) ** 2)


def roti_samples(obs, eph, h_km=ionosphere.IPP_HEIGHT_KM):
    """5-min ROTI and high-passed vertical TEC (dvtec, TECU; negative = depletion) per satellite arc and bin.

    Each sample sits at the pierce point of the bin's mid-time ray.  ``dvtec`` converts the slant TEC perturbation
    to vertical with the shell factor, so it describes the plasma at that point, not the line of sight.
    """
    rx_lat, rx_lon = (np.degrees(v) for v in _geodetic(obs['rx']))
    out = []
    for prn in np.unique(obs['prn']):
        s = np.flatnonzero(obs['prn'] == prn)
        t = obs['mjd'][s] * 86400.0
        stec = (C_LIGHT / F1 * obs['l1'][s] - C_LIGHT / F2 * obs['l2'][s]) * TECU_PER_M
        dt = np.diff(t)
        rot = np.diff(stec) / (dt / 60.0)
        ok = (np.abs(dt - 30.0) < 1.0) & ~obs['lli'][s][1:] & (np.abs(rot) < ROT_SLIP)
        dstec = _arc_highpass(t, stec, np.concatenate([[0], np.cumsum(~ok)]))[1:]
        tb = np.floor(t[1:] / BIN_S)
        for b in np.unique(tb[ok]):
            k = ok & (tb == b)
            if k.sum() < MIN_ROT:
                continue
            tm = (b + 0.5) * BIN_S
            pos = sat_ecef(eph, prn, ((tm / 86400.0 - _GPS0_MJD) % 7) * 86400.0)
            el, az = el_az(obs['rx'], pos)
            if not np.isfinite(el[0]) or el[0] < ELEV_MASK_DEG:
                continue
            lat, lon, _, _ = ionosphere.pierce_point(el[0], az[0], h_km, rx_lat, rx_lon)
            dv = np.nanmedian(dstec[k]) / mapping(el[0], h_km) if np.isfinite(dstec[k]).any() else np.nan
            out.append((tm / 86400.0, prn, el[0], lat, lon, np.std(rot[k]), dv))
    a = np.array(out, float).reshape(-1, len(SAMPLE_KEYS))
    return {n: a[:, i] for i, n in enumerate(SAMPLE_KEYS)}


def station_day(sta, mjd_day, download=True):
    """Cached ROTI samples of one station and UT day; an empty table with ``missing=True`` if the file is absent."""
    path = ROTI_DIR / f'{sta}_{int(mjd_day)}_v2.npz'  # v2 adds dvtec
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        import hatanaka

        y, doy = _day(mjd_day)
        try:
            raw = _get_first(OBS_URLS, y=y, doy=doy, sta=sta)
        except urllib.error.HTTPError as e:
            # A day can still be in flight on the archive: remember "no file" only once it is old enough.
            if e.code != 404 or time.time() / 86400.0 + 40587.0 - mjd_day < MISSING_AFTER_DAYS:
                raise
            res, missing = roti_samples_empty(), True
        else:
            eph = _nav(mjd_day)
            res, missing = roti_samples(parse_obs(hatanaka.decompress(raw).decode()), eph), False
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f'.{path.stem}.{os.getpid()}.part.npz')  # concurrent SkyML workers
        np.savez_compressed(tmp, missing=missing, **res)
        tmp.replace(path)
    return dict(np.load(path))


def roti_samples_empty():
    return {n: np.zeros(0) for n in SAMPLE_KEYS}


def _nav(mjd_day):
    path = ROTI_DIR / 'nav' / f'brdc_{int(mjd_day)}.rnx.gz'
    if not path.exists():
        y, doy = _day(mjd_day)
        raw = _get_nav(y, doy)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f'.{path.name}.{os.getpid()}.part')
        tmp.write_bytes(raw)
        tmp.replace(path)
    return parse_nav(gzip.decompress(path.read_bytes()).decode())


def night_samples(night, stations=STATIONS):
    """All stations' ROTI samples for the two UT days a night spans (evening MJD ``night``)."""
    parts = [station_day(s, d, download=False) for s in stations for d in (night, night + 1)]
    return {k: np.concatenate([p[k] for p in parts]) for k in SAMPLE_KEYS}


def _dist_km(lat1, lon1, lat2, lon2):
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    c = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dl)
    return (ionosphere.EARTH_RADIUS_KM + ionosphere.IPP_HEIGHT_KM) * np.arccos(np.clip(c, -1, 1))


def _aligned_dist_km(lat0, lon0, lat, lon):
    """Distance with the north-south (magnetic-meridian, declination ~-2 deg at LCO) offset divided by MER_SCALE."""
    rr = ionosphere.EARTH_RADIUS_KM + ionosphere.IPP_HEIGHT_KM
    dn = rr * np.radians(lat - lat0)
    de = rr * np.cos(np.radians(lat0)) * np.radians((lon - lon0 + 180.0) % 360.0 - 180.0)
    return np.hypot(de, dn / MER_SCALE)


def features_at(samples, mjd, lat, lon, radius_km=200.0, window_min=10.0, aligned=False):
    """Per query: max and median ROTI, sample count, min and median dvtec of samples within ``radius_km`` and
    ``window_min``; ``aligned`` measures distance along the field (north-south counts 1/MER_SCALE)."""
    out = np.full((np.size(mjd), 5), np.nan)
    out[:, 2] = 0
    order = np.argsort(samples['mjd'], kind='stable')
    t = samples['mjd'][order]
    s_lat, s_lon, s_roti = (samples[k][order] for k in ('ipp_lat', 'ipp_lon', 'roti'))
    dvtec = samples.get('dvtec', np.full(t.size, np.nan))[order]
    dist = _aligned_dist_km if aligned else _dist_km
    w = window_min / 1440.0
    for i, (m, la, lo) in enumerate(zip(np.atleast_1d(mjd), np.atleast_1d(lat), np.atleast_1d(lon), strict=True)):
        sl = slice(np.searchsorted(t, m - w, 'left'), np.searchsorted(t, m + w, 'right'))
        k = dist(la, lo, s_lat[sl], s_lon[sl]) <= radius_km
        if k.any():
            r, d = s_roti[sl][k], dvtec[sl][k]
            d = d[np.isfinite(d)]
            out[i] = (r.max(), np.median(r), k.sum(), d.min() if d.size else np.nan, np.median(d) if d.size else np.nan)
    return out


def dvtec_event(samples, mjd, lat, lon, window_min=5.0, min_el=45.0, max_km=400.0, drift_kms=0.0):
    """Per query, the single ray nearest the pierce point: its dvtec [TECU], field-aligned distance [km], elevation.

    "Nearest" uses ``_aligned_dist_km`` among samples within ``window_min`` and above ``min_el`` (less multipath and
    a smaller slant-to-vertical error); ``drift_kms`` moves each sample east by drift x (query time - sample time),
    for plasma that drifts with the eastward night-time flow.  Queries with no ray within ``max_km`` get NaN.
    """
    out = np.full((np.size(mjd), 3), np.nan)
    ok = (samples['el'] >= min_el) & np.isfinite(samples['dvtec'])
    order = np.flatnonzero(ok)[np.argsort(samples['mjd'][ok], kind='stable')]
    t, s_lat, s_lon, s_dv, s_el = (samples[k][order] for k in ('mjd', 'ipp_lat', 'ipp_lon', 'dvtec', 'el'))
    w = window_min / 1440.0
    for i, (m, la, lo) in enumerate(zip(np.atleast_1d(mjd), np.atleast_1d(lat), np.atleast_1d(lon), strict=True)):
        sl = slice(np.searchsorted(t, m - w, 'left'), np.searchsorted(t, m + w, 'right'))
        if sl.stop <= sl.start:
            continue
        shift = drift_kms * (m - t[sl]) * 86400.0 / (111.2 * np.cos(np.radians(la)))  # deg of longitude
        d = _aligned_dist_km(la, lo, s_lat[sl], s_lon[sl] + shift)
        j = int(np.argmin(d))
        if d[j] <= max_km:
            out[i] = (s_dv[sl][j], d[j], s_el[sl][j])
    return out


def exposure_flag(mjd, ipp, threshold=1.0, radius_km=200.0, window_min=10.0, stations=STATIONS, download=True):
    """Plasma-bubble QA for one exposure: max ROTI near each telescope's pierce point and a tri-state flag.

    ``ipp`` maps a telescope label to its (lat, lon) at 285 km.  ``bubble`` is 1 when any telescope has a ROTI
    sample >= ``threshold``, 0 when at least one telescope has samples and none reaches it, and -1 (unknown) when
    no telescope has samples or the data cannot be obtained; this never raises.
    """
    night = int(np.floor(mjd - 0.5))
    used, errors, parts = [], [], []
    for sta in stations:
        try:
            got = [station_day(sta, d, download=download) for d in (night, night + 1)]
        except Exception as e:  # noqa: BLE001 - QA must not stop the subtraction
            errors.append(f'{sta}: {e!r}')
            continue
        if not all(bool(g['missing']) for g in got):
            used.append(sta)
        parts += got
    out = {
        'stations': used,
        'errors': errors,
        'threshold': threshold,
        'radius_km': radius_km,
        'window_min': window_min,
        'telescopes': {},
    }
    smp = {
        k: np.concatenate([p.get(k, np.full(p['mjd'].size, np.nan)) for p in parts]) if parts else np.zeros(0)
        for k in ('mjd', 'ipp_lat', 'ipp_lon', 'roti', 'dvtec')
    }
    for label, (lat, lon) in ipp.items():
        f = features_at(smp, np.array([mjd]), np.array([lat]), np.array([lon]), radius_km, window_min)[0]
        g = features_at(smp, np.array([mjd]), np.array([lat]), np.array([lon]), radius_km, window_min, aligned=True)[0]
        out['telescopes'][label] = {
            'roti_max': float(f[0]),
            'n_samples': int(f[2]),
            'dvtec_min': float(g[3]),
            'dvtec_median': float(g[4]),
            'ipp_lat': float(lat),
            'ipp_lon': float(lon),
        }
    peaks = [v['roti_max'] for v in out['telescopes'].values() if v['n_samples'] > 0]
    out['bubble'] = -1 if not peaks else int(max(peaks) >= threshold)
    return out


TABLE_DIR = ionosphere.CONTEXT_TABLE_DIR / 'roti'
ROLES = ionosphere.ROLES


def roti_at(mjd, lat, lon, stations=STATIONS, download=False, radius_km=200.0, window_min=10.0, aligned=False):
    """Max and median ROTI [TECU/min] and sample count near each (mjd, lat, lon) query (see ``features_at``).

    Queries are grouped by night; a station-day that is absent (or cannot be fetched) simply contributes no samples,
    so a query with no samples gets count 0 and NaN ROTI.
    """
    mjd, lat, lon = (np.atleast_1d(np.asarray(a, np.float64)) for a in (mjd, lat, lon))
    out = np.full((mjd.size, 3), np.nan)
    out[:, 2] = 0
    night = np.floor(mjd - 0.5).astype(int)
    for n in np.unique(night):
        parts = []
        for sta in stations:
            for d in (n, n + 1):
                try:
                    parts.append(station_day(sta, d, download=download))
                except Exception:  # noqa: BLE001 - a missing station-day only removes samples
                    pass
        if not parts:
            continue
        smp = {k: np.concatenate([p.get(k, np.full(p['mjd'].size, np.nan)) for p in parts]) for k in SAMPLE_KEYS}
        rows = night == n
        out[rows] = features_at(smp, mjd[rows], lat[rows], lon[rows], radius_km, window_min, aligned)[:, :3]
    return out


def build_table(iono_table, out, stations=STATIONS):
    """Per-role ROTI for every row of a bundled ionosphere context table, from the cached station-days only.

    Field-aligned search (``aligned=True``, 200 km, 10 min): bubbles are elongated along the magnetic meridian, and
    the aligned distance lifts the share of rows with samples from ~60% to ~92%.
    """
    t = np.load(iono_table)
    table = {'obstime_mjd': t['obstime_mjd'].astype(np.float64)}
    for role in ROLES:
        f = roti_at(table['obstime_mjd'], t[f'{role}_ipp_lat'], t[f'{role}_ipp_lon'], stations, aligned=True)
        table |= {f'{role}_roti_max': f[:, 0].astype(np.float32), f'{role}_roti_med': f[:, 1].astype(np.float32),
                  f'{role}_roti_n': f[:, 2].astype(np.int32)}
        print(f'{role}: {np.mean(f[:, 2] > 0):.1%} of rows have samples', flush=True)
    np.savez_compressed(out, **table)
    return table


def bundled_table():
    """Bundled ROTI tables, row-aligned with ``ionosphere._bundled_tables()``."""
    tables = [dict(np.load(p, allow_pickle=False)) for p in sorted(TABLE_DIR.glob('*.npz'))]
    iono = ionosphere._bundled_tables()
    if not tables or iono is None:
        return None
    merged = {k: np.concatenate([t[k] for t in tables]) for k in tables[0]}
    order = np.argsort(merged['obstime_mjd'], kind='stable')
    merged = {k: v[order] for k, v in merged.items()}
    if not np.array_equal(merged['obstime_mjd'], iono['obstime_mjd']):
        raise ValueError('bundled ROTI and ionosphere tables are not row-aligned; rebuild the former')
    return merged


def _fetch_retry(job, tries=4):
    for k in range(tries):
        try:
            station_day(*job)
            return True
        except (OSError, http.client.HTTPException, ValueError, EOFError) as e:  # noqa: PERF203 - transient errors
            print(f'{job} try {k + 1}: {e!r}', flush=True)
            time.sleep(5 * (k + 1))
    return False


def fetch(nights, workers=4, stations=STATIONS):
    """Download and reduce the station-days of the given evening MJDs that are not cached; returns the failures."""
    days = sorted({int(d) for n in nights for d in (n, n + 1)})
    for d in days:  # navigation first, serially (shared by all stations)
        try:
            _nav(d)
        except (OSError, http.client.HTTPException) as e:
            print(f'navigation {d}: {e!r}', flush=True)
    jobs = [(s, d) for d in days for s in stations]
    with ThreadPoolExecutor(workers) as pool:
        failed = [j for j, ok in zip(jobs, pool.map(_fetch_retry, jobs), strict=True) if not ok]
    print(f'{len(jobs) - len(failed)}/{len(jobs)} station-days done; failed: {failed}')
    return failed


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=('fetch', 'table'))
    ap.add_argument('--nights', help='fetch: text file with one evening MJD per line')
    ap.add_argument('--iono-table', help='table: bundled ionosphere context table')
    ap.add_argument('--out', help='table: output .npz')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--stations', nargs='+', default=list(STATIONS))
    args = ap.parse_args()
    if args.step == 'table':
        build_table(args.iono_table, args.out, tuple(args.stations))
        return
    fetch(np.loadtxt(args.nights, dtype=int, ndmin=1), args.workers, tuple(args.stations))


if __name__ == '__main__':
    main()
