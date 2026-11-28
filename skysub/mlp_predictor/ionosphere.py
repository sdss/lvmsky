"""Ionospheric context for SkyML: IONEX vertical TEC at the F-region pierce point of each LVM telescope.

Geometry contract -- the van Rhijn factor is applied exactly once, elsewhere:
``data.airglow_geometry_scale`` divides every airglow coefficient by V(z, 285 km) (times extinction), so the
network models zenith-equivalent emission.  Here the line of sight enters ONLY through the horizontal
position of the ionospheric pierce point (IPP).  Every TEC value returned is VERTICAL TEC at that point,
exactly as published in the IONEX map; nothing is multiplied or divided by a mapping function or by V.

Maps, in order of preference: UPC quarter-hourly UQRG, UPC rapid 2-hourly UPRG, CODE final 1-hourly.  A day
whose file is missing, truncated or unparsable falls through to the next source (UPC serves a truncated
UQRG for 2025-11-26, for example); the source used is returned per row.  Daily files are cached in
``CACHE_DIR`` (override with the ``LVMSKY_IONEX_DIR`` environment variable) and never re-downloaded once
present.  A time no source covers raises instead of producing NaN.

Run ``python -m mlp_predictor.ionosphere --oi-table PATH --out PATH`` from ``skysub/`` to build the
per-exposure table for the Science, SkyE and SkyW telescopes.
"""

import argparse
import datetime as dt
import http.client
import os
import subprocess
import urllib.request
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

IPP_HEIGHT_KM = 285.0  # same F-region layer as data.LAYER_HEIGHTS_KM['ionospheric_red']
EARTH_RADIUS_KM = 6371.0
LCO_LAT_DEG, LCO_LON_DEG = -29.00889, -70.68992
TELESCOPES = ('sci', 'skye', 'skyw')
_SKYSUB = Path(__file__).resolve().parents[1]
CACHE_DIR = Path(os.environ.get('LVMSKY_IONEX_DIR', _SKYSUB / 'ionosphere_cache'))
GEOMAG_GRID = _SKYSUB / 'sky_decomp' / 'data' / 'geomag' / 'igrf_285km_lco.npz'
SOURCES = {  # tried in this order
    'UQRG': 'https://chapman.upc.es/tomion/rapid/{y}/{doy:03d}_{yymmdd}.15min/uqrg{doy:03d}0.{yy}i.Z',
    'UPRG': 'https://chapman.upc.es/tomion/rapid/{y}/{doy:03d}_{yymmdd}.15min/uprg{doy:03d}0.{yy}i.Z',
    'CODG': 'http://ftp.aiub.unibe.ch/CODE/{y}/COD0OPSFIN_{y}{doy:03d}0000_01D_01H_GIM.INX.gz',
}
_MJD0 = dt.datetime(1858, 11, 17)


def _mjd(y, mo, d, h=0, mi=0, s=0):
    """MJD of a calendar epoch; IONEX writes the end of day as hour 24."""
    return (dt.datetime(y, mo, d) - _MJD0).days + (h * 3600 + mi * 60 + s) / 86400.0


def ionex_file(mjd_day, source, download=True):
    """Cached path of one source's IONEX file for integer MJD ``mjd_day``, downloading it if missing."""
    day = _MJD0 + dt.timedelta(days=int(mjd_day))
    url = SOURCES[source].format(y=day.year, yy=f'{day.year % 100:02d}', doy=day.timetuple().tm_yday,
                                 yymmdd=day.strftime('%y%m%d'))
    path = CACHE_DIR / source / url.rsplit('/', 1)[1]
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={'User-Agent': 'lvmsky/1.0'})
        with urllib.request.urlopen(req, timeout=120) as resp:
            part = path.with_name(f'.{path.name}.{os.getpid()}.part')  # unique per worker: no shared temp file
            part.write_bytes(resp.read())
        part.replace(path)
    return path


def _axis(first, last, step):
    return np.arange(first, last + step / 2, step)


def _rows_to_grid(rows, nlat, nlon):
    """(nlat, nlon) grid from per-latitude rows; a row of the wrong length stays NaN."""
    grid = np.full((nlat, nlon), np.nan)
    good = [i for i, row in enumerate(rows[:nlat]) if len(row) == nlon]
    for i in good:
        grid[i] = rows[i]
    return grid, nlat - len(good)


def read_ionex(path):
    """Parse an IONEX file: epochs [MJD], native lat/lon axes [deg], TEC and RMS maps [TECU].

    Missing values (9999, or the non-standard ``nan`` token some UQRG files carry) and any latitude row
    whose length differs from the longitude axis are returned as NaN; ``bad_rows`` counts those rows.
    """
    path = Path(path)
    raw = (subprocess.run(['gzip', '-dc', str(path)], check=True, capture_output=True).stdout
           if path.suffix in ('.Z', '.gz') else path.read_bytes())
    head, exponent, epochs = {}, -1, []
    maps, bad_rows = {'TEC': [], 'RMS': []}, {'TEC': 0, 'RMS': 0}
    kind, rows, map_exp = None, [], exponent
    for line in raw.decode('ascii', 'replace').splitlines():
        label = line[60:].strip()
        if not any(c.isalpha() for c in label):  # data record
            if kind and rows:
                text = line.rstrip()
                rows[-1].extend(float(text[k:k + 5]) for k in range(0, len(text), 5))
            continue
        if label in ('LAT1 / LAT2 / DLAT', 'LON1 / LON2 / DLON', 'HGT1 / HGT2 / DHGT'):
            head[label.split()[0]] = [float(v) for v in line[:60].split()]
        elif label == 'MAPPING FUNCTION':
            head['mapping'] = line[:60].strip()
        elif label == 'EXPONENT':  # header default, or inside a map (then it applies to that map only)
            if kind:
                map_exp = int(line[:60].split()[0])
            else:
                exponent = int(line[:60].split()[0])
        elif label in ('START OF TEC MAP', 'START OF RMS MAP'):
            kind, rows, map_exp = label.split()[2], [], exponent
        elif label == 'EPOCH OF CURRENT MAP' and kind == 'TEC':
            epochs.append(_mjd(*(int(v) for v in line[:36].split())))
        elif label.startswith('LAT/LON1') and kind:
            rows.append([])
        elif label in ('END OF TEC MAP', 'END OF RMS MAP') and kind:
            grid, n_bad = _rows_to_grid(rows, _axis(*head['LAT1']).size, _axis(*head['LON1']).size)
            bad_rows[kind] += n_bad
            maps[kind].append(np.where(grid == 9999, np.nan, grid * 10.0 ** map_exp))
            kind = None
    lat, lon = _axis(*head['LAT1']), _axis(*head['LON1'])
    shape = (len(epochs), lat.size, lon.size)
    tec = np.asarray(maps['TEC']).reshape(shape)
    rms = np.asarray(maps['RMS']).reshape(shape) if maps['RMS'] else np.full(shape, np.nan)
    return dict(epoch=np.asarray(epochs), lat=lat, lon=lon, tec=tec, rms=rms, bad_rows=bad_rows,
                height_km=head['HGT1'][0], mapping=head.get('mapping', ''))


@lru_cache(maxsize=4)
def _day_maps(mjd_day, download=True):
    """Maps of the first source in SOURCES whose file for this day downloads and parses."""
    errors = []
    for source in SOURCES:
        try:
            path = ionex_file(mjd_day, source, download=download)
            maps = read_ionex(path)
            if maps['epoch'].size < 2:
                raise ValueError(f'{path.name}: fewer than two maps')
            return maps | {'source': source, 'file': path.name}
        except (OSError, ValueError, KeyError, subprocess.CalledProcessError, http.client.HTTPException) as exc:
            errors.append(f'{source}: {type(exc).__name__}: {exc}')
    day = _MJD0 + dt.timedelta(days=int(mjd_day))
    raise FileNotFoundError(f'No IONEX map for MJD {int(mjd_day)} ({day:%Y-%m-%d}); {"; ".join(errors)}. '
                            f'Place a daily file under {CACHE_DIR}/<SOURCE>/.')


def _sample(maps, field, k, lat, lon):
    """Bilinear interpolation of ``maps[field][k]`` at (lat, lon) (IONEX 1.0 manual, eq. 2)."""
    lat_ax, lon_ax = maps['lat'], maps['lon']
    fy = (lat - lat_ax[0]) / (lat_ax[1] - lat_ax[0])
    fx = (np.mod(lon - lon_ax[0], 360.0)) / (lon_ax[1] - lon_ax[0])
    y0 = np.clip(np.floor(fy).astype(int), 0, lat_ax.size - 2)
    x0 = np.clip(np.floor(fx).astype(int), 0, lon_ax.size - 2)
    wy, wx = fy - y0, fx - x0
    g = maps[field]
    return ((1 - wy) * (1 - wx) * g[k, y0, x0] + (1 - wy) * wx * g[k, y0, x0 + 1]
            + wy * (1 - wx) * g[k, y0 + 1, x0] + wy * wx * g[k, y0 + 1, x0 + 1])


def vtec_at(mjd, lat, lon, *, download=True):
    """Vertical TEC and its RMS [TECU] at (mjd, lat, lon), plus the map source of each row.

    Rotated-map time interpolation between the two bracketing maps (IONEX 1.0 manual, eq. 3): each map is
    sampled at longitude lon + (t - T_i) * 360 deg/day, so the Sun-fixed ionosphere is not smeared.
    Each day file is parsed once per call, so pass every query of a batch together.
    """
    mjd, lat, lon = (np.atleast_1d(np.asarray(a, np.float64)) for a in (mjd, lat, lon))
    out = {'tec': np.full(mjd.shape, np.nan), 'rms': np.full(mjd.shape, np.nan)}
    source = np.empty(mjd.shape, dtype='<U4')
    for day in np.unique(np.floor(mjd)).astype(int):
        rows = np.floor(mjd) == day
        maps = _day_maps(day, download=download)
        source[rows] = maps['source']
        if mjd[rows].max() > maps['epoch'][-1]:  # last quarter-hour of the day: append next day's maps
            nxt = _day_maps(day + 1, download=download)
            maps = {k: np.concatenate([maps[k], nxt[k]]) for k in ('epoch', 'tec', 'rms')} | {
                'lat': maps['lat'], 'lon': maps['lon']}
        t = mjd[rows]
        i = np.clip(np.searchsorted(maps['epoch'], t, side='right') - 1, 0, maps['epoch'].size - 2)
        t0, t1 = maps['epoch'][i], maps['epoch'][i + 1]
        w = (t - t0) / (t1 - t0)
        for field in ('tec', 'rms'):
            e0 = _sample(maps, field, i, lat[rows], lon[rows] + (t - t0) * 360.0)
            e1 = _sample(maps, field, i + 1, lat[rows], lon[rows] + (t - t1) * 360.0)
            out[field][rows] = (1 - w) * e0 + w * e1
    for key in ('tec', 'rms'):
        if not np.isfinite(out[key]).all():
            raise ValueError(f'{int((~np.isfinite(out[key])).sum())} {key.upper()} samples are not finite (map gaps)')
    return out['tec'], out['rms'], source


def pierce_point(alt_deg, az_deg, h_km=IPP_HEIGHT_KM, lat0=LCO_LAT_DEG, lon0=LCO_LON_DEG):
    """Geographic IPP of a line of sight from LCO, and its north/east offset [km] along the shell."""
    el, az, phi0 = np.radians(alt_deg), np.radians(az_deg), np.radians(lat0)
    psi = np.pi / 2 - el - np.arcsin(EARTH_RADIUS_KM / (EARTH_RADIUS_KM + h_km) * np.cos(el))
    lat = np.arcsin(np.sin(phi0) * np.cos(psi) + np.cos(phi0) * np.sin(psi) * np.cos(az))
    lon = lon0 + np.degrees(np.arcsin(np.sin(psi) * np.sin(az) / np.cos(lat)))
    s = (EARTH_RADIUS_KM + h_km) * psi
    return np.degrees(lat), lon, s * np.cos(az), s * np.sin(az)


@lru_cache(maxsize=1)
def _geomag_interp():
    g = np.load(GEOMAG_GRID)
    axes = (g['year'], g['lat'], g['lon'])
    return (RegularGridInterpolator(axes, g['declination_deg'].astype(np.float64)),
            RegularGridInterpolator(axes, g['inclination_deg'].astype(np.float64)))


def geomag_at(lat, lon, mjd):
    """IGRF declination, inclination and dip latitude [deg] at 285 km from the bundled grid (any IPP height)."""
    year = np.clip(2000.0 + (np.asarray(mjd, np.float64) - 51544.5) / 365.25, 2023.0, 2027.0)
    pts = np.column_stack(np.broadcast_arrays(year, lat, lon))
    dec_i, inc_i = _geomag_interp()
    dec, inc = dec_i(pts), inc_i(pts)
    return dec, inc, np.degrees(np.arctan(np.tan(np.radians(inc)) / 2.0))


def telescope_features(mjd, alt_deg, az_deg, h_km=IPP_HEIGHT_KM, *, download=True):
    """Ionospheric features of telescope pointings; all TEC values are vertical (see module docstring).

    ``dlog10vtec_dt`` is the centred +-30 min rate [dex/h]; ``dlog10vtec_dlat`` / ``dlog10vtec_dlon`` are the
    local map gradients at the pierce point [dex/deg], from +-half a grid cell (1.25 deg lat, 2.5 deg lon).
    """
    mjd = np.atleast_1d(np.asarray(mjd, np.float64))
    lat, lon, dn, de = pierce_point(alt_deg, az_deg, h_km)
    half, dlat, dlon = 0.5 / 24.0, 1.25, 2.5
    queries = [(mjd, lat, lon), (mjd + half, lat, lon), (mjd - half, lat, lon), (mjd, lat + dlat, lon),
               (mjd, lat - dlat, lon), (mjd, lat, lon + dlon), (mjd, lat, lon - dlon)]
    tec, rms, source = vtec_at(*(np.concatenate(q) for q in zip(*queries, strict=True)), download=download)
    n = mjd.size
    lt = np.log10(np.clip(tec, 0.1, None)).reshape(len(queries), n)
    dec, _, dip_lat = geomag_at(lat, lon, mjd)
    dec0 = np.radians(geomag_at(LCO_LAT_DEG, LCO_LON_DEG, mjd)[0])
    return {
        'ipp_lat': lat, 'ipp_lon': lon, 'vtec': tec[:n], 'vtec_rms': rms[:n],
        'dvtec_dt': tec[n:2 * n] - tec[2 * n:3 * n],  # TECU per hour, centred
        'dlog10vtec_dt': lt[1] - lt[2],
        'dlog10vtec_dlat': (lt[3] - lt[4]) / (2 * dlat),
        'dlog10vtec_dlon': (lt[5] - lt[6]) / (2 * dlon),
        'dip_lat': dip_lat, 'mag_dec': dec,
        'ipp_dx_mag_mer': dn * np.cos(dec0) + de * np.sin(dec0),  # along the magnetic meridian, +north
        'ipp_dx_mag_zon': -dn * np.sin(dec0) + de * np.cos(dec0),  # across it, +magnetic east
    }, source[:n]


def sun_alt_conjugate(mjd):
    """Solar altitude [deg] at the magnetic conjugate point of the LCO zenith pierce point (bundled IGRF trace).

    Conjugate photoelectrons heat the LCO F region once the conjugate point is sunlit, before local dawn.
    """
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, get_sun
    from astropy.time import Time

    g = np.load(GEOMAG_GRID)
    mjd = np.atleast_1d(np.asarray(mjd, np.float64))
    year = np.clip(2000.0 + (mjd - 51544.5) / 365.25, g['year'][0], g['year'][-1])
    lat, lon = (np.interp(year, g['year'], g[k]) for k in ('conjugate_lat_deg', 'conjugate_lon_deg'))
    time = Time(mjd, format='mjd', scale='utc')
    frame = AltAz(obstime=time, location=EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=0 * u.m))
    return get_sun(time).transform_to(frame).alt.to_value(u.deg)


def build_exposure_table(expnum, mjd, pointings, h_km=IPP_HEIGHT_KM):
    """Per-exposure table; ``pointings`` maps each of TELESCOPES to (alt_deg, az_deg) arrays."""
    mjd = np.asarray(mjd, np.float64)
    n = mjd.size
    alt = np.concatenate([np.asarray(pointings[t][0], np.float64) for t in TELESCOPES])
    az = np.concatenate([np.asarray(pointings[t][1], np.float64) for t in TELESCOPES])
    feats, source = telescope_features(np.tile(mjd, len(TELESCOPES)), alt, az, h_km)  # one pass over day files
    table = {'expnum': np.asarray(expnum), 'obstime_mjd': mjd, 'ipp_height_km': float(h_km),
             'sun_alt_conjugate': sun_alt_conjugate(mjd)}
    for k, tel in enumerate(TELESCOPES):
        part = slice(k * n, (k + 1) * n)
        table |= {f'{tel}_alt': alt[part], f'{tel}_az': az[part]}
        table |= {f'{tel}_{name}': v[part] for name, v in feats.items()}
        table[f'{tel}_ionex_source'] = source[part]
    return table


def _pointings_from_oi_table(d):
    """Science / SkyE / SkyW (alt, az) from an oi_table.npz; near/far are mapped by the ``ew`` column."""
    names = [str(n) for n in d['ctx_names']]
    col = {n: names.index(n) for n in ('alt', 'az_sin', 'az_cos', 'ew')}

    def altaz(ctx):
        return ctx[:, col['alt']].astype(np.float64), np.degrees(
            np.arctan2(ctx[:, col['az_sin']], ctx[:, col['az_cos']])) % 360.0

    sci, near, far = (altaz(np.asarray(d[f'ctx_{a}'])) for a in ('sci', 'near', 'far'))
    near_is_e = np.asarray(d['ctx_near'])[:, col['ew']] > 0
    if not np.array_equal(near_is_e, np.asarray(d['ctx_far'])[:, col['ew']] < 0):
        raise ValueError('near/far ew labels are not one SkyE and one SkyW per row')
    pick = lambda a, b: tuple(np.where(near_is_e, x, y) for x, y in zip(a, b, strict=True))  # noqa: E731
    return {'sci': sci, 'skye': pick(near, far), 'skyw': pick(far, near)}


# --- Bundled per-corpus context tables --------------------------------------------------------------------------
# One npz per corpus in CONTEXT_TABLE_DIR holds the features of every META row for the three ctx roles (sci, near,
# far = the SCI / SKY_NEAR / SKY_FAR pointings), keyed by obstime_mjd.  The training augment reads them first, so a
# known corpus needs no IONEX download; rows not found (new exposures, another corpus) are computed from IONEX.
CONTEXT_TABLE_DIR = _SKYSUB / 'sky_decomp' / 'data' / 'ionosphere'
ROLES = ('sci', 'near', 'far')
_ROLE_KIND = {'sci': 'sci', 'near': 'sky1', 'far': 'sky2'}
MIN_TABLE_ALT_DEG = 5.0
_TABLE_FEATURES = ('vtec', 'vtec_rms', 'dvtec_dt', 'dlog10vtec_dt', 'dlog10vtec_dlat', 'dlog10vtec_dlon', 'dip_lat',
                   'mag_dec', 'ipp_dx_mag_mer', 'ipp_dx_mag_zon', 'ipp_lat', 'ipp_lon')


def build_context_table(meta_fits, out):
    """Precompute the bundled context table of a corpus META (pointing geometry exactly as data.py computes it)."""
    from astropy.table import Table

    from . import data

    meta = Table.read(meta_fits, hdu='META')
    upper = {str(c).upper(): c for c in meta.colnames}
    mjd = np.asarray(data._extract_obstime_mjd(meta, upper), np.float64)
    geo = {r: data._compute_astropy_geometry(meta, upper, _ROLE_KIND[r]) for r in ROLES}
    # Rows with an arm below MIN_TABLE_ALT_DEG carry unphysical pointing metadata (some reach -45 deg) and never
    # enter training; they are left out, so the augment would recompute (and fail loudly) if one ever appeared.
    ok = np.all([np.asarray(geo[r]['alt'], np.float64) >= MIN_TABLE_ALT_DEG for r in ROLES], axis=0)
    mjd = mjd[ok]
    alt = np.concatenate([np.asarray(geo[r]['alt'], np.float64)[ok] for r in ROLES])
    az = np.concatenate([(np.degrees(np.arctan2(geo[r]['az_sin'], geo[r]['az_cos'])) % 360.0)[ok] for r in ROLES])
    feats, source = telescope_features(np.tile(mjd, len(ROLES)), alt, az)
    n = mjd.size
    print(f'context table: {n} of {ok.size} META rows ({int((~ok).sum())} with an arm below {MIN_TABLE_ALT_DEG} deg)')
    table = {'obstime_mjd': mjd, 'ipp_height_km': IPP_HEIGHT_KM, 'source_meta': Path(meta_fits).name,
             'sun_alt_conjugate': sun_alt_conjugate(mjd).astype(np.float32)}
    for k, role in enumerate(ROLES):
        part = slice(k * n, (k + 1) * n)
        table |= {f'{role}_alt': alt[part].astype(np.float32), f'{role}_az': az[part].astype(np.float32),
                  f'{role}_ionex_source': source[part]}
        table |= {f'{role}_{f}': np.asarray(feats[f][part], np.float32) for f in _TABLE_FEATURES}
    np.savez_compressed(out, **table)
    return table


@lru_cache(maxsize=1)
def _bundled_tables():
    tables = [dict(np.load(p, allow_pickle=False)) for p in sorted(CONTEXT_TABLE_DIR.glob('*.npz'))]
    if not tables:
        return None
    keys = [k for k in tables[0] if np.ndim(tables[0][k]) == 1 and k != 'source_meta']
    merged = {k: np.concatenate([t[k] for t in tables]) for k in keys}
    order = np.argsort(merged['obstime_mjd'], kind='stable')
    return {k: v[order] for k, v in merged.items()}


def context_from_tables(mjd, alt_by_role, az_by_role, tol_deg=1e-3):
    """Rows of the bundled tables matching each query row by obstime and by all three arms' alt/az.

    Returns ``(found, table_rows)``: a boolean mask over the queries and, for the found ones, their table row index.
    """
    tab = _bundled_tables()
    mjd = np.asarray(mjd, np.float64)
    if tab is None:
        return np.zeros(mjd.size, bool), np.zeros(0, int)
    pos = np.clip(np.searchsorted(tab['obstime_mjd'], mjd), 0, tab['obstime_mjd'].size - 1)
    found = np.abs(tab['obstime_mjd'][pos] - mjd) < 1e-9
    for role in ROLES:
        dalt = np.abs(tab[f'{role}_alt'][pos] - np.asarray(alt_by_role[role]))
        daz = np.abs((tab[f'{role}_az'][pos] - np.asarray(az_by_role[role]) + 180.0) % 360.0 - 180.0)
        found &= (dalt < tol_deg) & (daz < tol_deg)
    return found, pos[found]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--oi-table', help='oi_table.npz with ctx_{sci,near,far}, expnum, obstime_mjd')
    ap.add_argument('--meta-fits', help='corpus *_meta_only.fits: write its bundled context table to --out')
    ap.add_argument('--out', required=True)
    ap.add_argument('--height-km', type=float, default=IPP_HEIGHT_KM, help='pierce-point height for VTEC sampling')
    args = ap.parse_args()
    if args.meta_fits:
        table = build_context_table(args.meta_fits, args.out)
        print(f'wrote {args.out}: {table["obstime_mjd"].size} rows')
        return
    d = np.load(args.oi_table)
    mjd = np.asarray(d['obstime_mjd'], np.float64)
    table = build_exposure_table(d['expnum'], mjd, _pointings_from_oi_table(d), args.height_km)
    np.savez(args.out, **table)
    print(f'wrote {args.out}: {mjd.size} exposures, sources '
          f'{dict(zip(*np.unique(table["sci_ionex_source"], return_counts=True), strict=True))}')


if __name__ == '__main__':
    main()
