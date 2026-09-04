"""Data loading + geometry-feature pipeline.

Consolidates the following cells from the split-zodi notebook:
- ``f4b17802``        : decomposition + META loader, geometry feature builders,
                        solar-activity table, airglow scale, effective extinction.
- ``20eb0763``        : triplet builder (near / far -> sci).
- ``ed6334b4``        : reconstruction helpers (`reconstruct_with_lsf`).
- ``cd99f8df``        : row-level triplet filters (kappa-sigma, kappa-MAD,
                        LMC/SMC exclusion).
- ``ecliptic-ctx-augment``     : per-arm ecliptic coordinate augmentation.
- ``physics-priors-augment``   : Leinert zodi V-band prior features.

All functions preserve the notebook's public signatures.  A small
``DataPipeline`` class at the bottom composes them into a single
``load_and_prepare()`` entry point for notebook callers.
"""

from __future__ import annotations

import io
import os
import re
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from astropy.coordinates import (
    AltAz,
    BarycentricMeanEcliptic,
    EarthLocation,
    SkyCoord,
    get_sun,
)
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sky_decomp.fit import reconstruct_component_spectra
from sky_decomp.lsf_surface_iterative import (
    LSFSurfaceState,
    SkyDecompLSFSurfaceIterative,
    load_lsf_surface_state,
)

# Reused across `airglow_van_rhijn_matrix` and per-arm geometry augmentation.
LMC_EXCLUSION = {"ra_center_deg": 80.9, "dec_center_deg": -69.75, "radius_deg": 8.0}
SMC_EXCLUSION = {"ra_center_deg": 13.2, "dec_center_deg": -72.83, "radius_deg": 5.0}


# ==========================================================================
# Extracted from notebook cell id=f4b17802
# ==========================================================================

# Reused decomposition/context loading helpers
LAYER_HEIGHTS_KM = {
    'mesospheric_oh': 87.0,
    'mesospheric_atomic': 95.0,
    'ionospheric_red': 285.0,
}

TIME_FEATURE_NAMES = {
    'obstime_year_sin',
    'obstime_year_cos',
    'obstime_day_sin',
    'obstime_day_cos',
    'obstime_lunation_sin',
    'obstime_lunation_cos',
}

VAN_RHIJN_FEATURES = {
    'vanrhijn_87km': LAYER_HEIGHTS_KM['mesospheric_oh'],
    'vanrhijn_95km': LAYER_HEIGHTS_KM['mesospheric_atomic'],
    'vanrhijn_285km': LAYER_HEIGHTS_KM['ionospheric_red'],
}


def _as_array(x):
    arr = np.asarray(x)
    if arr.dtype.kind in ('U', 'S', 'O'):
        return None
    return arr.astype(np.float32)


def _coerce_coef_hdu_to_table(coef_hdu):
    data = coef_hdu.data
    if isinstance(coef_hdu, (fits.BinTableHDU, fits.TableHDU)):
        return Table(data)

    arr = np.asarray(data, dtype=np.float32)
    if arr.ndim != 2:
        raise ValueError(f'Expected 2D COEF image, got shape={arr.shape}')

    n_coef = arr.shape[1]
    names = []
    for i in range(n_coef):
        key = f'COEF{i:04d}'
        names.append(str(coef_hdu.header.get(key, f'coef_{i:04d}')))
    return Table({name: arr[:, i] for i, name in enumerate(names)})


def _select_context_from_labels(meta, meta_upper, labels, base_name):
    e_key = f'SKYE_{base_name.upper()}'
    w_key = f'SKYW_{base_name.upper()}'
    if e_key not in meta_upper or w_key not in meta_upper:
        return None

    arr_e = _as_array(meta[meta_upper[e_key]])
    arr_w = _as_array(meta[meta_upper[w_key]])
    if arr_e is None or arr_w is None:
        raise ValueError(f'Labeled context columns for {base_name} are non-numeric.')

    is_e = labels == 'SKYE'
    is_w = labels == 'SKYW'
    if not np.all(is_e | is_w):
        bad = np.unique(labels[~(is_e | is_w)])
        raise ValueError(f'Unexpected label values: {bad}')

    return np.where(is_e, arr_e, arr_w).astype(np.float32)


def _table_to_float32_matrix(tbl, value_name):
    names = list(tbl.colnames)
    cols = []
    numeric_names = []
    for name in names:
        arr = _as_array(tbl[name])
        if arr is not None:
            cols.append(arr)
            numeric_names.append(name)

    if len(cols) == 0:
        raise ValueError(f'No numeric {value_name} columns found.')

    return np.column_stack(cols).astype(np.float32), numeric_names


def _altitude_to_zenith_deg(alt_deg):
    alt = np.asarray(alt_deg, dtype=np.float64)
    return np.clip(90.0 - alt, 0.0, 89.9)


def _van_rhijn_factor(alt_deg, layer_height_km, earth_radius_km=6371.0):
    z_rad = np.deg2rad(_altitude_to_zenith_deg(alt_deg))
    shell_ratio = float(earth_radius_km) / (float(earth_radius_km) + float(layer_height_km))
    denom = 1.0 - (shell_ratio * np.sin(z_rad)) ** 2
    denom = np.clip(denom, 1e-6, None)
    return (1.0 / np.sqrt(denom)).astype(np.float32)


def _extract_obstime_mjd(meta, meta_upper):
    from astropy.time import Time

    candidates = ('OBSTIME', 'MJD', 'MJD_OBS', 'MJD-OBS', 'DATE_OBS', 'DATE-OBS')
    for key in candidates:
        if key not in meta_upper:
            continue

        raw = np.asarray(meta[meta_upper[key]])
        if raw.dtype.kind in ('i', 'u', 'f'):
            return raw.astype(np.float64)

        s = np.asarray(raw).astype(str)
        try:
            return Time(s, format='isot', scale='utc').mjd.astype(np.float64)
        except ValueError:
            return Time(s).mjd.astype(np.float64)

    raise KeyError('Could not find an OBSTIME/MJD/DATE-OBS-like META column for time features')


def _build_obstime_feature(meta, meta_upper, feature_name):
    mjd = _extract_obstime_mjd(meta, meta_upper)
    two_pi = 2.0 * np.pi

    if feature_name == 'obstime_year_sin':
        return np.sin(two_pi * mjd / 365.2422).astype(np.float32)

    if feature_name == 'obstime_year_cos':
        return np.cos(two_pi * mjd / 365.2422).astype(np.float32)

    frac_day = mjd - np.floor(mjd)
    if feature_name == 'obstime_day_sin':
        return np.sin(two_pi * frac_day).astype(np.float32)

    if feature_name == 'obstime_day_cos':
        return np.cos(two_pi * frac_day).astype(np.float32)

    if feature_name == 'obstime_lunation_sin':
        return np.sin(two_pi * mjd / 29.53058867).astype(np.float32)

    if feature_name == 'obstime_lunation_cos':
        return np.cos(two_pi * mjd / 29.53058867).astype(np.float32)

    raise KeyError(f'Unsupported obstime feature name: {feature_name}')


# --- Astropy-computed pointing geometry ------------------------------------
# alt / az / airmass and sun / moon positions are computed from RA, DEC and
# obstime via astropy rather than read directly from META. Centralising the
# geometry here means a per-column error in the input META (which happens for
# online-derived quantities) does not silently propagate through the pipeline.
# `moon_phase` is intentionally not in this set: illumination fraction is not
# a pure geometry feature and is still read from META.
ASTROPY_GEOMETRY_FEATURES = {
    'alt', 'az_sin', 'az_cos', 'airmass',
    'moon_alt', 'moon_az_sin', 'moon_az_cos', 'moon_sep',
    'sun_alt', 'sun_az_sin', 'sun_az_cos', 'sun_sep',
    'sci_sep',
}

# META columns stored as an angle in degrees that wraps at 0/360.
# `<stem>_sin` / `<stem>_cos` in context_cols pull the raw column and fold it.
CYCLIC_META_DEGREE_FEATURES = {'moon_phase'}

_LCO_EARTHLOCATION = None


def _lco_earth_location():
    """LCO EarthLocation, cached at first use."""
    global _LCO_EARTHLOCATION
    if _LCO_EARTHLOCATION is not None:
        return _LCO_EARTHLOCATION
    from astropy.coordinates import EarthLocation
    import astropy.units as u
    try:
        _LCO_EARTHLOCATION = EarthLocation.of_site('Las Campanas Observatory')
    except Exception:
        # Fallback to published LCO coordinates (Baade / du Pont site).
        _LCO_EARTHLOCATION = EarthLocation(
            lat=-29.00889 * u.deg,
            lon=-70.68992 * u.deg,
            height=2380.0 * u.m,
        )
    return _LCO_EARTHLOCATION


def _pointing_ra_dec_columns(meta_upper, kind):
    """Return the META column names for the pointing RA / DEC of one kind.

    Uses SKY_NEAR_RA / SKY_NEAR_DEC and SKY_FAR_RA / SKY_FAR_DEC for the two
    sky arms, so the near/far ordering (by separation from science) and the
    east/west assignment are picked up from the same columns that already
    drive the rest of the pipeline; no consultation of SKY_NEAR_LABEL /
    SKY_FAR_LABEL is needed for the pointing geometry.
    """
    if kind == 'sci':
        ra_key = next((meta_upper[k] for k in ('SCI_RA', 'RA') if k in meta_upper), None)
        dec_key = next((meta_upper[k] for k in ('SCI_DEC', 'DEC') if k in meta_upper), None)
    elif kind == 'sky1':
        ra_key = meta_upper.get('SKY_NEAR_RA')
        dec_key = meta_upper.get('SKY_NEAR_DEC')
    elif kind == 'sky2':
        ra_key = meta_upper.get('SKY_FAR_RA')
        dec_key = meta_upper.get('SKY_FAR_DEC')
    else:
        raise ValueError(f'unexpected kind {kind!r}')
    if ra_key is None or dec_key is None:
        raise KeyError(f'RA / DEC columns for kind={kind!r} not found in META')
    return ra_key, dec_key


def _compute_astropy_geometry(meta, meta_upper, kind):
    """Compute all astropy-derived geometry features for one pointing kind.

    Returns a dict keyed by feature name (subset of ASTROPY_GEOMETRY_FEATURES)
    with float32 arrays of length n_rows. Pointing RA / DEC come from the
    kind-specific columns; obstime comes from `_extract_obstime_mjd`. Sun and
    Moon topocentric positions are computed at LCO for every obstime. Airmass
    is sec(z) = 1 / sin(alt), NaN below the horizon.
    """
    from astropy.coordinates import AltAz, SkyCoord, get_body, get_sun
    from astropy.time import Time
    import astropy.units as u

    ra_key, dec_key = _pointing_ra_dec_columns(meta_upper, kind)
    ra = np.asarray(meta[ra_key], dtype=np.float64)
    dec = np.asarray(meta[dec_key], dtype=np.float64)

    # Science pointing coords are needed for the sky-to-sci angular separation feature.
    sci_ra_key, sci_dec_key = _pointing_ra_dec_columns(meta_upper, 'sci')
    sci_ra_col = np.asarray(meta[sci_ra_key], dtype=np.float64)
    sci_dec_col = np.asarray(meta[sci_dec_key], dtype=np.float64)

    mjd = _extract_obstime_mjd(meta, meta_upper)
    time = Time(mjd, format='mjd', scale='utc')

    lco = _lco_earth_location()
    frame = AltAz(obstime=time, location=lco)

    pointing = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    pointing_altaz = pointing.transform_to(frame)
    sci_pointing = SkyCoord(ra=sci_ra_col * u.deg, dec=sci_dec_col * u.deg, frame='icrs')
    sci_sep = pointing.separation(sci_pointing).to_value(u.deg)

    sun = get_sun(time)
    moon = get_body('moon', time, location=lco)
    sun_altaz = sun.transform_to(frame)
    moon_altaz = moon.transform_to(frame)

    sun_sep = pointing.separation(sun).to_value(u.deg)
    moon_sep = pointing.separation(moon).to_value(u.deg)

    alt_deg = pointing_altaz.alt.to_value(u.deg)
    sinalt = np.sin(np.deg2rad(alt_deg))
    airmass = np.where(sinalt > 1e-6, 1.0 / np.clip(sinalt, 1e-6, None), np.nan)

    az_rad = np.deg2rad(pointing_altaz.az.to_value(u.deg))
    moon_az_rad = np.deg2rad(moon_altaz.az.to_value(u.deg))
    sun_az_rad = np.deg2rad(sun_altaz.az.to_value(u.deg))

    return {
        'alt': alt_deg.astype(np.float32),
        'az_sin': np.sin(az_rad).astype(np.float32),
        'az_cos': np.cos(az_rad).astype(np.float32),
        'airmass': airmass.astype(np.float32),
        'moon_alt': moon_altaz.alt.to_value(u.deg).astype(np.float32),
        'moon_az_sin': np.sin(moon_az_rad).astype(np.float32),
        'moon_az_cos': np.cos(moon_az_rad).astype(np.float32),
        'moon_sep': moon_sep.astype(np.float32),
        'sun_alt': sun_altaz.alt.to_value(u.deg).astype(np.float32),
        'sun_az_sin': np.sin(sun_az_rad).astype(np.float32),
        'sun_az_cos': np.cos(sun_az_rad).astype(np.float32),
        'sun_sep': sun_sep.astype(np.float32),
        'sci_sep': sci_sep.astype(np.float32),
    }


# ---------------------------------------------------------------------------
# Solar activity proxies (F10.7 daily flux and Kp geomagnetic index).
#
# F10.7 is the standard EUV proxy driving thermospheric density and OI 6300 /
# OI 6364 recombination amplitude; Kp captures short-timescale auroral /
# geomagnetic driving that also modulates the 285 km group's brightness.
# Both are fetched from GFZ Potsdam's authoritative archive
#     https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt
# (definitive + quicklook, back to 1932, updated daily) and cached under
# solar_activity_cache/.  The raw file is re-downloaded if older than
# SOLAR_ACTIVITY_MAX_AGE_DAYS or missing.
# ---------------------------------------------------------------------------
SOLAR_ACTIVITY_FEATURES = {'f107', 'f107_81d', 'kp'}
SOLAR_ACTIVITY_SOURCE_URL = 'https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt'
SOLAR_ACTIVITY_CACHE_DIR = Path('solar_activity_cache')
SOLAR_ACTIVITY_RAW_PATH = SOLAR_ACTIVITY_CACHE_DIR / 'Kp_ap_Ap_SN_F107_since_1932.txt'
SOLAR_ACTIVITY_NPZ_PATH = SOLAR_ACTIVITY_CACHE_DIR / 'kp_f107.npz'
SOLAR_ACTIVITY_MAX_AGE_DAYS = 7

_SOLAR_ACTIVITY_TABLE = None  # module-level cache of parsed arrays


def _download_gfz_solar_activity(url, dest, timeout=60, verbose=True):
    """Fetch the GFZ Kp / ap / F10.7 archive to `dest` via a plain HTTP GET."""
    from urllib.request import Request, urlopen
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + '.part')
    req = Request(url, headers={'User-Agent': 'lvmsky-notebook/1.0'})
    with urlopen(req, timeout=timeout) as resp, open(tmp, 'wb') as fh:
        fh.write(resp.read())
    tmp.replace(dest)
    if verbose:
        print(f'  solar-activity archive downloaded: {dest} '
              f'({dest.stat().st_size / 1e6:.2f} MB)')


def _parse_gfz_solar_activity_raw(path):
    """Parse the GFZ text file. Returns dict of arrays: kp (3-hourly) + F10.7 (daily)."""
    df = pd.read_csv(path, comment='#', sep=r'\s+', header=None, engine='python')
    if df.shape[1] < 26:
        raise RuntimeError(
            f'Unexpected GFZ Kp/F10.7 file layout: got {df.shape[1]} columns, '
            f'expected at least 26. First data row: {df.iloc[0].tolist()!r}.'
        )
    year = df.iloc[:, 0].astype(int).to_numpy()
    month = df.iloc[:, 1].astype(int).to_numpy()
    day = df.iloc[:, 2].astype(int).to_numpy()
    # columns 3-6 = bookkeeping (days since 1932-01-01, Bartels rotation number, day of rotation)
    kp_cols = df.iloc[:, 7:15].astype(float).to_numpy()   # 8 three-hourly Kp values (Kp1..Kp8)
    # columns 15-22 = 8 ap values; 23 = Ap daily; 24 = SN; 25 = F10.7 obs; 26 = F10.7 adj
    f107_obs = df.iloc[:, 25].astype(float).to_numpy()

    from astropy.time import Time
    iso = [f'{y:04d}-{m:02d}-{d:02d}' for y, m, d in zip(year, month, day)]
    mjd_day0 = Time(iso, format='iso', scale='utc').mjd.astype(np.float64)
    # 3-hour Kp windows start at 0h, 3h, 6h, ..., 21h UT.
    hour_offsets = np.arange(8, dtype=np.float64) * (3.0 / 24.0)
    kp_mjd_start = (mjd_day0[:, None] + hour_offsets[None, :]).reshape(-1)
    kp_value = kp_cols.reshape(-1)
    # Daily F10.7 attached to noon UT so np.interp reads it as a midday sample.
    f107_mjd = mjd_day0 + 0.5
    good_kp = np.isfinite(kp_value) & (kp_value >= 0.0)
    good_f107 = np.isfinite(f107_obs) & (f107_obs > 0.0)
    return {
        'kp_mjd_start': kp_mjd_start[good_kp],
        'kp_value': kp_value[good_kp],
        'f107_mjd': f107_mjd[good_f107],
        'f107_obs': f107_obs[good_f107],
    }


def _f107_running_mean(mjd, f107, window_days=81.0):
    """Symmetric running mean of daily F10.7 over `window_days` (~3 solar rotations)."""
    order = np.argsort(mjd)
    mjd_s = mjd[order]
    f107_s = f107[order].astype(np.float64)
    n = mjd_s.size
    half = float(window_days) / 2.0
    lo = np.searchsorted(mjd_s, mjd_s - half, side='left')
    hi = np.searchsorted(mjd_s, mjd_s + half, side='right')
    csum = np.concatenate([[0.0], np.cumsum(f107_s)])
    sums = csum[hi] - csum[lo]
    counts = (hi - lo).astype(np.float64)
    with np.errstate(invalid='ignore', divide='ignore'):
        out_s = np.where(counts > 0, sums / counts, np.nan)
    inv = np.empty(n, dtype=int)
    inv[order] = np.arange(n)
    return out_s[inv]


def load_solar_activity_table(
    url=SOLAR_ACTIVITY_SOURCE_URL,
    raw_path=SOLAR_ACTIVITY_RAW_PATH,
    npz_path=SOLAR_ACTIVITY_NPZ_PATH,
    max_age_days=SOLAR_ACTIVITY_MAX_AGE_DAYS,
    force_refresh=False,
    verbose=True,
):
    """Return sorted numpy arrays for 3-hourly Kp and daily F10.7 (obs + 81-day mean).

    Downloads the GFZ archive if the raw cache is missing or older than
    `max_age_days` days.  Reuses a stale cache silently if the download fails.
    """
    global _SOLAR_ACTIVITY_TABLE
    if _SOLAR_ACTIVITY_TABLE is not None and not force_refresh:
        return _SOLAR_ACTIVITY_TABLE

    import time as _time_mod
    need_download = force_refresh or (not raw_path.exists())
    if not need_download:
        age_days = (_time_mod.time() - raw_path.stat().st_mtime) / 86400.0
        if age_days > float(max_age_days):
            if verbose:
                print(f'Solar-activity cache is {age_days:.1f} days old '
                      f'(> {max_age_days} d); refreshing.')
            need_download = True
    if need_download:
        try:
            _download_gfz_solar_activity(url, raw_path, verbose=verbose)
        except Exception as exc:
            if raw_path.exists():
                print(f'Solar-activity download failed '
                      f'({type(exc).__name__}: {exc}); reusing existing cache at {raw_path}.')
            else:
                raise RuntimeError(
                    f'Could not download solar-activity archive from {url} '
                    f'({type(exc).__name__}: {exc}); no local cache. '
                    f'Provide the file manually at {raw_path} to proceed.') from exc

    parsed = _parse_gfz_solar_activity_raw(raw_path)
    kp_order = np.argsort(parsed['kp_mjd_start'])
    kp_mjd_start = parsed['kp_mjd_start'][kp_order].astype(np.float64)
    kp_value = parsed['kp_value'][kp_order].astype(np.float32)
    f107_order = np.argsort(parsed['f107_mjd'])
    f107_mjd = parsed['f107_mjd'][f107_order].astype(np.float64)
    f107_obs = parsed['f107_obs'][f107_order].astype(np.float32)
    f107_81d = _f107_running_mean(f107_mjd, f107_obs, window_days=81.0).astype(np.float32)

    table = {
        'kp_mjd_start': kp_mjd_start,
        'kp_value': kp_value,
        'f107_mjd': f107_mjd,
        'f107_obs': f107_obs,
        'f107_81d': f107_81d,
    }
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(npz_path, **table)
    _SOLAR_ACTIVITY_TABLE = table
    if verbose:
        print(f'Solar-activity table ready: '
              f'{kp_mjd_start.size} Kp windows '
              f'(MJD {kp_mjd_start.min():.1f} .. {kp_mjd_start.max():.1f}), '
              f'{f107_mjd.size} daily F10.7 samples; cached to {npz_path}.')
    return table


def _solar_activity_lookup(mjd_query, feature):
    """Evaluate one solar-activity feature at each observation MJD."""
    table = load_solar_activity_table()
    mjd_query = np.asarray(mjd_query, dtype=np.float64)
    if feature == 'kp':
        # nearest 3-hour window that starts at or before the observation
        idx = np.searchsorted(table['kp_mjd_start'], mjd_query, side='right') - 1
        n_kp = table['kp_value'].size
        oob = (idx < 0) | (idx >= n_kp)
        idx = np.clip(idx, 0, n_kp - 1)
        val = table['kp_value'][idx].astype(np.float32).copy()
        val[oob] = np.nan
        return val
    if feature in ('f107', 'f107_81d'):
        src = table['f107_obs'] if feature == 'f107' else table['f107_81d']
        val = np.interp(mjd_query, table['f107_mjd'], src.astype(np.float64),
                        left=np.nan, right=np.nan)
        return val.astype(np.float32)
    raise KeyError(f'unknown solar-activity feature: {feature!r}')


def _decode_cyclic_context(ctx_names, ctx_matrix):
    """Fold `<name>_sin`/`<name>_cos` pairs back into a single `<name>` in degrees [0, 360).

    Purely for display: the model still consumes the sin/cos pair, but plots
    (histograms, coefficient-vs-context relationship panels) are much more
    legible with a single 0-360 axis than with a pair of [-1, 1] projections.
    """
    names = list(ctx_names)
    mat = np.asarray(ctx_matrix, dtype=np.float64)
    out_names = []
    out_cols = []
    handled = set()
    for i, name in enumerate(names):
        if name in handled:
            continue
        if name.endswith('_sin'):
            stem = name[:-4]
            cos_name = stem + '_cos'
            if cos_name in names:
                cos_idx = names.index(cos_name)
                rad = np.arctan2(mat[:, i], mat[:, cos_idx])
                out_names.append(stem)
                out_cols.append((np.rad2deg(rad) % 360.0).astype(np.float32))
                handled.add(name)
                handled.add(cos_name)
                continue
        if name.endswith('_cos') and name in handled:
            continue
        out_names.append(name)
        out_cols.append(mat[:, i].astype(np.float32))
    return out_names, (np.column_stack(out_cols).astype(np.float32)
                       if out_cols else np.empty((mat.shape[0], 0), dtype=np.float32))


def _resolve_base_context_array(meta, meta_upper, labels, kind, base_name):
    key = base_name.upper()

    if kind == 'sci':
        sci_key = f'SCI_{key}'
        if sci_key in meta_upper:
            arr = _as_array(meta[meta_upper[sci_key]])
            if arr is None:
                raise ValueError(f'Context column {sci_key} is non-numeric.')
            return arr

    if key in meta_upper:
        arr = _as_array(meta[meta_upper[key]])
        if arr is None:
            raise ValueError(f'Context column {base_name} is non-numeric.')
        return arr

    if labels is not None:
        arr = _select_context_from_labels(meta, meta_upper, labels, base_name)
        if arr is not None:
            return arr

    return None


def _resolve_context_feature(meta, meta_upper, labels, kind, feature_name, cache):
    if feature_name in cache:
        return cache[feature_name]

    if feature_name in ('obstime_mjd', 'obstime_mjd_z'):
        raise ValueError('Direct time features are excluded from model context; use only periodic year/day/lunation sine and cosine terms.')

    if feature_name in ASTROPY_GEOMETRY_FEATURES:
        # Compute all astropy features in one pass the first time any is asked for.
        if '_astropy_computed' not in cache:
            for _name, _arr in _compute_astropy_geometry(meta, meta_upper, kind).items():
                cache[_name] = np.asarray(_arr, dtype=np.float32)
            cache['_astropy_computed'] = True
        return cache[feature_name]

    if feature_name in SOLAR_ACTIVITY_FEATURES:
        if '_solar_activity_computed' not in cache:
            _sa_mjd = _extract_obstime_mjd(meta, meta_upper)
            for _sa_name in SOLAR_ACTIVITY_FEATURES:
                cache[_sa_name] = np.asarray(
                    _solar_activity_lookup(_sa_mjd, _sa_name), dtype=np.float32)
            cache['_solar_activity_computed'] = True
        return cache[feature_name]

    if feature_name == 'ew':
        if kind == 'sci':
            arr = np.zeros(len(meta), dtype=np.float32)
        else:
            arr = np.where(labels == 'SKYE', 1.0,
                           np.where(labels == 'SKYW', -1.0, 0.0)).astype(np.float32)
        cache[feature_name] = arr
        return arr

    if feature_name in TIME_FEATURE_NAMES:
        arr = _build_obstime_feature(meta, meta_upper, feature_name)
    elif feature_name == 'zenith_deg':
        alt = _resolve_context_feature(meta, meta_upper, labels, kind, 'alt', cache)
        arr = _altitude_to_zenith_deg(alt).astype(np.float32)
    elif (feature_name.endswith('_sin') or feature_name.endswith('_cos')) \
            and feature_name[:-4] in CYCLIC_META_DEGREE_FEATURES:
        raw = _resolve_base_context_array(meta, meta_upper, labels, kind, feature_name[:-4])
        if raw is None:
            raise KeyError(feature_name)
        rad = np.deg2rad(np.asarray(raw, dtype=np.float64))
        arr = (np.sin(rad) if feature_name.endswith('_sin') else np.cos(rad)).astype(np.float32)
    elif feature_name in VAN_RHIJN_FEATURES:
        alt = _resolve_context_feature(meta, meta_upper, labels, kind, 'alt', cache)
        arr = _van_rhijn_factor(alt, VAN_RHIJN_FEATURES[feature_name])
    else:
        arr = _resolve_base_context_array(meta, meta_upper, labels, kind, feature_name)
        if arr is None:
            raise KeyError(feature_name)

    cache[feature_name] = np.asarray(arr, dtype=np.float32)
    return cache[feature_name]


def _build_context_matrix(meta, context_columns, kind):
    meta_upper = {c.upper(): c for c in meta.colnames}
    labels = None

    if kind in ('sky1', 'sky2'):
        label_col = 'SKY_NEAR_LABEL' if kind == 'sky1' else 'SKY_FAR_LABEL'
        if label_col not in meta_upper:
            raise KeyError(f'Missing required META label column: {label_col}')
        labels = np.char.upper(np.char.strip(np.asarray(meta[meta_upper[label_col]]).astype(str)))

    ctx_names = []
    ctx_cols = []
    missing_cols = []
    cache = {}

    for cname in context_columns:
        try:
            arr = _resolve_context_feature(meta, meta_upper, labels, kind, cname, cache)
        except KeyError:
            missing_cols.append(cname)
            continue
        ctx_names.append(cname)
        ctx_cols.append(arr)

    if missing_cols:
        raise KeyError(f'Missing requested context columns: {missing_cols}')
    if len(ctx_cols) == 0:
        raise ValueError('No usable context columns were assembled.')

    return np.column_stack(ctx_cols).astype(np.float32), ctx_names


def _find_chi2_column(meta_tbl):
    names = {c.upper(): c for c in meta_tbl.colnames}
    for cand in ('REDUCED_CHI2', 'CHI2_REDUCED', 'CHI2', 'RCHI2'):
        if cand in names:
            return names[cand]
    raise KeyError('No chi2-like column found in decomposition META table')


def read_decomp_dataset(decomp_fits_path, input_fits_path, context_columns,
                        decomp_kind='sky1', return_chi2=False, return_err=False):
    if context_columns is None or len(context_columns) == 0:
        raise ValueError('context_columns must be a non-empty list.')

    kind = decomp_kind.lower()
    if kind not in ('sky1', 'sky2', 'sci'):
        raise ValueError("decomp_kind must be one of: 'sky1', 'sky2', 'sci'")

    with fits.open(decomp_fits_path) as hdul_dec, fits.open(input_fits_path) as hdul_in:
        coef_tbl = _coerce_coef_hdu_to_table(hdul_dec['COEF'])
        coef_mat, coef_names = _table_to_float32_matrix(coef_tbl, 'coefficient')

        # COEF_ERR is the active-set 1-sigma uncertainty on each coefficient
        # (see fit.SkyDecomp._coef_err_active_set); older FITS files omit it,
        # so we fall back to NaN and callers can decide whether to weight.
        coef_err_mat = None
        if return_err:
            if 'COEF_ERR' in hdul_dec:
                coef_err_tbl = _coerce_coef_hdu_to_table(hdul_dec['COEF_ERR'])
                coef_err_mat, coef_err_names = _table_to_float32_matrix(coef_err_tbl, 'coefficient error')
                if coef_err_names != coef_names:
                    raise ValueError('COEF_ERR column ordering does not match COEF')
            else:
                coef_err_mat = np.full_like(coef_mat, np.nan, dtype=np.float32)

        meta = Table(hdul_in['META'].data)
        ctx_mat, ctx_names = _build_context_matrix(meta, context_columns, kind)

        if coef_mat.shape[0] != ctx_mat.shape[0]:
            raise ValueError(
                f'Row count mismatch: COEF has {coef_mat.shape[0]} rows, META has {ctx_mat.shape[0]} rows'
            )

        good = np.isfinite(coef_mat).all(axis=1) & np.isfinite(ctx_mat).all(axis=1)
        coef_mat = coef_mat[good]
        ctx_mat = ctx_mat[good]
        if coef_err_mat is not None:
            coef_err_mat = coef_err_mat[good]

        if not return_chi2 and not return_err:
            return coef_mat, ctx_mat, coef_names, ctx_names

        chi2_used = None
        if return_chi2:
            dec_meta = Table(hdul_dec['META'].data)
            chi2_col = _find_chi2_column(dec_meta)
            chi2_full = np.asarray(dec_meta[chi2_col], dtype=np.float64)
            chi2_used = chi2_full[good]
            if chi2_used.shape[0] != coef_mat.shape[0]:
                raise ValueError(
                    f'Aligned chi2 rows ({chi2_used.shape[0]}) do not match coef rows ({coef_mat.shape[0]})'
                )

        result = [coef_mat, ctx_mat, coef_names, ctx_names]
        if return_chi2:
            result.append(chi2_used)
        if return_err:
            result.append(coef_err_mat)
        return tuple(result)


# ---------------------------------------------------------------------------
# Coefficient grouping by physical origin / emission layer.
# Single definition, used by the wavelength resolver, the structure-function
# diagnostics and the model.  Coefficients are consumed in order so each one
# lands in exactly one group.
# ---------------------------------------------------------------------------

def _contains_any(name, needles):
    return any(token in name for token in needles)


# Coefficient-name schema. One-to-one with the design_names produced by
# SkyDecomp._build_static_basis() in skysub/sky_decomp/fit.py. Every entry
# is (regex, group, wavelength_a, description). Any coefficient that does
# not match a pattern raises in _build_group_indices() -- silently
# absorbing unknown names into an 'other' group previously hid that
# HO2/FeO/O2Ac were routed into 'mesospheric' and that OI 6300/6364 +
# OI 7774/8446 (ATOM_Or, ATOM_Orc_*) were routed into 'atomic' (~95 km)
# instead of 'ionospheric' (~285 km).
#
# Group layers (see GROUP_HEIGHT_FEATURE):
#   moon         no thin-shell height (scattered moonlight + spline)
#   continuum    ~87 km  HO2 / FeO / O2Ac mesopause chemiluminescence
#   mesospheric  ~87 km  OH Meinel bands + O2 atmospheric band
#   atomic       ~95 km  K, Na, N I, OI 5577 (mesopause metal / metastable)
#   ionospheric  ~285 km OI 6300/6364 red + OI 7774/8446 F-region recomb
COEF_SCHEMA = [
    (re.compile(r'^OH_\d{3}$'),         'mesospheric', None,   'OH Meinel band group (wavelength from basis)'),
    (re.compile(r'^O2_b\d+$'),          'mesospheric', 8645.0, 'O2 (b1Sigma) atmospheric band'),
    (re.compile(r'^Moon_bs\d+$'),       'moon',        None,   'Moon scattered-light spline knot'),
    (re.compile(r'^Zodi_bs\d+$'),       'zodi',        None,   'Zodiacal-light color spline knot'),
    (re.compile(r'^HO2$'),              'continuum',   None,   'HO2 diffuse continuum'),
    (re.compile(r'^FeO$'),              'continuum',   None,   'FeO orange-arc diffuse continuum'),
    (re.compile(r'^O2Ac$'),             'continuum',   None,   'O2 Chamberlain diffuse continuum'),
    (re.compile(r'^ATOM_K$'),           'atomic',      7698.96, 'K I 7699'),
    (re.compile(r'^ATOM_N$'),           'ionospheric', 5199.0,  'N I 5199 [NI] 2D->4S F-region metastable'),
    (re.compile(r'^ATOM_Na$'),          'atomic',      5892.9,  'Na D'),
    (re.compile(r'^ATOM_Og$'),          'atomic',      5577.34, 'OI 5577 green line'),
    (re.compile(r'^ATOM_Or$'),          'ionospheric', 6300.3,  'OI 6300/6364 red doublet'),
    (re.compile(r'^ATOM_Orc_OI0777$'),  'ionospheric', 7774.2,  'OI 7774 F-region recombination'),
    (re.compile(r'^ATOM_Orc_OI0845$'),  'ionospheric', 8446.4,  'OI 8446 F-region recombination'),
]

# Canonical output order for group dicts.
_COEF_GROUP_ORDER = ('moon', 'zodi', 'continuum', 'mesospheric', 'ionospheric', 'atomic')


def _lookup_coef_schema(name):
    """Return (group, wavelength_a, description) for one coefficient, else None."""
    for pat, group, lam, desc in COEF_SCHEMA:
        if pat.match(str(name)):
            return group, lam, desc
    return None


def _build_group_indices(coef_names):
    """Route each coefficient name to exactly one group via COEF_SCHEMA.

    Raises RuntimeError on any unmatched name so a new coefficient family
    added to fit.py becomes a loud error here rather than a silent 'other'.
    """
    groups = {g: [] for g in _COEF_GROUP_ORDER}
    unmatched = []
    for j, name in enumerate(coef_names):
        hit = _lookup_coef_schema(name)
        if hit is None:
            unmatched.append((j, str(name)))
        else:
            groups[hit[0]].append(j)
    if unmatched:
        raise RuntimeError(
            'Unrecognised coefficient names (no COEF_SCHEMA match): '
            f'{unmatched[:8]}{"..." if len(unmatched) > 8 else ""}. '
            'Extend COEF_SCHEMA in the loaders cell and cross-check '
            'against SkyDecomp._build_static_basis() in '
            'skysub/sky_decomp/fit.py.'
        )
    return {g: np.asarray(idx, dtype=int) for g, idx in groups.items() if idx}


# Self-test: verify the schema on a canonical name list. Runs when the
# cell is executed, so a regression in COEF_SCHEMA fails immediately.
_SCHEMA_TEST_NAMES = (
    [f'OH_{i:03d}' for i in (0, 42, 401)]
    + [f'Moon_bs{i:02d}' for i in (0, 15, 28)]
    + ['HO2', 'FeO', 'O2Ac']
    + ['ATOM_K', 'ATOM_N', 'ATOM_Na', 'ATOM_Og', 'ATOM_Or',
       'ATOM_Orc_OI0777', 'ATOM_Orc_OI0845']
    + ['O2_b01']
)
_SCHEMA_TEST_EXPECTED = {
    'moon': 3, 'continuum': 3,
    'mesospheric': 4,  # 3x OH_### + O2_b01
    'atomic': 3,       # K, Na, Og
    'ionospheric': 4,  # N, Or, Orc_OI0777, Orc_OI0845
}
_schema_test = _build_group_indices(_SCHEMA_TEST_NAMES)
_schema_test_sizes = {g: int(idx.size) for g, idx in _schema_test.items()}
if _schema_test_sizes != _SCHEMA_TEST_EXPECTED:
    raise RuntimeError(
        f'COEF_SCHEMA self-test failed: got {_schema_test_sizes}, '
        f'expected {_SCHEMA_TEST_EXPECTED}'
    )
try:
    _build_group_indices(['this_is_not_a_valid_coef'])
except RuntimeError:
    pass
else:
    raise RuntimeError('COEF_SCHEMA self-test: unknown name did not raise')

# ---------------------------------------------------------------------------
# Geometric normalisation of airglow coefficients  (PHYSICAL units)
#
# Airglow amplitudes scale as
#       A_obs = A_intrinsic * V(z; h) * 10**(-0.4 * k(lambda) * (X - 1))
# where V is the van Rhijn slant-path enhancement through a thin emitting
# shell at height h, and the power of ten is the extinction suffered by that
# emission on the way down.  Dividing an observed amplitude by this factor
# recovers the intrinsic layer emissivity, which is the only quantity
# comparable between two different lines of sight.
#
# The normalisation therefore has to happen HERE, on physical coefficients,
# before the sqrt transform and before any robust scaling.  Doing it on
# scaled or sqrt-transformed values is not equivalent: the scaler subtracts a
# median, so (sqrt(c) - m) / V is not sqrt(c / V) - m, and the sqrt means the
# correct divisor would be sqrt(V) rather than V.
# ---------------------------------------------------------------------------

# HO2 (H+O2+M chemiluminescence, ~87 km), FeO (Fe+O3 orange arc, ~85-90 km)
# and O2Ac (O+O+M -> O2 Chamberlain bands, ~90-95 km) are mesopause
# thin-shell emitters, not aerosol/Rayleigh continuum -- their sky-to-sci
# transfer needs the same van Rhijn factor as the OH bands. They are kept
# in a separate coefficient GROUP (broadband basis vs OH line groups, so
# their compression pipeline stays sqrt-identity with no PCA) but the
# geometric normalisation is now shared with mesospheric at 87 km.
AIRGLOW_GROUPS = {'mesospheric', 'atomic', 'ionospheric', 'continuum'}

GROUP_HEIGHT_FEATURE = {
    'mesospheric': 'vanrhijn_87km',
    'atomic': 'vanrhijn_95km',
    'ionospheric': 'vanrhijn_285km',
    'continuum': 'vanrhijn_87km',
}

# Only used as a last-resort fallback when a per-coefficient wavelength cannot
# be determined; see resolve_coef_wavelengths_a().
GROUP_EFFECTIVE_WAVELENGTH_A = {
    'mesospheric': 8500.0,
    'continuum': 7500.0,   # HO2/FeO/O2Ac mid-band; basis route usually resolves per-coefficient
    'atomic': 5750.0,
    'ionospheric': 6300.0,
}


def _group_height_feature_name(group_name):
    """van Rhijn context feature (hence effective layer height) for a group.

    mesospheric -> 87 km   (OH Meinel bands + O2 atmospheric band)
    continuum   -> 87 km   (HO2 / FeO / O2Ac mesopause chemiluminescence)
    atomic      -> 95 km   (K / Na / N I / OI 5577 mesopause metal / metastable)
    ionospheric -> 285 km  (OI 6300/6364 red + OI 7774/8446 F-region recomb)
    moon / other -> None (not thin-shell emission)
    """
    return GROUP_HEIGHT_FEATURE.get(group_name)


def _load_lco_extinction_curve():
    lvmcore_dir = os.environ.get('LVMCORE_DIR')
    if not lvmcore_dir:
        raise EnvironmentError('LVMCORE_DIR is not set; cannot load lco_extinction.txt')
    curve_path = Path(lvmcore_dir) / 'etc' / 'lco_extinction.txt'
    curve = np.loadtxt(curve_path, dtype=np.float64)
    if curve.ndim != 2 or curve.shape[1] < 2:
        raise ValueError(f'Unexpected extinction-curve shape: {curve.shape}')
    wave_a = np.asarray(curve[:, 0], dtype=np.float64)
    ext_k = np.asarray(curve[:, 1], dtype=np.float64)
    order = np.argsort(wave_a)
    return wave_a[order], ext_k[order]


LCO_EXTINCTION_WAVE_A, LCO_EXTINCTION_K = _load_lco_extinction_curve()


def _lco_extinction_k(wavelength_a):
    """Extinction coefficient (mag/airmass) at one or many wavelengths (Angstrom)."""
    lam = np.asarray(wavelength_a, dtype=np.float64)
    k = np.interp(
        lam,
        LCO_EXTINCTION_WAVE_A,
        LCO_EXTINCTION_K,
        left=float(LCO_EXTINCTION_K[0]),
        right=float(LCO_EXTINCTION_K[-1]),
    )
    return float(k) if np.ndim(wavelength_a) == 0 else np.asarray(k, dtype=np.float64)


# --- per-coefficient effective wavelength -----------------------------------
#
# Extinction varies strongly across the LVM range (k ~ 0.11 at 5577 A vs
# ~0.02 in the z band), so a single wavelength per group is not good enough:
# the 'mesospheric' group alone spans OI 5577, Na D and the OH/O2 bands.
#
# The primary wavelength source is coef_wavelengths_from_basis (reads the
# reconstructed basis directly). The name-based fallback below consults
# COEF_SCHEMA only, so new coefficient names cannot silently pick up a
# default guessed from digit tokens or species substrings.


def _wavelength_from_name(name):
    """Wavelength (Angstrom) from COEF_SCHEMA, else NaN."""
    hit = _lookup_coef_schema(name)
    if hit is None:
        return float('nan')
    _group, lam, _desc = hit
    return float(lam) if lam is not None and np.isfinite(lam) else float('nan')


def coef_wavelengths_from_basis(
    coef_names,
    wave,
    lsf_sigma,
    base_dir=None,
    n_spline_knots=25,
    split_zodi=False,
    n_zodi_spline_knots=3,
    palace_oh_suffix=None,
    palace_diffuse_suffix=None,
    cache_path=None,
    only_indices=None,
    verbose=True,
    return_k_eff=False,
):
    """Per-coefficient wavelength centroid and B^2-weighted effective extinction.

    Exploits linearity of the decomposition: reconstructing with coef = e_j
    returns basis component j on its own, so we read two scalars per
    coefficient in one reconstruction pass:

      1. INTENSITY-WEIGHTED CENTROID (reporting / binning):
             lambda_eff(j) = sum(lambda |f_j|) / sum(|f_j|)

      2. B^2-WEIGHTED EFFECTIVE EXTINCTION under the LCO stellar curve:
             k_eff(j) = sum(|f_j|^2 k(lambda)) / sum(|f_j|^2)
         For a broadband basis component the least-squares decomposition
         amplitude picks up the B^2-weighted average of any wavelength-dependent
         factor multiplying the true emissivity, so <k>_{B^2} is the correct
         effective extinction, not k(centroid). For narrow-line components
         |f_j| is sharply peaked and both weightings agree; the correction only
         matters for HO2 / FeO / O2Ac (continuum group) and any other broadband
         basis functions.

    Costs one reconstruction call per coefficient. Both scalars are cached to
    `cache_path` (npz keys `wavelengths_a`, `k_eff_a`); an older cache without
    `k_eff_a` is auto-invalidated and recomputed. Components that reconstruct
    to zero (outside the wavelength range) leave NaN in both arrays and fall
    back to name parsing / group defaults downstream.

    Returns `wavelengths_a` by default, or `(wavelengths_a, k_eff_a)` when
    `return_k_eff=True`.
    """
    coef_names = [str(n) for n in coef_names]
    n_coef = len(coef_names)
    wave = np.asarray(wave, dtype=np.float64)

    # The cache key must include the wavelength grid: the same coefficient
    # names evaluated on a different grid give different centroids, so keying
    # on names alone would silently return stale values after a change of input
    # product.
    grid_key = float(np.nansum(wave.astype(np.float64) * np.arange(1, wave.size + 1)))
    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.exists():
            cached = np.load(cache_path, allow_pickle=True)
            same_names = [str(x) for x in cached['coef_names']] == coef_names
            same_grid = ('grid_key' in cached.files
                         and np.isclose(float(cached['grid_key']), grid_key, rtol=1e-12))
            has_k_eff = 'k_eff_a' in cached.files
            if same_names and same_grid and has_k_eff:
                if verbose:
                    print(f'Basis wavelengths + k_eff loaded from cache: {cache_path}')
                lam_cached = np.asarray(cached['wavelengths_a'], dtype=np.float64)
                k_cached = np.asarray(cached['k_eff_a'], dtype=np.float64)
                return (lam_cached, k_cached) if return_k_eff else lam_cached
            if verbose:
                if not same_names:
                    reason = 'coefficient names'
                elif not same_grid:
                    reason = 'wavelength grid'
                else:
                    reason = 'missing k_eff_a (older cache format)'
                print(f'Cache {cache_path} does not match on {reason}; recomputing.')

    if base_dir is None:
        if '_infer_base_dir_for_reconstruction' not in globals():
            raise RuntimeError(
                'base_dir not given and _infer_base_dir_for_reconstruction() is '
                'not defined yet; run the reconstruction-helpers cell first or '
                'pass base_dir= explicitly.'
            )
        base_dir = _infer_base_dir_for_reconstruction()

    if only_indices is None:
        targets = np.arange(n_coef)
    else:
        targets = np.unique(np.asarray(only_indices, dtype=int))
        if verbose:
            print(f'Basis wavelengths: evaluating {targets.size} of {n_coef} '
                  f'coefficients (one reconstruction call each).')

    k_wave = np.asarray(_lco_extinction_k(wave), dtype=np.float64)
    lam_eff = np.full(n_coef, np.nan, dtype=np.float64)
    k_eff = np.full(n_coef, np.nan, dtype=np.float64)
    for j in targets:
        unit = np.zeros(n_coef, dtype=np.float64)
        unit[j] = 1.0
        comps = reconstruct_component_spectra(
            wave=wave,
            coef=unit,
            lsf_sigma=lsf_sigma,
            n_spline_knots=n_spline_knots,
            split_zodi=split_zodi,
            n_zodi_spline_knots=n_zodi_spline_knots,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            base_dir=base_dir,
        )
        f = np.abs(np.asarray(comps['total'], dtype=np.float64))
        tot = float(np.nansum(f))
        if np.isfinite(tot) and tot > 0.0:
            lam_eff[j] = float(np.nansum(wave * f) / tot)
            f2 = f * f
            f2_sum = float(np.nansum(f2))
            if np.isfinite(f2_sum) and f2_sum > 0.0:
                k_eff[j] = float(np.nansum(f2 * k_wave) / f2_sum)

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache_path, coef_names=np.asarray(coef_names),
                 wavelengths_a=lam_eff, k_eff_a=k_eff,
                 grid_key=np.float64(grid_key))
        if verbose:
            print(f'Basis wavelengths + k_eff cached to: {cache_path}')

    if verbose:
        n_ok = int(np.isfinite(lam_eff).sum())
        n_ok_k = int(np.isfinite(k_eff).sum())
        print(f'Basis wavelengths derived for {n_ok}/{n_coef} coefficients; '
              f'basis B^2-weighted extinction k for {n_ok_k}/{n_coef}.')
    return (lam_eff, k_eff) if return_k_eff else lam_eff


def resolve_coef_wavelengths_a(
    coef_names,
    group_indices=None,
    basis_wavelengths_a=None,
    verbose=True,
):
    """Per-coefficient effective wavelength (Angstrom) plus its provenance.

    Priority, per coefficient:
      1. basis centroid from coef_wavelengths_from_basis(), where finite --
         reads the reconstructed basis directly, so it already carries the
         PMD-driven populations (including the mesopause rotational Boltzmann
         factor for OH), the LSF, and the blend structure across overlapping
         bands. This is the definitionally-correct weighting.
      2. explicit wavelength token or species lookup baked into COEF_SCHEMA.
      3. group default from GROUP_EFFECTIVE_WAVELENGTH_A.

    Returns (wavelengths_a, sources).  The printed table is deliberately
    verbose for the airglow groups: a silently wrong wavelength turns into a
    silently wrong extinction correction.
    """
    coef_names = [str(n) for n in coef_names]
    n_coef = len(coef_names)
    lam = np.full(n_coef, np.nan, dtype=np.float64)
    source = ['unset'] * n_coef

    if basis_wavelengths_a is not None:
        basis = np.asarray(basis_wavelengths_a, dtype=np.float64)
        if basis.size != n_coef:
            raise ValueError(
                f'basis_wavelengths_a has {basis.size} entries, expected {n_coef}'
            )
        ok = np.isfinite(basis) & (basis > 0.0)
        lam[ok] = basis[ok]
        for j in np.flatnonzero(ok):
            source[j] = 'basis'

    for j in range(n_coef):
        if np.isfinite(lam[j]):
            continue
        cand = _wavelength_from_name(coef_names[j])
        if np.isfinite(cand):
            lam[j] = cand
            source[j] = 'name'

    group_of = {}
    if group_indices is not None:
        for g, idx in group_indices.items():
            for j in np.asarray(idx, dtype=int):
                group_of[int(j)] = g
        for g, idx in group_indices.items():
            default = GROUP_EFFECTIVE_WAVELENGTH_A.get(g)
            if default is None:
                continue
            for j in np.asarray(idx, dtype=int):
                if not np.isfinite(lam[int(j)]):
                    lam[int(j)] = float(default)
                    source[int(j)] = f'group:{g}'

    if verbose:
        airglow_idx = [j for j in range(n_coef)
                       if group_of.get(j) in AIRGLOW_GROUPS]
        print(f'Per-coefficient wavelengths resolved for {n_coef} coefficients '
              f'({len(airglow_idx)} in airglow groups, which are the ones that '
              f'receive a geometry correction).')
        if airglow_idx:
            print(f"  {'coefficient':<28s} {'group':<13s} {'lambda[A]':>10s} "
                  f"{'k[mag/X]':>9s}  source")
            for j in airglow_idx:
                kj = _lco_extinction_k(lam[j]) if np.isfinite(lam[j]) else np.nan
                print(f'  {coef_names[j]:<28s} {group_of.get(j, "-"):<13s} '
                      f'{lam[j]:>10.1f} {kj:>9.4f}  {source[j]}')
            from collections import Counter
            print('  wavelength provenance:',
                  dict(Counter(source[j] for j in airglow_idx)))

    return lam, source


def assert_context_is_physical(ctx, ctx_names, rtol=1e-4, atol=1e-4):
    """Fail loudly if scaler-normalised context reaches the geometry code.

    Two independent checks:
      1. physical ranges (alt within +/-90 deg, airmass >= 1)
      2. any vanrhijn_* column must agree with recomputing van Rhijn from the
         'alt' column.  Those two are computed independently upstream, so
         disagreement means the context has been transformed on the way in.

    This is the regression test for the class of bug where geometry is
    evaluated on RobustScaler output: a scaled altitude of ~0.5 becomes a
    zenith angle of ~89.5 deg and every row gets a near-horizon van Rhijn
    factor of ~6.
    """
    ctx = np.asarray(ctx, dtype=np.float64)
    names = [str(n).strip().lower() for n in ctx_names]

    if ctx.ndim != 2:
        raise ValueError(f'context must be 2-D, got shape {ctx.shape}')
    if ctx.shape[1] != len(names):
        raise ValueError(
            f'context has {ctx.shape[1]} columns but {len(names)} names given'
        )
    for required in ('alt', 'airmass'):
        if required not in names:
            raise KeyError(
                f"geometry normalisation needs '{required}' among the context columns"
            )

    alt = ctx[:, names.index('alt')]
    airmass = ctx[:, names.index('airmass')]

    if np.nanmin(alt) < -90.0 or np.nanmax(alt) > 90.0:
        raise ValueError(
            f'alt outside [-90, 90] deg (min={np.nanmin(alt):.4g}, '
            f'max={np.nanmax(alt):.4g}); context looks scaled, not physical'
        )
    if np.nanmin(airmass) < 0.99:
        raise ValueError(
            f'airmass below 1 (min={np.nanmin(airmass):.4g}); '
            f'context looks scaled, not physical'
        )

    for feat, height_km in VAN_RHIJN_FEATURES.items():
        if feat not in names:
            continue
        stored = ctx[:, names.index(feat)]
        recomputed = np.asarray(_van_rhijn_factor(alt, float(height_km)), dtype=np.float64)
        if not np.allclose(stored, recomputed, rtol=rtol, atol=atol, equal_nan=True):
            worst = float(np.nanmax(np.abs(stored - recomputed)))
            raise ValueError(
                f"context column '{feat}' disagrees with van Rhijn recomputed "
                f"from 'alt' (max abs diff {worst:.4g}); context is not physical"
            )

    return True


def _airglow_ctx_columns(ctx_phys, ctx_names, check_physical=True):
    """Validate context and return (ctx, names, alt, airmass) in physical units."""
    ctx_phys = np.asarray(ctx_phys, dtype=np.float64)
    if ctx_phys.ndim == 1:
        ctx_phys = ctx_phys[None, :]
    names = [str(n).strip().lower() for n in ctx_names]
    if check_physical:
        assert_context_is_physical(ctx_phys, names)
    return (ctx_phys, names,
            ctx_phys[:, names.index('alt')],
            ctx_phys[:, names.index('airmass')])


def _coef_wavelengths_with_fallback(coef_wavelengths_a, group_indices, n_coef):
    """Per-coefficient wavelength array with group defaults filled in."""
    n_coef = int(n_coef)
    if coef_wavelengths_a is None:
        lam = np.full(n_coef, np.nan, dtype=np.float64)
    else:
        lam = np.asarray(coef_wavelengths_a, dtype=np.float64).astype(np.float64).copy()
        if lam.size != n_coef:
            raise ValueError(
                f'coef_wavelengths_a has {lam.size} entries, expected {n_coef}'
            )
    for group_name, idx in group_indices.items():
        fallback = GROUP_EFFECTIVE_WAVELENGTH_A.get(group_name)
        if fallback is None:
            continue
        idx = np.asarray(idx, dtype=int)
        if idx.size == 0:
            continue
        bad = ~(np.isfinite(lam[idx]) & (lam[idx] > 0.0))
        lam[idx[bad]] = float(fallback)
    return lam


def airglow_coef_extinction_k(coef_wavelengths_a, group_indices, n_coef,
                              coef_extinction_k=None):
    """Per-coefficient extinction coefficient (mag/airmass).

    If `coef_extinction_k` is supplied it is used verbatim -- that is the hook
    for the empirically fitted effective extinction (see
    fit_effective_extinction).  Otherwise the generic LCO stellar curve is
    interpolated at each coefficient's effective wavelength.
    """
    n_coef = int(n_coef)
    if coef_extinction_k is not None:
        k = np.asarray(coef_extinction_k, dtype=np.float64)
        if k.size != n_coef:
            raise ValueError(
                f'coef_extinction_k has {k.size} entries, expected {n_coef}'
            )
        return np.where(np.isfinite(k), k, 0.0)

    lam = _coef_wavelengths_with_fallback(coef_wavelengths_a, group_indices, n_coef)
    k = np.zeros(n_coef, dtype=np.float64)
    ok = np.isfinite(lam) & (lam > 0.0)
    if ok.any():
        k[ok] = np.atleast_1d(_lco_extinction_k(lam[ok])).astype(np.float64)
    return k


def airglow_van_rhijn_matrix(ctx_phys, ctx_names, group_indices, n_coef,
                             check_physical=True):
    """V(z; h) per row and coefficient; exactly 1.0 for non-airglow coefficients."""
    ctx_phys, _names, alt, _airmass = _airglow_ctx_columns(
        ctx_phys, ctx_names, check_physical)
    n_coef = int(n_coef)
    van_rhijn = np.ones((ctx_phys.shape[0], n_coef), dtype=np.float64)

    for group_name, idx in group_indices.items():
        if group_name not in AIRGLOW_GROUPS:
            continue
        feature_name = _group_height_feature_name(group_name)
        if feature_name is None or feature_name not in VAN_RHIJN_FEATURES:
            continue
        idx = np.asarray(idx, dtype=int)
        if idx.size == 0:
            continue
        height_km = float(VAN_RHIJN_FEATURES[feature_name])
        van_rhijn[:, idx] = np.asarray(
            _van_rhijn_factor(alt, height_km), dtype=np.float64)[:, None]

    return van_rhijn


def airglow_extinction_matrix(ctx_phys, ctx_names, group_indices, n_coef,
                              coef_wavelengths_a=None, coef_extinction_k=None,
                              check_physical=True):
    """10**(-0.4 k (X - 1)) per row and coefficient; 1.0 for non-airglow.

    Normalised to X = 1 rather than X = 0, so what this removes is the
    extinction *relative to zenith*.  The recovered quantity is therefore a
    zenith-equivalent amplitude, not an absolute layer emissivity.  That is
    self-consistent and cancels in the sky-to-science transfer.
    """
    ctx_phys, _names, _alt, airmass = _airglow_ctx_columns(
        ctx_phys, ctx_names, check_physical)
    n_coef = int(n_coef)
    k_all = airglow_coef_extinction_k(
        coef_wavelengths_a, group_indices, n_coef, coef_extinction_k)

    extinction = np.ones((ctx_phys.shape[0], n_coef), dtype=np.float64)
    for group_name, idx in group_indices.items():
        if group_name not in AIRGLOW_GROUPS:
            continue
        idx = np.asarray(idx, dtype=int)
        if idx.size == 0:
            continue
        extinction[:, idx] = 10.0 ** (
            -0.4 * k_all[idx][None, :] * (airmass[:, None] - 1.0))

    return extinction


def airglow_geometry_scale(
    ctx_phys,
    ctx_names,
    group_indices,
    n_coef,
    coef_wavelengths_a=None,
    coef_extinction_k=None,
    check_physical=True,
):
    """Geometry factor V(z; h) * 10**(-0.4 k (X - 1)) per row and coefficient.

    Returns an (n_rows, n_coef) array that is exactly 1.0 for every
    non-airglow coefficient (moon, other), since scattered moonlight needs
    scattering geometry rather than van Rhijn. The continuum group (HO2,
    FeO, O2Ac) IS airglow -- it's mesopause chemiluminescence, not aerosol
    or Rayleigh scattering -- and is treated at 87 km, same as OH.

    ctx_phys must be in physical units.  Divide physical coefficients by this
    to get zenith-equivalent amplitudes; multiply a prediction by it to get
    back to observed amplitudes at the target line of sight.

    Pass `coef_extinction_k` to use an empirically fitted effective extinction
    instead of the generic stellar curve.
    """
    van_rhijn = airglow_van_rhijn_matrix(
        ctx_phys, ctx_names, group_indices, n_coef, check_physical=check_physical)
    extinction = airglow_extinction_matrix(
        ctx_phys, ctx_names, group_indices, n_coef,
        coef_wavelengths_a=coef_wavelengths_a,
        coef_extinction_k=coef_extinction_k,
        check_physical=False)
    return np.clip(van_rhijn * extinction, 1e-6, None)


# ---------------------------------------------------------------------------
# Empirical effective airglow extinction
#
# The generic LCO curve is a *stellar* extinction curve.  Applying it to
# airglow over-corrects, because airglow is a quasi-uniform extended source:
# photons scattered out of the beam are largely replaced by photons scattered
# in from adjacent lines of sight, and scattering (Rayleigh + aerosol) is
# 70-100% of the total k everywhere airglow matters.  The effective airglow
# extinction is therefore well below the stellar value.
#
# Rather than model multiple scattering, fit the effective coefficient from
# the data.  The two sky pointings are simultaneous, so for coefficient j
#
#   ln(A_near/A_far) - ln(V_near/V_far) = -0.4 ln10 k_eff (X_near - X_far) + eps
#
# with eps the gravity-wave fluctuation between the two lines of sight (zero
# mean in log, uncorrelated with the airmass difference).  This is a Bouguer
# fit that uses airglow as its own source, and the recovered k_eff absorbs the
# multiple-scattering correction, the airglow-versus-stellar difference, the
# site aerosol level and any residual error in the assumed layer height.
#
# Two cautions, both documented in the methods section:
#   * h and k are nearly degenerate over the observed zenith range, so the
#     layer height is HELD FIXED here and only k_eff is fitted.  The product
#     V * 10^(-0.4 k (X-1)) is what the data constrain.
#   * a systematic horizontal gradient in layer brightness that correlates
#     with elevation at a fixed site would bias k_eff.  The returned table
#     includes split-half values so that can be checked.
# ---------------------------------------------------------------------------

def _airglow_height_per_coef(group_indices, n_coef):
    """Effective layer height per coefficient; NaN for non-airglow."""
    heights = np.full(int(n_coef), np.nan, dtype=np.float64)
    for group_name, idx in group_indices.items():
        if group_name not in AIRGLOW_GROUPS:
            continue
        feature_name = _group_height_feature_name(group_name)
        if feature_name is None or feature_name not in VAN_RHIJN_FEATURES:
            continue
        heights[np.asarray(idx, dtype=int)] = float(VAN_RHIJN_FEATURES[feature_name])
    return heights


def _masked_ols_rowconst(y, mask, x_row):
    """OLS of y_ij on [1, x_i] over masked entries, cluster-robust on rows.

    The airmass difference x is constant within a row, which collapses every
    sufficient statistic -- and the entire cluster-robust meat matrix -- to
    per-row counts and per-row sums:

        A_c' u_c = [sum_j u_ij, x_i sum_j u_ij] = g_i * [1, x_i]

    so no flattening and no per-cluster Python loop is needed.  The previous
    implementation looped over clusters (one iteration per row, via np.split)
    and was called about six times per wavelength bin, which dominated the
    runtime on large row counts.
    """
    y_masked = np.where(mask, y, 0.0)
    n_row = mask.sum(axis=1).astype(np.float64)
    s_row = y_masked.sum(axis=1)

    s1 = float(n_row.sum())
    sx = float((x_row * n_row).sum())
    sxx = float((x_row * x_row * n_row).sum())
    sy = float(s_row.sum())
    sxy = float((x_row * s_row).sum())

    xtx = np.array([[s1, sx], [sx, sxx]], dtype=np.float64)
    xtx_inv = np.linalg.pinv(xtx)
    beta = xtx_inv @ np.array([sy, sxy], dtype=np.float64)

    g = s_row - beta[0] * n_row - beta[1] * x_row * n_row
    g2 = g * g
    m01 = float((x_row * g2).sum())
    meat = np.array([[float(g2.sum()), m01],
                     [m01, float((x_row * x_row * g2).sum())]], dtype=np.float64)
    cov = xtx_inv @ meat @ xtx_inv
    return beta, cov, n_row


def fit_effective_extinction(
    coef_near,
    coef_far,
    ctx_near,
    ctx_far,
    ctx_names,
    group_indices,
    coef_wavelengths_a,
    n_wavelength_bins=8,
    wavelength_bin_edges=None,
    min_positive_fraction=0.80,
    clip_sigma=3.0,
    clip_iters=3,
    min_rows=100,
    min_pairs=500,
    clip_sample_rows=20000,
    seed=0,
    verbose=True,
):
    """Fit effective airglow extinction per wavelength bin from the sky pairs.

    Returns a DataFrame with one row per wavelength bin:
      lam_lo, lam_hi, lam_mid, n_coef, n_rows, n_pairs, retained_frac
      k_generic   -- generic LCO stellar curve at lam_mid, for comparison
      k_eff, k_eff_err  -- fitted effective extinction (cluster-robust error)
      intercept, intercept_err -- should be ~0; a significant value indicates a
                       relative throughput offset between the two sky channels
      k_eff_half1, k_eff_half2 -- split-half stability check
      ok          -- whether the bin is well enough constrained to be used

    Performance notes: van Rhijn is evaluated once per distinct layer height as
    an (n_rows,) vector rather than as an (n_rows, n_coef) matrix, the design
    statistics are accumulated from per-row sums instead of flattened
    (n_rows * n_coef) arrays, and the cluster-robust covariance is closed-form
    (see _masked_ols_rowconst).  Only the robust scale used for sigma clipping
    touches individual elements, and it is estimated on a row subsample.
    """
    rng = np.random.default_rng(seed)
    coef_near = np.asarray(coef_near, dtype=np.float64)
    coef_far = np.asarray(coef_far, dtype=np.float64)
    n_rows_all, n_coef = coef_near.shape

    ctx_near_arr, names, alt_near, airmass_near = _airglow_ctx_columns(
        ctx_near, ctx_names, check_physical=True)
    ctx_far_arr, _, alt_far, airmass_far = _airglow_ctx_columns(
        ctx_far, ctx_names, check_physical=False)
    delta_airmass = np.ascontiguousarray(airmass_near - airmass_far)

    height_per_coef = _airglow_height_per_coef(group_indices, n_coef)
    distinct_heights = np.unique(height_per_coef[np.isfinite(height_per_coef)])
    log_vr_ratio = {
        float(h): (np.log(_van_rhijn_factor(alt_near, float(h)))
                   - np.log(_van_rhijn_factor(alt_far, float(h))))
        for h in distinct_heights
    }

    lam = _coef_wavelengths_with_fallback(coef_wavelengths_a, group_indices, n_coef)
    airglow_cols = np.flatnonzero(np.isfinite(height_per_coef))
    if airglow_cols.size == 0:
        if verbose:
            print('fit_effective_extinction: no airglow coefficients; skipping.')
        return pd.DataFrame()

    if verbose:
        print(f'Effective-extinction fit: {airglow_cols.size} airglow coefficients, '
              f'{n_rows_all} rows, {len(distinct_heights)} distinct layer heights.')
        print(f'  airmass difference (near - far): '
              f'16-84% = {np.nanpercentile(delta_airmass, 16):+.3f}..'
              f'{np.nanpercentile(delta_airmass, 84):+.3f}, '
              f'|dX| median = {np.nanmedian(np.abs(delta_airmass)):.3f}')

    lam_ag = lam[airglow_cols]
    if wavelength_bin_edges is None:
        finite = np.isfinite(lam_ag)
        if finite.sum() < 2:
            return pd.DataFrame()
        wavelength_bin_edges = np.unique(np.nanpercentile(
            lam_ag[finite], np.linspace(0.0, 100.0, int(n_wavelength_bins) + 1)))
    wavelength_bin_edges = np.asarray(wavelength_bin_edges, dtype=np.float64)
    if wavelength_bin_edges.size < 2:
        return pd.DataFrame()

    half_assignment = rng.integers(0, 2, n_rows_all).astype(bool)
    slope_to_k = -1.0 / (0.4 * np.log(10.0))
    rows_out = []

    for b in range(wavelength_bin_edges.size - 1):
        lo, hi = wavelength_bin_edges[b], wavelength_bin_edges[b + 1]
        last = (b == wavelength_bin_edges.size - 2)
        in_bin = (lam_ag >= lo) & ((lam_ag <= hi) if last else (lam_ag < hi))
        cols = airglow_cols[in_bin]
        if cols.size == 0:
            continue

        near_block = coef_near[:, cols]
        far_block = coef_far[:, cols]

        # column-level selection only -- a per-row amplitude cut would be
        # selection on the outcome and biases the slope, because whichever arm
        # sits at lower airmass has the smaller van Rhijn factor and fails the
        # cut unless it carries a positive gravity-wave fluctuation
        positive = (near_block > 0.0) & (far_block > 0.0)
        keep_col = positive.mean(axis=0) >= float(min_positive_fraction)
        if not keep_col.any():
            continue
        cols = cols[keep_col]
        near_block = near_block[:, keep_col]
        far_block = far_block[:, keep_col]
        mask = positive[:, keep_col]

        # log amplitude ratio with van Rhijn removed, extinction left in.
        # van Rhijn is constant within a layer height, so it is subtracted as a
        # column vector per height rather than as a full matrix.
        with np.errstate(divide='ignore', invalid='ignore'):
            y_mat = np.log(np.where(mask, near_block, 1.0)) - np.log(
                np.where(mask, far_block, 1.0))
        for h in distinct_heights:
            sel = height_per_coef[cols] == h
            if sel.any():
                y_mat[:, sel] -= log_vr_ratio[float(h)][:, None]
        mask &= np.isfinite(y_mat)
        n_possible = int(mask.size)

        finite_x = np.isfinite(delta_airmass)
        if not finite_x.all():
            mask &= finite_x[:, None]

        for _ in range(int(clip_iters)):
            if mask.sum() < 10:
                break
            beta, _cov, _n_row = _masked_ols_rowconst(y_mat, mask, delta_airmass)
            resid = y_mat - (beta[0] + beta[1] * delta_airmass[:, None])
            if n_rows_all > int(clip_sample_rows):
                sample_rows = rng.choice(n_rows_all, int(clip_sample_rows), replace=False)
            else:
                sample_rows = np.arange(n_rows_all)
            vals = resid[sample_rows][mask[sample_rows]]
            if vals.size < 10:
                break
            med = float(np.median(vals))
            mad = 1.4826 * float(np.median(np.abs(vals - med)))
            if not np.isfinite(mad) or mad <= 0.0:
                break
            new_mask = mask & (np.abs(resid - med) < clip_sigma * mad)
            if new_mask.sum() == mask.sum():
                mask = new_mask
                break
            mask = new_mask

        n_pairs = int(mask.sum())
        lam_mid = float(np.nanmedian(lam[cols]))
        k_generic = float(_lco_extinction_k(lam_mid))
        if n_pairs < 10:
            continue

        beta, cov, n_row = _masked_ols_rowconst(y_mat, mask, delta_airmass)
        n_rows = int((n_row > 0).sum())
        retained = n_pairs / max(n_possible, 1)
        k_eff = float(beta[1] * slope_to_k)
        k_err = float(np.sqrt(max(cov[1, 1], 0.0)) * abs(slope_to_k))

        halves = []
        for which in (False, True):
            half_mask = mask & (half_assignment == which)[:, None]
            if half_mask.sum() >= 10:
                bh, _, _ = _masked_ols_rowconst(y_mat, half_mask, delta_airmass)
                halves.append(float(bh[1] * slope_to_k))
            else:
                halves.append(np.nan)

        rel_err = k_err / abs(k_eff) if k_eff != 0 else np.inf
        ok = bool(n_rows >= min_rows and n_pairs >= min_pairs
                  and np.isfinite(k_eff) and np.isfinite(k_err) and rel_err < 0.5)

        rows_out.append({
            'lam_lo': float(lo), 'lam_hi': float(hi), 'lam_mid': lam_mid,
            'n_coef': int(cols.size), 'n_rows': n_rows, 'n_pairs': n_pairs,
            'retained_frac': float(retained),
            'k_generic': k_generic, 'k_eff': k_eff, 'k_eff_err': k_err,
            'k_ratio': k_eff / k_generic if k_generic > 0 else np.nan,
            'intercept': float(beta[0]),
            'intercept_err': float(np.sqrt(max(cov[0, 0], 0.0))),
            'k_eff_half1': halves[0], 'k_eff_half2': halves[1],
            'ok': ok,
        })

    table = pd.DataFrame(rows_out)
    if verbose and not table.empty:
        print('  fitted effective extinction by wavelength bin:')
        cols_show = ['lam_mid', 'n_coef', 'n_rows', 'n_pairs', 'retained_frac',
                     'k_generic', 'k_eff', 'k_eff_err', 'k_ratio', 'intercept',
                     'k_eff_half1', 'k_eff_half2', 'ok']
        print(table[cols_show].to_string(
            index=False, float_format=lambda v: f'{v:.4g}'))
        low_ret = table[table['retained_frac'] < 0.5]
        if len(low_ret):
            print(f'  NOTE: {len(low_ret)} bin(s) retained under 50% of pairs. Heavy '
                  f'row-level loss reintroduces selection on the outcome and biases '
                  f'k_eff upward; consider raising min_positive_fraction so whole '
                  f'columns are dropped instead.')
        bad_int = table[np.abs(table['intercept'])
                        > 3.0 * table['intercept_err'].clip(lower=1e-12)]
        if len(bad_int):
            print(f'  NOTE: {len(bad_int)} bin(s) have an intercept >3 sigma from '
                  f'zero, which points to a relative throughput offset between '
                  f'the two sky channels rather than an extinction effect.')
    elif verbose:
        print('  effective-extinction fit produced no usable bins.')

    return table


def resolve_coef_extinction_k(
    coef_names,
    coef_wavelengths_a,
    group_indices,
    fit_table=None,
    clip_to_generic=True,
    verbose=True,
    coef_basis_k_generic=None,
):
    """Per-coefficient extinction, preferring the fitted k_eff over the generic curve.

    Well-constrained bins (`ok` in the fit table) are interpolated in
    wavelength; everything else falls back to the generic LCO curve.  With
    `clip_to_generic`, fitted values are restricted to [0, k_generic]: the
    multiple-scattering argument makes the effective airglow extinction a lower
    bound problem, so a fitted value above the stellar curve indicates noise or
    an unmodelled gradient rather than real physics.

    `coef_basis_k_generic`, when provided, is the per-coefficient B^2-weighted
    effective k from coef_wavelengths_from_basis(); where finite it overrides
    the point-sample `_lco_extinction_k(lambda_centroid)` as both the generic
    reference and the [0, k_generic] clip. For narrow-line coefficients this is
    a no-op (|f_j|^2 is peaked at the centroid); for broadband basis components
    (HO2, FeO, O2Ac) it correctly reflects the integrated behaviour of the LCO
    stellar curve across the basis function's support.
    """
    coef_names = [str(n) for n in coef_names]
    n_coef = len(coef_names)
    lam = _coef_wavelengths_with_fallback(coef_wavelengths_a, group_indices, n_coef)
    k_generic = airglow_coef_extinction_k(lam, group_indices, n_coef)
    if coef_basis_k_generic is not None:
        basis_k = np.asarray(coef_basis_k_generic, dtype=np.float64)
        if basis_k.size != n_coef:
            raise ValueError(
                f'coef_basis_k_generic has {basis_k.size} entries, expected {n_coef}'
            )
        replace = np.isfinite(basis_k) & (basis_k >= 0.0)
        if replace.any():
            k_generic = np.where(replace, basis_k, k_generic)
    k_out = k_generic.copy()
    source = ['generic'] * n_coef

    airglow_cols = np.concatenate(
        [np.asarray(idx, dtype=int) for g, idx in group_indices.items()
         if g in AIRGLOW_GROUPS and len(idx) > 0]
    ) if any(g in AIRGLOW_GROUPS and len(idx) > 0
             for g, idx in group_indices.items()) else np.zeros(0, dtype=int)

    usable = None
    if fit_table is not None and len(fit_table):
        usable = fit_table[fit_table['ok'].astype(bool)].sort_values('lam_mid')

    n_clipped = 0
    if usable is not None and len(usable) >= 1:
        lam_nodes = np.asarray(usable['lam_mid'], dtype=np.float64)
        k_nodes = np.asarray(usable['k_eff'], dtype=np.float64)
        for j in airglow_cols:
            if not np.isfinite(lam[j]):
                continue
            k_fit = float(np.interp(lam[j], lam_nodes, k_nodes,
                                    left=k_nodes[0], right=k_nodes[-1]))
            if clip_to_generic:
                k_clipped = float(np.clip(k_fit, 0.0, k_generic[j]))
                if k_clipped != k_fit:
                    n_clipped += 1
                k_fit = k_clipped
            k_out[j] = k_fit
            source[j] = 'fitted'

    if verbose:
        n_fit = sum(1 for s in source if s == 'fitted')
        print(f'Extinction resolved for {airglow_cols.size} airglow coefficients: '
              f'{n_fit} fitted, {airglow_cols.size - n_fit} generic.')
        if n_clipped:
            print(f'  {n_clipped} fitted value(s) clipped into [0, k_generic].')
        if airglow_cols.size:
            ratio = np.divide(k_out[airglow_cols], k_generic[airglow_cols],
                              out=np.full(airglow_cols.size, np.nan),
                              where=k_generic[airglow_cols] > 0)
            with np.errstate(invalid='ignore'):
                print(f'  k_eff / k_generic: median {np.nanmedian(ratio):.3f}, '
                      f'range {np.nanmin(ratio):.3f}..{np.nanmax(ratio):.3f}')

    return k_out, source



# ==========================================================================
# Extracted from notebook cell id=20eb0763
# ==========================================================================

# New helper: build simultaneous triplets (near, far -> sci) from decomposition files
def _load_decomp_with_row_index(
    decomp_fits_path,
    input_fits_path,
    context_columns,
    decomp_kind,
    return_chi2=False,
    return_err=False,
):
    """Load one decomp product and keep original source-row indices after finite filtering."""
    out = read_decomp_dataset(
        decomp_fits_path=decomp_fits_path,
        input_fits_path=input_fits_path,
        context_columns=context_columns,
        decomp_kind=decomp_kind,
        return_chi2=return_chi2,
        return_err=return_err,
    )

    # Unpack in the order used by read_decomp_dataset: base, then chi2, then err.
    _it = iter(out)
    coef_mat = next(_it); ctx_mat = next(_it); coef_names = next(_it); ctx_names = next(_it)
    chi2_used = next(_it) if return_chi2 else None
    coef_err_mat = next(_it) if return_err else None

    with fits.open(decomp_fits_path) as hdul_dec, fits.open(input_fits_path) as hdul_in:
        coef_tbl_full = _coerce_coef_hdu_to_table(hdul_dec['COEF'])
        coef_full, _ = _table_to_float32_matrix(coef_tbl_full, 'coefficient')
        meta_full = Table(hdul_in['META'].data)
        ctx_full, _ = _build_context_matrix(meta_full, context_columns, decomp_kind)
        meta_upper_full = {c.upper(): c for c in meta_full.colnames}
        obstime_mjd_full = _extract_obstime_mjd(meta_full, meta_upper_full)

    if coef_full.shape[0] != ctx_full.shape[0]:
        raise ValueError(
            f'Row count mismatch while building row_index: COEF has {coef_full.shape[0]} rows, META has {ctx_full.shape[0]} rows'
        )

    good = np.isfinite(coef_full).all(axis=1) & np.isfinite(ctx_full).all(axis=1)
    row_index = np.flatnonzero(good).astype(np.int64)
    obstime_mjd = np.asarray(obstime_mjd_full[good], dtype=np.float64)

    if row_index.size != coef_mat.shape[0]:
        raise RuntimeError(
            f'Internal row-index mismatch: row_index has {row_index.size}, data has {coef_mat.shape[0]}'
        )

    result = [coef_mat, ctx_mat, coef_names, ctx_names]
    if return_chi2:
        result.append(chi2_used)
    result.extend([row_index, obstime_mjd])
    if return_err:
        result.append(coef_err_mat)
    return tuple(result)


def build_triplet_coef_dataset(
    input_fits_path,
    sky_near_decomp_fits_path,
    sky_far_decomp_fits_path,
    sci_decomp_fits_path,
    context_columns,
    return_chi2=False,
    return_err=True,
):
    """Build aligned triplet arrays for coefficient-transfer experiments.

    ``return_err=True`` (default) also pulls the per-coefficient 1-sigma
    uncertainties from the ``COEF_ERR`` HDU (see fit.SkyDecomp._coef_err_active_set)
    and stores them as ``coef_err_near``, ``coef_err_far`` and ``coef_err_sci``.
    Older decomposition files without ``COEF_ERR`` yield all-NaN arrays.
    """

    def _load(kind, path):
        return _load_decomp_with_row_index(
            decomp_fits_path=path,
            input_fits_path=input_fits_path,
            context_columns=context_columns,
            decomp_kind=kind,
            return_chi2=return_chi2,
            return_err=return_err,
        )

    near = _load('sky1', sky_near_decomp_fits_path)
    far = _load('sky2', sky_far_decomp_fits_path)
    sci = _load('sci', sci_decomp_fits_path)

    def _unpack(out):
        it = iter(out)
        coef = next(it); ctx = next(it); cnames = next(it); xnames = next(it)
        chi2 = next(it) if return_chi2 else None
        row = next(it); mjd = next(it)
        err = next(it) if return_err else None
        return coef, ctx, cnames, xnames, chi2, row, mjd, err

    coef_near, ctx_near, coef_names_n, ctx_names_n, chi2_near, row_near, mjd_near, coef_err_near = _unpack(near)
    coef_far, ctx_far, coef_names_f, ctx_names_f, chi2_far, row_far, mjd_far, coef_err_far = _unpack(far)
    coef_sci, ctx_sci, coef_names_s, ctx_names_s, chi2_sci, row_sci, mjd_sci, coef_err_sci = _unpack(sci)

    if coef_names_n != coef_names_f or coef_names_n != coef_names_s:
        raise ValueError('Coefficient name mismatch across near/far/sci decomposition products.')
    if ctx_names_n != ctx_names_f or ctx_names_n != ctx_names_s:
        raise ValueError('Context name mismatch across near/far/sci products.')

    row_common = np.intersect1d(np.intersect1d(row_near, row_far), row_sci)
    if row_common.size == 0:
        raise ValueError('No shared source rows across near/far/sci after finite filtering.')

    idx_near = np.searchsorted(row_near, row_common)
    idx_far = np.searchsorted(row_far, row_common)
    idx_sci = np.searchsorted(row_sci, row_common)

    if (not np.array_equal(row_near[idx_near], row_common)
            or not np.array_equal(row_far[idx_far], row_common)
            or not np.array_equal(row_sci[idx_sci], row_common)):
        raise RuntimeError('Row-index alignment failure while building triplet dataset.')

    mjd_stack = np.column_stack([mjd_near[idx_near], mjd_far[idx_far], mjd_sci[idx_sci]])
    obstime_mjd = np.nanmedian(mjd_stack, axis=1).astype(np.float64)

    out = {
        'coef_near': coef_near[idx_near],
        'coef_far': coef_far[idx_far],
        'coef_sci': coef_sci[idx_sci],
        'ctx_near': ctx_near[idx_near],
        'ctx_far': ctx_far[idx_far],
        'ctx_sci': ctx_sci[idx_sci],
        'coef_names': coef_names_n,
        'ctx_names': ctx_names_n,
        'row_index': row_common.astype(np.int64),
        'obstime_mjd': obstime_mjd,
        'n_rows': int(row_common.size),
    }

    if return_chi2:
        out['chi2_near'] = chi2_near[idx_near]
        out['chi2_far'] = chi2_far[idx_far]
        out['chi2_sci'] = chi2_sci[idx_sci]

    if return_err:
        out['coef_err_near'] = coef_err_near[idx_near]
        out['coef_err_far'] = coef_err_far[idx_far]
        out['coef_err_sci'] = coef_err_sci[idx_sci]

    # Attach the science-pointing RA/Dec from META so downstream filters
    # can exclude tiles by sky region (see LMC/SMC exclusion in
    # apply_triplet_filters).  Try prefixed columns first, then unprefixed.
    try:
        with fits.open(input_fits_path) as _hdul_meta:
            _meta_tbl = Table(_hdul_meta['META'].data)
        _meta_upper = {c.upper(): c for c in _meta_tbl.colnames}
        _ra_col = next((_meta_upper[k] for k in ('SCI_RA', 'RA') if k in _meta_upper), None)
        _dec_col = next((_meta_upper[k] for k in ('SCI_DEC', 'DEC') if k in _meta_upper), None)
        if _ra_col is not None and _dec_col is not None:
            out['sci_ra'] = np.asarray(_meta_tbl[_ra_col], dtype=np.float64)[row_common]
            out['sci_dec'] = np.asarray(_meta_tbl[_dec_col], dtype=np.float64)[row_common]
            out['sci_radec_source'] = (_ra_col, _dec_col)
        else:
            print('  WARNING: sci RA/Dec columns not found in META '
                  '(looked for SCI_RA/SCI_DEC then RA/DEC); '
                  'field-level exclusion in apply_triplet_filters will be skipped.')
    except Exception as _exc:
        print(f'  WARNING: could not attach sci_ra/sci_dec '
              f'({type(_exc).__name__}: {_exc}).')

    print(
        f"Triplet dataset built: n_rows={out['n_rows']}, n_coef={out['coef_near'].shape[1]}, n_ctx={out['ctx_near'].shape[1]}"
        + (f" | sci pointing from META columns ({out['sci_radec_source'][0]}, {out['sci_radec_source'][1]})"
           if 'sci_radec_source' in out else "")
    )
    return out

# ==========================================================================
# Extracted from notebook cell id=ed6334b4
# ==========================================================================

# Reused reconstruction/prediction helpers for quick visual checks
def _meta_row_to_dict_upper(meta_row):
    names = list(meta_row.colnames) if hasattr(meta_row, 'colnames') else list(meta_row.dtype.names)
    return {str(k).upper(): k for k in names}

def _safe_float(x):
    arr = np.asarray(x)
    if arr.size == 0:
        raise ValueError('Empty value cannot be converted to float')
    if arr.shape != ():
        arr = arr.ravel()[0]
    return float(arr)


def _infer_base_dir_for_reconstruction():
    _cwd = Path.cwd().resolve()
    candidates = [_cwd / 'sky_decomp' / 'data',
                  _cwd.parent / 'skysub' / 'sky_decomp' / 'data',
                  _cwd, _cwd.parent]
    if 'PALACE_DIR' in globals():
        try:
            p = Path(PALACE_DIR).resolve()
            candidates.extend([p / 'sky_decomp' / 'data', p, p.parent])
        except Exception:
            pass

    # Prefer roots that have moon_zodi/ ROLO albedo (needed for split_zodi=True
    # reconstruction) alongside palace/PMD and a solar reference; fall back to
    # any root that at least has palace/PMD + solar for legacy no-split fits.
    def _has_solar(cand):
        return ((cand / 'Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt').exists()
                or (cand.parent / 'Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt').exists())

    with_moon_zodi = [c for c in candidates
                      if (c / 'palace' / 'PMD').exists()
                      and (c / 'moon_zodi' / 'eso_skycalc_rolo_moon_albedo.dat').exists()
                      and _has_solar(c)]
    if with_moon_zodi:
        return with_moon_zodi[0]
    for cand in candidates:
        if (cand / 'palace' / 'PMD').exists() and _has_solar(cand):
            return cand

    raise FileNotFoundError('Could not infer reconstruction base_dir containing palace/PMD, moon_zodi and solar reference file')

def load_lsf_state_if_available(decomp_fits_path, spectrum_index):
    """Return an LSFSurfaceState from a decomposition FITS row, or None if absent.

    The new sky_decomp.lsf_surface_iterative module writes LSF_COEF / LSF_KNOTS /
    LSF_META extensions alongside COEF and META.  Older decomposition outputs do
    not have those extensions; this returns None so callers can fall back to a
    Gaussian LSF from the input FITS.
    """
    path = Path(decomp_fits_path)
    if not path.exists():
        return None
    try:
        with fits.open(str(path)) as hdul:
            ext_names = {h.name for h in hdul}
            if 'LSF_COEF' not in ext_names:
                return None
        return load_lsf_surface_state(str(path), int(spectrum_index))
    except (KeyError, IndexError, ValueError) as exc:
        print(f'  LSF surface state unavailable in {path.name} row {spectrum_index}: '
              f'{type(exc).__name__}: {exc}')
        return None


def load_o2_vector_if_available(decomp_fits_path, spectrum_index):
    """Return the unit-integrated O2 template for one row, or None if absent.

    The reworked decomposition pipeline (fit.py + decompose_parallel.py) writes a
    VECTOR_O2 ImageHDU (n_rows x n_wave) alongside COEF / META. Older
    decomposition products predate the split of the O2 amplitude into the
    coefficient, so their VECTOR_O2 is missing and the O2 basis reconstructs to
    zero (matching the historical behaviour of reconstruct_component_spectra).
    """
    path = Path(decomp_fits_path)
    if not path.exists():
        return None
    try:
        with fits.open(str(path)) as hdul:
            if 'VECTOR_O2' not in {h.name for h in hdul}:
                return None
            data = np.asarray(hdul['VECTOR_O2'].data, dtype=np.float64)
    except (KeyError, IndexError, ValueError) as exc:
        print(f'  VECTOR_O2 unavailable in {path.name} row {spectrum_index}: '
              f'{type(exc).__name__}: {exc}')
        return None
    if data.ndim != 2 or int(spectrum_index) >= data.shape[0]:
        return None
    row = data[int(spectrum_index)]
    if not np.isfinite(row).any() or float(np.nansum(np.abs(row))) == 0.0:
        return None
    return row


def reconstruct_with_lsf(wave, coef, lsf, *, n_spline_knots=25, base_dir=None,
                         o2_vector=None, coef_err=None,
                         split_zodi=True, n_zodi_spline_knots=3,
                         palace_oh_suffix=None, palace_diffuse_suffix=None):
    """Reconstruct component spectra, dispatching on the LSF representation.

    ``lsf`` is either an ``LSFSurfaceState`` (uses the wavelength-dependent
    B-spline kernel via ``SkyDecompLSFSurfaceIterative._assemble_refined_matrices``)
    or a per-pixel Gaussian sigma vector / scalar (uses the parent-class
    ``reconstruct_component_spectra``).  Returns the same components dict as
    ``reconstruct_component_spectra``.

    If ``coef_err`` (per-coefficient 1σ, same shape as ``coef``) is
    provided, the returned dict additionally contains ``sigma`` (dict of
    per-component 1σ flux uncertainty per pixel) and ``sigma_total``
    (quadrature sum across independent components).  See
    ``SkyDecompBase._components_sigma_from_coef_err`` for the propagation.
    """
    if isinstance(lsf, LSFSurfaceState):
        model = SkyDecompLSFSurfaceIterative(
            wave,
            lsf_sigma=1.0,  # dummy; only stick matrices are used with the surface path
            n_spline_knots=n_spline_knots,
            base_dir=base_dir,
            palace_oh_suffix=palace_oh_suffix,
            palace_diffuse_suffix=palace_diffuse_suffix,
            split_zodi=bool(split_zodi),
            n_zodi_spline_knots=int(n_zodi_spline_knots),
        )
        coef_arr = np.asarray(coef, float).ravel()
        model._set_lsf_state(lsf)
        mats = model._assemble_refined_matrices()
        # VECTOR_O2 on disk is already convolved by the fitted LSF surface at fit
        # time; injecting it into matrix_o2_stick and letting _assemble_refined_matrices
        # re-convolve would double-broaden the O2 shape. Override the assembled O2
        # block verbatim instead.
        if o2_vector is not None:
            o2_vec = np.asarray(o2_vector, float).ravel()
            if o2_vec.shape != model.wave.shape:
                raise ValueError(
                    f'o2_vector shape mismatch: expected {model.wave.shape}, got {o2_vec.shape}'
                )
            mats['o2'] = o2_vec[None, :]
        n_expected = sum(m.shape[0] for m in mats.values())
        if coef_arr.size != n_expected:
            raise ValueError(
                f"reconstruct_with_lsf: coef size mismatch "
                f"(got {coef_arr.size}, expected {n_expected})")
        comps = model._components_from_coef(coef_arr, mats)
        comps['total'] = (comps['oh'] + comps['moon'] + comps.get('zodi', 0) + comps['diffuse']
                          + comps['atom'] + comps['orc'] + comps['o2'])
        if coef_err is not None:
            err_arr = np.asarray(coef_err, float).ravel()
            if err_arr.size != coef_arr.size:
                raise ValueError(
                    f'coef_err length mismatch: expected {coef_arr.size}, '
                    f'got {err_arr.size}'
                )
            sigma_comps = model._components_sigma_from_coef_err(err_arr, mats)
            comps['sigma'] = sigma_comps
            comps['sigma_total'] = np.sqrt(
                sigma_comps['oh'] ** 2
                + sigma_comps['moon'] ** 2
                + sigma_comps.get('zodi', np.zeros_like(sigma_comps['moon'])) ** 2
                + sigma_comps['diffuse'] ** 2
                + sigma_comps['atom'] ** 2
                + sigma_comps['orc'] ** 2
                + sigma_comps['o2'] ** 2
            )
        return comps
    return reconstruct_component_spectra(
        wave=wave, coef=coef, lsf_sigma=lsf,
        n_spline_knots=n_spline_knots, base_dir=base_dir, o2_vector=o2_vector,
        coef_err=coef_err,
        split_zodi=bool(split_zodi),
        n_zodi_spline_knots=int(n_zodi_spline_knots),
        palace_oh_suffix=palace_oh_suffix,
        palace_diffuse_suffix=palace_diffuse_suffix,
    )

# ==========================================================================
# Extracted from notebook cell id=cd99f8df
# ==========================================================================

# Filtering borrowed from the original notebook workflow, adapted to triplets
import plotly.express as px


def _angular_separation_deg_vec(ra_deg, dec_deg, ra_c_deg, dec_c_deg):
    """Great-circle angular separation (degrees), vectorised over (ra, dec)."""
    ra = np.deg2rad(np.asarray(ra_deg, dtype=np.float64))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=np.float64))
    ra_c = np.deg2rad(float(ra_c_deg))
    dec_c = np.deg2rad(float(dec_c_deg))
    cos_sep = (np.sin(dec) * np.sin(dec_c)
               + np.cos(dec) * np.cos(dec_c) * np.cos(ra - ra_c))
    return np.rad2deg(np.arccos(np.clip(cos_sep, -1.0, 1.0)))

LMC_EXCLUSION = {'name': 'LMC', 'ra_deg': 81, 'dec_deg': -69.7, 'radius_deg': 10.0}
SMC_EXCLUSION = {'name': 'SMC', 'ra_deg': 14, 'dec_deg': -73, 'radius_deg': 10}

def _kappa_sigma_row_mask(x, kappa=5.0, n_iter=3):
    x = np.asarray(x, dtype=np.float64)
    keep = np.isfinite(x).all(axis=1)
    if not np.any(keep):
        return keep

    for _ in range(n_iter):
        mu = np.nanmean(x[keep], axis=0)
        sig = np.nanstd(x[keep], axis=0)
        sig = np.where(np.isfinite(sig) & (sig > 0), sig, 1.0)
        within = np.all(np.abs(x - mu) <= (kappa * sig), axis=1)
        within &= np.isfinite(x).all(axis=1)
        new_keep = keep & within
        if new_keep.sum() == keep.sum() or new_keep.sum() == 0:
            break
        keep = new_keep

    return keep


def _kappa_mad_row_mask(x, kappa=4.0, n_iter=3):
    """Robust variant of _kappa_sigma_row_mask using per-column median + MAD.

    Median and MAD are not inflated by the outliers themselves, so at kappa=4
    the threshold reflects the actual bulk of the distribution.  Necessary for
    OH coefficients where a small fraction of rows carry decomposition failures
    at ~1e10-1e15 that would drag the mean/std of _kappa_sigma_row_mask so high
    that the outliers pass their own gate.  MAD is scaled by 1.4826 so kappa is
    directly comparable to sigma in the Gaussian limit.
    """
    x = np.asarray(x, dtype=np.float64)
    keep = np.isfinite(x).all(axis=1)
    if not np.any(keep):
        return keep

    for _ in range(n_iter):
        med = np.nanmedian(x[keep], axis=0)
        mad = np.nanmedian(np.abs(x[keep] - med), axis=0) * 1.4826
        mad = np.where(np.isfinite(mad) & (mad > 0), mad, 1.0)
        within = np.all(np.abs(x - med) <= (kappa * mad), axis=1)
        within &= np.isfinite(x).all(axis=1)
        new_keep = keep & within
        if new_keep.sum() == keep.sum() or new_keep.sum() == 0:
            break
        keep = new_keep

    return keep


def _loglog_slopes(component, log_wave, min_pixels=200):
    """Row-wise log-log slope of ``component`` against wavelength.

    Fits ``log f = a + b log(lambda)`` per row over finite, strictly positive
    samples only (the QP leaves exact zeros wherever a family is switched off,
    and those carry no colour information).  Returns NaN for rows with fewer
    than ``min_pixels`` usable samples or a degenerate lever arm.
    """
    flux = np.asarray(component, dtype=np.float64)
    x = np.asarray(log_wave, dtype=np.float64)[None, :]
    good = np.isfinite(flux) & (flux > 0.0)
    y = np.log(np.where(good, flux, 1.0))
    n = good.sum(axis=1).astype(np.float64)
    sx = np.where(good, x, 0.0).sum(axis=1)
    sy = np.where(good, y, 0.0).sum(axis=1)
    sxx = np.where(good, x * x, 0.0).sum(axis=1)
    sxy = np.where(good, x * y, 0.0).sum(axis=1)
    den = n * sxx - sx * sx
    ok = (n >= float(min_pixels)) & (den > 0.0)
    slope = np.full(flux.shape[0], np.nan, dtype=np.float64)
    np.divide(n * sxy - sx * sy, den, out=slope, where=ok)
    slope[~ok] = np.nan
    return slope


def split_zodi_reversal_diagnostics(
    decomp_fits_paths,
    row_index,
    wave,
    chunk_rows=256,
    min_component_frac=0.05,
):
    """Per-arm moon/zodi colours and role-reversal flags for the given rows.

    A *reversal* is a row where the fitted moon continuum is REDDER than the
    fitted zodi continuum -- ``moon_slope > zodi_slope`` -- i.e. the two
    families have swapped roles.  Physics puts the moon near -3.7 (scattered
    sunlight, Rayleigh-dominated) and the zodi near -0.3 (Leinert reddening),
    so the ordering is unambiguous when both families carry flux.

    Reversed rows are actively harmful to the ML stage rather than merely
    noisy: the moon coefficients of a reversed row describe zodiacal light and
    vice versa, so the network is asked to learn a moon-geometry mapping from
    zodi-shaped targets.  That is the same corruption the moon/zodi coupling
    head was quietly absorbing.

    Rows where one family is essentially switched off carry no ordering
    information and are NOT flagged: the test only applies where the moon's
    share of the moon+zodi continuum lies inside
    ``[min_component_frac, 1 - min_component_frac]``.  Dark-time rows sit at a
    moon share of ~0.02 under the current decomposition priors and are exempt
    by construction.

    Reads the refined ``COMP_MOON`` / ``COMP_ZODI`` component spectra written
    by ``decompose_parallel``, so no basis is rebuilt and no reconstruction
    convention has to be reproduced here.

    Returns ``{arm: {'moon_slope', 'zodi_slope', 'moon_frac', 'reversed',
    'testable'}}`` with one entry per requested row.
    """
    idx = np.asarray(row_index, dtype=np.int64)
    log_wave = np.log(np.asarray(wave, dtype=np.float64))
    out = {}
    for arm, path in dict(decomp_fits_paths).items():
        moon_slope = np.full(idx.size, np.nan)
        zodi_slope = np.full(idx.size, np.nan)
        moon_frac = np.full(idx.size, np.nan)
        try:
            with fits.open(path, memmap=True) as hdul:
                if 'COMP_MOON' not in hdul or 'COMP_ZODI' not in hdul:
                    print(f"  reversal check: {arm} has no COMP_MOON/COMP_ZODI HDU; skipping arm.")
                    continue
                for start in range(0, idx.size, int(chunk_rows)):
                    sel = idx[start:start + int(chunk_rows)]
                    moon = np.asarray(hdul['COMP_MOON'].data[sel], dtype=np.float64)
                    zodi = np.asarray(hdul['COMP_ZODI'].data[sel], dtype=np.float64)
                    sl = slice(start, start + sel.size)
                    moon_slope[sl] = _loglog_slopes(moon, log_wave)
                    zodi_slope[sl] = _loglog_slopes(zodi, log_wave)
                    m_tot = np.nansum(np.where(moon > 0, moon, 0.0), axis=1)
                    z_tot = np.nansum(np.where(zodi > 0, zodi, 0.0), axis=1)
                    total = m_tot + z_tot
                    moon_frac[sl] = np.where(total > 0, m_tot / np.maximum(total, 1e-300), np.nan)
        except FileNotFoundError:
            print(f"  reversal check: {arm} decomposition not found at {path}; skipping arm.")
            continue
        frac_lo = float(min_component_frac)
        testable = (
            np.isfinite(moon_slope) & np.isfinite(zodi_slope) & np.isfinite(moon_frac)
            & (moon_frac >= frac_lo) & (moon_frac <= 1.0 - frac_lo)
        )
        out[arm] = {
            'moon_slope': moon_slope,
            'zodi_slope': zodi_slope,
            'moon_frac': moon_frac,
            'separation': zodi_slope - moon_slope,
            'testable': testable,
        }
    return out


def split_zodi_decomp_paths(data_root, stem, suffix, every10=True):
    """Per-arm decomposition FITS paths, keyed the way the reversal helpers want."""
    base = f'{data_root}/{stem}' + ('_every10' if every10 else '')
    return {'near': f'{base}_decomp_sky1{suffix}.fits',
            'far':  f'{base}_decomp_sky2{suffix}.fits',
            'sci':  f'{base}_decomp_sci{suffix}.fits'}


def split_zodi_reversal_keep_mask(
    decomp_fits_paths,
    row_index,
    wave,
    min_component_frac=0.05,
    min_separation=0.0,
    label='sample',
    verbose=True,
):
    """Keep-mask (True = usable) dropping moon/zodi role-reversed rows.

    For the every10-derived DIAGNOSTIC samples.  These deliberately skip the
    training-time target filters (hard coefficient bounds, kappa-sigma) so the
    model's hardest genuine cases stay in the evaluation set -- but a reversal
    is not a hard case, it is a MISLABELLED one: the row's moon coefficients
    describe zodiacal light.  Scoring a prediction against swapped labels
    measures nothing, and such rows show up as spurious outliers in the moon
    panel while never having been trained on.  So this belongs with the chi2
    decomposition-failure gate, not with the target-outlier filters.

    Returns a boolean mask over ``row_index``; rows whose reversal state cannot
    be determined (one family switched off) are kept.
    """
    idx = np.asarray(row_index, dtype=np.int64)
    try:
        stats = split_zodi_reversal_diagnostics(
            decomp_fits_paths, idx, wave,
            min_component_frac=float(min_component_frac))
    except Exception as exc:  # missing COMP_* HDU, absent file, ...
        if verbose:
            print(f'  {label} reversal gate SKIPPED '
                  f'({type(exc).__name__}: {exc})')
        return np.ones(idx.size, dtype=bool)
    if not stats:
        if verbose:
            print(f'  {label} reversal gate SKIPPED (no arm had COMP_MOON/COMP_ZODI)')
        return np.ones(idx.size, dtype=bool)
    reversed_any = np.zeros(idx.size, dtype=bool)
    testable_any = np.zeros(idx.size, dtype=bool)
    for st in stats.values():
        reversed_any |= st['testable'] & (st['separation'] < float(min_separation))
        testable_any |= st['testable']
    keep = ~reversed_any
    if verbose:
        print(f'  {label} moon/zodi reversal gate: dropped '
              f'{int(reversed_any.sum())}/{idx.size} row(s) with swapped '
              f'moon/zodi colours in any arm '
              f'({int(testable_any.sum())} testable, '
              f'{int((~testable_any).sum())} with one family switched off)')
    return keep


def apply_triplet_filters(
    triplet_data,
    thin_every_n=1,
    chi2_qmax=90.0,
    chi2_min=0.0,
    chi2_max=10.0,
    hard_coef_bounds=None,
    kappa=6.0,
    kappa_iter=3,
    oh_kappa=4.0,
    oh_kappa_iter=3,
    exclude_field_regions=None,
    airmass_max=3.0,
    reversal_decomp_fits=None,
    reversal_wave=None,
    reversal_min_component_frac=0.05,
    reversal_min_separation=0.0,
):
    if hard_coef_bounds is None:
        hard_coef_bounds = {'feo': (0.0, 1.0), 'atom_k': (0.0, 1.0)}

    coef_names_local = [str(n) for n in triplet_data['coef_names']]
    coef_name_l = [n.lower() for n in coef_names_local]

    coef_near = np.asarray(triplet_data['coef_near'], dtype=np.float32)
    coef_far = np.asarray(triplet_data['coef_far'], dtype=np.float32)
    coef_sci = np.asarray(triplet_data['coef_sci'], dtype=np.float32)
    ctx_near = np.asarray(triplet_data['ctx_near'], dtype=np.float32)
    ctx_far = np.asarray(triplet_data['ctx_far'], dtype=np.float32)
    ctx_sci = np.asarray(triplet_data['ctx_sci'], dtype=np.float32)

    n0 = coef_near.shape[0]
    keep = np.ones(n0, dtype=bool)

    # Exclude tiles whose science pointing falls within a specified angular
    # radius of a listed sky region.  Defaults to LMC and SMC at 10 deg each --
    # both Clouds combine bright stellar populations, dense H II regions and
    # diffuse ionised gas at velocities close enough to airglow to blend with
    # the sky-component fit of coef_sci (methods sect 11.2).  Pass an empty
    # list to disable.
    if exclude_field_regions is None:
        exclude_field_regions = [LMC_EXCLUSION, SMC_EXCLUSION]
    if exclude_field_regions:
        if 'sci_ra' not in triplet_data or 'sci_dec' not in triplet_data:
            print("Science-field exclusion requested but triplet has no "
                  "sci_ra/sci_dec; skipping.  Re-run build_triplet_coef_dataset "
                  "so it attaches the science pointing coordinates.")
        else:
            sci_ra = np.asarray(triplet_data['sci_ra'], dtype=np.float64)
            sci_dec = np.asarray(triplet_data['sci_dec'], dtype=np.float64)
            field_mask = np.ones(n0, dtype=bool)
            for region in exclude_field_regions:
                reg_name = str(region.get('name', 'unnamed'))
                ra_c = float(region['ra_deg'])
                dec_c = float(region['dec_deg'])
                r_deg = float(region['radius_deg'])
                sep = _angular_separation_deg_vec(sci_ra, sci_dec, ra_c, dec_c)
                inside = np.isfinite(sep) & (sep <= r_deg)
                print(
                    f"Science-field exclusion around {reg_name} "
                    f"(ra={ra_c:.3f} deg, dec={dec_c:.3f} deg, radius={r_deg:.1f} deg): "
                    f"excluded {int(inside.sum())}/{inside.size} "
                    f"({100.0 * inside.mean():.1f}%)"
                )
                field_mask &= ~inside
            keep &= field_mask
            print(
                f"Combined science-field exclusion: kept {int(field_mask.sum())}/"
                f"{len(field_mask)} ({100.0 * field_mask.mean():.1f}%)"
            )

    ctx_names_l = [str(n).strip().lower() for n in triplet_data['ctx_names']]
    alt_idx = ctx_names_l.index('alt') if 'alt' in ctx_names_l else None
    airmass_idx = ctx_names_l.index('airmass') if 'airmass' in ctx_names_l else None
    vr87_idx = ctx_names_l.index('vanrhijn_87km') if 'vanrhijn_87km' in ctx_names_l else None
    vr95_idx = ctx_names_l.index('vanrhijn_95km') if 'vanrhijn_95km' in ctx_names_l else None
    vr285_idx = ctx_names_l.index('vanrhijn_285km') if 'vanrhijn_285km' in ctx_names_l else None

    if alt_idx is None:
        print("Altitude sanity filter: context column 'alt' not found; skipping altitude >= 0 check.")

    physical_mask = np.ones(n0, dtype=bool)
    if alt_idx is not None:
        alt_ok = (
            np.isfinite(ctx_near[:, alt_idx]) & (ctx_near[:, alt_idx] >= 0.0)
            & np.isfinite(ctx_far[:, alt_idx]) & (ctx_far[:, alt_idx] >= 0.0)
            & np.isfinite(ctx_sci[:, alt_idx]) & (ctx_sci[:, alt_idx] >= 0.0)
        )
        physical_mask &= alt_ok
        print(
            f"Altitude sanity filter (alt >= 0 in near/far/sci): kept {alt_ok.sum()}/{alt_ok.size} "
            f"({100.0 * alt_ok.mean():.1f}%)"
        )

    # airmass = sec(z) diverges near the horizon; anything above ~3 is not
    # science-usable and dominates the loss through the (X-1) extinction term.
    if airmass_idx is not None and airmass_max is not None:
        am_max = float(airmass_max)
        am_ok = (
            np.isfinite(ctx_near[:, airmass_idx]) & (ctx_near[:, airmass_idx] <= am_max)
            & np.isfinite(ctx_far[:, airmass_idx]) & (ctx_far[:, airmass_idx] <= am_max)
            & np.isfinite(ctx_sci[:, airmass_idx]) & (ctx_sci[:, airmass_idx] <= am_max)
        )
        physical_mask &= am_ok
        print(
            f"Airmass sanity filter (airmass <= {am_max:g} in near/far/sci): "
            f"kept {am_ok.sum()}/{am_ok.size} ({100.0 * am_ok.mean():.1f}%)"
        )
    elif airmass_idx is None and airmass_max is not None:
        print("Airmass sanity filter: context column 'airmass' not found; skipping.")

    for label, idx in (('vanrhijn_87km', vr87_idx), ('vanrhijn_95km', vr95_idx), ('vanrhijn_285km', vr285_idx)):
        if idx is None:
            continue
        vr_ok = (
            np.isfinite(ctx_near[:, idx]) & (ctx_near[:, idx] >= 1.0)
            & np.isfinite(ctx_far[:, idx]) & (ctx_far[:, idx] >= 1.0)
            & np.isfinite(ctx_sci[:, idx]) & (ctx_sci[:, idx] >= 1.0)
        )
        physical_mask &= vr_ok
        print(
            f"{label} sanity filter (finite and >= 1): kept {vr_ok.sum()}/{vr_ok.size} "
            f"({100.0 * vr_ok.mean():.1f}%)"
        )

    keep &= physical_mask
    print(
        f"Combined physical context filter: kept {physical_mask.sum()}/{len(physical_mask)} "
        f"({100.0 * physical_mask.mean():.1f}%)"
    )

    if all(k in triplet_data for k in ('chi2_near', 'chi2_far', 'chi2_sci')):
        chi2_stack = np.column_stack(
            [
                np.asarray(triplet_data['chi2_near'], dtype=np.float64),
                np.asarray(triplet_data['chi2_far'], dtype=np.float64),
                np.asarray(triplet_data['chi2_sci'], dtype=np.float64),
            ]
        )
        # Use max (not nanmax) so any-arm-NaN chi2 propagates NaN and
        # the row gets dropped by the isfinite gate below.  This
        # matches the intent that ANY per-arm decomposition failure
        # disqualifies the whole observation.
        chi2_combined = np.max(chi2_stack, axis=1)
        chi2_finite = chi2_combined[np.isfinite(chi2_combined)]
        chi2_hi = np.nanpercentile(chi2_finite, chi2_qmax)
        chi2_upper = min(float(chi2_max), float(chi2_hi)) if chi2_max is not None else float(chi2_hi)
        chi2_mask = np.isfinite(chi2_combined) & (chi2_combined >= chi2_min) & (chi2_combined <= chi2_upper)
        keep &= chi2_mask
        print(
            f"Triplet chi2 filter: min={chi2_min:.3g}, qmax={chi2_qmax:.1f}%=>{chi2_hi:.3g}, "
            f"upper={chi2_upper:.3g} | keep={chi2_mask.sum()}/{len(chi2_mask)} ({100.0*chi2_mask.mean():.1f}%)"
        )
        fig_chi2_triplet = px.histogram(
            x=chi2_combined[keep],
            nbins=80,
            title='Triplet combined reduced chi2 distribution (rows used for training)',
            labels={'x': 'max(reduced chi2 near/far/sci)', 'y': 'count'},
        )
        fig_chi2_triplet.update_layout(template='plotly_white', bargap=0.03)
        fig_chi2_triplet.show()
    else:
        print('Triplet chi2 columns not present; chi2 filtering skipped.')

    for cname, (lo, hi) in hard_coef_bounds.items():
        idxs = np.where(np.array(coef_name_l) == str(cname).lower())[0]
        if idxs.size == 0:
            print(f"Manual hard clip: coefficient {cname} not found; skipping.")
            continue

        j = int(idxs[0])
        within = (
            np.isfinite(coef_near[:, j]) & (coef_near[:, j] >= lo) & (coef_near[:, j] <= hi)
            & np.isfinite(coef_far[:, j]) & (coef_far[:, j] >= lo) & (coef_far[:, j] <= hi)
            & np.isfinite(coef_sci[:, j]) & (coef_sci[:, j] >= lo) & (coef_sci[:, j] <= hi)
        )
        keep &= within
        print(
            f"Manual hard clip {cname}: [{lo:.3g}, {hi:.3g}] | kept {within.sum()}/{within.size} ({100.0 * within.mean():.1f}%)"
        )

    # OH per-row-mean MAD filter.  Per-column filtering on OH doesn't work:
    # many OH columns have near-zero medians so their column-MAD is tiny and
    # legitimate rows look like outliers.  Instead we filter on the same
    # scalar the diagnostic histogram shows -- per-row mean(|OH|) taken as
    # the max across the three arms so any arm with a runaway decomposition
    # marks the row.  Robust median+MAD on this scalar puts pathological
    # rows (mean |OH| ~ 1e6 to 1e15) trillions of MADs above the bulk while
    # leaving rows with real airglow variability untouched.
    if oh_kappa is not None and float(oh_kappa) > 0:
        oh_idx_local = np.array(
            [j for j, n in enumerate(coef_name_l) if n.startswith('oh_')],
            dtype=int)
        if oh_idx_local.size > 0:
            _oh_row_scalar = np.maximum.reduce([
                np.nanmean(np.abs(coef_near[:, oh_idx_local]), axis=1),
                np.nanmean(np.abs(coef_far[:, oh_idx_local]), axis=1),
                np.nanmean(np.abs(coef_sci[:, oh_idx_local]), axis=1),
            ]).astype(np.float64)
            _finite = np.isfinite(_oh_row_scalar)
            _med = float(np.nanmedian(_oh_row_scalar[_finite]))
            _mad = float(np.nanmedian(np.abs(_oh_row_scalar[_finite] - _med))) * 1.4826
            _mad = max(_mad, 1e-30)
            _threshold = _med + float(oh_kappa) * _mad
            oh_mask = _finite & (_oh_row_scalar <= _threshold)
            keep &= oh_mask
            print(
                f"OH per-row-mean MAD filter (kappa={oh_kappa:.1f}, "
                f"n_oh_coefs={oh_idx_local.size} per arm; "
                f"median={_med:.3g}, MAD={_mad:.3g}, threshold={_threshold:.3g}): "
                f"kept {int(oh_mask.sum())}/{len(oh_mask)} "
                f"({100.0 * oh_mask.mean():.1f}%)"
            )

    # Moon/zodi role-reversal filter.  A reversed row hands the ML a
    # moon-labelled zodi spectrum, so it corrupts the geometry mapping the
    # network exists to learn rather than just adding noise -- one bad arm
    # disqualifies the observation, matching the chi2 filter's convention.
    reversal_stats = None
    if reversal_decomp_fits:
        if reversal_wave is None:
            print("Moon/zodi reversal filter: reversal_wave not given; skipping. "
                  "Pass the native wavelength grid (e.g. from "
                  "cfg.data.input_fits_for_basis).")
        else:
            reversal_stats = split_zodi_reversal_diagnostics(
                reversal_decomp_fits,
                np.asarray(triplet_data['row_index'], dtype=np.int64),
                reversal_wave,
                min_component_frac=float(reversal_min_component_frac),
            )
            rev_any = np.zeros(n0, dtype=bool)
            testable_any = np.zeros(n0, dtype=bool)
            moon_absent_all = np.ones(n0, dtype=bool)
            for arm, st in reversal_stats.items():
                # Which side is missing, for the untestable rows: a moon-absent
                # row (dark time, moon pinned near a 0.02 share) has no moon to
                # mislabel, whereas a zodi-absent row hides an undetermined
                # zodi colour.  They are not the same blind spot.
                moon_absent_all &= (
                    np.isfinite(st['moon_frac'])
                    & (st['moon_frac'] < float(reversal_min_component_frac))
                )
                arm_rev = st['testable'] & (st['separation'] < float(reversal_min_separation))
                st['reversed'] = arm_rev
                rev_any |= arm_rev
                testable_any |= st['testable']
                n_test = int(st['testable'].sum())
                if n_test == 0:
                    print(f"  {arm:>4s}: no rows with both families present")
                    continue
                t = st['testable']
                print(
                    f"  {arm:>4s}: reversed {int(arm_rev.sum())}/{n_test} testable "
                    f"({100.0 * arm_rev.sum() / n_test:.2f}%), "
                    f"median moon slope {np.nanmedian(st['moon_slope'][t]):+.2f}, "
                    f"zodi {np.nanmedian(st['zodi_slope'][t]):+.2f}, "
                    f"margin p10 {np.nanpercentile(st['separation'][t], 10):+.2f}"
                )
            rev_mask = ~rev_any
            keep &= rev_mask
            print(
                f"Moon/zodi reversal filter (moon_slope > zodi_slope in ANY arm, "
                f"tested where {reversal_min_component_frac:.2f} <= moon share <= "
                f"{1.0 - reversal_min_component_frac:.2f}): "
                f"kept {int(rev_mask.sum())}/{len(rev_mask)} "
                f"({100.0 * rev_mask.mean():.1f}%); "
                f"{int(rev_any.sum())} rows reversed, "
                f"{int((~testable_any).sum())} rows had no arm with both families present"
                f" (of those, {int((~testable_any & moon_absent_all).sum())} moon-absent"
                f" in every arm -- nothing to mislabel -- and"
                f" {int((~testable_any & ~moon_absent_all).sum())} with a"
                f" negligible zodi somewhere)"
            )
    else:
        print("Moon/zodi reversal filter: reversal_decomp_fits not given; skipping.")

    coef_concat = np.hstack([coef_near, coef_far, coef_sci]).astype(np.float32)
    kappa_mask = _kappa_sigma_row_mask(coef_concat, kappa=float(kappa), n_iter=int(kappa_iter))
    keep &= kappa_mask
    print(
        f"Kappa-sigma filter (kappa={kappa:.1f}): kept {kappa_mask.sum()}/{len(kappa_mask)} ({100.0 * kappa_mask.mean():.1f}%)"
    )

    # Thinning is applied LAST so the filter-fraction prints above reflect
    # counts against the full pre-thinning dataset.  It selects every N-th row
    # in the ORIGINAL row indexing; combined with the accumulated filter mask
    # via boolean AND, the result is the intersection.
    if int(thin_every_n) > 1:
        thin_mask = np.zeros(n0, dtype=bool)
        thin_mask[:: int(thin_every_n)] = True
        keep &= thin_mask
        print(f"Final thinning: every {int(thin_every_n)}-th of {n0} original rows; "
              f"combined with earlier filters, kept {int(keep.sum())} rows")
    else:
        print(f"Final thinning disabled: kept {int(keep.sum())}/{n0} rows after filters")

    if keep.sum() == 0:
        raise RuntimeError('Filtering removed all rows; relax thresholds.')

    out = {
        'coef_near': coef_near[keep],
        'coef_far': coef_far[keep],
        'coef_sci': coef_sci[keep],
        'ctx_near': ctx_near[keep],
        'ctx_far': ctx_far[keep],
        'ctx_sci': ctx_sci[keep],
        'coef_names': coef_names_local,
        'ctx_names': list(triplet_data['ctx_names']),
        'mask': keep,
        'row_index': np.asarray(triplet_data['row_index'])[keep],
    }
    if 'obstime_mjd' in triplet_data:
        out['obstime_mjd'] = np.asarray(triplet_data['obstime_mjd'], dtype=np.float64)[keep]
    for k in ('chi2_near', 'chi2_far', 'chi2_sci',
              'coef_err_near', 'coef_err_far', 'coef_err_sci',
              'sci_ra', 'sci_dec'):
        if k in triplet_data:
            out[k] = np.asarray(triplet_data[k])[keep]
    if 'sci_radec_source' in triplet_data:
        out['sci_radec_source'] = triplet_data['sci_radec_source']
    if reversal_stats:
        # Per-arm fitted colours, so the split can be audited without
        # re-reading the decomposition.  ``reversal`` is subset to the
        # surviving rows to line up with coef_*/ctx_*; ``reversal_all`` keeps
        # FULL-LENGTH arrays over the original rows, because the rows the
        # filter removed are the ones worth inspecting and they are by
        # definition absent from the subset.  ``mask`` selects survivors.
        out['reversal'] = {
            arm: {k: np.asarray(v)[keep] for k, v in st.items()}
            for arm, st in reversal_stats.items()
        }
        out['reversal_all'] = {
            arm: {k: np.asarray(v) for k, v in st.items()}
            for arm, st in reversal_stats.items()
        }

    print(
        f"Filtered triplet shapes: near={out['coef_near'].shape}, far={out['coef_far'].shape}, "
        f"sci={out['coef_sci'].shape}"
    )

    # Fold cyclic sin/cos pairs (az, moon_az, sun_az, moon_phase) back to
    # a single 0-360 degree axis per feature so histograms are readable.
    _display_names, _display_near = _decode_cyclic_context(out['ctx_names'], out['ctx_near'])
    _, _display_far = _decode_cyclic_context(out['ctx_names'], out['ctx_far'])
    _, _display_sci = _decode_cyclic_context(out['ctx_names'], out['ctx_sci'])

    # Distribution plot: one panel per context feature, three fields shown as
    # step lines (no fill).  Bin edges are computed per panel from that
    # feature's own min/max, so each panel has an independent x-range and
    # bin width -- which matters because the features span very different
    # ranges (degrees, airmass, sin/cos, van Rhijn factors).
    _hist_field_arrays = {
        'sky_near': _display_near,
        'sky_far': _display_far,
        'science': _display_sci,
    }
    _hist_field_colors = {
        'sky_near': '#1f77b4',
        'sky_far': '#9467bd',
        'science': '#2ca02c',
    }
    _hist_field_order = ['sky_near', 'sky_far', 'science']
    _hist_n_bins = 60

    context_order = list(_display_names)
    n_context = len(context_order)
    n_facet_cols = min(5, max(1, n_context))
    n_facet_rows = int(np.ceil(n_context / n_facet_cols))
    _n_panels_total = n_facet_rows * n_facet_cols

    fig_ctx_hist = make_subplots(
        rows=n_facet_rows,
        cols=n_facet_cols,
        subplot_titles=[c.replace('_', ' ') for c in context_order]
            + [''] * (_n_panels_total - n_context),
        horizontal_spacing=0.05,
        vertical_spacing=0.09,
    )

    for k, cname in enumerate(context_order):
        row = k // n_facet_cols + 1
        col = k % n_facet_cols + 1
        # sci_sep is 0 by construction for the science pointing; skip its delta at 0.
        panel_field_order = [f for f in _hist_field_order
                             if not (cname == 'sci_sep' and f == 'science')]
        per_field_vals = {}
        for f in panel_field_order:
            v = np.asarray(_hist_field_arrays[f][:, k], dtype=np.float64)
            per_field_vals[f] = v[np.isfinite(v)]
        combined = np.concatenate([per_field_vals[f] for f in panel_field_order])
        if combined.size < 2:
            continue
        lo = float(np.nanmin(combined))
        hi = float(np.nanmax(combined))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            continue
        edges = np.linspace(lo, hi, _hist_n_bins + 1)
        step_x = np.repeat(edges, 2)[1:-1]
        for f in panel_field_order:
            v = per_field_vals[f]
            if v.size == 0:
                continue
            counts, _ = np.histogram(v, bins=edges)
            step_y = np.repeat(counts, 2)
            fig_ctx_hist.add_trace(
                go.Scattergl(
                    x=step_x,
                    y=step_y,
                    mode='lines',
                    line=dict(color=_hist_field_colors[f], width=1.4),
                    name=f,
                    legendgroup=f,
                    showlegend=(k == 0),
                    hovertemplate=f + ': %{y}<extra></extra>',
                ),
                row=row,
                col=col,
            )

    fig_ctx_hist.for_each_annotation(lambda a: a.update(font=dict(size=11)))
    fig_ctx_hist.update_xaxes(showline=True, mirror=True, ticks='outside',
                              ticklen=4, showticklabels=True)
    fig_ctx_hist.update_yaxes(showline=True, mirror=True, ticks='outside',
                              ticklen=4, showticklabels=True)
    fig_ctx_hist.update_layout(
        template='plotly_white',
        title='Context parameter distributions by field (rows used after all filters)',
        height=n_facet_rows * 220 + 200,
        width=min(2000, n_facet_cols * 340 + 140),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0.0),
        margin=dict(l=60, r=20, t=90, b=50),
    )
    fig_ctx_hist.show()

    return out


# (Wiring moved to DataPipeline.load_and_prepare().)

# ==========================================================================
# Extracted from notebook cell id=ecliptic-ctx-augment
# ==========================================================================

# ECLIPTIC-CTX-V1: augment triplet ctx with PER-ARM ecliptic coordinates.
# Adds three features to each arm's ctx (following the same convention as
# `alt`, `az_sin`, `az_cos`: shared feature names, arm-specific values):
#   ecl_beta_deg -- raw ecliptic latitude in degrees (-90..+90). Zodi
#     amplitude peaks at beta=0 (ecliptic plane), so cos(beta) drives the
#     scaling; the MLP learns that from the raw latitude directly.
#   ecl_lon_sin/cos -- cyclic embedding of ecliptic longitude.
# `build_triplet_coef_dataset` in this notebook only attaches sci_ra/sci_dec;
# `_fill_arm_radec_from_meta_fits` reads SKY_NEAR/SKY_FAR RA/Dec directly from
# a meta FITS so per-arm ecliptic works even without changing the loader.
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
import astropy.units as u

ECLIPTIC_FEATURE_NAMES = ['ecl_beta_deg', 'ecl_lon_sin', 'ecl_lon_cos']
# Older names we strip on sight so re-running is idempotent.
_LEGACY_ECLIPTIC_NAMES = [
    'sci_ecl_lat_sin', 'sci_ecl_lat_cos',
    'sci_ecl_lon_sin', 'sci_ecl_lon_cos',
    'ecl_lat_sin', 'ecl_lat_cos',
]


def _ecliptic_features(ra_deg, dec_deg):
    _coords = SkyCoord(ra=np.asarray(ra_deg) * u.deg,
                       dec=np.asarray(dec_deg) * u.deg, frame='icrs')
    _ecl = _coords.transform_to(BarycentricMeanEcliptic())
    _lat_deg = np.asarray(_ecl.lat.degree, dtype=np.float32)
    _lon_rad = np.asarray(_ecl.lon.radian, dtype=np.float64)
    return np.stack([_lat_deg,
                     np.sin(_lon_rad).astype(np.float32),
                     np.cos(_lon_rad).astype(np.float32)],
                    axis=1).astype(np.float32)


def _fill_arm_radec_from_meta_fits(triplet, meta_fits_path=None):
    """Attach near/far RA/Dec to triplet by reading a meta FITS directly.

    Tries an explicit ``meta_fits_path`` first; otherwise falls back to
    ``{_DECOMP_STEM}_meta_only.fits`` then ``{_DECOMP_STEM}.fits``.
    """
    _row_idx = np.asarray(triplet['row_index'], dtype=int)
    _candidates = []
    if meta_fits_path is not None:
        _candidates.append(str(meta_fits_path))
    else:
        _stem = globals().get('_DECOMP_STEM', f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}')
        _candidates.extend([f'{_stem}_meta_only.fits', f'{_stem}.fits'])
    _meta = None
    _used_path = None
    for _p in _candidates:
        if not Path(_p).exists():
            continue
        try:
            with fits.open(_p) as _hdul:
                _meta = Table(_hdul['META'].data)
            _used_path = _p
            break
        except (KeyError, OSError) as _exc:
            print(f'  meta FITS candidate {_p} unusable: '
                  f'{type(_exc).__name__}: {_exc}')
    if _meta is None:
        print('  WARNING: no meta FITS available for per-arm RA/Dec fill; '
              'ecliptic will fall back to sci pointing for near/far arms')
        return
    _meta_up = {c.upper(): c for c in _meta.colnames}
    for _pref, _rakey, _deckey in (
        ('near', 'SKY_NEAR_RA', 'SKY_NEAR_DEC'),
        ('far', 'SKY_FAR_RA', 'SKY_FAR_DEC'),
    ):
        if f'{_pref}_ra' in triplet and f'{_pref}_dec' in triplet:
            continue
        if _rakey in _meta_up and _deckey in _meta_up:
            triplet[f'{_pref}_ra'] = np.asarray(
                _meta[_meta_up[_rakey]], dtype=np.float64)[_row_idx]
            triplet[f'{_pref}_dec'] = np.asarray(
                _meta[_meta_up[_deckey]], dtype=np.float64)[_row_idx]
    print(f'  filled per-arm RA/Dec from meta FITS ({_used_path})')


def _augment_triplet_with_ecliptic(triplet, force=True, meta_fits_path=None):
    """Append per-arm ecliptic features to ctx_near/ctx_far/ctx_sci in place.

    ``force=True`` (default): strip pre-existing ecliptic features and
    recompute from RA/Dec so re-running never leaves stale values in place.
    ``meta_fits_path``: optional override for the meta FITS used when the
    triplet is missing near/far RA/Dec (see ``_fill_arm_radec_from_meta_fits``).
    """
    _names = list(triplet.get('ctx_names', []))
    _to_strip = set(_LEGACY_ECLIPTIC_NAMES)
    if force:
        _to_strip.update(ECLIPTIC_FEATURE_NAMES)
    _strip_idx = [i for i, _n in enumerate(_names) if _n in _to_strip]
    if _strip_idx:
        _keep = [i for i in range(len(_names)) if i not in _strip_idx]
        for _arm in ('ctx_near', 'ctx_far', 'ctx_sci'):
            triplet[_arm] = np.asarray(triplet[_arm],
                                       dtype=np.float32)[:, _keep]
        _names = [_names[i] for i in _keep]
        triplet['ctx_names'] = _names
        print(f'  stripped {len(_strip_idx)} pre-existing ecliptic feature(s) '
              'so per-arm augment recomputes cleanly')
    if 'sci_ra' not in triplet or 'sci_dec' not in triplet:
        raise RuntimeError('triplet missing sci_ra/sci_dec; cannot compute '
                           'ecliptic features')
    if ('near_ra' not in triplet or 'far_ra' not in triplet
            or 'near_dec' not in triplet or 'far_dec' not in triplet):
        try:
            _fill_arm_radec_from_meta_fits(triplet, meta_fits_path=meta_fits_path)
        except Exception as _e:
            print(f'  meta FITS fallback failed: {type(_e).__name__}: {_e}')
    _sci_feats = _ecliptic_features(triplet['sci_ra'], triplet['sci_dec'])
    for _arm_prefix, _ctx_key in (('near', 'ctx_near'), ('far', 'ctx_far'),
                                  ('sci', 'ctx_sci')):
        if _arm_prefix == 'sci':
            _feats = _sci_feats
        elif f'{_arm_prefix}_ra' in triplet and f'{_arm_prefix}_dec' in triplet:
            _feats = _ecliptic_features(triplet[f'{_arm_prefix}_ra'],
                                        triplet[f'{_arm_prefix}_dec'])
        else:
            print(f'  WARNING: {_arm_prefix}_ra/{_arm_prefix}_dec unavailable; '
                  f'falling back to sci ecliptic for ctx_{_arm_prefix}')
            _feats = _sci_feats
        _prev = np.asarray(triplet[_ctx_key], dtype=np.float32)
        triplet[_ctx_key] = np.concatenate([_prev, _feats],
                                           axis=1).astype(np.float32)
    triplet['ctx_names'] = _names + ECLIPTIC_FEATURE_NAMES


# (Wiring moved to DataPipeline.load_and_prepare().)

# ==========================================================================
# Extracted from notebook cell id=physics-priors-augment
# ==========================================================================

# PHYSICS-PRIORS-CTX-V1: augment triplet ctx with the Leinert zodi V-band
# surface brightness predictor.
#
# The K-S moon-scattered-sky column that used to sit next to it was dropped
# 2026-08-21 (val+test |r|=0.043 vs 0.162 for zodi_log10_v: the K-S proxy is a
# smooth function of moon_alt / moon_sep / moon_phase / airmass that the
# encoder learns internally).  Its helpers and interpolation constant are
# removed from this variant.
import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u

PHYSICS_PRIOR_FEATURE_NAMES = [
    'zodi_log10_v',
    'moon_fli',           # (1 - cos(phase))/2, 0=new .. 1=full
    'moon_up_smooth',     # sigmoid(moon_alt/5), smooth horizon step
    'moon_airmass_up',    # moon LOS airmass, zeroed below horizon
    'moon_signal_proxy',  # rough moon flux at pointing (KS-like)
    # Phase A' (2026-08-26): explicit non-linear interaction features that
    # attack the residual_ctx_attribution RF rankings on continuum
    # (moon_fli * moon_phase_cos, moon_signal_proxy * ecl_lon_{cos,sin}).
    'moon_fli_x_phase_cos',   # 2nd-order moon phase (cos - cos^2)/2
    'moon_sig_x_lon_cos',     # moon-scatter x zodi anisotropy (cos)
    'moon_sig_x_lon_sin',     # moon-scatter x zodi anisotropy (sin)
]

# ---------------------------------------------------------------------------
# Leinert 1998 A&AS 127, 1, Table 17.  V-band zodi surface brightness in
# S10_sun/deg^2 units at 500 nm.  Grid: |helio-ecl-lon| (deg) x |ecl-lat| (deg).
# NaN corners = sun-avoidance zone; augment fills these with the ecliptic-pole
# floor.  Values are 15-20% accurate for physics-prior purposes.
# ---------------------------------------------------------------------------
_LEINERT_LON = np.array([
      0,   5,  10,  15,  20,  25,  30,  35,  40,  45,
     50,  55,  60,  65,  70,  75,  80,  85,  90,  95,
    100, 105, 110, 115, 120, 125, 130, 135, 140, 145,
    150, 155, 160, 165, 170, 175, 180,
], dtype=np.float64)
_LEINERT_BETA = np.array([
      0,   5,  10,  15,  20,  25,  30,  45,  60,  75,  90,
], dtype=np.float64)
_LEINERT_BRIGHTNESS = np.array([
    # cols: |beta| =  0,    5,   10,   15,   20,   25,   30,   45,   60,   75,   90
    [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,  610,  260,  130,   85],  # lon=  0
    [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1400,   580,  245,  128,   84],   # lon=  5
    [np.nan, np.nan, np.nan, np.nan, np.nan, 1600, 1350,   560,  240,  126,   84],     # lon= 10
    [np.nan, np.nan, np.nan, np.nan, 2200, 1550, 1300,   550,  240,  125,   83],       # lon= 15
    [np.nan, np.nan, np.nan, 2600, 2100, 1500, 1300,   540,  240,  125,   83],         # lon= 20
    [np.nan, np.nan, 2900, 2400, 2000, 1450, 1250,   535,  235,  122,   82],           # lon= 25
    [np.nan, 3200, 2600, 2200, 1900, 1400, 1200,   520,  225,  122,   82],             # lon= 30
    [ 2600, 2500, 2300, 2000, 1700, 1300, 1100,   490,  215,  120,   82],              # lon= 35
    [ 2000, 1900, 1800, 1600, 1400, 1200, 1000,   475,  210,  120,   82],              # lon= 40
    [ 1800, 1700, 1550, 1370, 1180, 1050,  900,   460,  210,  118,   81],              # lon= 45
    [ 1400, 1400, 1300, 1200, 1050,  950,  830,   440,  200,  115,   81],              # lon= 50
    [ 1100, 1100, 1050, 1000,  900,  830,  760,   420,  200,  115,   81],              # lon= 55
    [  850,  830,  790,  740,  680,  620,  580,   380,  195,  112,   80],              # lon= 60
    [  700,  690,  670,  640,  600,  560,  520,   360,  190,  111,   80],              # lon= 65
    [  600,  590,  570,  550,  530,  490,  460,   340,  185,  110,   80],              # lon= 70
    [  500,  495,  475,  445,  410,  385,  360,   300,  180,  108,   80],              # lon= 75
    [  440,  435,  420,  400,  380,  360,  340,   280,  178,  107,   80],              # lon= 80
    [  390,  385,  370,  355,  340,  325,  310,   260,  175,  106,   80],              # lon= 85
    [  340,  335,  325,  315,  300,  285,  275,   240,  170,  105,   80],              # lon= 90
    [  290,  287,  280,  270,  260,  252,  245,   225,  165,  103,   80],              # lon= 95
    [  260,  258,  253,  247,  238,  228,  225,   210,  160,  101,   80],              # lon=100
    [  245,  240,  237,  232,  225,  220,  215,   200,  155,   99,   80],              # lon=105
    [  225,  223,  221,  218,  213,  207,  205,   192,  150,   97,   80],              # lon=110
    [  215,  213,  211,  208,  205,  200,  195,   185,  148,   96,   80],              # lon=115
    [  200,  200,  198,  194,  190,  185,  180,   175,  145,   95,   80],              # lon=120
    [  195,  195,  193,  190,  185,  181,  178,   172,  143,   95,   80],              # lon=125
    [  190,  190,  188,  185,  180,  177,  175,   170,  140,   94,   80],              # lon=130
    [  185,  185,  183,  180,  175,  172,  170,   165,  138,   93,   80],              # lon=135
    [  186,  186,  184,  181,  176,  173,  170,   165,  138,   93,   80],              # lon=140
    [  188,  188,  186,  183,  178,  176,  173,   168,  138,   93,   80],              # lon=145
    [  195,  195,  192,  188,  185,  182,  180,   175,  138,   93,   80],              # lon=150
    [  200,  200,  197,  192,  188,  186,  184,   178,  140,   93,   80],              # lon=155
    [  210,  210,  205,  200,  194,  190,  188,   180,  142,   94,   80],              # lon=160
    [  225,  225,  220,  213,  205,  200,  195,   185,  145,   95,   80],              # lon=165
    [  240,  240,  233,  225,  213,  207,  200,   187,  147,   95,   80],              # lon=170
    [  250,  250,  245,  235,  222,  212,  205,   188,  148,   95,   80],              # lon=175
    [  260,  260,  253,  242,  228,  216,  208,   190,  148,   95,   80],              # lon=180 (gegenschein)
], dtype=np.float64)


def _leinert_zodi_log10_v(helio_lon_deg, beta_deg):
    """Bilinear interpolant of Leinert Table 17; returns log10(S10)."""
    _lon = np.abs(((np.asarray(helio_lon_deg, dtype=np.float64) + 180.0) % 360.0)
                   - 180.0)
    _beta = np.abs(np.asarray(beta_deg, dtype=np.float64))
    _log_tab = np.log10(_LEINERT_BRIGHTNESS)
    _lo_i = np.clip(np.searchsorted(_LEINERT_LON, _lon, side='right') - 1,
                     0, len(_LEINERT_LON) - 2)
    _lo_j = np.clip(np.searchsorted(_LEINERT_BETA, _beta, side='right') - 1,
                     0, len(_LEINERT_BETA) - 2)
    _fl = ((_lon - _LEINERT_LON[_lo_i])
           / (_LEINERT_LON[_lo_i + 1] - _LEINERT_LON[_lo_i]))
    _fb = ((_beta - _LEINERT_BETA[_lo_j])
           / (_LEINERT_BETA[_lo_j + 1] - _LEINERT_BETA[_lo_j]))
    _v00 = _log_tab[_lo_i, _lo_j]
    _v01 = _log_tab[_lo_i, _lo_j + 1]
    _v10 = _log_tab[_lo_i + 1, _lo_j]
    _v11 = _log_tab[_lo_i + 1, _lo_j + 1]
    return ((1 - _fl) * ((1 - _fb) * _v00 + _fb * _v01)
            + _fl * ((1 - _fb) * _v10 + _fb * _v11))


def _sun_ecliptic_longitude_deg(obstime_mjd):
    _t = Time(np.asarray(obstime_mjd, dtype=np.float64),
              format='mjd', scale='utc')
    _sun = get_sun(_t).transform_to(BarycentricMeanEcliptic())
    return np.asarray(_sun.lon.degree, dtype=np.float64)


def _augment_triplet_with_physics_priors(triplet, force=True):
    """Append Leinert V-band zodi surface brightness (one column per arm) to ctx.

    Requires ecl_beta_deg / ecl_lon_sin/cos (from ECLIPTIC-CTX-V1) and
    obstime_mjd for the sun ecliptic longitude used to project each arm's
    ecl_lon to helio-relative coordinates.
    """
    _names = list(triplet.get('ctx_names', []))
    _to_strip = set(PHYSICS_PRIOR_FEATURE_NAMES) if force else set()
    _strip_idx = [i for i, n in enumerate(_names) if n in _to_strip]
    if _strip_idx:
        _keep = [i for i in range(len(_names)) if i not in _strip_idx]
        for _arm in ('ctx_near', 'ctx_far', 'ctx_sci'):
            triplet[_arm] = np.asarray(
                triplet[_arm], dtype=np.float32)[:, _keep]
        _names = [_names[i] for i in _keep]
        triplet['ctx_names'] = _names
        print(f'  stripped {len(_strip_idx)} pre-existing physics-prior '
              f'feature(s) so augment recomputes cleanly')

    if 'obstime_mjd' not in triplet:
        raise RuntimeError('triplet missing obstime_mjd; cannot compute '
                            'Leinert zodi or K-S phase')

    _idx = {n: i for i, n in enumerate(_names)}
    _need = ('moon_alt', 'moon_sep', 'alt', 'moon_phase_cos',
             'ecl_beta_deg', 'ecl_lon_sin', 'ecl_lon_cos')
    _missing = [n for n in _need if n not in _idx]
    if _missing:
        raise RuntimeError(f'triplet ctx missing required features: '
                            f'{_missing}. Run ECLIPTIC-CTX-V1 first.')

    _sun_lon = _sun_ecliptic_longitude_deg(triplet['obstime_mjd'])

    _new_cols = {}
    for _arm_key in ('ctx_near', 'ctx_far', 'ctx_sci'):
        _arm_ctx = np.asarray(triplet[_arm_key], dtype=np.float64)
        _lon_s = _arm_ctx[:, _idx['ecl_lon_sin']]
        _lon_c = _arm_ctx[:, _idx['ecl_lon_cos']]
        _lon = np.rad2deg(np.arctan2(_lon_s, _lon_c)) % 360.0
        _helio_lon = ((_lon - _sun_lon + 180.0) % 360.0) - 180.0
        _zodi = _leinert_zodi_log10_v(
            _helio_lon, _arm_ctx[:, _idx['ecl_beta_deg']])
        # Sun-avoidance NaN (< ~30 deg from sun) -> ecliptic-pole floor.
        _zodi = np.where(np.isfinite(_zodi), _zodi, np.log10(80.0))
        # Moon-derived physical scaffolding (explicit nonlinear projections
        # of moon_alt / moon_sep / moon_phase so the encoder does not have to
        # reconstruct them internally from the raw angular ctx features).
        _phase_cos = _arm_ctx[:, _idx['moon_phase_cos']]
        _fli = 0.5 * (1.0 - _phase_cos)
        _alt = _arm_ctx[:, _idx['moon_alt']]
        _up_smooth = 1.0 / (1.0 + np.exp(-_alt / 5.0))
        _alt_rad = np.deg2rad(np.clip(_alt, 1.0, 89.0))
        _X_moon = 1.0 / np.sin(_alt_rad)
        _airmass_up = _up_smooth * _X_moon
        # 2026-08-24: V-band Krisciunas-Schaefer extinction on moon LOS.
        # k_V=0.15 mag/airmass. Fixes ~6x moon over-prediction at alt~4 deg.
        _moon_ext = 10.0 ** (-0.4 * 0.15 * _X_moon)
        _sep = _arm_ctx[:, _idx['moon_sep']]
        _signal_proxy = _fli * _up_smooth * _moon_ext / (1.0 + _sep / 45.0)**2
        _fli_x_phase_cos = _fli * _phase_cos
        _sig_x_lon_cos = _signal_proxy * _lon_c
        _sig_x_lon_sin = _signal_proxy * _lon_s
        _new_cols[_arm_key] = np.column_stack([
            _zodi, _fli, _up_smooth, _airmass_up, _signal_proxy,
            _fli_x_phase_cos, _sig_x_lon_cos, _sig_x_lon_sin,
        ]).astype(np.float32)

    for _arm_key, _new in _new_cols.items():
        _prev = np.asarray(triplet[_arm_key], dtype=np.float32)
        triplet[_arm_key] = np.concatenate(
            [_prev, _new], axis=1).astype(np.float32)
    triplet['ctx_names'] = _names + PHYSICS_PRIOR_FEATURE_NAMES

    print(f'  physics-prior augment applied: added '
          f'{PHYSICS_PRIOR_FEATURE_NAMES} per arm '
          f'(n_ctx now {len(triplet["ctx_names"])})')


# (Wiring moved to DataPipeline.load_and_prepare().)
