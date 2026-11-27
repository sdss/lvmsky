"""One table with every ionospheric and thermospheric quantity computed for the corpus, per exposure and telescope.

The corpus joins, row by row, the bundled ionosphere, thermosphere and ROTI tables (``sky_decomp/data/ionosphere``)
and adds what SkyML otherwise computes on the fly: the Swarm IBP bubble prior, the per-event dVTEC, the observed
F-region line coefficients of the decomposition, the sky-arm altimeter and q, exposure identifiers and solar
activity.  Columns are ``<arm>_<quantity>`` for the arms ``sci``, ``near``, ``far`` (``near`` is the sky telescope
named in ``sky_near_label``), or a plain name for per-exposure values.  ``COLUMNS`` gives unit, source and meaning;
it is also stored in the file as ``columns_json``.  All quantities are vertical or local except the observed line
coefficients, which are as seen along each line of sight (divide by ``<arm>_vanrhijn_<h>km`` for zenith-equivalent).

Run from ``skysub/`` (needs the corpus META and decomposition files and the ROTI station-day cache)::

    python -m mlp_predictor.iono_corpus --out sky_decomp/data/ionosphere/corpus/iono_corpus_<stem>.npz
    python -m mlp_predictor.iono_corpus --plot sky_decomp/data/ionosphere/corpus/iono_corpus_<stem>.npz
"""

import argparse
import json
from pathlib import Path

import numpy as np

from mlp_predictor import bubbles, ionosphere, roti, thermosphere

ARMS = ('sci', 'near', 'far')
LINES = {'ATOM_Or': '6300', 'ATOM_Orc_OI0777': '7774', 'ATOM_Orc_OI0845': '8446', 'ATOM_N': '5199', 'ATOM_Og': '5577'}
DECOMP_ROOT = '/Users/ik52/obs/sas/sdsswork/users/u6058164/medians_1.3.2_knox/sky_decomposition_1.3.2_knox'
DECOMP_STEM = 'lvmsframe_median_stack_1.3.2_knox_gaia1over100'

_ROW = {  # per-exposure columns: unit, source, meaning
    'obstime_mjd': ('d', 'META', 'mid-exposure time, UT MJD'),
    'expnum': ('', 'META', 'exposure number'),
    'night_mjd': ('d', 'derived', 'evening MJD of the night, floor(obstime_mjd - 0.5)'),
    'tile_id': ('', 'META', 'tile'),
    'near_is_skye': ('bool', 'META', 'True when SkyE is the near sky telescope (otherwise SkyW)'),
    'in_training_set': ('bool', 'notebook filters', 'row passes the training filters of the knox corpus (13,682 rows)'),
    'sun_alt': ('deg', 'META', 'Sun altitude at LCO'),
    'sun_alt_conjugate': ('deg', 'ionosphere', 'Sun altitude at the IGRF magnetic conjugate point of the pierce point'),
    'f107': ('sfu', 'solar table', 'daily F10.7 of the observing day'),
    'f107_prev': ('sfu', 'solar table', 'F10.7 of the previous day (NRLMSIS input)'),
    'f107a': ('sfu', 'solar table', '81-day mean F10.7'),
    'kp': ('', 'solar table', '3-hour Kp'),
    'ap': ('nT', 'solar table', '3-hour ap from Kp'),
    'sky_q': ('dex', 'decomposition + map', '|(log z7774_near - log z7774_far) - 2 (log VTEC_near - log VTEC_far)|: '
              'near/far mismatch beyond the map (bubble proxy)'),
}
_ARM = {  # per-arm columns
    'alt': ('deg', 'META', 'pointing altitude'),
    'az': ('deg', 'META', 'pointing azimuth (N through E)'),
    'ipp_lat': ('deg', 'ionosphere', 'latitude of the 285 km pierce point'),
    'ipp_lon': ('deg', 'ionosphere', 'longitude of the 285 km pierce point'),
    'local_time': ('h', 'derived', 'solar local time at the pierce point, -6..18 around midnight'),
    'vanrhijn_95km': ('', 'derived', 'thin-shell slant factor at 95 km (5577, Na layer)'),
    'vanrhijn_285km': ('', 'derived', 'thin-shell slant factor at 285 km (F-region lines)'),
    'ionex_source': ('', 'ionosphere', 'IONEX map used: UQRG, UPRG or CODE'),
    'vtec': ('TECU', 'IONEX', 'vertical TEC at the pierce point'),
    'vtec_rms': ('TECU', 'IONEX', 'map RMS at the pierce point'),
    'dvtec_dt': ('TECU/h', 'IONEX', 'VTEC rate, +-30 min'),
    'dlog10vtec_dt': ('dex/h', 'IONEX', 'relative VTEC rate, +-30 min'),
    'dlog10vtec_dlat': ('dex/deg', 'IONEX', 'map gradient north-south'),
    'dlog10vtec_dlon': ('dex/deg', 'IONEX', 'map gradient east-west'),
    'dip_lat': ('deg', 'IGRF', 'magnetic dip latitude atan(tan I / 2)'),
    'mag_dec': ('deg', 'IGRF', 'magnetic declination'),
    'ipp_dx_mag_mer': ('km', 'IGRF', 'pierce-point offset from LCO along the magnetic meridian'),
    'ipp_dx_mag_zon': ('km', 'IGRF', 'pierce-point offset from LCO across the magnetic meridian'),
    'hmf2': ('km', 'IRI-2020', 'F2 peak height'),
    'log10_nmf2': ('log10 m^-3', 'IRI-2020', 'F2 peak density, profile scaled to the map VTEC'),
    'b_bot': ('km', 'IRI-2020', 'bottomside thickness'),
    'log10_tec_scale': ('dex', 'IRI-2020', 'map VTEC / IRI VTEC'),
    'log10_o2': ('log10 cm^-3', 'NRLMSIS 2.1', 'O2 density at hmF2'),
    'log10_n2': ('log10 cm^-3', 'NRLMSIS 2.1', 'N2 density at hmF2'),
    'log10_o': ('log10 cm^-3', 'NRLMSIS 2.1', 'O density at hmF2'),
    'tn': ('K', 'NRLMSIS 2.1', 'neutral temperature at hmF2'),
    'log10_model6300': ('log10 R', 'Link & Cogger', 'model zenith 6300 brightness from the two profiles'),
    'h6300': ('km', 'Link & Cogger', 'emission-weighted mean height of the 6300 emission'),
    'log10_model7774': ('log10, rel.', 'IRI-2020', 'column of n_e^2 (7774 proxy)'),
    'wind_mer_mag': ('m/s', 'HWM14', 'magnetic-meridional wind at 285 km (positive north)'),
    'wind_w': ('m/s', 'HWM14', 'vertical drift driven by that wind, -U sin I cos I'),
    'bubble_p': ('', 'Swarm IBP + IBI', 'calibrated probability that a bubble reaches this field line'),
    'roti_max': ('TECU/min', 'GNSS RINEX', 'max ROTI of GPS rays within 200 km (field-aligned) and 10 min; NaN: none'),
    'roti_med': ('TECU/min', 'GNSS RINEX', 'median ROTI of those rays'),
    'roti_n': ('', 'GNSS RINEX', 'number of those ray samples'),
    'dvtec_event': ('TECU', 'GNSS RINEX', 'high-passed vertical TEC on the nearest ray (el >= 30 deg, +-5 min)'),
    'dvtec_event_dist_km': ('km', 'GNSS RINEX', 'field-aligned distance of that ray'),
    'dvtec_event_el': ('deg', 'GNSS RINEX', 'elevation of that ray'),
    'sky_altimeter': ('dex', 'decomposition', 'log z6300 - 0.5 log z7774 (zenith-equivalent); sky arms only'),
}
_ARM |= {f'coef_{v}': ('decomp. units', 'decomposition', f'observed {v} line coefficient (line of sight)')
         for v in LINES.values()}
COLUMNS = {k: v for k, v in _ROW.items()} | {f'{a}_{k}': v for a in ARMS for k, v in _ARM.items()}


def _match(mjd, query):
    """Index into sorted-unique ``mjd`` of the entry within 1e-7 d of each query time, -1 where absent."""
    pos = np.minimum(np.maximum(np.searchsorted(mjd, query), 1), mjd.size - 1)
    pos = pos - (np.abs(mjd[pos - 1] - query) < np.abs(mjd[pos] - query))  # the nearer of the two neighbours
    return np.where(np.abs(mjd[pos] - query) < 1e-7, pos, -1)


def build(out, decomp_root=DECOMP_ROOT, decomp_stem=DECOMP_STEM, training_mjd=None):
    from astropy.table import Table

    from mlp_predictor import data
    from mlp_predictor.config import DataConfig

    iono, thermo, rt = ionosphere._bundled_tables(), thermosphere.bundled_table(), roti.bundled_table()
    mjd = iono['obstime_mjd']
    n = mjd.size
    t = {'obstime_mjd': mjd, 'night_mjd': np.floor(mjd - 0.5).astype(np.int32),
         'sun_alt_conjugate': iono['sun_alt_conjugate']}
    for k in ('f107_prev', 'f107a', 'ap'):
        t[k] = thermo[k].astype(np.float32)
    t['f107'] = data._solar_activity_lookup(mjd, 'f107').astype(np.float32)
    t['kp'] = data._solar_activity_lookup(mjd, 'kp').astype(np.float32)

    cfg = DataConfig()
    cfg.decomp_data_root, cfg.decomp_stem = decomp_root, decomp_stem
    meta = Table.read(cfg.input_fits_meta, hdu='META')
    upper = {str(c).upper(): c for c in meta.colnames}
    m_mjd = np.asarray(data._extract_obstime_mjd(meta, upper), np.float64)
    order = np.argsort(m_mjd, kind='stable')
    j = order[_match(m_mjd[order], mjd)]
    if (_match(m_mjd[order], mjd) < 0).any():
        raise ValueError('some table rows have no META row')
    t['expnum'] = np.asarray(meta['expnum'])[j].astype(np.int64)
    t['tile_id'] = np.asarray(meta['tile_id'])[j].astype(np.int64)
    t['sun_alt'] = np.asarray(meta['sun_alt'], np.float32)[j]
    t['near_is_skye'] = np.char.lower(np.asarray(meta['sky_near_label']).astype(str))[j] == 'skye'

    for a in ARMS:
        for k in ('alt', 'az', 'ionex_source', *ionosphere._TABLE_FEATURES):
            t[f'{a}_{k}'] = iono[f'{a}_{k}']
        for k in thermosphere.FEATURES:
            t[f'{a}_{k}'] = thermo[f'{a}_{k}'].astype(np.float32)
        for k in ('roti_max', 'roti_med', 'roti_n'):
            t[f'{a}_{k}'] = rt[f'{a}_{k}']
        el = t[f'{a}_alt'].astype(np.float64)
        t[f'{a}_vanrhijn_95km'] = roti.mapping(el, 95.0).astype(np.float32)
        t[f'{a}_vanrhijn_285km'] = roti.mapping(el, ionosphere.IPP_HEIGHT_KM).astype(np.float32)
        t[f'{a}_local_time'] = bubbles.local_time_h(mjd, t[f'{a}_ipp_lon']).astype(np.float32)
        apex = bubbles.apex_height_km(t[f'{a}_dip_lat'].astype(np.float64), ionosphere.IPP_HEIGHT_KM)
        t[f'{a}_bubble_p'] = bubbles.reach_probability(mjd, t[f'{a}_ipp_lon'], apex, t['f107']).astype(np.float32)
        t[f'{a}_dvtec_event'], t[f'{a}_dvtec_event_dist_km'], t[f'{a}_dvtec_event_el'] = (
            np.full(n, np.nan, np.float32) for _ in range(3))
    for night in np.unique(t['night_mjd']):
        parts = []
        for sta in roti.STATIONS:
            for d in (night, night + 1):
                try:
                    parts.append(roti.station_day(sta, d, download=False))
                except FileNotFoundError:
                    pass
        if not parts:
            continue
        smp = {k: np.concatenate([p[k] for p in parts]) for k in roti.SAMPLE_KEYS}
        r = t['night_mjd'] == night
        for a in ARMS:
            f = roti.dvtec_event(smp, mjd[r], t[f'{a}_ipp_lat'][r], t[f'{a}_ipp_lon'][r], min_el=30.0)
            t[f'{a}_dvtec_event'][r], t[f'{a}_dvtec_event_dist_km'][r], t[f'{a}_dvtec_event_el'][r] = f.T

    tri = data.build_triplet_coef_dataset(
        input_fits_path=cfg.input_fits_meta, sky_near_decomp_fits_path=cfg.coef_fits('sky1'),
        sky_far_decomp_fits_path=cfg.coef_fits('sky2'), sci_decomp_fits_path=cfg.coef_fits('sci'),
        context_columns=list(cfg.context_columns), return_err=False)
    names = [str(c) for c in tri['coef_names']]
    k = _match(mjd, np.asarray(tri['obstime_mjd'], np.float64))
    keep = k >= 0
    for a in ARMS:
        for line, v in LINES.items():
            col = np.full(n, np.nan, np.float32)
            col[k[keep]] = np.asarray(tri[f'coef_{a}'], np.float32)[keep, names.index(line)]
            t[f'{a}_coef_{v}'] = col
    zlog = {a: {v: np.log10(np.where(t[f'{a}_coef_{v}'] > 0, t[f'{a}_coef_{v}'], np.nan) / t[f'{a}_vanrhijn_285km'])
                for v in ('6300', '7774')} for a in ('near', 'far')}
    t['sci_sky_altimeter'] = np.full(n, np.nan, np.float32)  # Science's own lines are the prediction target
    for a in ('near', 'far'):
        t[f'{a}_sky_altimeter'] = (zlog[a]['6300'] - 0.5 * zlog[a]['7774']).astype(np.float32)
    t['sky_q'] = np.abs(zlog['near']['7774'] - zlog['far']['7774']
                        - 2 * np.log10(t['near_vtec'] / t['far_vtec'])).astype(np.float32)
    t['in_training_set'] = (np.zeros(n, bool) if training_mjd is None
                            else _match(np.sort(np.asarray(training_mjd, np.float64)), mjd) >= 0)

    missing = set(COLUMNS) - set(t)
    if missing:
        raise RuntimeError(f'columns without values: {sorted(missing)}')
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, columns_json=json.dumps(COLUMNS), source_meta=Path(cfg.input_fits_meta).name,
                        models=json.dumps(thermosphere.model_versions()), **{c: t[c] for c in COLUMNS})
    print(f'wrote {out}: {n} rows, {len(COLUMNS)} columns')
    return t


def load(path):
    """The corpus as a dict of columns plus its column descriptions."""
    z = np.load(path, allow_pickle=False)
    cols = json.loads(str(z['columns_json']))
    return {c: z[c] for c in cols}, cols


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', help='build the corpus into this .npz')
    ap.add_argument('--training-triplet', help='pickled filtered triplet: marks in_training_set')
    ap.add_argument('--decomp-root', default=DECOMP_ROOT, help='directory of the median stack and its decompositions')
    ap.add_argument('--stem', default=DECOMP_STEM, help='median-stack file stem (file name without .fits)')
    ap.add_argument('--plot', help='write histograms of this corpus next to it (figs/)')
    args = ap.parse_args()
    if args.out:
        tm = None
        if args.training_triplet:
            import pickle

            tm = pickle.loads(Path(args.training_triplet).read_bytes())['obstime_mjd']
        build(args.out, decomp_root=args.decomp_root, decomp_stem=args.stem, training_mjd=tm)
    if args.plot:
        from mlp_predictor.diagnostics_cells import iono_corpus_hist

        iono_corpus_hist.plot(args.plot)


if __name__ == '__main__':
    main()
