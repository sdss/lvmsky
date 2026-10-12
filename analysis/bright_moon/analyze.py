#!/usr/bin/env python3
"""Select DRP All full-Moon exposures and measure their sky continuum."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import html
import json
import mmap
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, get_body
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
SOURCE = Path('/data/sas/sdsswork/lvm/spectro/redux/1.2.1/drpall-1.2.1.fits')
SAS = Path('/data/sas')
# Same fixed site as skysub/sky_decomp/moon_zodi_model.py.
LCO = EarthLocation.from_geodetic(-70.692 * u.deg, -29.0089 * u.deg, 2281 * u.m)
iers.conf.auto_download = False
iers.conf.auto_max_age = None


def lunar_angle(mjd):
    t = Time(mjd, format='mjd', scale='utc')
    frame = GeocentricTrueEcliptic(equinox=t)
    moon = get_body('moon', t, ephemeris='builtin').transform_to(frame)
    sun = get_body('sun', t, ephemeris='builtin').transform_to(frame)
    return (moon.lon - sun.lon).wrap_at(360 * u.deg).to_value(u.rad)


def full_moons(start, end):
    grid = np.arange(np.floor(start) - 33, np.ceil(end) + 34, 0.25)
    angle = lunar_angle(grid)
    brackets = np.flatnonzero((np.sin(angle[:-1]) > 0) & (np.sin(angle[1:]) <= 0))
    roots = [brentq(lambda x: np.sin(lunar_angle(x)), grid[i], grid[i + 1], xtol=1e-7)
             for i in brackets]
    return np.asarray([x for x in roots if start - 2 <= x <= end + 2])


def geometry(mjd):
    t = Time(mjd, format='mjd', scale='utc')
    frame = AltAz(obstime=t, location=LCO, pressure=0 * u.hPa)
    moon = get_body('moon', t, location=LCO, ephemeris='builtin')
    sun = get_body('sun', t, location=LCO, ephemeris='builtin')
    return moon.transform_to(frame).alt.deg, sun.transform_to(frame).alt.deg


def measure(row):
    """Median 5400–5500 A continuum per good sky fiber, then across fibers."""
    result = {'expnum': int(row.expnum), 'measurement_status': 'ok'}
    try:
        with fits.open(SAS / row.location, memmap=True) as h:
            if h['FLUX'].header.get('BUNIT') != 'erg / (Angstrom s cm2)':
                raise ValueError('Unexpected flux units')
            wave = h['WAVE'].data
            # Sparse sky-fiber windows: avoid reading entire 489 MB SFrames into the page cache.
            h._file._mmap.madvise(mmap.MADV_RANDOM)
            band = np.flatnonzero((wave >= 5400) & (wave <= 5500))
            if len(band) < 100:
                raise ValueError('Continuum window coverage')
            sl = slice(band[0], band[-1] + 1)
            fibers = Table(h['SLITMAP'].data)
            for arm in ['SkyE', 'SkyW']:
                ids = np.flatnonzero((fibers['telescope'] == arm) & (fibers['fibstatus'] == 0))
                # Read narrow wavelength slices, keeping the large SFrames memory mapped.
                flux = np.array(h['FLUX'].data[ids, sl], dtype=float)
                flux += h['SKY'].data[ids, sl]
                good = (h['MASK'].data[ids, sl] == 0) & (h['IVAR'].data[ids, sl] > 0)
                good &= np.isfinite(flux)
                flux[~good] = np.nan
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    per_fiber = np.nanmedian(flux, axis=1)
                per_fiber[np.mean(good, axis=1) < 0.8] = np.nan
                values = per_fiber[np.isfinite(per_fiber)]
                result[arm.lower() + '_n_fibers'] = len(values)
                result[arm.lower() + '_continuum'] = float(np.median(values)) if len(values) >= 8 else np.nan
            result['fluxcal_header'] = h[0].header.get('FLUXCAL', '')
    except (OSError, ValueError, KeyError) as exc:
        result['measurement_status'] = f'{type(exc).__name__}: {exc}'
    return result


def figure(fig):
    fig.update_layout(template='plotly_white', hovermode='closest', height=440,
                      margin=dict(t=70, b=100), font=dict(size=14),
                      legend=dict(orientation='h', y=-0.25, x=0))
    fig.update_xaxes(showspikes=False)
    fig.update_yaxes(showspikes=False)
    return fig.to_html(full_html=False, include_plotlyjs=False,
                       config={'displaylogo': False, 'responsive': True})


def main():
    table = Table.read(SOURCE)
    df = table.to_pandas()
    for col in df.select_dtypes(include='object'):
        df[col] = df[col].map(lambda x: x.decode() if isinstance(x, bytes) else x)
    assert df.expnum.is_unique, 'Exposure identifiers must be unique'
    start = Time(df.obstime.to_numpy(dtype=str), format='isot', scale='utc').mjd
    mid = start + df.exptime.to_numpy() / 172800
    df['midpoint_mjd_utc'] = mid
    df['midpoint_utc'] = Time(mid, format='mjd').isot
    moons = full_moons(mid.min(), mid.max())
    nearest = np.abs(mid[:, None] - moons).argmin(axis=1)
    df['full_moon_utc'] = Time(moons[nearest], format='mjd').isot
    df['days_from_full_moon'] = mid - moons[nearest]
    df['moon_alt_mid_deg'], df['sun_alt_mid_deg'] = geometry(mid)
    moon = get_body('moon', Time(mid, format='mjd'), ephemeris='builtin')
    sun = get_body('sun', Time(mid, format='mjd'), ephemeris='builtin')
    elong = moon.separation(sun).rad
    phase = np.arctan2(sun.distance.to_value(u.km) * np.sin(elong),
                      moon.distance.to_value(u.km) - sun.distance.to_value(u.km) * np.cos(elong))
    df['illumination_mid'] = (1 + np.cos(phase)) / 2
    df['possible_lunar_eclipse'] = np.rad2deg(np.pi - elong) < 1.8
    df['within_48h'] = df.days_from_full_moon.abs() <= 2
    df['moon_up'] = df.moon_alt_mid_deg > 0
    df['astronomical_night'] = df.sun_alt_mid_deg < -18
    df['selected'] = df.within_48h & df.moon_up & df.astronomical_night & (df.exptime > 0)
    df['dark_control'] = False
    # ponytail: deterministic 300-exposure control sample; match fields and weather for causal inference.
    dark = df[(df.moon_alt_mid_deg < -5) & df.astronomical_night & (df.exptime > 0)]
    control_ids = dark.sample(n=min(300, len(dark)), random_state=42).index
    df.loc[control_ids, 'dark_control'] = True
    targets = df[df.selected | df.dark_control]
    cache_path = OUT / 'continuum_measurements.csv'
    cached = pd.read_csv(cache_path) if cache_path.exists() else pd.DataFrame(columns=['expnum'])
    pending = targets[~targets.expnum.isin(cached.expnum)]
    print(f'All: {len(df)}; selected: {df.selected.sum()}; measure: {len(pending)}', flush=True)
    records = []
    with ThreadPoolExecutor(max_workers=24) as pool:
        for i, result in enumerate(pool.map(measure, pending.itertuples(index=False)), 1):
            records.append(result)
            if i % 100 == 0:
                print(f'Spectra measured: {i}/{len(pending)}', flush=True)
                pd.concat([cached, pd.DataFrame(records)], ignore_index=True).to_csv(cache_path, index=False)
    cached = pd.concat([cached, pd.DataFrame(records)], ignore_index=True)
    cached.to_csv(cache_path, index=False)
    df = df.merge(cached, on='expnum', how='left', validate='one_to_one')
    df['sky_continuum'] = df[['skye_continuum', 'skyw_continuum']].median(axis=1)
    df['sky_continuum_raw'] = df.sky_continuum
    log_flux = np.log10(df.sky_continuum.where(df.sky_continuum > 0))
    log_centre = log_flux.median()
    log_mad_sigma = 1.4826 * (log_flux - log_centre).abs().median()
    df['continuum_gross_outlier'] = log_flux - log_centre > 8 * log_mad_sigma
    df.loc[df.continuum_gross_outlier, 'sky_continuum'] = np.nan
    df['sky_continuum_1e15'] = df.sky_continuum * 1e15
    df['group'] = np.where(df.selected, 'Full-Moon nights', 'Moon-down control')
    sel = df[df.selected].copy()
    sel['brightness_rank_in_cycle'] = sel.groupby('full_moon_utc').sky_continuum.rank(ascending=False, method='first')
    ma0, _ = geometry(start[df.selected])
    ma1, _ = geometry(start[df.selected] + df.loc[df.selected, 'exptime'].to_numpy() / 86400)
    sel['moon_up_at_both_endpoints'] = (ma0 > 0) & (ma1 > 0)
    counts = pd.DataFrame({'full_moon_utc': Time(moons, format='mjd').isot})
    totals = df.groupby('full_moon_utc').agg(all_exposures=('expnum', 'size'),
                within_48h=('within_48h', 'sum'), selected=('selected', 'sum'))
    stats = sel.groupby('full_moon_utc').agg(nights=('mjd', 'nunique'),
                illumination_min=('illumination_mid', 'min'), moon_alt_median_deg=('moon_alt_mid_deg', 'median'),
                measured=('sky_continuum', 'count'), median_continuum_1e15=('sky_continuum_1e15', 'median'))
    counts = counts.merge(totals, on='full_moon_utc', how='left').merge(stats, on='full_moon_utc', how='left')
    for col in ['all_exposures', 'within_48h', 'selected', 'nights', 'measured']:
        counts[col] = counts[col].fillna(0).astype(int)
    counts['calendar_month'] = counts.full_moon_utc.str[:7]
    monthly = counts.groupby('calendar_month', as_index=False).agg(full_moons=('full_moon_utc', 'size'), selected=('selected', 'sum'))
    top = sel[sel.brightness_rank_in_cycle <= 10].sort_values(['full_moon_utc', 'brightness_rank_in_cycle'])
    for name, data in [('all_exposures', df), ('selected_exposures', sel), ('full_moon_cycles', counts),
                       ('calendar_months', monthly), ('brightest_per_cycle', top),
                       ('moon_up_within_48h', df[df.within_48h & df.moon_up]),
                       ('calibration_outliers', df[df.continuum_gross_outlier]),
                       ('one_brightest_per_cycle', top[top.brightness_rank_in_cycle == 1])]:
        data.to_csv(OUT / f'{name}.csv', index=False)
    plot_df = df[(df.selected | df.dark_control) & (df.sky_continuum > 0)].copy()
    bright_med = sel.sky_continuum.median()
    dark_med = df.loc[df.dark_control, 'sky_continuum'].median()
    summary = dict(source=str(SOURCE), sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        all_exposures=len(df), selected=int(len(sel)), nights=int(sel.mjd.nunique()),
        cycles_with_exposures=int((counts.selected > 0).sum()), total_cycles=len(counts),
        measured_selected=int(sel.sky_continuum.notna().sum()), dark_control=int(df.dark_control.sum()),
        raw_measured_selected=int(sel.sky_continuum_raw.notna().sum()),
        continuum_gross_outliers=int(df.continuum_gross_outlier.sum()),
        outlier_log10_median=float(log_centre), outlier_log10_mad_sigma=float(log_mad_sigma),
        measured_dark_control=int(df.loc[df.dark_control, 'sky_continuum'].notna().sum()),
        broader_moon_up_sample=int((df.within_48h & df.moon_up).sum()),
        possible_eclipse_selected=int(sel.possible_lunar_eclipse.sum()),
        illumination_min=float(sel.illumination_mid.min()), illumination_median=float(sel.illumination_mid.median()),
        median_bright_continuum=float(bright_med), median_dark_continuum=float(dark_med),
        ratio_of_medians=float(bright_med / dark_med),
        moon_up_both_endpoints=int(sel.moon_up_at_both_endpoints.sum()))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    assert sel.within_48h.all() and sel.moon_up.all() and sel.astronomical_night.all()
    assert counts.selected.sum() == len(sel) == monthly.selected.sum()
    assert top.groupby('full_moon_utc').size().max() <= 10
    plots = [figure(px.bar(counts, x='full_moon_utc', y='selected', title='Exposure count at every full Moon',
                          labels={'full_moon_utc': 'Full Moon (UTC)', 'selected': 'Selected exposures'})),
             figure(px.scatter(sel, x='days_from_full_moon', y='moon_alt_mid_deg', color='illumination_mid', render_mode='svg',
                          hover_data=['expnum'], title='Bright lunar disk above the horizon',
                          labels={'days_from_full_moon': 'Days from full Moon', 'moon_alt_mid_deg': 'Moon altitude (degrees)', 'illumination_mid': 'Illuminated fraction'})),
             figure(px.box(plot_df, x='group', y='sky_continuum_1e15', log_y=True, points=False,
                          title='Measured sky continuum: full-Moon selection and dark control',
                          labels={'group': '', 'sky_continuum_1e15': 'Continuum (10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹)'})),
             figure(px.scatter(plot_df, x='moon_alt_mid_deg', y='sky_continuum_1e15', color='group', log_y=True, render_mode='svg',
                          hover_data=['expnum'], title='Measured background versus Moon altitude',
                          labels={'moon_alt_mid_deg': 'Moon altitude (degrees)', 'sky_continuum_1e15': 'Continuum (10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹)', 'group': 'Sample'}))]
    import plotly.offline
    (OUT / 'plotly.min.js').write_text(plotly.offline.get_plotlyjs())
    display_cols = ['expnum', 'mjd', 'midpoint_utc', 'full_moon_utc', 'brightness_rank_in_cycle',
                    'days_from_full_moon', 'illumination_mid', 'moon_alt_mid_deg', 'sci_moon_sep',
                    'sky_continuum_1e15', 'skye_n_fibers', 'skyw_n_fibers']
    names = {'full_moon_utc': 'Full Moon (UTC)', 'all_exposures': 'All exposures assigned',
             'within_48h': 'Within ±48 h', 'selected': 'Night + Moon up', 'nights': 'Observing nights',
             'illumination_min': 'Min. illuminated fraction', 'moon_alt_median_deg': 'Median Moon altitude (°)',
             'measured': 'Measured exposures', 'median_continuum_1e15': 'Median continuum (10⁻¹⁵)',
             'calendar_month': 'Full-Moon month', 'expnum': 'Exposure ID', 'mjd': 'Night ID (MJD)',
             'midpoint_utc': 'Midpoint (UTC)', 'brightness_rank_in_cycle': 'Rank in cycle',
             'days_from_full_moon': 'Days from full Moon', 'illumination_mid': 'Illuminated fraction',
             'moon_alt_mid_deg': 'Moon altitude (°)', 'sci_moon_sep': 'Science–Moon separation (°)',
             'sky_continuum_1e15': 'Continuum (10⁻¹⁵)', 'skye_n_fibers': 'SkyE fibers', 'skyw_n_fibers': 'SkyW fibers'}
    winner_cols = ['full_moon_utc', 'expnum', 'midpoint_utc', 'sky_continuum_1e15', 'moon_alt_mid_deg']
    winners = top[top.brightness_rank_in_cycle == 1][winner_cols]
    top_display = top[display_cols].copy()
    for col in ['brightness_rank_in_cycle', 'skye_n_fibers', 'skyw_n_fibers']:
        top_display[col] = top_display[col].astype('Int64')
    content = f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Bright Moon exposures — LVM DRP All</title>
<script src="plotly.min.js"></script><style>body{{font:17px/1.6 system-ui;margin:32px auto;max-width:1200px;padding:0 20px;color:#172536}}h1,h2{{line-height:1.2}}.cards{{display:flex;flex-wrap:wrap;gap:16px}}.card{{padding:18px;background:#edf4fa;border-radius:10px;flex:1;min-width:160px}}.card strong{{font-size:30px;display:block}}table{{border-collapse:collapse;font-size:13px;white-space:nowrap}}td,th{{padding:7px 10px;border-bottom:1px solid #dce3ea;text-align:right}}th{{background:#edf4fa}}.scroll{{overflow:auto;max-height:540px}}a{{color:#145da0}}.note{{border-left:4px solid #5689b8;padding:10px 18px;background:#f4f7fa}}</style>
<style>thead th{{position:sticky;top:0;z-index:1}}code{{overflow-wrap:anywhere}}.modebar{{font-size:12px}}</style>
<h1>Bright Moon exposures in LVM DRP All</h1>
<p>All {len(df):,} unique exposures in the Data Reduction Pipeline exposure catalogue (DRP All) 1.2.1, from {min(df.obstime)[:10]} to {max(df.obstime)[:10]}, were evaluated. The primary sample below applies the astronomical-night condition in addition to the ±48-hour and Moon-up cuts.</p>
<div class="cards"><div class="card"><strong>{len(sel):,}</strong>astronomical-night exposures, Moon up, ±48 h</div><div class="card"><strong>{summary['broader_moon_up_sample']:,}</strong>±48 h + Moon up, including twilight</div><div class="card"><strong>{sel.mjd.nunique()}</strong>observing nights</div><div class="card"><strong>{summary['ratio_of_medians']:.1f}×</strong>ratio of median measured continua</div></div>
<h2>Selection and interpretation</h2><p>An exposure is selected when its UTC midpoint lies within ±48 hours of a geocentric full Moon, the geometric Moon altitude at Las Campanas is above 0°, the Sun altitude is below −18° (astronomical night), and integration time is positive. Counts include all DRP All exposures meeting these conditions, with one row per exposure and no additional reduction-quality cut. Every source row has DRPQUAL = 0.</p>
<p>Selected illuminated fractions range from {sel.illumination_mid.min():.3f} to {sel.illumination_mid.max():.3f}, with median {sel.illumination_mid.median():.3f}. The Moon is above the horizon at both exposure endpoints for {sel.moon_up_at_both_endpoints.sum():,} of {len(sel):,} selected exposures. The midpoint rule defines the primary count.</p>
<p>Full Moon is the opposition of apparent geocentric ecliptic longitudes, found with 6-hour brackets and a numerical root. Astropy's built-in ephemeris supplies positions; apparent illumination uses the Sun–Moon–Earth phase angle. LCO is fixed at longitude −70.692°, latitude −29.0089°, elevation 2281 m, matching this project's sky-decomposition model; refraction is disabled. The installed IERS Earth-orientation table is used offline with its predictive-age restriction relaxed. Each lunation is retained independently, including calendar months with two full Moons. The dates are rounded to milliseconds for bookkeeping, not a claim of ephemeris accuracy.</p>
<p>The {int((df.within_48h & df.moon_up).sum()):,} exposures within ±48 hours with the Moon up form the broader sample before the astronomical-night cut. <a href="all_exposures.csv">The full audit table</a> records all selection flags, including twilight. The nearest full Moon is assigned to every exposure; the ±48-hour flag determines membership in its window.</p>
{plots[0]}{plots[1]}
<h2>Measured evidence for a bright background</h2>
<p>An SFrame is the pipeline's calibrated, sky-subtracted spectral file. SkyE and SkyW are the two dedicated sky telescopes. The continuum is measured directly from each selected SFrame: reconstruct pre-subtraction spectra as FLUX + SKY, retain SkyE/SkyW fibers with FIBSTATUS = 0 and pixels with MASK = 0 and positive IVAR, take the median over 5400–5500 Å per fiber, then the median across fibers for each sky telescope. A fiber requires ≥80% valid pixels; an arm requires ≥8 valid fibers. The exposure statistic is the median of the available arm medians. Units are calibrated flux density per fiber; no surface-brightness conversion is applied. This window avoids the strong 5577 Å airglow line.</p>
<p>Raw continuum measurements cover {summary['raw_measured_selected']:,} / {len(sel):,} selected exposures; {summary['measured_selected']:,} remain for comparison and ranking after the calibration check. A deterministic random sample (seed 42) of {summary['dark_control']} astronomical-night exposures with the Moon below −5° provides a descriptive control. Median continuum is {bright_med * 1e15:.3f} × 10⁻¹⁵ in the full-Moon sample and {dark_med * 1e15:.3f} × 10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹ in the control; their ratio is {bright_med / dark_med:.2f}.</p>
<p class="note">Calibration check: positive log₁₀ continuum values more than eight robust standard deviations (1.4826 × median absolute deviation) above the combined measured-sample median receive a gross-outlier flag. Faint values are retained. This marks {summary['continuum_gross_outliers']} exposure(s). Exposure 23816 has FLUXCAL = STD and raw SkyE/SkyW continua of 0.319/0.370 erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹, roughly 12 orders above adjacent MOD-calibrated exposures. Its raw values and selection membership are preserved; its continuum is treated as suspect for ranking and comparison. See <a href="calibration_outliers.csv">the calibration-outlier audit</a>. This conservative screen identifies gross scale discontinuities; ordinary calibration variation remains part of the measurement uncertainty.</p>
<p>Measured dark-control coverage: {summary['measured_dark_control']} / {summary['dark_control']} sampled exposures. A conservative geometric eclipse flag (Moon within 1.8° of the antisolar direction) marks {summary['possible_eclipse_selected']} selected exposures; eclipses can dim the Moon despite a large illuminated fraction. This flag is retained in CSVs and does not alter the primary count or measured-background ranking.</p>
<p class="note">The measured continuum is total sky-field light. Lunar geometry and the control comparison support a lunar-background interpretation; stars, diffuse emission, atmospheric conditions and flux calibration also affect it. The control is not matched in pointing, airmass, date or weather, so the ratio is descriptive rather than an isolated measurement of the lunar component. The ranking uses measured background, and can include non-lunar contaminants. Positive values are shown on the logarithmic plots; CSV files retain all measurements and statuses.</p>
{plots[2]}{plots[3]}
<h2>Counts per full-Moon cycle</h2><p>Zero-count cycles are retained. “All exposures” counts the nearest-cycle assignment; “within_48h” applies the time window only. “Selected” adds Moon-up and astronomical-night conditions. Continuum columns use units of 10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹.</p>
<p>DRP MJD is the pipeline observing-night identifier; precise UTC midpoints are recorded separately. {summary['cycles_with_exposures']} of {len(counts)} full-Moon cycles contain selected exposures.</p>
<div class="scroll">{counts.rename(columns=names).to_html(index=False, na_rep='—', float_format=lambda x: f'{x:.3f}')}</div>
<h2>One brightest-background exposure per cycle</h2><p>Continuum units below: 10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹ fiber⁻¹. These are candidates selected by measured total sky background.</p>
<div class="scroll">{winners.rename(columns=names).to_html(index=False, na_rep='—', float_format=lambda x: f'{x:.3f}')}</div>
<h2>Brightest exposures in each cycle</h2><p>Up to ten selected exposures per cycle, ranked by the measured sky continuum. Rank 1 is the recommended brightest-background exposure for that lunation. Sky telescope directions can differ from the science pointing; sci_moon_sep is the DRP science–Moon separation in degrees. Illumination and altitude are recomputed at midpoint. Full lists retain both sky-arm measurements and their fiber counts.</p>
<p>An em dash denotes an unavailable summary value.</p><details><summary>Show up to ten brightest exposures per cycle</summary><div class="scroll">{top_display.rename(columns=names).to_html(index=False, na_rep='—', float_format=lambda x: f'{x:.3f}')}</div></details>
<h2>Data and reproduction</h2><p><a href="selected_exposures.csv">All selected exposures</a> · <a href="brightest_per_cycle.csv">Ten brightest per cycle</a> · <a href="full_moon_cycles.csv">Cycle counts</a> · <a href="calendar_months.csv">Calendar-month totals</a> · <a href="all_exposures.csv">All-exposure audit</a> · <a href="continuum_measurements.csv">Measured continua and statuses</a> · <a href="summary.json">Summary and source checksum</a> · <a href="analyze.py">Reproduction script</a></p>
<p><a href="moon_up_within_48h.csv">Broader Moon-up sample, including twilight</a> · <a href="one_brightest_per_cycle.csv">One brightest exposure per cycle</a> · <a href="test_analysis.py">Runnable scientific checks</a></p>
<p>Source: <code>{html.escape(str(SOURCE))}</code>. SHA-256: <code>{summary['sha256']}</code>. Run <code>/home/ik52/miniconda3/envs/astro3.12/bin/python analyze.py</code>. The continuum cache is reused on rerun; remove it when changing the input release or measurement method.</p>
<p>Methods reference: <a href="https://docs.astropy.org/en/stable/api/astropy.coordinates.get_body.html">Astropy get_body</a> and <a href="https://docs.astropy.org/en/stable/coordinates/solarsystem.html">Solar System ephemerides</a>.</p></html>'''
    (OUT / 'index.html').write_text(content)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
