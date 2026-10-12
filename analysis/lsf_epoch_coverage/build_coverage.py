"""Audit individual-fiber LSF coverage in the frozen 512-exposure SolarBlue run."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import hashlib
import os
import re
import sys
import html

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'skysub/medians_computation/src')]
from lvm_medians.gaia import _passband, _synthetic, _trapezoid_weights
from lvm_medians.stack import _flux_scale, gaia_ratio_clean_mask
from analysis.bright_moon.analyze import full_moons

RUN = Path('/data/sas/sdsswork/users/u6058164/lvm_solarblue_all_fibers_drp132_moonalt60_fli095')
DATA = Path('/data/sas/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2')
CACHE = DATA / 'gaia/fibers'
OUT = Path(os.environ.get('SAS_BASE_DIR', '/data/sas')) / 'sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132'
EPOCH_FILE = Path('/home/ik52/sdss5/lvmcore/calibrations/calibration-epochs.yaml')


def fit_files():
    original=sorted((RUN/'parameters').glob('solarblue-shard-*_sci_meta_coef_*.fits'))
    additional=sorted((OUT/'additional_clean_fits/parameters').glob('solarblue-shard-*_sci_meta_coef_*.fits'))
    return original+additional


def fitted_exposure_ids():
    ids=set()
    for path in fit_files():
        with fits.open(path,memmap=True) as h:ids.update(map(int,np.unique(h['INPUT_META'].data['expnum'])))
    return ids


def frame(data):
    df = Table(data).to_pandas()
    for c in df.select_dtypes('object'):
        df[c] = df[c].map(lambda v: v.decode().strip() if isinstance(v, bytes) else str(v).strip())
    return df


def gaia_path(expnum):
    path=CACHE / f'lvmGAIA-fibers-{int(expnum):08d}.fits'
    if not path.resolve().is_relative_to(DATA.resolve()):
        raise ValueError('Gaia source must be the user-specified gaia_fiber_matching_1.3.2 directory')
    return path


def photometry_provenance(selected):
    files=[]
    for row in selected.itertuples():
        path=gaia_path(row.expnum);digest=None
        if path.exists():
            with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
        stat=Path(row.source_path).stat()
        files.append(dict(expnum=int(row.expnum),gaia_path=str(path),gaia_sha256=digest,
            sframe_path=row.source_path,sframe_size=stat.st_size,sframe_mtime_ns=stat.st_mtime_ns))
    return dict(source_directory=str(CACHE.resolve()),selection_threshold=.01,validation_version=2,files=files)


def photometry(record):
    """Reuse coordinate-verified Gaia fluxes; measure the current DRP denominator."""
    exp = int(record['expnum'])
    cache = gaia_path(exp)
    with fits.open(record['source_path'], memmap=True) as h:
        slit = frame(h['SLITMAP'].data)
        slit = slit[slit.telescope.isin(['Sci', 'SkyE', 'SkyW'])].copy()
        slit['source_fiber_index'] = slit.index
        comments = ''.join(h[0].header.get('COMMENT', []))
        epochs = dict((cam.upper(), int(epoch)) for epoch, cam in
                      re.findall(r'/calib/(\d+)/lvm-mlsf-([brz][123])\.fits', comments))
        slit['exptime_seconds'] = float(h[0].header['EXPTIME'])
        slit['expnum'] = exp
        slit['gaia_verified'] = False
        slit['gaia_ratio'] = np.nan
        slit['synthetic_g'] = np.nan
        slit['synthetic_g_error'] = np.nan
        slit['g_valid_fraction'] = np.nan
        slit['gaia_g_flux'] = np.nan
        slit['n_gaia_sources'] = -1
        slit['gaia_clean_1pct'] = False
        slit['gaia_clean_10pct'] = False
        slit['gaia_cache'] = str(cache) if cache.exists() else 'unavailable'
        if cache.exists():
            table=frame(fits.getdata(cache,1))
            assert table.expnum.eq(exp).all() and table.fiberid.is_unique
            cached=table.set_index('fiberid').reindex(slit.fiberid)
            assert cached.index.is_unique
            current = SkyCoord(slit.ra.to_numpy()*u.deg, slit.dec.to_numpy()*u.deg)
            previous = SkyCoord(cached.fiber_ra.to_numpy()*u.deg, cached.fiber_dec.to_numpy()*u.deg)
            matched = current.separation(previous).arcsec < .05
            wave = np.asarray(h['WAVE'].data, dtype=float)
            fw, response, zero, pivot = _passband(None)
            assert np.allclose(cached.svo_g_zero_point_jy.dropna(), zero)
            assert np.allclose(cached.svo_g_pivot_angstrom.dropna(), pivot)
            weights = np.interp(wave, fw, response, left=0, right=0)*wave*_trapezoid_weights(wave)
            scale = _flux_scale(h[0].header, h['FLUX'])
            observed, errors, coverage = [], [], []
            for ids in np.array_split(slit.source_fiber_index.to_numpy(), 8):
                flux = (np.asarray(h['FLUX'].data[ids], dtype=float)
                        + np.asarray(h['SKY'].data[ids], dtype=float))*scale
                ivar = np.asarray(h['IVAR'].data[ids], dtype=float)/scale**2
                skyivar = np.asarray(h['SKY_IVAR'].data[ids], dtype=float)/scale**2
                variance = np.divide(1, ivar, out=np.full_like(ivar, np.nan), where=ivar>0)
                variance += np.divide(1, skyivar, out=np.full_like(skyivar, np.nan), where=skyivar>0)
                values, error, fraction = _synthetic(flux, variance, weights, h['MASK'].data[ids] == 0)
                observed.extend(values); errors.extend(error); coverage.extend(fraction)
            observed = np.asarray(observed)
            ratio = np.divide(cached.gaia_g_flux.to_numpy(), observed,
                              out=np.full(len(slit), np.nan), where=observed>0)
            phot = dict(n_gaia_sources=cached.n_gaia_sources.to_numpy(),
                        gaia_g_n_valid_sources=cached.gaia_g_n_valid_sources.to_numpy(),
                        lvm_flux_plus_sky_g=observed, lvm_flux_plus_sky_g_valid_fraction=coverage,
                        ratio_gaia_to_lvm_flux_plus_sky_g=ratio)
            verified = (matched & np.isfinite(observed) & (observed>0) & (np.asarray(coverage)>=.95)
                        & np.isfinite(cached.svo_g_zero_point_jy.to_numpy())
                        & np.isfinite(cached.svo_g_pivot_angstrom.to_numpy())
                        & (phot['n_gaia_sources']>=0)
                        & (phot['n_gaia_sources']==phot['gaia_g_n_valid_sources'])
                        & np.isfinite(cached.gaia_g_flux.to_numpy()) & (cached.gaia_g_flux.to_numpy()>=0))
            slit['gaia_verified'] = verified
            slit['gaia_ratio'] = ratio
            slit['synthetic_g'] = observed
            slit['synthetic_g_error'] = errors
            slit['g_valid_fraction'] = coverage
            slit['gaia_g_flux'] = cached.gaia_g_flux.to_numpy()
            slit['n_gaia_sources'] = cached.n_gaia_sources.to_numpy()
            for threshold, col in [(.01, 'gaia_clean_1pct'), (.1, 'gaia_clean_10pct')]:
                slit[col] = verified & gaia_ratio_clean_mask(phot, np.nextafter(threshold, 0))
    return slit[['expnum','fiberid','source_fiber_index','exptime_seconds','gaia_verified','gaia_ratio',
                 'synthetic_g','synthetic_g_error','g_valid_fraction','gaia_g_flux','n_gaia_sources',
                 'gaia_clean_1pct','gaia_clean_10pct','gaia_cache','telescope','spectrographid','fibstatus']], epochs


def summarize(group):
    good = group[group.accepted]
    nights = good.groupby('mjd').exptime_seconds.sum()
    seconds = float(nights.sum())
    times = np.sort(good.observation_mjd_utc.unique())
    return pd.Series(dict(n_selected=len(group), n_verified=int(group.gaia_verified.sum()),
        n_clean_1pct=int(group.gaia_clean_1pct.sum()), n_accepted=len(good),
        n_nights=len(nights), total_hours=seconds/3600,
        effective_nights=float(seconds**2/(nights**2).sum()) if seconds>0 else 0,
        largest_night_share=float(nights.max()/seconds) if seconds>0 else np.nan,
        time_span_days=float(np.ptp(times)) if len(times) else np.nan,
        largest_gap_days=float(np.max(np.diff(times))) if len(times)>1 else np.nan,
        min_days_from_full=good.days_from_full_moon.min(), max_days_from_full=good.days_from_full_moon.max(),
        min_fli=good.moon_fli.min(), max_fli=good.moon_fli.max(),
        median_g_snr=(good.synthetic_g/good.synthetic_g_error).median(),
        median_lsf_chi2=good.chi2_red.median(),
        measured=int(len(good)>0), repeat_exposures=int(len(good)>=3),
        repeat_nights=int(len(nights)>=2)))


def summarize_groups(rows, keys):
    groups = rows.groupby(keys, sort=True)
    base = groups.agg(n_selected=('accepted','size'), n_verified=('gaia_verified','sum'),
                      n_clean_1pct=('gaia_clean_1pct','sum'))
    good = rows[rows.accepted].copy()
    good['g_snr'] = good.synthetic_g/good.synthetic_g_error
    good = good.sort_values('observation_mjd_utc')
    good['gap'] = good.groupby(keys).observation_mjd_utc.diff()
    accepted = good.groupby(keys).agg(n_accepted=('accepted','size'),
        n_nights=('mjd','nunique'), total_seconds=('exptime_seconds','sum'),
        time_first=('observation_mjd_utc','min'),time_last=('observation_mjd_utc','max'),
        largest_gap_days=('gap','max'), min_days_from_full=('days_from_full_moon','min'),
        max_days_from_full=('days_from_full_moon','max'),min_fli=('moon_fli','min'),max_fli=('moon_fli','max'),
        median_g_snr=('g_snr','median'),median_lsf_chi2=('chi2_red','median'))
    nights = good.groupby(keys+['mjd']).exptime_seconds.sum().rename('seconds').reset_index()
    nights['squared_seconds'] = nights.seconds**2
    concentration = nights.groupby(keys).agg(squared_sum=('squared_seconds','sum'),max_seconds=('seconds','max'))
    out = base.join(accepted).join(concentration)
    count_cols=['n_accepted','n_nights','total_seconds']
    out[count_cols]=out[count_cols].fillna(0)
    out['total_hours']=out.total_seconds/3600
    out['effective_nights']=(out.total_seconds**2/out.squared_sum).fillna(0)
    out['largest_night_share']=out.max_seconds/out.total_seconds
    out['time_span_days']=out.time_last-out.time_first
    out['measured']=(out.n_accepted>0).astype(int)
    out['repeat_exposures']=(out.n_accepted>=3).astype(int)
    out['repeat_nights']=(out.n_nights>=2).astype(int)
    return out.drop(columns=['total_seconds','squared_sum','max_seconds','time_first','time_last']).reset_index()


def write_time_overview(selected):
    """Expose exact observing times against nominal epoch validity intervals."""
    definitions=yaml.safe_load(EPOCH_FILE.read_text())['epochs']
    starts=sorted(definitions)
    ends={e:starts[i+1] if i+1<len(starts) else np.nan for i,e in enumerate(starts)}
    mapping=selected.copy()
    mapping['epoch_start_mjd']=mapping.nominal_epoch
    mapping['epoch_end_mjd_exclusive']=mapping.nominal_epoch.map(ends)
    mapping['epoch_start_utc']=Time(mapping.epoch_start_mjd.to_numpy(),format='mjd').isot
    mapping['epoch_end_utc_exclusive']=[Time(v,format='mjd').isot if np.isfinite(v) else 'open' for v in mapping.epoch_end_mjd_exclusive]
    assert mapping.expnum.is_unique
    assert (mapping.mjd>=mapping.epoch_start_mjd).all()
    assert ((mapping.mjd<mapping.epoch_end_mjd_exclusive)|mapping.epoch_end_mjd_exclusive.isna()).all()
    mapping.to_csv(OUT/'exposure_epoch_mapping.csv',index=False)
    timing=mapping.groupby(['nominal_epoch','epoch_start_utc','epoch_end_utc_exclusive']).agg(
        exposures=('expnum','size'),observing_nights=('mjd','nunique'),
        first_exposure_utc=('obstime','min'),last_exposure_utc=('obstime','max'),
        full_moon_neighborhoods=('full_moon_utc','nunique')).reset_index()
    timing.to_csv(OUT/'epoch_time_coverage.csv',index=False)
    nights=mapping.groupby(['nominal_epoch','mjd','full_moon_utc']).agg(exposures=('expnum','size'),
        first_exposure_utc=('obstime','min'),last_exposure_utc=('obstime','max')).reset_index()
    nights.to_csv(OUT/'observing_night_epoch_mapping.csv',index=False)
    import matplotlib.dates as mdates
    fig,ax=plt.subplots(figsize=(14,8))
    stop=pd.to_datetime(mapping.obstime,utc=True).max()+pd.Timedelta(days=4)
    for y,row in enumerate(timing.itertuples()):
        left=pd.to_datetime(row.epoch_start_utc,utc=True)
        right=pd.to_datetime(row.epoch_end_utc_exclusive,utc=True) if row.epoch_end_utc_exclusive!='open' else stop
        ax.hlines(y,left,right,color='#c7c7c7',linewidth=7,zorder=1)
        times=pd.to_datetime(mapping.loc[mapping.nominal_epoch.eq(row.nominal_epoch),'obstime'],utc=True)
        ax.scatter(times,np.full(len(times),y),s=15,color='#145da0',alpha=.65,zorder=2)
    ax.set_yticks(range(len(timing)),[f'{r.nominal_epoch} | {r.epoch_start_utc[:10]}' for r in timing.itertuples()])
    ax.invert_yaxis();ax.set_ylabel('Calibration epoch MJD | start date (UTC)')
    ax.set_xlabel('Recorded exposure time (UTC)')
    ax.set_title(f'{len(mapping)} local exposures in time: calibration intervals and observed timestamps')
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3));ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.tick_params(axis='x',rotation=45);ax.grid(axis='x',alpha=.2)
    fig.tight_layout();fig.savefig(OUT/'exposure_epoch_timeline.png',dpi=150);plt.close(fig)
    display_timing=timing.copy()
    for column in ['epoch_start_utc','epoch_end_utc_exclusive']:
        display_timing[column]=display_timing[column].str[:10]
    labels={'nominal_epoch':'Calibration epoch start MJD','epoch_start_utc':'Epoch start date (UTC)',
        'epoch_end_utc_exclusive':'Next epoch start date (UTC)', 'exposures':'Moon/Sun-selected exposures (before Gaia)',
        'observing_nights':'Nights with Moon/Sun-selected exposures','first_exposure_utc':'First selected exposure (UTC)',
        'last_exposure_utc':'Last selected exposure (UTC)','full_moon_neighborhoods':'Distinct sampled lunations (nearest full Moon)',
        'mjd':'Observing-night MJD','full_moon_utc':'Nearest full Moon (UTC)'}
    return f"""<h2>Exposure times relative to calibration epochs</h2>
<p>Gray bars: calibration validity intervals [start, next start). Blue points: all selected exposure UTC timestamps, before Gaia filtering; nearby points can overlap.
Epochs are assigned by observing-night MJD. All nine camera calibration IDs agree with the configured assignments.</p>
<p><b>Selected exposures</b> pass Moon altitude &gt;60°, illumination &gt;95% and Sun altitude &lt;−18°; Gaia is checked later for each fiber.
<b>Sampled lunations</b> counts distinct nearest-full-Moon dates represented in an epoch. Epoch boundaries are shown as dates; the next epoch start is excluded from the current interval.</p>
<img src="exposure_epoch_timeline.png" alt="Actual exposure timestamps against calibration epoch intervals">
<p>Downloads: <a href="exposure_epoch_mapping.csv">every exposure and its epoch</a>, <a href="camera_epochs.json">camera calibration IDs</a>,
<a href="epoch_time_coverage.csv">epoch coverage</a>, <a href="observing_night_epoch_mapping.csv">observed nights</a>.
The tables show the {len(timing)} epochs with data; the inventory contains {len(starts)} configured epochs.</p>
<div class="scroll">{display_timing.rename(columns=labels).to_html(index=False)}</div>
<p>Supplementary phase bins group illumination and waxing/waning branch. Dates and calibration epochs are given by the timeline and night table.</p>
<h3>Every observed night</h3><div class="scroll">{nights.rename(columns=labels).to_html(index=False)}</div>
"""


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    frozen=pd.read_csv(RUN/'selected_exposures.csv')
    source_meta = frame(fits.getdata(DATA/'lvmsframe_median_stack_1.3.2_gaia1over100_meta_only.fits', 1))
    source_meta['source_path']=source_meta.path.map(lambda v:'/data/sas/sdsswork/'+str(v).split('/sdsswork/',1)[1])
    source_meta['geometry_selected']=(source_meta.moon_alt>60)&(source_meta.moon_fli>.95)&(source_meta.sun_alt<-18)
    source_meta['source_available']=pd.Series(pd.NA,index=source_meta.index,dtype='boolean')
    candidates=source_meta[source_meta.geometry_selected].copy()
    candidates['source_available']=candidates.source_path.map(lambda v:Path(v).is_file())
    source_meta.loc[candidates.index,'source_available']=candidates.source_available
    source_meta['completed_in_frozen_run']=source_meta.expnum.isin(frozen.expnum)
    fitted_ids=fitted_exposure_ids()
    source_meta['has_completed_fit']=source_meta.expnum.isin(fitted_ids)
    source_meta.to_csv(OUT/'exposure_selection_audit.csv',index=False)
    selected=candidates[candidates.source_available].copy().sort_values(['mjd','expnum'])
    selected['completed_in_frozen_run']=selected.expnum.isin(frozen.expnum)
    selected['has_completed_fit']=selected.expnum.isin(fitted_ids)
    assert selected.expnum.is_unique and set(frozen.expnum)<=set(selected.expnum)
    selection_counts=dict(catalog_exposures=len(source_meta),geometry_selected=len(candidates),
        local_sframe_exposures=len(selected),completed_individual_fit_exposures=len(fitted_ids),
        additional_local_expnums=selected.loc[~selected.completed_in_frozen_run,'expnum'].astype(int).tolist())
    (OUT/'exposure_selection_counts.json').write_text(json.dumps(selection_counts,indent=2))
    selected['observation_mjd_utc'] = Time(selected.obstime.to_numpy(dtype=str),format='isot').mjd
    # Manifest obstime is the recorded observation timestamp; retain that convention.
    t = selected.observation_mjd_utc.to_numpy()
    moons = full_moons(t.min()-35,t.max()+35)
    nearest = np.abs(t[:,None]-moons).argmin(axis=1)
    selected['full_moon_utc'] = Time(moons[nearest],format='mjd').isot
    selected['days_from_full_moon'] = t-moons[nearest]
    selected['phase_bin'] = np.where(selected.days_from_full_moon<0,'waxing','waning') + pd.cut(
        selected.moon_fli,bins=[.95,.97,.99,1.00001],labels=['_95-97%','_97-99%','_99-100%']).astype(str)
    definitions = yaml.safe_load(EPOCH_FILE.read_text())['epochs']
    starts = np.array(sorted(definitions))
    selected['nominal_epoch'] = starts[np.searchsorted(starts,selected.mjd,side='right')-1]
    photo_file = OUT/'photometry.parquet'
    camera_file = OUT/'camera_epochs.json'
    provenance_file=OUT/'gaia_photometry_provenance.json'
    provenance=photometry_provenance(selected)
    if (photo_file.exists() and camera_file.exists() and provenance_file.exists()
        and json.loads(provenance_file.read_text())==provenance):
        photos=pd.read_parquet(photo_file); cameras=json.loads(camera_file.read_text())
    else:
        tables=[]; cameras={}
        with ThreadPoolExecutor(max_workers=6) as pool:
            for i,(phot,epochs) in enumerate(pool.map(photometry,selected.to_dict('records'))):
                tables.append(phot); cameras[str(int(phot.expnum.iloc[0]))]=epochs
                if (i+1)%32==0:print(f'Photometry {i+1}/{len(selected)}',flush=True)
        photos=pd.concat(tables,ignore_index=True)
        photos.to_parquet(photo_file,index=False)
        camera_file.write_text(json.dumps(cameras,indent=2))
        provenance_file.write_text(json.dumps(provenance,indent=2))
    assert not photos.duplicated(['expnum','fiberid']).any()
    assert len(photos)==len(selected)*1920
    assert photos.loc[photos.gaia_verified,'gaia_cache'].map(lambda p:Path(p).resolve().is_relative_to(DATA.resolve())).all()
    pieces=[]
    for path in fit_files():
        with fits.open(path,memmap=True) as h:
            inp=frame(h['INPUT_META'].data)
            meta=frame(h['META'].data)
            quality=frame(h['LSF_META'].data)
            assert len(inp)==len(meta) and np.array_equal(meta.source_row,np.arange(len(inp)))
            inp['fit_status']=meta.fit_status.to_numpy()
            quality=quality.merge(inp[['expnum','fiberid','spectrographid','telescope','fibstatus',
                'fit_valid_fraction','fit_status']],left_on='spectrum_index',right_index=True,validate='many_to_one')
            flags=frame(Table(h['LSF_FAILURES'].data)[['spectrum_index','channel']])
            bad=set(zip(flags.spectrum_index,flags.channel))
            quality['lsf_flagged']=[(i,c) in bad for i,c in zip(quality.spectrum_index,quality.channel)]
            quality['quality_ok']=(quality.available & quality.converged & (quality.n_valid_pixels>0)
                & (quality.fit_status=='Solved') & (quality.fibstatus==0) & (quality.fit_valid_fraction>=.8)
                & ~quality.lsf_flagged & np.isfinite(quality.chi2_red) & (quality.chi2_red>=0)
                & np.isfinite(quality.sigma_pix_median) & (quality.sigma_pix_median>0))
            pieces.append(quality.drop(columns='spectrum_index'))
        print(f'Quality {path.name}',flush=True)
    quality=pd.concat(pieces,ignore_index=True).drop(columns=['telescope','spectrographid','fibstatus'])
    grid=photos.merge(pd.DataFrame({'channel':['B','R','Z']}),how='cross')
    rows=grid.merge(quality,on=['expnum','fiberid','channel'],how='left',validate='one_to_one')
    for col in ['quality_ok','available','converged','lsf_flagged']:rows[col]=rows[col].fillna(False).astype(bool)
    rows['fit_status']=rows.fit_status.fillna('not_fitted')
    rows=rows.merge(selected[['expnum','mjd','observation_mjd_utc','moon_fli','full_moon_utc',
        'days_from_full_moon','nominal_epoch','phase_bin']],on='expnum',validate='many_to_one')
    rows['camera']=rows.channel+rows.spectrographid.astype(str)
    rows['calibration_epoch']=[cameras[str(e)].get(c,-1) for e,c in zip(rows.expnum,rows.camera)]
    rows['accepted']=rows.gaia_clean_1pct & rows.quality_ok & (rows.calibration_epoch>=0)
    assert len(rows)==len(selected)*1920*3
    assert not rows.duplicated(['expnum','fiberid','channel']).any()
    assert (rows.loc[rows.accepted,'gaia_ratio']<.01).all()
    rows.to_parquet(OUT/'fiber_exposure_channel.parquet',index=False)
    keys=['calibration_epoch','full_moon_utc','fiberid','telescope','channel']
    coverage=summarize_groups(rows,keys)
    coverage.to_csv(OUT/'fiber_lunation_coverage.csv',index=False)
    phase=summarize_groups(rows,keys+['phase_bin'])
    phase.to_csv(OUT/'fiber_phase_coverage.csv',index=False)
    epochs=rows.groupby(['calibration_epoch','channel']).agg(selected_exposures=('expnum','nunique'),
        selected_fibers=('fiberid','nunique'),verified_rows=('gaia_verified','sum'),clean_rows=('gaia_clean_1pct','sum'),
        accepted_rows=('accepted','sum'),first_mjd=('mjd','min'),last_mjd=('mjd','max')).reset_index()
    good=rows[rows.accepted]
    support=good.groupby(['calibration_epoch','channel']).agg(accepted_exposures=('expnum','nunique'),
        accepted_fibers=('fiberid','nunique'),nights=('mjd','nunique'),lunations=('full_moon_utc','nunique')).reset_index()
    epochs=epochs.merge(support,how='left').fillna(0)
    repeated=coverage.groupby(['calibration_epoch','channel']).agg(repeat_exposure_fiber_lunations=('repeat_exposures','sum'),
        repeat_night_fiber_lunations=('repeat_nights','sum')).reset_index()
    epochs=epochs.merge(repeated)
    epochs['calibration_date_utc']=[Time(e,format='mjd').iso[:10] if e>=0 else 'unknown' for e in epochs.calibration_epoch]
    epochs['first_observation_utc']=[Time(e,format='mjd').iso[:10] for e in epochs.first_mjd]
    epochs['last_observation_utc']=[Time(e,format='mjd').iso[:10] for e in epochs.last_mjd]
    epochs.to_csv(OUT/'epoch_coverage.csv',index=False)
    selected['actual_lsf_epochs']=[','.join(map(str,sorted(set(cameras[str(e)].values())))) for e in selected.expnum]
    selected.to_csv(OUT/'exposure_epoch_mapping.csv',index=False)
    time_html=write_time_overview(selected)
    inventory=[]
    for i,e in enumerate(starts):
        subset=selected[selected.nominal_epoch==e]
        inventory.append(dict(epoch=int(e),date=Time(e,format='mjd').iso[:10],
            end_mjd=int(starts[i+1]) if i+1<len(starts) else None,trigger=definitions[int(e)]['trigger'],
            selected_exposures=len(subset),selected_lunations=subset.full_moon_utc.nunique()))
    inv=pd.DataFrame(inventory);inv.to_csv(OUT/'nominal_epoch_inventory.csv',index=False)
    payload=dict(exposures=len(selected),fiber_spectra=len(photos),date_start=selected.obstime.min(),
        date_end=selected.obstime.max(),gaia_cache_exposures=int(photos[photos.gaia_cache!='unavailable'].expnum.nunique()),
        verified_fiber_spectra=int(photos.gaia_verified.sum()),clean_1pct_fiber_spectra=int(photos.gaia_clean_1pct.sum()),
        clean_10pct_fiber_spectra=int(photos.gaia_clean_10pct.sum()),accepted_channel_rows=int(rows.accepted.sum()),
        full_moons=selected.full_moon_utc.nunique(),actual_epochs=rows.calibration_epoch.nunique(),
        accepted_epochs=good.calibration_epoch.nunique(),nominal_epochs_selected=selected.nominal_epoch.nunique(),
        total_nominal_epochs=len(starts),unknown_camera_rows=int((rows.calibration_epoch<0).sum()),
        gaia_source_directory=str(CACHE),geometry_selected_exposures=len(candidates),
        completed_fit_exposures=len(fitted_ids),phase_repeat_night_samples=int(phase.repeat_nights.sum()))
    (OUT/'summary.json').write_text(json.dumps(payload,indent=2))
    for channel in ['B','R','Z']:
        sub=coverage[coverage.channel==channel]
        grid=sub.pivot_table(index='fiberid',columns=['full_moon_utc','calibration_epoch'],values='n_nights',aggfunc='sum')
        verified_grid=sub.pivot_table(index='fiberid',columns=['full_moon_utc','calibration_epoch'],values='n_verified',aggfunc='sum').reindex_like(grid)
        colormap=plt.get_cmap('viridis').copy();colormap.set_bad('#cccccc')
        values=np.ma.masked_where(verified_grid.to_numpy()==0,grid.to_numpy())
        fig,ax=plt.subplots(figsize=(15,7))
        im=ax.imshow(values,aspect='auto',origin='lower',interpolation='nearest',vmin=0,cmap=colormap)
        ax.set_xticks(range(len(grid.columns)),[f'{d[:10]}\n{e}' for d,e in grid.columns],rotation=90,fontsize=8)
        ticks=np.linspace(0,len(grid)-1,9).astype(int)
        ax.set_yticks(ticks,grid.index[ticks]);ax.set_ylabel('Physical fiber ID')
        ax.set_title(f'{channel}: accepted observing nights per physical fiber, lunation and LSF calibration epoch')
        fig.colorbar(im,ax=ax,label='Accepted observing nights');fig.tight_layout()
        fig.savefig(OUT/f'coverage_{channel}.png',dpi=150);plt.close(fig)
    phase_summary=rows.groupby(['phase_bin','channel']).agg(selected_exposures=('expnum','nunique'),
        verified_fiber_spectra=('gaia_verified','sum'),accepted_channel_measurements=('accepted','sum')).reset_index()
    phase_good=good.groupby(['phase_bin','channel']).agg(accepted_exposures=('expnum','nunique'),
        accepted_fibers=('fiberid','nunique'),distinct_nights=('mjd','nunique')).reset_index()
    phase_summary=phase_summary.merge(phase_good,how='left').fillna(0)
    phase_counts=phase.groupby(['phase_bin','channel']).agg(measured_fiber_lunations=('measured','sum'),
        repeat_night_fiber_lunations=('repeat_nights','sum')).reset_index()
    phase_summary=phase_summary.merge(phase_counts)
    phase_summary.to_csv(OUT/'phase_summary.csv',index=False)
    phase_html=phase_summary.to_html(index=False,float_format=lambda x:f'{x:,.0f}')
    mismatches=int((rows.calibration_epoch!=rows.nominal_epoch).sum())
    table=epochs.to_html(index=False,float_format=lambda x:f'{x:,.0f}')
    content=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Individual-fiber LSF coverage</title>
<style>body{{font:16px/1.5 sans-serif;max-width:1400px;margin:30px auto;padding:0 20px;color:#17202a}}table{{border-collapse:collapse;font-size:12px}}td,th{{padding:5px;border:1px solid #ccd}}.scroll{{overflow:auto}}img{{width:100%}}code{{background:#eee}}</style>
<h1>Individual-fiber LSF coverage at bright Moon</h1><p>LSF is the line-spread function, the instrumental spectral response. B, R and Z are the three spectral channels.
SFrames are the reduced flux-calibrated exposure spectra. Physical fiber ID identifies a fixed instrument fiber across exposures.</p>
<p>Corrected Gaia audit: {payload['geometry_selected_exposures']} geometry-selected exposures, {payload['exposures']} local SFrames, {payload['completed_fit_exposures']} exposures with completed individual fits, {payload['fiber_spectra']:,} individual fiber spectra,
{payload['date_start'][:10]} to {payload['date_end'][:10]} UTC, sampling {payload['full_moons']} full-Moon neighborhoods. Science, SkyE and SkyW fibers are included; physical fiber IDs remain separate.</p>
{time_html}<h2>Corrected source and complete exposure selection</h2><p>This revision replaces the earlier Gaia audit that used a different directory.
The lunar/geometry conditions are applied to every exposure in the supplied DRP 1.3.2 metadata before the individual-fiber Gaia cut.
<a href="exposure_selection_audit.csv">Full exposure audit</a> records each step; all local opportunities are represented even when an individual fit or Gaia verification is unavailable.</p><h2>Measured support</h2><p>Eligibility counts below preserve the user-requested frozen Gaia &lt;1% audit. The chronology above includes every selected exposure regardless of this eligibility.</p><p>Accepted LSF measurements exist in {payload['accepted_epochs']} of {payload['actual_epochs']} actual epochs in this run.
Local Gaia aperture tables cover {payload['gaia_cache_exposures']}/{payload['exposures']} exposures.
The tables give the corrected per-fiber coverage, including the 2026 epochs.
Single-exposure samples measure an LSF with a single-night baseline; repetition requires multiple accepted measurements.</p>
<h2>Selection and interpretation</h2><p>The exposure selection is Moon altitude &gt;60°, solar altitude &lt;−18°, and illuminated fraction &gt;0.95.
Coverage therefore describes the vicinity of full Moon. Illuminated fraction alone does not identify waxing versus waning;
phase tables split these branches and use 95–97%, 97–99% and 99–100% illumination bins.</p>
<p>Gaia eligibility requires the summed catalog-source Gaia G flux in the 17.65 arcsec fiber aperture to be strictly less than
1/100 of the synthetic G flux of the current DRP 1.3.2 <code>FLUX+SKY</code> spectrum. Both use the stored Gaia G passband calibration.
The synthetic flux is a photon-weighted mean flux density, with inverse-variance and mask validity and at least 95% of the G-passband weight on the available SFrame wavelength grid.
This validity fraction describes masked-pixel coverage on that grid; it is distinct from the full passband wavelength coverage.
Gaia cache coordinates must agree within 0.05 arcsec; all matched catalog sources must have valid photometry.
A verified aperture containing zero catalog sources has zero summed catalog flux. Unverified apertures have unknown contamination.</p>
<p>The requested directory supplies the DRP 1.3.2 median-stack exposure metadata. Its median spectra do not identify individual-fiber contamination.
Detailed local Gaia aperture tables come exclusively from <code>{CACHE}</code>, the user-specified
<a href="http://lumen.abudhabi.nyu.edu/sas/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/">Lumen Gaia directory</a>.
The <a href="gaia_photometry_provenance.json">source manifest</a> records file paths and SHA-256 checksums. Their Gaia numerator is reused after coordinate and passband checks;
the denominator and its uncertainty are recomputed from each current SFrame, including FLUX and SKY variance.</p>
<p>Verified Gaia measurements: {payload['verified_fiber_spectra']:,}/{payload['fiber_spectra']:,} fiber spectra;
strict 1% selection: {payload['clean_1pct_fiber_spectra']:,}. A separate 10% flag ({payload['clean_10pct_fiber_spectra']:,}) is retained for reproducibility.</p>
<p>An accepted channel measurement additionally requires a good hardware fiber, main fit Solved, channel available and converged,
positive valid-pixel count and width, finite nonnegative recorded reduced chi-square, input-valid fraction ≥0.8,
and a known LSF calibration epoch. Every channel entry in LSF_FAILURES is filtered, including AlmostSolved and fallback.
These are technical quality checks; reduced chi-square has no calibrated scientific acceptance threshold here.
Accepted channel measurements: {payload['accepted_channel_rows']:,}. Each individual fiber spectrum can contribute up to three B/R/Z measurements;
fiber-spectrum counts and channel-measurement counts therefore have different units.</p>
<h2>Observation dates and calibration epochs</h2><p>The epoch in the tables is the actual <code>lvm-mlsf-[camera]</code> calibration directory recorded in each SFrame COMMENT.
It is read separately for B/R/Z and spectrographs 1/2/3. The nominal epoch is the latest start MJD not later than the observing-night MJD,
using the local lvmcore calibration-epochs.yaml. Actual and nominal epochs have {mismatches} mismatched channel rows across {len(selected)} local exposures and nine cameras. The exposure mapping preserves both. Calendar dates are UTC conversions of epoch MJD;
observing-night MJD and observation timestamp are retained separately.
The observation_mjd_utc field converts the recorded obstime timestamp directly; no exposure-midpoint adjustment is applied.
Distinct observing-night counts describe sampling, without asserting statistical independence.</p>
<p>A lunation is indexed by its nearest computed full Moon (built-in Astropy ephemeris). Signed elapsed days preserve the waxing/waning branch.
Each lunation is split at actual calibration changes: measurements from different epochs are never treated as one calibration sample.</p>
<h2>Coverage metrics and practical use</h2><p>For each epoch × lunation × physical fiber × channel, CSV tables give selected, verified, clean and accepted exposure counts;
distinct observing nights; exposure hours; time span and largest gap; illuminated-fraction and signed-day ranges;
median synthetic-G S/N and recorded LSF reduced chi-square. Medians use accepted rows only.</p>
<p>Effective nights = (Σ night exposure time)² / Σ (night exposure time)². Largest-night share measures concentration.
One accepted exposure establishes a measurement; ≥3 accepted exposures support within-sample repetition;
≥2 nights support a cross-night comparison. These count thresholds are descriptive, not precision guarantees.
S/N in broadband G is not the information content of the solar or sky lines. Measured LSF uncertainty and repeatability are required
before claiming a given resolution accuracy or transferring LSF to another phase.</p>
<p>Phase-table denominators are the selected exposures in each populated bin. Unobserved bins have zero coverage; they are absent from the sparse CSV.
The frozen run samples full-Moon conditions and supports no inference of coverage at lower illumination.
The nominal inventory lists all {len(starts)} configured epochs, including epochs with zero selected exposures.</p>
<h2>Supplementary illumination-bin coverage</h2><p>The corrected selection contains {payload['phase_repeat_night_samples']:,} fiber × lunation × illumination-bin × channel samples with accepted measurements on at least two observing nights.
Phase summaries cover all selected lunations together. Distinct nights count nights with at least one accepted fiber;
a given fiber may have fewer. Use fiber_phase_coverage.csv for the per-fiber, per-lunation values.</p><div class="scroll">{phase_html}</div>
<h2>Actual LSF epoch summary</h2><p>Counts below combine fibers for inventory only. The CSVs preserve individual fibers.
Repeat-night entries count fiber–lunation pairs; repeat-exposure entries also count fiber–lunation pairs.</p><div class="scroll">{table}</div>
<h2>Distinct-night coverage</h2><p>The first line of each column label is the full-Moon UTC date; the second is the actual calibration epoch MJD. Rows identify individual physical fibers.
Color is the number of accepted observing nights. Gray cells have selected exposures and zero verified Gaia measurements, so contamination is unknown.
Purple zero-valued cells have verified Gaia measurements and zero accepted nights.</p>
{''.join(f'<img src="coverage_{c}.png" alt="{c} channel distinct observing nights by physical fiber and lunation">' for c in ['B','R','Z'])}
<h2>Reproducible data</h2><ul>
{''.join(f'<li><a href="{p}">{p}</a></li>' for p in ['epoch_time_coverage.csv','observing_night_epoch_mapping.csv','epoch_coverage.csv','phase_summary.csv','exposure_epoch_mapping.csv','nominal_epoch_inventory.csv','fiber_lunation_coverage.csv','fiber_phase_coverage.csv','fiber_exposure_channel.parquet','photometry.parquet','camera_epochs.json','gaia_photometry_provenance.json','exposure_selection_counts.json','exposure_selection_audit.csv','summary.json'])}</ul>
<p>Build script: <code>analysis/lsf_epoch_coverage/build_coverage.py</code>. Production fits and source catalogs are read only.
Coverage is scoped to this completed run, rather than all exposures in the survey.</p></html>'''
    (OUT/'lsf_epoch_coverage.html').write_text(content)
    print(json.dumps(payload,indent=2),flush=True)


if __name__=='__main__':main()
