"""Measure LSF repeatability within each actual epoch, lunation and physical fiber."""
from pathlib import Path
from types import SimpleNamespace
import json
import hashlib
import sys

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from skysub.sky_decomp.lsf_spline2d import (
    mspline_basis, mspline_bin_integrals, _offset_knots, evaluate_lsf_diagnostics,
    evaluate_bspline_basis, LSF_SPLINE2D_REPRESENTATION)
from analysis.solarblue_lsf_serious_corpus.analyze_lsf import ANCHORS, CENTRAL_ANCHOR
from analysis.lsf_epoch_coverage.build_coverage import RUN, OUT, frame, fit_files

DEST=OUT/'stability'
KEYS=['calibration_epoch','full_moon_utc','fiberid','telescope','channel','wavelength_angstrom']
METRICS=['fwhm_angstrom','w50_angstrom','w90_angstrom','centroid_angstrom']


def diagnostics(masses, half_width):
    """Batch the existing 1201-point diagnostic definition without changing its grid."""
    count=masses.shape[1]
    delta=np.linspace(-half_width,half_width,1201)
    knots=_offset_knots(count,half_width)
    centers=.5*(knots[:-4]+knots[4:])
    variance=((knots[-1]-knots[0])/(len(knots)-1))**2/3
    cdf=masses@mspline_bin_integrals(np.full(delta.shape,-half_width),delta,
                                   n_basis=count,half_width=half_width).T
    quantiles=[]
    for q in [.05,.25,.5,.75,.95]:
        hi=np.argmax(cdf>=q,axis=1);lo=np.maximum(hi-1,0)
        lower=cdf[np.arange(len(cdf)),lo];upper=cdf[np.arange(len(cdf)),hi]
        fraction=np.divide(q-lower,upper-lower,out=np.zeros(len(cdf)),where=upper>lower)
        quantiles.append(delta[lo]+fraction*(delta[hi]-delta[lo]))
    q=np.asarray(quantiles).T
    density=masses@mspline_basis(delta,n_basis=count,half_width=half_width).T
    peak=np.argmax(density,axis=1);half=density[np.arange(len(density)),peak]/2
    indices=np.arange(len(delta))[None,:]
    left=np.max(np.where((indices<=peak[:,None])&(density<=half[:,None]),indices,-1),axis=1)
    right=np.min(np.where((indices>=peak[:,None])&(density<=half[:,None]),indices,len(delta)),axis=1)
    fwhm=np.full(len(masses),np.nan)
    ok=(left>=0)&(right<len(delta))
    rows=np.flatnonzero(ok);l=left[ok];r=right[ok]
    yl=density[rows,l];yh=density[rows,l+1]
    xl=delta[l]+(half[ok]-yl)/(yh-yl)*(delta[l+1]-delta[l])
    yl=density[rows,r];yh=density[rows,r-1]
    xr=delta[r]+(half[ok]-yl)/(yh-yl)*(delta[r-1]-delta[r])
    fwhm[ok]=xr-xl
    centroid=masses@centers
    return dict(fwhm_angstrom=fwhm,centroid_angstrom=centroid,
        sigma_angstrom=np.sqrt(np.maximum(masses@(centers**2+variance)-centroid**2,0)),
        w50_angstrom=q[:,3]-q[:,1],w90_angstrom=q[:,4]-q[:,0],
        asymmetry=(q[:,4]+q[:,0]-2*q[:,2])/(q[:,4]-q[:,0]))


def extract():
    selection=pd.read_parquet(OUT/'fiber_exposure_channel.parquet')
    selection=selection[selection.accepted].copy()
    wanted=selection.set_index(['expnum','fiberid','channel'])
    results=[];audit=[]
    for p in fit_files():
        with fits.open(p,memmap=True) as h:
            inp=frame(h['INPUT_META'].data)
            coef=h['LSF_COEF'].data;knots=h['LSF_KNOTS'].data
            config=h['LSF_CONFIG'].data
            for ci,channel in enumerate(['B','R','Z']):
                keys=pd.MultiIndex.from_arrays([inp.expnum,inp.fiberid,np.full(len(inp),channel)],
                                              names=['expnum','fiberid','channel'])
                take=np.flatnonzero(keys.isin(wanted.index))
                if not len(take):continue
                cf=config[ci]
                assert cf['channel'].strip()==channel and cf['legacy_kernel_representation'].strip()==LSF_SPLINE2D_REPRESENTATION
                basis_count=int(cf['n_basis']);degree=int(cf['degree']);half=float(cf['offset_half_width_angstrom'])
                knot=knots[take[0],ci,:int(cf['n_knots'])].astype(float)
                assert np.all(knots[take,ci,:len(knot)]==knot)
                wavelengths=np.array(ANCHORS[channel])
                assert wavelengths.min()>=float(cf['fit_lower']) and wavelengths.max()<=float(cf['fit_upper'])
                basis=evaluate_bspline_basis(wavelengths,knot,degree)
                for ids in np.array_split(take,max(1,(len(take)+1023)//1024)):
                    coefficients=np.asarray(coef[ids,ci,:,:basis_count],dtype=float)
                    masses=np.einsum('wb,ntb->nwt',basis,coefficients)
                    assert np.isfinite(masses).all() and masses.min()>=-1e-7
                    assert np.allclose(masses.sum(axis=2),1,atol=2e-5)
                    values=diagnostics(masses.reshape(-1,masses.shape[-1]),half)
                    meta=wanted.loc[keys[ids]].reset_index().loc[np.repeat(np.arange(len(ids)),len(wavelengths))].reset_index(drop=True)
                    meta['wavelength_angstrom']=np.tile(wavelengths,len(ids))
                    for name,value in values.items():meta[name]=value
                    flat_masses=masses.reshape(-1,masses.shape[-1])
                    for component in range(flat_masses.shape[1]):meta[f'lsf_mass_{component:02d}']=flat_masses[:,component]
                    meta['offset_half_width_angstrom']=half
                    meta['diagnostic_valid']=(np.isfinite(meta[METRICS]).all(axis=1)
                        & (meta[['fwhm_angstrom','w50_angstrom','w90_angstrom']]>0).all(axis=1))
                    results.append(meta)
                # Check the batch path against the existing scientific evaluator on stored states.
                for index in take[[0,len(take)//2,-1]]:
                    state=SimpleNamespace(coefficients={},knot_vectors={},degrees={},
                        config={'offset_half_width_angstrom':half},legacy_kernel_representation=LSF_SPLINE2D_REPRESENTATION)
                    for c,j in zip(['B','R','Z'],range(3)):
                        state.coefficients[c]=np.asarray(coef[index,j],dtype=float)
                        state.knot_vectors[c]=np.asarray(knots[index,j],dtype=float)
                        state.degrees[c]=int(config[j]['degree'])
                    reference=evaluate_lsf_diagnostics(state,wavelengths)
                    mass=basis@state.coefficients[channel].T
                    actual=diagnostics(mass,half)
                    error=max(float(np.max(np.abs(actual[k]-reference[k]))) for k in actual)
                    assert error<1e-9,(p.name,channel,error)
                    audit.append(dict(shard=p.name,channel=channel,row=int(index),max_abs_error=error))
        print(f'Diagnostics {p.name}',flush=True)
    metrics=pd.concat(results,ignore_index=True)
    assert len(metrics)==len(selection)*5
    metrics.to_parquet(DEST/'accepted_lsf_diagnostics.parquet',index=False)
    pd.DataFrame(audit).to_csv(DEST/'diagnostic_equivalence_checks.csv',index=False)
    return metrics


def scatter(data,keys):
    """One empirical scatter estimate per physical fiber and wavelength, no pooling of raw fibers."""
    base=data.groupby(keys).agg(n_exposures=('expnum','size'),n_nights=('mjd','nunique'),
        first_mjd=('observation_mjd_utc','min'),last_mjd=('observation_mjd_utc','max'),
        median_g_snr=('g_snr','median'),median_chi2=('chi2_red','median'))
    center=data.groupby(keys)[METRICS].median()
    dev=(data[METRICS]-data.groupby(keys)[METRICS].transform('median')).abs()
    joined=data[keys].copy()
    for c in METRICS:joined[c]=dev[c].to_numpy()
    meddev=joined.groupby(keys)[METRICS].median()*1.4826
    p90=joined.groupby(keys)[METRICS].quantile(.9)
    result=base.copy()
    for c in METRICS:
        result[c+'_median']=center[c]
        result[c+'_rmad']=meddev[c]
        result[c+'_p90_abs_deviation']=p90[c]
        if c!='centroid_angstrom':
            result[c+'_rmad_pct']=100*meddev[c]/center[c]
            result[c+'_p90_abs_deviation_pct']=100*p90[c]/center[c]
    result['time_span_days']=result.last_mjd-result.first_mjd
    return result.reset_index()


def main():
    DEST.mkdir(exist_ok=True)
    diagnostic_file=DEST/'accepted_lsf_diagnostics.parquet'
    with (OUT/'fiber_exposure_channel.parquet').open('rb') as stream:
        selection_digest=hashlib.file_digest(stream,'sha256').hexdigest()
    provenance=dict(selection_sha256=selection_digest,diagnostics_version='correct_gaia132_with_profile_masses',
        fitted_files=[dict(path=str(p),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns)
        for p in fit_files()])
    provenance_file=DEST/'diagnostic_provenance.json'
    if diagnostic_file.exists() and provenance_file.exists() and json.loads(provenance_file.read_text())==provenance:
        raw=pd.read_parquet(diagnostic_file)
    else:
        raw=extract();provenance_file.write_text(json.dumps(provenance,indent=2))
    valid=raw[raw.diagnostic_valid].copy()
    valid['g_snr']=valid.synthetic_g/valid.synthetic_g_error
    summary=scatter(valid,KEYS)
    summary.to_csv(DEST/'fiber_lunation_stability.csv',index=False)
    nightly=scatter(valid,KEYS+['mjd'])
    nightly.to_csv(DEST/'fiber_night_repeatability.csv',index=False)
    repeat=nightly[nightly.n_exposures>=3].copy()
    within=repeat.groupby(KEYS).agg(n_repeat_nights=('mjd','size'),
        within_night_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct','median'),
        within_night_w90_rmad_pct=('w90_angstrom_rmad_pct','median'),
        within_night_centroid_rmad_angstrom=('centroid_angstrom_rmad','median')).reset_index()
    # Equal weight for each sufficiently repeated night, regardless of exposure count.
    for c in METRICS:repeat[c]=repeat[c+'_median']
    dev=(repeat[METRICS]-repeat.groupby(KEYS)[METRICS].transform('median')).abs()
    repeated_nights=repeat[KEYS].copy()
    repeated_nights['between_night_fwhm_abs_pct']=100*dev.fwhm_angstrom.to_numpy()/repeat.groupby(KEYS).fwhm_angstrom.transform('median').to_numpy()
    repeated_nights['between_night_w90_abs_pct']=100*dev.w90_angstrom.to_numpy()/repeat.groupby(KEYS).w90_angstrom.transform('median').to_numpy()
    repeated_nights['between_night_centroid_abs_angstrom']=dev.centroid_angstrom.to_numpy()
    between=repeated_nights.groupby(KEYS).quantile(.9).reset_index()
    between=between.merge(within[KEYS+['n_repeat_nights']],on=KEYS)
    between=between[between.n_repeat_nights>=2]
    between.to_csv(DEST/'fiber_between_night_stability.csv',index=False)
    comparisons=repeat.sort_values('first_mjd').copy()
    grouped=comparisons.groupby(KEYS)
    comparisons['previous_night_mjd']=grouped.mjd.shift()
    comparisons['night_gap_days']=comparisons.mjd-comparisons.previous_night_mjd
    for c in METRICS:
        prev=grouped[c].shift()
        comparisons[c+'_night_change']=comparisons[c]-prev
        if c!='centroid_angstrom':comparisons[c+'_night_change_pct']=100*(comparisons[c]/prev-1)
    comparisons=comparisons[comparisons.previous_night_mjd.notna()]
    comparisons.to_csv(DEST/'consecutive_night_changes.csv',index=False)
    comparisons['abs_fwhm_night_change_pct']=comparisons.fwhm_angstrom_night_change_pct.abs()
    full_change=comparisons.groupby(KEYS).agg(n_night_pairs=('mjd','size'),
        median_abs_full_night_change_pct=('abs_fwhm_night_change_pct','median')).reset_index()
    full_change.to_csv(DEST/'fiber_full_night_changes.csv',index=False)
    samples=summary[summary.n_exposures>=3].merge(within,on=KEYS,how='left')
    samples=samples.merge(between.drop(columns='n_repeat_nights'),on=KEYS,how='left')
    samples=samples.merge(full_change,on=KEYS,how='left')
    central=samples[np.array([w==CENTRAL_ANCHOR[c] for w,c in zip(samples.wavelength_angstrom,samples.channel)])]
    epoch_keys=['calibration_epoch','full_moon_utc','channel']
    epoch=central.groupby(epoch_keys).agg(repeated_fibers=('fiberid','size'),
        median_exposures_per_fiber=('n_exposures','median'),
        median_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct','median'),
        p90_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct',lambda x:x.quantile(.9)),
        median_within_night_rmad_pct=('within_night_fwhm_rmad_pct','median'),
        cross_night_fibers=('between_night_fwhm_abs_pct','count'),
        median_between_night_p90_pct=('between_night_fwhm_abs_pct','median'),
        p90_between_night_p90_pct=('between_night_fwhm_abs_pct',lambda x:x.quantile(.9)),
        median_centroid_rmad_angstrom=('centroid_angstrom_rmad','median')).reset_index()
    epoch.to_csv(DEST/'epoch_lunation_stability.csv',index=False)
    spectral=samples.groupby(['channel','wavelength_angstrom']).agg(fiber_lunation_samples=('fiberid','size'),
        median_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct','median'),
        p90_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct',lambda x:x.quantile(.9)),
        median_w90_rmad_pct=('w90_angstrom_rmad_pct','median'),
        median_centroid_rmad_angstrom=('centroid_angstrom_rmad','median'),
        cross_night_samples=('between_night_fwhm_abs_pct','count'),
        median_between_night_p90_pct=('between_night_fwhm_abs_pct','median'),
        p90_between_night_p90_pct=('between_night_fwhm_abs_pct',lambda x:x.quantile(.9))).reset_index()
    spectral.to_csv(DEST/'wavelength_stability.csv',index=False)
    channel_summary=central.groupby('channel').agg(fiber_lunation_samples=('fiberid','size'),
        median_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct','median'),
        p90_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct',lambda x:x.quantile(.9)),
        median_within_night_rmad_pct=('within_night_fwhm_rmad_pct','median'),
        cross_night_samples=('median_abs_full_night_change_pct','count'),
        median_abs_full_night_change_pct=('median_abs_full_night_change_pct','median'),
        p90_abs_full_night_change_pct=('median_abs_full_night_change_pct',lambda x:x.quantile(.9))).reset_index()
    channel_summary.to_csv(DEST/'channel_summary.csv',index=False)
    telescope_summary=central.groupby(['telescope','channel']).agg(fiber_lunation_samples=('fiberid','size'),
        median_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct','median'),
        p90_fwhm_rmad_pct=('fwhm_angstrom_rmad_pct',lambda x:x.quantile(.9)),
        cross_night_samples=('median_abs_full_night_change_pct','count'),
        median_abs_full_night_change_pct=('median_abs_full_night_change_pct','median'),
        p90_abs_full_night_change_pct=('median_abs_full_night_change_pct',lambda x:x.quantile(.9))).reset_index()
    telescope_summary.to_csv(DEST/'telescope_channel_summary.csv',index=False)
    tails=central[central.between_night_fwhm_abs_pct>10].copy()
    tails.to_csv(DEST/'central_night_departure_review_gt10pct.csv',index=False)
    tail_nights=nightly.merge(tails[KEYS],on=KEYS,how='inner')
    tail_nights.to_csv(DEST/'review_flagged_night_diagnostics.csv',index=False)
    tail_counts=tails.groupby(['calibration_epoch','full_moon_utc','telescope','channel']).size().rename('review_fiber_samples').reset_index()
    tail_counts.to_csv(DEST/'review_flag_counts.csv',index=False)
    for channel in ['B','R','Z']:
        sub=central[central.channel.eq(channel)]
        grid=sub.pivot(index='fiberid',columns=['full_moon_utc','calibration_epoch'],values='fwhm_angstrom_rmad_pct')
        fig,ax=plt.subplots(figsize=(14,7))
        cmap=plt.get_cmap('magma').copy();cmap.set_bad('#ccc')
        values=grid.to_numpy();positive=values[values>0]
        from matplotlib.colors import LogNorm
        floor=min(.1,positive.min());display=np.where(values==0,floor,values)
        image=ax.imshow(display,aspect='auto',origin='lower',interpolation='nearest',cmap=cmap,
                        norm=LogNorm(vmin=min(.1,positive.min()),vmax=positive.max()))
        ax.set_xticks(range(len(grid.columns)),[f'{d[:10]}\n{e}' for d,e in grid.columns],rotation=90,fontsize=8)
        ticks=np.linspace(0,len(grid)-1,9).astype(int);ax.set_yticks(ticks,grid.index[ticks])
        ax.set_ylabel('Physical fiber ID');ax.set_title(f'{channel} at {CENTRAL_ANCHOR[channel]:.0f} Å: within-lunation FWHM robust scatter')
        fig.colorbar(image,ax=ax,label='1.4826 × MAD / median FWHM (%) — logarithmic color scale')
        fig.tight_layout();fig.savefig(DEST/f'stability_{channel}.png',dpi=150);plt.close(fig)
    # Show individual-fiber examples from the largest and smallest B scatter groups.
    b=epoch[epoch.channel.eq('B')].sort_values('median_fwhm_rmad_pct')
    choices=pd.concat([b.head(2),b.tail(2)]).drop_duplicates(epoch_keys)
    anomaly=epoch[epoch.calibration_epoch.eq(60380)&epoch.channel.eq('B')].copy()
    choices=pd.concat([choices,anomaly],ignore_index=True)
    example_html=[];example_records=[]
    for i,row in enumerate(choices.itertuples()):
        fibers=central[(central.calibration_epoch==row.calibration_epoch)&central.full_moon_utc.eq(row.full_moon_utc)
            & central.channel.eq('B')&central.telescope.eq('Sci')].copy()
        fibers['distance']=(fibers.fwhm_angstrom_rmad_pct-fibers.fwhm_angstrom_rmad_pct.median()).abs()
        fid=181 if i==len(choices)-1 else int(fibers.sort_values(['distance','fiberid']).fiberid.iloc[0])
        plotted=valid[(valid.calibration_epoch==row.calibration_epoch)&valid.full_moon_utc.eq(row.full_moon_utc)&valid.fiberid.eq(fid)]
        fig,axes=plt.subplots(4,3,figsize=(14,10),sharex=True,sharey='row')
        plotted_metrics=[('fwhm_angstrom','ΔFWHM (%)',1.),('w50_angstrom','ΔW50 (%)',1.),
                         ('w90_angstrom','ΔW90 (%)',1.),('centroid_angstrom','ΔCentroid (Å)',.01)]
        for col,channel in enumerate(['B','R','Z']):
            sub=plotted[plotted.channel.eq(channel)].sort_values('observation_mjd_utc')
            w=CENTRAL_ANCHOR[channel]
            points=sub[sub.wavelength_angstrom.eq(w)].copy()
            axes[0,col].set_title(f'{channel} at {w:.0f} Å')
            for metric_row,(metric,label,threshold) in enumerate(plotted_metrics):
                ax=axes[metric_row,col]
                normalization=points[metric].median()
                points['change']=points[metric]-normalization if metric=='centroid_angstrom' else 100*(points[metric]/normalization-1)
                artist=ax.plot(points.days_from_full_moon,points.change,'.',markersize=5)[0]
                for _,night in points.groupby('mjd'):
                    ax.plot(night.days_from_full_moon,night.change,'-',lw=.7,color=artist.get_color())
                ax.axhline(0,color='black',lw=.5);ax.set_yscale('symlog',linthresh=threshold)
                if col==0:ax.set_ylabel(label)
            axes[-1,col].set_xlabel('Days from full Moon (UTC)')
        count=plotted[plotted.channel.eq('B')&plotted.wavelength_angstrom.eq(CENTRAL_ANCHOR['B'])]
        fig.suptitle(f'Physical fiber {fid} • calibration epoch {row.calibration_epoch} • full Moon {row.full_moon_utc[:10]}\nB coverage: {count.expnum.nunique()} exposures across {count.mjd.nunique()} observing nights')
        fig.tight_layout();name=f'fiber_example_{i+1}.png';fig.savefig(DEST/name,dpi=150);plt.close(fig)
        example_html.append(f'<img src="{name}" alt="Separate B/R/Z panels for physical fiber {fid}">')
        example_records.append(dict(calibration_epoch=row.calibration_epoch,full_moon_utc=row.full_moon_utc,
            fiberid=fid,selection='SkyW fiber 181 illustrating a large inter-night LSF anomaly' if i==len(choices)-1 else 'Science fiber nearest median Science B robust scatter in this group'))
    pd.DataFrame(example_records).to_csv(DEST/'example_fibers.csv',index=False)
    payload=dict(accepted_channel_measurements=int(raw[['expnum','fiberid','channel']].drop_duplicates().shape[0]),
        diagnostic_rows=len(raw),valid_diagnostic_rows=len(valid),invalid_diagnostic_rows=int((~raw.diagnostic_valid).sum()),
        repeated_fiber_lunation_channels=len(central),cross_night_fiber_lunation_channels=int(central.between_night_fwhm_abs_pct.notna().sum()),
        epochs=int(epoch.calibration_epoch.nunique()),lunations=int(epoch.full_moon_utc.nunique()),
        diagnostic_equivalence_checks=len(pd.read_csv(DEST/'diagnostic_equivalence_checks.csv')))
    (DEST/'summary.json').write_text(json.dumps(payload,indent=2))
    labels={'calibration_epoch':'Calibration epoch (MJD)','full_moon_utc':'Full Moon (UTC)',
        'channel':'Channel','telescope':'Fiber telescope','wavelength_angstrom':'Wavelength (Å)',
        'repeated_fibers':'Repeated fiber samples','median_exposures_per_fiber':'Median exposures per fiber',
        'fiber_lunation_samples':'Fiber–lunation samples','median_fwhm_rmad_pct':'Median FWHM scatter (%)',
        'p90_fwhm_rmad_pct':'P90 FWHM scatter (%)','median_within_night_rmad_pct':'Median within-night scatter (%)',
        'cross_night_fibers':'Cross-night fiber samples','cross_night_samples':'Cross-night samples',
        'median_between_night_p90_pct':'Median night P90 departure (%)',
        'p90_between_night_p90_pct':'P90 across fibers: night P90 departure (%)',
        'median_centroid_rmad_angstrom':'Median centroid scatter (Å)',
        'median_w90_rmad_pct':'Median W90 scatter (%)',
        'median_abs_full_night_change_pct':'Median full sampled-night change (%)',
        'p90_abs_full_night_change_pct':'P90 full sampled-night change (%)',
        'review_fiber_samples':'Review-flagged fiber samples'}
    def table(df):return df.rename(columns=labels).to_html(index=False,float_format=lambda x:f'{x:.3f}')
    html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>LSF stability within calibration epochs and lunations</title>
<style>body{{font:16px/1.5 sans-serif;max-width:1400px;margin:30px auto;padding:0 20px;color:#17202a}}table{{border-collapse:collapse;font-size:12px}}td,th{{padding:5px;border:1px solid #ccc}}.scroll{{overflow:auto}}img{{width:100%}}</style>
<h1>Individual-fiber line-spread function (LSF): stability within epochs and lunations</h1>
<p>Repeatability is measured separately for each physical fiber, calibration epoch, lunation, channel and wavelength.
A lunation labels the sampled bright-Moon window. Selection: verified Gaia G-band flux / synthetic observed G-band flux from FLUX+SKY &lt;0.01 and the fit-quality checks in the <a href="../lsf_epoch_coverage.html">main report</a>.</p>
<p><b>{payload['accepted_channel_measurements']:,} accepted measurements</b>: one exposure × physical fiber × channel (B/R/Z), evaluated at five wavelengths per channel.
{payload['repeated_fiber_lunation_channels']:,} fiber–lunation–channel samples have at least three valid exposures;
{payload['cross_night_fiber_lunation_channels']:,} support a comparison of at least two nights with ≥3 exposures each.
The repeated samples cover {payload['epochs']} calibration epochs and {payload['lunations']} lunations. Acceptance establishes fit quality; precision requires uncertainty validation.</p>
<h2>Metrics and interpretation</h2><p>FWHM is the literal full width at half maximum; W50 and W90 enclose the central 50% and 90% of the LSF area.
The centroid is the mean wavelength offset in Å. The model's stored M-spline basis, wavelength knots and ±2.5 Å offset support are used directly,
with the same 1201-point diagnostic grid as the existing scientific evaluator. Batch results are checked against that evaluator on
{payload['diagnostic_equivalence_checks']} stored channel states to absolute numerical tolerance 10⁻⁹.
Diagnostic rows requiring finite positive widths: {payload['valid_diagnostic_rows']:,}/{payload['diagnostic_rows']:,}.</p>
<ul><li><b>Within-lunation scatter:</b> ≥3 exposures per fiber; widths use 100 × 1.4826 × MAD / median width (%), centroid uses 1.4826 × MAD (Å). MAD is median absolute deviation.</li>
<li><b>Within-night scatter:</b> the same calculation for nights with ≥3 exposures; summarize qualifying nights by their median for each fiber.</li>
<li><b>Between-night departure:</b> ≥2 qualifying nights; P90 over nights of 100 × |night median − median of night medians| / median of night medians.</li>
<li><b>Sampled-night change:</b> full percentage change between adjacent observed night medians, relative to the earlier night. Gaps can exceed one day.</li></ul>
<p>Summary medians and P90 (90th percentile) are calculated across individual-fiber statistics, with equal weight per fiber–lunation sample.
The between-night table's outer P90 summarizes each fiber's P90 departure; this differs from a full night-to-night change.
Centroid changes use Å. The factor 1.4826 scales MAD to a Gaussian-equivalent spread; with three exposures, MAD can underrepresent an outlier.</p>
<p>These statistics measure recovered-LSF repeatability, including fit noise and observing conditions. Their populations differ, so subtracting them does not isolate hardware drift.
This bright-Moon sample does not establish stability at lower illumination. G-band S/N and chi-square remain accompanying diagnostics.</p>
<h2>Central-wavelength overview</h2><p>Full change uses adjacent sampled qualifying night medians, relative to the earlier night.
For each fiber–lunation sample, the median absolute full change is computed over its adjacent sampled-night pairs (gaps can exceed one day); the table summarizes those per-fiber medians.
This full change differs from the between-night departure about a central median.</p><div class="scroll">{table(channel_summary)}</div>
<h2>Central-wavelength results by epoch and lunation</h2><p>Central anchors: B 4700 Å, R 6600 Å, Z 8500 Å.
Within-night columns need ≥3 exposures per night; between-night columns need ≥2 such nights.
NaN means insufficient repetition. Single-exposure periods are listed in the main report.</p>
<div class="scroll">{table(epoch)}</div>
<h2>Science and sky-fiber summaries</h2><p>Science, SkyE and SkyW samples are summarized separately; the same verified 1% selection is used throughout.</p>
<div class="scroll">{table(telescope_summary)}</div><h2>Large-departure review inventory</h2><p>These counts use each fiber’s P90 night departure &gt;10%, a descriptive inspection flag. All flagged measurements remain in the stability summaries.</p>
<div class="scroll">{table(tail_counts)}</div><h2>Dependence on wavelength</h2><p>These summaries combine per-fiber, per-lunation statistics across sampled epochs.
B, R and Z are separate spectral channels. W90 tests LSF wings as well as the core width.</p><div class="scroll">{table(spectral)}</div>
<h2>All-fiber stability plots</h2><p>The <a href="../lsf_epoch_coverage.html">main report</a> shows coverage and FWHM, W50, W90 and centroid scatter against physical fiber ID,
with separate Science, Sky East and Sky West columns for each epoch/lunation. Numerical tables retain exact per-fiber values.</p>
<h2>Individual-fiber time series</h2><p>Examples are drawn from the two lowest and two highest median B-scatter epoch–lunation groups.
Within each group, the Science fiber nearest the median Science-fiber B scatter is shown. Only one physical fiber appears in each figure,
with B/R/Z columns and FWHM, W50, W90 and centroid rows at central wavelengths. Width changes are percentages relative to the plotted median; centroid changes subtract that median, in Å.
Points show every accepted observation; lines connect observations within a night. Axes are symmetric-log, linear through ±1% for widths and ±0.01 Å for the centroid; all outliers are retained.
The first four examples illustrate typical Science fibers in selected groups; a fifth shows SkyW fiber 181 from the large-departure group.
Exposure/night counts in titles describe B coverage; the tables retain channel-specific counts. Sparse single-night examples illustrate exposure repeatability.</p>
{''.join(example_html)}
<h2>Reproducible tables</h2><ul>
{''.join(f'<li><a href="{p}">{p}</a></li>' for p in ['channel_summary.csv','telescope_channel_summary.csv','review_flag_counts.csv','central_night_departure_review_gt10pct.csv','review_flagged_night_diagnostics.csv','epoch_lunation_stability.csv','wavelength_stability.csv','fiber_full_night_changes.csv',
'fiber_lunation_stability.csv','fiber_night_repeatability.csv','fiber_between_night_stability.csv','consecutive_night_changes.csv',
'accepted_lsf_diagnostics.parquet','diagnostic_equivalence_checks.csv','example_fibers.csv','diagnostic_provenance.json','summary.json'])}</ul>
<p>Source fits are read only. Analysis script: analysis/lsf_epoch_coverage/analyze_stability.py.</p></html>'''
    (DEST/'lsf_epoch_lunation_stability.html').write_text(html)
    print(json.dumps(payload,indent=2),flush=True)
    print(spectral.to_string(index=False),flush=True)


if __name__=='__main__':main()
