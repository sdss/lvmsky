"""Publish a chronological, per-lunation audit with actual individual-fiber LSF profiles."""
from pathlib import Path
import json
import re
import html
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from analysis.lsf_epoch_coverage.build_coverage import OUT,fit_files,write_time_overview
from analysis.lsf_epoch_coverage.measurement_view import publish_measurements
from skysub.sky_decomp.lsf_spline2d import mspline_basis

CENTRAL={'B':4700.,'R':6600.,'Z':8500.}
TELESCOPES={'Sci':'Science','SkyE':'Sky East','SkyW':'Sky West'}
CHANNEL_COLORS={'B':'#0072B2','R':'#D55E00','Z':'#009E73'}
SCATTER_METRICS=[('fwhm_angstrom_rmad_pct','FWHM scatter (%)',1.),
                 ('w50_angstrom_rmad_pct','W50 scatter (%)',1.),
                 ('w90_angstrom_rmad_pct','W90 scatter (%)',1.),
                 ('centroid_angstrom_rmad','Centroid scatter (Å)',.01)]
ABSOLUTE_SCATTER_METRICS=[('fwhm_angstrom_rmad','FWHM scatter (Å)',.01),
                        ('w50_angstrom_rmad','W50 scatter (Å)',.01),
                        ('w90_angstrom_rmad','W90 scatter (Å)',.01),
                        ('centroid_angstrom_rmad','Centroid scatter (Å)',.01)]


def plot_fiber_metrics(rows,repeat,selected_exposures,label,epoch,absolute=False,levels=None):
    """Show all fibers; a symmetric-log scale retains zero and every outlier."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    assert not repeat.duplicated(['fiberid','channel']).any()
    plotted_metrics=ABSOLUTE_SCATTER_METRICS if absolute else SCATTER_METRICS
    spectra=rows.groupby(['telescope','fiberid','expnum']).agg(accepted=('accepted','max'),gaia_verified=('gaia_verified','max')).reset_index()
    coverage=spectra.groupby(['telescope','fiberid']).agg(
        accepted_exposures=('accepted','sum'),verified_exposures=('gaia_verified','sum')).reset_index()
    identity=rows[['fiberid','spectrographid']].drop_duplicates()
    assert identity.fiberid.is_unique
    ranges=identity.groupby('spectrographid').fiberid.agg(['min','max']).sort_index()
    assert (ranges['max'].to_numpy()[:-1]<ranges['min'].to_numpy()[1:]).all()
    boundaries=(ranges['max'].to_numpy()[:-1]+ranges['min'].to_numpy()[1:])/2
    scatter_rows=[2,4,6,8] if levels is not None else [1,2,3,4]
    value_rows=[1,3,5,7]
    fig,axes=plt.subplots(9 if levels is not None else 5,3,figsize=(14,21) if levels is not None else (14,12),
                         sharex=True,sharey='row',gridspec_kw={'width_ratios':[2,1,1]})
    level_axes={}
    for col,(telescope,title) in enumerate(TELESCOPES.items()):
        group=coverage[coverage.telescope.eq(telescope)]
        axes[0,col].set_title(f'{title} ({telescope}) • {group.fiberid.nunique()} fibers',pad=18)
        verified=group[group.verified_exposures>0]
        unknown=group[group.verified_exposures.eq(0)]
        axes[0,col].scatter(verified.fiberid,verified.accepted_exposures,s=7,color='#222',alpha=.65)
        axes[0,col].scatter(unknown.fiberid,np.zeros(len(unknown)),s=14,color='#777',marker='x')
        for sp,limits in ranges.iterrows():
            axes[0,col].text((limits['min']+limits['max'])/2,1.015,f'SP{sp}',ha='center',fontsize=8,
                             transform=axes[0,col].get_xaxis_transform())
        for channel,color in CHANNEL_COLORS.items():
            points=repeat[repeat.telescope.eq(telescope)&repeat.channel.eq(channel)]
            for row,(metric,_,_) in enumerate(plotted_metrics,1):
                values=points[metric]
                assert np.isfinite(values).all() and values.ge(0).all()
                axes[scatter_rows[row-1],col].scatter(points.fiberid,values,s=7,color=color,alpha=.65)
        if levels is not None:
            group_levels=levels[levels.telescope.eq(telescope)]
            for row,(scatter_metric,_,_) in enumerate(plotted_metrics,1):
                metric=scatter_metric.split('_rmad')[0]
                value_axis=axes[value_rows[row-1],col];level_axes[row,col]=value_axis
                for channel,color in CHANNEL_COLORS.items():
                    channel_levels=group_levels[group_levels.channel.eq(channel)]
                    for _,segment in channel_levels.groupby('spectrographid'):
                        segment=segment.sort_values('fiberid')
                        for statistic,style in [('median','-'),('mean','--')]:
                            values=segment[f'{metric}_{statistic}']
                            assert np.isfinite(values).all()
                            value_axis.plot(segment.fiberid,values,color=color,linestyle=style,linewidth=.8,
                                alpha=.85,marker='_' if len(segment)==1 else None,markersize=6)
                if col==0:value_axis.set_ylabel(f'{["FWHM","W50","W90","Centroid offset"][row-1]} value (Å)',fontsize=8)
        if not repeat.telescope.eq(telescope).any():
            for ax in axes[scatter_rows,col]:
                ax.text(.02,.97,'Scatter needs ≥3 exposures',ha='left',va='top',
                        fontsize=7,transform=ax.transAxes,bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1))
        axes[-1,col].set_xlabel('Physical fiber ID')
    axes[0,0].set_ylabel('Accepted exposures\n(≥1 channel)',fontsize=8)
    axes[0,0].set_ylim(-.5,max(3,selected_exposures)*1.06)
    for row,(metric,title,threshold) in enumerate(plotted_metrics,1):
        upper=max(threshold,float(repeat[metric].max()) if len(repeat) else threshold)*1.12
        ticks=np.array([0,.5,1,2,5,10,20,50,100]) if threshold==1 else np.array([0,.005,.01,.02,.05,.1,.2,.5,1,2,5])
        panel_row=scatter_rows[row-1]
        axes[panel_row,0].set_ylabel(title,fontsize=8)
        for ax in axes[panel_row]:
            ax.set_yscale('symlog',linthresh=threshold,linscale=1)
            ax.set_ylim(0,upper)
            ax.set_yticks(ticks[ticks<=upper])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'))
        if levels is not None and len(levels):
            prefix=metric.split('_rmad')[0]
            values=levels[[f'{prefix}_mean',f'{prefix}_median']].to_numpy()
            low,high=float(values.min()),float(values.max())
            pad=max((high-low)*.08,abs(high)*.01,1e-4)
            for col in range(3):level_axes[row,col].set_ylim(low-pad,high+pad)
    for ax in axes.flat:
        ax.set_xlim(0,1950);ax.set_xticks([1,500,1000,1500,1944])
        for boundary in boundaries:ax.axvline(boundary,color='#777',linestyle='--',linewidth=.7)
        ax.grid(axis='y',alpha=.22);ax.tick_params(axis='x',labelsize=9);ax.tick_params(axis='y',labelsize=8)
    for ax in axes[:,1:].flat:ax.set_xticks([1,648,1296,1944])
    handles=[Line2D([],[],marker='o',linestyle='',color=color,label=f'{channel}: {CENTRAL[channel]:.0f} Å')
             for channel,color in CHANNEL_COLORS.items()]
    handles.append(Line2D([],[],marker='o',linestyle='',color='#222',label='Accepted exposure counts'))
    handles.append(Line2D([],[],marker='x',linestyle='',color='#777',label='Gaia unverified (coverage row only)'))
    if levels is not None:
        handles.extend([Line2D([],[],color='#555',linestyle='-',linewidth=.8,label='Median value (Å)'),
                        Line2D([],[],color='#555',linestyle='--',linewidth=.8,label='Mean value (Å)'),
                        Line2D([],[],color='#555',marker='_',linestyle='',markersize=6,label='Isolated fiber value')])
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.965 if levels is not None else .952),ncol=3,fontsize=8)
    units='absolute scatter in Å' if absolute else 'relative width scatter; centroid in Å'
    if levels is not None:units+='; mean and median levels in Å'
    fig.suptitle(f'Full Moon {label} • calibration epoch {epoch}: {units}',fontsize=13)
    scale='0.01 Å for all metrics' if absolute else '1% (widths) / 0.01 Å (centroid)'
    fig.text(.5,.012,f'Scatter axes: linear through {scale}, logarithmic above. All outliers retained. '
             'No scatter point = fewer than three accepted exposures.',ha='center',fontsize=9)
    if levels is not None:
        fig.text(.5,.035,'Each metric has two separate rows: median/mean values in Å (linear), then scatter points. Values need ≥1 exposure; scatter needs ≥3.',ha='center',fontsize=8)
    if not absolute:
        fig.text(.5,.023,'Width % = 100 × (1.4826 MAD) / median of the same metric for this fiber, lunation and epoch.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.05,1,.93] if levels is not None else [0,.05,1,.92])
    return fig


def profile_summary(rows):
    cols=[f'lsf_mass_{i:02d}' for i in range(11)]
    mass=rows[cols].median().to_numpy(dtype=float)
    assert np.isfinite(mass).all() and mass.min()>=-1e-7 and mass.sum()>0
    half=float(rows.offset_half_width_angstrom.iloc[0])
    assert rows.offset_half_width_angstrom.eq(half).all()
    delta=np.linspace(-half,half,1201)
    return delta,mspline_basis(delta,n_basis=11,half_width=half)@(mass/mass.sum())


def show_table(df,labels=None):
    return df.rename(columns=labels or {}).to_html(index=False,na_rep='—',float_format=lambda v:f'{v:.2f}',
        formatters={column:lambda v:f'{v:.4f}' for column in ['FWHM (Å)','W50 (Å)','W90 (Å)','Centroid (Å)']})


def main():
    dest=OUT/'lunations';dest.mkdir(exist_ok=True)
    summary=json.loads((OUT/'summary.json').read_text())
    counts=json.loads((OUT/'exposure_selection_counts.json').read_text())
    exposures=pd.read_csv(OUT/'exposure_epoch_mapping.csv').sort_values('obstime')
    sample=pd.read_parquet(OUT/'fiber_exposure_channel.parquet',columns=['expnum','fiberid','channel','telescope',
        'spectrographid','calibration_epoch','full_moon_utc','mjd','accepted','gaia_verified','gaia_clean_1pct','fit_status'])
    diag_cols=['expnum','fiberid','channel','telescope','spectrographid','calibration_epoch','full_moon_utc',
        'mjd','observation_mjd_utc','wavelength_angstrom','fwhm_angstrom','w50_angstrom','w90_angstrom','centroid_angstrom','chi2_red','diagnostic_valid',
        'offset_half_width_angstrom']+[f'lsf_mass_{i:02d}' for i in range(11)]
    diagnostics=pd.read_parquet(OUT/'stability/accepted_lsf_diagnostics.parquet',columns=diag_cols)
    diagnostics=diagnostics[diagnostics.diagnostic_valid]
    central=diagnostics[np.array([w==CENTRAL[c] for w,c in zip(diagnostics.wavelength_angstrom,diagnostics.channel)])].copy()
    level_keys=['calibration_epoch','full_moon_utc','fiberid','telescope','spectrographid','channel']
    value_metrics=['fwhm_angstrom','w50_angstrom','w90_angstrom','centroid_angstrom']
    metric_levels=central.groupby(level_keys)[value_metrics].agg(['mean','median'])
    metric_levels.columns=[f'{metric}_{statistic}' for metric,statistic in metric_levels.columns]
    metric_levels['n_exposures']=central.groupby(level_keys).size()
    metric_levels=metric_levels.reset_index()
    metric_levels.to_csv(OUT/'fiber_lunation_metric_levels.csv',index=False)
    measurement_section=publish_measurements(central,exposures,OUT)
    metrics=pd.read_csv(OUT/'stability/epoch_lunation_stability.csv')
    within=pd.read_csv(OUT/'stability/fiber_night_repeatability.csv',usecols=['calibration_epoch','full_moon_utc',
        'fiberid','channel','mjd','wavelength_angstrom','n_exposures'])
    within=within[np.array([w==CENTRAL[c] for w,c in zip(within.wavelength_angstrom,within.channel)])]
    changes=pd.read_csv(OUT/'stability/fiber_full_night_changes.csv')
    changes=changes[np.array([w==CENTRAL[c] for w,c in zip(changes.wavelength_angstrom,changes.channel)])]
    per_fiber=pd.read_csv(OUT/'stability/fiber_lunation_stability.csv')
    per_fiber=per_fiber[np.array([w==CENTRAL[c] for w,c in zip(per_fiber.wavelength_angstrom,per_fiber.channel)])]
    pair_changes=pd.read_csv(OUT/'stability/consecutive_night_changes.csv',usecols=['calibration_epoch','full_moon_utc',
        'fiberid','channel','wavelength_angstrom','mjd','previous_night_mjd','fwhm_angstrom_night_change_pct'])
    pair_changes=pair_changes[np.array([w==CENTRAL[c] for w,c in zip(pair_changes.wavelength_angstrom,pair_changes.channel)])]
    pair_changes['absolute_change_pct']=pair_changes.fwhm_angstrom_night_change_pct.abs()
    from astropy.io import fits
    from astropy.table import Table
    from analysis.lsf_epoch_coverage.build_coverage import RUN
    indexes=[]
    for path in fit_files():
        with fits.open(path,memmap=True) as h:
            t=Table(h['INPUT_META'].data)[['expnum','fiberid']].to_pandas()
            t['fit_file']=str(path);t['fit_row']=np.arange(len(t));indexes.append(t)
    fit_index=pd.concat(indexes,ignore_index=True)
    inputs=pd.read_parquet(OUT/'photometry.parquet')
    working_fibers=inputs.loc[inputs.fibstatus.eq(0),'fiberid'].nunique()
    inputs=inputs[inputs.gaia_verified & inputs.gaia_clean_1pct & inputs.fibstatus.eq(0)]
    inputs=inputs.merge(exposures[['expnum','mjd','obstime','full_moon_utc','nominal_epoch','source_path','has_completed_fit']],on='expnum',validate='many_to_one')
    assert inputs.gaia_ratio.lt(.01).all()
    inputs.to_parquet(OUT/'gaia_clean_fiber_sample.parquet',index=False)
    clean=pd.read_parquet(OUT/'fiber_exposure_channel.parquet')
    clean=clean[clean.accepted].merge(fit_index,on=['expnum','fiberid'],validate='many_to_one')
    assert clean.gaia_ratio.lt(.01).all() and clean.gaia_verified.all() and clean.quality_ok.all()
    assert len(clean)==summary['accepted_channel_rows']
    good_fiber_count=clean.fiberid.nunique()
    accepted_spectrum_channels=clean.groupby(['expnum','fiberid']).size()
    accepted_spectra=len(accepted_spectrum_channels)
    complete_spectra=int(accepted_spectrum_channels.eq(3).sum())
    clean.to_parquet(OUT/'clean_lsf_sample.parquet',index=False)
    clean.groupby(['full_moon_utc','calibration_epoch','fiberid','telescope','channel']).agg(
        accepted_exposures=('expnum','size'),observing_nights=('mjd','nunique'),
        first_timestamp_mjd_utc=('observation_mjd_utc','min'),last_timestamp_mjd_utc=('observation_mjd_utc','max')).reset_index().to_csv(
            OUT/'clean_lsf_fiber_inventory.csv',index=False)
    del clean
    inventory=[];sections=[];examples=[]
    for i,((moon,epoch),exp) in enumerate(exposures.groupby(['full_moon_utc','nominal_epoch'],sort=True)):
        identifier=f'lunation-{i+1:02d}';label=moon[:10]
        rows=sample[sample.full_moon_utc.eq(moon)&sample.calibration_epoch.eq(epoch)]
        good=rows[rows.accepted]
        data=central[central.full_moon_utc.eq(moon)&central.calibration_epoch.eq(epoch)]
        repeat=per_fiber[per_fiber.full_moon_utc.eq(moon)&per_fiber.calibration_epoch.eq(epoch)&(per_fiber.n_exposures>=3)]
        levels=metric_levels[metric_levels.full_moon_utc.eq(moon)&metric_levels.calibration_epoch.eq(epoch)]
        multi=changes[changes.full_moon_utc.eq(moon)&changes.calibration_epoch.eq(epoch)]
        available=exp.has_completed_fit if 'has_completed_fit' in exp else exp.expnum.isin(good.expnum)
        info=dict(full_moon_utc=moon,epoch_ids=','.join(map(str,sorted(exp.nominal_epoch.unique()))),
            first_exposure_utc=exp.obstime.min(),last_exposure_utc=exp.obstime.max(),
            selected_exposures=len(exp),selected_nights=exp.mjd.nunique(),completed_fit_exposures=int(available.sum()),
            gaia_verified_exposures=rows[rows.gaia_verified].expnum.nunique(),
            accepted_exposures=good.expnum.nunique(),accepted_physical_fibers=good.fiberid.nunique(),
            repeated_fiber_channels=len(repeat),cross_night_fiber_channels=len(multi))
        inventory.append(info)
        section=f'<section id="{identifier}"><h2>{label}: observed {info["first_exposure_utc"]} to {info["last_exposure_utc"]} UTC</h2>'
        section+=f'<p>Calibration epoch: <b>{info["epoch_ids"]}</b>. {len(exp)} Moon/Sun-selected exposures (before Gaia) on {info["selected_nights"]} nights; '
        section+=f'{info["completed_fit_exposures"]} have individual fits; {info["gaia_verified_exposures"]} have verified Gaia measurements. '
        section+=f'<b>{info["accepted_physical_fibers"]} physical fibers</b> have at least one accepted LSF channel measurement.</p>'
        table=[];scatter=[]
        for channel in ['B','R','Z']:
            accepted=good[good.channel.eq(channel)]
            rs=repeat[repeat.channel.eq(channel)];ms=multi[multi.channel.eq(channel)]
            table.append(dict(Channel=channel,**{'Wavelength (Å)':CENTRAL[channel],'Accepted exposures':accepted.expnum.nunique(),
                'Accepted fibers':accepted.fiberid.nunique(),'Fibers with ≥3 exposures':len(rs),
                'Fibers with ≥2 qualifying nights':len(ms),'Median within-lunation scatter (%)':rs.fwhm_angstrom_rmad_pct.median(),
                'P90 within-lunation scatter (%)':rs.fwhm_angstrom_rmad_pct.quantile(.9),
                'Median sampled-night change (%)':ms.median_abs_full_night_change_pct.median(),
                'P90 sampled-night change (%)':ms.median_abs_full_night_change_pct.quantile(.9)}))
            scatter.append(dict(Channel=channel,**{'Fibers with ≥3 exposures':len(rs),
                'FWHM (%)':rs.fwhm_angstrom_rmad_pct.median(),
                'FWHM (Å)':rs.fwhm_angstrom_rmad.median(),
                'W50 (%)':rs.w50_angstrom_rmad_pct.median(),
                'W50 (Å)':rs.w50_angstrom_rmad.median(),
                'W90 (%)':rs.w90_angstrom_rmad_pct.median(),
                'W90 (Å)':rs.w90_angstrom_rmad.median(),
                'Centroid (Å)':rs.centroid_angstrom_rmad.median()}))
        if any(row['Fibers with ≥2 qualifying nights'] for row in table):
            section+='<p>Repeated measurements support within-lunation and cross-night comparisons.</p>'
        elif any(row['Fibers with ≥3 exposures'] for row in table):
            section+='<p>Repeated measurements support within-lunation scatter. Cross-night comparison needs two nights with ≥3 accepted exposures each.</p>'
        elif len(good):
            section+='<p>Individual LSF measurements are available; fewer than three per fiber leaves repeatability unestimated.</p>'
        else:
            section+='<p>Selected exposures are listed; accepted LSF sampling is shown below.</p>'
        section+='<h3>Coverage and sampled-night FWHM change</h3><div class="scroll">'+show_table(pd.DataFrame(table).drop(columns=['Median within-lunation scatter (%)','P90 within-lunation scatter (%)']))+'</div>'
        section+='<h3>Within-lunation scatter: percent and Å</h3><p>Median of per-fiber robust scatter at the wavelengths above. '
        section+='Absolute scatter: 1.4826 × MAD in Å; relative width scatter: 100 × absolute scatter / median width. '
        section+='Each is calculated per fiber before summarizing. These are spreads, not the widths themselves. — means fewer than three accepted exposures.</p><div class="scroll">'+show_table(pd.DataFrame(scatter))+'</div>'
        fig=plot_fiber_metrics(rows,repeat,len(exp),label,epoch,levels=levels)
        name=f'{identifier}_coverage.png';fig.savefig(dest/name,dpi=150);plt.close(fig)
        section+='<h3>Coverage and stability by physical fiber</h3><p>Columns: Science, Sky East (SkyE), Sky West (SkyW). '
        section+='X: physical fiber ID. Top: exposures with ≥1 accepted LSF channel (black); each following metric has a value row, then a separate scatter row. '
        section+='B/R/Z are colored points. All outliers are retained; scatter axes are linear near zero and logarithmic above the stated threshold.</p>'
        section+='<p><b>Separate panels:</b> each metric has a value row (Å, linear axis) followed by a scatter row. Thin solid lines show each fiber’s median and dashed lines its arithmetic mean. '
        section+='They use all accepted exposures of that fiber/channel in this period (≥1); scatter points require ≥3. '
        section+='Lines join measured fiber values within each spectrograph, without smoothing or crossing spectrograph boundaries.</p>'
        section+='<p>For one accepted exposure, mean and median coincide. An isolated fiber value is shown as a short horizontal tick in the value panel.</p>'
        section+='<p><b>Percent denominator:</b> the median of the same width metric across accepted exposures of this physical fiber in this lunation and calibration epoch. '
        section+='Science occupies half the panel width; each sky telescope occupies one quarter. Dashed lines separate spectrographs SP1/SP2/SP3.</p>'
        section+=f'<img src="lunations/{name}" alt="{label}: coverage and four stability metrics against fiber ID, separated by telescope">'
        fig=plot_fiber_metrics(rows,repeat,len(exp),label,epoch,absolute=True,levels=levels)
        name=f'{identifier}_absolute_scatter.png';fig.savefig(dest/name,dpi=150);plt.close(fig)
        section+='<h3>Absolute scatter in Å</h3><p>The same fibers and accepted exposures, with FWHM, W50, W90 and centroid spread in Å. '
        section+='Axes are linear through 0.01 Å and logarithmic above; all outliers are retained.</p>'
        section+='<p>Value panels show solid median and dashed mean lines in Å on linear axes. Separate scatter panels show points for the parameter’s spread; each panel has one Y-axis.</p>'
        section+=f'<img src="lunations/{name}" alt="{label}: absolute FWHM, W50, W90 and centroid scatter in Angstrom against fiber ID">'
        if not data.empty:
            science=data[data.telescope.eq('Sci')]
            candidates=science.groupby('fiberid').agg(channels=('channel','nunique'),nights=('mjd','nunique'),exposures=('expnum','nunique'),sp=('spectrographid','first'))
            candidates=candidates[candidates.channels.eq(3)].copy()
            # Coverage, not measured scatter, determines the illustrative fiber.
            if len(candidates):
                candidates['center_distance']=abs(candidates.index.to_numpy()-np.median(candidates.index))
                chosen=candidates.sort_values(['nights','exposures','center_distance'],ascending=[False,False,True]).index[0]
                example=data[data.fiberid.eq(chosen)]
                fig,axes=plt.subplots(2,3,figsize=(14,8))
                nightly_records=[]
                for j,channel in enumerate(['B','R','Z']):
                    subset=example[example.channel.eq(channel)].sort_values('observation_mjd_utc')
                    axes[1,j].plot(pd.to_datetime(subset.observation_mjd_utc,unit='D',origin='1858-11-17'),subset.fwhm_angstrom,
                                   '.',color='#8b8b8b',markersize=4)
                    for night,measurements in subset.groupby('mjd'):
                        if len(measurements)<3:continue
                        delta,density=profile_summary(measurements)
                        date=pd.to_datetime(float(night),unit='D',origin='1858-11-17').strftime('%m-%d')
                        axes[0,j].plot(delta,density,lw=1.5,label=f'{date}, N={len(measurements)}')
                        axes[1,j].plot(pd.to_datetime(measurements.observation_mjd_utc.median(),unit='D',origin='1858-11-17'),
                                       measurements.fwhm_angstrom.median(),'s',color='black',markersize=6)
                        nightly_records.append(dict(full_moon_utc=moon,fiberid=int(chosen),channel=channel,mjd=int(night),
                            n_exposures=len(measurements),median_fwhm_angstrom=measurements.fwhm_angstrom.median(),
                            median_chi2=measurements.chi2_red.median()))
                    axes[0,j].set_title(f'{channel} at {CENTRAL[channel]:.0f} Å');axes[0,j].set_xlabel('Wavelength offset (Å)')
                    axes[0,j].set_ylabel('Normalized LSF density (Å⁻¹)')
                    if axes[0,j].lines:axes[0,j].legend(fontsize=8)
                    else:
                        axes[0,j].set_axis_off()
                        axes[0,j].text(.5,.5,f'{channel} at {CENTRAL[channel]:.0f} Å\nNo night with ≥3 accepted exposures',ha='center',va='center',transform=axes[0,j].transAxes,fontsize=10)
                    axes[1,j].set_ylabel('FWHM (Å)');axes[1,j].set_xlabel('Observation date (UTC)');axes[1,j].tick_params(axis='x',rotation=45)
                    import matplotlib.dates as mdates
                    axes[1,j].xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3,maxticks=5))
                    axes[1,j].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d\n%H:%M'))
                    if subset.observation_mjd_utc.nunique()==1:
                        point=pd.to_datetime(subset.observation_mjd_utc.iloc[0],unit='D',origin='1858-11-17')
                        axes[1,j].set_xlim(point-pd.Timedelta(minutes=30),point+pd.Timedelta(minutes=30))
                fig.suptitle(f'Physical Science fiber {int(chosen)} • full Moon {label}: nightly median reconstructed profiles and individual widths')
                fig.tight_layout();name=f'{identifier}_profiles.png';fig.savefig(dest/name,dpi=150);plt.close(fig)
                section+=f'<p><a href="lunations/{name}">Supplementary profile illustration: Science fiber {int(chosen)}</a>. '
                section+='Chosen by coverage; normalized median M-spline profiles need ≥3 exposures per night. Median individual FWHM can differ from reconstructed profile FWHM.</p>'
                examples.extend(nightly_records)
        else:
            section+='<p>These exposure times are represented in the selection audit. Individual-fiber stability is unmeasured until verified Gaia eligibility and completed LSF fits are available.</p>'
        outliers=pair_changes[pair_changes.full_moon_utc.eq(moon)&pair_changes.calibration_epoch.eq(epoch)].sort_values('absolute_change_pct',ascending=False).head(5)
        if len(outliers):
            section+='<h3>Largest measured changes between sampled nights</h3><p>These rows remain in the audit. Large changes require fit inspection before constructing a calibration reference.</p>'
            section+='<div class="scroll">'+show_table(outliers[['fiberid','channel','previous_night_mjd','mjd','absolute_change_pct']],
                {'fiberid':'Physical fiber','channel':'Channel','previous_night_mjd':'Earlier night MJD','mjd':'Later night MJD','absolute_change_pct':'Absolute FWHM change (%)'})+'</div>'
        section+='</section>';sections.append(section)
        print(f'Published section {i+1}: {moon}',flush=True)
    overview=pd.DataFrame(inventory);overview.to_csv(OUT/'lunation_inventory.csv',index=False)
    pd.DataFrame(examples).to_csv(OUT/'lunation_profile_examples.csv',index=False)
    # The sample explicitly covers every selected exposure, fiber and channel opportunity.
    assert overview.selected_exposures.sum()==len(exposures)
    assert not sample.duplicated(['expnum','fiberid','channel']).any()
    chronology=write_time_overview(exposures)
    navigation=' · '.join(f'<a href="#lunation-{i+1:02d}">{r.full_moon_utc[:10]} / {r.epoch_ids}</a>' for i,r in enumerate(overview.itertuples()))
    reference_link = ('<p><a href="lsf_calibration/lsf_reference_surfaces.html">Dispersion-weighted LSF reference FITS, statistics and methods</a> (verified Gaia flux ratio &lt;0.1; separate epoch and lunation products).</p>'
                      if (OUT/'lsf_calibration/lsf_reference_surfaces.html').exists() else '')
    content=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Per-fiber LSF samples and stability by lunation</title>
<style>body{{font:16px/1.5 sans-serif;max-width:1450px;margin:30px auto;padding:0 20px;color:#17202a}}h2{{margin-top:45px}}table{{border-collapse:collapse;font-size:13px}}td,th{{padding:6px;border:1px solid #ccc}}.scroll{{overflow:auto}}img{{width:100%}}.note{{padding:16px;background:#edf4fa;border-left:5px solid #145da0}}section{{border-top:2px solid #ddd}}</style>
<h1>Per-fiber line-spread function (LSF): coverage and stability</h1>
{reference_link}
<div class="note"><p><b>{counts['geometry_selected']} Moon/Sun-selected exposures · {summary['full_moons']} sampled lunations · {summary['actual_epochs']} calibration epochs.</b>
Selection: Moon altitude &gt;60°, illumination &gt;95%, Sun altitude &lt;−18°. Gaia tables and individual fits cover all {summary['completed_fit_exposures']} selected exposures.</p>
<p><b>{summary['accepted_channel_rows']:,} accepted LSF measurements.</b> One measurement is one exposure × physical fiber × channel (B, R or Z), after Gaia contamination &lt;1% and fit-quality checks. One fiber in one exposure can contribute three measurements.</p>
<p><b>{accepted_spectra:,} distinct fiber spectra</b> have an accepted LSF; {complete_spectra:,} pass in all three channels.
Each spectrum has one common iterative fit, with separate B/R/Z LSF substeps. Gaia, overall fit status and convergence are checked at spectrum level; channel-specific LSF diagnostics are then checked separately.</p>
<p><b>{good_fiber_count:,} physical fibers have an accepted measurement in at least one lunation.</b> This total spans the complete data set; each period has its own coverage below.</p></div>
<p>Gaia source: only the <a href="http://lumen.abudhabi.nyu.edu/sas/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/">requested Lumen data set</a>, <code>{summary['gaia_source_directory']}</code>.
Contamination = summed aperture-source Gaia G-band flux / synthetic G-band flux from observed FLUX+SKY. The ratio must be &lt;0.01; unknown contamination remains unknown.
See the <a href="gaia_photometry_provenance.json">source manifest</a>.</p>
<h2>Dates and coverage</h2><p>A <b>calibration epoch</b> is the validity interval of the calibration recorded in each reduced exposure (SFrame); its ID is the epoch start MJD.
A <b>lunation</b> is labelled by the nearest full-Moon UTC date. These observations cover its bright-Moon window.
The timeline gives exact exposure UTC dates and calibration boundaries. MJD means Modified Julian Date; integer observing-night MJD groups nights, while continuous timestamp MJD positions exposures.</p>
<p>The default plots place <b>physical fiber ID on X</b>, with separate Science, Sky East and Sky West columns.
Science has half the plot width; Sky East and Sky West each have one quarter. Dashed vertical lines mark the verified spectrograph boundaries (between fiber IDs 648/649 and 1296/1297).
They show exposure counts with ≥1 accepted LSF channel (black) and all four scatter metrics (B/R/Z colors); zero scatter is shown, while insufficient repetition has no scatter point.
Coverage-row gray crosses identify unverified Gaia eligibility; black zero counts are verified opportunities with no accepted fit.
Each metric shares a Y scale across telescopes within a figure; scales vary between lunations. The <a href="fiber_lunation_coverage.csv">coverage table</a> also gives hours, time span, largest gap and concentration among nights.
Effective nights = (Σ t<sub>night</sub>)² / Σ t<sub>night</sub>², where t<sub>night</sub> is total accepted exposure time that night. Night counts describe sampling, not independence.</p>
<h2>What stability means</h2><ul>
<li><b>FWHM:</b> full width at half maximum.</li>
<li><b>W50 / W90:</b> widths enclosing the central 50% / 90% of LSF area.</li>
<li><b>Centroid:</b> mean wavelength offset of the LSF, in Å.</li></ul>
<p><b>Within-lunation scatter</b> is computed separately for each fiber, calibration epoch, lunation, channel and wavelength, using ≥3 accepted exposures.
Absolute scatter for FWHM, W50, W90 and centroid: 1.4826 × MAD (Å). Relative width scatter: 100 × absolute scatter / median width (%). MAD is median absolute deviation.
<b>The denominator is the median of that same metric for that physical fiber, channel, wavelength, lunation and calibration epoch, using the plotted accepted exposures.</b>
The tables below give the median of these per-fiber scatters for B/R/Z at 4700/6600/8500 Å. Detailed tables retain five wavelengths per channel.</p>
<p><b>Sampled-night FWHM change:</b> absolute percentage change between adjacent observed night medians, relative to the earlier night; each night needs ≥3 exposures.
Tables summarize each fiber's median change with a median and P90 (90th percentile) across fibers. Night gaps can exceed one day.
— means insufficient sampling. One measurement gives an LSF, but repetition is needed to measure stability. These statistics include fitting noise and observing conditions.</p>
<h2>Using scatter to assess LSF transfer</h2><p>Absolute scatter quantifies recovered-LSF variation when reusing a reference for the <b>same physical fiber, calibration epoch and wavelength</b> under the sampled bright-Moon conditions.
Smaller FWHM, W50, W90 and centroid spreads support more repeatable parameters. These empirical exposure-to-exposure spreads include fit noise.
Formal fit errors, uncertainty on a median and transfer prediction errors require separate validation.</p>
<p>Direct transfer validation requires a reference built without the target exposure, followed by comparison on that exposure; withholding a whole night tests transfer between nights.
Full-profile residuals are needed alongside width and centroid metrics to test shape. The current report measures repeatability; it does not yet perform this held-out reference test.</p>
{measurement_section}
<h2>All sampled lunations and calibration epochs</h2><p>{navigation}</p><div class="scroll">{show_table(overview,{'full_moon_utc':'Nearest full Moon (UTC; lunation label)','epoch_ids':'Calibration epoch start MJD','first_exposure_utc':'First selected exposure UTC','last_exposure_utc':'Last selected exposure UTC','selected_exposures':'Moon/Sun-selected exposures (before Gaia)','selected_nights':'Observed nights (before Gaia)','completed_fit_exposures':'Exposures with completed fits','gaia_verified_exposures':'Exposures with ≥1 Gaia-verified fiber','accepted_exposures':'Exposures with ≥1 accepted fiber-channel','accepted_physical_fibers':'Fibers with ≥1 accepted LSF','repeated_fiber_channels':'Fiber-channel samples with ≥3 exposures','cross_night_fiber_channels':'Fiber-channel samples with ≥2 qualifying nights'})}</div>
{chronology}
{''.join(sections)}
<h2>Data and quality checks</h2><p>Accepted measurements require a working fiber, a Solved fit, channel convergence, positive valid-pixel count and widths,
finite nonnegative chi-square, input-valid fraction ≥0.8 and no channel failure flag.
Gaia verification requires matching aperture coordinates, complete source photometry, positive observed G and ≥95% valid passband weight on the available wavelength grid.
Observed G and errors are recomputed from DRP 1.3.2 FLUX+SKY. Verified apertures with zero catalog sources have zero Gaia flux.</p>
<p><a href="fiber_lunation_metric_levels.csv">Mean and median parameter values for each fiber, channel, lunation and epoch</a>: computed directly from the accepted measurements, in Å, with supporting exposure counts.</p>
<ul>{''.join(f'<li><a href="{p}">{p}</a></li>' for p in ['gaia_clean_fiber_sample.parquet','clean_lsf_sample.parquet','clean_lsf_fiber_inventory.csv','exposure_selection_audit.csv','exposure_epoch_mapping.csv','lunation_inventory.csv','fiber_exposure_channel.parquet','photometry.parquet','fiber_lunation_coverage.csv','gaia_photometry_provenance.json','lunation_profile_examples.csv','stability/fiber_lunation_stability.csv','stability/fiber_night_repeatability.csv','stability/fiber_full_night_changes.csv','stability/lsf_stability_diagnostics.html'])}</ul>
<p>Diagnostics use saved LSF coefficients, checked against the production evaluator. G-band S/N and fit chi-square accompany the tables.
Inspect large departures and validate uncertainties before constructing calibration-reference profiles.</p>
<h2>Code and reproducibility</h2>
<p>Scripts: <code>/home/ik52/sdss5/lvmsky/analysis/lsf_epoch_coverage/</code>.</p>
<p>The production code is <code>analysis/solarblue_all_fibers/run_all_fibers.py</code>; the LSF evaluator is <code>skysub/sky_decomp/lsf_spline2d.py</code>.
Coefficient products contain the individual fit and LSF for every physical fiber. Full diagnostic vectors are additionally retained for nine designated reference fibers per exposure.
For ordinary fibers, native wavelength grids, inverse variance (IVAR) and masks are retrieved from the source SFrame using the saved exposure and source-fiber row.</p>
<ul><li><code>build_coverage.py</code>: select exposure times, verify per-fiber Gaia photometry, read calibration epochs, apply channel-quality checks and count coverage.</li>
<li><code>complete_additional_fits.py</code>: fit 46 additional exposures using Gaia-clean working fibers.</li>
<li><code>analyze_stability.py</code>: evaluate saved LSF coefficients and calculate separate per-fiber within-night, within-lunation and sampled-night statistics.</li>
<li><code>publish_lunation_report.py</code>: export the clean samples and assemble the chronology, coverage maps, profile examples and numerical findings.</li>
<li><code>measurement_view.py</code>: publish individual-fiber measurements, ordered exposure positions within observed lunations, and the interactive Å plots.</li>
<li><code>test_coverage.py</code>, <code>test_stability.py</code>, <code>test_lunation_report.py</code>: check selection boundaries, sampling metrics, diagnostic equivalence and profile normalization.</li>
<li><code>README.md</code>: reproduction commands, estimator definitions, input provenance and report reviews.</li></ul>
<p><code>gaia_clean_fiber_sample.parquet</code>: {len(inputs):,} Gaia-clean working-fiber spectra before fit-quality checks.
<code>clean_lsf_sample.parquet</code>: {summary['accepted_channel_rows']:,} accepted channel measurements with source fit files and coefficient rows.
<code>fiber_exposure_channel.parquet</code>: all selected opportunities and their selection flags.</p></html>'''
    detail=OUT/'stability/lsf_epoch_lunation_stability.html'
    if 'Central-wavelength overview' in detail.read_text():
        (OUT/'stability/lsf_stability_diagnostics.html').write_text(detail.read_text())
    companion=re.sub(r'(href|src)="([^"#][^"]*)"',
        lambda m:m.group(1)+'="'+(m.group(2) if m.group(2).startswith(('http:','https:')) else '../'+m.group(2))+'"',content)
    detail.write_text(companion)
    (OUT/'lsf_epoch_coverage.html').write_text(content)
    (OUT/'lsf_epoch_coverage.html').chmod(0o644)
    print('Per-lunation inventory and complete sample checks passed',flush=True)


if __name__=='__main__':main()
