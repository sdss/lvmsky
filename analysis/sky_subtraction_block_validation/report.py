"""Continuous English report comparing sky-only blocks with native spectra and measured LSFs."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import json,sys,html
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'analysis/sky_subtraction_validation'))
from report_fibers import trace,finish,page,UNIT
from report_spectral_residuals import offset_residual
from report_lsf import surface_map,profile_comparison
from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
from skysub.sky_subtraction.pipeline import render,FACTOR,FAMILIES
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics,mspline_basis
from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
HERE=Path(__file__).resolve().parent
MEASUREMENTS=Path('/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/individual_lsf_measurements.parquet')
METRICS={'fwhm_angstrom':'FWHM (Å)','w50_angstrom':'W50 (Å)','w90_angstrom':'W90 (Å)','centroid_angstrom':'Centroid (Å)'}


def spectrum(title,wave,observed,model,components,sigma,valid,source=False,fragment=False):
    good=valid&np.isfinite(sigma)&(sigma>0)&np.isfinite(model)
    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,row_heights=[.70,.30],vertical_spacing=.12)
    labels={'oh':'OH airglow','moon':'Scattered moonlight','zodi':'Zodiacal light','diffuse':'Diffuse airglow continuum','atom':'Atomic airglow','orc':'Oxygen recombination','o2':'O₂ airglow'}
    additive={k:v for k,v in components.items() if k in FAMILIES or k in ('OH correction','B residual transfer')}
    np.testing.assert_allclose(sum(additive.values())*FACTOR,model*FACTOR,rtol=1e-8,atol=1e-8)
    for key,value in additive.items():trace(fig,wave,np.where(valid,value*FACTOR,np.nan),labels.get(key,key),1,width=.65)
    trace(fig,wave,np.where(valid,observed*FACTOR,np.nan),'Observed',1,color='#333333',width=1)
    trace(fig,wave,np.where(valid,model*FACTOR,np.nan),'Applied sky' if source else 'Best fit',1,color='#dc4a2c',width=1.4)
    offset=offset_residual(fig,wave,observed*FACTOR,model*FACTOR,sigma*FACTOR,good)
    if source:fig.data[-1].name='Source + sky residual (−5 RMS)'
    for level in [-1,1]:
        trace(fig,wave,np.where(good,level,np.nan),'−1σ' if level<0 else '+1σ',2,color='#9ebbd2',width=.6)
    fig.data[-1].fill='tonexty';fig.data[-1].fillcolor='rgba(135,175,210,.17)'
    for level in [-1,1]:
        fig.add_annotation(x=1,xref='paper',y=level,yref='y2',text='−1σ' if level<0 else '+1σ',showarrow=False,xanchor='left',xshift=3,yshift=level*8,font=dict(size=10))
    trace(fig,wave,np.where(good,0,np.nan),'Zero',2,color='#777777',width=.7,showlegend=False)
    residual=np.divide(observed-model,sigma,out=np.full(len(wave),np.nan),where=good)
    trace(fig,wave,np.where(good,residual,np.nan),'Source + sky residual / σ' if source else '(Data − model) / σ',2,color='#34526e',width=.7)
    fig.update_xaxes(title_text='Air wavelength (Å)',row=2,col=1,range=[float(wave[0]),float(wave[-1])])
    fig.update_yaxes(title_text=UNIT,row=1,col=1)
    fig.update_yaxes(title_text='Residual / σ',row=2,col=1)
    fig.add_annotation(x=0,y=1.045,xref='paper',yref='paper',text='Observed spectrum and components; flux residual at −5 RMS',showarrow=False,xanchor='left',font=dict(size=12))
    fig.add_annotation(x=0,y=.30,xref='paper',yref='paper',text='Input-noise units; shaded ±1σ with reference lines',showarrow=False,xanchor='left',font=dict(size=12))
    body='<h1>'+html.escape(title)+'</h1><p>σ = 1/√IVAR from the input CFrame (median IVAR for median diagnostics). The thin black curves show input error around the −5 RMS offset; RMS is unweighted over displayed valid pixels. '
    if source:body+='Remaining FLUX contains source emission plus sky residual; its normalized values are not a zero-target fit test. '
    body+='All wavelengths and outliers remain available; drag to zoom, use the toolbar to reset. Input masks apply; protected-line masks also apply to Sky/median fits. Plots retain the channel-edge pixels trimmed from the QP. Units: 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹.</p>'
    content=body+finish(fig,title,860)
    if fragment:return content.replace('<h1>','<h3>').replace('</h1>','</h3>')
    return page(title,content+'<p><a href="00_MAIN_sky_block_report.html">Main report</a></p>')


def maps(r,row):
    e=row['expnum'];science=np.asarray(r.h['LSF_FIBERS'].data['input_row'],int);slit=r.h['SLITMAP'].data
    with fits.open(row['drp'],memmap=True) as drp:
        np.testing.assert_array_equal(drp['WAVE'].data,r.wave)
        np.testing.assert_array_equal(drp['SLITMAP'].data['fiberid'],slit['fiberid'])
        from astropy import units as u
        scale=u.Unit(drp['FLUX'].header['BUNIT']).to(u.erg/u.s/u.cm**2/u.AA)
        a=np.asarray(drp['FLUX'].data[science],float)*scale*FACTOR
        b=np.asarray(r.h['FLUX'].data[science],float)*FACTOR
        valid=(drp['MASK'].data[science]==0)&(r.h['MASK'].data[science]==0)&(drp['IVAR'].data[science]>0)&(r.h['IVAR'].data[science]>0)
    fig,axes=plt.subplots(3,2,figsize=(15,11),layout='constrained')
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#888888')
    map_rows=[]
    for j,(ch,lo,hi) in enumerate([('B',3600,5787),('R',5787,7454),('Z',7454,9801)]):
        take=(r.wave>=lo)&(r.wave<hi);common=valid[:,take]
        limit=float(np.nanpercentile(np.abs(np.r_[a[:,take][common],b[:,take][common]]),98))
        for k,values in enumerate([a,b]):
            image=np.full((len(slit),take.sum()),np.nan)
            image[science]=np.where(common,values[:,take],np.nan)
            ax=axes[j,k];im=ax.imshow(image,aspect='auto',origin='lower',extent=[r.wave[take][0],r.wave[take][-1],slit['fiberid'].min()-.5,slit['fiberid'].max()+.5],cmap=cmap,vmin=-limit,vmax=limit,interpolation='nearest')
            ax.set_title(ch+' · '+('DRP SFrame FLUX' if k==0 else 'Block SFrame FLUX'));ax.set_ylabel('Physical fiber ID');ax.set_xlabel('Air wavelength (Å)')
            for sp in [1,2]:ax.axhline(np.max(slit['fiberid'][slit['spectrographid']==sp])+.5,color='black',lw=.6)
            fig.colorbar(im,ax=ax,label='Flux (10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹)',shrink=.82)
        map_rows.append(dict(expnum=e,channel=ch,shared_abs_limit=limit,valid_pixels=int(common.sum()),
            drp_median_abs_flux=float(np.median(abs(a[:,take][common]))),block_median_abs_flux=float(np.median(abs(b[:,take][common])))))
    fig.suptitle(f'Exposure {e}: remaining source + sky residual',fontsize=17)
    fig.savefig(HERE/f'flux_maps_{e}.png',dpi=125);plt.close(fig)
    return map_rows


def lsf_comparison(r,row,measured):
    """Match exact exposure, physical fiber and channel to accepted individual measurements."""
    e=row['expnum'];m=measured[(measured.expnum==e)&(measured.telescope=='Sci')&measured.channel.isin(['R','Z'])]
    records=[]
    for f,group in m.groupby('fiberid'):
        state=r.global_state(int(f))
        for _,point in group.iterrows():
            values=evaluate_lsf_diagnostics(state,np.array([point.wavelength_angstrom]))
            record=dict(expnum=e,fiberid=int(f),spectrographid=int(point.spectrographid),channel=point.channel,
                wavelength_A=float(point.wavelength_angstrom),lunation=int(point.lunation_number))
            for key in METRICS:
                reference=float(point[key]);block=float(values[key][0]);record[key+'_measured']=reference;record[key+'_block']=block;record[key+'_difference']=block-reference
            records.append(record)
    if not records:return records,[]
    frame=pd.DataFrame(records);summaries=[]
    fig,axes=plt.subplots(4,2,figsize=(11,12),sharex=True,layout='constrained')
    for j,ch in enumerate('RZ'):
        c=frame[frame.channel==ch].sort_values('fiberid')
        for i,(key,label) in enumerate(METRICS.items()):
            ax=axes[i,j];ax.scatter(c.fiberid,c[key+'_measured'],s=5,label='Individual measured LSF',alpha=.7)
            ax.scatter(c.fiberid,c[key+'_block'],s=4,label='Sky-only block',alpha=.8)
            ax.set_ylabel(label);ax.grid(alpha=.2)
            if i==0:ax.set_title(ch+f' · {len(c)} matched fibers')
            if i==3:ax.set_xlabel('Physical fiber ID')
            if len(c):summaries.append(dict(expnum=e,channel=ch,metric=key,n_fibers=len(c),median_difference=float(c[key+'_difference'].median()),median_abs_difference=float(c[key+'_difference'].abs().median())))
    axes[0,1].legend(loc='upper center',bbox_to_anchor=(-.1,1.36),ncol=2,fontsize=10)
    fig.suptitle(f'Exposure {e}: R/Z block LSF and per-fiber Science fits at one wavelength per channel\nOriginal centroids retained',fontsize=16)
    fig.savefig(HERE/f'lsf_comparison_{e}.png',dpi=125);plt.close(fig)
    return records,summaries


def exposure(row):
    e=row['expnum'];print('Reporting',e,flush=True)
    measured=pd.read_parquet(MEASUREMENTS)
    with SkyBlockReconstructor(row['output'],cframe=row['input']) as r:
        science=r.h['LSF_FIBERS'].data;blocks=Table(r.h['BLOCK_LSF_META'].data).to_pandas()
        section=f'<section id="exp{e}"><h2>Exposure {e}</h2><p>{html.escape(row.get("purpose","Independent Science-LSF comparison").split(";")[0])} · UTC {html.escape(str(r.h[0].header["OBSTIME"]))} · {len(science)} Science fibers. <a href="products/{Path(row["output"]).name}">Compact SFrame</a></p>'
        selected_rows=r.h['SKY_FIBER_COEF'].data
        links=[];roundtrips=[]
        for sp in [1,2,3]:
            rows=selected_rows[selected_rows['spectrograph']==sp];skyrow=rows[len(rows)//2]
            f=int(skyrow['fiberid']);d=r.sky_fiber(f)
            name=f'sky_fit_{e}_{f}.html';(HERE/name).write_text(spectrum(f'Exposure {e} · SP{sp} · {d["telescope"]} fiber {f}',r.wave,d['observed'],d['model'],d['components'],d['sigma'],d['valid']))
            links.append(f'<a href="{name}">SP{sp} Sky fiber {f}</a>')
        for label in r.labels:
            j=r.labels.index(label);observed=r.h['MED_FLUX'].data[j];weight=r.h['MED_IVAR'].data[j]
            sigma=np.divide(1,np.sqrt(np.maximum(weight,0)),out=np.full(len(weight),np.nan),where=weight>0)
            model,components=r.medians[label];valid=(weight>0)&np.isfinite(observed)&~r.h['MED_FIT_MASK'].data.astype(bool)
            name=f'median_{e}_{label}.html';(HERE/name).write_text(spectrum(f'Exposure {e} · {label} preparation median',r.wave,observed,model,components,sigma,valid))
            links.append(f'<a href="{name}">{label} median</a>')
        ids=[680,1447,971] if e==59750 else [int(science['fiberid'][np.flatnonzero(r.h['SLITMAP'].data['spectrographid'][science['input_row']]==sp)[len(np.flatnonzero(r.h['SLITMAP'].data['spectrographid'][science['input_row']]==sp))//2]]) for sp in [1,2,3]]
        for f in ids:
            values=r.sky(f);index=int(np.flatnonzero(r.h['SLITMAP'].data['fiberid']==f)[0])
            state=r.global_state(f);state_index=int(np.flatnonzero(science['fiberid']==f)[0])
            for ch in 'BRZ':np.testing.assert_allclose(state.coefficients[ch],r.state('',state_index).coefficients[ch],atol=1e-13)
            _,components,_=render(r.model_sci,r.coef,state)
            components={k:v/FACTOR for k,v in components.items()};components['OH correction']=values['SKY_OH'];components['B residual transfer']=values['SKY_RES']
            sigma=np.divide(1,np.sqrt(np.maximum(r.ivar[index],0)),out=np.full(len(r.wave),np.nan),where=r.ivar[index]>0)
            valid=(r.mask[index]==0)&np.isfinite(r.observed[index])&(r.ivar[index]>0)
            name=f'science_{e}_{f}.html';(HERE/name).write_text(spectrum(f'Exposure {e} · SP{int(r.h["SLITMAP"].data["spectrographid"][index])} · Science fiber {f}',r.wave,r.observed[index],values['SKY'],components,sigma,valid,source=True))
            links.append(f'<a href="{name}">Science fiber {f}</a>');roundtrips.append(f)
        section+='<p>'+' · '.join(links)+'</p><p>Sky examples: middle selected physical fiber in each spectrograph. Science examples: middle processed fiber in each spectrograph; exposure 59750 retains previously inspected fibers 680, 1447, 971. Preparation medians keep their original LSFs.</p>'
        section+='<table><tr><th>Block (R/Z only)</th><th>Sky fibers</th><th>χ² / pixel</th><th>Slit rank</th><th>Selected fiber range</th><th>Normalization / centroid error</th></tr>'
        for _,b in blocks.iterrows():section+=f'<tr><td>SP{b.spectrograph} {b.channel}</td><td>{b.n_fibers}</td><td>{b.chi2_per_pixel:.3f}</td><td>{b.slit_rank}</td><td>{b.selected_min:g}–{b.selected_max:g}</td><td>{b.normalization_error:.1e} / {b.centroid_error:.1e} Å</td></tr>'
        section+='</table>'
        surface_map(r,e)
        section+=f'<img src="global_lsf_maps_{e}.png" alt="Full wavelength centered LSF maps for exposure {e}"><p>Global surfaces on a regular sampling of the full wavelength and pseudo-slit domains. Colors show FWHM, equal-tail W50/W90, centroid and quantile asymmetry; each panel uses its own labeled scale. Scales are independent across panels and exposures. The centroid panel shows numerical centering error, using a symmetric half-range of at least 10⁻⁶ Å.</p>'
        section+=profile_comparison(r,e,measured)
        flux_rows=maps(r,row);records,summaries=lsf_comparison(r,row,measured)
        section+=f'<img src="flux_maps_{e}.png" alt="Full wavelength Science flux maps for exposure {e}"><p>Each channel pair shares a symmetric scale at the pooled 98th absolute percentile of common valid pixels; saturated values remain in the FITS. Gray gaps are masked or non-Science rows. Physical fiber order and spectrograph boundaries follow SLITMAP. These maps show remaining source plus sky residual, not a pure sky-fit residual.</p>'
        if records:
            section+=f'<img src="lsf_comparison_{e}.png" alt="Measured and block LSF metrics for exposure {e}"><p>Accepted individual Science profiles from the lunation analysis, matched by exposure, physical fiber and wavelength. W50/W90 enclose 50%/90% of profile probability. Numerical widths and the original centroid are compared in Å; these independent fits are empirical references rather than ground truth. <a href="http://lumen.abudhabi.nyu.edu/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/lsf_epoch_coverage.html">Reference selection and acceptance methodology</a>.</p>'
        section+='</section>'
        integrity=dict(expnum=e,roundtrip_fiberids=roundtrips,global_blocks=len(blocks),science_fibers=len(science),
            normalization_max=float(blocks.normalization_error.max()),centroid_max_A=float(blocks.centroid_error.max()),
            flux_identity_max=float(np.nanmax(abs(r.h['FLUX'].data[science['input_row']]+r.h['SKY'].data[science['input_row']]-r.observed[science['input_row']]))))
    (HERE/f'section_{e}.html').write_text(section)
    (HERE/f'metrics_{e}.json').write_text(json.dumps(dict(integrity=integrity,flux=flux_rows,lsf_summary=summaries,lsf_records=records),indent=2))
    return integrity


def assemble(rows):
    body='<h1>Sky-only block LSF reconstruction</h1><p>LSF = line-spread function; ML = machine learning; QP = quadratic program; IVAR = inverse observed-noise variance. One centered LSF surface per spectrograph for each of R and Z, inferred jointly from clean SkyE/SkyW fibers. Science spectra contribute to median preparation and OH correction; the block LSF fit uses Sky fibers exclusively. B is not fitted: each Science fiber keeps its survey LSF there.</p>'
    body+='<p><a href="README_method.md">Method and CLI</a> · <a href="run_manifest.json">Input/output manifest</a> · <a href="lsf_differences.csv">Matched LSF measurements</a> · <a href="product_validation.json">All-fiber checks</a></p>'
    body+='<p><b>Model:</b> normalized cubic M-splines in offset × four cubic wavelength controls × four cubic pseudo-slit controls. Nonnegative coefficients, unit integral and zero centroid hold continuously. Each Sky fiber has its own frozen decomposition vector. The input wavelength solution is retained. No additional width penalty or spatial regularization is imposed.</p>'
    body+='<p><b>Comparison:</b> the five previous exposures reuse exactly their saved medians, ML predictions and B/OH corrections. Three additional exposures come from different sampled lunations with dense accepted Science-LSF coverage. Individual Science fits test transfer from Sky positions; profile centroids remain in their original wavelength convention. Fixed component/calibration errors and input wavelength offsets can contribute to disagreements; solver success is a feasibility check.</p>'
    body+='<p><b>Noise and selection:</b> Gaia stellar flux / observed total fiber flux in the Gaia G passband &lt;0.01, hardware-good Sky fibers; native input IVAR, input/protected-line masks, and 5 Å channel-edge trimming in QP. σ = 1/√input IVAR; uncertainties are observed noise, with no propagated coefficient/LSF covariance. Selected-range endpoints in each table delimit data support; Science positions outside that interval test extrapolation. Slit rank four verifies algebraic coverage, not uniform information. χ²/pixel uses trimmed Sky fitting pixels and frozen coefficients estimated from the same spectra; it is neither reduced χ² nor held-out accuracy. Exposure 16819 Z has χ²/pixel ≈11–13: substantial spectral mismatch persists despite feasible moment constraints.</p>'
    decomposition_section=HERE/'fiber_decomposition_section.html'
    if decomposition_section.exists():body+=decomposition_section.read_text()
    fiber_section=HERE/'fiber_profiles_section.html'
    if fiber_section.exists():body+=fiber_section.read_text()
    comparison_section=HERE/'individual_comparison_section.html'
    if comparison_section.exists():body+=comparison_section.read_text()
    sky_section=HERE/'sky_diagnostics_section.html'
    if sky_section.exists():body+=sky_section.read_text()
    residual_section=HERE/'residual_budget_section.html'
    if residual_section.exists():body+=residual_section.read_text()
    body+='<p>'+' · '.join(f'<a href="#exp{r["expnum"]}">{r["expnum"]}</a>' for r in rows)+'</p>'
    allmetrics=[json.loads((HERE/f'metrics_{r["expnum"]}.json').read_text()) for r in rows]
    summaries=[x for m in allmetrics for x in m['lsf_summary']]
    if summaries:
        body+='<h2>Matched Science LSF comparison</h2><p>Signed difference = block − individual fit. W50 = q75−q25 and W90 = q95−q05 are equal-tail profile widths. Median absolute differences use matched accepted fibers; counts are per channel, at each measurement’s recorded wavelength.</p><table><tr><th>Exposure</th><th>Channel</th><th>Fibers</th><th>|ΔFWHM| median (Å)</th><th>|ΔW50| median (Å)</th><th>|ΔW90| median (Å)</th><th>|Δcentroid| median (Å)</th></tr>'
        for e,ch in sorted(set((s['expnum'],s['channel']) for s in summaries)):
            group={s['metric']:s for s in summaries if s['expnum']==e and s['channel']==ch}
            body+=f'<tr><td>{e}</td><td>{ch}</td><td>{group["fwhm_angstrom"]["n_fibers"]}</td>'+''.join(f'<td>{group[key]["median_abs_difference"]:.4f}</td>' for key in METRICS)+'</tr>'
        body+='</table>'
    for row in rows:body+=(HERE/f'section_{row["expnum"]}.html').read_text()
    assert '<details' not in body and '<summary' not in body
    assert body.count('<section id=')==len(rows)
    (HERE/'00_MAIN_sky_block_report.html').write_text(page('Sky-only block LSF reconstruction',body))
    (HERE/'integrity.json').write_text(json.dumps([m['integrity'] for m in allmetrics],indent=2))
    pd.DataFrame([x for m in allmetrics for x in m['lsf_records']]).to_csv(HERE/'lsf_differences.csv',index=False)
    pd.DataFrame(summaries).to_csv(HERE/'lsf_summary.csv',index=False)
    pd.DataFrame([x for m in allmetrics for x in m['flux']]).to_csv(HERE/'flux_summary.csv',index=False)
    (HERE/'plotly.min.js').write_text(get_plotlyjs())
    (HERE/'README_method.md').write_text((ROOT/'skysub/sky_subtraction_block/README.md').read_text())


if __name__=='__main__':
    rows=json.loads((HERE/'run_manifest.json').read_text())
    if '--assemble-only' not in sys.argv:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(4) as pool:list(pool.map(exposure,rows))
    assemble(rows)
