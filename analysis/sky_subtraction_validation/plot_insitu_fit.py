"""Separate interactive spectral diagnostics with their own residuals and legends."""
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from report_fibers import ROOT,page,trace,finish,UNIT
from skysub.sky_subtraction.reconstruction import diagnostic_hdul
from report_spectral_residuals import offset_residual
from report_contrast import contrast_rows,contrast_html



def write_fiber_diagnostic(h,drp,path,fiberid):
    """Write one page from an already reconstructed exposure view."""
    expnum=int(Path(path).stem.rsplit('-',1)[1])
    wave=np.asarray(h['WAVE'].data,float)
    row=int(np.flatnonzero(h['SLITMAP'].data['fiberid']==fiberid)[0])
    assert int(drp['SLITMAP'].data['fiberid'][row])==fiberid
    np.testing.assert_array_equal(wave,drp['WAVE'].data)
    values={name:np.asarray(h[name].data[row],float)*1e14 for name in ['OBSFLUX','FIT_MODEL','FIT_SKY','FIT_LINES','FIT_SOLAR','FIT_CONT','SKY','SKY_ML','SKY_OH','SKY_RES','FLUX']}
    np.testing.assert_allclose(values['FIT_MODEL'],values['FIT_SKY']+values['FIT_CONT'],rtol=2e-6,atol=1e-6,equal_nan=True)
    np.testing.assert_allclose(values['SKY'],values['SKY_ML']+values['SKY_OH']+values['SKY_RES'],rtol=2e-6,atol=1e-6,equal_nan=True)
    np.testing.assert_allclose(values['FLUX'],values['OBSFLUX']-values['SKY'],rtol=2e-6,atol=1e-6,equal_nan=True)
    valid=(h['MASK'].data[row]==0)&(h['IVAR'].data[row]>0)
    used=valid&~h['FIT_MASK'].data[row].astype(bool)&np.isfinite(values['FIT_MODEL'])
    sigma=np.divide(1e14,np.sqrt(np.maximum(h['IVAR'].data[row],0)),out=np.full(len(wave),np.nan),where=valid)
    qa=h['LSF_QA'].data;qa=qa[qa['input_row']==row]
    decision=' · '.join(f'{q["channel"]}: '+('in situ' if q['accepted'] else 'survey') for q in qa)
    body=f'<p><a href="00_MAIN_sky_subtraction_report.html">All exposures</a></p><h1>{expnum} · fiber {fiberid}</h1><p>Adopted LSF — B: survey · {decision}. σ(λ) = observed CFrame IVAR(λ)⁻¹ᐟ²; not sky-model error. Main residual offset: −5 RMS; RMS = √mean((data − fit)²) over valid fitted pixels.</p>'
    phot=h['GAIA'].data[row]
    body+=f'<p>Gaia G / observed G: {float(phot["ratio_gaia_to_lvm_flux_plus_sky_g"]):.3g}.</p>'
    if expnum==16819 and fiberid==489:
        blue=(wave<5787)&valid
        spike=int(np.nanargmax(np.where(blue,values['OBSFLUX'],np.nan)))
        body+=f'<p>Fiber-specific note: input Blue spike at {wave[spike]:.1f} Å; input MASK={int(h["MASK"].data[row,spike])}. This is not an exposure-wide label.</p>'
    if expnum==44808 and fiberid==1834:
        body+='<p>Fiber-specific note: previously inspected Blue stellar-contamination case; not an exposure-wide classification.</p>'
    if 'quality_bin_counts' in qa.names:
        body+='<h2>Full-channel line information</h2><p>Four equal wavelength regions per channel. Each needs ≥2 features with contrast ≥1 and peak S/N ≥5; features separated by ≥10 Å. Contrast = A/max(C, σA); S/N = A/σA. A: local fitted line amplitude; σA: input-noise amplitude error. C: larger observed or fitted baseline, using 75th percentiles in 4–10 Å sidebands. Medians include all measured features, not only passing features.</p><table><tr><th>Channel / region (Å)</th><th>Passing features</th><th>Median contrast</th><th>Median peak S/N</th><th>Channel decision</th></tr>'
        for q in qa:
            lo,hi=(5787.,7454.) if q['channel']=='R' else (7454.,9801.)
            hi=min(hi,float(wave[-1]+np.median(np.diff(wave))))
            edges=np.linspace(lo,hi,5)
            reasons=dict(insufficient_lines="insufficient numerical support",continuum_solver="continuum solver",line_information="insufficient line information",not_converged="iteration convergence",profile_guard="profile guard",insitu="accepted",unusable_survey_fallback="unusable survey reference")
            adopted="In situ" if q["accepted"] else ("Invalid sky" if q["reason"]=="unusable_survey_fallback" else "Survey")
            outcome=adopted+" — "+reasons.get(str(q["reason"]),str(q["reason"]))
            for j in range(4):
                contrast=q['quality_bin_contrast'][j];snr=q['quality_bin_snr'][j]
                contrast=f'{contrast:.2g}' if np.isfinite(contrast) else '—'
                snr=f'{snr:.2g}' if np.isfinite(snr) else '—'
                body+=f'<tr><td>{q["channel"]}: {edges[j]:.0f}–{edges[j+1]:.0f}</td><td>{int(q["quality_bin_counts"][j])}</td><td>{contrast}</td><td>{snr}</td><td>{outcome}</td></tr>'
        body+='</table>'
    body+='<p>The fitted sky combines emission lines with a solar/diffuse spectrum containing absorption features. The free continuum can absorb source light or sky-continuum mismatch; it is not a physical stellar-spectrum estimate. Cubic spline knots are spaced by 20 Å; flexibility can produce local bends. R/Z continua are fitted independently, so channel edges can jump.</p>'
    body+='<p>Strongest-peak summary below supplements the full-channel gate.</p>'
    body+=contrast_html(contrast_rows(wave,np.asarray(h['FIT_LINES'].data[row],float)*1e14,values['FIT_CONT'],np.asarray(h['FIT_SOLAR'].data[row],float)*1e14,used,values['OBSFLUX']))
    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=.10,row_heights=[.7,.3],
        subplot_titles=['Observed spectrum = fitted sky + nuisance continuum', '(Data − fit) / input flux error σ'])
    for key,name,color in [('FIT_LINES','Fitted sky emission lines','#1674aa'),('FIT_SOLAR','Fitted solar/diffuse spectrum','#6d5ca5'),('FIT_CONT','Nuisance continuum spline','#dc991d'),('OBSFLUX','Observed spectrum','#333333'),('FIT_MODEL','Sky + continuum fit','#c33b38')]:
        trace(fig,wave,np.where(valid,values[key],np.nan),name,1,color=color,width=.55 if key in ('FIT_LINES','FIT_SOLAR','FIT_CONT') else (1.5 if key=='FIT_MODEL' else 1.1))
    offset_residual(fig,wave,values['OBSFLUX'],values['FIT_MODEL'],sigma,used)
    residual=np.divide(values['OBSFLUX']-values['FIT_MODEL'],sigma)
    trace(fig,wave,np.where(used,residual,np.nan),'Fit residual / σ',2,color='#333333')
    fig.add_hrect(y0=-3,y1=3,fillcolor='#dce8ee',opacity=.5,line_width=0,row=2,col=1,layer='below')
    fig.add_hrect(y0=-1,y1=1,fillcolor='#8fb7cd',opacity=.55,line_width=0,row=2,col=1,layer='below')
    for level,label in [(1,'+1σ'),(-1,'−1σ')]:
        fig.add_hline(y=level,line_color='#35657e',line_width=1.4,line_dash='dash',row=2,col=1)
        fig.add_annotation(x=.98,xref='x2 domain',y=level,yref='y2',text=label,showarrow=False,xanchor='right',
            yanchor='bottom' if level>0 else 'top',yshift=7 if level>0 else -7,font=dict(color='#35657e',size=12),bgcolor='rgba(255,255,255,.85)')
    fig.add_hline(y=0,line_color='#555555',line_width=.8,row=2,col=1)
    fig.update_yaxes(title_text=UNIT,row=1,col=1);fig.update_yaxes(title_text='(Data − fit) / σ',row=2,col=1)
    fig.update_xaxes(title_text='Air wavelength (Å)',row=2,col=1)
    fig.update_xaxes(range=[float(wave[0]),float(wave[-1])])
    body+=finish(fig,'In-situ LSF fitting in R and Z channels',820)
    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=.10,row_heights=[.65,.35],
        subplot_titles=['Observed spectrum and actual SKY subtracted', 'Remaining SFrame FLUX: source + sky residual'])
    for key,name,color,dash in [('SKY_ML','ML sky before corrections','#1674aa','dash'),('SKY_OH','OH correction','#168257',None),('SKY_RES','Sky-arm correction','#dc991d',None),('OBSFLUX','Observed spectrum','#333333',None),('SKY','Final SKY','#993366',None)]:
        trace(fig,wave,np.where(valid,values[key],np.nan),name,1,color=color,dash=dash,width=.55 if key in ('SKY_ML','SKY_OH','SKY_RES') else (1.5 if key=='SKY' else 1.1))
    common=valid&(drp['MASK'].data[row]==0)&(drp['IVAR'].data[row]>0)
    error=np.where(common,sigma,np.nan)
    for y,name,fill in [(-error,'−1σ',None),(error,'Observed ±1σ','tonexty')]:
        fig.add_trace(go.Scatter(x=wave,y=y,mode='lines',name=name,showlegend=fill is not None,
            line=dict(color='#7ca2b8',width=.7),fill=fill,fillcolor='rgba(100,150,180,.25)',hoverinfo='skip'),row=2,col=1)
    for y,name,color in [(drp['FLUX'].data[row]*1e14,'DRP FLUX','#888888'),(values['FLUX'],'SKYML FLUX','#0078b7')]:
        trace(fig,wave,np.where(common,y,np.nan),name,2,color=color)
    fig.add_hline(y=0,line_color='#555555',line_width=.8,row=2,col=1)
    fig.update_yaxes(title_text=UNIT)
    fig.update_xaxes(title_text='Air wavelength (Å)',row=2,col=1)
    fig.update_xaxes(range=[float(wave[0]),float(wave[-1])])
    body+=finish(fig,'Final ML sky subtraction (B/R/Z)',820)
    output=ROOT/f'insitu_fit_{expnum}_{fiberid}.html'
    output.write_text(page(f'{expnum} / {fiberid}: fit and sky subtraction',body))
    assert '<details' not in body and '<summary' not in body
    print(output.name,flush=True)


def interactive(path,fiberid):
    expnum=int(Path(path).stem.rsplit("-",1)[1])
    with diagnostic_hdul(path,[fiberid]) as h,fits.open(Path(path).with_name(f"lvmSFrame-{expnum:08d}.fits")) as drp:
        write_fiber_diagnostic(h,drp,path,fiberid)

if __name__=='__main__':
    outputs={e['expnum']:e['output'] for e in json.loads((ROOT/'run_manifest_v2.json').read_text())}
    selected={(c['expnum'],c['fiberid']) for c in json.loads((ROOT/'fiber_cases.json').read_text())}
    selected.update(tuple(map(int,p.stem.split('_')[-2:])) for p in ROOT.glob('insitu_fit_*.html'))
    for expnum,fiberid in sorted(selected):
        interactive(outputs[expnum],fiberid)
