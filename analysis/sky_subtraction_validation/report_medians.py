"""Reconstructed Science-median decomposition and matched-fiber subtraction stages."""
import json
import numpy as np
from plotly.subplots import make_subplots
from report_spectral_residuals import offset_residual


def science_median_data(h,drp):
    from skysub.sky_subtraction.io import median_spectrum
    labels=list(h['MED_LABELS'].data['telescope']);index=labels.index('Sci')
    rows=np.flatnonzero(h['GAIA'].data['selected']&(h['GAIA'].data['telescope']=='Sci'))
    assert len(rows)>0
    np.testing.assert_array_equal(h['SLITMAP'].data['fiberid'],drp['SLITMAP'].data['fiberid'])
    np.testing.assert_array_equal(h['WAVE'].data,drp['WAVE'].data)
    observed=np.asarray(h['MED_FLUX'].data[index],float)
    rebuilt,_,_=median_spectrum(h['OBSFLUX'].data,h['IVAR'].data,h['MASK'].data,rows)
    np.testing.assert_allclose(rebuilt,observed,rtol=2e-6,atol=1e-21,equal_nan=True)
    components=np.asarray(h['MED_COMPONENTS'].data[index],float)
    model=np.asarray(h['MED_MODEL'].data[index],float)
    np.testing.assert_allclose(components.sum(axis=0),model,rtol=2e-6,atol=1e-21,equal_nan=True)
    np.testing.assert_allclose(h['OBSFLUX'].data[rows]-h['SKY'].data[rows],h['FLUX'].data[rows],rtol=2e-6,atol=1e-21,equal_nan=True)
    data=dict(wave=np.asarray(h['WAVE'].data,float),n_fibers=len(rows),observed=observed*1e14,model=model*1e14,
              components=components*1e14,names=list(h['COMPONENT_NAMES'].data['family']),
              fit_mask=h['MED_FIT_MASK'].data.astype(bool),fit_ivar=np.asarray(h['MED_IVAR'].data[index],float)/1e28)
    valid=(h['MASK'].data==0)&(drp['MASK'].data==0)&(h['IVAR'].data>0)&(drp['IVAR'].data>0)
    stages={'DRP':drp['FLUX'].data,'ML before corrections':h['OBSFLUX'].data-h['SKY_ML'].data,
            'ML + OH':h['OBSFLUX'].data-h['SKY_ML'].data-h['SKY_OH'].data,'Final SKYML':h['FLUX'].data}
    for values in stages.values():valid &= np.isfinite(values)
    mask=(~valid).astype('u1')
    data['stages']={}
    for name,values in stages.items():
        median,ivar,count=median_spectrum(values,h['IVAR'].data,mask,rows)
        data['stages'][name]=median*1e14
    data['sigma']=np.divide(1e14,np.sqrt(ivar),out=np.full_like(ivar,np.nan),where=ivar>0)
    if 'SKYML_META' in h:
        data['oh_correction']=json.loads(bytes(h['SKYML_META'].data).decode())['oh_correction']
    return data


def oh_status(info):
    from html import escape
    decision='applied' if info['accepted'] else 'ML retained — '+str(info.get('reason','correction quality guard'))
    return '<p><b>OH correction:</b> '+escape(decision)+'.</p>'


def science_median_metrics(data):
    """Native-noise fit residuals and matching-pixel subtraction-stage scatter."""
    from validate import highpass
    wave=data['wave'];residual=data['observed']-data['model']
    valid=~data['fit_mask']&(data['fit_ivar']>0)&np.isfinite(residual)
    common=valid.copy()
    for spectrum in data['stages'].values():common &= np.isfinite(spectrum)
    filtered={name:highpass(spectrum[None,:],common[None,:])[0] for name,spectrum in data['stages'].items()}
    output=[]
    for channel,lo,hi in [('B',3600,5787),('R',5787,7454),('Z',7454,9801)]:
        band=(wave>=lo)&(wave<hi);used=valid&band
        if not used.any():continue
        take=common&band
        output.append(dict(channel=channel,rms=float(np.sqrt(np.mean(residual[used]**2))),
            normalized_rms=float(np.sqrt(np.mean(residual[used]**2*data['fit_ivar'][used]))),
            median_sigma=float(np.median(1/np.sqrt(data['fit_ivar'][used]))),
            stages={name:float(np.sqrt(np.mean(values[take]**2))) if take.any() else None for name,values in filtered.items()}))
    return output


def median_metrics_table(data):
    body='<p>RMS in 10⁻¹⁴ flux units. Fit RMS = √mean[(data−fit)²] over valid fitted pixels. Fit σ = MED_IVAR⁻¹ᐟ²: propagated input error of the median, not residual scatter or model uncertainty. Subtraction RMS: 51-pixel masked mean removed, same valid pixels; includes source structure.</p>'
    body+='<table><tr><th>Channel</th><th>Fit RMS</th><th>Median input σ</th><th>√mean[(data−fit)²/σ²]</th><th>DRP RMS</th><th>ML RMS</th><th>+OH RMS</th><th>Final RMS</th></tr>'
    for row in science_median_metrics(data):
        values=[row['rms'],row['median_sigma'],row['normalized_rms'],*row['stages'].values()]
        body+='<tr><td>'+row['channel']+'</td>'+''.join('<td>'+('—' if x is None else f'{x:.3g}')+'</td>' for x in values)+'</tr>'
    return body+'</table>'


def science_median_plots(data,expnum):
    from report_fibers import trace,finish,UNIT
    wave=data['wave']
    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=.10,row_heights=[.7,.3],
        subplot_titles=['Observed Science median and reconstructed decomposition','Decomposition residual / approximate input median σ'])
    colors=['#9955aa','#dc991d','#997433','#777777','#2585ba','#71a43d','#cb6485']
    for name,values,color in zip(data['names'],data['components'],colors):trace(fig,wave,values,name,1,color=color,width=.55)
    trace(fig,wave,data['observed'],'Observed median',1,color='#222222',width=1.1)
    trace(fig,wave,data['model'],'Science decomposition',1,color='#c33b38',width=1.5)
    used=~data['fit_mask']&(data['fit_ivar']>0)
    sigma=np.divide(1,np.sqrt(data['fit_ivar']),out=np.full_like(data['fit_ivar'],np.nan),where=data['fit_ivar']>0)
    offset_residual(fig,wave,data['observed'],data['model'],sigma,used)
    residual=(data['observed']-data['model'])*np.sqrt(np.maximum(data['fit_ivar'],0))
    trace(fig,wave,np.where(used,residual,np.nan),'Fit residual / σ',2,color='#333333')
    fig.add_hrect(y0=-1,y1=1,fillcolor='#dce8ee',opacity=.6,line_width=0,row=2,col=1,layer='below')
    for level,label in [(1,'+1σ'),(-1,'−1σ')]:
        fig.add_hline(y=level,line_color='black',line_width=.8,row=2,col=1)
        fig.add_annotation(x=1,xref='x2 domain',y=level,yref='y2',text=label,showarrow=False,xanchor='right',
            yanchor='bottom' if level>0 else 'top',yshift=7 if level>0 else -7)
    fig.add_hline(y=0,line_color='#777777',line_width=.5,row=2,col=1)
    fig.update_yaxes(title_text=UNIT,row=1,col=1);fig.update_yaxes(title_text='(Data − model) / σ',row=2,col=1)
    fig.update_xaxes(title_text='Air wavelength (Å)',row=2,col=1)
    fig.update_xaxes(range=[float(wave[0]),float(wave[-1])])
    body=median_metrics_table(data)
    body+='<p>Recorded MED_FLUX and coefficient-reconstructed model; residuals use MED_FIT_MASK. σ = MED_IVAR⁻¹ᐟ², approximate input-median error. Main residual offset: −5 RMS of fitted data−model.</p>'
    body+=finish(fig,f'{expnum} · Science median decomposition · {data["n_fibers"]} selected fibers',850)
    fig=make_subplots(rows=1,cols=1)
    for name,color in [('DRP','#888888'),('ML before corrections','#dc991d'),('ML + OH','#168257'),('Final SKYML','#0078b7')]:
        trace(fig,wave,data['stages'][name],name,1,color=color)
    for sign in [-1,1]:trace(fig,wave,sign*data['sigma'],'Approx. input median ±1σ',1,color='black',showlegend=sign==1)
    fig.add_hline(y=0,line_color='#777777',line_width=.5,row=1,col=1)
    fig.update_yaxes(title_text=UNIT);fig.update_xaxes(title_text='Air wavelength (Å)',range=[float(wave[0]),float(wave[-1])])
    body+='<p>ML → +SKY_OH → +SKY_RES (final): per-fiber subtraction before the common-valid-pixel median. σ: independent Gaussian input-median approximation; source scatter and sky-model errors uncalibrated.</p>'
    body+=finish(fig,f'{expnum} · Science median after subtraction: source + sky residual',630)
    if 'oh_correction' in data:body+=oh_status(data['oh_correction'])
    from audit_median_fits import median_line_probes,audit_summary
    lines=sum(values for name,values in zip(data['names'],data['components']) if name in ('oh','atom','orc','o2'))
    probes=median_line_probes(wave,data['observed'],lines,data['model'],np.where(used,data['fit_ivar'],0))
    body+=audit_summary({'Sci':{'line_probes':probes}},expnum)
    return body
