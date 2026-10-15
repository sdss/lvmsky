"""Whole-exposure FLUX maps and explicit LSF iteration profiles."""
from copy import deepcopy
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from astropy.io import fits
from skysub.sky_subtraction.reconstruction import diagnostic_hdul
from validate import ROOT, highpass
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_density,evaluate_lsf_diagnostics
from report_medians import science_median_data,science_median_plots
from report_cases import case_table,diagnostic_cases
from report_lsf_summary import fallback_summaries,fallback_table

CHANNELS=[('B',3600,5787),('R',5787,7454),('Z',7454,9801)]
UNIT='Flux (10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹)'
WINDOW_GROUPS={
    'blue':[('[O II] 3727',3600,4050),('Ca II H/K',3800,4300),('[O I] 5577',5350,5750)],
    'red':[('[O I] 6300/6364',6100,6500),('Hα → R-channel end',6500,7454),('OH + telluric B',6650,7454)],
    'far_red':[('Telluric A + OH',7454,8250),('8300 Å + Ca II triplet',8100,8900),('Red OH + telluric H₂O',8950,9800)],
}


def lsf_selection_diagram():
    """Show the production decision, independently for each Science fiber/channel."""
    return '''<section id="lsf-selection" class="card" style="text-align:center">
<h2>Which line-spread function (LSF) is used for this fiber?</h2>
<p><b>B → survey LSF</b> · <b>R and Z → independent decisions below</b></p>
<div class="notice"><b>Prepare inputs — no fitting yet</b><br>Read observed CFrame spectrum, input errors and masks; load survey LSF as the initial profile.</div>
<div aria-hidden="true">↓</div>
<div><b>Check numerical support</b><br>Valid pixels and wavelength-basis coverage; this is a solver precondition, not a quality test.</div>
<div aria-hidden="true">↓ yes</div>
<div class="notice"><b>Fit starts here: alternate continuum/amplitudes ↔ LSF</b><br>Model: sky emission + solar/diffuse spectrum + continuum spline ≥0 (20 Å knot spacing).<br>Up to 12 iterations. Convergence: W50/W90 (50%/90% enclosed-flux widths) step &lt;1%; centroid step &lt;0.015 Å.</div>
<div aria-hidden="true">↓</div>
<div class="notice"><b>Do observed sky lines constrain the full channel?</b><br>Four equal wavelength regions; each needs ≥2 independent sky features (≥10 Å apart).<br>Each counted feature: local line-amplitude contrast ≥1 and line-amplitude S/N ≥5.<br>Local weighted line + continuum + slope fit gives amplitude A and input-noise error σA; S/N = A/σA.<br>Contrast = A/max(C, σA); C is the larger observed or fitted baseline in 4–10 Å sidebands (75th percentiles).</div>
<div aria-hidden="true">↓</div>
<div><b>Converged + profile guards satisfied?</b><br>At 9 wavelength probes: FWHM 0.9–3 Å; W50/W90 change ≤35% from survey, ≤50% from Science median; centroid change ≤0.25 Å from survey.</div>
<div aria-hidden="true">↙　　　　　　　　　↘</div>
<div style="display:grid;grid-template-columns:1fr 1fr;gap:12px">
<div style="background:#e9f4ed;padding:12px"><b>All checks pass → in-situ LSF</b></div>
<div style="background:#edf1f6;padding:12px"><b>Insufficient line information, solver issue,<br>or any guard fails → measured survey LSF</b><br>Survey reference unusable → channel marked invalid.</div>
</div>
<p>Adopted LSF → common Science-IFU ML coefficients + median-derived OH correction + Blue sky-arm residual correction → SKY; FLUX = observed − SKY.<br>Local spline and local fit amplitudes are used only to estimate LSF.</p>
<p style="font-size:14px">Conservative contrast/S/N operating points; synthetic tests check clean, bright-continuum, low-S/N and incomplete-coverage cases. Passing these guards does not establish independently calibrated LSF accuracy.</p>
</section>'''


def flux_arrays(h, drp):
    """Keep every science row and use identical valid pixels in both images."""
    slit=h['SLITMAP'].data
    rows=np.flatnonzero(np.char.strip(slit['telescope'].astype(str))=='Sci')
    rows=rows[np.argsort(slit['fiberid'][rows])]
    assert np.array_equal(slit['fiberid'],drp['SLITMAP'].data['fiberid'])
    np.testing.assert_array_equal(h['WAVE'].data,drp['WAVE'].data)
    good=(h['MASK'].data[rows]==0)&(drp['MASK'].data[rows]==0)&(h['IVAR'].data[rows]>0)&(drp['IVAR'].data[rows]>0)
    arrays=[np.asarray(x['FLUX'].data[rows],float)*1e14 for x in (drp,h)]
    good &= np.isfinite(arrays[0])&np.isfinite(arrays[1])
    return rows,good,[np.where(good,a,np.nan) for a in arrays]


def flux_map(h,drp,expnum,median,window_group=None):
    rows,good,arrays=flux_arrays(h,drp)
    wave=np.asarray(h['WAVE'].data,float)
    bands=WINDOW_GROUPS[window_group] if window_group else CHANNELS
    fig=plt.figure(figsize=(16,15),layout='constrained')
    grid=fig.add_gridspec(6,2,height_ratios=[3,1,3,1,3,1])
    axes=np.array([[fig.add_subplot(grid[2*j,col]) for col in range(2)] for j in range(3)])
    limits=[]
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#b8bcc2')
    for j,(ch,lo,hi) in enumerate(bands):
        band=(wave>=lo)&(wave<hi)
        pair=[a[:,band] for a in arrays]
        values=np.concatenate([a[np.isfinite(a)] for a in pair])
        limit=max(float(np.percentile(np.abs(values),98)),1e-12)
        limits.append(limit)
        for col,a in enumerate(pair):
            image=axes[j,col].imshow(a,origin='lower',aspect='auto',interpolation='nearest',cmap=cmap,vmin=-limit,vmax=limit,
                extent=(wave[band][0]-.25,wave[band][-1]+.25,-.5,len(rows)-.5),rasterized=True)
            axes[j,col].set_title(f'{ch} · '+('Standard DRP SFrame FLUX' if col==0 else 'SKYML SFrame FLUX'))
            ticks=np.linspace(0,len(rows)-1,6,dtype=int)
            axes[j,col].set_yticks(ticks,labels=h['SLITMAP'].data['fiberid'][rows[ticks]])
            axes[j,col].set_ylabel('Science fiber ID (ordered rows)')
            axes[j,col].set_xlabel('Air wavelength (Å)')
        fig.colorbar(image,ax=axes[j,:],label=UNIT,extend='both',fraction=.025,pad=.02)
        spectrum=fig.add_subplot(grid[2*j+1,:])
        for label,color in [('DRP','#777777'),('Final SKYML','#0078b7')]:
            spectrum.plot(wave[band],median['stages'][label][band],color=color,lw=.7,label=label)
        spectrum.axhline(0,color='black',lw=.5);spectrum.set_xlim(lo,hi)
        spectrum.set_ylabel('Median FLUX\n(10⁻¹⁴ units)');spectrum.set_xlabel('Air wavelength (Å)')
        spectrum.set_title(f'Same {median["n_fibers"]} Gaia-selected Science fibers',fontsize=10)
        for sign in [-1,1]:spectrum.plot(wave[band],sign*median['sigma'][band],color='black',lw=.5)
        for y,label,color in [(.78,'DRP','#777777'),(.48,'SKYML','#0078b7'),(.18,'±1σ','black')]:
            spectrum.text(1.005,y,label,color=color,transform=spectrum.transAxes,fontsize=8,va='center')
        spectrum.grid(alpha=.15)
    fig.suptitle(f'Exposure {expnum} · all {len(rows)} science fibers · '+(window_group.replace('_',' ')+' diagnostics' if window_group else 'full channels'),fontsize=16)
    filename=f'map_{expnum}'+('_'+window_group+'_windows' if window_group else '')+'.png'
    fig.savefig(ROOT/filename,dpi=190);plt.close(fig)
    return filename,limits,len(rows)


def quality(h,drp,expnum):
    rows,good,arrays=flux_arrays(h,drp)
    wave=np.asarray(h['WAVE'].data,float)
    good &= h['SLITMAP'].data['skyml_applied'][rows,None]
    good &= ~h['MED_FIT_MASK'].data.astype(bool)[None,:]
    hp=[highpass(a,good) for a in arrays]
    raw=highpass(np.asarray(h['SKY_ML'].data[rows],float)*1e14,good)
    ivar=np.asarray(h['IVAR'].data[rows],float)/1e28
    result=[]
    for ch,lo,hi in CHANNELS:
        band=(wave>=lo)&(wave<hi)
        threshold=np.nanpercentile(np.maximum(raw[:,band],0),99,axis=1)*.04
        pixels=good&band[None,:]&(raw>threshold[:,None])
        count=pixels.sum(axis=1)
        scores=[np.sqrt(np.divide(np.sum(np.where(pixels,a*a*ivar,0),axis=1),count,out=np.full(len(rows),np.nan),where=count>=10)) for a in hp]
        valid=np.isfinite(scores[0])&np.isfinite(scores[1])
        result.append(dict(expnum=expnum,channel=ch,n_fibers=int(valid.sum()),drp=float(np.median(scores[0][valid])),skyml=float(np.median(scores[1][valid])),fraction_lower=float(np.mean(scores[1][valid]<scores[0][valid]))))
    return result


def iteration_image(h,input_row,survey,qa,expnum):
    fiberid=int(h['SLITMAP'].data['fiberid'][input_row])
    wave=np.asarray(h['WAVE'].data,float)
    history=h['FIT_LSF_HISTORY'].data
    fig,axes=plt.subplots(2,3,figsize=(16,8),layout='constrained')
    delta=np.linspace(-2.5,2.5,201)
    for j,(ch,lo,hi) in enumerate(CHANNELS[1:]):
        records=history[(history['input_row']==input_row)&(history['channel']==ch)]
        if not len(records):continue
        available=(wave>=lo)&(wave<hi)&~h['FIT_MASK'].data[input_row].astype(bool)
        line=np.where(available,h['FIT_LINES'].data[input_row],-np.inf)
        centre=float(wave[np.argmax(line)])
        initial=evaluate_lsf_density(survey,np.array([centre]),delta)[0]
        row=qa[qa['channel']==ch][0]
        axes[j,0].plot(delta,initial,color='black',ls='--',lw=2,label='Start: survey' if row['accepted'] else 'Adopted: survey start')
        widths=[]
        norm=Normalize(1,max(2,len(records)));cmap=plt.get_cmap('viridis')
        for record in records:
            state=deepcopy(survey)
            state.coefficients[ch]=np.asarray(record['coefficients']).reshape(survey.coefficients[ch].shape)
            profile=evaluate_lsf_density(state,np.array([centre]),delta)[0]
            color=cmap(norm(record['iteration']))
            axes[j,0].plot(delta,profile,color=color,lw=1.2)
            axes[j,1].plot(delta,profile-initial,color=color,lw=1.2)
            widths.append(record['w50'])
        last=records[-1]
        state.coefficients[ch]=np.asarray(last['coefficients']).reshape(survey.coefficients[ch].shape)
        final=evaluate_lsf_density(state,np.array([centre]),delta)[0]
        axes[j,0].plot(delta,final,color='#c02f63',lw=2,label=f'Last trial: iteration {int(last["iteration"])}')
        axes[j,0].set_title(f'{ch}: LSF at {centre:.1f} Å · '+('adopted in situ' if row['accepted'] else 'survey adopted'))
        axes[j,0].set_ylabel('Normalized LSF (Å⁻¹)');axes[j,0].legend(fontsize=9)
        axes[j,1].axhline(0,color='black',lw=.7)
        axes[j,1].set_title('Each iteration minus survey start')
        axes[j,1].set_ylabel('LSF density change (Å⁻¹)')
        for col in [0,1]:axes[j,col].set_xlabel('Wavelength offset (Å)');axes[j,col].grid(alpha=.2)
        initial_width=evaluate_lsf_diagnostics(survey,last['probe_wave'])['w50_angstrom']
        change=(np.asarray(widths)/initial_width-1)*100
        scale=max(float(np.max(np.abs(change))),.1)
        image=axes[j,2].imshow(change,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-scale,vmax=scale,
            extent=(last['probe_wave'][0],last['probe_wave'][-1],.5,len(records)+.5),interpolation='nearest')
        axes[j,2].set_title('Central-50% width change from start')
        axes[j,2].set_xlabel('Air wavelength (Å)');axes[j,2].set_ylabel('Fit iteration')
        axes[j,2].set_yticks(np.arange(1,len(records)+1))
        fig.colorbar(image,ax=axes[j,2],label='Width change (%)',fraction=.045)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes[j,:2],label='Fit iteration',ticks=np.arange(1,len(records)+1),fraction=.02)
    fig.suptitle(f'Exposure {expnum} · fiber {fiberid} · what changes during LSF recovery',fontsize=15)
    name=f'iterations_{expnum}_{fiberid}.png';fig.savefig(ROOT/name,dpi=160);plt.close(fig)
    return name


def overview(exposures,cases,page):
    body='<h1>Sky subtraction: whole-exposure comparison</h1><p class="notice">Left: DRP FLUX. Right: SKYML FLUX. All science fibers; shared colorbar within each pair. Click to enlarge.</p>'
    body+='<p>Units: 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹. Gray: masked. Pair limits: joint 98th percentile of |FLUX|; colors saturate beyond the limits.</p>'
    summary=fallback_summaries(exposures)
    (ROOT/'lsf_fallback_summary.json').write_text(json.dumps(summary,indent=2))
    body+=fallback_table(summary)
    body+=lsf_selection_diagram()
    scores=[];updated_cases=[]
    for e in exposures:
        with diagnostic_hdul(e['output']) as h,fits.open(e['output'].replace('lvmSFrame-skyml-','lvmSFrame-'),memmap=True) as drp:
            current_cases=diagnostic_cases(h,e['expnum']);updated_cases.extend(current_cases)
            qa=h['LSF_QA'].data
            adoption={ch:(int(np.sum(qa['accepted'][qa['channel']==ch])),int(np.sum(qa['channel']==ch))) for ch in 'RZ'}
            from plot_insitu_fit import write_fiber_diagnostic
            ids={c['fiberid'] for c in current_cases}
            ids.update(int(p.stem.split('_')[-1]) for p in ROOT.glob(f"insitu_fit_{e['expnum']}_*.html"))
            for fid in sorted(ids):write_fiber_diagnostic(h,drp,e['output'],fid)
            median=science_median_data(h,drp)
            full,_,n=flux_map(h,drp,e['expnum'],median)
            windows=[flux_map(h,drp,e['expnum'],median,group)[0] for group in WINDOW_GROUPS]
            median_plots=science_median_plots(median,e['expnum'])
            current=quality(h,drp,e['expnum']);scores.extend(current)
        exp=e['expnum']
        body+=f'<h2 id="exp{exp}">Exposure {exp} · {e["purpose"]}</h2><a href="{full}"><img src="{full}" width="3040" height="2850" loading="eager" alt="All {n} science fibers: direct DRP and SKYML FLUX, B R Z"></a>'
        body+=median_plots
        body+='<p>RMS of high-pass FLUX / input σ (CFrame IVAR⁻¹ᐟ²); source structure included. Same ≥10 valid line pixels per fiber; <a href="#metric">pixel selection</a>.</p><table><tr><th>Channel</th><th>DRP line-scale RMS / noise</th><th>SKYML line-scale RMS / noise ↓ better</th><th>Fibers with lower SKYML RMS</th></tr>'
        for s in current:body+=f'<tr><td>{s["channel"]}</td><td>{s["drp"]:.2f}</td><td>{s["skyml"]:.2f}</td><td>{s["fraction_lower"]:.0%} of {s["n_fibers"]}</td></tr>'
        body+='</table>'
        for zoom in windows:body+=f'<a href="{zoom}"><img src="{zoom}" width="3040" height="2850" loading="lazy" alt="Matched Science diagnostic windows"></a>'
        body+=case_table(current_cases,exp,median["wave"])
        fractions='; '.join(f'{ch}: {n/total:.1%} ({n}/{total})' for ch,(n,total) in adoption.items())
        body+=f'<p>All processed Science fibers — in-situ LSF: {fractions}; {e["seconds"]:.0f} s.</p>'
    body+='<h3 id="metric">Metric</h3><p>Median per-fiber RMS of high-pass FLUX / observed σ. High-pass: masked 51-pixel mean removed. Line pixels: >4% of each fiber/channel’s 99th-percentile positive high-pass ML sky; protected science lines masked. This includes source structure and sky residuals; it is not calibrated χ². σ uses observed IVAR; full sky uncertainty remains uncalibrated.</p>'
    body+='<p><a href="MEDIAN_FIT_AUDIT.md">Median-fit investigation</a> · <a href="median_fit_audit.json">Fit measurements</a> · <a href="../../PLOTTING_AND_REPORT_RULES.md">Plotting rules</a> · <a href="quality_v2.json">Quality measurements</a> · <a href="run_manifest_v2.json">Five output FITS</a> · <a href="FITS_STRUCTURE.md">FITS structure and sizes</a> · <a href="pipeline_README.md">Method and CLI</a> · <a href="integrity_v2.json">FITS checks</a> · <a href="reviews_diagnostic_clarity.md">Independent presentation reviews</a></p>'
    (ROOT/'fiber_cases.json').write_text(json.dumps(updated_cases,indent=2,allow_nan=False))
    (ROOT/'quality_v2.json').write_text(json.dumps(scores,indent=2))
    (ROOT/'00_MAIN_sky_subtraction_report.html').write_text(page('Whole-exposure sky subtraction quality',body))
