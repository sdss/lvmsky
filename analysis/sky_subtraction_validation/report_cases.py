"""Representative reproducible per-exposure examples with measured selection evidence."""
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter,minimum_filter1d
from report_contrast import contrast_rows


def diagnostic_cases(h,expnum):
    rows=np.asarray(h['LSF_FIBERS'].data['input_row'],int)
    wave=np.asarray(h['WAVE'].data,float);phot=h['GAIA'].data;qa=h['LSF_QA'].data
    lookup={(int(q['input_row']),str(q['channel'])):q for q in qa}
    rms={ch:np.full(len(rows),np.nan) for ch in 'RZ'}
    spike=np.full(len(rows),np.nan);spike_wave=spike.copy()
    blue=(wave<5787)
    for start in range(0,len(rows),64):
        ids=rows[start:start+64]
        obs=np.asarray(h['OBSFLUX'].data[ids],float);model=np.asarray(h['FIT_MODEL'].data[ids],float)
        ivar=np.asarray(h['IVAR'].data[ids],float)
        good=~h['FIT_MASK'].data[ids].astype(bool)&np.isfinite(obs)&np.isfinite(model)&(ivar>0)
        for ch,lo,hi in [('R',5787,7454),('Z',7454,9801)]:
            used=good&(wave>=lo)&(wave<hi)
            count=used.sum(axis=1)
            rms[ch][start:start+len(ids)]=np.sqrt(np.divide(np.sum(np.where(used,(obs-model)**2*ivar,0),axis=1),count,out=np.full(len(ids),np.nan),where=count>0))
        valid=(h['MASK'].data[ids][:,blue]==0)&(ivar[:,blue]>0)&np.isfinite(obs[:,blue])
        neighbors=minimum_filter1d(valid.astype('u1'),size=9,axis=1,mode='constant',cval=0).astype(bool)
        sharp=(obs[:,blue]-median_filter(np.where(valid,obs[:,blue],0),size=(1,9),mode='nearest'))*np.sqrt(np.maximum(ivar[:,blue],0))
        available=neighbors&~h['MED_FIT_MASK'].data[blue].astype(bool)[None,:]
        sharp=np.where(available,sharp,-np.inf)
        peak=np.argmax(sharp,axis=1)
        spike[start:start+len(ids)]=sharp[np.arange(len(ids)),peak]
        spike_wave[start:start+len(ids)]=wave[blue][peak]
    ratio=np.asarray(phot['ratio_gaia_to_lvm_flux_plus_sky_g'][rows],float)
    stellar=np.asarray(phot['gaia_g_flux'][rows],float)
    accepted=np.array([all(bool(lookup[(int(row),ch)]['accepted']) for ch in 'RZ') for row in rows])
    score=np.maximum(rms['R'],rms['Z']);control=accepted&np.isfinite(ratio)&(ratio>=0)&(ratio<.01)&np.isfinite(score)
    selected=[];taken=set()
    def choose(order,label):
        for i in order:
            if int(i) not in taken:
                selected.append((int(i),label));taken.add(int(i));return
        raise ValueError('Distinct diagnostic Science fibers are required')
    pool=np.flatnonzero(control)
    if len(pool):choose(pool[np.argsort(score[pool],kind='stable')],'Low-residual control')
    else:choose(np.argsort(np.where(np.isfinite(score),score,np.inf),kind='stable'),'Lowest available fit residual')
    choose(np.argsort(np.where(np.isfinite(stellar),-stellar,np.inf),kind='stable'),'Highest summed Gaia stellar flux')
    choose(np.argsort(np.where(np.isfinite(spike),-spike,np.inf),kind='stable'),'Sharp-input-feature stress test')
    choose(np.argsort(np.where(np.isfinite(score),-score,np.inf),kind='stable'),'Large-fit-residual stress test')
    fallback=np.flatnonzero(~accepted)
    if len(fallback):choose(fallback[np.argsort(np.where(np.isfinite(score[fallback]),-score[fallback],np.inf),kind='stable')],'Survey-fallback example')
    else:choose(np.argsort(np.where(np.isfinite(ratio),-ratio,np.inf),kind='stable'),'High Gaia contamination example')
    # Show several adopted fits even when only Z has sufficient line contrast.
    n_accepted=np.array([sum(bool(lookup[(int(row),ch)]['accepted']) for ch in 'RZ') for row in rows])
    pool=np.flatnonzero((n_accepted>0)&np.isfinite(score))
    order=pool[np.lexsort((score[pool],-n_accepted[pool]))]
    for i in [i for i in order if int(i) not in taken][:3]:
        choose([i],'Additional adopted in-situ example')
    def finite(value):return float(value) if np.isfinite(value) else None
    output=[]
    for i,label in selected:
        row=int(rows[i]);fiberid=int(h['SLITMAP'].data['fiberid'][row])
        contrast=contrast_rows(wave,np.asarray(h['FIT_LINES'].data[row],float),np.asarray(h['FIT_CONT'].data[row],float),np.asarray(h['FIT_SOLAR'].data[row],float),~h['FIT_MASK'].data[row].astype(bool),np.asarray(h['OBSFLUX'].data[row],float))
        output.append(dict(expnum=expnum,fiberid=fiberid,input_row=row,selection=label,ratio=finite(ratio[i]),gaia_g_flux=finite(stellar[i]),
            fit_rms={ch:finite(rms[ch][i]) for ch in 'RZ'},sharpness=finite(spike[i]),sharp_wave=finite(spike_wave[i]),
            source_contrast={c['channel']:finite(c['source_ratio']) for c in contrast},observed_contrast={c['channel']:finite(c['observed_ratio']) for c in contrast},
            qa={ch:dict(accepted=bool(lookup[(row,ch)]['accepted']),reason=str(lookup[(row,ch)]['reason']),counts=np.asarray(lookup[(row,ch)]['quality_bin_counts'],int).tolist()) for ch in 'RZ'},
            page=f'insitu_fit_{expnum}_{fiberid}.html'))
    return output


def case_table(cases,expnum,wave):
    def value(x):return '—' if x is None else f'{x:.2g}'
    body='<h3>Diagnostic fibers</h3><div class="notice"><b>Sky-line coverage:</b> Q1–Q4 are consecutive wavelength regions, blue → red. Each needs ≥2 usable sky features; each feature needs contrast ≥1 and line-amplitude S/N ≥5.<br>'
    for channel,lo,hi in [('R',5787,7454),('Z',7454,9801)]:
        native=wave[(wave>=lo)&(wave<hi)]
        edges=np.linspace(native[0],native[-1]+np.median(np.diff(native)),5)
        regions='; '.join(f'Q{i+1}: {a:.1f}–{b:.1f}' for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])))
        body+=f'<div><b>{channel} (Å):</b> {regions}.</div>'
    body+='In-situ adoption requires all four regions, convergence and profile checks to pass (<a href="#lsf-selection">full LSF rule</a>). B uses survey. Fit RMS and Gaia ratio are diagnostics, not acceptance thresholds.</div><table><caption style="text-align:left;padding:10px 8px"><b>How to read:</b> numbers are usable sky-line counts in Q1–Q4. Each region needs ≥2 lines. For example, [1, 0, 2, 6] passes 2/4 regions → survey LSF. All 4/4 must pass to allow in-situ LSF; convergence and profile checks also apply.</caption><tr><th>Fiber</th><th>Selection / evidence</th><th>Gaia G / observed G</th><th>Fit RMS / σ</th><th>Usable sky-feature counts · ≥2 per region</th><th>Adopted LSF</th></tr>'
    for c in [c for c in cases if c['expnum']==expnum]:
        label=c['selection']
        if label.startswith('Sharp'):label+=f': {value(c["sharpness"])} × input σ at {c["sharp_wave"]:.1f} Å'
        adopted=[];counts=[];fit_rms=[]
        for channel in 'RZ':
            q=c['qa'][channel]
            reason='insufficient line coverage' if q['reason']=='line_information' else q['reason'].replace('_',' ')
            decision='in situ' if q['accepted'] else 'survey — '+reason
            adopted.append(f'<div><b>{channel}:</b> {decision}</div>')
            values=', '.join(f'Q{i+1}: {n}' for i,n in enumerate(q['counts']))
            covered=sum(n>=2 for n in q['counts'])
            counts.append(f'<div><b>{channel}:</b> {values}<br>{covered}/4 regions have ≥2 lines → {"in situ" if q["accepted"] else "survey"}</div>')
            fit_rms.append(f'<div><b>{channel}:</b> {value(c["fit_rms"][channel])}</div>')
        body+=f'<tr><td><a href="{c["page"]}">{c["fiberid"]}</a></td><td>{label}</td><td>{value(c["ratio"])}</td><td>{"".join(fit_rms)}</td><td>{"".join(counts)}</td><td>{"".join(adopted)}</td></tr>'
    return body+'</table><p>σ = observed CFrame IVAR⁻¹ᐟ². Fit RMS/σ = √mean[(data−fit)²/σ²]. Sharpness is a 9-pixel positive deviation / σ. These auxiliary fits do not independently establish LSF accuracy.</p>'


if __name__=='__main__':
    from astropy.io import fits
    from skysub.sky_subtraction.reconstruction import diagnostic_hdul
    root=Path(__file__).resolve().parent;cases=[]
    for e in json.loads((root/'run_manifest_v2.json').read_text()):
        with diagnostic_hdul(e['output']) as h:cases.extend(diagnostic_cases(h,e['expnum']))
    assert len(cases)==40 and len({(c['expnum'],c['fiberid']) for c in cases})==40
    assert all(sum(c['expnum']==e['expnum'] for c in cases)==8 for e in json.loads((root/'run_manifest_v2.json').read_text()))
    (root/'fiber_cases.json').write_text(json.dumps(cases,indent=2,allow_nan=False))
    print(f'Selected {len(cases)} diagnostic fibers.')
