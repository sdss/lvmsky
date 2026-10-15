"""Native-noise diagnostics and controlled coefficient/shift trials for selected Sky fibers."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import json,re,sys,time
import numpy as np
import pandas as pd
from astropy.time import Time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from report import HERE,spectrum,page,assemble
from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
from skysub.sky_subtraction_block.lsf import FIT_CHANNELS,individual_coefficients,channel_tensors,fit_block
from skysub.sky_subtraction.pipeline import render,FACTOR
CHANNELS={'B':(3600.,5787.),'R':(5787.,7454.),'Z':(7454.,9801.)}


def residual_metrics(observed,model,sigma,good):
    z=(observed[good]-model[good])/sigma[good]
    if not len(z):raise ValueError('No valid pixels for residual metrics')
    return dict(n_pixels=len(z),chi2_per_pixel=float(np.mean(z*z)),median_z=float(np.median(z)),
        within_1sigma=float(np.mean(abs(z)<=1)),median_abs_z=float(np.median(abs(z))))


def shift_trial(wave,observed,model,sigma,good,limit=.25):
    """First-order model(lambda-delta), with fixed coefficients and profile shape."""
    derivative=np.gradient(model,wave)
    x=derivative[good]/sigma[good];y=(observed[good]-model[good])/sigma[good]
    denominator=float(x@x)
    delta=float(np.clip(-(x@y)/denominator,-limit,limit)) if denominator>0 else 0.
    return model-delta*derivative,delta


def logged_timings():
    records=[]
    for p in sorted(HERE.glob('run_*.log')):
        start=end=None
        for line in p.read_text().splitlines():
            match=re.match(r'(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d+) INFO (.*)',line)
            if not match:continue
            date=datetime.strptime(match[1],'%Y-%m-%d %H:%M:%S,%f');message=match[2]
            if message.startswith('Decomposing ') and 'selected Sky fibers' in message:start=date
            if message.startswith('Rendering '):end=date
        if start and end:records.append(dict(expnum=int(p.stem.split('_')[1]),logged_stage_seconds=(end-start).total_seconds()))
    assert len(records)==8
    return records


def analyze(row):
    e=row['expnum'];print('Sky diagnostics',e,flush=True)
    records=[];examples={};scores={};trials=[];benchmark_samples=defaultdict(list);tensor_seconds=0.
    with SkyBlockReconstructor(row['output'],cframe=row['input']) as r:
        mjd=float(Time(r.h[0].header['OBSTIME'],scale='utc').mjd)
        rows=r.h['SKY_FIBER_COEF'].data
        for n,rowfit in enumerate(rows):
            f=int(rowfit['fiberid']);i=int(rowfit['input_row']);label=str(rowfit['telescope']);sp=int(rowfit['spectrograph'])
            model=r.model(label);median_state=r.state('MED_',r.labels.index(label));state=r.global_state(f)
            coef=np.asarray(rowfit['coef'],float);observed=r.observed[i]
            sigma=np.divide(1,np.sqrt(np.maximum(r.ivar[i],0)),out=np.full(len(r.wave),np.nan),where=r.ivar[i]>0)
            mask=r.mask[i].copy();mask[r.h['MED_FIT_MASK'].data.astype(bool)]|=1
            valid=(mask==0)&np.isfinite(observed)&np.isfinite(sigma)&(sigma>0)
            initial,_,_=render(model,coef,median_state)
            current,components,_=render(model,coef,state)
            refit_coef,status=individual_coefficients(model,state,observed*FACTOR,r.ivar[i]/FACTOR**2,mask)
            refit,_,_=render(model,refit_coef,state)
            initial=initial/FACTOR;current=current/FACTOR;refit=refit/FACTOR
            trials.append(dict(fiberid=f,coef=refit_coef.tolist(),status=status))
            chi_total=[]
            for ch,(lo,hi) in CHANNELS.items():
                take=(r.wave>=lo)&(r.wave<hi);x=r.wave[take]
                good=valid&take&(r.wave>=x[0]+5)&(r.wave<=x[-1]-5)
                shifted,delta=shift_trial(r.wave,observed,current,sigma,good)
                for variant,values in [('median_lsf',initial),('block',current),('coef_refit',refit),('linear_shift',shifted)]:
                    info=residual_metrics(observed,values,sigma,good)
                    records.append(dict(expnum=e,fiberid=f,telescope=label,spectrograph=sp,channel=ch,variant=variant,
                        diagnostic_shift_A=delta if variant=='linear_shift' else 0.,**info))
                    if variant=='block':chi_total.append(info['chi2_per_pixel'])
            score=float(np.mean(chi_total));scores[f]=dict(score=score,telescope=label,spectrograph=sp)
            examples[f]=dict(observed=observed.copy(),model=current,components={k:v/FACTOR for k,v in components.items()},sigma=sigma,valid=valid)
            if e==59750:
                start=time.perf_counter();designs=channel_tensors(model,coef,median_state);tensor_seconds+=time.perf_counter()-start
                benchmark_samples[sp].append(dict(position=float(f),wave=r.wave,flux=observed*FACTOR,ivar=r.ivar[i]/FACTOR**2,mask=mask,designs=designs))
            if (n+1)%25==0:print(e,n+1,'Sky fibers analyzed',flush=True)
        frame=pd.DataFrame(records);frame.to_csv(HERE/f'sky_metrics_{e}.csv',index=False)
        (HERE/f'sky_refit_coefficients_{e}.json').write_text(json.dumps(dict(exposure=e,description='Diagnostic coefficient QP at frozen block LSF; products unchanged',fibers=trials)))
        selected=[]
        for label in ('SkyE','SkyW'):
            for sp in [1,2,3]:
                candidates=[f for f,s in scores.items() if s['telescope']==label and s['spectrograph']==sp]
                if candidates:selected.append((max(candidates,key=lambda f:scores[f]['score']),'Largest mean B/R/Z χ²/pixel in this telescope/spectrograph'))
        control=min(scores,key=lambda f:scores[f]['score'])
        if control not in [f for f,_ in selected]:selected.append((control,'Lowest mean B/R/Z χ²/pixel in this exposure (optimistically selected low-residual example)'))
        intro=f'<h1>Exposure {e}: individual Sky-fiber fits</h1><p>Actual block LSF with the frozen individual component coefficients used in the block solve. σ is the native CFrame input error. Each figure is an independent full-range spectrum with flux residual at −5 RMS and an aligned residual/σ panel. Metrics trim 5 Å at channel edges and use the same pixels for every diagnostic trial. χ²/pixel is the mean squared normalized residual, not reduced χ². σ describes input noise without fitted-model covariance.</p>'
        body=intro;body+=example_table(frame,selected)
        for f,reason in selected:
            info=scores[f];d=examples[f]
            title=f'Exposure {e} · MJD {mjd:.5f} · {info["telescope"]} · SP{info["spectrograph"]} · fiber {f}'
            body+='<p>'+reason+'. Pseudo-slit coordinate: physical fiber ID '+str(f)+'.</p>'+spectrum(title,r.wave,d['observed'],d['model'],d['components'],d['sigma'],d['valid'],fragment=True)
        (HERE/f'sky_fibers_{e}.html').write_text(page(f'Exposure {e}: Sky fiber fits',body+'<p><a href="00_MAIN_sky_block_report.html">Main report</a></p>'))
        # The front-page examples show the challenged exposure and a within-exposure control.
        if e==16819:
            front=[]
            for f,reason in [(max(scores,key=lambda f:scores[f]['score']),'Largest mean B/R/Z χ²/pixel'),(control,'Lowest mean B/R/Z χ²/pixel control')]:
                d=examples[f];info=scores[f]
                title=f'Exposure {e} · {info["telescope"]} · SP{info["spectrograph"]} · fiber {f} · MJD {mjd:.5f}'
                front.append('<p>'+reason+'.</p>'+spectrum(title,r.wave,d['observed'],d['model'],d['components'],d['sigma'],d['valid'],fragment=True))
            (HERE/'sky_front_examples.html').write_text(''.join(front))
        if e==59750:
            benchmark=[]
            template=r.state('MED_',r.labels.index('Sci'));slit=r.h['SLITMAP'].data
            for sp,samples in benchmark_samples.items():
                domain=slit['fiberid'][slit['spectrographid']==sp]
                for ch in FIT_CHANNELS:
                    coefficient,_,qa=fit_block(samples,template,ch,float(domain.min()),float(domain.max()))
                    np.testing.assert_allclose(coefficient,r.blocks[(sp,ch)]['coef'],atol=1e-7,rtol=1e-7)
                    benchmark.append(dict(spectrograph=sp,channel=ch,normal_seconds=qa['normal_seconds'],solve_seconds=qa['solve_seconds']))
            (HERE/'block_timing_replay.json').write_text(json.dumps(dict(expnum=e,tensor_seconds=tensor_seconds,blocks=benchmark,description='Replay with saved individual coefficients, one native thread; concurrent exposure diagnostic processes'),indent=2))
    return e


def example_table(frame,selected):
    body='<table><tr><th>Fiber / selection</th><th>B χ²/pixel</th><th>R χ²/pixel</th><th>Z χ²/pixel</th><th>Pixels within ±1σ (B / R / Z)</th></tr>'
    for f,reason in selected:
        group=frame[(frame.fiberid==f)&(frame.variant=='block')].set_index('channel')
        identity=frame[(frame.fiberid==f)&(frame.variant=='block')].iloc[0]
        body+=f'<tr><td>{f} · {identity.telescope} · SP{int(identity.spectrograph)} · {reason}</td>'+''.join(f'<td>{group.loc[ch,"chi2_per_pixel"]:.2f}</td>' for ch in 'BRZ')+'<td>'+' / '.join(f'{group.loc[ch,"within_1sigma"]:.0%}' for ch in 'BRZ')+'</td></tr>'
    return body+'</table><p>B: 3600–5787 Å; R: 5787–7454 Å; Z: 7454–9801 Å. For independent Gaussian noise and an adequate fixed model, approximately 68% of residuals fall within ±1σ; coefficient estimation and shared calibration introduce correlations.</p>'


def assemble_diagnostics(rows):
    data=pd.concat([pd.read_csv(HERE/f'sky_metrics_{row["expnum"]}.csv') for row in rows],ignore_index=True)
    data.to_csv(HERE/'sky_fiber_residual_metrics.csv',index=False)
    summary=data.groupby(['expnum','channel','variant']).agg(fibers=('fiberid','nunique'),median_chi2=('chi2_per_pixel','median'),median_within_1sigma=('within_1sigma','median'),median_abs_shift_A=('diagnostic_shift_A',lambda x:float(np.median(abs(x))))).reset_index()
    summary.to_csv(HERE/'sky_fiber_trial_summary.csv',index=False)
    timings=logged_timings();(HERE/'logged_block_timings.json').write_text(json.dumps(timings,indent=2))
    replay=json.loads((HERE/'block_timing_replay.json').read_text())
    body='<div id="sky-fibers"><h2>Individual Sky fibers: fit quality and improvement trials</h2><p>These are actual Sky-fiber fits, independent of the remaining source emission in Science. Full-range diagnostics show the observed data, seven additive components, total model, flux residual at −5 RMS with thin input-error curves, and residual/σ with visible ±1σ references.</p>'
    body+='<p>'+' · '.join(f'<a href="sky_fibers_{row["expnum"]}.html">{row["expnum"]}: SkyE/SkyW, SP1–3</a>' for row in rows)+'</p>'
    body+='<h3>Actual Sky fits: challenged and typical examples</h3>'+(HERE/'sky_front_examples.html').read_text()
    body+='<h3 id="block-runtime">Runtime</h3><p>The logged pre-render coefficient → tensor → nine-QP stage takes '+f'{min(t["logged_stage_seconds"] for t in timings):.1f}–{max(t["logged_stage_seconds"] for t in timings):.1f} s per exposure on Lumen, with 8–12 worker processes. For 59750 it takes 24.0 s. A replay of its nine QPs from saved individual coefficients takes '+f'{sum(b["solve_seconds"] for b in replay["blocks"]):.3f} s for QP setup, solve and validation, '+f'{sum(b["normal_seconds"] for b in replay["blocks"]):.3f} s assembling normal equations, and {replay["tensor_seconds"]:.2f} s building spectral tensors. Replay uses one native thread with concurrent exposure diagnostics; it excludes median/ML preparation, Science rendering and FITS writing. <a href="block_timing_replay.json">Replay measurements</a> · <a href="logged_block_timings.json">Original stage timings</a>.</p>'
    body+='<h3>Controlled diagnostic trials</h3><p>Median LSF: the initial telescope-median LSF with the saved individual coefficients. Block: the delivered global LSF with those coefficients frozen. Coefficient refit: a new component QP at the frozen block LSF. Shift trial: a first-order model(λ−δ) approximation, one scalar δ per fiber/channel, bounded to ±0.25 Å, with the original frozen coefficients. Trials use identical native-noise pixels. The coefficient QP fits the full valid spectrum, while the table reports trimmed channels separately, so an individual channel statistic may increase. These conditional fits diagnose model flexibility; they are not held-out recovery tests and do not alter any SFrame.</p>'
    body+='<table><tr><th>Exposure / channel</th><th>Fibers</th><th>Median-LSF χ²/pixel</th><th>Block χ²/pixel</th><th>Coefficient-refit χ²/pixel</th><th>Shift-trial χ²/pixel</th><th>Median |δ| (Å)</th></tr>'
    for e,ch in sorted(set(zip(summary.expnum,summary.channel))):
        if ch=='B':continue
        group=summary[(summary.expnum==e)&(summary.channel==ch)].set_index('variant')
        body+=f'<tr><td>{e} / {ch}</td><td>{int(group.loc["block","fibers"])}</td>'+''.join(f'<td>{group.loc[v,"median_chi2"]:.2f}</td>' for v in ['median_lsf','block','coef_refit','linear_shift'])+f'<td>{group.loc["linear_shift","median_abs_shift_A"]:.4f}</td></tr>'
    body+='</table><p><a href="sky_fiber_residual_metrics.csv">Every selected fiber and channel: residuals versus input noise</a> · <a href="sky_fiber_trial_summary.csv">Trial summary</a></p>'
    body+='<h3>What the trials establish</h3><p>Exposure 16819 Z already has median χ²/pixel 11.66 with the initial median LSF; the delivered block gives 11.88, coefficient refitting 11.75, and the constant-shift trial 11.86. These two simple changes do not bring this exposure to input-noise residuals. In 12630 Z, the shift trial improves 1.55 to 1.33: the effect depends on the exposure.</p>'
    trial_path=HERE/'sky_independent_lsf_trials.csv'
    if trial_path.exists():
        trial=pd.read_csv(trial_path)
        body+='<p>A separate centered LSF for each of two 16819 example fibers tests the effect of spatial coupling while retaining the same frozen populations and wavelength basis. This is an in-sample flexibility trial. <a href="sky_fibers_16819_independent_trial.html">Full independent-LSF trial spectra</a>.</p><table><tr><th>Fiber / channel</th><th>Block χ²/pixel</th><th>Independent centered LSF χ²/pixel</th></tr>'
        for _,r in trial.iterrows():body+=f'<tr><td>{int(r.fiberid)} / {r.channel}</td><td>{r.block_chi2:.2f}</td><td>{r.individual_centered_chi2:.2f}</td></tr>'
        body+='</table><p>For these two fibers, removing spatial coupling alone changes Z χ²/pixel only from 23.43 to 23.37 and from 11.60 to 11.58. Flexible individual profiles with the same centering and frozen populations therefore do not resolve the large residuals in these examples.</p>'
    body+='<p><b>Next tests:</b> use the residual shapes to separate line-profile mismatch from continuum/template, transmission, flux-calibration and noise-model effects; test wavelength-dependent individual shifts while keeping the intrinsic LSF centered; then choose slit flexibility and smoothness using held-out Sky fibers. Component and LSF updates can alternate as separate linear QPs, but the combined coefficient × LSF problem is not one linear QP. The current constant-shift trial cannot diagnose wavelength-dependent offsets.</p>'
    body+='</div>'
    assert '<details' not in body and data.fiberid.notna().all()
    assert len(data)==sum(pd.read_csv(HERE/f'sky_metrics_{row["expnum"]}.csv').shape[0] for row in rows)
    (HERE/'sky_diagnostics_section.html').write_text(body)
    assemble(rows)
    print(summary.to_string(index=False))


if __name__=='__main__':
    rows=json.loads((HERE/'run_manifest.json').read_text())
    if '--assemble-only' not in sys.argv:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(4) as pool:list(pool.map(analyze,rows))
    assemble_diagnostics(rows)
