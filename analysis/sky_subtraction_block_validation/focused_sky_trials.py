"""Test spatial coupling on two Sky examples while retaining centering and frozen populations."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from copy import deepcopy
import json,sys
import numpy as np
import pandas as pd
from astropy.time import Time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from sky_diagnostics import HERE,CHANNELS,residual_metrics,example_table,assemble_diagnostics
from report import spectrum,page
from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
from skysub.sky_subtraction_block.lsf import channel_tensors,solve_normal,offset_centers
from skysub.sky_subtraction.pipeline import render,FACTOR


def main():
    rows=json.loads((HERE/'run_manifest.json').read_text());row=next(r for r in rows if r['expnum']==16819)
    data=pd.read_csv(HERE/'sky_metrics_16819.csv');block=data[data.variant=='block']
    scores=block.groupby('fiberid').chi2_per_pixel.mean()
    worst=int(scores.idxmax());typical=int((scores-scores.median()).abs().idxmin())
    selected=[(worst,'Largest mean B/R/Z χ²/pixel'),(typical,'Typical comparison: nearest the median mean B/R/Z χ²/pixel')]
    records=[];front=example_table(data,selected);trial_body='';saved=[]
    with SkyBlockReconstructor(row['output'],cframe=row['input']) as r:
        mjd=float(Time(r.h[0].header['OBSTIME'],scale='utc').mjd)
        for f,reason in selected:
            rowfit=r.h['SKY_FIBER_COEF'].data[r.h['SKY_FIBER_COEF'].data['fiberid']==f][0]
            label=str(rowfit['telescope']);sp=int(rowfit['spectrograph']);i=int(rowfit['input_row'])
            d=r.sky_fiber(f);state=r.global_state(f);candidate=deepcopy(state)
            designs=channel_tensors(r.model(label),rowfit['coef'],r.state('MED_',r.labels.index(label)))
            for ch,(take,design,fixed) in designs.items():
                wave=r.wave[take];good=d['valid'][take]&(wave>=wave[0]+5)&(wave<=wave[-1]-5)
                target=d['observed'][take]*FACTOR-fixed
                weight=r.ivar[i,take]/FACTOR**2
                a=design[good]*np.sqrt(weight[good,None]);y=target[good]*np.sqrt(weight[good])
                no,nw=state.coefficients[ch].shape
                coef,status=solve_normal(a.T@a,-a.T@y,(no,nw,1),offset_centers(no,state.config['offset_half_width_angstrom']))
                candidate.coefficients[ch]=coef[:,:,0]
            trial,components,_=render(r.model(label),rowfit['coef'],candidate)
            trial/=FACTOR
            for ch,(lo,hi) in CHANNELS.items():
                if ch=='B':continue
                take=(r.wave>=lo)&(r.wave<hi);wave=r.wave[take]
                good=d['valid']&take&(r.wave>=wave[0]+5)&(r.wave<=wave[-1]-5)
                original=residual_metrics(d['observed'],d['model'],d['sigma'],good)
                changed=residual_metrics(d['observed'],trial,d['sigma'],good)
                assert changed['chi2_per_pixel']<=original['chi2_per_pixel']+1e-5
                records.append(dict(expnum=16819,fiberid=f,telescope=label,spectrograph=sp,channel=ch,
                    block_chi2=original['chi2_per_pixel'],individual_centered_chi2=changed['chi2_per_pixel']))
            title=f'Exposure 16819 · {label} · SP{sp} · fiber {f} · MJD {mjd:.5f}'
            front+='<p>'+reason+'. Actual delivered block fit.</p>'+spectrum(title,r.wave,d['observed'],d['model'],d['components'],d['sigma'],d['valid'],fragment=True)
            trial_body+='<p>Diagnostic trial: a separate centered LSF per fiber/channel with frozen component coefficients, the same wavelength controls, the same noise pixels and no spatial coupling. This is a flexible in-sample comparison; the delivered SFrame remains unchanged.</p>'+spectrum(title+' · independent centered LSF trial',r.wave,d['observed'],trial,{k:v/FACTOR for k,v in components.items()},d['sigma'],d['valid'],fragment=True)
            saved.append(dict(fiberid=f,coefficients={ch:candidate.coefficients[ch].tolist() for ch in 'BRZ'}))
    (HERE/'sky_front_examples.html').write_text(front)
    (HERE/'sky_fibers_16819_independent_trial.html').write_text(page('Exposure 16819: independent centered LSF trials',trial_body+'<p><a href="00_MAIN_sky_block_report.html">Main report</a></p>'))
    pd.DataFrame(records).to_csv(HERE/'sky_independent_lsf_trials.csv',index=False)
    (HERE/'sky_independent_lsf_coefficients.json').write_text(json.dumps(saved,indent=2))
    assemble_diagnostics(rows)
    print(pd.DataFrame(records).to_string(index=False))


if __name__=='__main__':main()
