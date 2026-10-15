"""Measure saved median fits and test local line displacement without changing FITS."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))


def local_line_probe(wave,observed,lines,model,ivar,center):
    """Diagnostic only: shared template amplitude, baseline, slope and displacement."""
    used=(np.abs(wave-center)<4)&np.isfinite(observed)&np.isfinite(ivar)&(ivar>0)
    if used.sum()<8:return None
    x=wave[used]-center;y=observed[used];weight=np.sqrt(ivar[used])
    design=np.column_stack([lines[used],np.ones(len(x)),x])
    start=np.linalg.lstsq(design*weight[:,None],y*weight,rcond=None)[0]
    def residual(p):
        return (p[0]*np.interp(wave[used]-p[3],wave,lines)+p[1]+p[2]*x-y)*weight
    fit=least_squares(residual,[np.clip(start[0],1e-6,3-1e-6),start[1],start[2],0],
                      bounds=([0,-np.inf,-np.inf,-.5],[3,np.inf,np.inf,.5]))
    zero=least_squares(lambda p:residual([*p,0]),fit.x[:3],
                       bounds=([0,-np.inf,-np.inf],[3,np.inf,np.inf]))
    if not fit.success or not zero.success:return None
    return dict(wavelength=float(center),shift=float(fit.x[3]),amplitude=float(fit.x[0]),
                rms_saved=float(np.sqrt(np.mean((y-model[used])**2))),
                rms_local_zero=float(np.sqrt(np.mean((residual([*zero.x,0])/weight)**2))),
                rms_local=float(np.sqrt(np.mean((residual(fit.x)/weight)**2))))


def median_line_probes(wave,observed,lines,model,ivar):
    peaks,_=find_peaks(lines,distance=max(1,int(np.ceil(30/np.median(np.diff(wave))))))
    peaks=sorted((j for j in peaks if 7500<wave[j]<9700),key=lambda j:-lines[j])[:12]
    probes=[local_line_probe(wave,observed,lines,model,ivar,wave[j]) for j in peaks]
    return [p for p in probes if p is not None]


def audit_product(path):
    from skysub.sky_subtraction.reconstruction import SkyMLReconstructor
    from skysub.sky_subtraction.pipeline import render
    from validate import highpass
    output={}
    with SkyMLReconstructor(path) as r:
        for i,label in enumerate(r.labels):
            observed=r.h['MED_FLUX'].data[i]*1e14;model,components=r.medians[label];model=model*1e14
            ivar=r.h['MED_IVAR'].data[i]/1e28
            used=(ivar>0)&np.isfinite(observed)&np.isfinite(model)&~r.h['MED_FIT_MASK'].data.astype(bool)
            channels={}
            for ch,lo,hi in [('B',3600,5787),('R',5787,7454),('Z',7454,9801)]:
                take=used&(r.wave>=lo)&(r.wave<hi);residual=observed[take]-model[take]
                channels[ch]=dict(rms=float(np.sqrt(np.mean(residual**2))),
                    normalized_rms=float(np.sqrt(np.mean(residual**2*ivar[take]))),
                    median_sigma=float(np.median(1/np.sqrt(ivar[take]))))
            lines=sum(components[k] for k in ('oh','atom','orc','o2'))*1e14
            probes=median_line_probes(r.wave,observed,lines,model,np.where(used,ivar,0))
            output[label]=dict(channels=channels,line_probes=probes)
            if label=='Sci':
                fixed={}
                for name,coef in [('decomposition',r.h['MED_COEF'].data[i]),('ML',r.coef),('ML_OH',r.oh_coef)]:
                    sky,_,_=render(r.model_sci,coef,r.state('MED_',i))
                    residual=highpass((observed-sky)[None,:],used[None,:])[0]
                    fixed[name]={ch:float(np.sqrt(np.mean(residual[used&(r.wave>=lo)&(r.wave<hi)]**2)))
                                 for ch,lo,hi in [('B',3600,5787),('R',5787,7454),('Z',7454,9801)]}
                output[label]['fixed_median_lsf_rms']=fixed
                provenance=json.loads(bytes(r.h['SKYML_META'].data).decode())
                output[label]['oh_correction']=provenance['oh_correction']
    return output


def audit_summary(audit,expnum):
    probes=sorted(audit['Sci']['line_probes'],key=lambda p:-p['rms_saved'])[:3]
    body='<h3>Local line-position check</h3><p>Three largest saved-fit RMS among ≤12 strongest modeled Z peaks, ≥30 Å apart. Local ±4 Å refits: amplitude + baseline + slope, with shift fixed at 0 or free within ±0.5 Å. Positive shift → longer wavelength. Same input IVAR; RMS in 10⁻¹⁴ units. Diagnostic only.</p><table><tr><th>Line region (Å)</th><th>Shift (Å)</th><th>Saved-fit RMS</th><th>Zero-shift local RMS</th><th>Free-shift local RMS</th></tr>'
    for p in probes:
        body+=f'<tr><td>{p["wavelength"]:.1f}</td><td>{p["shift"]:+.3f}</td><td>{p["rms_saved"]:.3g}</td><td>{p["rms_local_zero"]:.3g}</td><td>{p["rms_local"]:.3g}</td></tr>'
    return body+'</table>'


if __name__=='__main__':
    manifest=json.loads((ROOT/'run_manifest_v2.json').read_text())
    audit={str(e['expnum']):audit_product(e['output']) for e in manifest}
    (ROOT/'median_fit_audit.json').write_text(json.dumps(audit,indent=2,allow_nan=False))
    print('Audited',len(audit),'exposures')
