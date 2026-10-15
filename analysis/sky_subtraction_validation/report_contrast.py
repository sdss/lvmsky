"""Trial-based sky-line contrast, for inspection rather than LSF acceptance."""
import numpy as np
from scipy.signal import find_peaks


def contrast_rows(wave,lines,source,solar,used,observed):
    peaks=find_peaks(np.where(np.isfinite(lines),lines,0))[0]
    output=[]
    for channel,lo,hi in [('R',5787,7454),('Z',7454,9801)]:
        candidates=peaks[used[peaks]&(wave[peaks]>=lo)&(wave[peaks]<hi)&(lines[peaks]>0)&np.isfinite(source[peaks])&np.isfinite(solar[peaks])]
        selected=[]
        for i in candidates[np.argsort(lines[candidates])[::-1]]:
            if all(abs(wave[i]-wave[j])>=50 for j in selected):selected.append(i)
            if len(selected)==5:break
        ix=np.asarray(selected,int)
        def ratio(continuum):
            return float(np.median(np.divide(lines[ix],continuum[ix],out=np.full(len(ix),np.inf),where=continuum[ix]>0))) if len(ix) else np.nan
        star=ratio(source);total=ratio(source+solar)
        local=[]
        for i in ix:
            side=used&(np.abs(wave-wave[i])>=2)&(np.abs(wave-wave[i])<=8)&np.isfinite(observed)
            continuum=float(np.median(observed[side])) if np.any(side) else np.nan
            if np.isfinite(continuum) and continuum>0:local.append(float(lines[i]/continuum))
        measured=float(np.median(local)) if local else np.nan
        output.append(dict(channel=channel,source_ratio=star,total_ratio=total,observed_ratio=measured,peaks=wave[ix].tolist(),
            source_dominated=bool(np.isfinite(star) and star<=1),low_contrast=bool(np.isfinite(measured) and measured<=1)))
    return output


def contrast_html(rows):
    def number(x):return '∞' if np.isinf(x) else ('—' if not np.isfinite(x) else f'{x:.2g}')
    body='<p>Contrast diagnostic: median fitted emission-peak / observed local continuum, using up to five peaks separated by ≥50 Å (2–8 Å sidebands). This strongest-peak summary is supplemental; adoption follows the separate full-channel gate. Stellar absorption can bias individual local measurements.</p><table><tr><th>Channel</th><th>Sky peak / observed continuum</th><th>Sky peak / fitted nuisance continuum</th><th>Peak wavelengths (Å)</th><th>Review</th></tr>'
    for r in rows:
        body+=f'<tr><td>{r["channel"]}</td><td>{number(r["observed_ratio"])}</td><td>{number(r["source_ratio"])}</td><td>{", ".join(f"{w:.1f}" for w in r["peaks"])}</td><td>{"Low observed contrast: inspect fit" if r["low_contrast"] else "Inspect spectrum"}</td></tr>'
    return body+'</table>'
