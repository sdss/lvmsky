"""Flux-space residuals and input error curves on the fitted spectrum axis."""
import numpy as np


def offset_residual(fig,wave,observed,model,sigma,used):
    from report_fibers import trace
    residual=np.asarray(observed)-np.asarray(model)
    good=np.asarray(used,bool)&np.isfinite(residual)&np.isfinite(sigma)&(sigma>0)
    if not np.any(good):return np.nan
    rms=float(np.sqrt(np.mean(residual[good]**2)));offset=-5*rms
    for sign in [-1,1]:
        trace(fig,wave,np.where(good,offset+sign*sigma,np.nan),'Offset ±1σ',1,color='black',showlegend=sign==1,width=.65)
    trace(fig,wave,np.where(good,offset,np.nan),'Residual baseline',1,color='#888888',dash='dot',showlegend=False,width=.6)
    trace(fig,wave,np.where(good,offset+residual,np.nan),'Data − model (−5 RMS)',1,color='#34526e',width=.9)
    fig.data[-1].customdata=np.where(good,residual,np.nan)
    fig.data[-1].hovertemplate=f'λ %{{x:.2f}} Å<br>Data − model %{{customdata:.4g}}<br>Offset {offset:.4g}<extra></extra>'
    return offset
