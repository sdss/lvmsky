"""Repeat numerical recovery checks and record their acceptance/error outcomes."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import os
os.environ.setdefault('MPLBACKEND','Agg')
import sys
import numpy as np
from astropy.table import Table

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))
path=ROOT.parents[1]/'skysub/sky_subtraction/tests/test_core.py'
spec=importlib.util.spec_from_file_location('injection_tests',path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
original=module.fit_channel
records=[]
case={}


def record_fit(*args,**kwargs):
    fit=original(*args,**kwargs)
    candidate,info=fit.candidate,fit.qa
    wave,_,_,_,reference,channel=args[:6]
    truth=deepcopy(reference)
    if case['kind']=='known_broader':
        c=reference.coefficients[channel]
        c=.5*c+.25*np.roll(c,1,axis=0)+.25*np.roll(c,-1,axis=0)
        truth.coefficients[channel]=c/c.sum(axis=0)
    probes=np.linspace(wave[0]+1,wave[-1]-1,9)
    d=module.evaluate_lsf_diagnostics(candidate,probes)
    t=module.evaluate_lsf_diagnostics(truth,probes)
    error=float(np.max(np.abs(d['w50_angstrom']/t['w50_angstrom']-1)))
    records.append(dict(**case,channel=channel,accepted=info['accepted'],reason=info['reason'],max_w50_error=error))
    return fit


module.fit_channel=record_fit
for channel,lo,hi in [('R',5787.,7454.),('Z',7454.,9800.)]:
    for contrast in [0.,1.,10.]:
        case=dict(kind='known_broader',contrast=contrast)
        module.test_known_lsf_under_stellar_structure(channel,lo,hi,contrast)
    for contrast in [.3,3.,10.]:
        case=dict(kind='empirical_M5V',contrast=contrast)
        module.test_empirical_m5v_bands(channel,lo,hi,contrast)
assert sum(r['accepted'] for r in records if r['contrast']==0)==2
Table(rows=records).write(ROOT/'injection_metrics.csv',format='ascii.csv',overwrite=True)
print(Table(rows=records))
