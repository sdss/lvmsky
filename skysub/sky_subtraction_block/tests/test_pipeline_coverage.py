"""Unavailable individual Sky fits are explicit data outcomes."""
from types import SimpleNamespace
import numpy as np
from skysub.sky_subtraction_block import pipeline


def test_individual_unavailable_fit_is_recorded(monkeypatch):
    slit=np.array([('SkyE',36,1)],dtype=[('telescope','U8'),('fiberid','i4'),('spectrographid','i4')])
    model=SimpleNamespace(design_names=['lines','moon','zodi'])
    renderer=SimpleNamespace(model=lambda label:model,h={'MED_FIT_MASK':SimpleNamespace(data=np.zeros(8,bool))})
    work=dict(renderer=renderer,slit=slit,templates={'SkyE':object()},flux=np.ones((1,8)),
              ivar=np.ones((1,8)),mask=np.zeros((1,8),dtype='i4'),matrices={'SkyE':object()})
    monkeypatch.setattr(pipeline,'_WORK',work)
    def unavailable(*args,**kwargs):raise ValueError('insufficient usable pixels')
    monkeypatch.setattr(pipeline,'individual_coefficients',unavailable)
    result=pipeline._fit_sky(0)
    assert not result['usable'] and result['fiberid']==36 and np.isnan(result['coef']).all()
    assert result['reason']=='insufficient usable pixels' and result['status']=='FiberFitUnavailable'
