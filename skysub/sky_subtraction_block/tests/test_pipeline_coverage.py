"""Median coefficients are the default; individual Sky fits are opt-in."""
from types import SimpleNamespace
import numpy as np
from skysub.sky_subtraction_block import pipeline


def test_individual_unavailable_fit_is_recorded(monkeypatch):
    slit=np.array([('SkyE',36,1)],dtype=[('telescope','U8'),('fiberid','i4'),('spectrographid','i4')])
    model=SimpleNamespace(design_names=['lines','moon','zodi'])
    renderer=SimpleNamespace(model=lambda label:model,h={'MED_FIT_MASK':SimpleNamespace(data=np.zeros(8,bool))})
    work=dict(renderer=renderer,slit=slit,templates={'SkyE':object()},flux=np.ones((1,8)),
              ivar=np.ones((1,8)),mask=np.zeros((1,8),dtype='i4'),matrices={'SkyE':object()},fit_sky_fibers=True)
    monkeypatch.setattr(pipeline,'_WORK',work)
    def unavailable(*args,**kwargs):raise ValueError('insufficient usable pixels')
    monkeypatch.setattr(pipeline,'individual_coefficients',unavailable)
    result=pipeline._fit_sky(0)
    assert not result['usable'] and result['fiberid']==36 and np.isnan(result['coef']).all()
    assert result['reason']=='insufficient usable pixels' and result['status']=='FiberFitUnavailable'


def test_default_uses_own_telescope_median_without_fitting(monkeypatch):
    from inspect import signature
    assert signature(pipeline.run).parameters['fit_sky_fibers'].default is False
    slit=np.array([('SkyE',36,1),('SkyW',38,1),('SkyE',40,1)],
        dtype=[('telescope','U8'),('fiberid','i4'),('spectrographid','i4')])
    coefficients=np.array([[10.,20.,30.],[1.,2.,3.],[4.,5.,6.]])
    renderer=SimpleNamespace(labels=['Sci','SkyW','SkyE'],h={'MED_COEF':SimpleNamespace(data=coefficients)})
    monkeypatch.setattr(pipeline,'_WORK',dict(renderer=renderer,slit=slit,fit_sky_fibers=False))
    def forbidden(*args,**kwargs):raise AssertionError('Default must not fit individual Sky fibers')
    monkeypatch.setattr(pipeline,'individual_coefficients',forbidden)
    rows=[pipeline._fit_sky(i) for i in range(3)]
    for row,expected in zip(rows,[coefficients[2],coefficients[1],coefficients[2]]):
        np.testing.assert_array_equal(row['coef'],expected)
        assert row['usable'] and row['status']=='MedianCoefficients'
        assert not np.shares_memory(row['coef'],coefficients)
