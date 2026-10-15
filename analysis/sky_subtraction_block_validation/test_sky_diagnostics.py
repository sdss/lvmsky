import numpy as np
import pytest
from analysis.sky_subtraction_block_validation.sky_diagnostics import shift_trial,residual_metrics


def test_linear_shift_recovers_signed_small_shift_and_improves_residual():
    wave=np.linspace(6000,6100,1001)
    model=2+np.exp(-.5*((wave-6050)/1.1)**2)
    derivative=np.gradient(model,wave)
    observed=model-.08*derivative
    sigma=np.full_like(wave,.01);good=np.ones(len(wave),bool)
    trial,delta=shift_trial(wave,observed,model,sigma,good)
    assert delta==pytest.approx(.08)
    np.testing.assert_allclose(trial,observed,atol=1e-14)
    assert residual_metrics(observed,trial,sigma,good)['chi2_per_pixel']<1e-20
    limited,delta=shift_trial(wave,model-derivative,model,sigma,good)
    assert delta==.25


def test_input_error_normalization_and_empty_selection():
    observed=np.array([2.,-2.,100.]);model=np.zeros(3);sigma=np.array([2.,2.,1.]);good=np.array([True,True,False])
    metrics=residual_metrics(observed,model,sigma,good)
    assert metrics['n_pixels']==2 and metrics['chi2_per_pixel']==1 and metrics['within_1sigma']==1
    with pytest.raises(ValueError):residual_metrics(observed,model,sigma,np.zeros(3,bool))
