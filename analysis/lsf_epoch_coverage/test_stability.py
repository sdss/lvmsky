"""Check LSF diagnostic interpolation and separation of physical fibers."""
from types import SimpleNamespace
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_stability import diagnostics, scatter, METRICS, evaluate_lsf_diagnostics, LSF_SPLINE2D_REPRESENTATION
import numpy as np
import pandas as pd


def test_diagnostics_equivalence():
    rng=np.random.default_rng(52)
    coef=rng.dirichlet(np.ones(11),size=4).T
    state=SimpleNamespace(coefficients={c:coef for c in ['B','R','Z']},
        knot_vectors={c:np.array([3600.]*4+[5786.5]*4) for c in ['B','R','Z']},
        degrees={c:3 for c in ['B','R','Z']},config={'offset_half_width_angstrom':2.5},
        legacy_kernel_representation=LSF_SPLINE2D_REPRESENTATION)
    wave=np.array([3800.,4700.,5650.])
    from analyze_stability import evaluate_bspline_basis
    mass=evaluate_bspline_basis(wave,state.knot_vectors['B'],3)@coef.T
    reference=evaluate_lsf_diagnostics(state,wave)
    actual=diagnostics(mass,2.5)
    for k in actual:assert np.allclose(actual[k],reference[k],rtol=0,atol=1e-9),k


def test_fiber_specific_centering():
    df=pd.DataFrame(dict(fiberid=[1,1,1,2,2,2],expnum=[1,2,3]*2,mjd=[60000]*6,
        observation_mjd_utc=[60000.,60000.01,60000.02]*2,g_snr=[100.]*6,chi2_red=[1.]*6))
    for c in METRICS:df[c]=[1.,2.,3.,10.,20.,30.]
    result=scatter(df,['fiberid']).set_index('fiberid')
    assert result.loc[1,'fwhm_angstrom_median']==2 and result.loc[2,'fwhm_angstrom_median']==20
    assert np.isclose(result.loc[1,'fwhm_angstrom_rmad_pct'],74.13)
    assert np.isclose(result.loc[1,'fwhm_angstrom_rmad_pct'],result.loc[2,'fwhm_angstrom_rmad_pct'])
    assert result.n_exposures.eq(3).all() and result.n_nights.eq(1).all()


if __name__=='__main__':
    test_diagnostics_equivalence();test_fiber_specific_centering()
    print('Stability checks passed')
