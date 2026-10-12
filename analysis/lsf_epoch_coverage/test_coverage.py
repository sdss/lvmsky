"""Small checks for strict contamination cuts and independent-night metrics."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from build_coverage import summarize, summarize_groups, gaia_ratio_clean_mask, gaia_path
import numpy as np
import pandas as pd


def test_strict_gaia_selection():
    data=dict(n_gaia_sources=np.array([1,1,0,-1,1]),
        gaia_g_n_valid_sources=np.array([1,1,0,-1,1]),
        lvm_flux_plus_sky_g=np.ones(5),lvm_flux_plus_sky_g_valid_fraction=np.array([1,1,1,1,.94]),
        ratio_gaia_to_lvm_flux_plus_sky_g=np.array([.009,.01,0,np.nan,.001]))
    assert gaia_ratio_clean_mask(data,np.nextafter(.01,0)).tolist()==[True,False,True,False,False]


def test_night_concentration_and_empty_selection():
    df=pd.DataFrame(dict(accepted=[True,True,True],gaia_verified=[True]*3,gaia_clean_1pct=[True]*3,
        mjd=[60000,60000,60001],exptime_seconds=[900.,900.,1800.],
        observation_mjd_utc=[60000.,60000.01,60001.],days_from_full_moon=[-1.,-.99,0.],
        moon_fli=[.99, .99, 1.],synthetic_g=[1.,1.,1.],synthetic_g_error=[.1,.1,.1],chi2_red=[1.,2.,3.]))
    result=summarize(df)
    df['key']=1
    fast=summarize_groups(df,['key']).iloc[0]
    for col in result.index:
        assert np.isclose(result[col],fast[col],equal_nan=True), col
    assert result.n_nights==2 and result.effective_nights==2
    assert result.largest_night_share==.5 and result.total_hours==1
    assert result.repeat_nights==1 and result.median_lsf_chi2==2
    df['accepted']=False
    result=summarize(df)
    fast=summarize_groups(df,['key']).iloc[0]
    for col in result.index:
        assert np.isclose(result[col],fast[col],equal_nan=True), col
    assert result.n_accepted==0 and result.n_nights==0 and result.effective_nights==0
    assert np.isnan(result.median_lsf_chi2)


def test_gaia_source_guard():
    import build_coverage
    original=build_coverage.CACHE
    try:
        build_coverage.CACHE=Path('/data/sas/sdsswork/lvm/sandbox/calib/gaia/exposure/fibers')
        try:gaia_path(3983)
        except ValueError:pass
        else:raise AssertionError('An unrelated Gaia source must be rejected')
    finally:build_coverage.CACHE=original
    assert gaia_path(3983).is_relative_to(build_coverage.DATA)

if __name__=='__main__':
    test_strict_gaia_selection();test_night_concentration_and_empty_selection()
    test_gaia_source_guard()
    print('Coverage checks passed')
