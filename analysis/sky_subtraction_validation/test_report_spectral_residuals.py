"""Check offset/error arithmetic, contrast diagnostics, and requested windows."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from plotly.subplots import make_subplots
from report_spectral_residuals import offset_residual
from report_contrast import contrast_rows
from report_overview import WINDOW_GROUPS


def test_offset_is_five_fitted_residual_rms_and_errors_stay_in_flux_units():
    fig=make_subplots(rows=2,cols=1)
    wave=np.arange(4.);data=np.array([2.,4.,999.,8.]);model=np.array([1.,2.,3.,4.])
    sigma=np.array([.1,.2,.3,.4]);used=np.array([True,True,False,True])
    offset=offset_residual(fig,wave,data,model,sigma,used)
    np.testing.assert_allclose(offset,-5*np.sqrt(7))
    np.testing.assert_allclose(np.asarray(fig.data[0].y)[used],offset-sigma[used])
    np.testing.assert_allclose(np.asarray(fig.data[1].y)[used],offset+sigma[used])
    np.testing.assert_allclose(np.asarray(fig.data[-1].y)[used],offset+(data-model)[used])
    assert np.isnan(fig.data[-1].y[2]) and fig.data[-1].yaxis=='y'
    assert fig.data[0].line.color=='black' and fig.data[0].line.width<1


def test_contrast_uses_unmasked_separated_sky_peaks_not_total_line_plus_continuum():
    wave=np.arange(5787,9801,.5)
    centers=[5900,6100,6400,6900,7300,7600,8000,8500,9000,9500]
    lines=sum(10*np.exp(-((wave-c)/.7)**2) for c in centers)
    source=np.full_like(wave,20.);solar=np.full_like(wave,5.)
    used=np.ones(len(wave),bool);used[wave==5900]=False
    rows=contrast_rows(wave,lines,source,solar,used,source+solar)
    assert 5900 not in rows[0]['peaks']
    for row in rows:
        np.testing.assert_allclose(row['source_ratio'],.5)
        np.testing.assert_allclose(row['total_ratio'],.4)
        assert row['source_dominated']
        np.testing.assert_allclose(row['observed_ratio'],.4)
    assert not contrast_rows(wave,lines,np.ones_like(wave),solar,used,source+solar)[0]['source_dominated']


def test_requested_windows_include_ca_triplet_and_diagnostic_lines():
    def covered(group,wavelength):return any(lo<=wavelength<=hi for _,lo,hi in WINDOW_GROUPS[group])
    assert all(covered('blue',w) for w in [3727,3934,3968,5577])
    assert all(covered('red',w) for w in [6300,6563,6716,6731])
    assert all(covered('far_red',w) for w in [7600,8300,8498,8542,8662,9400])


def test_broad_windows_cover_red_chip_and_telluric_a_context():
    assert any(lo<=6500 and hi>=7454 for _,lo,hi in WINDOW_GROUPS['red'])
    assert any(lo<=7500 and hi>=8000 for _,lo,hi in WINDOW_GROUPS['far_red'])


def test_lsf_selection_diagram_states_actual_gates_and_scope():
    from report_overview import lsf_selection_diagram
    diagram=lsf_selection_diagram()
    assert all(text in diagram for text in ['B → survey LSF','R and Z → independent','35%','50%','0.25 Å','12 iterations','line-amplitude contrast ≥1'])
    assert '<details' not in diagram


def test_scheme_distinguishes_preparation_fitting_and_solver_support():
    from report_overview import lsf_selection_diagram
    d=lsf_selection_diagram()
    assert d.index('no fitting yet') < d.index('Check numerical support') < d.index('Fit starts here')
    assert 'solver precondition, not a quality test' in d
    assert 'each needs ≥2 independent sky features' in d
    assert '20 Å knot spacing' in d


def test_context_windows_are_at_least_400_angstrom_wide():
    for windows in WINDOW_GROUPS.values():
        assert all(hi-lo>=400 for _,lo,hi in windows)
