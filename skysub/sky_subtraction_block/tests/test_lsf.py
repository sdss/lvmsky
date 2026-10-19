import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.interpolate import BSpline
from skysub.sky_subtraction_block.lsf import slit_basis,offset_centers,solve_normal,fit_block
from skysub.sky_decomp.lsf_spline2d import mspline_basis
from types import SimpleNamespace


def test_exact_continuous_constraints_and_injected_recovery():
    rng=np.random.default_rng(87);no,nw,ns=11,4,4
    centers=offset_centers(no,2.5)
    u=np.linspace(-2.5,2.5,20001)
    np.testing.assert_allclose(trapezoid(u[:,None]*mspline_basis(u,n_basis=no,half_width=2.5),u,axis=0),centers,atol=1e-11)
    # Symmetric but variable truth; all slit controls differ, even across the two telescope groups.
    truth=rng.uniform(.02,1,(no,nw,ns));truth=(truth+truth[::-1])/2;truth/=truth.sum(axis=0)
    positions=[1,30,100,200,350,500,600,648]
    basis,knots=slit_basis(positions,1,648,ns)
    samples=[];wave=np.linspace(6000,7000,240)
    for position,b in zip(positions,basis):
        design=rng.uniform(0,2,(len(wave),no*nw))
        flux=design@np.einsum('kqs,s->kq',truth,b).ravel()
        samples.append(dict(position=position,wave=wave,flux=flux,ivar=np.ones(len(wave)),mask=np.zeros(len(wave)),
            designs={'R':(np.ones(len(wave),bool),design,np.zeros(len(wave)))}))
    state=SimpleNamespace(coefficients={'R':truth[:,:,0]},config={'offset_half_width_angstrom':2.5})
    fitted,_,qa=fit_block(samples,state,'R',1,648,ns,centered=True)
    np.testing.assert_allclose(fitted,truth,atol=3e-5)
    assert qa['slit_rank']==4 and qa['n_fibers']==8
    wb=BSpline.design_matrix(np.linspace(0,1,31),[0]*4+[1]*4,3).toarray()
    sb=BSpline.design_matrix(np.linspace(1,648,47),knots,3).toarray()
    profiles=np.einsum('kqs,aq,bs->abk',fitted,wb,sb)
    assert profiles.min()>0
    np.testing.assert_allclose(profiles.sum(axis=-1),1,atol=1e-8)
    np.testing.assert_allclose(profiles@centers,0,atol=1e-8)


def test_centroid_constraint_is_optional_and_off_by_default():
    rng=np.random.default_rng(11);no,nw,ns=11,4,4;centers=offset_centers(no,2.5)
    # Truth with a nonzero centroid (shifted LSF): only the unconstrained default can recover it.
    truth=np.exp(-.5*((np.arange(no)-6.5)/1.5)**2)[:,None,None]*rng.uniform(.9,1.1,(no,nw,ns));truth/=truth.sum(axis=0)
    assert np.abs(np.einsum('k,kqs->qs',centers,truth)).min()>.05
    positions=[1,30,100,200,350,500,600,648];basis,_=slit_basis(positions,1,648,ns);wave=np.linspace(6000,7000,240)
    samples=[]
    for position,b in zip(positions,basis):
        design=rng.uniform(0,2,(len(wave),no*nw))
        samples.append(dict(position=position,wave=wave,flux=design@np.einsum('kqs,s->kq',truth,b).ravel(),ivar=np.ones(len(wave)),
            mask=np.zeros(len(wave)),designs={'R':(np.ones(len(wave),bool),design,np.zeros(len(wave)))}))
    state=SimpleNamespace(coefficients={'R':truth[:,:,0]},config={'offset_half_width_angstrom':2.5})
    free,_,qa=fit_block(samples,state,'R',1,648,ns)
    np.testing.assert_allclose(free,truth,atol=1e-4)
    assert not qa['centered'] and qa['centroid_error']>.05
    forced,_,qa=fit_block(samples,state,'R',1,648,ns,centered=True)
    assert qa['centered'] and qa['centroid_error']<1e-7 and np.abs(forced-truth).max()>1e-3


def test_input_and_solver_guards():
    with pytest.raises(ValueError):slit_basis([1],1,1)
    with pytest.raises(ValueError):slit_basis([0,2],0,2,3)
    with pytest.raises(ValueError):slit_basis([-1],0,2)
    # All-positive centroids cannot satisfy normalization and zero first moment.
    with pytest.raises(RuntimeError):solve_normal(np.eye(8),np.zeros(8),(2,2,2),np.ones(2))


def test_blue_channel_is_never_replaced_by_a_block():
    from skysub.sky_subtraction_block.lsf import evaluate_block, FIT_CHANNELS
    assert FIT_CHANNELS == ('R', 'Z')
    rng = np.random.default_rng(3)
    base = SimpleNamespace(coefficients={c: rng.uniform(size=(11, 4)) for c in 'BRZ'},
                           knot_vectors={c: np.r_[[6000.] * 4, [7000.] * 4] for c in 'BRZ'},
                           degrees={c: 3 for c in 'BRZ'}, metrics={})
    blue = base.coefficients['B'].copy()
    blue_knots = base.knot_vectors['B'].copy()
    knots = np.r_[[1.] * 4, [648.] * 4]
    fit_knots = np.r_[[5500.] * 4, [7500.] * 4]
    coef = np.zeros((11, 4, 4))
    coef[5] = 1.
    blocks = {(1, c): dict(coef=coef, slit_knots=knots, wave_knots=fit_knots, wave_degree=3, qa=dict(status='Solved'))
              for c in 'RZ'}
    state = evaluate_block(base, blocks, 1, 300.)
    np.testing.assert_array_equal(state.coefficients['B'], blue)
    np.testing.assert_array_equal(state.knot_vectors['B'], blue_knots)
    np.testing.assert_allclose(state.coefficients['R'].sum(axis=0), 1)
    # Block coefficients are interpreted only on the wavelength basis they were fitted on.
    for c in 'RZ':
        np.testing.assert_array_equal(state.knot_vectors[c], fit_knots)
    np.testing.assert_array_equal(base.knot_vectors['R'], np.r_[[6000.] * 4, [7000.] * 4])
    with pytest.raises(KeyError):
        evaluate_block(base, {(1, 'R'): blocks[(1, 'R')]}, 1, 300.)
    with pytest.raises(KeyError):
        evaluate_block(base, {(1, c): dict(coef=coef, slit_knots=knots, qa=dict(status='Solved')) for c in 'RZ'}, 1, 300.)


def test_fwhm_batch_matches_profile_diagnostics():
    from skysub.sky_subtraction_block.lsf import offset_centers
    from skysub.sky_subtraction_block.science import fwhm_batch
    from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics
    from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
    rng = np.random.default_rng(5)
    knots = np.r_[[5000.] * 4, [9500.] * 4]
    centers = offset_centers(11, 2.5)
    wave = np.array([6000., 7000., 8000.])
    basis = evaluate_bspline_basis(wave, knots, 3)
    basis = basis.toarray() if hasattr(basis, 'toarray') else basis
    sigma = rng.uniform(.4, 1.2, 4)
    coef = np.stack([np.exp(-.5 * (centers / s) ** 2) * rng.uniform(.8, 1.2, 11) for s in sigma])
    for c in coef:
        c = np.repeat(c[:, None], 4, axis=1)
        state = SimpleNamespace(coefficients={ch: c for ch in 'BRZ'}, knot_vectors={ch: knots for ch in 'BRZ'},
                                degrees={ch: 3 for ch in 'BRZ'}, config={'offset_half_width_angstrom': 2.5},
                                legacy_kernel_representation='lsf-spline2d')
        expected = evaluate_lsf_diagnostics(state, wave)['fwhm_angstrom']
        np.testing.assert_allclose(fwhm_batch(basis @ c.T, 2.5), expected, rtol=0, atol=1e-12)
