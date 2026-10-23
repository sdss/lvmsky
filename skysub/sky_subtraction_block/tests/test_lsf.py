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


def test_regularization_preserves_data_objective_and_shrinks_roughness():
    from skysub.sky_subtraction_block.lsf import regularized_normal
    shape = (7, 4, 4)
    rng = np.random.default_rng(8)
    truth = rng.uniform(.01, 1, shape); truth /= truth.sum(axis=0)
    hessian = np.eye(truth.size); linear = -truth.ravel()
    reference = np.ones(shape[:2])/shape[0]
    zero, vector = regularized_normal(hessian, linear, shape, reference)
    np.testing.assert_array_equal(zero, hessian)
    np.testing.assert_array_equal(vector, linear)
    penalized, vector = regularized_normal(hessian, linear, shape, reference, 100., 10., 10., 0.)
    fitted, _ = solve_normal(penalized, vector, shape)
    assert np.linalg.norm(np.diff(fitted, n=2, axis=0)) < np.linalg.norm(np.diff(truth, n=2, axis=0))
    np.testing.assert_array_equal(hessian, np.eye(truth.size))
    np.testing.assert_allclose(fitted.sum(axis=0), 1., atol=1e-8)


def _coverage_samples():
    rng=np.random.default_rng(52); no,nw,ns=7,4,4
    truth=rng.uniform(.05,1,(no,nw,ns));truth/=truth.sum(axis=0)
    positions=np.linspace(1,648,16);basis,_=slit_basis(positions,1,648,ns)
    wave=np.linspace(6000,7000,100);samples=[]
    for position,b in zip(positions,basis):
        design=rng.uniform(.1,2,(len(wave),no*nw))
        flux=design@np.einsum('kqs,s->kq',truth,b).ravel()
        samples.append(dict(position=position,wave=wave,flux=flux,ivar=np.ones(len(wave)),mask=np.zeros(len(wave)),
            designs={c:(np.ones(len(wave),bool),design,np.zeros(len(wave))) for c in 'RZ'}))
    state=SimpleNamespace(coefficients={c:truth[:,:,0] for c in 'RZ'},config={'offset_half_width_angstrom':2.5})
    return samples,state,truth,positions


def test_local_effective_df_handles_smoothing_equalities_and_active_faces():
    from skysub.sky_subtraction_block.lsf import local_effective_df
    coef=np.full((3,2,2),1/3);hessian=np.eye(coef.size)
    for penalty in [0.,.1,1.,100.]:
        df,active,_=local_effective_df(hessian,hessian*(1+penalty),coef)
        assert active==12
        np.testing.assert_allclose(df,8/(1+penalty),atol=1e-12)
    df,_,_=local_effective_df(hessian,hessian,coef,np.array([-1.,0.,1.]))
    np.testing.assert_allclose(df,4,atol=1e-12)
    coef[:]=0;coef[1]=1
    df,active,_=local_effective_df(hessian,hessian,coef)
    assert df==0 and active==4


@pytest.mark.parametrize('case,reason',[
    ('empty','insufficient_independent_sky_positions'),
    ('few','insufficient_independent_sky_positions'),
    ('boundary','unsupported_slit_boundary'),
    ('gap','large_sky_coverage_gap'),
    ('bad_channel','insufficient_independent_sky_positions')])
def test_sparse_or_masked_sky_uses_explicit_reference(case,reason):
    from skysub.sky_subtraction_block.lsf import fit_supported_block
    samples,state,truth,positions=_coverage_samples()
    if case=='empty':samples=[]
    elif case=='few':samples=samples[:3]
    elif case=='boundary':samples=samples[1:]
    elif case=='gap':samples=samples[:4]+samples[-4:]
    elif case=='bad_channel':
        for s in samples:s['ivar'][:]=0
    result,_,qa=fit_supported_block(samples,state,'R',1,648,expected_positions=positions)
    np.testing.assert_array_equal(result,np.repeat(state.coefficients['R'][:,:,None],4,axis=2))
    assert not qa['adopted'] and not qa['centered'] and qa['fallback_reason']==reason
    assert qa['status']=='ReferenceFallback' and np.isnan(qa['chi2_per_pixel'])
    np.testing.assert_allclose(result.sum(axis=0),1,atol=1e-12)


def test_supported_sky_fit_matches_original_and_survives_one_fiber_loss():
    from skysub.sky_subtraction_block.lsf import fit_supported_block
    samples,state,truth,positions=_coverage_samples()
    supported,_,qa=fit_supported_block(samples,state,'R',1,648,expected_positions=positions)
    original,_,_=fit_block(samples,state,'R',1,648)
    np.testing.assert_array_equal(supported,original)
    assert qa['adopted'] and qa['fallback_reason']==''
    samples.pop(7)
    partial,_,qa=fit_supported_block(samples,state,'R',1,648,expected_positions=positions)
    assert qa['adopted']
    np.testing.assert_allclose(partial,truth,atol=3e-4)


def test_channel_coverage_and_solver_fallback_are_independent(monkeypatch):
    import skysub.sky_subtraction_block.lsf as module
    samples,state,_,positions=_coverage_samples()
    # A fully masked R design is unusable while Z keeps the same physical fiber.
    take,design,fixed=samples[0]['designs']['R']
    samples[0]['designs']['R']=(take,design,np.full_like(fixed,np.nan))
    _,_,red=module.fit_supported_block(samples,state,'R',1,648,expected_positions=positions)
    _,_,zed=module.fit_supported_block(samples,state,'Z',1,648,expected_positions=positions)
    assert red['fallback_reason']=='unsupported_slit_boundary' and red['rejected_channel_fibers']==1
    assert zed['adopted'] and zed['rejected_channel_fibers']==0
    def unavailable(*args,**kwargs):raise RuntimeError('injected solver failure')
    monkeypatch.setattr(module,'fit_block',unavailable)
    _,_,qa=module.fit_supported_block(samples,state,'Z',1,648,expected_positions=positions)
    assert not qa['adopted'] and qa['fallback_reason']=='block_solver_failure'


@pytest.mark.parametrize('offset_count', [7, 15])
def test_mixed_basis_roundtrip_and_exact_native_rendering(offset_count):
    from pathlib import Path
    from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
    from skysub.sky_decomp.result_io import build_lsf_hdus, lsf_surface_state_from_cubes
    from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics, evaluate_lsf_density
    from skysub.sky_subtraction_block.lsf import resize_state
    from skysub.sky_subtraction.tests.test_product import _small_model
    from skysub.sky_subtraction.rendering import render
    from skysub.sky_subtraction.fiber_lsf import sky_designs
    base = load_survey_ascii_lsf(Path(__file__).parents[2]/'sky_decomp/data/lsf', 500)
    state = resize_state(base, offset_count, 6)
    hdus = build_lsf_hdus([state])
    restored = lsf_surface_state_from_cubes(tuple(h.data for h in hdus), 0)
    for ch in 'BRZ':
        np.testing.assert_array_equal(restored.coefficients[ch], state.coefficients[ch])
    np.testing.assert_array_equal(state.coefficients['B'], base.coefficients['B'])
    probes = np.array([4500., 6500., 8500.])
    diagnostics = evaluate_lsf_diagnostics(state, probes)
    density = evaluate_lsf_density(state, probes, np.linspace(-2.5, 2.5, 1201))
    np.testing.assert_allclose(trapezoid(density, np.linspace(-2.5, 2.5, 1201)), 1., atol=1e-8)
    assert np.isfinite(diagnostics['fwhm_angstrom']).all()
    model = _small_model(); coef = np.full(len(model.design_names), .001)
    expected, _, matrices = render(model, coef, restored)
    for ch, (take, parts) in sky_designs(model, coef, restored, matrices).items():
        actual = (parts['lines']+parts['solar'])@restored.coefficients[ch].ravel()+parts['fixed']
        np.testing.assert_allclose(actual, expected[take], rtol=1e-10, atol=1e-12)


def test_failed_arms_from_the_flag_or_from_older_preparations_r2():
    from skysub.sky_subtraction_block.pipeline import failed_arms

    summary=lambda r2:f'status=Solved | npar=388 | ngood=12307 | chi2_red=19.7 | R2={r2} | line_weight=0.0005'
    # Recorded flag wins over the summary.
    assert failed_arms({'median_fits':{'SkyE':dict(failed=True,summary=summary(0.99)),
                                       'SkyW':dict(failed=False,summary=summary(-5.0))}})==('SkyE',)
    # Preparations written before the flag: the R2 of the fit summary (13190: SkyE -16.35033).
    old={'median_fits':{'SkyE':dict(summary=summary(-16.35033)),'SkyW':dict(summary=summary(0.99628))}}
    assert failed_arms(old)==('SkyE',)
    assert failed_arms({'median_fits':{'SkyE':dict(summary=summary(0.9)),'SkyW':dict(summary=summary('nan'))}})==()
    with pytest.raises(ValueError,match='Both'):
        failed_arms({'median_fits':{l:dict(failed=True) for l in ('SkyE','SkyW')}})


def test_older_preparations_get_the_moon_glare_test_from_their_medians():
    from skysub.sky_subtraction_block.pipeline import failed_arms

    wave=np.linspace(3600,9800,2000); sky=(wave/4500)**-1.5; glare=sky*(wave/4500)**2.5
    summary=lambda r2:f'status=Solved | chi2_red=2 | R2={r2} | line_weight=0.0005'
    prov={'median_fits':{'SkyE':dict(summary=summary(0.96)),'SkyW':dict(summary=summary(0.997))}}   # exposure 13248
    meta={'skye_moon_sep':6.0,'skyw_moon_sep':62.3,'moon_alt':47.6,'moon_fli':0.94}
    assert failed_arms(prov)==()                                            # R2 alone misses it
    assert failed_arms(prov,(wave,{'SkyE':glare,'SkyW':sky},meta))==('SkyE',)
    assert failed_arms(prov,(wave,{'SkyE':glare,'SkyW':sky},dict(meta,moon_alt=-5.0)))==()


def test_shared_metadata_checks_only_the_preserved_blue_channel():
    from skysub.sky_subtraction_block.science import uniform_metadata
    from astropy.io import fits
    from astropy.table import Table

    rows = [dict(spectrum_index=i, channel=ch, status="survey" if ch == "B" else f"block_{i}",
                 value=float("nan")) for i in range(2) for ch in "BRZ"]
    hdus = fits.HDUList([fits.PrimaryHDU(),
                        fits.BinTableHDU(Table(rows=rows), name="LSF_META"),
                        fits.BinTableHDU(Table(dict(fiberid=[1, 2])), name="LSF_FIBERS")])
    reader = SimpleNamespace(h=hdus)
    assert uniform_metadata(reader)
    hdus["LSF_META"].data["status"][3] = "changed"
    assert not uniform_metadata(reader)
