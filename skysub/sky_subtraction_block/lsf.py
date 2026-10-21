"""One convex LSF fit per spectrograph/channel, with frozen fiber populations."""
from copy import deepcopy
import time
import clarabel
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh, null_space
from scipy.interpolate import BSpline
from skysub.sky_decomp.lsf_spline2d import _offset_knots, mspline_basis, evaluate_lsf_density

FIT_CHANNELS = ('R', 'Z')


def slit_basis(position, lower, upper, count=4):
    if count < 4 or not upper > lower:
        raise ValueError('Cubic slit basis needs count >= 4 and a finite nonempty domain')
    knots = np.r_[np.repeat(lower, 4), np.linspace(lower, upper, count-2)[1:-1], np.repeat(upper, 4)]
    return BSpline.design_matrix(position, knots, 3).toarray(), knots


def offset_centers(count, half_width):
    knots = _offset_knots(count, half_width)
    # A unit-integral cubic M-spline's first moment is its five-knot mean.
    return np.array([np.mean(knots[k:k+5]) for k in range(count)])


def solve_normal(hessian, linear, shape, centers=None):
    """Enforce unit integral on every wavelength/slit control pair, and zero centroid when ``centers`` is given."""
    no, nw, ns = shape
    size = no*nw*ns
    identity = sp.eye(nw*ns, format='csc')
    rows = [sp.kron(np.ones((1,no)), identity)]
    if centers is not None:
        rows.append(sp.kron(centers[None,:], identity))
    equality = sp.vstack(rows, format='csc')
    constraint = sp.vstack([equality, -sp.eye(size)], format='csc')
    rhs = np.r_[np.ones(nw*ns), np.zeros((len(rows)-1)*nw*ns+size)]
    scale = max(float(np.trace(hessian)/size), 1e-30)
    settings = clarabel.DefaultSettings(); settings.verbose = False
    settings.tol_gap_abs = settings.tol_feas = settings.tol_gap_rel = 1e-9
    solution = clarabel.DefaultSolver(sp.triu(sp.csc_matrix(hessian/scale)).tocsc(), linear/scale,
        constraint, rhs, [clarabel.ZeroConeT(len(rows)*nw*ns), clarabel.NonnegativeConeT(size)], settings).solve()
    if str(solution.status) not in ('Solved', 'AlmostSolved'):
        raise RuntimeError(f'Block LSF QP: {solution.status}')
    coefficient = np.asarray(solution.x).reshape(shape)
    # Keep the equality solution intact; clipping and renormalizing alter the first moment.
    if (not np.isfinite(coefficient).all() or coefficient.min() < -1e-8
        or np.max(np.abs(coefficient.sum(axis=0)-1)) > 1e-7
        or (centers is not None and np.max(np.abs(np.einsum('k,kqs->qs',centers,coefficient))) > 1e-7)):
        raise ValueError('Block LSF solution violates positivity or moment constraints')
    return coefficient, str(solution.status)


def component_matrices(model, state):
    """Component matrices of a model and LSF state; independent of the coefficients, so shared by all fibers of a telescope."""
    from skysub.sky_subtraction.pipeline import render
    return render(model, np.zeros(len(model.design_names)), state)[2]


def channel_tensors(model, coef, state, matrices=None):
    """Affine native-pixel designs ``{channel: (take, design, fixed)}`` of the shared sky renderer, R and Z only."""
    from skysub.sky_subtraction.fiber_lsf import sky_designs
    if matrices is None:
        matrices = component_matrices(model, state)
    return {ch: (take, parts['lines']+parts['solar'], parts['fixed'])
            for ch, (take, parts) in sky_designs(model, coef, state, matrices).items()}


def individual_coefficients(model, state, flux, ivar, mask, matrices=None):
    """One full component QP with fixed median LSF and median O2 temperature."""
    from skysub.sky_subtraction.pipeline import FAMILIES
    if matrices is None:
        matrices = component_matrices(model, state)
    slices = model._component_slices(matrices)
    good = (mask==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)
    if good.sum() <= len(model.design_names):
        raise ValueError('Sky fiber has insufficient valid pixels for decomposition')
    result = model._fit_design(np.vstack([matrices[k] for k in FAMILIES]),
        np.where(good,flux,0),np.where(good,ivar,0),
        moon_slice=slices['moon'],zodi_slice=slices['zodi'],diffuse_slice=slices['diffuse'])
    if str(result['status']) not in ('Solved','AlmostSolved') or not np.isfinite(result['coef']).all():
        raise RuntimeError(f"Sky-fiber coefficient QP: {result['status']}")
    return np.asarray(result['coef']), str(result['status'])


def resize_state(base, offset_count, wave_count):
    """Project median R/Z profiles onto the requested basis; preserve B exactly."""
    state = deepcopy(base)
    for channel in FIT_CHANNELS:
        if base.coefficients[channel].shape == (offset_count, wave_count):
            continue
        old = base.knot_vectors[channel]; degree = base.degrees[channel]
        lower, upper = old[degree], old[-degree-1]
        wave = np.linspace(lower, upper, max(32, wave_count*4))
        basis, knots = slit_basis(wave, lower, upper, wave_count)
        delta = np.linspace(-base.config['offset_half_width_angstrom'], base.config['offset_half_width_angstrom'], 241)
        offset = mspline_basis(delta, n_basis=offset_count, half_width=base.config['offset_half_width_angstrom'])
        design = np.einsum('dk,wq->wdkq', offset, basis).reshape(-1, offset_count*wave_count)
        target = evaluate_lsf_density(base, wave, delta).ravel()
        coefficient, _ = solve_normal(design.T@design, -design.T@target, (offset_count, wave_count, 1))
        state.coefficients[channel] = coefficient[:, :, 0]
        state.knot_vectors[channel] = knots
        state.degrees[channel] = 3
    count = max(c.shape[0] for c in state.coefficients.values())
    state.tap_offsets = np.arange(-(count//2), count//2+1)
    return state


def regularized_normal(hessian, linear, shape, reference, offset=0., wave=0., slit=0., prior=0.):
    """Dimensionless second differences and median-profile prior, scaled by mean data curvature."""
    scale = max(float(np.trace(hessian)/len(linear)), 1e-30)
    penalty = np.zeros_like(hessian)
    for axis, strength in enumerate((offset, wave, slit)):
        if strength:
            factors = [np.eye(n) for n in shape]
            factors[axis] = np.diff(factors[axis], n=2, axis=0)
            operator = np.kron(np.kron(factors[0], factors[1]), factors[2])
            penalty += strength * (operator.T@operator)
    reference = np.broadcast_to(reference[:, :, None], shape).ravel()
    if prior:
        penalty += prior*np.eye(len(linear))
    return hessian+scale*penalty, linear-scale*prior*reference


def _sample_pixels(sample, channel):
    take, design, fixed = sample['designs'][channel]
    wave = sample['wave'][take]
    y = sample['flux'][take]-fixed
    ivar = sample['ivar'][take]
    good = (sample['mask'][take]==0)&np.isfinite(y)&np.isfinite(ivar)&(ivar>0)
    good &= (wave>=wave[0]+5)&(wave<=wave[-1]-5)
    return design, y, ivar, good


def local_effective_df(hessian, penalized, coefficient, centers=None):
    """Local smoother trace on the positive-coefficient face and equality tangent.

    This conditions on the active face and fixed nuisance coefficients. It does
    not include uncertainty from selecting a face, regularizer or preparation.
    """
    no,nw,ns=coefficient.shape
    threshold=max(1e-8,1e-6*float(coefficient.max()))
    active=coefficient.ravel()>threshold
    identity=np.eye(nw*ns)
    equality=np.kron(np.ones((1,no)),identity)
    if centers is not None:equality=np.vstack([equality,np.kron(centers[None,:],identity)])
    tangent=null_space(equality[:,active])
    if not tangent.shape[1]:return 0.,int(active.sum()),threshold
    data=tangent.T@hessian[np.ix_(active,active)]@tangent
    total=tangent.T@penalized[np.ix_(active,active)]@tangent
    values,vectors=eigh((total+total.T)/2)
    take=values>max(float(values.max())*1e-10,1e-30)
    vectors=vectors[:,take]
    trace=float(np.sum(np.sum(vectors*(data@vectors),axis=0)/values[take]))
    return float(np.clip(trace,0,tangent.shape[1])),int(active.sum()),threshold


def fit_supported_block(samples, state, channel, lower, upper, slit_count=4, centered=False,
                        offset_roughness=0., wave_roughness=0., slit_roughness=0., reference_prior=0.,
                        expected_positions=None):
    """Fit a covered block, or adopt its projected Science-median reference.

    Coverage is checked before solving, using channel-specific valid samples.
    Entire-block fallback avoids discontinuities across unsupported slit intervals.
    These conservative geometry checks establish support, not physical LSF accuracy.
    """
    no, nw = state.coefficients[channel].shape
    usable = [s for s in samples if _sample_pixels(s, channel)[3].sum() > no*nw]
    positions = np.unique([s['position'] for s in usable])
    full_basis, knots = slit_basis(positions if len(positions) else [lower], lower, upper, slit_count)
    basis = full_basis if len(positions) else full_basis[:0]
    rank = int(np.linalg.matrix_rank(basis)) if len(positions) else 0
    condition = float(np.linalg.cond(basis)) if rank == slit_count else np.inf
    expected = np.unique(expected_positions if expected_positions is not None else [lower, upper])
    nominal_gap = float(np.diff(expected).max()) if len(expected)>1 else upper-lower
    max_gap = float(np.diff(positions).max()) if len(positions)>1 else upper-lower
    minimum = float(positions.min()) if len(positions) else np.nan
    maximum = float(positions.max()) if len(positions) else np.nan
    reason = ''
    if len(positions) < 2*slit_count:
        reason = 'insufficient_independent_sky_positions'
    elif rank < slit_count or not np.isfinite(condition) or condition > 1e3:
        reason = 'unidentified_slit_basis'
    elif minimum > lower or maximum < upper:
        reason = 'unsupported_slit_boundary'
    elif max_gap > 2*nominal_gap:
        reason = 'large_sky_coverage_gap'
    coverage = dict(adopted=not bool(reason), fallback_reason=reason, requested_centered=bool(centered),
        input_fibers=len(samples), rejected_channel_fibers=len(samples)-len(usable),
        slit_condition=condition, expected_max_gap=nominal_gap, selected_max_gap=max_gap)
    if not reason:
        try:
            coefficient, knots, qa = fit_block(usable, state, channel, lower, upper, slit_count, centered,
                offset_roughness, wave_roughness, slit_roughness, reference_prior)
            qa.update(coverage)
            qa.update(data_chi2=qa['chi2'],data_n_pixels=qa['n_pixels'])
            return coefficient, knots, qa
        except (RuntimeError, ValueError):
            # Numerical solver/constraint failure has the same explicit reference path.
            reason = 'block_solver_failure'
            coverage.update(adopted=False, fallback_reason=reason)
    coefficient = np.repeat(state.coefficients[channel][:, :, None], slit_count, axis=2)
    reference_chi2=0.;reference_pixels=0
    for sample in usable:
        design,y,ivar,good=_sample_pixels(sample,channel)
        residual=y[good]-design[good]@state.coefficients[channel].ravel()
        reference_chi2+=float(residual@(ivar[good]*residual));reference_pixels+=int(good.sum())
    centers = offset_centers(no, state.config['offset_half_width_angstrom'])
    qa = dict(status='ReferenceFallback', centered=False, n_fibers=len(usable),
        n_pixels=sum(int(_sample_pixels(s, channel)[3].sum()) for s in usable),
        normal_seconds=0., solve_seconds=0., chi2=np.nan, chi2_per_pixel=np.nan,
        slit_rank=rank, selected_min=minimum, selected_max=maximum,
        normalization_error=float(np.max(abs(coefficient.sum(axis=0)-1))),
        centroid_error=float(np.max(abs(np.einsum('k,kqs->qs', centers, coefficient)))),
        coefficient_min=float(coefficient.min()),data_chi2=reference_chi2,data_n_pixels=reference_pixels,
        nominal_df=0,effective_df=0.,active_coefficients=0,active_threshold=0.,df_seconds=0.,
        hessian_bytes=0,penalized_bytes=0,coefficient_bytes=coefficient.nbytes,
        sample_design_bytes=sum(s['designs'][channel][1].nbytes for s in usable),**coverage)
    return coefficient, knots, qa


def fit_block(samples, state, channel, lower, upper, slit_count=4, centered=False,
              offset_roughness=0., wave_roughness=0., slit_roughness=0., reference_prior=0.):
    """Accumulate small normal equations without a full multi-fiber design array.

    ``centered`` adds the zero-first-moment equality; ``centroid_error`` is reported either way."""
    started=time.perf_counter()
    no,nw = state.coefficients[channel].shape
    shape = (no,nw,slit_count); size = int(np.prod(shape))
    hessian = np.zeros((size,size)); linear = np.zeros(size)
    basis,knots = slit_basis([s['position'] for s in samples],lower,upper,slit_count)
    counts=[]; chi_constant=0
    for sample,weights in zip(samples,basis):
        design,y,ivar,good = _sample_pixels(sample,channel)
        if good.sum() <= no*nw:
            raise ValueError('Sky fiber has insufficient valid channel pixels')
        root = np.sqrt(ivar[good]); a=design[good]*root[:,None]; target=y[good]*root
        normal = a.T@a; vector = a.T@target
        hessian += np.kron(normal, np.outer(weights,weights))
        linear -= np.kron(vector,weights)
        chi_constant += target@target; counts.append(int(good.sum()))
    centers=offset_centers(no,state.config['offset_half_width_angstrom'])
    normal_seconds=time.perf_counter()-started
    solve_started=time.perf_counter()
    penalized, objective = regularized_normal(hessian, linear, shape, state.coefficients[channel],
        offset_roughness, wave_roughness, slit_roughness, reference_prior)
    coefficient,status = solve_normal(penalized,objective,shape,centers if centered else None)
    solve_seconds=time.perf_counter()-solve_started
    vector=coefficient.ravel()
    chi2=float(vector@hessian@vector+2*linear@vector+chi_constant)
    df_started=time.perf_counter()
    effective,active,threshold=local_effective_df(hessian,penalized,coefficient,centers if centered else None)
    df_seconds=time.perf_counter()-df_started
    return coefficient,knots,dict(status=status,centered=bool(centered),n_fibers=len(samples),n_pixels=sum(counts),
        normal_seconds=normal_seconds,solve_seconds=solve_seconds,chi2=chi2,chi2_per_pixel=chi2/sum(counts),slit_rank=int(np.linalg.matrix_rank(basis)),
        selected_min=float(min(s['position'] for s in samples)),selected_max=float(max(s['position'] for s in samples)),
        normalization_error=float(np.max(abs(coefficient.sum(axis=0)-1))),
        centroid_error=float(np.max(abs(np.einsum('k,kqs->qs',centers,coefficient)))),coefficient_min=float(coefficient.min()),
        nominal_df=(no-1-int(centered))*nw*slit_count,effective_df=effective,active_coefficients=active,active_threshold=threshold,
        df_seconds=df_seconds,hessian_bytes=hessian.nbytes,penalized_bytes=penalized.nbytes,
        coefficient_bytes=coefficient.nbytes,sample_design_bytes=sum(s['designs'][channel][1].nbytes for s in samples))


def evaluate_block(base, blocks, spectrograph, position):
    """``base`` state with R and Z replaced by the block surfaces; B keeps the base (survey or median) LSF."""
    state=deepcopy(base)
    for channel in FIT_CHANNELS:
        block=blocks[(int(spectrograph),channel)]
        basis=BSpline.design_matrix([position],block['slit_knots'],3).toarray()[0]
        state.coefficients[channel]=np.einsum('kqs,s->kq',block['coef'],basis)
        state.knot_vectors[channel]=np.array(block['wave_knots'],float)
        state.degrees[channel]=int(block['wave_degree'])
        adopted=block['qa'].get('adopted',True)
        state.metrics[channel]=dict(status=block['qa']['status'],model='sky_block_centered' if block['qa'].get('centered',True) else 'sky_block',
            quality='measured' if adopted else 'reference_fallback',converged=adopted,
            reason=block['qa'].get('fallback_reason',''))
    count = max(c.shape[0] for c in state.coefficients.values())
    state.tap_offsets = np.arange(-(count//2), count//2+1)
    return state
