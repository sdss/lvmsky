"""One convex LSF fit per spectrograph/channel, with frozen fiber populations."""
from copy import deepcopy
import time
import clarabel
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import BSpline
from skysub.sky_decomp.lsf_spline2d import _offset_knots

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


def solve_normal(hessian, linear, shape, centers):
    """Enforce unit integral and zero centroid on every wavelength/slit control pair."""
    no, nw, ns = shape
    size = no*nw*ns
    identity = sp.eye(nw*ns, format='csc')
    equality = sp.vstack([sp.kron(np.ones((1,no)), identity), sp.kron(centers[None,:], identity)], format='csc')
    constraint = sp.vstack([equality, -sp.eye(size)], format='csc')
    rhs = np.r_[np.ones(nw*ns), np.zeros(nw*ns+size)]
    scale = max(float(np.trace(hessian)/size), 1e-30)
    settings = clarabel.DefaultSettings(); settings.verbose = False
    settings.tol_gap_abs = settings.tol_feas = settings.tol_gap_rel = 1e-9
    solution = clarabel.DefaultSolver(sp.triu(sp.csc_matrix(hessian/scale)).tocsc(), linear/scale,
        constraint, rhs, [clarabel.ZeroConeT(2*nw*ns), clarabel.NonnegativeConeT(size)], settings).solve()
    if str(solution.status) not in ('Solved', 'AlmostSolved'):
        raise RuntimeError(f'Block LSF QP: {solution.status}')
    coefficient = np.asarray(solution.x).reshape(shape)
    # Keep the equality solution intact; clipping and renormalizing alter the first moment.
    if (not np.isfinite(coefficient).all() or coefficient.min() < -1e-8
        or np.max(np.abs(coefficient.sum(axis=0)-1)) > 1e-7
        or np.max(np.abs(np.einsum('k,kqs->qs',centers,coefficient))) > 1e-7):
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


def fit_block(samples, state, channel, lower, upper, slit_count=4):
    """Accumulate small normal equations without a full multi-fiber design array."""
    started=time.perf_counter()
    no,nw = state.coefficients[channel].shape
    shape = (no,nw,slit_count); size = int(np.prod(shape))
    hessian = np.zeros((size,size)); linear = np.zeros(size)
    basis,knots = slit_basis([s['position'] for s in samples],lower,upper,slit_count)
    counts=[]; chi_constant=0
    for sample,weights in zip(samples,basis):
        take,design,fixed = sample['designs'][channel]
        y = sample['flux'][take]-fixed
        ivar = sample['ivar'][take]
        good = (sample['mask'][take]==0)&np.isfinite(y)&np.isfinite(ivar)&(ivar>0)
        # Neighboring-channel solar leakage is fixed in each affine design; use interior pixels.
        wave = sample['wave'][take]
        good &= (wave>=wave[0]+5)&(wave<=wave[-1]-5)
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
    coefficient,status = solve_normal(hessian,linear,shape,centers)
    solve_seconds=time.perf_counter()-solve_started
    vector=coefficient.ravel()
    chi2=float(vector@hessian@vector+2*linear@vector+chi_constant)
    return coefficient,knots,dict(status=status,n_fibers=len(samples),n_pixels=sum(counts),
        normal_seconds=normal_seconds,solve_seconds=solve_seconds,chi2=chi2,chi2_per_pixel=chi2/sum(counts),slit_rank=int(np.linalg.matrix_rank(basis)),
        selected_min=float(min(s['position'] for s in samples)),selected_max=float(max(s['position'] for s in samples)),
        normalization_error=float(np.max(abs(coefficient.sum(axis=0)-1))),
        centroid_error=float(np.max(abs(np.einsum('k,kqs->qs',centers,coefficient)))),coefficient_min=float(coefficient.min()))


def evaluate_block(base, blocks, spectrograph, position):
    """``base`` state with R and Z replaced by the block surfaces; B keeps the base (survey or median) LSF."""
    state=deepcopy(base)
    for channel in FIT_CHANNELS:
        block=blocks[(int(spectrograph),channel)]
        basis=BSpline.design_matrix([position],block['slit_knots'],3).toarray()[0]
        state.coefficients[channel]=np.einsum('kqs,s->kq',block['coef'],basis)
        state.knot_vectors[channel]=np.array(block['wave_knots'],float)
        state.degrees[channel]=int(block['wave_degree'])
        state.metrics[channel]=dict(status=block['qa']['status'],model='sky_block_centered',quality='measured')
    return state
