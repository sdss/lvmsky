# Wavelength-resolved residual atlas (physics-space failure map).
#
# Sample n_sample test rows, reconstruct both the true and the ML-predicted
# sci-arm spectra, and aggregate the flux-space residual per pixel over the
# sample.  The atlas shows *where in wavelength* the model fails: a coherent
# blue tilt = zodi mis-calibration; sharp bumps at OH band centres = OH group
# mis-scaling; ROLO mineral bands = moon mis-projection.
#
# What to look for:
#   * flat, near-zero median residual with narrow 68% band -> healthy across
#     the whole range.
#   * coherent blue lift -> zodi head is under-predicting; check the
#     naive_baseline() close_zodi regime.
#   * OH band spikes -> the OH head is biased at specific transitions;
#     probably worth adding a per-group loss floor.
#   * moon band residuals (700-1000 nm mineral features) -> moon head is
#     picking the wrong knot amplitudes.
#
# NB: this reuses the SAME reconstruction path full_spectrum_batch_rmse uses,
# on the every10 corpus files (fast: ~1 s per row after basis build).
import time as _time
from pathlib import Path

required = ['filtered_triplet', 'mlp_artifacts', 'predict_sci_coefficients_default',
            'context_cols', 'build_triplet_coef_dataset',
            'reconstruct_component_spectra', 'reconstruct_with_lsf',
            'load_lsf_state_if_available', 'load_o2_vector_if_available',
            '_infer_base_dir_for_reconstruction', 'FACTOR']
_missing = [k for k in required if k not in globals()]
if _missing:
    raise RuntimeError('Missing kernel state: ' + ', '.join(_missing))

n_sample_atlas = 200
rng_seed_atlas = 42

_e10_stem = f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10'
_e10_suffix = _DECOMP_SUFFIX
EVERY10_INPUT = f'{_e10_stem}.fits'
EVERY10_NEAR  = f'{_e10_stem}_decomp_sky1{_e10_suffix}.fits'
EVERY10_FAR   = f'{_e10_stem}_decomp_sky2{_e10_suffix}.fits'
EVERY10_SCI   = f'{_e10_stem}_decomp_sci{_e10_suffix}.fits'

# Build aligned e10 triplet (matches full_spectrum_batch_rmse).
e10_triplet = build_triplet_coef_dataset(
    input_fits_path=EVERY10_INPUT,
    sky_near_decomp_fits_path=EVERY10_NEAR,
    sky_far_decomp_fits_path=EVERY10_FAR,
    sci_decomp_fits_path=EVERY10_SCI,
    context_columns=context_cols,
    return_chi2=False,
)
if '_augment_triplet_with_ecliptic' in globals():
    _augment_triplet_with_ecliptic(e10_triplet, meta_fits_path=EVERY10_INPUT)
if '_augment_triplet_with_physics_priors' in globals():
    _augment_triplet_with_physics_priors(e10_triplet)

_e10_n0 = int(e10_triplet['n_rows'])
if _e10_n0 == 0:
    raise RuntimeError('No aligned rows in e10_triplet')

# Apply LMC/SMC field exclusion (do NOT apply the kappa filters -- the atlas
# wants a realistic sample, tail included).
_keep = np.ones(_e10_n0, dtype=bool)
if 'sci_ra' in e10_triplet and 'sci_dec' in e10_triplet:
    _sci_ra = np.asarray(e10_triplet['sci_ra'], dtype=np.float64)
    _sci_dec = np.asarray(e10_triplet['sci_dec'], dtype=np.float64)
    for _region in (LMC_EXCLUSION, SMC_EXCLUSION):
        _sep = _angular_separation_deg_vec(
            _sci_ra, _sci_dec, _region['ra_deg'], _region['dec_deg'])
        _keep &= ~(np.isfinite(_sep) & (_sep <= float(_region['radius_deg'])))

_avail = np.where(_keep)[0]
_rng = np.random.default_rng(rng_seed_atlas)
_n_pick = int(min(n_sample_atlas, _avail.size))
sel_pos = _rng.choice(_avail, size=_n_pick, replace=False)
sel_rows = np.asarray(e10_triplet['row_index'], dtype=np.int64)[sel_pos]
print(f'  atlas sample: n = {_n_pick} of {_avail.size} available every10 rows')

# ML predictions for these rows.
coef_sci_pred_atlas = predict_sci_coefficients_default(
    mlp_artifacts,
    coef_near_phys=np.asarray(e10_triplet['coef_near'][sel_pos], dtype=np.float32),
    coef_far_phys=np.asarray( e10_triplet['coef_far'][sel_pos],  dtype=np.float32),
    ctx_near_phys=np.asarray( e10_triplet['ctx_near'][sel_pos],  dtype=np.float32),
    ctx_far_phys=np.asarray(  e10_triplet['ctx_far'][sel_pos],   dtype=np.float32),
    ctx_sci_phys=np.asarray(  e10_triplet['ctx_sci'][sel_pos],   dtype=np.float32),
).astype(np.float32)
coef_sci_true_atlas = np.asarray(e10_triplet['coef_sci'][sel_pos], dtype=np.float32)

# Load wavelength grid + LSF ref from the input FITS (matches the recon cell).
with fits.open(EVERY10_INPUT) as hdul:
    wave_arr = np.asarray(hdul['WAVE'].data, dtype=np.float64)
    lsf_sci_arr = np.asarray(hdul['LSF_SCI'].data, dtype=np.float64)

_wave_ref_recon = wave_arr if wave_arr.ndim == 1 else wave_arr[int(sel_rows[0])]
base_dir_guess = _infer_base_dir_for_reconstruction()

_lsf_model = SkyDecompLSFSurfaceIterative(
    _wave_ref_recon, lsf_sigma=1.0, n_spline_knots=N_MOON_KNOTS,
    base_dir=base_dir_guess,
    palace_oh_suffix='_joint_v2_updated',
    palace_diffuse_suffix='_joint_native_adam_invsky_p2_10000iter',
    split_zodi=SPLIT_ZODI, n_zodi_spline_knots=N_ZODI_KNOTS,
)


def _precache_decomp_state(decomp_path):
    _state = {'path': Path(decomp_path), 'has_lsf': False, 'o2_cube': None}
    if not _state['path'].exists():
        return _state
    try:
        with fits.open(str(_state['path']), memmap=False) as _hdul_dec:
            _names = {h.name for h in _hdul_dec}
            _state['has_lsf'] = all(_e in _names for _e in ('LSF_COEF', 'LSF_KNOTS', 'LSF_META'))
            if 'VECTOR_O2' in _names:
                _data = np.asarray(_hdul_dec['VECTOR_O2'].data, dtype=np.float64)
                if _data.ndim == 2:
                    _state['o2_cube'] = _data
    except (KeyError, IndexError, ValueError) as _exc:
        print(f'  precache failed for {_state["path"].name}: {_exc}')
    return _state


_state_sci = _precache_decomp_state(EVERY10_SCI)


def _reconstruct_sci_total(coef_row, row_idx, lsf_sigma_fallback):
    _lsf_state = None
    if _state_sci['has_lsf']:
        try:
            _lsf_state = load_lsf_surface_state(str(_state_sci['path']), int(row_idx))
        except (KeyError, IndexError, ValueError):
            _lsf_state = None
    _o2 = None
    if _state_sci['o2_cube'] is not None and int(row_idx) < _state_sci['o2_cube'].shape[0]:
        _o2_row = _state_sci['o2_cube'][int(row_idx)]
        if np.isfinite(_o2_row).any() and float(np.nansum(np.abs(_o2_row))) > 0.0:
            _o2 = _o2_row
    if isinstance(_lsf_state, LSFSurfaceState):
        _lsf_model._set_lsf_state(_lsf_state)
        _mats = _lsf_model._assemble_refined_matrices()
        if _o2 is not None:
            _mats['o2'] = np.asarray(_o2, float).ravel()[None, :]
        _comps = _lsf_model._components_from_coef(np.asarray(coef_row, float).ravel(), _mats)
        _total = (_comps['oh'] + _comps['moon'] + _comps.get('zodi', 0)
                  + _comps['diffuse'] + _comps['atom'] + _comps['orc'] + _comps['o2'])
        return np.asarray(_total, dtype=np.float64)
    # Rare fallback.
    _comps = reconstruct_with_lsf(
        wave=_lsf_model.wave, coef=coef_row, lsf=lsf_sigma_fallback,
        n_spline_knots=N_MOON_KNOTS, base_dir=base_dir_guess, o2_vector=_o2,
        split_zodi=SPLIT_ZODI, n_zodi_spline_knots=N_ZODI_KNOTS)
    _total = (_comps['oh'] + _comps['moon'] + _comps.get('zodi', 0)
              + _comps['diffuse'] + _comps['atom'] + _comps['orc'] + _comps['o2'])
    return np.asarray(_total, dtype=np.float64)


# Reconstruct N rows: shape (n_pick, n_pix).
_n_pix = _wave_ref_recon.size
_resid = np.full((_n_pick, _n_pix), np.nan, dtype=np.float64)
_truth = np.full((_n_pick, _n_pix), np.nan, dtype=np.float64)

_t0 = _time.perf_counter()
for _i, _rr in enumerate(sel_rows):
    _lsf_row = lsf_sci_arr if lsf_sci_arr.ndim == 1 else lsf_sci_arr[int(_rr)]
    _lsf_sigma_fb = _lsf_row / 2.35
    _y_true = _reconstruct_sci_total(coef_sci_true_atlas[_i], _rr, _lsf_sigma_fb) / FACTOR
    _y_pred = _reconstruct_sci_total(coef_sci_pred_atlas[_i], _rr, _lsf_sigma_fb) / FACTOR
    _resid[_i, :] = _y_pred - _y_true
    _truth[_i, :] = _y_true
    if _i and _i % 50 == 0:
        print(f'  reconstructed {_i}/{_n_pick} rows ({_time.perf_counter() - _t0:.1f}s)')
print(f'  atlas reconstruction: {_time.perf_counter() - _t0:.1f}s '
      f'({_n_pick} rows x 2 spectra)')

# Aggregate residuals per pixel.
_med   = np.nanmedian(_resid, axis=0)
_p16   = np.nanpercentile(_resid, 16, axis=0)
_p84   = np.nanpercentile(_resid, 84, axis=0)
_iqr   = _p84 - _p16
_truth_med = np.nanmedian(_truth, axis=0)

# Normalise residual by the local truth median so a common frac_residual axis
# is meaningful across the range.
_denom = np.where(np.abs(_truth_med) > 1e-30, _truth_med, np.nan)
_med_frac = _med   / _denom
_p16_frac = _p16   / _denom
_p84_frac = _p84   / _denom

# Plot: three panels (absolute residual with band, fractional, and truth flux for context).
_fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
    subplot_titles=('Median pred - true (absolute) with 68% band',
                    'Median pred - true / truth median (fractional)',
                    'Median truth flux (reference)'),
    vertical_spacing=0.06)

_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_p84,
    mode='lines', line=dict(color='rgba(0,120,220,0.0)'),
    showlegend=False), row=1, col=1)
_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_p16,
    mode='lines', line=dict(color='rgba(0,120,220,0.0)'),
    fill='tonexty', fillcolor='rgba(0,120,220,0.20)',
    showlegend=True, name='68% band (p16..p84)'), row=1, col=1)
_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_med,
    mode='lines', line=dict(color='steelblue', width=1.5),
    name='median residual'), row=1, col=1)
_fig.add_shape(type='line', x0=float(_wave_ref_recon[0]),
               x1=float(_wave_ref_recon[-1]), y0=0.0, y1=0.0,
               line=dict(color='black', dash='dot', width=1), row=1, col=1)

_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_p84_frac,
    mode='lines', line=dict(color='rgba(180,50,50,0.0)'),
    showlegend=False), row=2, col=1)
_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_p16_frac,
    mode='lines', line=dict(color='rgba(180,50,50,0.0)'),
    fill='tonexty', fillcolor='rgba(180,50,50,0.20)',
    showlegend=True, name='fractional 68% band'), row=2, col=1)
_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_med_frac,
    mode='lines', line=dict(color='firebrick', width=1.5),
    name='median fractional residual'), row=2, col=1)
_fig.add_shape(type='line', x0=float(_wave_ref_recon[0]),
               x1=float(_wave_ref_recon[-1]), y0=0.0, y1=0.0,
               line=dict(color='black', dash='dot', width=1), row=2, col=1)

_fig.add_trace(go.Scatter(x=_wave_ref_recon, y=_truth_med,
    mode='lines', line=dict(color='gray', width=1),
    showlegend=False), row=3, col=1)

_fig.update_yaxes(title_text='pred - true (flux units)', row=1, col=1)
_fig.update_yaxes(title_text='(pred - true) / truth_median', row=2, col=1,
                  range=[-0.5, 0.5])
_fig.update_yaxes(title_text='truth median flux', row=3, col=1, type='log')
_fig.update_xaxes(title_text='wavelength [Å]', row=3, col=1)
_fig.update_layout(height=850, width=1200,
                   title=f'Wavelength residual atlas (n = {_n_pick} rows, '
                         f'test-representative every10 sample)')
_fig.show()

# Compact summary: RMS of the median residual by band.
_bands = [('blue    3600-5500', 3600.0, 5500.0),
          ('mid     5500-7500', 5500.0, 7500.0),
          ('NIR     7500-9800', 7500.0, 9800.0)]
print()
print('Per-band RMS of the median residual and median fractional residual:')
for _name, _lo, _hi in _bands:
    _mask_lam = (_wave_ref_recon >= _lo) & (_wave_ref_recon < _hi)
    if not _mask_lam.any():
        continue
    _rms_abs = float(np.sqrt(np.nanmean(_med[_mask_lam] ** 2)))
    _rms_frac = float(np.sqrt(np.nanmean(_med_frac[_mask_lam] ** 2)))
    _mean_bias_frac = float(np.nanmean(_med_frac[_mask_lam]))
    print(f'  {_name:<20s} RMS|abs| = {_rms_abs:.3g}   '
          f'RMS|frac| = {_rms_frac:.3g}   mean_bias_frac = {_mean_bias_frac:+.3g}')
