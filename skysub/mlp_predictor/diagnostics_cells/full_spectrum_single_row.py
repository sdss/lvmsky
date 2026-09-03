# Full-spectrum reconstruction test for a single requested row using default coefficients
import plotly.graph_objects as go

_e10_stem = f"{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10"
_e10_suffix = _DECOMP_SUFFIX  # inherited from cell 6
EVERY10_INPUT = f"{_e10_stem}.fits"
EVERY10_NEAR = f"{_e10_stem}_decomp_sky1{_e10_suffix}.fits"
EVERY10_FAR = f"{_e10_stem}_decomp_sky2{_e10_suffix}.fits"
EVERY10_SCI = f"{_e10_stem}_decomp_sci{_e10_suffix}.fits"

# Set the row to reconstruct and inspect.
# REQUESTED_ROW = 612
# REQUESTED_ROW = 500
REQUESTED_ROW = 978
# REQUESTED_ROW = 1481
# REQUESTED_ROW = 741 #OK
# REQUESTED_ROW = 742
# REQUESTED_ROW = 830

# Overlay the frozen physical Moon/Zodi model on the three flux panels
# (see 5b below).  Costs ~0.25 s per arm; set False to skip it entirely.
SHOW_MOON_ZODI_MODEL = True

required = [
    "mlp_artifacts",
    "predict_sci_coefficients_default",
    "context_cols",
    "build_triplet_coef_dataset",
    "reconstruct_component_spectra",
    "reconstruct_with_lsf",
    "load_lsf_state_if_available",
    "load_o2_vector_if_available",
    "_infer_base_dir_for_reconstruction",
    "_moon_bs_indices_from_names",
    "_row_spline_roughness",
]
missing = [k for k in required if k not in globals()]
if missing:
    raise RuntimeError("Run the training + residual-correction cells first. Missing: " + ", ".join(missing))

# 1) Load coefficients/context from every10 decomposition products.
e10_triplet = build_triplet_coef_dataset(
    input_fits_path=EVERY10_INPUT,
    sky_near_decomp_fits_path=EVERY10_NEAR,
    sky_far_decomp_fits_path=EVERY10_FAR,
    sci_decomp_fits_path=EVERY10_SCI,
    context_columns=context_cols,
    return_chi2=False,
)
# ECLIPTIC-CTX-V1: match training-time ctx layout on the e10 triplet.
if '_augment_triplet_with_ecliptic' in globals():
    _augment_triplet_with_ecliptic(e10_triplet, meta_fits_path=EVERY10_INPUT)
# PHYSICS-PRIORS-CTX-V1: same augment on the e10 triplet.
if '_augment_triplet_with_physics_priors' in globals():
    _augment_triplet_with_physics_priors(e10_triplet)
n_e10 = int(e10_triplet["n_rows"])
row_index_e10 = np.asarray(e10_triplet["row_index"], dtype=np.int64)
coef_names_e10 = [str(n) for n in e10_triplet["coef_names"]]
moon_idx = _moon_bs_indices_from_names(coef_names_e10)
zodi_idx = np.array(
    [i for i, n in enumerate(coef_names_e10) if n.startswith('Zodi_bs')],
    dtype=np.int64,
)

if moon_idx.size < 3:
    raise RuntimeError("Expected at least 3 Moon_bs coefficients for spline diagnostics")

# 2) Load observed spectra, wavelength grid, and LSF from every10 input.
with fits.open(EVERY10_INPUT) as hdul:
    for ext in ("FLUX_SKY_NEAR", "FLUX_SKY_FAR", "FLUX_SCI", "WAVE", "LSF_SCI"):
        if ext not in hdul:
            raise KeyError(f"Missing extension {ext} in {EVERY10_INPUT}")

    wave_arr = np.asarray(hdul["WAVE"].data, dtype=np.float64)
    flux_near_all = np.asarray(hdul["FLUX_SKY_NEAR"].data, dtype=np.float64)
    flux_far_all = np.asarray(hdul["FLUX_SKY_FAR"].data, dtype=np.float64)
    flux_sci_true_all = np.asarray(hdul["FLUX_SCI"].data, dtype=np.float64)
    lsf_sci_arr = np.asarray(hdul["LSF_SCI"].data, dtype=np.float64)

n_spec, n_wave = flux_sci_true_all.shape
if row_index_e10.size != n_e10:
    raise ValueError(f"Triplet row_index length mismatch: {row_index_e10.size} vs n_rows={n_e10}")
if np.any(row_index_e10 < 0) or np.any(row_index_e10 >= n_spec):
    raise ValueError(
        f"Triplet row_index contains values outside [0, {n_spec - 1}] for {EVERY10_INPUT}"
    )

idx_row = int(REQUESTED_ROW)
triplet_pos = np.flatnonzero(row_index_e10 == idx_row)
if triplet_pos.size == 0:
    raise IndexError(
        f"REQUESTED_ROW={idx_row} is not available in aligned triplet rows. "
        f"Choose one of e10_triplet['row_index'] (size={row_index_e10.size})."
    )
triplet_pos = int(triplet_pos[0])

# Normalize WAVE/LSF arrays to per-row vectors, then select requested row.
wave_row = wave_arr if wave_arr.ndim == 1 else wave_arr[idx_row]
lsf_row = lsf_sci_arr if lsf_sci_arr.ndim == 1 else lsf_sci_arr[idx_row]

flux_near_row = flux_near_all[idx_row]
flux_far_row = flux_far_all[idx_row]
flux_sci_true_row = flux_sci_true_all[idx_row]

# 3) Predict SCI coefficients for the requested row using global default path.
coef_pred_row_batch = predict_sci_coefficients_default(
    mlp_artifacts,
    coef_near_phys=e10_triplet["coef_near"][triplet_pos: triplet_pos + 1],
    coef_far_phys=e10_triplet["coef_far"][triplet_pos: triplet_pos + 1],
    ctx_near_phys=e10_triplet["ctx_near"][triplet_pos: triplet_pos + 1],
    ctx_far_phys=e10_triplet["ctx_far"][triplet_pos: triplet_pos + 1],
    ctx_sci_phys=e10_triplet["ctx_sci"][triplet_pos: triplet_pos + 1],
)
coef_pred_row = np.asarray(coef_pred_row_batch[0], dtype=np.float64)

# 3a) Context values for the near-sky, far-sky and science pointings at this row.
#     Cyclic sin/cos pairs are folded back to 0-360 degree axes for readability.
_ctx_names_e10 = list(e10_triplet["ctx_names"])
_ctx_row_stack = np.stack([
    e10_triplet["ctx_near"][triplet_pos],
    e10_triplet["ctx_far"][triplet_pos],
    e10_triplet["ctx_sci"][triplet_pos],
], axis=0)
_ctx_disp_names, _ctx_disp_stack = _decode_cyclic_context(_ctx_names_e10, _ctx_row_stack)
_ctx_row_df = pd.DataFrame({
    "feature": _ctx_disp_names,
    "sky_near": _ctx_disp_stack[0],
    "sky_far": _ctx_disp_stack[1],
    "science": _ctx_disp_stack[2],
})
print(f"Context values at row {idx_row} (sin/cos pairs decoded to degrees):")
print(_ctx_row_df.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
print()

# 3b) Moon_bs coefficient diagnostics (global prior already applied).
coef_near_row = np.asarray(e10_triplet["coef_near"][triplet_pos], dtype=np.float64)
coef_far_row = np.asarray(e10_triplet["coef_far"][triplet_pos], dtype=np.float64)
coef_sci_row = np.asarray(e10_triplet["coef_sci"][triplet_pos], dtype=np.float64)

# Per-row coef_err arrays (LSF-propagated sigma feeds off these when the
# decomposition FITS lacks a FLUX_SIGMA_TOTAL HDU).  Missing HDU -> None,
# and the downstream WRMSE degrades to the median-floor path.
def _row_coef_err(triplet_key):
    _arr = e10_triplet.get(triplet_key)
    if _arr is None:
        return None
    _row = np.asarray(_arr[triplet_pos], dtype=np.float64)
    return _row if np.any(np.isfinite(_row)) else None

coef_err_near_row = _row_coef_err("coef_err_near")
coef_err_far_row  = _row_coef_err("coef_err_far")
coef_err_sci_row  = _row_coef_err("coef_err_sci")
moon_pred = coef_pred_row[moon_idx]
zodi_pred = coef_pred_row[zodi_idx] if zodi_idx.size else np.zeros(0)
zodi_near = coef_near_row[zodi_idx] if zodi_idx.size else np.zeros(0)
zodi_far  = coef_far_row [zodi_idx] if zodi_idx.size else np.zeros(0)
zodi_true = coef_sci_row [zodi_idx] if zodi_idx.size else np.zeros(0)
moon_near = coef_near_row[moon_idx]
moon_far = coef_far_row[moon_idx]
moon_true = coef_sci_row[moon_idx]

# 4) Reconstruct this row for SCI prediction and near/far self-consistency checks.
base_dir_guess = _infer_base_dir_for_reconstruction()

# Prefer the fitted wavelength-dependent LSF surface from each decomp file's
# LSF_COEF/LSF_KNOTS/LSF_META extensions (written by sky_decomp.lsf_surface_iterative);
# fall back to the input FITS's Gaussian LSF_SCI if that isn't present.
_lsf_state_near = load_lsf_state_if_available(EVERY10_NEAR, idx_row)
_lsf_state_far  = load_lsf_state_if_available(EVERY10_FAR,  idx_row)
_lsf_state_sci  = load_lsf_state_if_available(EVERY10_SCI,  idx_row)
_lsf_sigma_fallback = lsf_row / 2.35
print(f"  LSF source per arm (row {idx_row}): "
      f"near={'surface' if _lsf_state_near is not None else 'gaussian (LSF_SCI)'}, "
      f"far={'surface' if _lsf_state_far is not None else 'gaussian (LSF_SCI)'}, "
      f"sci={'surface' if _lsf_state_sci is not None else 'gaussian (LSF_SCI)'}")

# Per-arm unit-integrated O2 templates from the decomposition FITS
# (VECTOR_O2 extension). If absent (older decomp), the O2 basis stays at
# zero -- matching pre-2026-08-10 behaviour. For the predicted-sci
# reconstruction, use the sci-arm template so the shape is anchored on the
# same layer temperature the science pointing had; the amplitude comes
# from the predicted coef['O2_b01'].
_o2_vec_near = load_o2_vector_if_available(EVERY10_NEAR, idx_row)
_o2_vec_far  = load_o2_vector_if_available(EVERY10_FAR,  idx_row)
_o2_vec_sci  = load_o2_vector_if_available(EVERY10_SCI,  idx_row)
print(f"  O2 template per arm (row {idx_row}): "
      f"near={'VECTOR_O2' if _o2_vec_near is not None else 'zero'}, "
      f"far={'VECTOR_O2' if _o2_vec_far is not None else 'zero'}, "
      f"sci={'VECTOR_O2' if _o2_vec_sci is not None else 'zero'}")

comps_sci = reconstruct_with_lsf(
    wave=wave_row,
    coef=coef_pred_row,
    lsf=_lsf_state_sci if _lsf_state_sci is not None else _lsf_sigma_fallback,
    n_spline_knots=N_MOON_KNOTS,
    n_zodi_spline_knots=N_ZODI_KNOTS,
    base_dir=base_dir_guess,
    o2_vector=_o2_vec_sci,
)
comps_near_from_near = reconstruct_with_lsf(
    wave=wave_row,
    coef=coef_near_row,
    lsf=_lsf_state_near if _lsf_state_near is not None else _lsf_sigma_fallback,
    n_spline_knots=N_MOON_KNOTS,
    n_zodi_spline_knots=N_ZODI_KNOTS,
    base_dir=base_dir_guess,
    o2_vector=_o2_vec_near,
    coef_err=coef_err_near_row,
)
comps_far_from_far = reconstruct_with_lsf(
    wave=wave_row,
    coef=coef_far_row,
    lsf=_lsf_state_far if _lsf_state_far is not None else _lsf_sigma_fallback,
    n_spline_knots=N_MOON_KNOTS,
    n_zodi_spline_knots=N_ZODI_KNOTS,
    base_dir=base_dir_guess,
    o2_vector=_o2_vec_far,
    coef_err=coef_err_far_row,
)
# Reconstruction of the observed sci spectrum from the fitted sci coefs, so panel 3
# can separate the sky-decomposition fit residual (obs vs recon-from-sci-coef)
# from our transfer-model error (recon-from-pred vs recon-from-sci-coef).
comps_sci_true = reconstruct_with_lsf(
    wave=wave_row,
    coef=coef_sci_row,
    lsf=_lsf_state_sci if _lsf_state_sci is not None else _lsf_sigma_fallback,
    n_spline_knots=N_MOON_KNOTS,
    n_zodi_spline_knots=N_ZODI_KNOTS,
    base_dir=base_dir_guess,
    o2_vector=_o2_vec_sci,
    coef_err=coef_err_sci_row,
)

flux_sci_pred_row = np.asarray(comps_sci["total"], dtype=np.float64) / FACTOR
flux_sci_true_recon_row = np.asarray(comps_sci_true["total"], dtype=np.float64) / FACTOR
flux_near_recon_row = np.asarray(comps_near_from_near["total"], dtype=np.float64) / FACTOR
flux_far_recon_row = np.asarray(comps_far_from_far["total"], dtype=np.float64) / FACTOR

# 5) Single-row metrics.
resid_row = flux_sci_pred_row - flux_sci_true_row
rmse_row = float(np.sqrt(np.mean(resid_row ** 2)))
rmse_row_display = float(rmse_row * FACTOR)
mae_row = float(np.mean(np.abs(resid_row)))
rel_resid_row = resid_row / np.where(flux_sci_true_row != 0, flux_sci_true_row, np.nan)

rmse_near_recon = float(np.sqrt(np.mean((flux_near_recon_row - flux_near_row) ** 2)))
rmse_far_recon = float(np.sqrt(np.mean((flux_far_recon_row - flux_far_row) ** 2)))

# Pixel-space WRMSE for the same three arms.  Sigma source order:
#   1. FLUX_SIGMA_TOTAL HDU in the decomposition FITS (future LSF-aware
#      propagator output);
#   2. sigma_total returned by `reconstruct_component_spectra(coef_err=...)`
#      when the corresponding COEF_ERR array is available on this row;
#   3. None -> pixel_wrmse_per_row degrades to a per-row median-floor RMSE.
def _sigma_for_single_row(fits_path, hdu_row_idx, comps_dict, comps_true_dict=None):
    """Prefer FITS HDU sigma; else use comps sigma_total (from coef_err)."""
    _fits_sigma = load_pixel_sigma_if_available(fits_path,
                                                row_indices=[int(hdu_row_idx)])
    if _fits_sigma is not None:
        return np.asarray(_fits_sigma[0], dtype=np.float64) / FACTOR
    # Fallback: use the sigma_total attached to the reconstructed comps dict.
    _sig = comps_dict.get("sigma_total") if isinstance(comps_dict, dict) else None
    if _sig is not None:
        return np.asarray(_sig, dtype=np.float64) / FACTOR
    # As a last resort, try the "true" recomposition's sigma (only meaningful
    # for the sci arm where the true decomposition sigma is more directly
    # comparable to the observation noise floor).
    _sig_alt = (comps_true_dict.get("sigma_total")
                if isinstance(comps_true_dict, dict) else None)
    if _sig_alt is not None:
        return np.asarray(_sig_alt, dtype=np.float64) / FACTOR
    return None

_sig_near_pix = _sigma_for_single_row(EVERY10_NEAR, idx_row, comps_near_from_near)
_sig_far_pix  = _sigma_for_single_row(EVERY10_FAR,  idx_row, comps_far_from_far)
_sig_sci_pix  = _sigma_for_single_row(EVERY10_SCI,  idx_row, comps_sci, comps_sci_true)

wrmse_near_recon = float(pixel_wrmse_per_row(
    flux_near_recon_row, flux_near_row, _sig_near_pix)[0])
wrmse_far_recon  = float(pixel_wrmse_per_row(
    flux_far_recon_row,  flux_far_row,  _sig_far_pix)[0])
wrmse_row_pix    = float(pixel_wrmse_per_row(
    flux_sci_pred_row,   flux_sci_true_row, _sig_sci_pix)[0])

print("Single-row reconstruction summary (every10, default coefficients)")
print(f"  row index (input file) = {idx_row}")
print(f"  row index (triplet pos) = {triplet_pos}")
print(f"  n_wave     = {n_wave}")
print("  predictor  = deep group-head MLP")
print(
    "  Moon_bs roughness: "
    f"pred={_row_spline_roughness(moon_pred):.4g}, "
    f"near={_row_spline_roughness(moon_near):.4g}, "
    f"far={_row_spline_roughness(moon_far):.4g}, "
    f"sci_true={_row_spline_roughness(moon_true):.4g}"
)
print(f"  near self-recon pRMSE  = {rmse_near_recon:.6g}")
print(f"  near self-recon pWRMSE = {wrmse_near_recon:.6g}")
print(f"  far  self-recon pRMSE  = {rmse_far_recon:.6g}")
print(f"  far  self-recon pWRMSE = {wrmse_far_recon:.6g}")
print(f"  sci row pRMSE          = {rmse_row:.6g}")
print(f"  sci row pWRMSE         = {wrmse_row_pix:.6g}")
print(f"  sci row pRMSE (x{FACTOR:.3g} display units) = {rmse_row_display:.6g}")
print(f"  sci row MAE            = {mae_row:.6g}")

# 5b) Training-density audit in regime space. This replaces the earlier
# coef_err_sci percentile audit -- which was mostly a brightness proxy --
# with a direct measure of whether this row sits in a sparsely-covered part
# of the training distribution on the physical failure axes. That is the
# real cause of tail errors on rows the decomposition itself fits fine.
_REGIME_AXES = ('moon_alt', 'moon_fli', 'abs_ecl_beta', 'airmass', 'sun_alt')
_REGIME_K = 10
_REGIME_SPARSE_PERCENTILE = 90.0

def _regime_axis_key(axis_name):
    return 'ecl_beta_deg' if axis_name == 'abs_ecl_beta' else axis_name

def _regime_axis_value(vec, axis_name):
    return np.abs(vec) if axis_name == 'abs_ecl_beta' else vec

_pop_ctx_ss = np.asarray(filtered_triplet['ctx_sci'], dtype=np.float64)
_pop_ctx_names_ss = list(filtered_triplet['ctx_names'])
_train_idx_ss = np.asarray(
    mlp_artifacts.get('train_idx', np.arange(_pop_ctx_ss.shape[0])), dtype=int)
_train_idx_ss = _train_idx_ss[(_train_idx_ss >= 0) & (_train_idx_ss < _pop_ctx_ss.shape[0])]

_axis_names_present = []
_train_cols = []
for _a in _REGIME_AXES:
    _key = _regime_axis_key(_a)
    if _key not in _pop_ctx_names_ss:
        continue
    _train_cols.append(_regime_axis_value(
        _pop_ctx_ss[_train_idx_ss, _pop_ctx_names_ss.index(_key)], _a))
    _axis_names_present.append(_a)

sci_row_regime_audit = {}
sci_row_sparse_regime_flag = False
if not _axis_names_present:
    print()
    print('[regime audit skipped: no regime axes present in ctx_names]')
else:
    _train_stack = np.stack(_train_cols, axis=1)
    _finite_train = np.all(np.isfinite(_train_stack), axis=1)
    _train_stack = _train_stack[_finite_train]
    _mu = _train_stack.mean(axis=0)
    _sd = _train_stack.std(axis=0)
    _sd = np.where(_sd > 0, _sd, 1.0)
    _train_z = (_train_stack - _mu) / _sd

    # Cache the train-side kNN work per (axes, n_train, k) across cell re-runs.
    if '_regime_train_cache' not in globals():
        globals()['_regime_train_cache'] = {}
    _cache_key = tuple(_axis_names_present) + (int(_train_z.shape[0]), _REGIME_K)
    if _cache_key not in _regime_train_cache:
        try:
            from scipy.spatial import cKDTree as _cKDTree
            _tree = _cKDTree(_train_z)
            _train_d, _ = _tree.query(_train_z, k=_REGIME_K + 1)
            _train_knn_dist = _train_d[:, _REGIME_K]
        except ImportError:
            _tree = None
            _dsq_tr = ((_train_z[:, None, :] - _train_z[None, :, :]) ** 2).sum(axis=2)
            _dsq_tr.sort(axis=1)
            _train_knn_dist = np.sqrt(_dsq_tr[:, _REGIME_K])
        _regime_train_cache[_cache_key] = dict(
            tree=_tree, mu=_mu, sd=_sd,
            train_z=_train_z,
            train_stack=_train_stack,
            train_p90=float(np.percentile(
                _train_knn_dist, _REGIME_SPARSE_PERCENTILE)),
            train_knn_dist=_train_knn_dist,
        )
    _cache = _regime_train_cache[_cache_key]

    # Build the query row from the sci arm of e10_triplet for this row.
    _q_ctx = np.asarray(e10_triplet['ctx_sci'][triplet_pos], dtype=np.float64)
    _q_names = list(e10_triplet['ctx_names'])
    _q_vals = []
    _per_axis_pct = {}
    for _a in _axis_names_present:
        _key = _regime_axis_key(_a)
        if _key not in _q_names:
            _q_vals.append(np.nan)
            continue
        _raw = float(_q_ctx[_q_names.index(_key)])
        _val = float(np.abs(_raw)) if _a == 'abs_ecl_beta' else _raw
        _q_vals.append(_val)
        _tr_axis = _cache['train_stack'][:, _axis_names_present.index(_a)]
        _per_axis_pct[_a] = dict(
            value=_val,
            pop_percentile=float(100.0 * (_tr_axis < _val).mean()),
        )
    _q_arr = np.asarray(_q_vals, dtype=np.float64)
    if not np.all(np.isfinite(_q_arr)):
        print()
        print('[regime audit skipped: query row has NaN on a regime axis]')
    else:
        _q_z = (_q_arr - _cache['mu']) / _cache['sd']
        if _cache['tree'] is not None:
            _d, _ = _cache['tree'].query(_q_z[None, :], k=_REGIME_K + 1)
            _q_knn_all = _d[0]
        else:
            _dsq = ((_cache['train_z'] - _q_z[None, :]) ** 2).sum(axis=1)
            _q_knn_all = np.sqrt(np.sort(_dsq))
        # A near-zero first distance means the query is a training-set twin.
        _shift = 1 if _q_knn_all[0] < 1e-6 else 0
        _q_knn = float(_q_knn_all[_REGIME_K - 1 + _shift])
        _train_p90 = _cache['train_p90']
        _density_ratio = _q_knn / _train_p90 if _train_p90 > 0 else np.inf
        sci_row_sparse_regime_flag = bool(_q_knn > _train_p90)
        sci_row_regime_audit = dict(
            axis_names=list(_axis_names_present),
            axis_values=list(_q_vals),
            per_axis=_per_axis_pct,
            knn_distance=_q_knn,
            train_p90=_train_p90,
            density_ratio=_density_ratio,
            k=_REGIME_K,
            is_sparse=sci_row_sparse_regime_flag,
        )

        print()
        print('Regime-space training-density audit '
              '(this row vs training corpus):')
        print(f"  axes = {', '.join(_axis_names_present)}   "
              f"(k={_REGIME_K} nearest, {int(_cache['train_z'].shape[0])} train rows)")
        print(f"  {'axis':<16}{'row_value':>14}{'train_pct':>14}")
        for _a in _axis_names_present:
            _info = _per_axis_pct[_a]
            _tag = ' *' if (_info['pop_percentile'] > 95.0
                            or _info['pop_percentile'] < 5.0) else ''
            print(f"  {_a:<16}{_info['value']:>14.3f}"
                  f"{_info['pop_percentile']:>13.1f}%{_tag}")
        print(f"  kNN dist (standardized)     = {_q_knn:.3f}")
        print(f"  train p{_REGIME_SPARSE_PERCENTILE:.0f} threshold          "
              f"= {_train_p90:.3f}")
        print(f"  density ratio               = {_density_ratio:.2f}x")
        if sci_row_sparse_regime_flag:
            print(f'  [SPARSE_REGIME] row is beyond the p{_REGIME_SPARSE_PERCENTILE:.0f} '
                  f'training-density envelope;')
            print('                  the model is extrapolating in this regime; '
                  'residuals expected.')
        else:
            print(f'  [ok] row is inside the p{_REGIME_SPARSE_PERCENTILE:.0f} '
                  'training-density envelope.')

# Compact flag prefix + regime density subline reused in every plot title
# below, so the row's coverage status is always visible on the figure itself.
_flag_prefix = f"[{'SPARSE_REGIME' if sci_row_sparse_regime_flag else 'ok'}]"
if sci_row_regime_audit:
    _axis_parts = []
    for _a in sci_row_regime_audit['axis_names']:
        _info = sci_row_regime_audit['per_axis'][_a]
        _mark = '!' if (_info['pop_percentile'] > 95.0
                        or _info['pop_percentile'] < 5.0) else ''
        _axis_parts.append(f"{_a}={_info['value']:.1f}({_info['pop_percentile']:.0f}%{_mark})")
    _pct_subline = (f"regime kNN(k={sci_row_regime_audit['k']})"
                    f"={sci_row_regime_audit['knn_distance']:.2f} "
                    f"vs train_p{_REGIME_SPARSE_PERCENTILE:.0f}="
                    f"{sci_row_regime_audit['train_p90']:.2f} "
                    f"({sci_row_regime_audit['density_ratio']:.2f}x) · "
                    + ' · '.join(_axis_parts))
else:
    _pct_subline = ''


# 5b) Frozen physical Moon + Zodi model overlay (sky_decomp/moon_zodi_model.py).
#     This is independent of both the decomposition and the ML: it predicts the
#     scattered-moonlight and zodiacal continua from exposure-midpoint ephemeris
#     geometry alone (ROLO albedo x solar SED x Rayleigh/HG scattering for the
#     moon, Leinert B500 for the zodi).  Overlaying it on the three flux panels
#     gives a physical reference for the amplitudes the QP assigned to Moon_bs /
#     Zodi_bs, and for what the ML predicted at the sci pointing.  Each arm uses
#     its own pointing and its own LSF.
import warnings as _warnings


class _MzOverlayDisabled(Exception):
    """Internal sentinel: overlay switched off, not a failure."""


_mz_overlay = {}
if not SHOW_MOON_ZODI_MODEL:
    print("  moon/zodi physical model overlay: disabled "
          "(SHOW_MOON_ZODI_MODEL = False)")
try:
    if not SHOW_MOON_ZODI_MODEL:
        raise _MzOverlayDisabled
    from sky_decomp.moon_zodi_model import (
        MoonZodiInvalidObservationError,
        MoonZodiObservation,
        MoonZodiPhysicalModel,
    )

    with fits.open(EVERY10_INPUT) as _hdul_mz:
        _mz_meta = _hdul_mz["META"].data[idx_row]
        _mz_meta_names = set(_hdul_mz["META"].columns.names or ())
        _mz_lsf = {}
        for _arm, _ext in (("near", "LSF_SKY_NEAR"),
                           ("far", "LSF_SKY_FAR"),
                           ("sci", "LSF_SCI")):
            _a = (np.asarray(_hdul_mz[_ext].data, dtype=np.float64)
                  if _ext in _hdul_mz else lsf_sci_arr)
            _mz_lsf[_arm] = np.asarray(_a if _a.ndim == 1 else _a[idx_row],
                                       dtype=np.float64)

    # Exposure length: prefer a metadata column, else the pipeline's 900 s
    # default (decompose_parallel._WORKER_EXPOSURE_SECONDS).
    _mz_exp, _mz_exp_src = 900.0, "assumed_900s"
    for _c in ("exposure_seconds", "exptime"):
        if _c in _mz_meta_names:
            _v = float(_mz_meta[_c])
            if np.isfinite(_v) and _v > 0.0:
                _mz_exp, _mz_exp_src = _v, "metadata"
                break

    _mz_date_obs = (_mz_meta["date_obs"].decode().strip()
                    if isinstance(_mz_meta["date_obs"], bytes)
                    else str(_mz_meta["date_obs"]).strip())
    _mz_model = MoonZodiPhysicalModel()
    _mz_roles = {
        "near": ("sky_near", "sky_near_ra", "sky_near_dec"),
        "far":  ("sky_far",  "sky_far_ra",  "sky_far_dec"),
        "sci":  ("sci",      "sci_ra",      "sci_dec"),
    }
    for _arm, (_role, _rac, _decc) in _mz_roles.items():
        _obs = MoonZodiObservation(
            expnum=int(_mz_meta["expnum"]),
            date_obs=_mz_date_obs,
            role=_role,
            target_ra_deg=float(_mz_meta[_rac]),
            target_dec_deg=float(_mz_meta[_decc]),
            exposure_seconds=_mz_exp,
            exposure_seconds_source=_mz_exp_src,
        )
        try:
            # astropy warns about IERS coverage for these epochs; the model
            # deliberately pins the packaged table (compute_midpoint_geometry
            # sets iers.conf.auto_download=False), so the warning is expected.
            with _warnings.catch_warnings():
                _warnings.simplefilter("ignore")
                _pr = _mz_model.predict(wave_row, _mz_lsf[_arm], _obs,
                                        physical_to_fit_flux_scale=FACTOR)
            # predict() returns fit-flux units (scaled by FACTOR); divide back to
            # physical so these match every other trace in this cell, which is
            # stored physical and multiplied by FACTOR at plot time.
            _mz_overlay[_arm] = {
                "moon": np.asarray(_pr.moon, dtype=np.float64) / FACTOR,
                "zodi": np.asarray(_pr.zodi, dtype=np.float64) / FACTOR,
                "state": _pr.state,
            }
        except MoonZodiInvalidObservationError as _exc:
            print(f"  moon/zodi model: {_arm} arm not modellable ({_exc.reason})")
except _MzOverlayDisabled:
    pass
except Exception as _exc:  # missing data bundle, ephemeris, META columns, ...
    print(f"  moon/zodi model overlay unavailable: "
          f"{type(_exc).__name__}: {_exc}")

if _mz_overlay:
    # Compare the physical model against the amplitudes the QP actually fitted
    # (and, for sci, against what the ML predicted).  Ratios are band-integrated
    # so they are insensitive to per-pixel noise.
    _mz_decomp = {
        "near": comps_near_from_near,
        "far": comps_far_from_far,
        "sci": comps_sci_true,
    }
    _mz_rows = []
    for _arm, _ov in _mz_overlay.items():
        _geo = _ov["state"].geometry
        _row = {
            "arm": _arm,
            "moon_alt": _geo.moon_altitude_deg,
            "moon_sep": _geo.moon_separation_deg,
            "phase": _geo.signed_phase_deg,
            "zodi_b500": _geo.zodi_b500,
        }
        _mod_sum = _fit_sum = 0.0
        for _fam in ("moon", "zodi"):
            _mod = float(np.nansum(_ov[_fam]))
            _fit = float(np.nansum(np.asarray(
                _mz_decomp[_arm].get(_fam, 0.0), dtype=np.float64) / FACTOR))
            _row[f"{_fam}_fit/model"] = (_fit / _mod if abs(_mod) > 0 else np.nan)
            _mod_sum += _mod
            _fit_sum += _fit
        # moon and zodi are both reddened solar continua, so the QP can trade
        # amplitude between them almost freely (§1.2.2).  The combined ratio is
        # the identifiable quantity: if it sits near 1 while the two individual
        # ratios are far off, the disagreement is a *split* problem, not an
        # amplitude problem -- and only the split is degenerate.
        _row["(moon+zodi)_fit/model"] = (_fit_sum / _mod_sum
                                         if abs(_mod_sum) > 0 else np.nan)
        _mz_rows.append(_row)
    # Same ratio for the ML prediction at the sci pointing.
    if "sci" in _mz_overlay:
        _row = {"arm": "sci (ML pred)", "moon_alt": np.nan, "moon_sep": np.nan,
                "phase": np.nan, "zodi_b500": np.nan}
        _mod_sum = _pred_sum = 0.0
        for _fam in ("moon", "zodi"):
            _mod = float(np.nansum(_mz_overlay["sci"][_fam]))
            _pred = float(np.nansum(np.asarray(
                comps_sci.get(_fam, 0.0), dtype=np.float64) / FACTOR))
            _row[f"{_fam}_fit/model"] = (_pred / _mod if abs(_mod) > 0 else np.nan)
            _mod_sum += _mod
            _pred_sum += _pred
        _row["(moon+zodi)_fit/model"] = (_pred_sum / _mod_sum
                                         if abs(_mod_sum) > 0 else np.nan)
        _mz_rows.append(_row)
    _mz_state0 = next(iter(_mz_overlay.values()))["state"]
    print(f"  moon/zodi physical model {_mz_state0.model_id} "
          f"({_mz_state0.formula_version}), exposure {_mz_exp:.0f}s "
          f"[{_mz_exp_src}], flags={_mz_state0.flags}")
    print(f"    scientific_status = {_mz_state0.scientific_status!r}")
    print(f"    correction_scope  = {_mz_state0.correction_scope!r}")
    print(pd.DataFrame(_mz_rows).to_string(
        index=False, float_format=lambda v: f'{v:.3g}', na_rep='-'))
    # Read the ratios with the model's own scope in mind.  The fitted
    # correction is applied to the moon+zodi SUM (CORRECTION_SCOPE =
    # 'moon_plus_zodi'), so the combined column is the only one the model is
    # calibrated to reproduce.  The individual moon and zodi columns compare
    # against vectors the fit never constrained separately, so a large split
    # discrepancy there is NOT evidence that the QP mis-assigned the families
    # -- the two are degenerate in the model exactly as they are in the QP.
    print("    (moon+zodi)_fit/model is the calibrated comparison: ~1 means the "
          "total continuum")
    print("    amplitude agrees with the frozen physical prediction.  The "
          "per-family columns are")
    print("    indicative only -- the model's correction scope is the sum, so it "
          "does not claim to")
    print("    split moon from zodi any better than the decomposition does.")
    print("    NB scientific_status marks this model diagnostic-only; use it as "
          "a sanity reference,")
    print("    not as truth.")
    print()


# 6) Four-panel diagnostic plot:
#    row1: near observed vs reconstruction from near coefficients
#    row2: far observed vs reconstruction from far coefficients
#    row3: science true vs science prediction
#    row4: science relative residual
fig = make_subplots(
    rows=4,
    cols=1,
    shared_xaxes=True,
    vertical_spacing=0.04,
    subplot_titles=(
        "Near: observed vs reconstructed from near coefficients",
        "Far: observed vs reconstructed from far coefficients",
        "Science: observed / recon(sci coef) / recon(pred)",
        "Science residual: (pred - true) / true",
    ),
    row_heights=[0.24, 0.24, 0.34, 0.18],
)

# Physical Moon/Zodi model overlays, one pair per flux panel.  Dotted so they
# read as an external reference rather than as data or reconstruction.
def _add_mz_traces(_arm, _row):
    _ov = _mz_overlay.get(_arm)
    if _ov is None:
        return
    _mz_curves = (
        ("moon", _ov["moon"], "#9467bd", "dot", 1.2),
        ("zodi", _ov["zodi"], "#17becf", "dot", 1.2),
        # The model's correction scope is moon_plus_zodi, so the sum is the only
        # calibrated curve here -- drawn heavier than its two parts.
        ("moon+zodi", _ov["moon"] + _ov["zodi"], "#8c564b", "dashdot", 1.7),
    )
    for _fam, _y, _color, _dash, _w in _mz_curves:
        fig.add_trace(
            go.Scattergl(
                x=wave_row,
                y=_y * FACTOR,
                mode="lines",
                name=f"{_fam} physical model",
                legendgroup=f"mz_{_fam}",
                showlegend=(_row == 1),
                line=dict(color=_color, width=_w, dash=_dash),
            ),
            row=_row,
            col=1,
        )


fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_near_row * FACTOR,
        mode="lines",
        name="near true",
        line=dict(color="#7f7f7f", width=1.0),
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_near_recon_row * FACTOR,
        mode="lines",
        name="near recon(from near coef)",
        line=dict(color="#e41a1c", width=1.4),
    ),
    row=1,
    col=1,
)
_add_mz_traces("near", 1)

fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_far_row * FACTOR,
        mode="lines",
        name="far true",
        line=dict(color="#7f7f7f", width=1.0),
    ),
    row=2,
    col=1,
)
fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_far_recon_row * FACTOR,
        mode="lines",
        name="far recon(from far coef)",
        line=dict(color="#ff7f00", width=1.4),
    ),
    row=2,
    col=1,
)
_add_mz_traces("far", 2)

fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_sci_true_row * FACTOR,
        mode="lines",
        name="science observed",
        line=dict(color="#7f7f7f", width=1.0),
    ),
    row=3,
    col=1,
)
fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_sci_true_recon_row * FACTOR,
        mode="lines",
        name="science recon(from sci coef)",
        line=dict(color="#2ca02c", width=1.2, dash="dash"),
    ),
    row=3,
    col=1,
)
fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=flux_sci_pred_row * FACTOR,
        mode="lines",
        name="science recon(pred)",
        line=dict(color="#1f78b4", width=1.4),
    ),
    row=3,
    col=1,
)
_add_mz_traces("sci", 3)

fig.add_trace(
    go.Scattergl(
        x=wave_row,
        y=rel_resid_row,
        mode="lines",
        name="science residual",
        line=dict(color="#d62728", width=1.0),
    ),
    row=4,
    col=1,
)
fig.add_hline(y=0, line=dict(color="black", width=0.8, dash="dash"), row=4, col=1)

fig.update_yaxes(type="log", title_text="Near flux", row=1, col=1)
fig.update_yaxes(type="log", title_text="Far flux", row=2, col=1)
fig.update_yaxes(type="log", title_text="Science flux", row=3, col=1)
fig.update_yaxes(type="linear", title_text="(pred-true)/true", row=4, col=1)
fig.update_xaxes(title_text="Wavelength [A]", row=4, col=1)

fig.update_layout(
    template="plotly_white",
    title=dict(
        text=(
            f"{_flag_prefix} Every10 row {idx_row}<br>"
            f"<sub>pRMSE near / far / sci = {rmse_near_recon:.3g} / "
            f"{rmse_far_recon:.3g} / {rmse_row:.3g}  ·  pWRMSE = "
            f"{wrmse_near_recon:.3g} / {wrmse_far_recon:.3g} / "
            f"{wrmse_row_pix:.3g}  ·  sci display pRMSE = "
            f"{rmse_row_display:.3g}</sub><br>"
            f"<sub>{_pct_subline}</sub>"
        ),
        font=dict(size=13),
        x=0.02, xanchor='left',
        y=0.995, yanchor='top',
    ),
    height=1220,
    margin=dict(t=110, r=20, l=70, b=90),
    legend=dict(
        orientation="h",
        yanchor="top", y=-0.05,
        xanchor="left", x=0.0,
        font=dict(size=10),
    ),
)
fig.show()

# 7) Moon spline coefficient diagnostic figure (global-prior result).
moon_axis = np.arange(moon_idx.size)
fig_moon = go.Figure()
fig_moon.add_trace(
    go.Scatter(
        x=moon_axis,
        y=moon_near,
        mode="lines+markers",
        name="near",
        line=dict(color="#7f7f7f"),
    )
)
fig_moon.add_trace(
    go.Scatter(
        x=moon_axis,
        y=moon_far,
        mode="lines+markers",
        name="far",
        line=dict(color="#bdbdbd"),
    )
)
fig_moon.add_trace(
    go.Scatter(
        x=moon_axis,
        y=moon_true,
        mode="lines+markers",
        name="sci true",
        line=dict(color="#1f78b4"),
    )
)
fig_moon.add_trace(
    go.Scatter(
        x=moon_axis,
        y=moon_pred,
        mode="lines+markers",
        name="pred default",
        line=dict(color="#e41a1c"),
    )
)
fig_moon.update_layout(
    template="plotly_white",
    title=dict(
        text=(f"{_flag_prefix} Moon spline coefficients — row {idx_row}<br>"
              f"<sub>{_pct_subline}</sub>"),
        font=dict(size=12),
        x=0.02, xanchor='left',
    ),
    xaxis_title="Moon_bs coefficient index",
    yaxis_title="coefficient value",
    height=440,
    margin=dict(t=80),
    legend=dict(font=dict(size=10)),
)
fig_moon.show()

# 7b) Zodi spline coefficient diagnostic figure (split_zodi=True corpus).
if zodi_idx.size:
    zodi_axis = np.arange(zodi_idx.size)
    fig_zodi = go.Figure()
    for _y, _name, _color in (
        (zodi_near, 'near',        '#7f7f7f'),
        (zodi_far,  'far',         '#bdbdbd'),
        (zodi_true, 'sci true',    '#1f78b4'),
        (zodi_pred, 'pred default','#e41a1c'),
    ):
        fig_zodi.add_trace(go.Scatter(
            x=zodi_axis, y=_y, mode='lines+markers',
            name=_name, line=dict(color=_color),
        ))
    fig_zodi.update_layout(
        template='plotly_white',
        title=dict(
            text=(f"{_flag_prefix} Zodi spline coefficients — row {idx_row}<br>"
                  f"<sub>{_pct_subline}</sub>"),
            font=dict(size=12),
            x=0.02, xanchor='left',
        ),
        xaxis_title='coefficient index',
        yaxis_title='coefficient value',
        height=380,
        margin=dict(t=80),
        legend=dict(font=dict(size=10)),
    )
    fig_zodi.show()

# 8) Per-component reconstructions for all four arms, mirroring the
#    four-trace layout of the Moon spline diagnostic above but as
#    spectra over wavelength rather than coefficients over index.
#    comps_sci_true was reconstructed in section 4 so row 3 can use it.
_comps_by_arm = {
    "near": comps_near_from_near,
    "far": comps_far_from_far,
    "sci true": comps_sci_true,
    "pred default": comps_sci,
}
_arm_colors = {
    "near": "#7f7f7f",
    "far": "#bdbdbd",
    "sci true": "#1f78b4",
    "pred default": "#e41a1c",
}

# comps["*"] is in the same display scale as comps["total"], so no *FACTOR here.
def _non_moon_continuum(comps):
    return np.asarray(comps["diffuse"], dtype=np.float64)

def _line_component(comps):
    return (np.asarray(comps["oh"], dtype=np.float64)
            + np.asarray(comps["atom"], dtype=np.float64)
            + np.asarray(comps["orc"], dtype=np.float64)
            + np.asarray(comps["o2"], dtype=np.float64))

def _moon_spectrum(comps):
    return np.asarray(comps["moon"], dtype=np.float64)

def _zodi_spectrum(comps):
    if "zodi" not in comps:
        return np.zeros_like(comps["moon"])
    return np.asarray(comps["zodi"], dtype=np.float64)

# Sanity: total is defined by reconstruct_component_spectra as
#   oh + moon + diffuse + atom + orc + o2
# so lines + moon_spectrum + non_moon_continuum must equal it.
_total_pred = np.asarray(comps_sci["total"], dtype=np.float64)
_sum_pred = (_line_component(comps_sci)
             + _moon_spectrum(comps_sci)
             + _non_moon_continuum(comps_sci))
_max_diff = float(np.nanmax(np.abs(_total_pred - _sum_pred)))
_max_rel = float(np.nanmax(np.abs(_total_pred - _sum_pred)
                          / np.clip(np.abs(_total_pred), 1e-30, None)))
print(f"Component-sum check (pred): max abs diff = {_max_diff:.3g}, "
      f"max rel diff = {_max_rel:.3g}")

fig_continuum = go.Figure()
for _arm, _comps in _comps_by_arm.items():
    fig_continuum.add_trace(
        go.Scattergl(
            x=wave_row,
            y=_non_moon_continuum(_comps),
            mode="lines",
            name=_arm,
            line=dict(color=_arm_colors[_arm], width=1.2),
        )
    )
fig_continuum.update_layout(
    template="plotly_white",
    title=dict(
        text=(f"{_flag_prefix} Reconstructed non-moon continuum "
              f"(diffuse = HO2 + FeO + O2ac) — row {idx_row}"),
        font=dict(size=12),
        x=0.02, xanchor='left',
    ),
    xaxis_title="Wavelength [A]",
    yaxis_title=f"Flux (display units x{FACTOR:.3g})",
    height=440,
    margin=dict(t=60),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0, font=dict(size=10)),
)
fig_continuum.show()

# Physical-model overlay for the single-family spline-spectrum panels below.
# One trace per physical pointing: "sci true" and "pred default" share the sci
# pointing, so the model contributes three curves, not four.
def _add_mz_family_traces(_target_fig, _fam):
    for _arm, _color in (("near", "#7f7f7f"), ("far", "#bdbdbd"), ("sci", "#1f78b4")):
        _ov = _mz_overlay.get(_arm)
        if _ov is None:
            continue
        _target_fig.add_trace(
            go.Scattergl(
                x=wave_row,
                # comps[...] in these panels is plotted in fit units without a
                # *FACTOR, while _mz_overlay is stored physical -- hence *FACTOR.
                y=_ov[_fam] * FACTOR,
                mode="lines",
                name=f"{_arm} physical model",
                line=dict(color=_color, width=1.6, dash="dot"),
            )
        )


fig_moon_spectrum = go.Figure()
for _arm, _comps in _comps_by_arm.items():
    fig_moon_spectrum.add_trace(
        go.Scattergl(
            x=wave_row,
            y=_moon_spectrum(_comps),
            mode="lines",
            name=_arm,
            line=dict(color=_arm_colors[_arm], width=1.2),
        )
    )
_add_mz_family_traces(fig_moon_spectrum, "moon")
fig_moon_spectrum.update_layout(
    template="plotly_white",
    title=dict(
        text=(f"{_flag_prefix} Reconstructed moon spline spectrum "
              f"(comps['moon']) — row {idx_row}"),
        font=dict(size=12),
        x=0.02, xanchor='left',
    ),
    xaxis_title="Wavelength [A]",
    yaxis_title=f"Flux (display units x{FACTOR:.3g})",
    height=440,
    margin=dict(t=60),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0, font=dict(size=10)),
)
fig_moon_spectrum.show()

# Reconstructed zodi spline spectrum (split_zodi Zodi_bs family).
if 'zodi' in comps_sci:
    fig_zodi_spectrum = go.Figure()
    for _arm, _comps in _comps_by_arm.items():
        fig_zodi_spectrum.add_trace(
            go.Scattergl(
                x=wave_row,
                y=_zodi_spectrum(_comps),
                mode='lines',
                name=_arm,
                line=dict(width=1.4),
            )
        )
    _add_mz_family_traces(fig_zodi_spectrum, 'zodi')
    fig_zodi_spectrum.update_layout(
        template='plotly_white',
        title=dict(
            text=(f"{_flag_prefix} Reconstructed zodi spline spectrum "
                  f"(comps['zodi']) — row {idx_row}"),
            font=dict(size=12),
            x=0.02, xanchor='left',
        ),
        xaxis_title='Wavelength (A)',
        yaxis_title='flux (fit units)',
        height=400,
        margin=dict(t=60),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0.0, font=dict(size=10)),
    )
    fig_zodi_spectrum.show()

fig_lines = go.Figure()
for _arm, _comps in _comps_by_arm.items():
    fig_lines.add_trace(
        go.Scattergl(
            x=wave_row,
            y=_line_component(_comps),
            mode="lines",
            name=_arm,
            line=dict(color=_arm_colors[_arm], width=1.2),
        )
    )
fig_lines.update_layout(
    template="plotly_white",
    title=dict(
        text=(f"{_flag_prefix} Reconstructed line emission "
              f"(OH + atom + ORC + O2) — row {idx_row}"),
        font=dict(size=12),
        x=0.02, xanchor='left',
    ),
    xaxis_title="Wavelength [A]",
    yaxis_title=f"Flux (display units x{FACTOR:.3g})",
    height=480,
    margin=dict(t=60),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0, font=dict(size=10)),
)
fig_lines.show()

# 9) Per-component (pred - sci_true) residual spectrum. Row 3 shows only
#    the total residual; this decomposes it into moon / diffuse / lines
#    so a small broadband deficit in a component that spans the whole
#    wavelength range is visible even when the component panels above
#    make it look "close" on a linear-y axis. By construction the three
#    traces must sum to the blue-minus-green curve of row 3, and that
#    sum is drawn as a black dashed reference.
_delta_moon = _moon_spectrum(comps_sci) - _moon_spectrum(comps_sci_true)
_delta_zodi = _zodi_spectrum(comps_sci) - _zodi_spectrum(comps_sci_true)
_delta_diffuse = _non_moon_continuum(comps_sci) - _non_moon_continuum(comps_sci_true)
_delta_lines = _line_component(comps_sci) - _line_component(comps_sci_true)
_delta_total = _delta_moon + _delta_zodi + _delta_diffuse + _delta_lines

fig_deltas = go.Figure()
for _label, _y, _color in (
    ("moon (pred - sci recon)", _delta_moon, "#e41a1c"),
    ("diffuse (pred - sci recon)", _delta_diffuse, "#377eb8"),
    ("lines (pred - sci recon)", _delta_lines, "#4daf4a"),
    ("total (pred - sci recon)", _delta_total, "#000000"),
):
    fig_deltas.add_trace(
        go.Scattergl(
            x=wave_row,
            y=_y,
            mode="lines",
            name=_label,
            line=dict(color=_color, width=1.2,
                      dash="dash" if _label.startswith("total") else "solid"),
        )
    )
fig_deltas.add_hline(y=0, line=dict(color="rgba(0,0,0,0.4)", width=0.8, dash="dot"))
fig_deltas.update_layout(
    template="plotly_white",
    title=dict(
        text=(f"{_flag_prefix} Per-component prediction minus sci-arm "
              f"reconstruction (linear) — row {idx_row}"),
        font=dict(size=12),
        x=0.02, xanchor='left',
    ),
    xaxis_title="Wavelength [A]",
    yaxis_title=f"Delta flux (display units x{FACTOR:.3g})",
    height=480,
    margin=dict(t=60),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0, font=dict(size=10)),
)
fig_deltas.show()

# 10) Numeric integrated deltas per component in three wavelength bands, so
#     the sign and magnitude of each contribution to the red deficit is
#     visible even where the plot traces are noisy line-by-line.
_bands = [
    ("blue  (< 5500 A)", wave_row < 5500.0),
    ("green (5500-7500)", (wave_row >= 5500.0) & (wave_row < 7500.0)),
    ("red   (>= 7500 A)", wave_row >= 7500.0),
]
print()
print("Integrated (pred - sci recon) per component and wavelength band:")
print(f"  {'band':<18s} {'moon':>12s} {'zodi':>12s} {'diffuse':>12s} {'lines':>12s} {'total':>12s}")
for _bname, _mask in _bands:
    if not _mask.any():
        continue
    _sm = float(np.nansum(_delta_moon[_mask]))
    _sd = float(np.nansum(_delta_diffuse[_mask]))
    _sl = float(np.nansum(_delta_lines[_mask]))
    _st = float(np.nansum(_delta_total[_mask]))
    print(f"  {_bname:<18s} {_sm:>+12.4g} {_sd:>+12.4g} {_sl:>+12.4g} {_st:>+12.4g}")
