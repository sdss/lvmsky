# Sky-arm disagreement floor (irreducible-noise detector).
#
# For each test row r and each group g, compute the intrinsic-emissivity
# disagreement between the two sky arms and the ML absolute error:
#     Delta_g(r) = sqrt( mean_{j in g} (c_near_g[j]/G_near_g[j]
#                                       - c_far_g[j]/G_far_g[j])**2 ) / 2
#     err_g(r)   = sqrt( mean_{j in g} (c_pred_g[j] - c_true_g[j])**2 )
# Delta measures the gravity-wave / systematic noise seen simultaneously by
# the two sky arms; err is what the ML head produced.  If err ~ Delta, the
# model has reached the physically achievable floor -- pushing further is
# limited by intrinsic sky variability, not by the network.
#
# Per group we report:
#   median(err/Delta), P50(err), P50(Delta), and the fraction of rows where
#   err <= 1.5 * Delta (on-floor).
# Higher err/Delta = more headroom for ML improvement.

required = ['filtered_triplet', 'compress_geom_kwargs', '_group_indices_compress',
            'coef_near_all', 'coef_far_all', 'coef_sci_all',
            'ctx_near_all', 'ctx_far_all', 'ctx_sci_all',
            'test_idx', 'coef_pred_det']
_missing = [k for k in required if k not in globals()]
if _missing:
    raise RuntimeError('Run prerequisite cells first. Missing: ' + ', '.join(_missing))

_te = np.asarray(test_idx, dtype=int)
_ctx_near_te = np.asarray(ctx_near_all[_te], dtype=np.float64)
_ctx_far_te  = np.asarray(ctx_far_all[_te],  dtype=np.float64)
_ctx_sci_te  = np.asarray(ctx_sci_all[_te],  dtype=np.float64)
_coef_near_te = np.asarray(coef_near_all[_te], dtype=np.float64)
_coef_far_te  = np.asarray(coef_far_all[_te],  dtype=np.float64)
_coef_sci_te  = np.asarray(coef_sci_all[_te],  dtype=np.float64)
_coef_pred_te = np.asarray(coef_pred_det, dtype=np.float64)

_G_near = airglow_geometry_scale(_ctx_near_te, **compress_geom_kwargs)
_G_far  = airglow_geometry_scale(_ctx_far_te,  **compress_geom_kwargs)
_em_near = _coef_near_te / _G_near
_em_far  = _coef_far_te  / _G_far

# Row-subset masks on ctx_sci.
_ctx_names_ll = list(filtered_triplet['ctx_names'])
_ma_col = _ctx_names_ll.index('moon_alt') if 'moon_alt' in _ctx_names_ll else None
_eb_col = _ctx_names_ll.index('ecl_beta_deg') if 'ecl_beta_deg' in _ctx_names_ll else None

_masks = {'all': np.ones(_te.size, dtype=bool)}
if _ma_col is not None:
    _moon_alt_te = _ctx_sci_te[:, _ma_col]
    _masks['moon_up']   = _moon_alt_te > 0.0
    _masks['moon_down'] = _moon_alt_te <= 0.0
if _eb_col is not None:
    _abs_eb_te = np.abs(_ctx_sci_te[:, _eb_col])
    _masks['close_zodi'] = _abs_eb_te < 20.0

_group_names = list(_group_indices_compress.keys())

# --- Per-group scalar Delta and err ---
_delta_per_group = {}   # (n_te,) per group
_err_per_group   = {}
for gname, idx in _group_indices_compress.items():
    idx = np.asarray(idx, dtype=int)
    _diff = _em_near[:, idx] - _em_far[:, idx]
    _delta_per_group[gname] = 0.5 * np.sqrt(np.mean(_diff * _diff, axis=1))
    _resid = _coef_pred_te[:, idx] - _coef_sci_te[:, idx]
    _err_per_group[gname]   = np.sqrt(np.mean(_resid * _resid, axis=1))

print('=' * 90)
print(f'Sky-arm disagreement floor vs ML error ({_te.size} rows, {len(_group_names)} groups)')
print('=' * 90)

for regime, mask in _masks.items():
    _n = int(mask.sum())
    if _n < 10:
        continue
    print()
    print(f'--- regime: {regime}  (n_rows = {_n}) ---')
    _rows_out = []
    for gname in _group_names:
        _d = _delta_per_group[gname][mask]
        _e = _err_per_group[gname][mask]
        _finite = np.isfinite(_d) & np.isfinite(_e) & (_d > 1e-30)
        if _finite.sum() < 5:
            continue
        _dd = _d[_finite]
        _ee = _e[_finite]
        _ratio = _ee / _dd
        _rows_out.append({
            'group':          gname,
            'p50_delta':      float(np.median(_dd)),
            'p50_err':        float(np.median(_ee)),
            'p50_err/delta':  float(np.median(_ratio)),
            'p16_err/delta':  float(np.percentile(_ratio, 16)),
            'p84_err/delta':  float(np.percentile(_ratio, 84)),
            'on-floor_frac':  float(np.mean(_ratio <= 1.5)),
        })
    if _rows_out:
        print(pd.DataFrame(_rows_out).to_string(index=False,
              float_format=lambda v: f'{v:.4g}'))

print()
print('Interpretation:')
print('  p50_err/delta ~ 1  -> ML is at the sky-arm noise floor for this group/regime.')
print('  p50_err/delta > 2  -> substantial ML headroom above the floor.')
print('  p50_err/delta < 1  -> ML uses ctx information the arms do not carry (a '
      'good sign: the head is genuinely predictive).')
print('  on-floor_frac      -> fraction of rows where ML error <= 1.5 * Delta.')

# --- Scatter plot: ML err vs Delta, one panel per group, on the "all" regime ---
_n_g = len(_group_names)
_ncols = min(3, _n_g)
_nrows = int(np.ceil(_n_g / _ncols))
_fig = make_subplots(rows=_nrows, cols=_ncols,
                    subplot_titles=[f'{g}: err vs sky-arm Delta' for g in _group_names])
for _i, gname in enumerate(_group_names):
    _row = _i // _ncols + 1
    _col = _i %  _ncols + 1
    _d = _delta_per_group[gname]
    _e = _err_per_group[gname]
    _finite = np.isfinite(_d) & np.isfinite(_e) & (_d > 1e-30)
    _dd = _d[_finite]
    _ee = _e[_finite]
    _fig.add_trace(go.Scatter(
        x=_dd, y=_ee, mode='markers',
        marker=dict(size=4, opacity=0.5, color='steelblue'),
        showlegend=False, hoverinfo='skip'), row=_row, col=_col)
    # Percentile-based axis limits so a handful of near-zero Delta values
    # (solver-noise-level sky-arm agreement) do not stretch the log axes.
    _lo = max(float(np.percentile(_dd, 1)),
              float(np.percentile(_ee, 1)), 1e-10)
    _hi = max(float(np.percentile(_dd, 99)),
              float(np.percentile(_ee, 99)), _lo * 10)
    _grid = np.array([_lo, _hi])
    _fig.add_trace(go.Scatter(x=_grid, y=_grid, mode='lines',
                              line=dict(color='black', dash='dot', width=1),
                              showlegend=False, hoverinfo='skip'),
                   row=_row, col=_col)
    _log_lo = float(np.log10(_lo))
    _log_hi = float(np.log10(_hi))
    _fig.update_xaxes(title_text='Delta (sky-arm intrinsic disagreement)',
                      type='log', range=[_log_lo, _log_hi], row=_row, col=_col)
    _fig.update_yaxes(title_text='|ML pred - true|',
                      type='log', range=[_log_lo, _log_hi], row=_row, col=_col)
_fig.update_layout(height=280 * _nrows, width=340 * _ncols,
                   title='ML per-row RMS error vs sky-arm disagreement floor '
                         '(diagonal = at the floor)')
_fig.show()
