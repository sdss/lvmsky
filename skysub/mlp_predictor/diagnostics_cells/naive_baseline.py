# Naive baselines vs the ML model on the night-held-out test split.
#
# Metric:
#   sRMSE_g(r) = sqrt( mean_{j in g} (pred[r, j] - true[r, j])**2 )
#   sRMSE_g    = mean_r sRMSE_g(r)   (also report median_r for outlier robustness)
#
# sRMSE is a per-spectrum per-group figure of merit: it asks "how well does
# variant X reproduce group g on each row" and then averages over rows.
# Unlike eRMSE (per-column RMSE, averaged over columns), sRMSE gives each
# spectrum equal weight regardless of how many coefficients the group has.
# That kills the pathology where the ~400 OH columns drown out the moon+zodi
# groups the ML is actually predicting.
#
# Baselines:
#   B0  copy_near   hat{c} = coef_near                             (no physics)
#   B1  near_geo    hat{c} = em_near * G_sci                       (near arm re-projected onto sci geometry)
#   B2  mean_geo    hat{c} = 0.5*(em_near + em_far) * G_sci        (geometry-corrected symmetric average)
#   ML_default      coef_pred_det (ensemble mean)
#
# em_arm = coef_arm / G_arm; G_arm = airglow_geometry_scale(ctx_arm). G_arm is
# exactly 1 for moon/zodi, so B1's moon prediction reduces to coef_near_moon
# and B2's moon prediction is the plain arithmetic mean of the two sky arms.
#
# Row-subset regimes (masks on ctx_sci):
#   all           entire test set
#   moon_up       moon_alt > 0
#   moon_down     moon_alt <= 0
#   close_zodi    |ecl_beta_deg| < 20  (near the ecliptic plane -> zodi bright)

required = ['filtered_triplet', 'compress_geom_kwargs', '_group_indices_compress',
            'coef_near_all', 'coef_far_all', 'coef_sci_all',
            'ctx_near_all', 'ctx_far_all', 'ctx_sci_all',
            'test_idx', 'coef_pred_det']
_missing = [k for k in required if k not in globals()]
if _missing:
    raise RuntimeError('Run prerequisite cells first. Missing: ' + ', '.join(_missing))

_te_bl = np.asarray(test_idx, dtype=int)
_ctx_near_te = np.asarray(ctx_near_all[_te_bl], dtype=np.float64)
_ctx_far_te  = np.asarray(ctx_far_all[_te_bl],  dtype=np.float64)
_ctx_sci_te  = np.asarray(ctx_sci_all[_te_bl],  dtype=np.float64)
_coef_near_te = np.asarray(coef_near_all[_te_bl], dtype=np.float64)
_coef_far_te  = np.asarray(coef_far_all[_te_bl],  dtype=np.float64)
_coef_sci_te  = np.asarray(coef_sci_all[_te_bl],  dtype=np.float64)

_G_near = airglow_geometry_scale(_ctx_near_te, **compress_geom_kwargs)
_G_far  = airglow_geometry_scale(_ctx_far_te,  **compress_geom_kwargs)
_G_sci  = airglow_geometry_scale(_ctx_sci_te,  **compress_geom_kwargs)

_em_near = _coef_near_te / _G_near
_em_far  = _coef_far_te  / _G_far

_all_preds = {
    'B0_copy_near': np.clip(_coef_near_te,                       0.0, None).astype(np.float64),
    'B1_near_geo':  np.clip(_em_near * _G_sci,                   0.0, None).astype(np.float64),
    'B2_mean_geo':  np.clip(0.5 * (_em_near + _em_far) * _G_sci, 0.0, None).astype(np.float64),
    'ML_default':   np.asarray(coef_pred_det, dtype=np.float64),
}
_y_te = _coef_sci_te.astype(np.float64)

# --- Row-subset masks (built from sci-arm context) ---
_ctx_names_ll = list(filtered_triplet['ctx_names'])
_ma_col = _ctx_names_ll.index('moon_alt') if 'moon_alt' in _ctx_names_ll else None
_eb_col = _ctx_names_ll.index('ecl_beta_deg') if 'ecl_beta_deg' in _ctx_names_ll else None

_masks = {'all': np.ones(_te_bl.size, dtype=bool)}
if _ma_col is not None:
    _moon_alt_te = _ctx_sci_te[:, _ma_col]
    _masks['moon_up']   = _moon_alt_te > 0.0
    _masks['moon_down'] = _moon_alt_te <= 0.0
else:
    print('  (moon_alt not in ctx_names -> skipping moon_up/moon_down regimes)')
if _eb_col is not None:
    _abs_eb_te = np.abs(_ctx_sci_te[:, _eb_col])
    _masks['close_zodi'] = _abs_eb_te < 20.0
else:
    print('  (ecl_beta_deg not in ctx_names -> skipping close_zodi regime)')


def _srmse_per_group(y_true, y_pred, group_indices, row_mask):
    """Return {group: (mean_sRMSE, median_sRMSE, n_row)} on the masked rows."""
    y_true = np.asarray(y_true, dtype=np.float64)
    y_pred = np.asarray(y_pred, dtype=np.float64)
    _r = np.asarray(row_mask, dtype=bool)
    out = {}
    for gname, idx in group_indices.items():
        idx = np.asarray(idx, dtype=int)
        _resid_g = (y_pred[_r][:, idx] - y_true[_r][:, idx]) ** 2
        _row_rms = np.sqrt(np.mean(_resid_g, axis=1))
        out[gname] = (float(np.mean(_row_rms)),
                      float(np.median(_row_rms)),
                      int(_r.sum()))
    return out


_group_names = list(_group_indices_compress.keys())

print('=' * 90)
print(f'Naive baselines vs ML on test split ({_te_bl.size} rows, '
      f'{_y_te.shape[1]} coefficients, {len(_group_names)} groups)')
print('sRMSE_g(r) = sqrt(mean_{j in g} (pred[r,j] - true[r,j])**2); reported = mean_r.')
print('=' * 90)

for regime, mask in _masks.items():
    _n = int(mask.sum())
    print()
    print(f'--- regime: {regime}  (n_rows = {_n}) ---')
    if _n == 0:
        print('  (empty subset)')
        continue
    _rows_out = []
    for name, pred in _all_preds.items():
        _per_g = _srmse_per_group(_y_te, pred, _group_indices_compress, mask)
        row = {'variant': name}
        for g in _group_names:
            row[g] = _per_g[g][0]
        _rows_out.append(row)
    _df = pd.DataFrame(_rows_out).set_index('variant')
    print(_df.to_string(float_format=lambda v: f'{v:.4g}'))

    print('  ML vs best naive baseline per group '
          '(pct gain positive = ML wins over the best non-ML variant):')
    for g in _group_names:
        _col = _df[g]
        _ml = float(_col.loc['ML_default'])
        _bl_col = _col.drop('ML_default')
        _best_bl_name = _bl_col.idxmin()
        _best_bl_val = float(_bl_col.min())
        _winner = 'ML' if _ml <= _best_bl_val else _best_bl_name
        _gain_pct = 100.0 * (_best_bl_val - _ml) / max(_best_bl_val, 1e-30)
        print(f'    {g:<12s} winner={_winner:<12s}  ML={_ml:.4g}  '
              f'best_naive={_best_bl_name}:{_best_bl_val:.4g}  gain={_gain_pct:+.1f}%')

# Group-equal aggregate on the full test set: mean over groups of the per-variant
# per-group mean sRMSE.  Each group counts once, so moon+zodi are not drowned out
# by the 400-column OH block.
print()
print('=' * 90)
print('Group-equal aggregate on the full test set '
      '(mean over groups of the per-variant per-group mean sRMSE):')
print('=' * 90)
_agg_rows = []
for name, pred in _all_preds.items():
    _per_g = _srmse_per_group(_y_te, pred, _group_indices_compress, _masks['all'])
    _agg_rows.append({
        'variant': name,
        'group_equal_sRMSE': float(np.mean([_per_g[g][0] for g in _group_names])),
    })
_all_df = pd.DataFrame(_agg_rows).set_index('variant')
print(_all_df.to_string(float_format=lambda v: f'{v:.5g}'))

_bl_only = _all_df.drop('ML_default')
_best = _bl_only['group_equal_sRMSE'].idxmin()
_best_v = float(_bl_only['group_equal_sRMSE'].loc[_best])
_ml_v = float(_all_df['group_equal_sRMSE'].loc['ML_default'])
_gain_pct = 100.0 * (_best_v - _ml_v) / max(_best_v, 1e-30)
print()
print(f"Best naive baseline (group-equal aggregate): {_best!r} = {_best_v:.5g}")
print(f"ML_default (group-equal aggregate):          {_ml_v:.5g}")
print(f"ML improvement over best naive: {_gain_pct:+.1f}%")
if _ml_v < _best_v:
    print('Verdict (group-equal sRMSE): the network earns its complexity.')
elif abs(_gain_pct) < 2.0:
    print('Verdict (group-equal sRMSE): ML and best baseline within noise; '
          'inspect per-group / per-regime tables above.')
else:
    print('Verdict (group-equal sRMSE): the best naive baseline beats the ML model. '
          'Check moon/zodi rows in the per-regime tables to see whether the loss '
          'is defeated in the physically relevant regimes.')
