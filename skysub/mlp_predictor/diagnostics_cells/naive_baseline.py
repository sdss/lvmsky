# Naive baselines vs the ML model on the same night-held-out test split (§11.1).
#
#   B0  copy_near        hat{c} = coef_near                                (no physics)
#   B1  near_geo         hat{c} = em_near * G_sci                          (near arm re-projected onto sci geometry)
#   B2  mean_geo         hat{c} = 0.5 * (em_near + em_far) * G_sci         (geometry-corrected symmetric average)
#   ML_default           coef_pred_det (single-seed model)
#   ML_ensemble          arithmetic-mean of the 4-seed ensemble (if multi-seed cell has run)
#
# em_arm = coef_arm / G_arm; G_arm = airglow_geometry_scale(ctx_arm).  G_arm is exactly 1 for
# the moon and other non-airglow groups, so B1's moon prediction reduces to `coef_near` and
# B2's moon prediction is the plain arithmetic mean of the two sky arms.

required = ['filtered_triplet', 'compress_geom_kwargs', '_group_indices_compress',
            'mlp_artifacts', 'cmp_df', '_metric_row',
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
    'B0_copy_near': np.clip(_coef_near_te,                       0.0, None).astype(np.float32),
    'B1_near_geo':  np.clip(_em_near * _G_sci,                   0.0, None).astype(np.float32),
    'B2_mean_geo':  np.clip(0.5 * (_em_near + _em_far) * _G_sci, 0.0, None).astype(np.float32),
    # ML_default is the 4-seed ensemble prediction (see trainer cell + §12, 2026-08-11).
    'ML_default':   np.asarray(coef_pred_det, dtype=np.float32),
}

_y_te = _coef_sci_te.astype(np.float32)

# 2026-08-16: sigma-aware _metric_row -> mean/median/total WRMSE columns.
_sig_te = (coef_err_sci_all[test_idx].astype(np.float32)
           if 'coef_err_sci_all' in globals() else None)
_floor_te = (dict(DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP)
              if 'DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP' in globals() else None)
_rows = [{'variant': name,
          **_metric_row(_y_te, pred, name,
                        sigma=_sig_te,
                        group_indices=_group_indices_compress,
                        floor_by_group=_floor_te)}
         for name, pred in _all_preds.items()]
_summary_df = pd.DataFrame(_rows)

print('=' * 78)
print(f'Naive baselines vs ML on test split ({_te_bl.size} rows, {_y_te.shape[1]} coefficients)')
print('=' * 78)
print(_summary_df.to_string(index=False, float_format=lambda v: f'{v:.6g}'))

# Per-group mean coefficient bias.
_bias_rows = []
for gname, idx in _group_indices_compress.items():
    idx = np.asarray(idx, dtype=int)
    _mean_true = float(np.mean(_y_te[:, idx]))
    row = {'group': gname, 'n': int(idx.size), 'mean_true': _mean_true}
    for name, pred in _all_preds.items():
        row[f'bias_{name}_%'] = (
            100.0 * (float(np.mean(pred[:, idx])) - _mean_true)
            / max(abs(_mean_true), 1e-30))
    _bias_rows.append(row)
print()
print('Per-group mean coefficient bias on test rows (positive = over-predict, %):')
print(pd.DataFrame(_bias_rows).to_string(index=False, float_format=lambda v: f'{v:.3g}'))

# Per-group mean per-coefficient RMSE.
_rmse_group_rows = []
for gname, idx in _group_indices_compress.items():
    idx = np.asarray(idx, dtype=int)
    row = {'group': gname, 'n': int(idx.size)}
    for name, pred in _all_preds.items():
        _rmse_j = np.sqrt(np.mean(
            (pred[:, idx].astype(np.float64) - _y_te[:, idx].astype(np.float64)) ** 2,
            axis=0))
        row[f'eRMSE_{name}'] = float(np.mean(_rmse_j))
    _rmse_group_rows.append(row)
print()
print('Per-group mean per-coefficient eRMSE (lower is better):')
print(pd.DataFrame(_rmse_group_rows).to_string(index=False, float_format=lambda v: f'{v:.4g}'))

# Verdict against the best naive baseline.
# 2026-08-14 update: use median_eRMSE as the primary merit metric because
# mean_eRMSE is dominated by a small coefficient tail.
_bl_rows = _summary_df[~_summary_df['variant'].str.startswith('ML_')].copy()
_ml_row  = _summary_df[_summary_df['variant'] == 'ML_default'].iloc[0]
_best_bl = _bl_rows.sort_values(['median_eRMSE', 'mean_eRMSE']).iloc[0]

_ml_med_rmse = float(_ml_row['median_eRMSE'])
_best_bl_med_rmse = float(_best_bl['median_eRMSE'])
_gain_med_pct = 100.0 * (_best_bl_med_rmse - _ml_med_rmse) / max(_best_bl_med_rmse, 1e-30)

_ml_mean_rmse = float(_ml_row['mean_eRMSE'])
_best_bl_mean_rmse = float(_best_bl['mean_eRMSE'])
_gain_mean_pct = 100.0 * (_best_bl_mean_rmse - _ml_mean_rmse) / max(_best_bl_mean_rmse, 1e-30)

print()
print(f"Best naive baseline (by median_eRMSE): {_best_bl['variant']!r}")
print(f"  median_eRMSE={_best_bl_med_rmse:.4g}, mean_eRMSE={_best_bl_mean_rmse:.4g}, "
      f"median_corr={float(_best_bl['median_corr']):.4f}")
print("ML default:")
print(f"  median_eRMSE={_ml_med_rmse:.4g}, mean_eRMSE={_ml_mean_rmse:.4g}, "
      f"median_corr={float(_ml_row['median_corr']):.4f}")
print(f'ML improvement over best baseline: median_eRMSE {_gain_med_pct:+.1f}%  '
      f'(mean_eRMSE context: {_gain_mean_pct:+.1f}%)')

if _ml_med_rmse < _best_bl_med_rmse:
    print('Verdict (median_eRMSE-primary): the network earns its complexity on this test set.')
elif abs(_gain_med_pct) < 2.0:
    print('Verdict (median_eRMSE-primary): ML and best baseline are within noise; '
          'check per-group eRMSE and outlier-tail diagnostics for practical differences.')
else:
    print('Verdict (median_eRMSE-primary): the best naive baseline beats the ML model. '
          'Either the model overfits, the loss overfocuses a difficult tail, or '
          'the compressor discards directions the baseline preserves.')
