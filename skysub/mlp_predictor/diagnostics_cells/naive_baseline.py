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
    if regime == 'all':
        naive_baseline_per_group = {}
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
        if regime == 'all':
            naive_baseline_per_group[g] = dict(
                ml=_ml, best_naive=_best_bl_val,
                best_naive_name=str(_best_bl_name), gain_pct=_gain_pct)

# --- FLUX-SPACE companion for the mesospheric group -----------------------
# The coefficient sRMSE above weights all 357 OH sticks equally, but the OH
# block is internally DEGENERATE: neighbouring sticks trade amplitude with
# almost no change to the convolved spectrum, so a large coefficient error can
# be invisible in the delivered sky and a small one can dominate it.
#
# This is not hypothetical.  Across the three 2026-09 corpora the mesospheric
# COEFFICIENT gain read +1.2% (no diffuse/OH cap), -3.9% (cap W=0.15) and
# -3.6% (W=0.30) -- twice flagged as "LOSES" -- while exactly the same errors
# projected through the basis gave +5.5%, +4.9% and +5.0%.  The sign flips.
# The flux number is the one that describes the subtracted spectrum, and acting
# on the coefficient number alone cost a decomposition re-run.
#
# NOTE the basis: OH coefficients live on the CONVOLVED STICK matrix, not on
# `matrix_oh` (whose row integrals are 5.01x larger).  Getting that wrong
# silently rescales the residual.
_MES_FLUX_KEYS = ('SkyDecompLSFSurfaceIterative', '_infer_base_dir_for_reconstruction',
                  'N_MOON_KNOTS', 'SPLIT_ZODI', 'N_ZODI_KNOTS',
                  'DECOMP_DATA_ROOT', 'DECOMP_STEM')
naive_baseline_mesospheric_flux = None
if all(k in globals() for k in _MES_FLUX_KEYS):
    try:
        # Explicit: `fits` is in the shared exec namespace for some cells but not
        # reliably for this one, and a NameError here would be reported as the
        # companion being "unavailable" rather than as the trivial fix it is.
        from astropy.io import fits as _fits_bl
        _cn_bl = [str(n) for n in filtered_triplet['coef_names']]
        _oh_bl = np.array([j for j, n in enumerate(_cn_bl)
                           if n.startswith('OH_')], dtype=int)
        if _oh_bl.size:
            with _fits_bl.open(f'{DECOMP_DATA_ROOT}/{DECOMP_STEM}_every10.fits') as _hbl:
                _wbl = np.asarray(_hbl['WAVE'].data, dtype=np.float64)
                _wbl = _wbl if _wbl.ndim == 1 else _wbl[0]
            # Must be THIS corpus's basis: the telluric variant groups OH
            # differently, so projecting its coefficient errors through the
            # split-zodi stick matrix would measure the wrong thing.
            from mlp_predictor.data import (
                make_reconstruction_decomposer as _mk_decomp_bl)
            _mbl = _mk_decomp_bl(
                _wbl, n_spline_knots=N_MOON_KNOTS,
                base_dir=_infer_base_dir_for_reconstruction(),
                split_zodi=SPLIT_ZODI, n_zodi_spline_knots=N_ZODI_KNOTS,
                telluric=globals().get('TELLURIC_BASIS_KW'))
            _Mbl = np.asarray(_mbl._convolve_matrix_channelwise(_mbl.matrix_oh_stick),
                              dtype=np.float64)
            if _Mbl.shape[0] != _oh_bl.size:
                raise RuntimeError(f'OH stick basis has {_Mbl.shape[0]} rows but '
                                   f'there are {_oh_bl.size} OH coefficients')
            _flux_rmse = {}
            for _nm_bl, _pd_bl in _all_preds.items():
                _D = (np.asarray(_pd_bl, dtype=np.float64)[:, _oh_bl]
                      - _y_te[:, _oh_bl])
                _flux_rmse[_nm_bl] = float(np.sqrt(np.mean((_D @ _Mbl) ** 2)))
            _ml_f = _flux_rmse['ML_default']
            _bl_f = {k: v for k, v in _flux_rmse.items() if k != 'ML_default'}
            _bn_f = min(_bl_f, key=_bl_f.get)
            _gain_f = 100.0 * (_bl_f[_bn_f] - _ml_f) / max(_bl_f[_bn_f], 1e-30)
            _gain_c = (naive_baseline_per_group.get('mesospheric', {})
                       .get('gain_pct', float('nan')))
            print('\n  mesospheric in FLUX space (OH coefficient error projected '
                  'through the convolved stick basis):')
            print(f'    ML={_ml_f:.5g}  best_naive={_bn_f}:{_bl_f[_bn_f]:.5g}  '
                  f'gain={_gain_f:+.1f}%   [coefficient-space gain was '
                  f'{_gain_c:+.1f}%]')
            if np.isfinite(_gain_c) and (_gain_c < 0) and (_gain_f > 0):
                print('    NB the two DISAGREE IN SIGN.  The OH block is internally '
                      'degenerate, so trust\n    the flux number: the coefficient '
                      'metric is scoring a direction the spectrum\n    cannot see.')
            naive_baseline_mesospheric_flux = dict(
                ml=_ml_f, best_naive=_bl_f[_bn_f], best_naive_name=str(_bn_f),
                gain_pct=_gain_f, gain_pct_coef=float(_gain_c))
    except Exception as _exc_bl:
        print(f'  (mesospheric flux-space companion unavailable: '
              f'{type(_exc_bl).__name__}: {_exc_bl})')
else:
    print('  (mesospheric flux-space companion skipped: missing '
          + ', '.join(k for k in _MES_FLUX_KEYS if k not in globals()) + ')')

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


# Persisted so `headline_summary` can synthesise without recomputing.  Named
# without a leading underscore on purpose: the diagnostics cells share one
# exec-globals dict and underscore names are routinely clobbered by later cells.
naive_baseline_result = dict(
    per_group=naive_baseline_per_group,
    # Flux-space mesospheric gain; the coefficient-space one above is degenerate.
    mesospheric_flux=naive_baseline_mesospheric_flux,
    group_equal=dict(ml=_ml_v, best_naive=_best_v, best_naive_name=str(_best),
                     gain_pct=_gain_pct),
    n_test=int(_y_te.shape[0]),
)
