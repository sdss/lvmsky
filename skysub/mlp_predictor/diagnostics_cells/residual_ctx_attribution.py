# Residual -> context regression per group (missing-feature detector).
#
# For each group g:
#   1. Compute per-row per-group residual RMSE on the test split:
#          r_g(row) = sqrt(mean_{j in g} (pred[row, j] - true[row, j])**2)
#   2. Fit a linear regression on r_g(row) ~ ctx_sci features.
#      Report 5-fold CV R^2 and the top-|coef| features (standardised X).
#   3. Fit a random forest for non-linear structure; report CV R^2 and top
#      permutation importances.
#
# What the numbers mean:
#   linear R^2 > 0.1  -> ctx features carry recoverable signal the head is
#                        not modelling; broaden its input set or capacity.
#   RF R^2 - lin R^2 > 0.1 -> the residual is non-linearly related to
#                        an existing feature; the head is under-parameterised
#                        for that variable.
#   both near 0        -> residual is either irreducible physics (gravity
#                        waves, per-night atmospheric variability) or living
#                        in a non-context space (row-to-row time correlation,
#                        input noise).
#
# Requires sklearn (usually available in this env).

from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler

required = ['filtered_triplet', '_group_indices_compress',
            'coef_sci_all', 'ctx_sci_all', 'test_idx', 'coef_pred_det']
_missing = [k for k in required if k not in globals()]
if _missing:
    raise RuntimeError('Run prerequisite cells first. Missing: ' + ', '.join(_missing))

_te = np.asarray(test_idx, dtype=int)
_ctx_names = list(filtered_triplet['ctx_names'])
_X_full = np.asarray(ctx_sci_all[_te], dtype=np.float64)
_y_true = np.asarray(coef_sci_all[_te], dtype=np.float64)
_y_pred = np.asarray(coef_pred_det, dtype=np.float64)

# Drop constant / near-constant columns (would trip the standardiser).
_std_x = np.nanstd(_X_full, axis=0)
_keep_cols = np.where(_std_x > 1e-12)[0]
if len(_keep_cols) < len(_std_x):
    print(f'  dropping {len(_std_x) - len(_keep_cols)} constant ctx columns')
_X = _X_full[:, _keep_cols]
_ctx_names_k = [_ctx_names[i] for i in _keep_cols]

_kf = KFold(n_splits=5, shuffle=True, random_state=42)


def _cv_r2(_model, X, y):
    return float(np.mean(cross_val_score(_model, X, y, cv=_kf, scoring='r2', n_jobs=1)))


print('=' * 90)
print('Residual -> context regression per group (test split, 5-fold CV R^2)')
print('=' * 90)

_scaler = StandardScaler()
_X_std = _scaler.fit_transform(_X)

_summary = []
for gname, idx in _group_indices_compress.items():
    idx = np.asarray(idx, dtype=int)
    _r_g = np.sqrt(np.mean((_y_pred[:, idx] - _y_true[:, idx]) ** 2, axis=1))
    _finite = np.isfinite(_r_g)
    if _finite.sum() < 20:
        print(f'  {gname}: too few finite residuals; skipping')
        continue
    _yy = _r_g[_finite]
    _XX = _X_std[_finite]

    _lin = LinearRegression()
    _lin_r2 = _cv_r2(_lin, _XX, _yy)
    _lin.fit(_XX, _yy)
    _lin_coefs = _lin.coef_
    _lin_top = np.argsort(np.abs(_lin_coefs))[::-1][:5]

    _rf = RandomForestRegressor(n_estimators=200, max_depth=8, random_state=42, n_jobs=-1)
    _rf_r2 = _cv_r2(_rf, _XX, _yy)
    _rf.fit(_XX, _yy)
    _rf_top = np.argsort(_rf.feature_importances_)[::-1][:5]

    print()
    print(f'--- group: {gname}  (n_valid_rows = {int(_finite.sum())})  '
          f'residual mean = {float(np.mean(_yy)):.4g}  median = {float(np.median(_yy)):.4g}')
    print(f'  LinearRegression  CV R^2 = {_lin_r2:+.3f}')
    print('    top |coef| (standardised X):')
    for _i in _lin_top:
        print(f'      {_ctx_names_k[_i]:<22s} coef = {_lin_coefs[_i]:+.3g}')
    print(f'  RandomForest      CV R^2 = {_rf_r2:+.3f}   '
          f'(non-linear gain = {_rf_r2 - _lin_r2:+.3f})')
    print('    top importances:')
    for _i in _rf_top:
        print(f'      {_ctx_names_k[_i]:<22s} imp  = {_rf.feature_importances_[_i]:.3g}')

    _summary.append({
        'group': gname,
        'n': int(_finite.sum()),
        'residual_mean': float(np.mean(_yy)),
        'lin_R2': _lin_r2,
        'rf_R2':  _rf_r2,
        'nonlin_gain': _rf_r2 - _lin_r2,
    })

if _summary:
    print()
    print('Summary:')
    print(pd.DataFrame(_summary).to_string(index=False,
          float_format=lambda v: f'{v:.4g}'))
    print()
    print('Interpretation:')
    print('  lin_R2 > 0.10        -> ctx features linearly predict residual; head '
          'is under-using them.  Retrain with a wider trunk or head, or '
          'reconsider the compressor for that group.')
    print('  nonlin_gain > 0.10   -> RF finds structure that the linear model '
          'does not; head is under-parameterised for that feature (or the '
          'feature enters non-linearly, e.g. sign or threshold).')
    print('  both < 0.02          -> residual is irreducible from ctx.  Look at '
          'sky_arm_disagreement_floor() to check the noise floor and '
          'wavelength_residual_atlas() to see where the failure lives.')
