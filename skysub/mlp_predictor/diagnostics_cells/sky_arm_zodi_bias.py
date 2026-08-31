# Sky-arm zodi bias diagnostic per slice.
# Answers: is the residual Q1/Q4 zodi bias an ML tuning issue or a decomp-truth
# artifact? Compares geometry-corrected sky1 (B1) and sky2 (B_far) predictions
# vs sci truth on the same slices used by the ML per-slice metrics cell.
#   - If bias_B1 / bias_B2 in a slice matches bias_ML in sign and magnitude
#     -> ML is faithfully reproducing what the sky arms tell it; residual comes
#        from the decomp truth itself. No ML tuning can fix that.
#   - If bias_B1 / bias_B2 are near zero but bias_ML is high in the same slice
#     -> ML is genuinely over- or under-predicting; row-weight boost or
#        regime lift extension can help.

import numpy as np
import pandas as pd

required = ['filtered_triplet', 'test_idx',
            'coef_near_all', 'coef_far_all', 'coef_sci_all',
            'ctx_near_all', 'ctx_far_all', 'ctx_sci_all',
            'coef_pred_det', '_group_indices_compress',
            'compress_geom_kwargs', 'airglow_geometry_scale',
            '_moon_phase_deg_from_ctx']
_missing = [k for k in required if k not in globals()]
if _missing:
    raise RuntimeError('Run prerequisite cells first. Missing: ' + ', '.join(_missing))

_te = np.asarray(test_idx, dtype=int)
_ctx_names = list(filtered_triplet['ctx_names'])

_ctx_near_te = np.asarray(ctx_near_all[_te], dtype=np.float64)
_ctx_far_te  = np.asarray(ctx_far_all[_te],  dtype=np.float64)
_ctx_sci_te  = np.asarray(ctx_sci_all[_te],  dtype=np.float64)

_coef_near_te = np.asarray(coef_near_all[_te], dtype=np.float64)
_coef_far_te  = np.asarray(coef_far_all[_te],  dtype=np.float64)
_coef_sci_te  = np.asarray(coef_sci_all[_te],  dtype=np.float64)
_coef_ml_te   = np.asarray(coef_pred_det,       dtype=np.float64)

# Geometry-corrected sky arms projected onto sci geometry (matches B1/B2 in cell 34).
_G_near = airglow_geometry_scale(_ctx_near_te, **compress_geom_kwargs)
_G_far  = airglow_geometry_scale(_ctx_far_te,  **compress_geom_kwargs)
_G_sci  = airglow_geometry_scale(_ctx_sci_te,  **compress_geom_kwargs)
_em_near_sci = (_coef_near_te / _G_near) * _G_sci
_em_far_sci  = (_coef_far_te  / _G_far)  * _G_sci
_em_mean_sci = 0.5 * (_em_near_sci + _em_far_sci)

_zodi_idx = np.asarray(_group_indices_compress['zodi'], dtype=int)

_moon_phase_te = _moon_phase_deg_from_ctx(
    {'ctx_sci': _ctx_sci_te, 'ctx_names': _ctx_names})
_moon_alt_te = _ctx_sci_te[:, _ctx_names.index('moon_alt')]
_airmass_te = _ctx_sci_te[:, _ctx_names.index('airmass')]

_slices = {
    'all_test':                     np.ones(len(_te), dtype=bool),
    'moon_phase Q1 [2-110 deg]':    (_moon_phase_te >= 2)   & (_moon_phase_te < 110),
    'moon_phase Q2 [110-172 deg]':  (_moon_phase_te >= 110) & (_moon_phase_te < 172),
    'moon_phase Q3 [172-265 deg]':  (_moon_phase_te >= 172) & (_moon_phase_te < 265),
    'moon_phase Q4 [265-358 deg]':  (_moon_phase_te >= 265) & (_moon_phase_te < 358),
    'moon_alt > 0 (moon up)':       _moon_alt_te > 0,
    'moon_alt <= 0 (moon down)':    _moon_alt_te <= 0,
    'airmass <= 1.5 (low)':         _airmass_te <= 1.5,
    'airmass > 1.5 (high)':         _airmass_te > 1.5,
}


def _zodi_bias_pct(pred, truth, mask, idx):
    _t = float(np.mean(truth[mask][:, idx]))
    _p = float(np.mean(pred[mask][:, idx]))
    return 100.0 * (_p - _t) / max(abs(_t), 1e-30)


_rows = []
for _name, _mask in _slices.items():
    _n = int(_mask.sum())
    if _n < 20:
        continue
    _rows.append({
        'slice': _name,
        'n_rows': _n,
        'zodi_sci_mean': float(np.mean(_coef_sci_te[_mask][:, _zodi_idx])),
        'bias_B1_sky1_geo_%': _zodi_bias_pct(_em_near_sci, _coef_sci_te, _mask, _zodi_idx),
        'bias_B_sky2_geo_%':  _zodi_bias_pct(_em_far_sci,  _coef_sci_te, _mask, _zodi_idx),
        'bias_B_mean_geo_%':  _zodi_bias_pct(_em_mean_sci, _coef_sci_te, _mask, _zodi_idx),
        'bias_ML_default_%':  _zodi_bias_pct(_coef_ml_te,  _coef_sci_te, _mask, _zodi_idx),
    })
_sky_zodi_bias_df = pd.DataFrame(_rows)

print('=' * 110)
print(f'Sky-arm zodi bias diagnostic per slice (test rows, n={len(_te)}).')
print('Compares geometry-corrected sky1/sky2/mean and ML against sci truth on the zodi group.')
print('=' * 110)
print(_sky_zodi_bias_df.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
print()

print('-' * 110)
print('Per-slice verdict (delta = ML_bias - B_mean_bias): |delta| > 1.5%% flagged as ML-tuning-fixable.')
print('-' * 110)
for _, _r in _sky_zodi_bias_df.iterrows():
    _delta = _r['bias_ML_default_%'] - _r['bias_B_mean_geo_%']
    if abs(_delta) > 1.5:
        _tag = ' <-- ML residual is NOT truth-driven; tuning can help'
    elif abs(_delta) < 0.5:
        _tag = ' <-- ML bias tracks decomp truth; tuning cannot help'
    else:
        _tag = ' (borderline)'
    print(f"  {_r['slice']:<32} B_mean={_r['bias_B_mean_geo_%']:+6.2f}%  "
          f"ML={_r['bias_ML_default_%']:+6.2f}%  delta={_delta:+6.2f}%{_tag}")
