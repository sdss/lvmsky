# Coefficient-space relative residual vs true coefficient value, on log-y.
# Metric = |pred - true| / |true|.  Central line = median of the metric,
# 1σ band = median..p68.27, 2σ band = p68.27..p95.45.  The 68.27 / 95.45
# cuts are the enclosed-mass fractions of a zero-mean Gaussian for
# |z| <= 1 / <= 2, so the 1σ / 2σ labels keep their standard meaning.
# Samples with |true| == 0 are dropped (relative error undefined).  Moon
# and zodi are split into blue / green / red thirds by spline-knot
# wavelength.  Continuum, ionospheric, atomic are one panel each with
# one curve per coefficient.  OH (mesospheric) collapses the 403 coefs
# into a single per-row median|resid| / median|coef| vs per-row median
# |coef|.
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

required = ['coef_pred_det', 'coef_sci_all', 'test_idx',
            'coef_wavelengths_a', '_group_indices_compress',
            'filtered_triplet']
_missing_resid = [k for k in required if k not in globals()]
if _missing_resid:
    raise RuntimeError('Run trainer first. Missing: ' + ', '.join(_missing_resid))

_coef_names_resid = list(filtered_triplet['coef_names'])
_te_resid = np.asarray(test_idx, dtype=int)
_y_true_resid = np.asarray(coef_sci_all[_te_resid], dtype=np.float64)
_y_pred_resid = np.asarray(coef_pred_det, dtype=np.float64)
_resid_all = _y_pred_resid - _y_true_resid
_lambda_A = np.asarray(coef_wavelengths_a, dtype=np.float64)


def _rgba_from_hex(hex_color, alpha):
    _h = hex_color.lstrip('#')
    _r, _g, _b = int(_h[0:2], 16), int(_h[2:4], 16), int(_h[4:6], 16)
    return f'rgba({_r},{_g},{_b},{alpha:.2f})'


# 1σ/2σ enclosed-mass fractions of a zero-mean Gaussian for |z| <= 1 / <= 2.
_GAUSS_ENCLOSED_1S = 68.26894921370858
_GAUSS_ENCLOSED_2S = 95.44997361036416


def _binned_abs_percentiles(x, y, n_bins=18, min_per_bin=15):
    """Median, p68.27, p95.45 of |y| in `n_bins` quantile bins of x."""
    _ok = np.isfinite(x) & np.isfinite(y)
    if _ok.sum() < min_per_bin * 2:
        return None
    x, y = np.asarray(x)[_ok], np.asarray(y)[_ok]
    _abs_y = np.abs(y)
    _edges = np.unique(np.quantile(x, np.linspace(0.0, 1.0, n_bins + 1)))
    if _edges.size < 3:
        return None
    _idx = np.digitize(x, _edges[1:-1], right=False)
    _n_bin = _edges.size - 1
    _out = {k: np.full(_n_bin, np.nan) for k in ('x_ctr', 'p50', 'p1s', 'p2s')}
    _cnt = np.zeros(_n_bin, dtype=int)
    for b in range(_n_bin):
        _m = _idx == b
        _n_b = int(_m.sum())
        if _n_b < min_per_bin:
            continue
        _ays, _xs = _abs_y[_m], x[_m]
        _out['p50'][b] = float(np.nanmedian(_ays))
        _out['p1s'][b] = float(np.nanpercentile(_ays, _GAUSS_ENCLOSED_1S))
        _out['p2s'][b] = float(np.nanpercentile(_ays, _GAUSS_ENCLOSED_2S))
        _out['x_ctr'][b] = float(np.nanmedian(_xs))
        _cnt[b] = _n_b
    _out['count'] = _cnt
    return _out


def _add_percentile_band(fig, row, col, stats, *, color, name,
                         legendgroup, showlegend):
    _valid = (np.isfinite(stats['x_ctr']) & np.isfinite(stats['p50'])
              & (stats['p50'] > 0.0) & (stats['p2s'] > 0.0))
    if not _valid.any():
        return
    _x = stats['x_ctr'][_valid]
    _p50 = stats['p50'][_valid]
    _p1s = stats['p1s'][_valid]
    _p2s = stats['p2s'][_valid]
    fig.add_trace(go.Scatter(
        x=np.concatenate([_x, _x[::-1]]),
        y=np.concatenate([_p2s, _p1s[::-1]]),
        fill='toself', fillcolor=_rgba_from_hex(color, 0.12),
        line=dict(width=0), name=f'{name} 2σ',
        legendgroup=legendgroup, showlegend=False, hoverinfo='skip',
    ), row=row, col=col)
    fig.add_trace(go.Scatter(
        x=np.concatenate([_x, _x[::-1]]),
        y=np.concatenate([_p1s, _p50[::-1]]),
        fill='toself', fillcolor=_rgba_from_hex(color, 0.28),
        line=dict(width=0), name=f'{name} 1σ',
        legendgroup=legendgroup, showlegend=False, hoverinfo='skip',
    ), row=row, col=col)
    fig.add_trace(go.Scatter(
        x=_x, y=_p50, mode='lines',
        line=dict(color=color, width=2.2),
        name=name, legendgroup=legendgroup, showlegend=showlegend,
    ), row=row, col=col)


def _wave_third_slices(gidx, lam):
    """Split coef indices into three equal-wavelength thirds."""
    _lam_g = lam[gidx]
    _fin = np.isfinite(_lam_g)
    if not _fin.any():
        return None
    _lo, _hi = float(_lam_g[_fin].min()), float(_lam_g[_fin].max())
    _edges = np.linspace(_lo, _hi, 4)
    _out = []
    for b in range(3):
        _mask = (_lam_g >= _edges[b]) & (
            _lam_g <= _edges[b + 1] if b == 2 else _lam_g < _edges[b + 1])
        _out.append((gidx[_mask], (_edges[b], _edges[b + 1])))
    return _out


def _choose_xaxis_type(vals, min_span_decades=1.5):
    """Return 'log' if the positive tail spans > 1.5 decades, else 'linear'."""
    _v = np.asarray(vals)
    _v = _v[np.isfinite(_v) & (_v > 0.0)]
    if _v.size < 8:
        return 'linear'
    _p5, _p95 = np.percentile(_v, [5.0, 95.0])
    if _p5 <= 0.0 or _p95 <= 0.0:
        return 'linear'
    return 'log' if (np.log10(_p95 / _p5) > min_span_decades) else 'linear'


_subplot_titles_resid = [
    'Moon — blue third', 'Moon — green third', 'Moon — red third',
    'Zodi — blue third', 'Zodi — green third', 'Zodi — red third',
    'continuum', 'ionospheric', 'atomic',
    'mesospheric (OH lines): per-row median residual vs per-row median |OH coef|',
]
fig_coef_resid = make_subplots(
    rows=4, cols=3,
    subplot_titles=_subplot_titles_resid,
    specs=[
        [{}, {}, {}],
        [{}, {}, {}],
        [{}, {}, {}],
        [{'colspan': 3}, None, None],
    ],
    horizontal_spacing=0.07, vertical_spacing=0.09,
)

_MOON_COLOR = '#d62728'
_ZODI_COLOR = '#ff7f0e'
_OTHER_PALETTE = ['#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#17becf']
_OH_COLOR = '#17becf'

# Moon and zodi: three panels each, pooled over (row x coef) samples per third.
for _row_i, (_gname, _color) in enumerate((('moon', _MOON_COLOR),
                                            ('zodi', _ZODI_COLOR)), start=1):
    _idx_g = np.asarray(_group_indices_compress.get(_gname, []), dtype=int)
    if _idx_g.size == 0:
        continue
    _thirds = _wave_third_slices(_idx_g, _lambda_A)
    if _thirds is None:
        continue
    for _col_i, (_idx_third, (_lo_A, _hi_A)) in enumerate(_thirds, start=1):
        if _idx_third.size == 0:
            continue
        _x_flat = _y_true_resid[:, _idx_third].ravel()
        _y_flat = _resid_all[:, _idx_third].ravel()
        _pos_flat = np.abs(_x_flat) > 0.0
        _x_flat = _x_flat[_pos_flat]
        _rel_flat = np.abs(_y_flat[_pos_flat]) / np.abs(_x_flat)
        _stats = _binned_abs_percentiles(_x_flat, _rel_flat)
        if _stats is None:
            continue
        _label = f'{_gname} {_lo_A:.0f}-{_hi_A:.0f} Å'
        _add_percentile_band(
            fig_coef_resid, row=_row_i, col=_col_i,
            stats=_stats, color=_color, name=_label,
            legendgroup=_gname, showlegend=False,
        )
        _title_idx = 3 * (_row_i - 1) + (_col_i - 1)
        fig_coef_resid.layout.annotations[_title_idx].text = (
            f'{_gname.capitalize()} — {_lo_A:.0f}-{_hi_A:.0f} Å '
            f'(n_coef={_idx_third.size})'
        )
        _xt = _choose_xaxis_type(_x_flat)
        fig_coef_resid.update_xaxes(type=_xt, row=_row_i, col=_col_i)

# Row 3: continuum, ionospheric, atomic — one line per coef inside each panel.
_ROW3_GROUPS = ('continuum', 'ionospheric', 'atomic')
for _col_i, _gname in enumerate(_ROW3_GROUPS, start=1):
    _idx_g = np.asarray(_group_indices_compress.get(_gname, []), dtype=int)
    if _idx_g.size == 0:
        continue
    _x_all_panel = _y_true_resid[:, _idx_g].ravel()
    for _k_i, _k in enumerate(_idx_g):
        _x_col = _y_true_resid[:, _k]
        _y_col = _resid_all[:, _k]
        _pos_col = np.abs(_x_col) > 0.0
        _x_col = _x_col[_pos_col]
        _rel_col = np.abs(_y_col[_pos_col]) / np.abs(_x_col)
        _stats = _binned_abs_percentiles(_x_col, _rel_col)
        if _stats is None:
            continue
        _cname = _coef_names_resid[_k]
        _color = _OTHER_PALETTE[_k_i % len(_OTHER_PALETTE)]
        _add_percentile_band(
            fig_coef_resid, row=3, col=_col_i,
            stats=_stats, color=_color, name=_cname,
            legendgroup=f'{_gname}/{_cname}', showlegend=True,
        )
    _xt = _choose_xaxis_type(_x_all_panel)
    fig_coef_resid.update_xaxes(type=_xt, row=3, col=_col_i)

# Row 4: OH (mesospheric) — per-row median vs per-row median |coef|.
_oh_indices = np.asarray(_group_indices_compress.get('mesospheric', []), dtype=int)
if _oh_indices.size > 0:
    _oh_true = _y_true_resid[:, _oh_indices]
    _oh_resid = _resid_all[:, _oh_indices]
    _oh_med_true = np.nanmedian(np.abs(_oh_true), axis=1)
    _oh_med_abs_resid = np.nanmedian(np.abs(_oh_resid), axis=1)
    _oh_rel = _oh_med_abs_resid / np.where(_oh_med_true > 0.0,
                                            _oh_med_true, np.nan)
    _stats_oh = _binned_abs_percentiles(_oh_med_true, _oh_rel,
                                         n_bins=14, min_per_bin=25)
    if _stats_oh is not None:
        _add_percentile_band(
            fig_coef_resid, row=4, col=1,
            stats=_stats_oh, color=_OH_COLOR, name='OH per-row median',
            legendgroup='oh', showlegend=False,
        )
        _xt_oh = _choose_xaxis_type(_oh_med_true)
        fig_coef_resid.update_xaxes(type=_xt_oh, row=4, col=1)

for _r in range(1, 4):
    for _c in range(1, 4):
        fig_coef_resid.update_xaxes(title_text='coef true value',
                                     row=_r, col=_c)
        fig_coef_resid.update_yaxes(title_text='|pred − true| / |true|',
                                     type='log', row=_r, col=_c)
fig_coef_resid.update_xaxes(
    title_text='per-row median |OH coef|', row=4, col=1)
fig_coef_resid.update_yaxes(
    title_text='median|resid| / median|true|', type='log', row=4, col=1)

fig_coef_resid.update_layout(
    template='plotly_white',
    height=1150, width=1250,
    title=('Coefficient relative residual vs coefficient value '
           '— median |pred − true| / |true| + 1σ (p68.27) + 2σ (p95.45) '
           f'percentile bands, test split n_row={_te_resid.size}'),
    margin=dict(l=60, r=30, t=90, b=60),
    showlegend=True,
    legend=dict(orientation='v', x=1.02, y=1.0, font=dict(size=10)),
)
fig_coef_resid.show()
