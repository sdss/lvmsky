"""Signed O I / N I amplitude error against the context that drives each line, one figure per line.

Same construction as the parent notebook's amplitude-error-vs-context panels (mlp_predictor
diagnostics_cells/amplitude_error_vs_ctx.py): Delta = log10(A_pred / A_true), median per quantile bin of
the Science-pointing context, band = 16-84 percentile range of the rows in the bin.  For a single line
A = c * sum(template), so Delta = log10(c_pred / c_true).  Several methods are overlaid.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

N_BINS = 8
_CLOCKS = ('obstime_day_sin', 'obstime_day_cos', 'obstime_year_sin', 'obstime_year_cos')
_SOLAR = ('sun_alt', 'f107', 'f107_81d', 'kp')
CTX = {  # each line against its own layer's van Rhijn factor, geometry, solar activity and clocks
    'ATOM_Or': ('vanrhijn_285km', 'airmass', 'alt', *_SOLAR, *_CLOCKS),
    'ATOM_Orc_OI0777': ('vanrhijn_285km', 'airmass', 'alt', *_SOLAR, *_CLOCKS),
    'ATOM_Orc_OI0845': ('vanrhijn_285km', 'airmass', 'alt', *_SOLAR, *_CLOCKS),
    'ATOM_N': ('vanrhijn_285km', 'airmass', 'alt', *_SOLAR, *_CLOCKS),
    'ATOM_Og': ('vanrhijn_95km', 'airmass', 'alt', *_SOLAR, *_CLOCKS),
}
TITLE = {'ATOM_Or': '[O I] 6300/6364', 'ATOM_Orc_OI0777': 'O I 7774', 'ATOM_Orc_OI0845': 'O I 8446',
         'ATOM_N': '[N I] 5199', 'ATOM_Og': '[O I] 5577'}
STYLE = {'deployed': ('Deployed ensemble', '#d62728'), 'iono_branches': ('Isolated branches', '#9467bd'),
         'iono_branches_alpha1': ('Isolated branches, alpha = 1', '#ff7f0e'),
         'dedicated': ('Dedicated [O I] MLPs', '#1f77b4')}


def _rgba(hex_color, a):
    h = hex_color.lstrip('#')
    return f'rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})'


def error_vs_ctx(triplet, rows, preds, lines, line, keep=None):
    """Figure for one line; ``preds`` maps method -> (n_rows, len(lines)) predicted coefficients."""
    j = list(lines).index(line)
    names = [str(n) for n in triplet['coef_names']]
    truth = np.asarray(triplet['coef_sci'], np.float64)[rows, names.index(line)]
    ctx_names = [str(n).strip().lower() for n in triplet['ctx_names']]
    ctx = np.asarray(triplet['ctx_sci'], np.float64)[rows]
    keep = np.ones(len(rows), bool) if keep is None else keep
    feats = [f for f in CTX[line] if f in ctx_names]
    nrow = int(np.ceil(len(feats) / 3))
    fig = make_subplots(rows=nrow, cols=3, subplot_titles=feats, vertical_spacing=0.09, horizontal_spacing=0.07)
    summary = []
    for m, (label, color) in STYLE.items():
        if m not in preds or not np.isfinite(preds[m][:, j]).any():
            continue
        p = preds[m][:, j]
        ok = keep & (p > 0) & (truth > 0)
        d = np.where(ok, np.log10(np.where(ok, p, 1) / np.where(ok, truth, 1)), np.nan)
        summary.append(f'{label}: median {np.nanmedian(d):+.4f} dex, MAD {np.nanmedian(np.abs(d)):.4f}')
        for k, f in enumerate(feats):
            x = ctx[:, ctx_names.index(f)]
            ok_x = np.isfinite(x) & np.isfinite(d)
            qs = np.unique(np.nanpercentile(x[ok_x], np.linspace(0, 100, N_BINS + 1)))
            cx, md, lo, hi = [], [], [], []
            for b in range(len(qs) - 1):
                sel = ok_x & (x >= qs[b]) & (x <= qs[b + 1])
                if sel.sum() >= 8:
                    cx.append(float(np.median(x[sel])))
                    md.append(float(np.median(d[sel])))
                    lo.append(float(np.percentile(d[sel], 15.865)))
                    hi.append(float(np.percentile(d[sel], 84.135)))
            r, c = k // 3 + 1, k % 3 + 1
            fig.add_trace(go.Scatter(x=cx + cx[::-1], y=hi + lo[::-1], fill='toself', fillcolor=_rgba(color, 0.1),
                                     line=dict(width=0), mode='lines', showlegend=False, hoverinfo='skip'), r, c)
            fig.add_trace(go.Scatter(x=cx, y=md, mode='lines+markers', name=label, legendgroup=m,
                                     showlegend=k == 0, line=dict(color=color, width=1.6), marker=dict(size=5)),
                          r, c)
    for k in range(len(feats)):
        fig.add_hline(y=0.0, line=dict(color='black', width=0.8, dash='dash'), row=k // 3 + 1, col=k % 3 + 1)
    fig.update_layout(
        height=240 * nrow + 150,
        legend=dict(orientation='h', yanchor='bottom', y=1.03, xanchor='left', x=0.0, font=dict(size=10)),
        title=dict(text=(f'{TITLE[line]} ({line}): signed amplitude error vs its own context '
                         f'(median per quantile bin, {N_BINS} bins; n={int(keep.sum())} rows)<br>'
                         f'<sub>Delta = log10(c_pred/c_true), + = over-predicted. Band: 16-84 pct of rows '
                         f'in the bin. {"; ".join(summary)}</sub>'),
                   font=dict(size=13), x=0.02, xanchor='left'),
        margin=dict(t=150))
    fig.update_yaxes(title_text='Delta [dex]', row=1, col=1)
    return fig
