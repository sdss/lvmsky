"""[O I] 6300/6364 on held-out frames: near-arm and Science decompositions against the predictions.

``select_frames`` picks test frames spread evenly over the lunation (moon phase) and, within each phase bin,
over time from twilight (Sun altitude) and airmass.  ``frames_table`` lists them; ``frames_figure`` draws
their amplitudes; ``correlation_figure`` scatters Science [O I] against the near arm and each prediction
on all clean test rows, with the selected frames marked.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from astropy.io import fits
from plotly.subplots import make_subplots

from mlp_predictor.ml_utils import moon_phase_deg_from_ctx

LINE = 'ATOM_Or'
STYLE = {'deployed': ('Deployed ensemble', '#d62728', 'circle'),
         'iono_branches': ('Isolated branches', '#9467bd', 'diamond'),
         'dedicated': ('Dedicated [O I] MLP', '#1f77b4', 'square')}


def _ctx(triplet, rows, name):
    names = [str(n).strip().lower() for n in triplet['ctx_names']]
    return np.asarray(triplet['ctx_sci'], np.float64)[rows, names.index(name)]


def _coef(triplet, rows, arm):
    names = [str(n) for n in triplet['coef_names']]
    return np.asarray(triplet[f'coef_{arm}'], np.float64)[rows, names.index(LINE)]


def select_frames(triplet, rows, keep, n_phase=8, per_bin=3):
    """Positions (into ``rows``) of ``n_phase * per_bin`` clean frames, one night each where possible."""
    phase = moon_phase_deg_from_ctx(triplet)[rows]
    sun, air = _ctx(triplet, rows, 'sun_alt'), _ctx(triplet, rows, 'airmass')
    night = np.floor(np.asarray(triplet['obstime_mjd'], np.float64)[rows] - 0.5).astype(int)
    edges = np.linspace(0, 360, n_phase + 1)
    picked, used = [], set()
    for b in range(n_phase):
        cand = np.flatnonzero(keep & (phase >= edges[b]) & (phase < edges[b + 1]))
        if cand.size == 0:
            continue
        # spread over Sun altitude, then prefer a mid-range airmass and an unused night
        order = cand[np.argsort(sun[cand])]
        for chunk in np.array_split(order, min(per_bin, order.size)):
            fresh = [i for i in chunk if night[i] not in used] or list(chunk)
            i = min(fresh, key=lambda k: abs(air[k] - np.median(air[cand])))
            picked.append(i)
            used.add(night[i])
    return np.array(picked)


def frames_table(triplet, rows, preds, lines, frames, corpus_fits):
    """One row per selected frame: identity, context and the [O I] 6300 amplitudes (internal 1e14 units)."""
    j = list(lines).index(LINE)
    r = np.asarray(rows)[frames]
    with fits.open(corpus_fits, memmap=True) as h:
        meta = h['META'].data
        col = next(c for c in meta.columns.names if c.upper() in ('EXPNUM', 'EXP_NUM', 'EXPOSURE'))
        expnum = np.asarray(meta[col])[np.asarray(triplet['row_index'], int)[r]]
    tab = pd.DataFrame({
        'expnum': expnum.astype(int),
        'mjd': np.asarray(triplet['obstime_mjd'], np.float64)[r].round(3),
        'moon_phase_deg': moon_phase_deg_from_ctx(triplet)[r].round(0),
        'moon_fli': _ctx(triplet, r, 'moon_fli').round(2),
        'sun_alt': _ctx(triplet, r, 'sun_alt').round(1),
        'airmass': _ctx(triplet, r, 'airmass').round(2),
        'near': _coef(triplet, r, 'near'), 'far': _coef(triplet, r, 'far'), 'sci (truth)': _coef(triplet, r, 'sci'),
    })
    for m in STYLE:
        if m in preds:
            tab[m] = preds[m][frames, j]
    return tab.round(2)


def frames_figure(tab):
    """[O I] 6300 amplitude per selected frame: near/far decompositions, Science truth, predictions."""
    x = [f"{e}<br>{p:.0f}°" for e, p in zip(tab['expnum'], tab['moon_phase_deg'])]
    truth = tab['sci (truth)'].to_numpy()
    fig = make_subplots(2, 1, shared_xaxes=True, row_heights=[0.65, 0.35], vertical_spacing=0.05)
    fig.add_trace(go.Scatter(x=x, y=tab['near'], mode='markers', name='Near arm (decomposition)',
                             marker=dict(color='gray', symbol='triangle-up', size=9)), 1, 1)
    fig.add_trace(go.Scatter(x=x, y=tab['far'], mode='markers', name='Far arm (decomposition)',
                             marker=dict(color='lightgray', symbol='triangle-down', size=9)), 1, 1)
    fig.add_trace(go.Scatter(x=x, y=truth, mode='markers', name='Science (decomposition, truth)',
                             marker=dict(color='black', symbol='x-thin', size=12, line=dict(width=2))), 1, 1)
    for m, (label, color, symbol) in STYLE.items():
        if m in tab:
            fig.add_trace(go.Scatter(x=x, y=tab[m], mode='markers', name=label,
                                     marker=dict(color=color, symbol=symbol, size=7)), 1, 1)
            fig.add_trace(go.Scatter(x=x, y=100 * (tab[m] - truth) / truth, mode='markers', showlegend=False,
                                     marker=dict(color=color, symbol=symbol, size=7)), 2, 1)
    fig.add_trace(go.Scatter(x=x, y=100 * (tab['near'] - truth) / truth, mode='markers', showlegend=False,
                             marker=dict(color='gray', symbol='triangle-up', size=8)), 2, 1)
    fig.add_hline(y=0, line=dict(color='black', width=0.8, dash='dash'), row=2, col=1)
    fig.update_yaxes(type='log', title_text='[O I] 6300 amplitude (1e-14)', row=1, col=1)
    fig.update_yaxes(title_text='(x - truth) / truth [%]', row=2, col=1)
    fig.update_xaxes(title_text='exposure, moon phase (frames ordered by lunation)', tickangle=0,
                     tickfont=dict(size=9), row=2, col=1)
    fig.update_layout(height=700, width=1500, legend=dict(orientation='h', y=1.08),
                      title='[O I] 6300/6364 on selected held-out frames: sky arms, Science truth and predictions')
    return fig


def correlation_figure(triplet, rows, preds, lines, keep, frames):
    """Science [O I] 6300 truth against the near arm and each prediction, all clean test rows."""
    j = list(lines).index(LINE)
    truth = _coef(triplet, rows, 'sci')
    xs = {'near': ('Near arm (decomposition)', 'gray', _coef(triplet, rows, 'near'))}
    xs.update({m: (label, color, preds[m][:, j]) for m, (label, color, _) in STYLE.items() if m in preds})
    fig = make_subplots(1, len(xs), shared_yaxes=True, horizontal_spacing=0.03,
                        subplot_titles=[v[0] for v in xs.values()])
    ok_all = keep & (truth > 0)
    lo, hi = np.percentile(truth[ok_all], [0.5, 99.5])
    sel = np.zeros(len(rows), bool)
    sel[frames] = True
    for c, (label, color, x) in enumerate(xs.values(), start=1):
        ok = ok_all & (x > 0)
        d = np.log10(x[ok] / truth[ok])
        r = np.corrcoef(np.log10(x[ok]), np.log10(truth[ok]))[0, 1]
        fig.add_trace(go.Scattergl(x=x[ok & ~sel], y=truth[ok & ~sel], mode='markers', showlegend=False,
                                   marker=dict(size=3, color=color, opacity=0.35)), 1, c)
        fig.add_trace(go.Scatter(x=x[ok & sel], y=truth[ok & sel], mode='markers', showlegend=c == 1,
                                 name='selected frames', marker=dict(size=8, color='black', symbol='x-thin',
                                                                    line=dict(width=1.5))), 1, c)
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode='lines', showlegend=False,
                                 line=dict(color='black', dash='dash', width=1)), 1, c)
        fig.add_annotation(text=f'r(log) = {r:.4f}<br>median |Δ| = {np.median(np.abs(d)):.4f} dex<br>'
                                f'bias = {np.median(d):+.4f} dex<br>n = {ok.sum()}',
                           xref=f'x{c if c > 1 else ""} domain', yref='y domain', x=0.03, y=0.97,
                           showarrow=False, align='left', font=dict(size=10))
        fig.update_xaxes(type='log', range=np.log10([lo, hi]), title_text='x (1e-14)', row=1, col=c)
    fig.update_yaxes(type='log', range=np.log10([lo, hi]))
    fig.update_yaxes(title_text='Science [O I] 6300, decomposition (1e-14)', row=1, col=1)
    fig.update_layout(height=480, width=1500,
                      title='Science [O I] 6300/6364 against the near arm and the predictions (clean test rows)')
    return fig
