"""Plotly diagnostics for the [O I] experiment notebook; each function returns a figure or a table."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from astropy.io import fits
from plotly.subplots import make_subplots

from .train import log10

METHOD_STYLE = {
    'near_copy': ('Near-arm copy', 'gray'),
    'sky_geomean': ('Sky-arm geometric mean', 'darkgoldenrod'),
    'deployed_ensemble': ('Deployed ensemble', 'crimson'),
    'oi_mlp_or': ('[O I] MLP: 6300/6364 only', 'royalblue'),
    'oi_mlp_orc': ('[O I] MLP: 7774+8446 only', 'royalblue'),
    'oi_mlp_og': ('[O I] MLP: 5577 only', 'royalblue'),
    'oi_mlp_all': ('[O I] MLP: all O I lines', 'seagreen'),
    'oi_mlp_groups': ('[O I] MLPs per mechanism (or, orc, og)', 'royalblue'),
    'split_heads': ('Split ionospheric heads', 'darkorange'),
    'control': ('Parent model retrained (control)', 'gray'),
    'iono_branches': ('Isolated ionospheric + 5577 branches', 'purple'),
    'iono_branches_alpha1': ('Isolated branches, alpha = 1', 'seagreen'),
}
CORE_HALF_WIDTH = 3.0


def journal_style():
    """Journal-style Plotly axes, as in the deployed training notebook."""
    for ax in (pio.templates['plotly_white'].layout.xaxis, pio.templates['plotly_white'].layout.yaxis):
        ax.update(showgrid=False, showline=True, mirror=True, linecolor='black', linewidth=1,
                  ticks='inside', zeroline=False)
    pio.templates.default = 'plotly_white'


def or_ratio(t):
    i = list(t['oi_names']).index('ATOM_Or')
    return log10(t['coef_sci'][:, i]) - 0.5 * (log10(t['coef_near'][:, i]) + log10(t['coef_far'][:, i]))


def outlier_cut(t, clean, outlier_dex):
    r = or_ratio(t)
    i = list(t['oi_names']).index('ATOM_Or')
    sky = 0.5 * (log10(t['coef_near'][:, i]) + log10(t['coef_far'][:, i]))
    fig = make_subplots(1, 2, subplot_titles=(
        'Science / sky-arm ATOM_Or', 'Science vs sky-arm ATOM_Or'))
    fig.add_trace(go.Histogram(x=r, nbinsx=200, marker_color='gray', name='All rows'), 1, 1)
    for s in (-outlier_dex, outlier_dex):
        fig.add_vline(x=s, line_dash='dash', line_color='crimson', row=1, col=1)
    for keep, name, color in ((True, 'Kept', 'black'), (False, 'Removed', 'crimson')):
        m = clean == keep
        fig.add_trace(go.Scattergl(x=sky[m], y=sky[m] + r[m], mode='markers', name=f'{name} ({m.sum()})',
                                   marker=dict(size=3, color=color, opacity=0.5)), 1, 2)
    fig.update_xaxes(title_text='log10(sci / geometric mean of sky arms)', row=1, col=1)
    fig.update_yaxes(title_text='Rows', type='log', row=1, col=1)
    fig.update_xaxes(title_text='log10 sky-arm ATOM_Or', row=1, col=2)
    fig.update_yaxes(title_text='log10 science ATOM_Or', row=1, col=2)
    fig.update_layout(height=420, width=1100)
    return fig


def oi_correlations(t, arm='near'):
    names = list(t['oi_names'])
    rho = pd.DataFrame(log10(t[f'coef_{arm}']), columns=names).corr(method='spearman')
    fig = go.Figure(go.Heatmap(z=rho.values, x=names, y=names, zmin=0, zmax=1, colorscale='Viridis',
                               text=np.round(rho.values, 2), texttemplate='%{text}'))
    fig.update_layout(title=f'Spearman correlation of O I coefficients, sky {arm} arm',
                      height=420, width=520)
    return fig


def results_table(summary, subset='test_clean'):
    rows = []
    for line, methods in summary['results'][subset].items():
        for method, s in methods.items():
            rows.append(dict(line=line, method=METHOD_STYLE[method][0], **s))
    return pd.DataFrame(rows).set_index(['line', 'method']).round(4)


def pred_vs_true(t, preds, rows, line='ATOM_Or', methods=('deployed_ensemble', 'oi_mlp_or')):
    k = list(t['oi_names']).index(line)
    truth = log10(t['coef_sci'][rows, k])
    fig = make_subplots(1, len(methods), subplot_titles=[METHOD_STYLE[m][0] for m in methods],
                        shared_yaxes=True)
    lo, hi = np.percentile(truth, [0.5, 99.5])
    for j, m in enumerate(methods, 1):
        fig.add_trace(go.Scattergl(x=truth, y=preds[m][rows, k], mode='markers', showlegend=False,
                                   marker=dict(size=3, color=METHOD_STYLE[m][1], opacity=0.5)), 1, j)
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode='lines', showlegend=False,
                                 line=dict(color='black', dash='dash', width=1)), 1, j)
        fig.update_xaxes(title_text=f'log10 science {line} (decomposition)', row=1, col=j)
    fig.update_yaxes(title_text=f'log10 predicted {line}', row=1, col=1)
    fig.update_layout(height=430, width=1000)
    return fig


def error_vs_context(t, preds, rows, line='ATOM_Or',
                     features=('sun_alt', 'vanrhijn_285km', 'f107', 'kp'),
                     methods=('near_copy', 'deployed_ensemble', 'oi_mlp_or', 'oi_mlp_all'), n_bins=12):
    k = list(t['oi_names']).index(line)
    names = list(t['ctx_names'])
    truth = log10(t['coef_sci'][rows, k])
    fig = make_subplots(1, len(features), shared_yaxes=True, subplot_titles=features)
    for j, f in enumerate(features, 1):
        x = t['ctx_sci'][rows, names.index(f)]
        edges = np.unique(np.percentile(x, np.linspace(0, 100, n_bins + 1)))
        centre = 0.5 * (edges[1:] + edges[:-1])
        which = np.clip(np.digitize(x, edges) - 1, 0, len(centre) - 1)
        for m in methods:
            err = np.abs(preds[m][rows, k] - truth)
            med = [np.median(err[which == b]) if np.any(which == b) else np.nan for b in range(len(centre))]
            fig.add_trace(go.Scatter(x=centre, y=med, mode='lines+markers', name=METHOD_STYLE[m][0],
                                     legendgroup=m, showlegend=j == 1,
                                     line=dict(color=METHOD_STYLE[m][1],
                                               dash='dash' if m == 'oi_mlp_all' else 'solid')), 1, j)
        fig.update_xaxes(title_text=f, row=1, col=j)
    fig.update_yaxes(title_text=f'median |log10 error| {line} (dex)', row=1, col=1)
    fig.update_layout(height=420, width=1300, legend=dict(orientation='h', y=-0.25))
    return fig


def split_overlap_check(t, preds, line='ATOM_Og', seen_mad=0.0025):
    """Per-night deployed-ensemble error on each split: nights it trained on show a much smaller error."""
    k = list(t['oi_names']).index(line)
    err = np.abs(preds['deployed_ensemble'][:, k] - log10(t['coef_sci'][:, k]))
    night = np.floor(t['obstime_mjd'] - 0.5).astype(int)
    rows = []
    for split in ('train', 'val', 'test'):
        idx = t[f'{split}_idx']
        mad = np.array([np.median(err[idx][night[idx] == n]) for n in np.unique(night[idx])])
        p10, p50, p90 = np.percentile(mad, [10, 50, 90])
        rows.append({'split': split, 'nights': mad.size, 'per-night MAD p10': p10, 'p50': p50, 'p90': p90,
                     f'nights with MAD < {seen_mad}': np.mean(mad < seen_mad)})
    return pd.DataFrame(rows).set_index('split').round(4)


def every10_rows(every10_stack, expnum_by_row, rows):
    """every10 row numbers and table rows for the requested table rows present in every10."""
    with fits.open(every10_stack, memmap=True) as h:
        e10 = np.asarray(h['META'].data['expnum'])
    where = {int(e): j for j, e in enumerate(e10)}
    pairs = [(where[int(expnum_by_row[r])], r) for r in rows if int(expnum_by_row[r]) in where]
    return np.array([p[0] for p in pairs], int), np.array([p[1] for p in pairs], int)


def oi_window(every10_stack, every10_sci_decomp, t, rows, lo=6290.0, hi=6312.0):
    """Science spectrum around [O I] 6300 for the table rows that are in every10.

    Returns wave, observed, the decomposition's own [O I] line (COMP_ATOM, the only atomic
    line here), the rest of its sky model, and the matching table rows.
    """
    e10, rows = every10_rows(every10_stack, t['expnum'], rows)
    with fits.open(every10_stack, memmap=True) as hs, fits.open(every10_sci_decomp, memmap=True) as hd:
        wave = np.asarray(hs['WAVE'].data, float)
        w = (wave >= lo) & (wave <= hi)
        observed = np.asarray(hs['FLUX_SCI'].data[e10][:, w], float) * 1e14
        oi = np.asarray(hd['COMP_ATOM'].data[e10][:, w], float)
        # The COMP_* components sum to BESTFIT_LSF.
        rest = np.asarray(hd['BESTFIT_LSF'].data[e10][:, w], float) - oi
    return wave[w], observed, oi, rest, rows


def oi_scales(t, preds, rows, methods):
    """Predicted / decomposed science ATOM_Or per row: how much each method scales the true line."""
    k = list(t['oi_names']).index('ATOM_Or')
    truth = np.maximum(t['coef_sci'][rows, k], 1e-4)
    return {m: 10 ** preds[m][rows, k] / truth for m in methods}


def oi_example_spectra(window, t, preds, methods=('deployed_ensemble', 'oi_mlp_or'), percentiles=(50, 90, 99)):
    """Top: observed [O I] 6300 and the line each method subtracts. Bottom: what is left after subtraction."""
    wave, observed, oi, rest, rows = window
    scale = oi_scales(t, preds, rows, methods)
    brightness = oi.sum(axis=1)
    picks = [int(np.argmin(np.abs(brightness - np.percentile(brightness, q)))) for q in percentiles]
    titles = [f'Exposure {int(t["expnum"][rows[j]])} ([O I] brightness p{q})'
              for j, q in zip(picks, percentiles, strict=True)]
    fig = make_subplots(2, len(picks), subplot_titles=titles + [''] * len(picks), shared_xaxes=True,
                        row_heights=(0.6, 0.4), vertical_spacing=0.06)
    lines = {'[O I] fitted on the science spectrum (best possible)': (np.ones(len(rows)), 'gray', 'dot')}
    lines.update({f'[O I] predicted: {METHOD_STYLE[m][0]}': (scale[m], METHOD_STYLE[m][1], 'solid') for m in methods})
    for col, j in enumerate(picks, 1):
        line = observed[j] - rest[j]
        fig.add_trace(go.Scatter(x=wave, y=line, name='Observed [O I] (observed minus sky model without [O I])',
                                 legendgroup='obs', showlegend=col == 1, line=dict(color='black')), 1, col)
        for name, (s, color, dash) in lines.items():
            style = dict(color=color, dash=dash)
            fig.add_trace(go.Scatter(x=wave, y=s[j] * oi[j], name=name, legendgroup=name, showlegend=col == 1,
                                     line=style), 1, col)
            fig.add_trace(go.Scatter(x=wave, y=line - s[j] * oi[j], name=name, legendgroup=name,
                                     showlegend=False, line=style), 2, col)
        fig.add_hline(y=0, line_color='black', line_width=1, row=2, col=col)
        fig.update_xaxes(title_text='Air wavelength (Angstrom)', row=2, col=col)
    fig.update_yaxes(title_text='Flux (1e-14 erg/s/cm2/A)', row=1, col=1)
    fig.update_yaxes(title_text='Left after subtraction', row=2, col=1)
    fig.update_layout(height=650, width=1300, legend=dict(orientation='h', y=-0.15))
    return fig


def oi_leftover(window, t, preds, methods=('deployed_ensemble', 'oi_mlp_or', 'oi_mlp_all')):
    """[O I] 6300 flux left after subtraction, in % of the line flux, within +/-3 A of the line.

    Positive: the sky line was under-subtracted; negative: over-subtracted.
    """
    wave, observed, oi, rest, rows = window
    core = np.abs(wave - 6300.30) <= CORE_HALF_WIDTH
    line = oi[:, core].sum(axis=1)
    scale = oi_scales(t, preds, rows, methods)
    left = {'Fitted on the science spectrum (best possible)': observed - rest - oi}
    left.update({METHOD_STYLE[m][0]: observed - rest - scale[m][:, None] * oi for m in methods})
    pct = {name: 100 * r[:, core].sum(axis=1) / line for name, r in left.items()}
    colors = ['gray'] + [METHOD_STYLE[m][1] for m in methods]
    fig = make_subplots(1, 2, subplot_titles=(
        f'Leftover per exposure ({len(rows)} exposures)', 'Typical size of the leftover'))
    for (name, v), c in zip(pct.items(), colors, strict=True):
        fig.add_trace(go.Box(y=v, name=name, marker_color=c, boxpoints=False, legendgroup=name), 1, 1)
        fig.add_trace(go.Bar(x=[name], y=[np.median(np.abs(v))], marker_color=c, showlegend=False,
                             text=[f'{np.median(np.abs(v)):.1f}%'], textposition='outside'), 1, 2)
    fig.add_hline(y=0, line_color='black', line_width=1, row=1, col=1)
    fig.update_yaxes(title_text='[O I] 6300 left after subtraction (% of line flux)', row=1, col=1)
    fig.update_yaxes(title_text='median |leftover| (% of line flux)', row=1, col=2)
    fig.update_xaxes(showticklabels=False)
    fig.update_layout(height=480, width=1200, legend=dict(orientation='h', y=-0.15))
    table = pd.DataFrame({name: dict(median_pct=np.median(v), median_abs_pct=np.median(np.abs(v)),
                                     p90_abs_pct=np.percentile(np.abs(v), 90)) for name, v in pct.items()}).T
    return fig, table.round(2)


EXPOSURE_LINES = {'[O I] 5577': 5577.34, '[O I] 6300': 6300.30, '[O I] 6364': 6363.78,
                  'O I 7774': 7774.2, 'O I 8446': 8446.4}


def exposure_coefficients(paths, oi_names):
    """O I amplitudes of one exposure: the three median fits and the ML values each product subtracted."""
    rows = {}
    with fits.open(next(iter(paths.values()))) as h:
        labels = [str(x).strip() for x in h['MED_LABELS'].data[h['MED_LABELS'].columns.names[0]]]
        names = [str(n).strip() for n in h['ML_COEF'].data['name']]
        for label, coef in zip(labels, h['MED_COEF'].data, strict=True):
            rows[f'{label} median fit'] = coef
    for method, path in paths.items():
        with fits.open(path) as h:
            rows[f'{METHOD_STYLE[method][0]} (subtracted)'] = h['ML_COEF'].data['coef']
    return pd.DataFrame({k: [float(v[names.index(n)]) for n in oi_names] for k, v in rows.items()},
                        index=list(oi_names)).T


def _science_rows(h):
    slit = h['SLITMAP'].data
    tel = np.char.strip(np.asarray(slit['telescope']).astype(str))
    return np.flatnonzero((tel == 'Sci') & np.asarray(slit['skyml_applied'], bool))


def exposure_oi_spectra(paths, half_width=12.0):
    """Median over Science fibers of the delivered FLUX (top) and of SKY(method) - SKY(first method) (bottom)."""
    cube = {}
    for method, path in paths.items():
        with fits.open(path, memmap=True) as h:
            wave = np.asarray(h['WAVE'].data, float)
            rows = _science_rows(h)
            cols = {k: np.flatnonzero(np.abs(wave - w) <= half_width) for k, w in EXPOSURE_LINES.items()}
            pick = np.concatenate(list(cols.values()))
            flux = np.asarray(h['FLUX'].data[rows][:, pick], float) * 1e14
            ivar = np.asarray(h['IVAR'].data[rows][:, pick], float) / 1e28
            sky = np.asarray(h['SKY'].data[rows][:, pick], float) * 1e14
        flux[ivar <= 0] = np.nan
        err = 1.2533 * np.sqrt(np.nanmedian(np.where(ivar > 0, 1 / ivar, np.nan), 0) / np.sum(ivar > 0, 0))
        cube[method] = dict(flux=np.nanmedian(flux, 0), err=err, sky=np.nanmedian(sky, 0))
    methods = list(paths)
    fig = make_subplots(2, len(EXPOSURE_LINES), shared_xaxes='columns', vertical_spacing=0.08,
                        subplot_titles=list(EXPOSURE_LINES), row_heights=[0.6, 0.4])
    start = 0
    for c, (key, center) in enumerate(EXPOSURE_LINES.items(), start=1):
        sl = slice(start, start + cols[key].size)
        start += cols[key].size
        w = wave[cols[key]]
        for method in methods:
            label, color = METHOD_STYLE[method]
            fig.add_trace(go.Scatter(x=w, y=cube[method]['flux'][sl], name=label, legendgroup=method,
                                     showlegend=c == 1, line=dict(color=color, width=1.2, shape='hvh')), 1, c)
        base = cube[methods[0]]
        b = float(np.nanmedian(base['flux'][sl]))
        for sign in (1, -1):
            fig.add_trace(go.Scatter(x=w, y=b + sign * base['err'][sl], showlegend=False, hoverinfo='skip',
                                     line=dict(color='gray', dash='dash', width=0.8)), 1, c)
        for method in methods[1:]:
            label, color = METHOD_STYLE[method]
            diff = cube[method]['sky'][sl] - base['sky'][sl]
            fig.add_trace(go.Scatter(x=w, y=diff, name=f'SKY - SKY(deployed): {label}', showlegend=c == 1,
                                     line=dict(color=color, width=1.2, shape='hvh')), 2, c)
        fig.add_hline(y=0, line_color='lightgray', row=2, col=c)
        for r in (1, 2):
            fig.add_vline(x=center, line_dash='dot', line_color='lightgray', row=r, col=c)
        fig.update_xaxes(range=[center - half_width, center + half_width], row='all', col=c)
        fig.update_xaxes(title_text='Wavelength (A)', row=2, col=c)
    fig.update_yaxes(title_text='Median science FLUX (1e-14)', row=1, col=1)
    fig.update_yaxes(title_text='Sky difference (1e-14)', row=2, col=1)
    fig.update_layout(height=620, width=1500, legend=dict(orientation='h', y=-0.15))
    return fig


def exposure_leftover(paths, core=3.0, side=(6.0, 12.0)):
    """Median Science-fiber line flux left after subtraction, in % of the subtracted sky line (negative: too much)."""
    out = {}
    for method, path in paths.items():
        with fits.open(path, memmap=True) as h:
            wave = np.asarray(h['WAVE'].data, float)
            rows = _science_rows(h)
            out[METHOD_STYLE[method][0]] = {}
            for key, center in EXPOSURE_LINES.items():
                d = np.abs(wave - center)
                c, s = d <= core, (d > side[0]) & (d <= side[1])
                flux = np.nanmedian(np.asarray(h['FLUX'].data[rows][:, c], float), 0)
                base = np.nanmedian(np.nanmedian(np.asarray(h['FLUX'].data[rows][:, s], float), 0))
                sky = np.nanmedian(np.asarray(h['SKY'].data[rows][:, c], float), 0)
                out[METHOD_STYLE[method][0]][key] = 100 * np.sum(flux - base) / np.sum(sky)
    return pd.DataFrame(out).T
