"""Report for the nine-fiber test: spectra, residuals against input errors, and the in-fiber versus block comparison (R and Z)."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import html
import json
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / 'analysis/sky_subtraction_validation'))
from report_fibers import trace, finish, page, UNIT  # noqa: E402
from fiber_decomposition_test import EXPOSURES, VARIANTS, GROUPS, WORK, SOLARBLUE, EXTRA, selected_record, reference_file  # noqa: E402

LABEL = {'fiber': 'Per-fiber fit (own LSF, stored)', 'block_same': 'Block LSF, same coefficients', 'fiber_refit': 'Per-fiber LSF, coefficients refit',
         'block_refit': 'Block LSF, coefficients refit'}
SHORT = {'fiber': 'Per-fiber (stored)', 'block_same': 'Block, same coef.', 'fiber_refit': 'Per-fiber LSF, refit', 'block_refit': 'Block, refit'}
COLOR = {'fiber': '#1f77b4', 'block_same': '#e0a000', 'fiber_refit': '#6baed6', 'block_refit': '#d62728'}
FAMILY_LABEL = {'oh': 'OH airglow', 'moon': 'Scattered moonlight', 'zodi': 'Zodiacal light', 'diffuse': 'Diffuse airglow continuum', 'atom': 'Atomic airglow',
                'orc': 'Oxygen recombination', 'o2': 'O₂ airglow'}
FLUX_NOTE = 'Units: 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹.'


def clean_flags():
    flags = {}
    for e in EXPOSURES:
        path, rows = reference_file(e)
        with fits.open(path, memmap=True) as h:
            meta = h['INPUT_META'].data
            for row in rows:
                flags[(e, int(meta['fiberid'][row]))] = (bool(meta['gaia_clean_1pct'][row]), float(meta['gaia_ratio'][row]))
    return flags


def spectrum_page(e, k, data, flags, record):
    wave = data['wave']
    f, sp = int(data['fiberid'][k]), int(data['spectrograph'][k])
    fit = data['fit_pixels'][k]
    flux, sigma = data['flux'][k].astype(float), data['sigma'][k].astype(float)
    model = {v: data[f'model_{v}'][k].astype(float) for v in VARIANTS}
    good = fit & np.isfinite(sigma)
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[.52, .24, .24], vertical_spacing=.07)
    for key, label in FAMILY_LABEL.items():
        if key in ('moon', 'zodi', 'diffuse', 'oh', 'atom', 'orc', 'o2'):
            trace(fig, wave, np.where(good, data[f'comp_fiber_{key}'][k], np.nan), label + ' (per-fiber fit)', 1, width=.6)
    trace(fig, wave, np.where(good, flux, np.nan), 'Observed', 1, color='#333333', width=1)
    trace(fig, wave, np.where(good, model['fiber'], np.nan), LABEL['fiber'], 1, color=COLOR['fiber'], width=1.3)
    trace(fig, wave, np.where(good, model['block_refit'], np.nan), LABEL['block_refit'], 1, color=COLOR['block_refit'], width=1.3)
    trace(fig, wave, np.where(good, model['block_same'], np.nan), LABEL['block_same'], 1, color=COLOR['block_same'], width=1.1, dash='dot')
    fig.data[-1].visible = 'legendonly'
    rms = float(np.sqrt(np.mean((flux - model['fiber'])[good] ** 2)))
    for variant, offset in (('fiber', -5 * rms), ('block_refit', -10 * rms)):
        for sign in (-1, 1):
            trace(fig, wave, np.where(good, offset + sign * sigma, np.nan), 'Offset ±1σ', 1, color='black', showlegend=(variant == 'fiber' and sign == 1), width=.6)
        trace(fig, wave, np.where(good, offset + flux - model[variant], np.nan),
              ('Data − per-fiber model (−5 RMS)' if variant == 'fiber' else 'Data − block-refit model (−10 RMS)'), 1, color=COLOR[variant], width=.8)
    for variant in ('fiber', 'block_refit'):
        trace(fig, wave, np.where(good, (flux - model[variant]) / sigma, np.nan), f'(Data − model)/σ: {SHORT[variant]}', 2, color=COLOR[variant], width=.8)
    for level in (-1, 1):
        trace(fig, wave, np.where(good, level, np.nan), '−1σ' if level < 0 else '+1σ', 2, color='#9ebbd2', width=.6)
    fig.data[-1].fill = 'tonexty'
    fig.data[-1].fillcolor = 'rgba(135,175,210,.17)'
    trace(fig, wave, np.where(good, 0, np.nan), 'Zero', 2, color='#777777', width=.6, showlegend=False)
    for variant in ('block_same', 'block_refit'):
        trace(fig, wave, np.where(good, (model[variant] - model['fiber']) / sigma, np.nan), f'(Model − per-fiber model)/σ: {SHORT[variant]}', 3, color=COLOR[variant], width=.8)
    trace(fig, wave, np.where(good, 0, np.nan), 'Zero', 3, color='#777777', width=.6, showlegend=False)
    for level in (-1, 1):
        fig.add_annotation(x=1, xref='paper', y=level, yref='y2', text='−1σ' if level < 0 else '+1σ', showarrow=False, xanchor='left', xshift=3, yshift=level * 8, font=dict(size=10))
    for xb, name in ((5787, 'B | R'), (7454, 'R | Z')):
        for row in (1, 2, 3):
            fig.add_vline(x=xb, line=dict(color='#bbbbbb', width=.8, dash='dash'), row=row, col=1)
    fig.update_xaxes(title_text='Air wavelength (Å)', row=3, col=1, range=[float(wave[0]), float(wave[-1])])
    fig.update_yaxes(title_text=UNIT, row=1, col=1)
    fig.update_yaxes(title_text='Residual / σ', row=2, col=1)
    fig.update_yaxes(title_text='Model difference / σ', row=3, col=1)
    fig.add_annotation(x=0, y=1.035, xref='paper', yref='paper', text='Observed spectrum, components, models and flux residuals (offsets stated in the legend)', showarrow=False, xanchor='left', font=dict(size=12))
    fig.add_annotation(x=0, y=.455, xref='paper', yref='paper', text='Residuals in input-noise units; shaded ±1σ', showarrow=False, xanchor='left', font=dict(size=12))
    fig.add_annotation(x=0, y=.19, xref='paper', yref='paper', text='Difference of each block model from the per-fiber model, in input-noise units', showarrow=False, xanchor='left', font=dict(size=12))
    clean, ratio = flags[(e, f)]
    title = (f'Exposure {e} · MJD {record["mjd"]:.0f} · Science SP{sp} · fiber {f} · Gaia/observed flux {ratio:.3f}'
             + ('' if clean else ' (above 1%: stellar contamination)'))
    body = (f'<h1>{html.escape(title)}</h1><p>σ = 1/√IVAR: input observed noise (including the exact fit mask: input mask and the Hα-centred science-line mask). '
            'The per-fiber model is the stored SolarBlue fit of this spectrum with its own fitted R/Z LSF. The block models use the global block R/Z LSF evaluated at this fiber; B keeps the fiber’s own LSF in all models. '
            '“Coefficients refit” re-solves the production decomposition QP with the stated LSF held fixed. '
            'RMS is unweighted over fit pixels. Everything is shown over the native range; drag to zoom. ' + FLUX_NOTE + '</p>')
    name = f'fiber_decomp_spectrum_{e}_{f}.html'
    (HERE / name).write_text(page(title, body + finish(fig, title, 1050) + '<p><a href="00_MAIN_sky_block_report.html">Main report</a></p>'))
    return name


def summary_figures(metrics, categories):
    plt.rcParams.update({'font.size': 10})
    fibers = metrics.drop_duplicates(['expnum', 'fiberid']).sort_values(['expnum', 'spectrograph', 'fiberid'])
    fig, axes = plt.subplots(2, len(EXPOSURES), figsize=(4.2 * len(EXPOSURES), 7.5), layout='constrained', sharey='row')
    for j, e in enumerate(EXPOSURES):
        order = list(fibers[fibers.expnum == e].fiberid)
        for i, ch in enumerate(('R', 'Z')):
            ax = axes[i, j]
            for variant in VARIANTS:
                sel = metrics[(metrics.expnum == e) & (metrics.channel == ch) & (metrics.variant == variant)].set_index('fiberid').loc[order]
                ax.plot(range(len(order)), sel.chi2_per_pixel, marker='o', ms=4.5, lw=0, color=COLOR[variant], label=SHORT[variant], alpha=.85)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels([str(x) for x in order], rotation=90, fontsize=8)
            ax.set_yscale('log')
            ax.grid(alpha=.25)
            for edge in (2.5, 5.5):
                ax.axvline(edge, color='#999999', lw=.7)
            ax.set_title(f'{e} · {ch}')
            if j == 0:
                ax.set_ylabel('χ² per pixel (input σ)')
            if i == 1:
                ax.set_xlabel('Physical fiber ID (SP1 | SP2 | SP3)')
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='outside lower center', ncol=4, fontsize=10)
    fig.suptitle('Nine reference Science fibers per exposure: χ² per pixel on identical fit pixels', fontsize=14)
    fig.savefig(HERE / 'fiber_decomp_chi2.png', dpi=90)
    plt.close(fig)
    names = list(GROUPS)
    fig, axes = plt.subplots(2, len(EXPOSURES), figsize=(4.2 * len(EXPOSURES), 7.5), layout='constrained', sharey='row')
    for j, e in enumerate(EXPOSURES):
        for i, ch in enumerate(('R', 'Z')):
            ax = axes[i, j]
            for v, variant in enumerate(VARIANTS):
                sel = categories[(categories.expnum == e) & (categories.channel == ch) & (categories.variant == variant)]
                grouped = sel.groupby('category').apply(lambda g: g.chi2_sum.sum() / g.n_pixels.sum(), include_groups=False).reindex(names)
                ax.bar(np.arange(len(names)) + (v - 1.5) * .2, grouped.values, width=.2, color=COLOR[variant], label=SHORT[variant])
            ax.set_yscale('log')
            ax.set_xticks(range(len(names)))
            ax.set_xticklabels(names, rotation=75, ha='right', fontsize=8)
            ax.grid(alpha=.25, axis='y')
            ax.set_title(f'{e} · {ch}')
            if j == 0:
                ax.set_ylabel('χ² per pixel in category')
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='outside lower center', ncol=4, fontsize=10)
    fig.suptitle('χ² per pixel by pixel category (pooled over the nine fibers; category from distance to the nearest strong model line)', fontsize=14)
    fig.savefig(HERE / 'fiber_decomp_categories.png', dpi=90)
    plt.close(fig)


def block_facts():
    manifest = {r['expnum']: r for r in json.loads((HERE / 'run_manifest.json').read_text())}
    rows = []
    for e in EXPOSURES:
        with fits.open(manifest[e]['output'], memmap=True) as h:
            for b in h['BLOCK_LSF_META'].data:
                rows.append(dict(expnum=e, spectrograph=int(b['spectrograph']), channel=str(b['channel']).strip(), n_fibers=int(b['n_fibers']), n_pixels=int(b['n_pixels']),
                                 shape=tuple(h['BLOCK_LSF_COEF'].data[0].shape)))
    return pd.DataFrame(rows)


def section(metrics, categories, pages, flags, checks):
    block = block_facts()
    shape = block['shape'].iloc[0]
    sky_fibers = f'{block.n_fibers.min()}–{block.n_fibers.max()}'
    sky_pixels = f'{block.n_pixels.min():,}–{block.n_pixels.max():,}'
    body = '<div id="fiber-decomposition"><h2>Nine reference Science fibers: in-fiber Sky Decomposition versus the block LSF (R and Z)</h2>'
    body += ('<p>The nine reference fibers of the SolarBlue runs (three per spectrograph: 65, 329, 584 | 714, 976, 1234 | 1362, 1618, 1881), in each of 4108, 16819, 40620, 44808 and 60919. '
             'All four variants fit the same spectrum, with the same fit pixels, input σ, science-line mask, O₂ transmission and decomposition model; only the R and Z LSF and the coefficients differ. B keeps the fiber’s own LSF everywhere. '
             '<b>Per-fiber (stored)</b>: the stored SolarBlue decomposition of that spectrum, re-rendered from its recorded coefficients and LSF. '
             '<b>Block, same coef.</b>: the block R/Z LSF evaluated at the fiber, with those stored coefficients unchanged. '
             '<b>Per-fiber LSF, refit</b> and <b>Block, refit</b>: the production decomposition QP re-solved with the fiber’s own, or the block, R/Z LSF held fixed. '
             'The refit pair differs only in the LSF; the stored fit additionally reflects the fiber’s iterative LSF refinement. χ² per pixel is the mean squared (data − model)/σ over R or Z fit pixels, with 5 Å trimmed at the channel edges; it is not a reduced χ².</p>')
    worst = checks.max_abs_render_minus_stored_bestfit.max()
    body += (f'<p>Check: re-rendering the stored coefficients with the stored LSF reproduces the file’s BESTFIT_LSF for all 45 spectra (largest difference {worst:.1e} in units of 10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹, '
             f'against a median fitted flux of {checks.median_abs_stored_bestfit.median():.1f}).</p>')
    body += '<h3>Are the two fits comparable in dimension?</h3><table><tr><th>Item</th><th>In-fiber fit</th><th>Block fit</th></tr>'
    body += (f'<tr><td>LSF unknowns per channel</td><td>11 offset × 4 wavelength = 44 per fiber</td><td>11 × 4 × 4 slit = {np.prod(shape[1:]) if False else 11*4*4} per spectrograph (shared by all its fibers; {len(shape)}-index array {shape} stored per block)</td></tr>'
             f'<tr><td>Data constraining them</td><td>one spectrum (~1.7–2.3 thousand R or Z pixels)</td><td>{sky_fibers} Sky fibers, {sky_pixels} pixels per block fit</td></tr>'
             '<tr><td>Where the LSF is estimated</td><td>from the same Science spectrum it is tested on (in-sample)</td><td>from Sky fibers only; these Science fibers are out-of-sample</td></tr>'
             '<tr><td>Centroid</td><td>free: absorbs wavelength-solution error</td><td>zero by construction at every wavelength and slit position</td></tr>'
             '<tr><td>Other constraints</td><td>nonnegative, unit integral, coefficient monotonicity; offset second-difference penalty α = 0.01</td><td>nonnegative, unit integral, zero first moment; no penalty</td></tr>'
             '<tr><td>Component coefficients</td><td>refit together with the LSF over up to 30 cycles</td><td>fixed per Sky fiber from a fit at the telescope-median LSF; not iterated with the LSF</td></tr>'
             '<tr><td>Science fiber at test time</td><td>its own coefficients</td><td>ML-predicted coefficients in the product; this test also shows stored and refit coefficients</td></tr></table>')
    body += '<p>The block has four times as many unknowns per channel as a single fiber, but is constrained by two orders of magnitude more data; it carries no free centroid, so wavelength-solution errors cannot be absorbed by the LSF.</p>'
    cl = pd.DataFrame([dict(expnum=e, fiberid=f, clean=c, ratio=r) for (e, f), (c, r) in flags.items()])
    body += '<h3>Pooled result</h3><table><tr><th>Exposure</th><th>Channel</th><th>Gaia-clean fibers</th>' + ''.join(f'<th>χ²/pixel: {SHORT[v]}</th>' for v in VARIANTS) + '<th>RMS(block refit − per-fiber)/σ</th><th>Pixels within ±1σ: per-fiber / block refit</th></tr>'
    for e in EXPOSURES:
        n_clean = int(cl[cl.expnum == e].clean.sum())
        for ch in ('R', 'Z'):
            sel = metrics[(metrics.expnum == e) & (metrics.channel == ch)]
            med = sel.groupby('variant').chi2_per_pixel.median()
            within = sel.groupby('variant').within_1sigma.median()
            diff = sel[sel.variant == 'block_refit'].rms_model_difference_sigma.median()
            body += f'<tr><td>{e}</td><td>{ch}</td><td>{n_clean} / 9</td>' + ''.join(f'<td>{med[v]:.2f}</td>' for v in VARIANTS) + f'<td>{diff:.2f}</td><td>{within["fiber"]:.0%} / {within["block_refit"]:.0%}</td></tr>'
    body += '</table><p>Medians over the nine fibers. Stellar contamination (Gaia flux above 1% of the observed flux) is flagged in the fiber table; 16819 and 44808 have few clean fibers, so conclusions there are limited.</p>'
    body += '<h3>χ² per pixel for every fiber</h3><img src="fiber_decomp_chi2.png" alt="Chi-square per pixel of four model variants for nine fibers in each exposure, R and Z"><p>Log scale; vertical lines separate spectrographs; each point is one fiber and variant on identical pixels.</p>'
    body += '<h3>Where the residual sits</h3><img src="fiber_decomp_categories.png" alt="Chi-square per pixel by pixel category for four model variants"><p>Categories are exclusive: OH and other-line cores (&lt;1 Å from a line carrying at least 1% of the strongest line flux), wings (1–2.5 Å), the zone beyond the LSF support (2.5–4 Å), and continuum/solar structure (&gt;4 Å). Pooled over the nine fibers.</p>'
    body += '<h3>Interactive spectra</h3><p>Each page covers B, R and Z: components, observed spectrum, models, flux residuals at −5 RMS (per-fiber) and −10 RMS (block refit) with thin ±σ curves, normalized residuals, and each block model’s difference from the per-fiber model in σ.</p>'
    body += '<table><tr><th>Exposure</th><th>SP</th><th>Fiber</th><th>Gaia flux / observed</th>' + ''.join(f'<th>R / Z χ²/px: {SHORT[v]}</th>' for v in ('fiber', 'block_same', 'block_refit')) + '<th>Spectrum</th></tr>'
    for e in EXPOSURES:
        sel = metrics[metrics.expnum == e]
        for f in sorted(sel.fiberid.unique(), key=lambda x: (int(sel[sel.fiberid == x].spectrograph.iloc[0]), x)):
            one = sel[sel.fiberid == f]
            clean, ratio = flags[(e, f)]
            cells = ''.join('<td>' + ' / '.join(f'{one[(one.channel == ch) & (one.variant == v)].chi2_per_pixel.iloc[0]:.2f}' for ch in ('R', 'Z')) + '</td>' for v in ('fiber', 'block_same', 'block_refit'))
            body += f'<tr><td>{e}</td><td>SP{int(one.spectrograph.iloc[0])}</td><td>{f}</td><td>{ratio:.3f}{"" if clean else " (contaminated)"}</td>{cells}<td><a href="{pages[(e, f)]}">open</a></td></tr>'
    body += '</table><p><a href="fiber_decomp_metrics_all.csv">All metrics</a> · <a href="fiber_decomp_categories_all.csv">Category budget</a> · <a href="fiber_decomp_render_check_all.csv">Render check</a></p></div>'
    return body


def main():
    flags = clean_flags()
    metrics = pd.concat([pd.read_csv(HERE / f'fiber_decomp_metrics_{e}.csv') for e in EXPOSURES], ignore_index=True)
    categories = pd.concat([pd.read_csv(HERE / f'fiber_decomp_categories_{e}.csv') for e in EXPOSURES], ignore_index=True)
    checks = pd.concat([pd.read_csv(HERE / f'fiber_decomp_render_check_{e}.csv') for e in EXPOSURES], ignore_index=True)
    metrics.to_csv(HERE / 'fiber_decomp_metrics_all.csv', index=False)
    categories.to_csv(HERE / 'fiber_decomp_categories_all.csv', index=False)
    checks.to_csv(HERE / 'fiber_decomp_render_check_all.csv', index=False)
    pages = {}
    for e in EXPOSURES:
        data = np.load(WORK / f'fiber_decomposition_{e}.npz')
        record = selected_record(e)
        for k in range(len(data['fiberid'])):
            pages[(e, int(data['fiberid'][k]))] = spectrum_page(e, k, data, flags, record)
        print('Spectra', e, flush=True)
    summary_figures(metrics, categories)
    (HERE / 'fiber_decomposition_section.html').write_text(section(metrics, categories, pages, flags, checks))


if __name__ == '__main__':
    main()
