"""Few-fiber test: per-fiber Sky Decomposition LSF versus the global block LSF, as 2D surfaces and 1D profiles (R and Z)."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from copy import deepcopy
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fiber_decomposition_test import reference_file, EXPOSURES, FIT, CHANNEL_INDEX  # noqa: E402

WAVE_FRACTIONS = (.15, .5, .85)
HALF_WIDTH = 2.5
OFFSETS = np.linspace(-HALF_WIDTH, HALF_WIDTH, 501)
METRICS = {'fwhm_angstrom': 'FWHM (Å)', 'w50_angstrom': 'W50 (Å)', 'w90_angstrom': 'W90 (Å)', 'centroid_angstrom': 'Centroid (Å)', 'asymmetry': 'Quantile asymmetry'}
CLEAN = {}
LABELS = {'individual': 'Per-fiber Sky Decomposition', 'block': 'Global block'}


def reference_profiles(expnum):
    """The nine reference fibers' own SolarBlue LSF surfaces, exactly those that produced the stored BESTFIT_LSF."""
    from astropy.io import fits
    path, rows = reference_file(expnum)
    profiles, spectrograph, clean = {}, {}, {}
    with fits.open(path, memmap=True) as h:
        meta = h['INPUT_META'].data
        config = {str(c['channel']).strip(): c for c in h['LSF_CONFIG'].data}
        for row in rows:
            f = int(meta['fiberid'][row])
            spectrograph[f], clean[f] = int(meta['spectrographid'][row]), bool(meta['gaia_clean_1pct'][row])
            profiles[f] = {ch: (np.asarray(h['LSF_COEF'].data[row, CHANNEL_INDEX[ch]], float),
                                np.asarray(h['LSF_KNOTS'].data[row, CHANNEL_INDEX[ch], :int(config[ch]['n_knots'])], float), int(config[ch]['degree'])) for ch in FIT}
    chosen = sorted((spectrograph[f], f) for f in profiles)
    return profiles, chosen, clean


def density(coefficient, knots, degree, wave):
    from skysub.sky_decomp.lsf_spline2d import mspline_basis
    from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
    basis = evaluate_bspline_basis(np.clip(wave, knots[degree], knots[-degree - 1]), knots, degree)
    return basis @ coefficient.T @ mspline_basis(OFFSETS, n_basis=coefficient.shape[0], half_width=HALF_WIDTH).T


def profile_set(row):
    from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics
    from skysub.sky_subtraction_block.reconstruction import SkyBlockReconstructor
    e = row['expnum']
    profiles, chosen, clean = reference_profiles(e)
    out = []
    with SkyBlockReconstructor(row['output'], cframe=row['input']) as r:
        template = r.state('MED_', r.labels.index('Sci'))
        for sp, f in chosen:
            state = r.global_state(f)
            for ch in FIT:
                lower, upper = template.knot_vectors[ch][template.degrees[ch]], template.knot_vectors[ch][-template.degrees[ch] - 1]
                grid = np.linspace(lower + .001, upper - .001, 121)
                samples = np.array([lower + q * (upper - lower) for q in WAVE_FRACTIONS])
                coefficient, knots, degree = profiles[f][ch]
                item = dict(expnum=e, spectrograph=sp, fiberid=f, gaia_clean=clean[f], channel=ch, grid=grid, samples=samples,
                            individual2d=density(coefficient, knots, degree, grid), block2d=density(state.coefficients[ch], state.knot_vectors[ch], state.degrees[ch], grid),
                            individual1d=density(coefficient, knots, degree, samples), block1d=density(state.coefficients[ch], state.knot_vectors[ch], state.degrees[ch], samples))
                item['metrics_block'] = evaluate_lsf_diagnostics(state, samples)
                item['curve_block'] = evaluate_lsf_diagnostics(state, grid)
                alt = deepcopy(state)
                alt.coefficients[ch], alt.knot_vectors[ch], alt.degrees[ch] = coefficient, knots, degree
                item['metrics_individual'] = evaluate_lsf_diagnostics(alt, samples)
                item['curve_individual'] = evaluate_lsf_diagnostics(alt, grid)
                out.append(item)
    return out


def figure_metrics(items, e, ch):
    items = [i for i in items if i['channel'] == ch]
    colors = ('#1b9e77', '#7570b3', '#d95f02')
    fig, axes = plt.subplots(2 * len(METRICS), 3, figsize=(15, 3.1 * len(METRICS)), layout='constrained',
                             gridspec_kw=dict(height_ratios=[3, 1.2] * len(METRICS)), sharex='col')
    for col, sp in enumerate((1, 2, 3)):
        group = [i for i in items if i['spectrograph'] == sp]
        for row, (key, label) in enumerate(METRICS.items()):
            level, diff = axes[2 * row, col], axes[2 * row + 1, col]
            for n, item in enumerate(group):
                name = ('low', 'middle', 'high')[n] + ' fiber'
                level.plot(item['grid'], item['curve_individual'][key], color=colors[n], lw=1.3, label=f'{name}: per-fiber')
                level.plot(item['grid'], item['curve_block'][key], color=colors[n], lw=1.3, ls='--', label=f'{name}: block')
                diff.plot(item['grid'], item['curve_block'][key] - item['curve_individual'][key], color=colors[n], lw=1)
            diff.axhline(0, color='#666666', lw=.6)
            level.grid(alpha=.2)
            diff.grid(alpha=.2)
            if row == 0:
                level.set_title(f'SP{sp} · fibers ' + ' / '.join(str(i['fiberid']) for i in group))
            if col == 0:
                level.set_ylabel(label)
                diff.set_ylabel('Block − fiber')
            if row == len(METRICS) - 1:
                diff.set_xlabel('Air wavelength (Å)')
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='outside lower center', ncol=3, fontsize=10)
    fig.suptitle(f'Exposure {e}, {ch}: LSF metrics versus wavelength, per-fiber Sky Decomposition (solid) and global block (dashed); three fibers per spectrograph', fontsize=14)
    name = f'fiber_lsf_metrics_{ch}_{e}.png'
    fig.savefig(HERE / name, dpi=80)
    plt.close(fig)
    return name


def curve_rows(items):
    rows = []
    for item in items:
        for n, wave in enumerate(item['grid']):
            record = dict(expnum=item['expnum'], spectrograph=item['spectrograph'], fiberid=item['fiberid'], channel=item['channel'], wavelength_A=float(wave))
            for key in METRICS:
                record[key + '_individual'] = float(item['curve_individual'][key][n])
                record[key + '_block'] = float(item['curve_block'][key][n])
            rows.append(record)
    return rows


def figure_2d(items, e, ch):
    items = [i for i in items if i['channel'] == ch]
    fig, axes = plt.subplots(len(items), 3, figsize=(15, 2.5 * len(items)), layout='constrained', sharex='col')
    top = float(max(np.percentile(i['individual2d'], 99.5) for i in items))
    for row, item in enumerate(items):
        d = item['block2d'] - item['individual2d']
        spread = float(np.percentile(np.abs(d), 99.5))
        for col, (values, title, options) in enumerate((
                (item['individual2d'], 'Per-fiber Sky Decomposition', dict(cmap='viridis', vmin=0, vmax=top)),
                (item['block2d'], 'Global block', dict(cmap='viridis', vmin=0, vmax=top)),
                (d, 'Block − per-fiber', dict(cmap='RdBu_r', vmin=-spread, vmax=spread)))):
            ax = axes[row, col]
            mesh = ax.pcolormesh(item['grid'], OFFSETS, values.T, shading='nearest', rasterized=True, **options)
            if row == 0:
                ax.set_title(title)
            ax.axhline(0, color='white' if col < 2 else 'black', lw=.5, alpha=.6)
            if col == 0:
                ax.set_ylabel(f'SP{item["spectrograph"]} · fiber {item["fiberid"]}\nOffset (Å)')
            if row == len(items) - 1:
                ax.set_xlabel('Air wavelength (Å)')
            fig.colorbar(mesh, ax=ax, label='LSF density (Å⁻¹)' if col < 2 else 'Δ density (Å⁻¹)')
    fig.suptitle(f'Exposure {e}, {ch}: LSF as a function of wavelength and offset (shared density scale in the first two columns; symmetric difference scale per row)', fontsize=14)
    name = f'fiber_lsf_2d_{ch}_{e}.png'
    fig.savefig(HERE / name, dpi=80)
    plt.close(fig)
    return name


def figure_1d(items, e, ch):
    items = [i for i in items if i['channel'] == ch]
    rows = len(items)
    fig, axes = plt.subplots(2 * rows, 3, figsize=(15, 3.2 * rows), layout='constrained',
                             gridspec_kw=dict(height_ratios=[3, 1.1] * rows), sharex=True)
    for r_, item in enumerate(items):
        for c_, wave in enumerate(item['samples']):
            a, b = item['individual1d'][c_], item['block1d'][c_]
            ax = axes[2 * r_, c_]
            ax.plot(OFFSETS, a, color='#1f77b4', lw=1.4, label=LABELS['individual'])
            ax.plot(OFFSETS, b, color='#e07b00', lw=1.4, label=LABELS['block'])
            ax.set_title(f'SP{item["spectrograph"]} · fiber {item["fiberid"]} · {wave:.0f} Å', fontsize=10)
            ax.grid(alpha=.2)
            if c_ == 0:
                ax.set_ylabel('LSF density (Å⁻¹)')
            diff = axes[2 * r_ + 1, c_]
            diff.plot(OFFSETS, b - a, color='#a34242', lw=1)
            diff.axhline(0, color='#666666', lw=.6)
            diff.grid(alpha=.2)
            if c_ == 0:
                diff.set_ylabel('Block − fiber')
            if r_ == rows - 1:
                diff.set_xlabel('Offset from input wavelength (Å)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=2, fontsize=11)
    fig.suptitle(f'Exposure {e}, {ch}: one-dimensional LSF profiles and signed differences', fontsize=14)
    name = f'fiber_lsf_1d_{ch}_{e}.png'
    fig.savefig(HERE / name, dpi=80)
    plt.close(fig)
    return name


def table_rows(items):
    rows = []
    for item in items:
        for c_, wave in enumerate(item['samples']):
            a, b = item['individual1d'][c_], item['block1d'][c_]
            record = dict(expnum=item['expnum'], spectrograph=item['spectrograph'], fiberid=item['fiberid'], channel=item['channel'], wavelength_A=float(wave),
                          l1_distance=float(np.trapz(np.abs(b - a), OFFSETS)), peak_individual=float(a.max()), peak_block=float(b.max()))
            for key in ('fwhm_angstrom', 'w50_angstrom', 'w90_angstrom', 'centroid_angstrom'):
                record[key + '_individual'] = float(item['metrics_individual'][key][c_])
                record[key + '_block'] = float(item['metrics_block'][key][c_])
            rows.append(record)
    return rows


def section(frame, names, selection):
    body = '''<div id="fiber-profiles"><h2>R and Z: per-fiber Sky Decomposition versus the global block LSF</h2>
<p>The same nine reference Science fibers as in the spectral test above (three per spectrograph; the fibers listed below) in each of the five exposures. <b>Per-fiber Sky Decomposition</b>: the LSF of the SolarBlue decomposition of that fiber's own spectrum in that exposure, taken from the reference model that produced the displayed fit. <b>Global block</b>: the sky-only surface of the spectrograph evaluated at that fiber; these Science fibers were not used in the block fit. W50 = q75 − q25 and W90 = q95 − q05 of the profile; FWHM is the literal full width at half maximum; asymmetry = (q95 + q05 − 2·q50)/W90. B is not fitted. LSF density is in Å⁻¹ with unit integral over the ±2.5 Å offset support; both surfaces keep their own centroid (the block centroid is zero by construction). The Gaia column of the selection table flags fibers whose Gaia-to-observed flux ratio exceeds 1%: for those, stellar flux contaminates the Science spectrum and the per-fiber LSF is less reliable.</p>'''
    body += '<h3>Selected fibers</h3><table><tr><th>Exposure</th><th>Spectrograph</th><th>Physical fibers (low / middle / high)</th><th>Gaia-clean (&lt;1%) flags</th></tr>'
    for (e, sp), ids in selection.items():
        body += f'<tr><td>{e}</td><td>SP{sp}</td><td>{" / ".join(str(i) for i in ids)}</td><td>{" / ".join("clean" if CLEAN[(e, i)] else "Gaia-contaminated" for i in ids)}</td></tr>'
    body += '</table>'
    body += '<h3>Profile differences at the plotted wavelengths</h3><p>L1 distance = ∫|block − fiber| dΔλ (0: identical; 2: no overlap). FWHM and W50/W90 in Å; the table lists medians over the nine fibers for each exposure, channel and wavelength position (low/middle/high of the fitted channel range).</p>'
    body += '<table><tr><th>Exposure</th><th>Channel</th><th>Wavelength position</th><th>Median L1</th><th>Median ΔFWHM (block − fiber, Å)</th><th>Median ΔW50 (Å)</th><th>Median ΔW90 (Å)</th><th>Median Δcentroid (Å)</th></tr>'
    for (e, ch), g in frame.groupby(['expnum', 'channel']):
        g = g.copy()
        g['position'] = g.groupby(['spectrograph', 'fiberid']).cumcount()
        for position, label in enumerate(('low', 'middle', 'high')):
            s = g[g.position == position]
            body += f'<tr><td>{e}</td><td>{ch}</td><td>{label} ({s.wavelength_A.median():.0f} Å)</td><td>{s.l1_distance.median():.3f}</td>' + ''.join(
                f'<td>{(s[k + "_block"] - s[k + "_individual"]).median():+.3f}</td>' for k in ('fwhm_angstrom', 'w50_angstrom', 'w90_angstrom', 'centroid_angstrom')) + '</tr>'
    body += '</table>'
    for e in sorted(frame.expnum.unique()):
        body += f'<h3>Exposure {e}</h3>'
        for ch in FIT:
            body += (f'<h4>{ch}: metrics versus wavelength</h4><img src="{names[(e, ch, "metrics")]}" alt="LSF metrics versus wavelength of nine fibers, exposure {e} channel {ch}">'
                     f'<h4>{ch}: two-dimensional surfaces</h4><img src="{names[(e, ch, "2d")]}" alt="Two-dimensional LSF of nine fibers, exposure {e} channel {ch}">'
                     f'<h4>{ch}: profile shape at three wavelengths</h4><img src="{names[(e, ch, "1d")]}" alt="One-dimensional LSF profiles of nine fibers, exposure {e} channel {ch}">')
    body += '<p><a href="fiber_lsf_profile_metrics.csv">Values at the three plotted wavelengths</a> · <a href="fiber_lsf_metric_curves.csv">Metric curves (121 wavelengths)</a></p></div>'
    return body


if __name__ == '__main__':
    manifest = {r['expnum']: r for r in json.loads((HERE / 'run_manifest.json').read_text())}
    records, curves, names, selection = [], [], {}, {}
    for e in EXPOSURES:
        print('Exposure', e, flush=True)
        items = profile_set(manifest[e])
        for ch in FIT:
            names[(e, ch, '2d')] = figure_2d(items, e, ch)
            names[(e, ch, '1d')] = figure_1d(items, e, ch)
            names[(e, ch, 'metrics')] = figure_metrics(items, e, ch)
        records += table_rows(items)
        curves += curve_rows(items)
        for sp in (1, 2, 3):
            selection[(e, sp)] = sorted({i['fiberid'] for i in items if i['spectrograph'] == sp})
        CLEAN.update({(e, i['fiberid']): i['gaia_clean'] for i in items})
    frame = pd.DataFrame(records)
    frame.to_csv(HERE / 'fiber_lsf_profile_metrics.csv', index=False)
    pd.DataFrame(curves).to_csv(HERE / 'fiber_lsf_metric_curves.csv', index=False)
    (HERE / 'fiber_profiles_section.html').write_text(section(frame, names, selection))
    print(frame.groupby(['expnum', 'channel']).l1_distance.median().round(3).to_string())
