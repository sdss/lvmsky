"""Figures and HTML section for the Sky-fiber residual budget."""
from pathlib import Path
import html
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from residual_budget import CATEGORIES, WORK, TELLURIC  # noqa: E402

GROUPS = [('Sky-line cores and wings (<2.5 Å)', lambda c: 'core' in c or 'wing' in c),
          ('Outside LSF support (2.5–4 Å)', lambda c: 'Outside' in c),
          ('Solar structure', lambda c: 'Solar' in c),
          ('Smooth continuum', lambda c: 'Smooth' in c)]
COLORS = ['#c0392b', '#e69f00', '#3b6fb6', '#4a9d6b']
CONTROL = 3469
FOCUS = [16819, 60919, CONTROL]


def group_table(summary):
    rows = []
    for (e, ch), frame in summary.groupby(['expnum', 'channel']):
        total_chi, total_n = frame.chi2_sum.sum(), frame.n_pixels.sum()
        for name, test in GROUPS:
            part = frame[frame.category.map(test)]
            rows.append(dict(expnum=e, channel=ch, group=name, n_pixels=int(part.n_pixels.sum()),
                             chi2_per_pixel=part.chi2_sum.sum() / part.n_pixels.sum() if part.n_pixels.sum() else np.nan,
                             chi2_share=part.chi2_sum.sum() / total_chi, pixel_share=part.n_pixels.sum() / total_n))
    return pd.DataFrame(rows)


def budget_figure(table, exposures):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True, layout='constrained')
    width = .2
    for ax, ch in zip(axes, 'BRZ'):
        for k, (name, _) in enumerate(GROUPS):
            values = [table[(table.expnum == e) & (table.channel == ch) & (table.group == name)].chi2_per_pixel.iloc[0] for e in exposures]
            ax.bar(np.arange(len(exposures)) + (k - 1.5) * width, values, width, color=COLORS[k], label=name)
        ax.axhline(1, color='black', lw=.8)
        ax.set_yscale('log')
        ax.set_xticks(np.arange(len(exposures)), [str(e) for e in exposures], rotation=45)
        ax.set_title(f'{ch} channel')
        ax.set_xlabel('Exposure')
        ax.grid(axis='y', alpha=.25)
    axes[0].set_ylabel('Mean (data − model)² / σ² per pixel')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4)
    fig.suptitle('Sky-fiber residual by pixel category (all selected Sky fibers pooled; black line = 1)', fontsize=15)
    fig.savefig(HERE / 'residual_budget_categories.png', dpi=120)
    plt.close(fig)


def telluric_figure(telluric, exposures):
    windows = list(TELLURIC)
    values = np.array([[telluric[(telluric.expnum == e) & (telluric.window == w)].chi2_per_pixel.median() for w in windows] for e in exposures])
    fig, ax = plt.subplots(figsize=(11, 5), layout='constrained')
    image = ax.imshow(np.log10(values), cmap='viridis', aspect='auto', vmin=0, vmax=np.log10(70))
    for i in range(len(exposures)):
        for j in range(len(windows)):
            ax.text(j, i, f'{values[i, j]:.1f}', ha='center', va='center', color='white' if np.log10(values[i, j]) < 1 else 'black', fontsize=11)
    ax.set_xticks(range(len(windows)), windows, rotation=30, ha='right')
    ax.set_yticks(range(len(exposures)), [str(e) for e in exposures])
    ax.set_ylabel('Exposure')
    cb = fig.colorbar(image, ax=ax)
    cb.set_label('log₁₀ median χ² per pixel across Sky fibers')
    ax.set_title('Residual in telluric windows (numbers: median χ² per pixel)')
    fig.savefig(HERE / 'residual_budget_telluric.png', dpi=120)
    plt.close(fig)


def lowpass_map(exposures):
    fig, axes = plt.subplots(len(exposures), 1, figsize=(15, 4.3 * len(exposures)), layout='constrained')
    for ax, e in zip(np.atleast_1d(axes), exposures):
        data = np.load(WORK / f'residual_arrays_{e}.npz', allow_pickle=True)
        wave, z = data['wave'], data['z_low']
        image = ax.imshow(z, aspect='auto', origin='lower', extent=[wave[0], wave[-1], -.5, len(z) - .5],
                          vmin=-3, vmax=3, cmap='RdBu_r', interpolation='nearest')
        order = data['spectrograph']
        for k in np.flatnonzero(np.diff(order) != 0):
            ax.axhline(k + .5, color='black', lw=.6)
        for lo, hi in ((5787, 5787), (7454, 7454)):
            ax.axvline(lo, color='black', lw=.5, ls=':')
        ax.set_title(f'Exposure {e}')
        ax.set_ylabel('Selected Sky fiber (stored order)')
        ax.set_xlabel('Air wavelength (Å)')
        cb = fig.colorbar(image, ax=ax, shrink=.9)
        cb.set_label('Running-median (data − model) / σ')
    fig.suptitle('Smooth residual (≈20 Å running median) for every selected Sky fiber; shared scale ±3σ, saturated beyond', fontsize=15)
    fig.savefig(HERE / 'residual_budget_smooth_maps.png', dpi=110)
    plt.close(fig)


def table_html(frame, columns, formats):
    body = '<table><tr>' + ''.join(f'<th>{html.escape(c)}</th>' for c in columns) + '</tr>'
    for _, row in frame.iterrows():
        body += '<tr>' + ''.join(f'<td>{formats[c](row[c])}</td>' for c in columns) + '</tr>'
    return body + '</table>'


def main():
    summary = pd.read_csv(HERE / 'residual_budget_summary.csv')
    telluric = pd.read_csv(HERE / 'residual_telluric_metrics.csv')
    scale = pd.read_csv(HERE / 'residual_scale_metrics.csv')
    exposures = sorted(summary.expnum.unique())
    table = group_table(summary)
    table.to_csv(HERE / 'residual_group_summary.csv', index=False)
    budget_figure(table, exposures)
    telluric_figure(telluric, exposures)
    lowpass_map(FOCUS)

    medians = scale.groupby(['expnum', 'channel'])[['rms_z', 'high_frequency_rms_z', 'low_frequency_rms_z', 'low_frequency_chi2_share']].median().reset_index()
    rows = []
    for e in exposures:
        item = {'Exposure': e}
        for ch in 'BRZ':
            m = medians[(medians.expnum == e) & (medians.channel == ch)].iloc[0]
            item[f'{ch}: total / noise-scale / smooth'] = f'{m.rms_z:.2f} / {m.high_frequency_rms_z:.2f} / {m.low_frequency_rms_z:.2f}'
        rows.append(item)
    scale_table = pd.DataFrame(rows)
    scale_html = table_html(scale_table, list(scale_table.columns), {c: str for c in scale_table.columns})

    pivot = table.pivot_table(index='expnum', columns=['channel', 'group'], values='chi2_per_pixel')
    share = table.pivot_table(index='expnum', columns=['channel', 'group'], values='pixel_share')
    detail = '<table><tr><th>Exposure</th>' + ''.join(f'<th>{ch}: {html.escape(name)}</th>' for ch in 'BRZ' for name, _ in GROUPS) + '</tr>'
    for e in exposures:
        detail += f'<tr><td>{e}</td>' + ''.join(f'<td>{pivot.loc[e, (ch, name)]:.2f} <small>({share.loc[e, (ch, name)]:.0%} of pixels)</small></td>' for ch in 'BRZ' for name, _ in GROUPS) + '</tr>'
    detail += '</table>'

    section = f'''<div id="residual-budget"><h2>Where the Sky-fiber residual comes from</h2>
<p>Every selected Sky fiber is rendered from its stored coefficients and the delivered block LSF. Each valid pixel (5 Å channel edges trimmed) is assigned to exactly one category by its distance to, and the family of, the nearest strong model line (flux ≥1% of the fiber's strongest line). Pixels farther than 4 Å from strong lines are split into <i>solar structure</i> (model deviates from its own ≈20 Å running median by more than 1σ) and <i>smooth continuum</i>. The residual is r = (data − model)/σ with σ = 1/√IVAR from the input CFrame. The statistic is the mean r² per pixel; 1 would mean the residuals match the input noise.</p>
<img src="residual_budget_categories.png" alt="Mean squared normalized residual by pixel category, per channel and exposure">
<p>Values and pixel fractions (pooled over all selected Sky fibers; telescope groups SkyE and SkyW are in <a href="residual_category_metrics.csv">the per-fiber CSV</a>):</p>
{detail}
<h3>Noise scale versus smooth structure</h3>
<p>Median over Sky fibers of the root-mean-square r: total; noise scale (robust width of r after removing the running median, ≈ pixel-to-pixel noise); and the rms of the ≈20 Å running median of r (smooth model mismatch). An accurate model with correct IVAR gives about 1 / 1 / 0.</p>
{scale_html}
<p>Exposure 16819 B has noise scale 0.93 (consistent with IVAR) and a smooth component of 1.25σ; in Z the noise scale itself rises to 1.59 because sky-line and telluric structure is included in the pixel-to-pixel term. Other exposures have smooth components of 0.2–0.9σ. <a href="residual_scale_metrics.csv">Per-fiber values</a>.</p>
<h3>Telluric windows</h3>
<img src="residual_budget_telluric.png" alt="Median squared normalized residual in telluric windows per exposure">
<p>Windows are fixed wavelength ranges and include any sky emission inside them; they locate the problem but do not separate transmission errors from line errors. <a href="residual_telluric_metrics.csv">Per-fiber values</a>.</p>
<h3>Smooth residual for each Sky fiber</h3>
<img src="residual_budget_smooth_maps.png" alt="Running-median normalized residual versus wavelength and Sky fiber for three exposures">
<p>16819 and 60919 are the large-residual cases; 3469 is a low-residual control. Horizontal lines separate spectrographs; vertical dotted lines mark B/R and R/Z boundaries. Rows follow the stored order of the Sky-fiber coefficient table within each spectrograph.</p>
<p><a href="residual_budget_summary.csv">Pooled category table</a> · <a href="residual_group_summary.csv">Grouped table</a> · generator <code>residual_budget.py</code> in the repository.</p></div>'''
    (HERE / 'residual_budget_section.html').write_text(section)


if __name__ == '__main__':
    main()
