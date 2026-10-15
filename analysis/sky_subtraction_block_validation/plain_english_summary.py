"""Plain-language summary of the in-fiber versus block LSF comparison; every number is computed from the published CSV files."""
import html
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / 'analysis/sky_subtraction_validation'))
from report_fibers import page  # noqa: E402
from fiber_decomposition_test import EXPOSURES, gaia_status  # noqa: E402

OUTLIER = (40620, 714)
MAIN = '00_MAIN_sky_block_report.html'


def table(header, rows):
    return '<table><tr>' + ''.join(f'<th>{h}</th>' for h in header) + '</tr>' + ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in rows) + '</table>'


def main():
    m = pd.read_csv(HERE / 'fiber_decomp_metrics_all.csv')
    c = pd.read_csv(HERE / 'fiber_decomp_categories_all.csv')
    lsf = pd.read_csv(HERE / 'fiber_lsf_profile_metrics.csv')
    keep = ~((m.expnum == OUTLIER[0]) & (m.fiberid == OUTLIER[1]))
    paired = m[keep].pivot_table(index=['expnum', 'channel', 'fiberid'], columns='variant', values='chi2_per_pixel')
    paired['increase'] = paired.block_refit / paired.fiber - 1
    status = {e: gaia_status(e) for e in EXPOSURES}

    rows = []
    for e in EXPOSURES:
        fibers = sorted({f for (ee, _, f) in paired.index if ee == e})
        clean = sum(status[e][f][0] for f in fibers)
        for ch in ('R', 'Z'):
            g = paired.loc[(e, ch)]
            rows.append([e, ch, f'{clean} of {len(fibers)}', f'{g.fiber.median():.2f}', f'{g.block_refit.median():.2f}', f'{100 * g.increase.median():+.1f}%',
                         f'{int((g.increase > 0).sum())} of {len(g)}'])
    chi_table = table(['Exposure', 'Channel', 'Fibers without star contamination', 'Typical χ² per pixel, own LSF', 'Typical χ² per pixel, block LSF', 'Change', 'Fibers where block is worse'], rows)

    cat = c[~((c.expnum == OUTLIER[0]) & (c.fiberid == OUTLIER[1])) & c.variant.isin(['fiber', 'block_refit'])]
    cat_rows = []
    for ch in ('R', 'Z'):
        grouped = cat[cat.channel == ch].groupby(['category', 'variant']).agg(n=('n_pixels', 'sum'), chi=('chi2_sum', 'sum')).reset_index()
        grouped['pp'] = grouped.chi / grouped.n
        wide = grouped.pivot(index='category', columns='variant', values='pp')
        count = grouped[grouped.variant == 'fiber'].set_index('category').n
        for name in ('OH core (<1 Å)', 'OH wing (1–2.5 Å)', 'Other line cores (<1 Å)', 'Other line wings (1–2.5 Å)', 'Outside LSF support (2.5–4 Å)', 'Continuum (>4 Å from lines)'):
            cat_rows.append([ch, html.escape(name), f'{int(count[name]):,}', f'{wide.loc[name, "fiber"]:.2f}', f'{wide.loc[name, "block_refit"]:.2f}', f'{100 * (wide.loc[name, "block_refit"] / wide.loc[name, "fiber"] - 1):+.1f}%'])
    cat_table = table(['Channel', 'Pixels near', 'Pixels (all exposures)', 'χ² per pixel, own LSF', 'χ² per pixel, block LSF', 'Change'], cat_rows)

    lsf['d_fwhm'] = lsf.fwhm_angstrom_block - lsf.fwhm_angstrom_individual
    shape_rows = []
    for ch in ('R', 'Z'):
        s = lsf[lsf.channel == ch]
        shape_rows.append([ch, f'{s.l1_distance.median():.3f}', f'{100 * s.l1_distance.median() / 2:.1f}%', f'{s.fwhm_angstrom_individual.median():.2f}', f'{s.d_fwhm.median():+.3f}',
                           f'{s.centroid_angstrom_individual.abs().median():.3f}', f'{len(s)}'])
    shape_table = table(['Channel', 'Median L1 distance', 'Share of the profile that differs', 'Typical FWHM, own LSF (Å)', 'Median FWHM change with block (Å)', 'Typical off-center shift of the own LSF (Å)', 'Profiles compared'], shape_rows)

    r = paired.groupby('channel').increase
    worse = {ch: (int((paired.xs(ch, level='channel').increase > 0).sum()), len(paired.xs(ch, level='channel'))) for ch in ('R', 'Z')}
    same = float(np.nanmax(np.abs(paired.fiber - paired.fiber_refit)))
    outlier = m[(m.expnum == OUTLIER[0]) & (m.fiberid == OUTLIER[1]) & (m.variant == 'fiber')].set_index('channel').chi2_per_pixel
    ratio = status[OUTLIER[0]][OUTLIER[1]][1]
    contaminated = {e: sum(not status[e][f][0] for f in status[e] if f in set(paired.index.get_level_values('fiberid'))) for e in EXPOSURES}
    body = f'''<h1>Sky subtraction with a block LSF: plain-language summary</h1>
<p><a href="../{MAIN}">Full technical report</a> · all numbers below are computed from the files in <a href="../data/fiber_decomp_metrics_all.csv">data/</a>.</p>

<h2>What we tested</h2>
<p>The spectrograph blurs every spectral line by a small amount. This blur is called the line-spread function (LSF). To remove the sky from a science spectrum, the sky model has to be blurred in the same way.</p>
<p>The standard method (<b>own LSF</b>) measures the LSF separately for every fiber, using that fiber's own spectrum. The new method (<b>block LSF</b>) measures one smooth LSF for each spectrograph from the sky-only fibers. It then looks up the LSF at each science fiber's position. Both are used in the red (R) and infrared (Z) channels. The blue channel (B) is not fitted and keeps the survey LSF.</p>
<p>We compared the two on nine science fibers (three per spectrograph) in each of five exposures: 45 spectra. For each spectrum we fitted the same sky model to the same pixels with the same noise estimate. Only the LSF differs. The science fibers were not used to build the block LSF.</p>

<h2>Main findings</h2>
<ul>
<li><b>The two LSFs are very similar.</b> Typically about {100 * lsf[lsf.channel == "R"].l1_distance.median() / 2:.0f}% of the R profile and {100 * lsf[lsf.channel == "Z"].l1_distance.median() / 2:.0f}% of the Z profile differs in shape (table 3).</li>
<li><b>In R, the sky fit is practically unchanged.</b> The typical fit quality gets worse by {100 * r.median()["R"]:.1f}% and never by more than {100 * paired.xs("R", level="channel").increase.max():.1f}%.</li>
<li><b>In Z, the block LSF is slightly worse.</b> The typical change is {100 * r.median()["Z"]:+.1f}%, and the largest is {100 * paired.xs("Z", level="channel").increase.max():+.0f}%. The extra misfit sits at the sky-line positions. The continuum between the lines is as good as before (table 2).</li>
<li><b>The block LSF is not better on any spectrum.</b> Its fit quality is worse for {worse["R"][0]} of {worse["R"][1]} spectra in R and {worse["Z"][0]} of {worse["Z"][1]} in Z. Some loss is expected: each own-LSF fit was tuned on the very spectrum it is judged on, while the block LSF was not.</li>
<li><b>Fit quality is limited by something other than the LSF.</b> Exposure 16819 fits poorly with both methods, and 60919 moderately so (table 1). Re-solving the sky-component amounts for each LSF changes nothing (largest change in χ² per pixel {same:.0e}), so the amounts are not the limit either.</li>
</ul>

<h2>Table 1. How well the sky model fits the spectra</h2>
{chi_table}
<p>“χ² per pixel” is the average squared difference between data and model, in units of the measurement noise. A value of 1 means the leftover is as small as the noise; larger values mean visible misfit. “Typical” is the median over the fibers. Exposure {OUTLIER[0]} fiber {OUTLIER[1]} is left out: a star is in the fiber (Gaia flux is {100 * ratio:.0f}% of the observed flux), and its χ² per pixel is above 10⁴ in both channels ({outlier["R"]:.0f} in R, {outlier["Z"]:.0f} in Z) even with its own LSF. Star contamination means Gaia flux above 1% of the observed flux. Contaminated fibers in the table: {", ".join(f"{e}: {n}" for e, n in contaminated.items() if n)}.</p>
<img src="../figures/fiber_decomp_chi2.png" alt="Chi-square per pixel for every fiber and model variant">
<p>Each point is one fiber. Colors show four model variants: the stored fit with its own LSF, the block LSF with the stored component amounts, the own LSF with re-solved amounts, and the block LSF with re-solved amounts. Variants that share the same LSF lie almost on top of each other. The figure includes the star-contaminated fiber that the text leaves out.</p>

<h2>Table 2. Where the extra misfit comes from</h2>
{cat_table}
<p>Pixels are sorted by their distance to the nearest strong model line. All nine fibers of all five exposures are pooled, without the star-contaminated fiber above. The misfit grows near OH lines in Z (core {cat_rows[6][5]}, wing {cat_rows[7][5]}). In R the OH pixels change by {cat_rows[0][5]} (core) and {cat_rows[1][5]} (wing), and lines that are not OH change more: cores {cat_rows[2][5]} and wings {cat_rows[3][5]}. Those pixels are few, so the pooled value moves easily. The continuum does not change.</p>

<h2>Table 3. How different the LSF profiles are</h2>
{shape_table}
<p>L1 distance: 0 means identical profiles, 2 means no overlap, so half of it is the share of the profile that differs. FWHM is the full width at half maximum. The block LSF is centered by construction; the own LSF is free to be off-center, which absorbs small wavelength errors. Part of the L1 distance therefore comes from this shift and not from the shape. Each fiber is compared at three wavelengths per channel.</p>

<h2>Why the block LSF is a little worse</h2>
<p>Measured: the extra misfit is at the sky lines and not in the continuum, and it is larger in Z. We have not tested the causes. Possible contributors are:</p>
<ul>
<li>The own LSF can follow the real line shape of that spectrum, including errors in the wavelength solution. The block LSF cannot, because it is smooth along the slit and exactly centered.</li>
<li>The sky-fiber amounts used to build the block are fixed from a fit with the median LSF, so remaining errors in them reach the block LSF.</li>
</ul>

<h2>What this does and does not show</h2>
<ul>
<li>It is a test on five exposures with nine fibers each, not on whole exposures.</li>
<li>Exposures 16819 and 44808 have star contamination in all nine fibers, so their numbers are less reliable.</li>
<li>The comparison is fair for the sky-line pixels but cannot show that either LSF is the true one. Neither is ground truth.</li>
<li>Interactive spectra for all 45 fibers, the LSF maps and the profile plots are in the <a href="../{MAIN}#fiber-decomposition">full report</a>.</li>
</ul>'''
    (HERE / 'summary_plain_english.html').write_text(page('Block LSF sky subtraction: plain-language summary', body))
    print('written', len(body))


if __name__ == '__main__':
    main()
