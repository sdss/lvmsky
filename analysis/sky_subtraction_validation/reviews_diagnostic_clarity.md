# Diagnostic clarity reviews — 2026-10-15

Two independent fresh-context readers reviewed the auxiliary LSF-fitting diagnostic.

- Methodology/terminology review: separate the fitted emission spectrum, solar/diffuse spectrum, nuisance continuum, and actual final SKY; describe absorption as part of the solar spectrum. Define contrast aggregation, peak separation and variable peak count. Label low contrast for inspection. Applied in generators and linked pages.
- Rendered first-reader review: checked the fiber 59750/680 screenshot and source. Titles and legends occupy separate areas; full-range interaction and input-error normalization are clear. Requested an explicit RMS definition and inward ±1σ labels. Both applied.

Presentation review does not establish independent LSF-recovery accuracy. Flux products and LSF acceptance rules remain those of the existing v2 run.

## Selection scheme and broader windows

Two new independent fresh-context reviews checked the selection scheme against fiber_lsf.py/pipeline.py and a 1400-pixel rendered main-page screenshot. Applied: exact basis-weighted coverage wording, measured-survey fallback with an invalid-channel branch, explicit distinction between local nuisance fits and final ML rendering, and definitions of LSF/W50/W90. Layout was readable with visible branches and separate text areas. Hα-to-R-end and telluric-A context windows cover 6500–7454 Å and 7454–8050 Å.

## Observed full-channel information gate

Two independent fresh-context reviewers checked the new gate and report source. The first verified local observed line-amplitude fits, input-error S/N, per-quarter coverage, 10 Å peak separation and conservative thresholds. It identified a stale hard-coded overview count; that sentence was removed. The second requested exact contrast definitions, readable adoption reasons and an explicit distinction from the older five-peak summary; all were applied. The preparation block explicitly performs no fit, and 20 Å denotes continuum knot spacing. The updated gate is a conservative operating point, not an independently calibrated recovery probability.

## Compact parameter-based products

Two independent fresh-context reviewers inspected compact.py, reconstruction.py, writer integration and report callers. Applied findings: empty iteration histories produce unavailable fit curves with an all-excluded mask; reconstructed validity includes finite input IVAR; captions identify reconstructed median models; current callers/fixtures use SKY_ML; duplicate static fiber-report code was removed in favor of the interactive diagnostics. Native robust-weight reconstruction retains the pre-LSF spline and Huber scale. Real-product checks compare all three median models and representative individual SKY, corrections, fit components and weights to the archived full-diagnostic products. Additional 59750 fibers 680/1447/971 pass strict <2e-6 peak-relative checks. Full SKY can render with an unavailable CFrame path. Final browser layout is inspected separately.

Final rendered review: both independent reviewers inspected browser_compact_main.png and browser_compact_fiber.png. Titles, legends, axes, input-error references, tables and continuous visible layout passed; the contract accurately distinguishes reconstruction verification from independent scientific accuracy.


## Median-fit investigation and local decision evidence

Two independent fresh-context reviewers checked the median RMS/error definitions, local LSF decisions, example selection, broad windows, and rendered median/table/blue-window images at 1400/1050 pixels. Applied: explicit unweighted fit RMS; adjacent input-noise metric definition; line-amplitude S/N terminology; Gaia units clarification; zero-shift local control with identical amplitude/baseline/slope freedom; displacement sign and peak selection. A rounding issue made 99.54% appear as 100%; adoption now comes directly from current LSF_QA with one decimal and the full processed-fiber numerator/denominator. Rendered titles, legends, axes, colorbars and continuous visible content are readable. The investigation identifies actual systematic mismatch; these checks do not establish independently calibrated LSF accuracy or a scientific repair.

Final review confirmed fixed-LSF coefficient comparisons, matched-fiber subtraction and actual OH applied/retained provenance are distinguished. Causal wording is limited to the repeated exposure; the possible benefit of retained OH brightness terms requires a constrained test. The final local-probe screenshot preserves readable zero/free-shift controls, selection, units and displacement sign.


## Explicit sky-line coverage counts

Two independent fresh-context reviewers checked the revised table source and its 1150-pixel browser rendering. Each channel now has its own labeled line, Q1–Q4 feature counts, passing-region total, and adopted-LSF decision; the adjacent legend gives native-grid wavelength ranges. Both reviewers verified region edges against the production information gate and found the layout readable without overlap. All 80 displayed channel records have measured coverage. The wording distinguishes coverage from convergence/profile checks and from independently demonstrated LSF accuracy.

Caption follow-up: two fresh independent reviewers checked the visible table caption, worked count example, explicit per-cell threshold/adoption summaries, and browser_counts_caption.png. Both found the wording and rendering immediately readable; the caption retains convergence/profile requirements beyond coverage. The targeted table test passes.


## Exposure-level survey LSF statistics

Two fresh independent reviewers checked report_lsf_summary.py, its two tests, current FITS QA records and the rendered table. Both verified processed Science selection, separate R/Z totals, primary reason counts, unique-fiber union and invalid-channel separation; all five FITS summaries match the exported JSON. Applied wording: R and Z are assessed separately, avoiding a statistical-independence claim. The table was visually checked at 1400 and 1050 pixels with fully visible content and readable columns. Counts summarize adoption decisions, not independently demonstrated LSF recovery accuracy.

Channel-first summary: two fresh independent reviewers checked separate R-channel and Z-channel tables, processed denominators, in-situ/survey counts and percentages, reason totals and the 1400-pixel rendering. Both found the layout clear and readable; the two targeted aggregation/rendering tests pass. Unique R-or-Z counts remain in the linked JSON.
