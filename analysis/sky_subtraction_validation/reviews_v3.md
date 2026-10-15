# Independent reviews: whole-exposure maps and LSF iterations

Two independent fresh-context reviewers inspected the generated report, source, paired FLUX maps, and iteration images. One focused on first-reader comprehension and visual readability; the other examined scientific methodology and the quality measurements.

Both found the direct FLUX comparisons, shared paired scales, common masks, all-science row coverage, and iteration difference plots clear.

Applied refinements:

- Magenta iteration profiles identify the last fitted trial. When the trial is not adopted, the survey-start curve explicitly identifies the adopted profile.
- The report explains that maps include every science row, while quality statistics use fibers with SKYML applied and at least ten common eligible line pixels.
- LSF is expanded as line-spread function.
- The quality statistic is explicitly RMS of high-pass FLUX normalized by shared per-pixel observed noise, with its formula and per-fiber/channel sky-line threshold definition.

The report states that SKYML has lower median line-scale RMS in 12 of 15 exposure/channel comparisons, with higher values in 3469 R/Z and 59750 Z. This statistic measures residual fluctuations together with source structure and noise; it does not establish absolute sky accuracy or isolated LSF accuracy.

Paired-mask and quality-score tests passed. Browser rendering was checked separately from the reviews. The pipeline and output FITS are unchanged by this presentation revision.
