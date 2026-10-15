# CFrame sky subtraction: implementation plan

Status: implemented and exercised on five DRP 1.3.2 exposures, 2026-10-14. Code, CLI, bundled ASCII references, tests, and the reviewed validation report are present. Keep this plan uncommitted. See `analysis/sky_subtraction_validation/index.html` for measured results and scientific limits.

## Goal and working rules

Build a small Python package and CLI in `skysub/sky_subtraction/` that reads an arbitrary supported `lvmCFrame`, models the sky using SkyE/SkyW and the existing trained predictor, subtracts it from the science fibers, and writes an SFrame-compatible FITS beside the input. Preserve the data needed to reproduce and assess every stage. Validate on approximately five exposures, including difficult stellar contamination.

Use Ponytail throughout implementation and testing: reuse existing numerical code, add only the CFrame adapter, orchestration, per-fiber R/Z fitter, and output/QA needed for this workflow. Keep one readable pipeline, ordinary functions, and the existing dependencies. Add focused tests for new numerical behavior. Scientific documentation, labels, tables, and CLI help are in English. Do not retrain the predictor as part of this task.

Provisional output name: `lvmSFrame-skyml-<exposure>.fits`. Keep the suffix in one place so it can be finalized later. Always protect the input and standard DRP SFrame; refuse an existing destination by default and publish the completed FITS atomically.

## Verified starting points

| Responsibility | Existing code to reuse |
| --- | --- |
| Gaia ARI TAP query and fiber photometry | `skysub/medians_computation/src/lvm_medians/gaia.py` |
| Gaia selection, masks, medians, fiber identities, metadata | `skysub/medians_computation/src/lvm_medians/stack.py`, `metadata.py` |
| Continuous spline LSF and iterative decomposition | `skysub/sky_decomp/lsf_spline2d.py`, `lsf_surface_iterative.py` |
| Calibration loading and selected-fiber reference averages | `skysub/sky_decomp/lsf_reference.py` |
| Fit serialization | `skysub/sky_decomp/result_io.py` |
| Trained ensemble loading and prediction | `skysub/mlp_predictor/serialization.py`, `inference.py` |
| Sky-arm residual transfer | `skysub/mlp_predictor/sky_arm_correction.py` |
| OH amplitude and rotational-temperature correction | `skysub/mlp_predictor/sky_line_scaling.py` |
| Integration example and deployed ensemble selection | `skysub/notebook_example_predict_sky.ipynb` |

The Gaia code currently assumes an SFrame: observed total flux is `FLUX + SKY`. For CFrames, inspect the actual HDUs and use the observed pre-subtraction flux directly. Adapt the input boundary without adding sky twice or requiring a companion SFrame. Preserve existing SFrame callers and their tests.

### Survey LSF dependency

The related chat `01a13ac6-5bfb-7221-af48-fab3a913d581` is actively producing quality-screened ASCII references. On inspection, its worktree at `/home/ik52/.codex/worktrees/blue-lsf-quality/lvmsky/` already contains:

- `skysub/sky_decomp/data/lsf/survey_lsf_{B,R,Z}.csv`;
- `wavelength_knots.csv`, `coverage.csv`, `manifest.json`, and `README.md`;
- an ASCII loader being integrated and checked on real spectra.

The user selected `skysub/sky_decomp/data/lsf/` in that chat. Use that location. These files are not yet present in this checkout. Recheck completion and bring the validated files and necessary loader changes into this checkout before implementation; do not duplicate the active conversion. If that work is unavailable, convert the survey product under `/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/lsf_calibration` with the same full-profile representation, provenance, quality flags, and numerical round-trip checks.

“Survey-wise” means one wavelength-dependent reference **per physical fiber and channel**, pooled over the survey; it does not mean one shared profile for all fibers. Preserve knots, normalization, displacement support, units, native grid, and quality metadata. Read the manifest rather than hardcoding defaults: the new reference support differs from the default in the current spline module. Z fibers 84/85 currently have diagnostic warning profiles rejected by strict loading; define an explicit flagged outcome for an unusable fallback.

## Pipeline and intermediate contracts

### 1. Read and validate one exposure

Read FLUX, IVAR, MASK, WAVE, SLITMAP, observation metadata, and available DRP provenance. Verify dimensions, physical units, wavelength convention and grid, fiber IDs, and telescope labels. Join all tables by physical `fiberid`, never by accidental row position. Use the DRP fiber-status and pixel-mask rules. A usable pixel has valid status, finite flux, and finite positive IVAR.

Inspect a real CFrame and matching DRP SFrame first to fix the output HDU contract, arm-overlap handling, and any differences in uncertainty semantics. Support the requested input without relying on a pre-existing SFrame. Report unsupported grids or incomplete essential metadata explicitly.

### 2. Resolve Gaia and select fibers

Use the supplied Gaia path when present. The reference cache is `/data/sas/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/gaia/fibers/`. Follow the requested behavior: if no path is supplied, or the exposure is absent/incompatible, query ARI TAP using the existing Gaia DR3 query, aperture radius 17.65 arcsec, source columns, and retry/cache behavior. Record the query and cache provenance. Distinguish a successful query with zero matches from a query failure.

Reuse the existing G-passband metric: summed Gaia G flux divided by observed total synthetic G flux. Match `gaia_ratio_clean_mask(..., threshold=0.01)`, including positive observed flux, at least 95% valid passband coverage, complete source photometry, and the zero-source case. Cached photometry must correspond to the same exposure and spectral convention; recompute the denominator from the CFrame when needed.

Select eligible SkyE, SkyW, and Science fibers independently. For Science, if fewer than 50 satisfy the threshold, take the 50 eligible fibers with the smallest valid ratio, using a stable fiber-ID tie break. Successfully queried zero-source fibers have ratio zero for ranking. Unknown photometry is not a clean fiber. If fewer than 50 eligible fibers exist, use all eligible fibers with an explicit reduced-sample flag; require enough usable data for the downstream fit. SkyE and SkyW use only the Gaia-clean selection; the 50-fiber fallback applies only to Science (confirmed by the user).

Store every fiber's metric, eligibility, selection, selection reason, and Gaia status. The threshold is a stellar-contamination screen, not proof that the Science median contains only sky.

### 3. Build three masked medians

Produce SkyE, SkyW, and Science median fluxes on the input grid, pixel contributor counts, masks, uncertainty estimates, and selected fiber IDs. Exclude invalid pixels before taking medians; mark unsupported wavelengths. Use an explicitly documented uncertainty estimate for the median, checked against resampling/noise simulations; do not substitute the variance of an arithmetic mean.

Construct the effective survey LSF from the selected physical fibers, reusing reference averaging. An average LSF is only an approximation for a median of different profiles and wavelength-dependent contributors: test this approximation with synthetic medians and store its selection/weight convention. Keep the final per-fiber Blue LSF exact.

### 4. Decompose the medians

Fit all three medians with the existing sky decomposition. Fix Blue to its selected-fiber survey reference; initialize R/Z from the corresponding survey references and allow them to refine. Explicitly wire R/Z survey initialization: the related chat's existing fixed-Blue integration currently describes DRP initialization for R/Z.

Use the deployed model's coefficient names, ordering, spline basis, normalization, and moon/zodiacal conventions. Verify these against the ensemble archive, not just coefficient count. Preserve fitted components, coefficients, residuals, LSF states, masks, fit statistics, and convergence flags. Science decomposition provides diagnostics, the median LSF reference, and correction inputs; its coefficients do not replace the ML prediction.

### 5. Predict the Science sky and estimate corrections

Load the existing ensemble once. Use `predict_sky_from_minimal_inputs` with the SkyE/SkyW coefficients and required geometry, time, LSF, and exposure metadata. Preserve its near/far association and model-context conventions. Save raw predicted coefficients, ensemble spread, and reliability diagnostics.

Apply the existing sky-arm residual transfer using each Sky telescope's residual against its own fit. Then derive the existing OH band-amplitude and temperature-tilt correction from the Science median after that transfer. Reuse the science-line protection, robust weighting, and velocity handling already implemented; do not replace them with an unrestricted residual fit. Save correction parameters, masks, gates, and diagnostics separately from the raw prediction.

Predict once at the Science IFU center, with common median-derived OH parameters and a separate rendering LSF for each science fiber (confirmed by the user). Use the Science median LSF for the initial prediction context and record it. The subsequent per-fiber LSF fit refines the rendering; avoid an unrequested iterative ML/refitting loop.

### 6. Develop the per-fiber R/Z LSF fit

This is the main new numerical component. For each eligible science fiber, fit R and Z using the same continuous LSF surface representation, constraints, and line rendering as decomposition. Initialize with that fiber's survey reference. Use the predicted and OH-corrected sky-line spectrum as the starting line model and a single flexible spline continuum per fitted channel as a nuisance component.

Fit the LSF jointly with the nuisance continuum and only the minimal line-amplitude freedom needed to avoid interpreting amplitude mismatch as width. Choose continuum knot spacing through injection tests so it follows broad stellar structure without absorbing instrumental sky-line profiles. Preserve science-line masks, handle masked/cosmic pixels robustly, and assess usable sky-line information. The nuisance continuum is never subtracted as sky. Any local amplitude adjustments used solely to estimate LSF remain diagnostic unless explicitly promoted by a later scientific decision.

Assess each channel using convergence, supported line coverage, physical profile constraints, residual behavior, bound saturation, and profile metrics versus both its fiber survey reference and telescope median fits. Compare evaluated widths, wings, and centroids across wavelength, rather than raw coefficient distance alone. Calibrate numerical thresholds with controlled tests and freeze them before the final exposure validation; genuine inter-fiber differences must remain possible.

Accept the fitted profile or fall back to that fiber's survey profile independently in R and Z. Store candidate and adopted coefficients, quality statistics, and reason. If the survey profile is itself unusable, preserve the affected data with an explicit invalid-model flag rather than inventing a calibrated profile.

### 7. Render and subtract once

For each science fiber, render the predicted sky with fixed survey Blue and adopted R/Z profiles. Apply the median OH amplitude/temperature parameters to line templates regenerated at that fiber's LSF; do not reuse a median-resolution correction spectrum for narrow lines.

Transfer the sky-arm residual correction consistently to the target model. Inspect the existing correction's wavelength support and smoothing: an empirical residual already convolved with a Sky LSF cannot simply be reconvolved as if it were an intrinsic spectrum. Validate the current transfer at differing resolutions; handle any narrow-line contribution explicitly and document the supported transfer rule.

Assemble `sky = predicted components + residual correction + OH correction` using a clearly recorded order, then `science_output = observed - sky`. Each correction is applied once. Store the actual applied arrays so the result can be reconstructed exactly.

### 8. Write an auditable SFrame-compatible result

Preserve the standard spatial layout, metadata, wavelength grid, and necessary DRP HDUs, updating sky-subtraction semantics consistently. Define non-science-row behavior from the inspected DRP contract and record it. Preserve input flux in a dedicated extension if needed for reconstruction.

Additional extensions should contain, with concise documented names:

- per-fiber Gaia photometry/matches and selection flags; retain source matches when available and record aggregate-only cache provenance;
- median spectra, uncertainties, contributor counts, masks, fit coefficients/components/residuals, and median LSF states;
- raw ML coefficients, spread, reliability, feature/model provenance, and correction parameters;
- fitted and adopted per-fiber LSF coefficients, common knots/configuration, fallback reasons, and fit statistics;
- applied sky and correction arrays, calibration/model hashes, software version, input identity, and settings.

Reuse existing LSF/result serialization where compatible. Verify round-trip reconstruction. Specify IVAR and SKY_IVAR semantics: keep observational noise and estimated sky-model uncertainty distinguishable, avoid double-counting the original sky noise, and do not label ensemble spread as full predictive variance. Mark approximation and unmodeled covariance explicitly.

## Small implementation sequence

1. Inspect one real CFrame/SFrame pair, adopt the finalized ASCII references, and lock input/output and model compatibility contracts.
2. Implement CFrame/Gaia adaptation, selection, medians, and their focused tests.
3. Connect existing median decomposition, ML inference, residual transfer, and OH correction in one exposure pipeline.
4. Implement and validate the R/Z fiber fitter and channel-specific fallback.
5. Complete per-fiber rendering, FITS writing, CLI, and round-trip checks.
6. Run the five-exposure validation and deliver readable QA with reproducible commands.

Start with `pipeline.py`, `fiber_lsf.py`, `__main__.py`, and a small `tests/` directory; split I/O only if it materially improves readability. Review callers and tests whenever a shared helper changes. Keep existing notebook workflows operational.

Proposed CLI contract: `python -m skysub.sky_subtraction INPUT --ensemble MODEL [--gaia-dir DIR] [--lsf-calibration DIR] [--output PATH]`. The packaged ASCII reference is the calibration default. Add only operational options actually needed for reproducible execution. Log stages, counts, elapsed time, and fallback summaries.

## Validation and completion criteria

Focused tests cover shuffled fiber IDs; cached/TAP Gaia equivalence; zero matches versus unavailable photometry; the 0.01 threshold and 50-fiber fallback; masked medians and uncertainty; fixed Blue; survey R/Z initialization; model-basis compatibility; correction ordering; per-channel LSF fallback; FITS round trips; and overwrite protection. Re-run affected existing tests.

Controlled LSF tests use known R/Z surfaces, realistic noise/masks, and injected smooth continua, emission lines, and M-star absorption at several source/sky contrasts. Measure profile/width recovery and sky subtraction separately. Require hard-case contamination to lead to accurate recovery or an explicit fallback, and measure science-feature flux preservation. Include a case where identical raw coefficient distances correspond to different profile behavior. Set quantitative acceptance limits from these tests before assessing the final sample.

Choose approximately five exposures spanning a clean dark field, bright moon, strong/diffuse science emission, a bright stellar field including M-star contamination, and a challenging LSF epoch or low-S/N case. Use local data first; fetch needed inputs through the existing `getsas` workflow. Record exposure IDs and the reason for each choice. Separate threshold tuning from the final validation sample and identify any overlap with predictor training.

Compare against the matching DRP SFrame when available and against the survey-LSF-only version of this pipeline. Report per-channel sky-line residual bias/RMS, continuum residuals, science-line preservation, LSF diagnostics, fallback fraction, runtime, and memory. Include ablations for residual and OH corrections where they help attribute improvement. Five exposures establish an initial demonstration, not survey-wide generalization. Use compact closest-point Plotly hover, no unified hover, and no axis spikes.

Deliver the code/CLI, calibration data, tests, five output files, run manifest, and concise QA documentation with exact reproduction commands. If a notebook or published HTML report is produced, obtain two independent fresh-context reviews of comprehension, terminology, methodology, and visual readability, and incorporate feedback before delivery.

## Confirmed decisions and deferred naming

1. Confirmed: the 50-fiber fallback applies only to Science. SkyE and SkyW retain the Gaia-clean selection and explicit fit-quality gates.
2. Confirmed: one ML prediction at the Science IFU center, with common median-derived OH parameters, rendered with each fiber's own LSF.
3. Final product filename can be settled later; the provisional `skyml` suffix guarantees separation from DRP products.

Selection of continuum flexibility, LSF acceptance thresholds, uncertainty approximations, and residual-resolution handling are implementation research tasks with tests, rather than arbitrary constants to ask the user to invent.
