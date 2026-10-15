# Plotting and scientific report rules

This is the project reference for LVM Sky Decomposition figures and reports. Read it before creating or changing a diagnostic. Apply the user's latest explicit request first; a permission granted for one report does not change defaults for other reports.

## 1. Report entry point and directory structure

Publish one clearly named main HTML at the report directory root, for example:

```text
00_MAIN_sky_subtraction_report.html
spectra/       # individual interactive spectral diagnostics
figures/       # PNG/SVG figures
assets/        # local Plotly, scripts, styles
data/         # tables, measurements, manifests
docs/         # supporting methodology and review records
archive/       # earlier report versions and superseded diagnostics
```

- The root contains exactly one HTML entry point. Put supplementary HTML and assets in the appropriate subdirectories; use descriptive filenames throughout.
- Keep the directory URL as a server file listing. Use no default index document, including `index.html` or `index.htm`; share the explicit main-report URL.
- Publish HTML and linked products under `$SAS_BASE_DIR/sdsswork/users/u6058164/<descriptive-analysis-directory>/`. Keep generators in the repository.
- Make every relative link work from its actual page directory, including links back from individual diagnostics. Use local plotting assets and verify HTTP access to the main page and representative dependencies.
- Keep an earlier version in `archive/` when preserving it. The published root is not an accumulation of drafts, screenshots, or loose assets.

## 2. Visible content and concise explanations

- Write scientific text, labels, captions, and tables in English unless explicitly requested otherwise.
- Keep text minimal: a short measured outcome, understandable figure labels, and only the methodology needed to interpret the result. Place evidence near the conclusion.
- Display every section directly. Use ordinary headings and anchor links; no `details`/`summary`, collapsed accordions, hidden tabs, or selectors that hide report sections.
- Show telescope and spectrograph groups as separate visible plots. A physical-fiber selector is allowed only when explicitly authorized for that report; it does not authorize telescope/spectrograph selectors.
- Label multi-value diagnostics by channel and region, define wavelength ranges next to the table, and show how the values determine the decision; avoid unlabeled numeric sequences.
- Use informative table headers. Distinguish exposure selection, Gaia verification, attempted fits, accepted fits, and usable measurements. Use date-only labels when dates suffice; preserve exact times in provenance or hover.

## 3. Spectral diagnostics

- Use interactive Plotly figures covering the full native wavelength range, normally B/R/Z together. Preserve masks and the unsmoothed data; add zoom capability rather than replacing full spectra with selected line windows.
- Give each diagnostic its own Plotly figure, its own associated residual view, and its own legend. A spectrum plus its residual is one diagnostic; unrelated fits and subtraction stages belong in separate figures.
- Use flux display units of **10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹**. Convert flux and errors numerically, and inverse variance inversely to the square of that conversion. Keep normalized residuals dimensionless.
- Draw thin line/component traces first, the observed data next, and the complete best fit last; keep modeled sky-line curves thin so they do not obscure the data.
- Identify the observed spectrum, complete fitted model, and relevant physical components. Include exposure and physical fiber; include telescope, spectrograph, slit position, and MJD when comparing fibers or epochs.
- State the residual sign and noise source. For `r = (observed − model) / σ`, identify whether σ comes from input IVAR or another variance model; distinguish observed noise from full sky-model uncertainty.
- Plot actual SFrame FLUX directly: it already contains sky-subtracted spectra. Show standard DRP and new SFrame FLUX on the same axes and common valid pixels when comparing subtraction.

### Residual conventions

Show the primary data−model residual directly on the fitted spectrum axis in flux units, shifted to **−5 RMS**, where RMS is the unweighted root-mean-square residual over valid fitted pixels. Draw thin black curves at that offset ± the wavelength-dependent input flux error, without shading. Retain a separate aligned (data−model)/σ panel as an additional diagnostic; it does not replace the primary flux-space residual.
- **Normalized residual panel:** show zero and clearly labeled −1σ/+1σ reference lines and the ±1σ envelope. An additional ±3σ envelope is optional. Keep labels readable with pixel offsets when the displayed residual range is large. If shading is used, it supplements the visible reference lines.
- **Flux-offset residual view:** plot `observed − model + offset` on the spectrum's flux axis. Draw two thin black curves at `offset − input_error` and `offset + input_error`, with a stated offset and no uncertainty shading. This is the original main-project SolarBlue comparison convention.
- A panel showing remaining object FLUX after sky subtraction is not a zero-target fit residual. Label it as source plus sky residual; its input-error envelope follows wavelength-dependent σ in the displayed flux units.
- Preserve outliers. Any display clipping must be labeled and must retain an accessible full-range view; distinguish masked points from zero-valued measurements.

## 4. Titles, legends, interaction, and readability

- Reserve separate space for figure titles, subplot headings, axes, legends, and colorbars. They must not overlap. Each Plotly figure owns its legend; legends may sit in reserved space below the diagnostic.
- Verify the rendered layout at the actual browser width and a narrower desktop width. Source configuration alone does not establish visual correctness.
- Use compact axis labels and readable, slightly smaller tick fonts when needed. Increase figure height to accommodate required panels rather than layering information on conflicting axes.
- Use clickable legends and compact closest-point hover. Set `hovermode="closest"` and `showspikes=False` on both axes; never use unified hover, vertical hover lines, spikes, or large hover blocks covering data.
- Set Plotly `scrollZoom=False`: the mouse wheel scrolls the page even over a plot. Keep drag-to-zoom, pan, reset, and toolbar controls.
- Provide full-range zoom/pan and reset. Do not depend on WebGL for basic spectral traces when SVG traces suffice.

## 5. Whole-exposure maps and color scales

- Compare complete science-fiber SFrame FLUX as wavelength × physical-fiber maps, with standard and new subtraction side by side. Keep B/R/Z clearly separated and use the same valid-pixel mask, row order, and color limits within each comparison pair.
- Give line-window overview maps broad spectral context, normally 400–1000 Å; retain the full-channel view and interactive zoom for individual lines.
- Preserve the physical order from SLITMAP/pseudo-slit metadata. Do not concatenate or reorder telescope groups arbitrarily. Mark telescope/spectrograph boundaries from verified metadata.
- Use an informative multi-color sequential scale for positive LSF quantities such as widths. Use a zero-centered diverging scale for signed flux/residuals and signed LSF differences.
- Label colorbars with quantity and units. State whether scales are shared across panels/epochs and state any percentile limits or saturation. Masked or unknown values must be visually distinct from zero.
- For LSF maps, show the requested wavelength-dependent metrics, including FWHM, W50, W90, centroid, and asymmetry where relevant. A compact overview of one metric across all channels/spectrographs is preferable to many enormous, unreadable 1D plots.
- For side-by-side Science/SkyE/SkyW summaries, give Science substantially more width; the established coverage layout uses approximately 2:1:1. Mark spectrograph boundaries using the actual mapping.

## 6. LSF trends, levels, and scatter

- Show physical fibers separately in time-trend figures. Do not overlay different fibers as if they were repeated observations of one profile.
- Show metric level (mean/median) and scatter in separate paired panels, each with one Y-axis. Do not put level and scatter on dual Y-axes.
- Show absolute values/scatter in Å as well as requested percentages. Define every percentage denominator and compute levels, scatter, and denominators from the same stated plotted selection.
- Show individual measurements, observing-night summaries where appropriate, and actual coverage. A single observation is a measurement, not a stability estimate. Unknown coverage or Gaia status is not zero.
- Preserve true dates, full MJD values, sampling gaps, and calibration epochs. Calendar dates may accompany MJD axes. If lunations have equal-width blocks or exposure-rank coordinates, label that coordinate explicitly; it is not elapsed time or lunar phase.
- Use local observing windows for sparse series when appropriate; show every retained point and state window selection. Keep comparisons on common scales where their purpose requires it.
- For signed pairwise LSF changes, use later minus earlier, the same offset grid and normalized profiles, and a common symmetric scale. State any alignment; do not silently recenter profiles.

## 7. Evidence for in-situ LSF quality

- Show the observed spectrum, full candidate fit, components, and spectral residuals. Kernel shapes and iteration-width plots supplement these diagnostics; they do not establish recovery quality.
- Identify survey start, fitted candidate, and adopted LSF separately. Label the last trial as a candidate, especially when survey LSF is adopted.
- Distinguish the auxiliary LSF-fitting sky and source spline from the actual final ML SKY subtracted. Show actual additive corrections and final remaining FLUX when diagnosing subtraction.
- Use explicit labels such as “LSF fitting in R and Z channels”, “Fitted sky”, “Nuisance continuum spline”, and “Sky + continuum fit”; explain a candidate only when its adopted status matters. State that solar absorption is part of the sky continuum and that a nuisance spline is not a physical stellar-spectrum estimate.
- Give selected individual-fiber examples a compact table with selection reason, measured evidence, and an ordinary low-residual control. Identify a suspected artifact by fiber and wavelength; keep whole-exposure titles descriptive of observing conditions.
- Separate solver success, numerical convergence, physical-profile guards, spectral fit quality, and independently demonstrated recovery accuracy. Runtime and acceptance percentage are not accuracy measures.
- A claim about fit evolution requires the per-iteration spectral model, continuum, masks, and weights, or an exact reproducible reconstruction. Saved LSF widths alone do not establish spectral improvement.
- Compare initial/reference and refined fits under matching nuisance-parameter and weighting conventions. State held-out-line or injected-truth evidence when claiming LSF recovery accuracy.

## 8. Before delivery

- Verify array units, residual arithmetic, masks, physical ordering, complete wavelength coverage, and the actual compared products.
- Inspect rendered spectra, residual ±1σ references, legends, titles, annotations, and colorbars; inspect both an adopted and a survey-fallback example when applicable.
- Obtain at least two independent fresh-context reviews of a scientific notebook or HTML report covering first-reader comprehension, terminology, methodology, and visual readability. Incorporate concrete feedback before delivery.
- Verify the directory listing, the single named main HTML, relative links and representative HTTP assets. Share the explicit main-report link.
- Keep the final message short: main report/rules link, concrete changes, and material scientific limitations.

## Sources and precedence

Consolidated on 2026-10-15 from the user's corrections in the current SFrame chat, main-project `AGENTS.md`, `analysis/solarblue_lsf_serious_corpus/README.md`, `analysis/lsf_epoch_coverage/README.md`, and the reviewed main-project LSF coverage, survey-quality and stability chats. The latest SFrame request requires a primary flux-space residual at −5 RMS with input-error curves on the spectrum axis, with a separate normalized-residual panel as an additional view. Component traces are drawn first and thin, then data, then best fit. Mouse-wheel zoom is disabled. Report-specific scientific selection thresholds remain in their analysis documentation rather than becoming universal plotting cuts.
