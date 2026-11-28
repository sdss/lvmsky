# Changelog

## Unreleased

### Changed

- Block `single`, `batch` and the legacy entry point accept `--no-bubble-qa` (2026-11-28). The plasma-bubble QA
  flag ran in every block run with no way to turn it off, and offline each uncached GNSS request waits for its
  timeout. The switch reaches `sky_subtraction.pipeline.run(bubble_qa=...)`, gives `BUBBLE = -1`, leaves the sky
  unchanged and enters the product fingerprint only when given. Checked on exposure 47284.

- O I integration review (2026-11-28): offline exposure context now passes `download=False` through to IONEX;
  replacing `ATOM_Or` without a valid test-night scatter records an unknown uncertainty rather than the old
  MLP spread; standard block `single` records unavailable CPU affinity on macOS without crashing. Added
  regression checks for these paths and coefficient replacement/fallback. Product code hashes now include
  the ionosphere, thermosphere, bubble-prior and ROTI implementations.

- SkyML reads the dedicated Science [O I] 6300 model, and the model is rebuilt by one command.
  - `mlp_predictor.oi6300_build` now starts from the decomposed median stack: it builds what is missing of the
    ionospheric context (`iono_corpus.build_context`: IONEX map table, GNSS station-days and ROTI table,
    thermosphere table, corpus; a rebuilt file forces the ones after it, `--rebuild-context` rebuilds all), then
    selects the rows, validates by nights and saves the model. Without options it builds the default model
    (`--sci-lines` defaults to `oh o2 5577 7774 8446`, the night field is on; `--no-sky-field` and `--sci-lines`
    without values give the other variants). The same command is `python -m skysub.sky_subtraction_block.cli
    oi-model`. `iono_corpus.build` reads the three tables of its own stem instead of the merged tables of every
    stem; the rebuilt knox corpus equals the stored one in all 146 data columns (`in_training_set` aside).
  - A saved model holds `corpus`: for every corpus exposure the night field and the model's prediction, the
    test-night prediction for the rows it was trained on (`OI6300Model.corpus_row`, `corpus_coefficient`).
  - `mlp_predictor.oi6300_exposure` (new) applies a saved model to single exposures: it builds the input table from
    the three median decompositions, the observing context and `iono_corpus.exposure_columns` (the corpus columns
    computed for any exposure). On 360 knox exposures it gives the prediction of the build table within 1e-4 dex for
    95% of them (largest difference 0.0015 dex). A night-field model gives no value outside its corpus.
  - SkyML and the block: `--oi-model FILE` (`pipeline.run(oi_model=...)`). The model runs after the median fits; its
    `ATOM_Or` replaces the ensemble coefficient (uncertainty: the model's test-night scatter), and the ensemble
    coefficient is kept, with the reason recorded, when the model gives no value. New keywords `OI63SRC`, `OI63R`
    and `SKYML_META.oi6300` (`FITS_STRUCTURE.md`). The block fingerprint includes the model file only when the
    option is given. Checked on exposure 47284 (DRP 1.3.2): no-field model 138.2 and night-field model 145.4
    against the ensemble's 143.8 and the Science median's own 147.4; the no-field value agrees with the corpus
    prediction of that exposure within 0.0004 dex.
  - Sky-transfer ensemble inputs: `data.augment_ionospheric_context` appends the ionospheric context the model reads
    and, with `oi_model`, the columns `iono_sqrt_oi6300_model` and `iono_oi6300_model_avail`
    (`data.OI6300_MODEL_NAMES`); `inference.predict_sky_from_minimal_inputs(oi6300_coef=...)` builds them for an
    exposure. No ensemble has been trained with them yet.
  - The default and the no-field models were rebuilt with the new command; their test-night scores are unchanged
    (0.0055 and 0.0062 dex). The other four saved models have no `corpus` table.
  - Not done: the ionospheric context of a new stack was not built from scratch in one run (on knox the tables
    existed). The later integration review fixed the macOS affinity-reporting failure and ran standard
    `skysub.sky_subtraction_block.cli single` successfully on exposure 47284 with the night-field model.

- `mlp_predictor.sky_field` (new) and `OI6300Model(sky_field=True)`: optional inputs from a night-long field of the
  zenith [O I] 6300 measured by SkyNear and SkyFar (Gaussian process over all exposures of the night, read at the
  Science and SkyNear pierce points; 4 field inputs plus the number of samples within 1 h). Ported from
  `experiments_external_data_v1/WP8_lvm_internal/skyfield.py`; the new module reproduces its values to 3e-5 dex
  (one gradient to 7e-4 in the worst row) when given its kernel parameters. The kernel parameters are learned from
  the sky telescopes of the corpus at build time (`sky_field.fit`, 80 random nights) and stored on the model as
  `field_hyper` and in its JSON, so a rebuild on another DRP version learns them again.
  Post-processing only: the field uses later exposures of the night. Default off; saved models load unchanged.
  `oi6300_build --sky-field` builds the inputs. Three more models in `skysub/oi6300_models/` (`*_field_ens5`),
  test nights, 4,834 rows: no Science lines 0.0083 -> 0.0073 dex; OH, O2, 5577, 7774, 8446 0.0062 -> 0.0056; plus
  5199 0.0056 -> 0.0051 (0.0072, 0.0056, 0.0051 with the experiment's parameters). The default model is now OH, O2, 5577, 7774, 8446 with the night field (`oi6300_model_oh_o2_5577_7774_8446_field_ens5`); the description is rewritten around it. Test added
  (`sky_decomp/tests/test_oi6300.py`, 8 pass).
- Science [O I] 6300 model, second review before the first commit (against the independent review in the Obsidian
  note `OI 6300 Uncommitted Review 2026-11-27`); all six models rebuilt. Fixes: exposures with a pointing at
  ra = dec = 0 no longer give samples to the night field; a missing ROTI stays missing (it was read as 0.005
  TECU/min); a sky telescope of zero weight is no longer read by the base (a missing SkyFar line gave NaN at
  `a = 1`); the cross-validation learns the OH groups, the kernel parameters of the night field and the input
  ranges inside each run from the other nights only; `sky_field.fit` picks nights by usable exposures, records
  nights, samples, likelihood and convergence, and treats overflowing parameters as a rejected point; `build_table`
  compares every corpus line with the coefficient files, rejects repeated exposure numbers, and checks the rows and
  grid of the decomposition files through residual + best fit = stack flux; repeated `sci_lines`, a non-positive
  guard or ensemble size and too few nights are rejected before the fit; `iono_corpus._match` takes the nearer of
  the two neighbours (a time 5e-8 d above an entry was not matched). New: `OI6300Model.lines_used`, `quality`,
  `INPUTS_VERSION` (`load` refuses a model built for other inputs: models saved before this entry must be rebuilt);
  each build also writes `<name>_cv.csv` (test-night prediction of every training row) and records input files
  with checksums, git commit, module checksums, settings and package versions in `<name>.json`. Test-night sigma
  after the rebuild [dex]: default (Science lines + night field) 0.0055 (was 0.0056), + [N I] 5199 0.0051 (0.0051),
  Science lines without field 0.0062 (0.0062), + 5199 0.0056 (0.0056), field without Science lines 0.0071 (0.0073),
  neither 0.0083 (0.0083). `OI6300_MODEL.md` rewritten as a step-by-step description of the data selection and the
  build. The corpus label of `h6300` is corrected to the emission-weighted mean height (values unchanged).
- `OI6300Model.predict`: a guarded Science line that is missing for an exposure is now handled like a line above
  the guard, with the model trained without it; before, it was passed to the full model as a missing value. No
  retraining: saved models pick the rule up on load. On the knox 1.3.2 corpus no clean row has a missing guarded
  line, so predictions on the 10,426 clean rows are bit-identical and the test-night scores stand; 1,523 of 18,640
  rows change (1,969 with [N I] 5199), median 0.002 dex, maximum 0.02-0.03 dex.
- `mlp_predictor.oi6300_build` (new): one module that rebuilds the Science [O I] 6300 model from a decomposed median
  stack: `python -m mlp_predictor.oi6300_build --decomp-root ROOT --stem STEM --out oi6300_models --sci-lines oh o2
  5577 7774 8446`. Steps are separate functions (`build_table`, `row_sets`, `control_lines`, `cross_validate`,
  `scores`); it replaces the table code that was repeated in the model notebook and in two experiment scripts. It
  reproduces the earlier table exactly (same row sets, night groups, OH groups, control ratios), with three
  differences: rows with ra = dec = 0 in any telescope pointing are excluded (32 corpus rows, 4 of the 4,838
  training rows), the corpus must match the decomposition coefficients (checked), and the ROTI columns are the
  corrected ones. Each run writes `<name>.pkl` and `<name>.json` (input files, corpus checksum, row counts, OH
  groups, test-night scores, package versions); the same record is on the model as `meta`, the OH groups as
  `oh_groups`. Models rebuilt in `skysub/oi6300_models/` (pickles git-ignored, JSON tracked), nested
  cross-validation by nights on 4,834 rows, ensemble of 5: OH, O2, 5577, 7774, 8446 (current) 0.0062 dex, 10.5% of
  rows with a dropped line; plus 5199 0.0056, 13.0%; no Science lines 0.0083. `OI6300Model`: `fit` rejects training
  rows without a valid target, boosting early stopping is switched off explicitly (it was already off below 10,000
  rows). `iono_corpus` CLI takes `--decomp-root` and `--stem`. Description of the current model:
  `skysub/mlp_predictor/OI6300_MODEL.md`. Tests: `sky_decomp/tests/test_oi6300.py` (7 pass). The earlier pickles in
  the experiment `outputs/` directories are kept as history. Not yet read by SkyML.

- `mlp_predictor.roti`: the broadcast navigation file is now the first one with at least 24 GPS satellites
  (`_get_nav`, `MIN_NAV_GPS`), not the first that exists. `BRDC00WRD_S` holds only 8-11 GPS satellites on 28 days
  (MJD 60181-60213), so ROTI on those days came from a third of the rays. The cache for those days was refetched
  (32 satellites; old files kept in `ionosphere_cache/roti/backup_incomplete_nav/`), the bundled ROTI table
  `sky_decomp/data/ionosphere/roti/roti_ctx_...npz` was rebuilt, and the `<arm>_roti_*` and `<arm>_dvtec_event*`
  columns of the ionosphere corpus were replaced on the 31 affected nights (546 rows; rows with a ROTI value
  66% -> 98%, median ray samples 4 -> 13; night 60284 also gained samples from station-days cached after the first
  build). No other column changed. The corpus column `log10_nmf2` is in log10 m^-3: the label was wrong (cm^-3) in
  `iono_corpus.COLUMNS`, the corpus README and the file's `columns_json`; values are unchanged. Models trained
  before this fix used the old ROTI values on those rows: the five `OI6300Model` pickles in
  `experiments_ionospheric_oi_model/outputs/` and the ionospheric-context MLP candidates were not retrained; the two
  Science-line ensembles in `experiments_external_data_v1/science_lines/outputs/` were. Test:
  `sky_decomp/tests/test_roti.py::test_navigation_download_skips_an_incomplete_file`.

- `mlp_predictor.oi6300`: optional Science-spectrum airglow lines as inputs of `OI6300Model` (`sci_lines`, default
  none, so existing models and saved pickles behave as before). Allowed by the user on 2026-11-27: the mesospheric
  control lines `'oh'` and `'o2'`, and `'5577'`, `'7774'`, `'8446'`; `'5199'` is an option to leave out where [N I]
  may be nebular. All enter as Science/SkyNear zenith log ratios (`control_ratio`, `science_lines`); [O I] 6300 itself
  is never an input. Lines a nebula can emit are guarded: the ratio of each is predicted from the sky telescopes and
  the control lines, `line_levels` returns the excess in units of its scatter on nights left out of the fit, and
  `predict` drops a line above `guard` (3) sigma for that exposure and uses the model trained without it (one model
  per subset of the guarded lines; `guard=None` trains the full model only). Test-night sigma on the 4,838 model C
  rows (`experiments_external_data_v1/science_lines/verify_model.py`): no Science lines 0.0084 dex; OH, O2 0.0074;
  plus 5577, 7774, 8446 0.0063 (0.0062 unguarded); plus 5199 0.0057. The guard drops a line on 11-13% of clean rows.
  With 1 dex of emission injected into one Science line the unguarded model is biased by up to 0.010 dex and the
  guarded one is not (`science_lines/REPORT.md`). Ensembles of 5 saved by
  `experiments_external_data_v1/science_lines/train_models.py` in `science_lines/outputs/` (`..._scilines_ens5.pkl`,
  `..._scilines_ni_ens5.pkl` with 5199); the earlier models are unchanged. Not yet read by SkyML. Tests:
  `sky_decomp/tests/test_oi6300.py` (5 pass).

- `mlp_predictor.oi6300` (new): `OI6300Model` predicts the Science [O I] 6300 (slant, R; or the `ATOM_Or`
  coefficient) from the two sky telescopes, log I_sci = log V_sci + a log z_near + (1 - a) log z_far + delta(x), with
  delta a gradient-boosting correction on the sky lines, ctx and ionosphere of all three telescopes and their
  geometry; no Science-spectrum input. `fit(..., members=5)` builds an ensemble of 5 models over random night groups
  (predictions averaged). Inputs are kept only if the Clouds training rows cover the non-Clouds Science
  range (>= 0.9); `predict` returns a quality flag. `clean_rows` (>= 10 fibres, log10 EW > 1.2, valid lines) and
  `equivalent_width` define the training cuts. Trained and validated in
  `experiments_ionospheric_oi_model/notebook_feature_histograms.ipynb` ("Final model") with cross-validation by
  nights for three training sets: LMC+SMC Science, gas-poor Science elsewhere, and both. Trained on both, test-night
  sigma is 0.0075 dex on LMC+SMC and 0.0086 dex on gas-poor Science elsewhere (SkyNear rescaled alone: 0.0129,
  0.0122). Optional extra `oi6300` (scikit-learn). Not yet read by SkyML. Test: `sky_decomp/tests/test_oi6300.py`.

- Ionosphere corpus (`mlp_predictor.iono_corpus`, new): one table with every ionospheric and thermospheric quantity per
  exposure and telescope for the knox 1.3.2 corpus, `sky_decomp/data/ionosphere/corpus/iono_corpus_...npz` (18,640
  rows, 146 columns, 8.7 MB): map and IGRF, IRI/NRLMSIS/HWM14 and model lines, bubble prior, ROTI, per-event dVTEC,
  observed F-region and 5577 coefficients, sky-arm altimeter and q, solar activity, exposure ids, training flag.
  `README.md` there lists each column (unit, source, spread); `diagnostics_cells/iono_corpus_hist.py` draws the
  histograms in `figs/`. Not read by SkyML.

- `mlp_predictor.roti.dvtec_event`: per-event electron deficit on the single ray nearest each pierce point
  (field-aligned, elevation >= 30 deg, +-5 min); `diagnostics_cells/dvtec_events.py` tests it against the F-region
  lines. Inside bubbles the sign and the line ordering are right (7774 slope +0.53 [+0.23, +0.88], 5577 control ~0)
  but the slope is ~1/4 of the physical value because the nearest ray is a median 116 km away, about a bubble's
  width; a correction of the Science prediction does not reduce held-out error. Not used by SkyML.

- ROTI as a context option: `mlp_predictor.roti.roti_at`, `build_table` (`python -m mlp_predictor.roti table`) and
  `bundled_table`; bundled per-exposure table `sky_decomp/data/ionosphere/roti/roti_ctx_...npz` (0.5 MB, 18,640 knox
  rows, row-aligned with the ionosphere table; field-aligned 200 km / 10 min search, rays at 92% of rows).
  `data._augment_triplet_with_roti` appends `iono_log10_roti_max`, `iono_log10_roti_med`, `iono_roti_avail`; inference
  builds them for ensembles that read them, and `IONCTX` becomes `map+thermosphere+roti`. `features_at` uses a sorted
  time window (same result, ~100x faster). Held-out test: no gain over the thermosphere variant (6300 +1.7%
  [-5.2, +9.7]), so no ensemble uses it; the default is unchanged.

- `mlp_predictor.thermosphere` (new): per-telescope F-layer and thermosphere context at the 285 km pierce point.
  IRI-2020 profile (PyIRI) scaled to the UQRG VTEC (hmF2, NmF2, bottom thickness, VTEC/IRI ratio); NRLMSIS 2.1
  (pymsis) O2, N2, O and T at hmF2; HWM14 magnetic-meridional wind and the vertical drift it drives; a Link & Cogger
  (1988) model 6300 brightness, its emission height and a 7774 proxy. New `thermosphere` extra (pymsis, PyIRI);
  HWM14 is compiled by `python -m mlp_predictor.thermosphere build-hwm` from a pinned, SHA256-checked upstream
  commit (needs gfortran). Bundled table `sky_decomp/data/ionosphere/thermosphere/thermo_ctx_...npz` (18,640 knox
  rows, row-aligned with the ionosphere table). `data._augment_triplet_with_thermosphere` appends six `iono_*`
  columns (hmF2, log O2, T, log model 6300, h6300, wind drift); inference builds them for ensembles that use them.
  Notebook variant `IONO_CTX=thermo` (`run_notebook.py --ionoctx-thermo`). Held-out clean MAD vs the ionospheric
  context variant: 6300 -6.3% [-12.9, -0.4], 7774 -3.5%, 8446 -3.4% (-32/-28/-33% vs deployed);
  `diagnostics_cells/thermosphere_signal.py` makes the report figures. `mlp_predictor.roti` now also keeps a
  60-min high-passed vertical TEC perturbation (dvtec) per sample (cache files `*_v2.npz`).

- SkyML runs ensembles that read the ionospheric context: `predict_ifu` builds the nine map and six thermosphere
  `iono_*` columns per telescope (bundled tables first, otherwise IONEX, PyIRI, NRLMSIS and HWM14) and stops with an
  explanatory error when they cannot be built, instead of substituting values. New primary keys `IONCTX`
  (none / map / map+thermosphere) and `IONXSRC`; the per-telescope values, `thermosphere_source` and model versions
  are in `SKYML_META.prediction_diagnostics.triplet`. Exposure 14771 with the thermosphere ensemble: F-line
  coefficients move by <= 0.007 dex; the computed path matches the bundled values to 9e-6. The isolated-branches
  ensemble reproduces 14771 bit for bit. The default ensemble is unchanged.

- Full ionospheric context variant: `data.THERMOSPHERE_EXTRA_NAMES` (NmF2, bottom thickness, VTEC/IRI, N2, O, model
  7774, meridional wind), `BUBBLE_PRIOR_NAMES` (`iono_bubble_p`, Swarm IBP reach probability) and `SKY_LINE_NAMES`
  (`iono_sky_altimeter`, `iono_sky_q`, from the SkyNear/SkyFar amplitudes only; Science carries SkyNear's) with
  their augments, built at inference for ensembles that read them. Notebook `IONO_CTX=full`
  (`run_notebook.py --ionoctx-full`): vs `thermo` 7774 -5.5% [-12.9, +1.8], 8446 -8.8% [-15.7, -1.1], 6300 +2.4%
  (in noise); vs deployed -30 / -32 / -38% (6300 / 7774 / 8446). Review fixes: a 404 for a RINEX day younger than
  4 days is no longer cached as missing; IONEX downloads use per-process temp files; non-finite map RMS raises;
  in-map IONEX EXPONENT records are honoured; incomplete downloads fall through to the next map source; HWM14
  downloads are written atomically. The thermosphere table is now tracked (`.gitignore`) and packaged.

- `mlp_predictor.ionosphere` (new, not yet wired into the context or SkyML): IONEX vertical TEC
  (UPC UQRG 15 min, CODE fallback, daily files cached in `skysub/ionosphere_cache/`) at the 285 km
  pierce point of each telescope (Science, SkyE, SkyW). Also RMS, a +-30 min VTEC trend, IGRF
  dip latitude and declination from the bundled `sky_decomp/data/geomag/igrf_285km_lco.npz`, and
  pierce-point offsets along and across the magnetic meridian. All TEC values are vertical: the
  van Rhijn factor stays applied only once, in `data.airglow_geometry_scale`.
  `python -m mlp_predictor.ionosphere --oi-table ... --out ...` builds the per-exposure table.
  `pyproject.toml` gains `ionosphere` and `geomag` (ppigrf, dev-only) extras.
  Later additions: d log10 VTEC/dt, local VTEC gradient at the pierce point, and sun altitude at the IGRF-traced
  magnetic conjugate point (+8.1 N, -74.5 W).

- `mlp_predictor`: `predict_sky_from_minimal_inputs` now feeds the context columns in the order of the
  ensemble's `ctx_names`, so new columns no longer shift an older model's inputs. The deployed ensemble's order is
  unchanged. `data._augment_triplet_with_ionosphere` adds nine `iono_*` columns per arm, built automatically at
  inference for ensembles trained with them. The new optional model and trainer key `iono_ctx_groups` routes them
  to the listed isolated branches only; it is off by default, and older checkpoints load unchanged.
  `pyproject.toml` gains a dev-only `notebook` extra.
  The bundled `sky_decomp/data/ionosphere/iono_ctx_lvmsframe_median_stack_1.3.2_knox_gaia1over100.npz` (2.9 MB) holds
  the features of 18,640 knox META rows; the augment reads it first, so training on that corpus needs no IONEX
  download. The isolated-branches notebook gains an `IONO_CTX` variant (`run_notebook.py --ionoctx-branches` or
  `--ionoctx-all`) and an `ionoctx-compare` cell. Its corpus root and the 14771 CFrame now follow
  `LVMSKY_DECOMP_DATA_ROOT` and `SAS_BASE_DIR`.

- `mlp_predictor.bubbles` (new, not wired into the context or SkyML): plasma-bubble reach probability per
  telescope. The ESA Swarm IBP climatology (`ibpmodel`, new `bubbles` extra) is calibrated in field-line apex height
  on 2,538 Swarm A/B in-situ bubble-index passes over the LCO sector on the 771 corpus nights (VirES HAPI, cached in
  `skysub/ionosphere_cache/swarm_ibi/`). Calibration: `sky_decomp/data/ionosphere/bubble_reach_calibration.json`.
  `diagnostics_cells/bubble_signal.py` is the go/no-go: no gain on the deployed residual for any line (all 95%
  night-bootstrap intervals contain zero), and Swarm-confirmed bubble exposures do not differ in the internal
  near/far proxy q. Kept as infrastructure for a local measurement (GNSS ROTI); not adopted as a feature.

- `mlp_predictor.roti` (new, not wired into SkyML): GNSS ROTI (5 min scatter of the 30 s geometry-free phase TEC
  rate, GPS L1C/L2W, elevation >= 30 deg) at 285 km pierce points from public BKG 30 s RINEX 3 and broadcast
  navigation (new `roti` extra, `hatanaka==2.8.1`); station-day results cached in `skysub/ionosphere_cache/roti/`,
  raw RINEX not kept. `diagnostics_cells/roti_signal.py` pilot on 201 nights (OAFA, SANT, CORD; 3,568 clean
  exposures): the near/far proxy q doubles when ROTI >= 1 TECU/min at a sky telescope (x2.0, 95% [1.2, 2.9],
  37 nights), but such exposures are 3.8% of the total and adding ROTI gives no cross-validated gain on the
  deployed residual. Candidate use: per-exposure bubble flag, not a model input.

- SkyML plasma-bubble QA: `pipeline.plasma_bubble_qa` computes GNSS ROTI near each telescope's 285 km pierce point
  (`mlp_predictor.roti.exposure_flag`, five IGS stations, RINEX fetched from BKG and cached) and records `BUBBLE`
  (1/0/-1 unknown), `ROTISCI`/`ROTISKE`/`ROTISKW`, the new warning bit `plasma_bubble` (`1 << 26`, appended to
  `sky_decomp.reliability`) and `SKYML_META.plasma_bubble`. Thresholds live in `Settings` (1 TECU/min, 200 km,
  10 min). QA only: never raises and never changes the sky; `--no-bubble-qa` / `bubble_qa=False` skips it.
  Exposure 14771 reproduces `ML_COEF`, `SKY`, `FLUX` and `MED_COEF` bit for bit with `BUBBLE = 0`.

- `mlp_predictor`: optional split of the `ionospheric` group into one group per
  emission mechanism (`data.IONOSPHERIC_SPLIT`, `data.split_ionospheric_group()`),
  a per-group `head_dim` mapping in the model, and a `group_weight_overrides`
  trainer key. Off by default; the deployed ensemble loads and predicts as before.
  Used by `experiments_ionospheric_split_heads_v1` (not adopted).

- `mlp_predictor`: optional isolated airglow branches (`isolated_branches`,
  `isolated_branch_inputs` model/trainer keys): small per-group MLPs that bypass the
  trunk, read the full sci/near/far context plus chosen groups' near/far scores, and
  add to the group's alpha blend. Also `data.split_atomic_group()` moving [O I] 5577
  into its own `oi_green` group (95 km). Off by default; the deployed ensemble loads
  and predicts as before. Used by `experiments_ionospheric_branches_v1`.

- `mlp_predictor`: optional `fixed_blend_alpha` model/trainer key holding a group's
  near/far blend weight fixed instead of learned. Off by default.

- Extend the least-contaminated-fiber median fallback to SkyE/SkyW, targeting
  10 eligible fibers per Sky telescope (Science keeps 50). Persist fallback headers, selection counts
  and a prediction reliability warning. Retain finite `InsufficientProgress`
  median coefficients with their original status and a propagated fit-failure
  bit; invalid inputs and nonfinite coefficients remain fatal.
- Queue block-CLI path lists without scanning CFrame headers or contents first.
  Workers resolve each input and output path independently, with batch progress
  active during preparation. Missing/invalid inputs follow per-worker policies.
- Remove CFrame content hashing from ordinary block calculations and use
  size/mtime metadata for CLI reuse fingerprints. Legacy preparation reuse keeps
  its hash check when no metadata identity exists. Old CLI fingerprints require
  explicit recomputation; metadata identity does not detect same-size/mtime edits.
- Set a three-second default socket timeout when importing `sky_decomp`, covering
  decomposition and Sky Subtraction workers and existing-product validation
  before DRP/Astroquery imports.
- Batch decomposition and Sky Subtraction share CSV LSF loading, physical-fiber
  selection averaging, and fixed-Blue validation. Blue is held fixed; R/Z start
  from the same selected-fiber CSV means and are refined on the target spectrum.
- Reference mode is the batch default and uses the bundled `sky_decomp/data/lsf`
  calibration. CSV median inputs require `FIBER_ROLE`; explicit solar fitting
  and explicitly selected legacy FITS calibration directories remain supported.
- Reference provenance records the CSV manifest SHA256, averaging convention,
  membership use, and per-channel contributor counts in LSF metadata. Previous
  compact caches are invalidated by source and reference fingerprints.

### Validation

- 128 block-CLI, scheduling, block-LSF and compact-product tests passed in the
  Python 3.11 uv runtime, including deferred subprocess workers with successful,
  missing and malformed FITS inputs. Scientific fitting was stubbed in the new
  scheduling integration test; no full Utah exposure was recomputed.
- 35 targeted tests passed in the Python 3.11 uv runtime: shared reference and
  worker/API/compact contracts (27), plus numerical fixed-Blue/R/Z refinement,
  profile-mean, strong-continuum and SkyML compact checks (8).
- No new full-corpus calculation or neural-network training was performed.
