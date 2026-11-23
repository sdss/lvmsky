# Experiment changelog and retention review

Last reviewed: 2026-08-20

## Purpose and scope

This document is the historical index for experimental work under `skysub`.
It records when an experiment was made, what scientific or numerical question
it addressed, what the latest local evidence says, and whether its generated
data should be retained.

The primary inventory contains all 60 top-level directories named
`experiment_*`, `experiments_*`, or `sky_model_search_framework*` present on
2026-08-20. The moonless-continuum notebook/output series is included as one
additional experiment because its products live under the shared `output/`
directory rather than an experiment directory. Loose legacy notebooks are
listed separately when they do not have enough provenance for a scientific
status decision.

Most experiment directories are currently untracked by Git. Dates below are
therefore taken, in descending order of authority, from dated run IDs and
manifests, experiment logs, or the oldest/newest file timestamps. A date range
means that the local artifact tree was developed or updated over that range;
it is not a Git commit range. Only
`experiments_moon_zodi_baseline_decomp_jax_v1` and part of
`experiments_moon_zodi_lsf_surface_iterative_v1` have tracked Git history in
the current checkout.

## Status and retention vocabulary

- **Retain -- accepted/current**: an accepted result, a current dependency, or
  the best available reference dataset. Do not delete its source or canonical
  artifacts.
- **Retain -- diagnostic**: scientifically unaccepted but still active or
  uniquely informative. Keep until its successor is reviewed.
- **Archive compact evidence**: keep source, README/notebook, manifests,
  summaries, metrics, and representative figures; large reproducible arrays,
  checkpoints, and caches may then be removed.
- **Superseded / closed**: no longer a candidate for continuation. Preserve the
  compact historical record; generated data are deletion candidates unless a
  live dependency is documented.
- **Undocumented**: no authoritative result statement was found. Do not delete
  scientific data solely on directory age; first write a result summary or
  confirm that the run was abandoned.

Optimizer convergence, numerical integrity, or lower pooled loss is not
scientific acceptance. The status below follows the experiment's declared
local, per-channel, held-out, and physical gates where they exist.

## Executive review

The current experiment tree occupies about 37.2 GiB, excluding the shared
`output/` tree (about 1.31 GiB). The review supports three immediate decisions:

1. Keep all source code and compact result evidence for history.
2. Keep current dependency datasets, especially the independent-line OH
   reference, corrected-decoder corpus, accepted positive full-sky products,
   selected Moon/Zodi inputs, and the accepted/current LSF configuration.
3. Remove only explicitly listed reproducible or superseded binary products
   after a dry-run inventory. No deletion was performed during this review.

Approximately 10.3 GiB of generated data is listed below as immediate or
conditional cleanup candidates. About 4.8 GiB is low-risk reproducible cache,
environment, or superseded-input cleanup; the remaining approximately 5.5 GiB
requires compact-evidence archiving first. Totals are rounded independently.

## General, infrastructure, and early experiments

| Date | Experiment | Purpose | Result and current status | Data decision |
|---|---|---|---|---|
| 2026-01-29--02-14 | `experiment_gaia_templates` | Compare full-spectrum Gaia-band flux with Gaia TAP sources per LVM fiber. | Notebook prototype; no scientific result or acceptance record found. **Undocumented.** | Keep source/notebook; local Gaia passband tables are small. |
| 2026-01-30 | `experiments_gaia` | Store Gaia source/fiber cross-match material for exposure 9387. | Two data files only; no model or result summary. **Undocumented.** | Keep until the Gaia-template notebook is closed or reproduced. |
| 2026-03-04--03-12 | `sky_model_search_framework2` | Baseline and grouped Stage-1 sky-model search with alternative QP solvers. | Historical framework retained as a reference by its README; superseded architecturally by framework 3 and later decomposition experiments. **Superseded / closed.** | Archive README and run manifests; generated runs may be removed after confirming no notebook reads them. |
| 2026-03-04--03-12 | `sky_model_search_framework3` | Sparse regional Stage-1 emission/continuum model with shared line amplitudes and local continua. | Later framework prototype, but no final result summary; current production work follows other paths. **Undocumented historical prototype.** | Keep source/tests; archive or remove demo runs only after recording solver outcomes. |
| 2026-06-24--06-25 | `experiments_sky_mixing_v1` | Compare robust SkyE/SkyW mixing against DRP sky subtraction on exposures 24986 and 19943. | Executed diagnostic with FITS/HTML outputs; no promotion or cross-exposure validation. **Archive compact evidence.** | Keep experiment log, notebook, metrics, and representative HTML; binary model/output copies are deletion candidates. |
| 2026-08-06 | `experiments_jax_cpu_mps_benchmark_v1` | Establish JAX CPU/MPS backend and dtype capability. | CPU64, CPU32, and MPS32 gates passed; MPS64 was a capability probe. Later end-to-end solver benchmarks supersede this environment test. **Superseded benchmark.** | Keep lock/results; delete the reproducible `.venv` (about 595 MiB). |
| 2026-08-06--08-07 | `experiments_jax_native_lbfgs_benchmark_v1` | Test device-resident L-BFGS/L-BFGS-B on the full joint model. | L-BFGS-B float32 exact-bound diagnostics failed because negative subspace boundary ratios reversed a descent direction; no 1,000-step production candidate was authorized. **Closed diagnostic.** | Archive selection and exact-bound summaries/notebook; large repeated run/timing products are conditional deletion candidates. |
| 2026-08-07 | `experiments_jax_reduced_lbfgs_benchmark_v1` | Compare reduced exact-inner JAXopt, Optax, and Optimistix solvers on CPU32 and MPS32. | Preflight and physical gates passed; Optax was selected, and MPS was about 8x faster end-to-end. Both 1,000-update endpoints were still improving and remain diagnostic, not final sky models. **Retain -- numerical reference.** | Keep selected long-run summaries and reproducibility artifacts; repeated timing replicas may be compacted. |

## Continuum, diffuse, and full-sky experiments

| Date | Experiment | Purpose | Result and current status | Data decision |
|---|---|---|---|---|
| 2026-07-20 | `experiments_blue_continuum_joint_v2` | Fit 3650--5550 A Moon, Zodi, diffuse, and learned signed continuum-transfer modes. | Stored row-wise 5% gate passed, but later grouped review found train/holdout pointing leakage. Empirical signed modes also weaken physical interpretation. **Superseded by v3.** | Keep compact v2 comparison evidence; generated basis/notebook outputs may be archived. |
| 2026-07-20--08-04 | `experiments_blue_continuum_positive_v3` | Replace signed modes with a constructively nonnegative full-grid continuum model. | Positive by construction, but independent Moon-lit worst-bin values of 6.47% and 6.39% fail a universal 5% requirement. **Retain -- latest blue-continuum diagnostic.** | Keep source, deterministic 72-spectrum asset, notebook, and gate evidence. |
| 2026-07-20--08-11 | `experiments_positive_sky_decomposition_v1` | Build positive continuum states and a full-stack LSF-aware line decomposition. | Five-state-per-family continuum reached 0.969% validation and 0.898% test median error on 20 A summaries; all 14,026 full-sample group solves passed and the 20 A full-sample median was 0.923%. **Retain -- accepted/current reference.** | Do not delete canonical selection, accepted state basis, merged line products, or validation evidence. |
| 2026-08-02--08-06 | `notebook_experiments_moonless_far_continuum_variability.ipynb` and `output/moonless_far_*` | Explore Moon-free continuum, native-grid PCA/nearly-NMF bases, and OH-blend residual regularization. | Native v5 preserves all 12,401 pixels. The N0/N1/N2 blend-null pilot failed its predeclared skyline and matched-line gates, so no mass retraining was authorized. Earlier v2/v3/non-native-v4 products are superseded. **Retain native v5; close failed pilot.** | Keep `output/moonless_far_free_bspline_o2_anchor_native_v5`; delete v2, v3, and `o2_anchor_v4` after compact evidence is retained (about 138 MiB). |
| 2026-08-05 | `experiments_joint_diffuse_einstein_jax_v1` | Jointly learn diffuse-vector and within-group OH corrections on 10 spectra. | Training WRMS improved strongly, but the run is in-sample and has active-bound/degeneracy warnings. Superseded by the 50-spectrum v2 design. **Archive compact evidence.** | Keep final summary, selected endpoint, diagnostics, and parameter maps; prune intermediate checkpoints/model cubes. |
| 2026-08-05--08-06 | `experiments_joint_diffuse_einstein_jax_v2` | Extend joint diffuse/OH fitting to 50 spectra and test backend precision. | Long fit improved its initialization but final mean WRMS 0.09746 was 37.2% worse than the comparison fit and won 0/50; status is REJECT/diagnostic-only. Backend runs also stopped at the iteration cap. **Superseded / closed scientifically.** | Preserve export provenance, final summaries, and notebook; `outputs/full_50spec` is a conditional deletion candidate (about 1.04 GiB). |
| 2026-08-12--08-13 | `experiments_diffuse_refinement_v1` | Test alternating QP, B-spline, adaptive-weight, and joint-native Adam refinements for diffuse components. | The tree contains many `completed_diagnostic_only` runs and one explicit `completed_diagnostic_rejected` static-gate run, but no experiment-level README or selected final candidate. **Undocumented active/diagnostic collection.** | Do not delete `dark_fresh_decomposition_v1` or a possible selected run yet. First add a branch-level decision; then remove redundant smoke, 100/1,000-step, and dominated schedule runs. |

## Moon, zodiacal-light, OH, and LSF experiment chain

| Date | Experiment | Purpose | Result and current status | Data decision |
|---|---|---|---|---|
| 2026-08-11 | `experiments_moon_differential_pair_study_v1` | Build a native-grid, same-night Moon differential-pair and bright-anchor corpus with editable visual review. | Corpus build completed; `reviewed_corpus_v1` preserves the role-scoped manual decisions, and `visual_audit_v2` is the current readable audit. This is a selection/reference dataset, not a fitted Moon model. **Retain -- current corpus/reference.** | Keep manifests, reviewed corpus, visual decisions, summary, and v2 audit; the superseded v1 visual gallery is a deletion candidate. |
| 2026-08-05--08-11 | `experiments_moon_zodi_jax_v1` | Earliest JAX Moon/Zodi single/demo fit. | Prototype superseded by corpus and reviewed-clean branches; no authoritative result summary. **Superseded / undocumented.** | Keep minimal source/notebook; generated demo products may be removed. |
| 2026-07-17--08-04 | `experiments_moon_zodi_lmfit_v1` | Establish the detailed physical Moon and zodiacal-light forward model and packaged assets. | Foundational physical implementation reused by later blue-continuum work; not itself a production sky-subtraction model. **Retain -- dependency/reference.** | Keep source and checksum-locked physical assets. |
| 2026-08-05 | `experiments_moon_zodi_corpus_jax_v1` | First corpus-level analytic Moon/Zodi fit. | Explicitly identified as rejected by v2. **Superseded / closed.** | Keep source history; image gallery and derived products are deletion candidates. |
| 2026-08-05--08-06 | `experiments_moon_zodi_corpus_jax_v2` | Analytic corpus model on the strict bright-Moon cohort. | Median error reached 4.98%, but exposure p90 was 18.67%; failed the 5% contract. **Closed diagnostic.** | Archive coefficients, metrics, provenance, and representative plots; bulk gallery data may be removed. |
| 2026-08-07 | `experiments_moon_zodi_palace_anchored_jax_v2` | Add metadata-driven Moon/Zodi carriers and PALACE diffuse anchors. | Pilot branch superseded by diffuse-physical v3. **Archive compact evidence.** | Preserve summaries/manifests; duplicated pilot arrays are deletion candidates. |
| 2026-08-07 | `experiments_moon_zodi_diffuse_physical_jax_v3` | Fit global season/time/F10.7/van-Rhijn diffuse amplitudes with Moon and Zodi. | Executed pilot, but no experiment-level acceptance statement was found. Later reviewed-clean branches became the primary Moon/Zodi path. **Undocumented diagnostic.** | Keep source and compact pilot evidence; do not delete the only result before adding a decision. |
| 2026-08-07 | `experiments_moon_zodi_anchor_same_night_jax_v4` | Combine bright absolute anchors with same-night differential constraints. | Converged numerically but failed bright-anchor and all-band differential gates; status `diagnostic_only`. **Retain -- dependency.** | Do not delete its approximately 2.31 GiB run yet: v5 reads this exact run as an input. |
| 2026-08-07--08-15 | `experiments_moon_zodi_anchor_pairs_jax_v5` | Curate low/high-horizon pairs and refit anchor plus differential roles. | Combined solve stalled above the gradient threshold; B held-out pairs regressed and the 5% anchor gate failed. **Closed diagnostic, but downstream provenance reference.** | Keep source, manifest, summary, and selected pair corpus; compact duplicate figures/arrays only after decoupling from v4. |
| 2026-08-06--08-07 | `experiments_moon_zodi_baseline_decomp_jax_v1` | Fit shared Moon/Zodi parameters to baseline line-subtracted native spectra. | `fit_results.json` reports `failed_5pct_gate`. **Superseded / closed.** | Archive result/provenance; derived decomposition parts may be removed if no successor imports them. |
| 2026-08-11 | `experiments_moon_zodi_reviewed_clean_jax_v6` | Fit a reviewed-clean physical Moon/Zodi model to native-grid anchors and pairs. | P1 rejected; status `diagnostic_only`. **Superseded / closed.** | Keep compact report and parameters; generated fit arrays are deletion candidates. |
| 2026-08-11--08-17 | `experiments_moon_zodi_reviewed_clean_jax_v8` | Validate a nine-parameter geometric model, high-resolution carriers, and native-grid prepared inputs. | Numerical preflight passed; the 15-update fit was explicitly preliminary. Later analytic and differential branches import v8 code/data contracts. **Retain -- dependency/reference.** | Keep source, authoritative v5 inputs, and joint-input cache; delete superseded v3/v4 input caches and `decomposition_lsf_surface_n5_v2` (about 1.91 GiB). |
| 2026-08-12 | `experiments_moon_zodi_reviewed_clean_fast_analytic_jax_v9` | Establish a fast 31-parameter analytic Moon/Zodi pilot. | Fifteen-step diagnostic, not converged. Superseded by v10--v13. **Superseded.** | Keep source and compact metrics only. |
| 2026-08-12 | `experiments_moon_zodi_reviewed_clean_fast_analytic_cv_v10` | Run five-fold whole-MJD CV for the fast analytic model. | `moon_linear_zodi_quadratic31` was selected over v11, but the family still missed scientific accuracy gates and was finally closed by v13. **Retain -- historical selected baseline.** | Keep runtime, coefficients, and OOF summary; remove duplicate per-fold intermediates if reproducible. |
| 2026-08-12 | `experiments_moon_zodi_reviewed_clean_fast_analytic_cv_v11` | Change only the loss to training-only q05-relative normalization. | Improved pairs/equal-role metric but worsened anchors; explicitly not retained. **Superseded / closed.** | Archive notebook and gate table; bulk fold outputs may be removed. |
| 2026-08-12 | `experiments_moon_zodi_reviewed_clean_fast_analytic_cv_v12` | Continue the 31-parameter model to deterministic 3,000-update multistart convergence. | Stable nominal analytic fit used as the final frozen model in the v13 review. **Retain -- final analytic baseline.** | Keep contract, selected checkpoint/parameters, and fitted outputs needed by v13. |
| 2026-08-12 | `experiments_moon_zodi_reviewed_clean_fast_analytic_micro_v13` | Screen additional physics/geometry terms and close the analytic branch. | Final negative result: held fold-0 anchor/pair errors 16.33%/39.07%; no tested term transferred without regression or sign instability. **Retain -- final closure report.** | Keep final notebook, v12 nominal result, ledgers, and selected diagnostics; dominated screen cubes may be compacted. |
| 2026-08-12 | `experiments_moon_zodi_baseline_decomp_differential_jax_v1` | Apply A3 to absolute anchors, signed pairs, and non-anchor bright endpoints. | `summary.json` reports `diagnostic_failed_relative_gate`. Superseded by v2. **Closed diagnostic.** | Keep source/report/summary; `outputs/run_v1` and duplicate basis products are conditional cleanup candidates. |
| 2026-08-12 | `experiments_moon_zodi_baseline_decomp_differential_jax_v2` | Add an exact-audited detector-ratio approximation and physical contrast terms. | `summary.json` again reports `diagnostic_failed_relative_gate`. **Closed diagnostic.** | Archive compact evidence; large exact-projection/prediction arrays may be removed. |
| 2026-08-12 | `experiments_moon_zodi_pysr_v1` | Search a compact positive grey geometry correction. | Discovery-only protocol; no promoted expression. Superseded by v2. **Closed diagnostic.** | Keep search summary/expression audit; derived search caches may be removed. |
| 2026-08-12 | `experiments_moon_zodi_pysr_v2` | Search wavelength-dependent symbolic Rayleigh/aerosol corrections. | Best candidate improved micro loss but only 3/5 pairs and had 12.155% worst-band regression; status `diagnostic_failed_micro_gate`. **Closed.** | Preserve rejected-expression record and report; search/capacity caches are deletion candidates. |
| 2026-08-07 | `experiments_moon_zodi_noll_lsf_jax_v1` | Fit native-grid Noll Zodi and analytic Moon carriers with profiled diffuse terms. | Deployable Moon+Zodi failed all-band validation; profiled diffuse absorbed the missing continuum. Full-corpus/test run was not opened. **Closed diagnostic.** | README explicitly identifies bases/decomposition as reproducible intermediates; remove them after retaining splits, parameters, metrics, and figures (about 1.37 GiB). |
| 2026-08-15--08-17 | `experiments_moon_zodi_lsf_surface_iterative_v1` | Assess the opt-in iterative Moon/Zodi LSF-surface decomposition on a frozen holdout. | Immutable assessment runner; baseline remains default regardless of outcome. Current production code/tests depend on this contract. **Retain -- current validation path.** | Keep source, executed report, manifests, and final holdout evidence. |
| 2026-08-18 | `experiments_moon_zodi_lsf_spline_grid_v1` | Search R/Z B-spline capacity and knot placement while B remains constant. | 1,037 fits completed; validation retained information-adaptive cubic R6/Z6. Higher-capacity and arc-gap grids failed fixed bright-line gates. **Retain -- accepted/current configuration.** | Keep complete selection/gate evidence and winner/current holdout products. Rejected candidate spectra may be archived separately. |
| 2026-08-13--08-14 | `experiments_moon_zodi_full_sky_joint_adam_300_v1` | Jointly fit Moon/Zodi geometry plus all local line/diffuse amplitudes with fixed LSF. | All long trajectories remain diagnostic; balanced runs failed continuum/decomposition non-regression and the 30,000-update continuation did not establish scientific acceptance. **Superseded by later flexible-LSF/OH work.** | Archive manifests, final summaries, selected endpoints, and representative diagnostics; the 1.88 GiB output tree is a conditional deletion candidate. |
| 2026-08-18 | `experiments_moon_zodi_oh_lsf_joint_adam_v1` | Joint complete Moon/Zodi/diffuse/OH/non-OH model with shared OH weights and flexible LSF. | Aggregate loss improved, but B-channel and [O I] 5577 non-regression failed; status `diagnostic_failed_channel_partition_non_regression`. **Closed diagnostic.** | Preserve final gate/ablation summaries and notebook; approximately 877 MiB of run data are a conditional deletion candidate. |
| 2026-08-18--08-19 | `experiments_moon_zodi_oh_lsf_fluxnorm_adam_v2` | Add physical OH grouping, flux-normalized within-group corrections, and protected non-OH blends. | Objective fell 91.43%, but B core worsened 6.765x and [O I] 5577 SSE rose 118.26x in all 100 spectra. **Scientifically rejected.** | Keep README, final summary, ablation, provenance, and notebook; approximately 608 MiB of run/checkpoint data are a conditional deletion candidate. |
| 2026-08-19 | `experiments_moon_zodi_oh_lsf_fluxnorm_adam_v3` | Repeat flux-normalized OH learning from canonical unsuffixed PALACE with six-basis LSF. | Numerical gates passed, but B core/wing ratios reached 45.87/6.15 and [O I] 5577 SSE worsened 2,158.9x. **Scientifically rejected.** | Preserve compact audit and notebook; approximately 570 MiB of run data are a conditional deletion candidate. |
| 2026-08-19 | `experiments_moon_zodi_oh_hierarchical_lsf_chisq_adam_v1` | Use hierarchical PALACE-centered OH corrections, pseudo-chi-square, fixed B LSF, and flexible R/Z LSF. | The 10,000-step snapshot passed execution and local channel/line gates, including improved [O I] 5577, but extreme corrections and lack of MJD-disjoint validation prevent physical acceptance. **Retain -- active diagnostic candidate.** | Keep the corrected v2 trajectory and complete audit; the rejected free-B v1 safety run can be compacted after archiving. |
| 2026-08-19 | `experiments_moon_zodi_oh_palace_phase_lsf_adam_v1` | Test phase/PCHIP PALACE hierarchy with a less flexible B-channel LSF. | The production-linear trajectory reduced total loss, but [O I] 5577 RMS worsened 4.43x. Later B/R/Z six-basis work supersedes it. **Superseded / closed.** | Archive final diagnostic/physical-audit summaries; approximately 456 MiB of outputs are a deletion candidate. |
| 2026-08-19 | `experiments_moon_zodi_oh_palace_brz_lsf_adam_v2` | Fit PALACE OH hierarchy with six-basis surfaces in all B/R/Z channels. | 10,000-iteration in-sample candidate passed configured critical windows, including [O I] 5577 ratio 0.898, but no held-out linear-decomposition gate has passed. **Retain -- current diagnostic candidate.** | Keep full candidate until held-out assessment and export decision are complete. |
| 2026-08-19 | `experiments_moon_zodi_oh_palace_brz_lsf_transfer_adam_v3` | Transfer v2 OH state to disjoint exposure ranks and refit the target corpus. | Exposure IDs are disjoint but 50 MJDs overlap. Final/warm [O I] 5577 median RMS is 1.019 and only 44% improve; OH 9387 also slightly regresses. **Retain -- diagnostic transfer evidence, not acceptance.** | Keep while judging transferability; compact only the 100-step smoke run after final evidence is frozen. |
| 2026-08-19--08-20 | `experiments_moon_zodi_oh_flat_family_adam_v1` | Tie unresolved OH transitions into flat observable families and monitor checkpoint-local gates. | Iteration 1,000 passed all local gates; iteration 2,000 reduced pooled loss further but worsened [O I] 5577 by 50.49%. Status `diagnostic_failed_checkpoint_local_non_regression`; no 10,000-step continuation recommended. **Retain iteration 1,000 diagnostic; close continuation.** | Keep the iteration-1,000 checkpoint/evidence. Archive or delete the rejected iteration-2,000 state after preserving its gate summary. |

## OH grouping, decoder, and reference experiments

| Date | Experiment | Purpose | Result and current status | Data decision |
|---|---|---|---|---|
| 2026-06-25 | `experiments_oh_variability_v1` | Explore OH variability, ridge line-weight corrections, and a two-spectrum free-line pilot. | Early diagnostic with no formal accepted export; later independent-line and grouping experiments supersede its scientific role. Its 12-spectrum sample is still used by the decoder A/B test. **Archive/reference dependency.** | Keep selection and compact metrics; large per-spectrum tables/FITS may be archived after decoder A/B is self-contained. |
| 2026-06-26--08-05 | `experiments_oh_data_driven_v1` | Develop data-driven OH grouping/weights and full-model JAX prototypes. | Current README identifies `random_weights_vn217_noadd_fixedlsf_10spec_s250` as the main diagnostic run. It is an experimental reference, not a promoted model. Its local JAX environment is imported by later joint experiments. **Retain -- dependency/reference.** | Keep source, current run, input FITS, and `.venv_jax` while imports remain; older smoke runs can be pruned after a dependency audit. |
| 2026-07-27--08-05 | `experiments_oh_upper_state_v1` | Produce the nearly unregularized per-exposure independent-line OH reference. | Complete: 14,026 successful spectra, 11,393 independent lines, no failed spectra. This is the high-flexibility reference used by grouping work. **Retain -- canonical reference.** | Do not delete dataset manifests, sparse measurements, line/spectrum tables, or authoritative dense checkpoints without an external archived copy. |
| 2026-07-29--07-30 | `experiments_oh_group_search_v1` | Search model atoms plus promoted singleton lines and open lockboxes once. | `production_v1` K=3625 failed condition lockbox; pre-existing K=4500 `production_v2` passed the temporal lockbox and is the accepted artifact. **Retain -- accepted grouping.** | Keep accepted line table, full decision record, and lockbox evidence; rejected development candidates may be compacted. |
| 2026-07-30--07-31 | `experiments_oh_group_search_v2` | Compress the independent-line reference with Ward blends and later recurrent searches. | Current frozen H11 K=100 candidate passes global spectral gates but fails strict local-feature historical gates. The older `frozen_candidate` is explicitly superseded H8 evidence. **Retain current search; reject production promotion.** | Keep H11 baseline/current candidate/final reports. Delete the superseded `frozen_candidate` and dominated search scratch after preserving manifests. |
| 2026-07-30 | `experiments_oh_membership_jax_v1` | Test differentiable hard OH membership plus shared relative weights. | Feasibility run improved 55/77 condition spectra but used the condition set and cannot be tuned further; temporal lockbox remains reserved. **Archive compact evidence.** | Keep selected matrix, condition result, split/provenance, and summary; training-state duplicates are deletion candidates. |
| 2026-07-30--08-03 | `experiments_oh_einstein_jax_v1` | Learn 7,174 relative Einstein weights while profiling group amplitudes. | Fold-4 catastrophes fell from 2,072 to 277, but the result still failed catastrophe gates, was worse than B6 K350, and produced implausible factors up to 4,679. **Closed diagnostic.** | Preserve selected checkpoint, held-out comparison, and physical audit; bulky profiles/checkpoints are deletion candidates. |
| 2026-07-31 | `experiments_oh_group_comparison_v1` | Compare physical JAX, learned-membership JAX, and F1 groupings on one frozen condition sample. | Descriptive comparison only; F1 had lower OH WRMS in the selected low/high-lunar subset, but the temporal lockbox was untouched. **Retain -- compact comparison evidence.** | Keep complete small result tree. |
| 2026-08-15 | `experiments_oh_relative_weights_joint_adam_v1` | Learn shared positive within-group OH corrections jointly with local amplitudes. | 10,000 updates improved spectra but produced non-identifiable factors up to 4,015,810; physical result rejected. Warm restarts also failed their optimizer criterion. The earlier 100-step cache is still an input to later experiments. **Mixed: retain dependency, close long run.** | Keep `corpus300_adam100_v1`; remove rejected `corpus300_adam10000_v1` and warm-restart outputs after compact evidence (about 491 MiB). |
| 2026-08-15--08-17 | `experiments_oh_lsf_flexibility_closure_v1` | Test six-, twelve-, and twenty-four-basis R/Z LSF flexibility with OH/continuum fixed. | Predeclared diagnostic experiment; no authoritative completed result statement was found in the README. **Undocumented/incomplete.** | Keep source and current result tree until a final assessment is written. |
| 2026-08-18 | `experiments_oh_decoder_ab_v1` | Isolate the effect of corrected HITRAN upper-state decoding on 12 spectra. | Completed A/B diagnostic; superseded as a corpus input by the 300-spectrum corrected-decoder run. **Archive compact evidence.** | Keep summaries and representative viewer; duplicate full decompositions may be removed. |
| 2026-08-19 | `experiments_oh_decoder_corpus_v1` | Recompute the frozen 300-spectrum corpus with corrected decoding and canonical unsuffixed PALACE `Aij * gi`. | Complete corpus with status `diagnostic_only`; it is the direct input to current PALACE B/R/Z and transfer experiments. **Retain -- current dependency/reference.** | Do not delete selection, cache, catalogs, parts, manifest, or QC until all dependent experiments are self-contained. |
| 2026-08-13 | `experiments_palace_suffix_comparison_v1` | Compare canonical and suffixed PALACE tables on a common five-spectrum fit. | Compact contract/provenance comparison; still relevant to `palace_oh_suffix` selection history. **Retain -- contract evidence.** | Keep source and v3 summary/viewer; older v1/v2 viewers may be removed if v3 is self-contained. |

## ML transfer experiments

| Date | Experiment | Purpose | Result and current status | Data decision |
|---|---|---|---|---|
| 2026-11-19--11-20 | `experiments_oi_dedicated_mlp_v1` | Test small dedicated MLPs for the O I lines (6300/6364, 7774+8446, 5577, and all together) against the deployed ensemble on the 1.3.2 knox corpus, with rows whose science [O I] differs from the sky arms by more than 0.3 dex removed. | Held-out nights (n = 1,361): [O I] 6300/6364 MAD 0.0104 -> 0.0080 dex (-22.7%, night-bootstrap 95% [-28.4%, -13.6%]); similar gains for 7774, 8446 and 5577. [O I] 6300 left after subtraction (133 every10 test exposures, +/-3 A): median 2.60% -> 1.80% of the line flux (1.06% if fitted on the science spectrum). Moon/zodi reversal gate applied on the full decomp files (removes 0 rows). Exposure 14771 through SkyML with the O I values replaced (`exposure.py`): [O I] 6300 left in the Science fibers -2.9% (deployed) vs +0.2%; one exposure from a training night. Not measured on sky-subtracted CFrames; not wired into SkyML. **Retain -- diagnostic.** | Keep `notebook_oi_dedicated_mlp.ipynb`, `outputs/run_v1/results.json`, `outputs/run_v1/oi_models.pt`, `outputs/exposure_14771/`; `outputs/oi_table.npz` and `predictions.npz` are reproducible from `prepare.py` / `train.py`. |
| 2026-11-20 | `experiments_ionospheric_split_heads_v1` | Replace the deployed model's single `ionospheric` head (N I 5199, [O I] 6300/6364, O I 7774/8446) by three heads, one per emission mechanism (widths 96/64/32, same total 192), and retrain the 10-seed ensemble with Niv's notebook and config, plus a retrained unsplit control. | Held-out test, clean rows (n = 1,361), MAD dex, control -> split: [O I] 6300 0.0103 -> 0.0102, 7774 0.0210 -> 0.0211, 8446 0.0207 -> 0.0196, [N I] 0.0153 -> 0.0154; every night-bootstrap interval contains zero, and control vs deployed differs by as much. Dedicated [O I] MLPs remain 14-21% better. Exposure 14771: [O I] 6300 left in Science fibers -3.0% (split) vs -2.9% (deployed). **Not adopted.** | Keep `notebook_ionospheric_split_heads.ipynb`, `outputs/control_run.ipynb`; ensembles and the 14771 product in `outputs/` are reproducible with `run_notebook.py`. |
| 2026-11-20 | `experiments_ionospheric_branches_v1` | Move [O I] 6300/6364, O I 7774+8446, [N I] 5199 and [O I] 5577 off the shared trunk into four isolated branches (one per mechanism, 136 inputs: full sci/near/far context plus the five lines' near/far scores; hidden 48/24, 64/32, 32/16, 48/24), each added to its own scalar near/far alpha blend, as the parent does for zodi. 31,973 parameters vs 31,684 in the parent ionospheric head; per-line loss weights unchanged; Niv's notebook and config, 10 seeds. | Held-out test, clean rows (n = 1,361), MAD dex, control -> branches: [O I] 6300 0.0103 -> 0.0084 (-19% [-27, -14]), 7774 0.0210 -> 0.0179 (-15%), 8446 0.0207 -> 0.0166 (-20%), 5577 0.0049 -> 0.0039 (-20%), [N I] unchanged; night-bootstrap intervals exclude zero. Dedicated [O I] MLPs no longer better (-4% to +1%, intervals contain zero). Moon, zodi, continuum, OH unchanged. Exposure 14771 residual in Science fibers, deployed -> branches: 5577 -1.8% -> 0.0%, 6300 -2.9% -> -1.9%, 7774 -3.7% -> -2.4%, 8446 -3.6% -> -2.3%. An alpha = 1 variant (near arm only for the three F-region groups) is equal or slightly better on held-out rows (7774 -6%, 8446 -9%, 6300 -4% within noise) but over-subtracts [O I] 6300 on exposure 14771 (-4.4% vs -1.9%). **Promising; adoption pending review.** | Keep `notebook_ionospheric_branches.ipynb`; the ensemble and the 14771 product in `outputs/` are reproducible with `run_notebook.py`. |

## Loose legacy notebooks without experiment-level provenance

The following top-level notebooks are experimental in character but do not
currently have an experiment README, immutable manifest, or declared result
gate. They are not assigned a scientific acceptance status by this review:

- `manual_experiment_short.ipynb` (2026-03-17)
- `notebook_experiments_bspline_lsf.ipynb` (2026-08-07)
- `notebook_experiments_lsf_surface_iterative_brightness.ipynb` (2026-08-11)
- `notebook_nd_decomp_experiment.ipynb` (2026-08-11)

Before deleting any associated data, either move the notebook into the
relevant experiment directory and add a short result record, or mark it as an
abandoned prototype with its input/output paths.

## Proposed data deletion plan

No files have been deleted. Paths below are proposals based on the current
dependency search and local disk usage.

### Tier A -- low-risk reproducible or superseded data

These items can be removed after one final path-existence dry run and after
confirming that the compact evidence named above is present:

| Approx. size | Path | Reason |
|---:|---|---|
| 595 MiB | `experiments_jax_cpu_mps_benchmark_v1/.venv` | Reproducible environment; lock and benchmark results remain. |
| 106 MiB | experiment `__pycache__`, `.ruff_cache`, and `.ipynb_checkpoints` directories | Interpreter/editor caches, not scientific results. |
| 360 MiB | `experiments_moon_zodi_reviewed_clean_jax_v8/outputs/spherical_model_inputs_q56_drp_nine_parameter_v3` | Superseded by v5. |
| 720 MiB | `experiments_moon_zodi_reviewed_clean_jax_v8/outputs/spherical_model_inputs_q56_drp_nine_parameter_v4` | Superseded by v5. |
| 878 MiB | `experiments_moon_zodi_reviewed_clean_jax_v8/outputs/decomposition_lsf_surface_n5_v2` | Superseded by the joint-input v3 decomposition cache. |
| 745 MiB | `experiments_moon_zodi_noll_lsf_jax_v1/outputs/bases` | README labels this a reproducible intermediate; scientific gate failed. |
| 657 MiB | `experiments_moon_zodi_noll_lsf_jax_v1/outputs/decomposition` | README labels this a reproducible intermediate; retain splits/metrics/fit/figures. |
| 234 MiB | `experiments_moon_differential_pair_study_v1/outputs/visual_audit_v1` | Superseded by `visual_audit_v2`; retain the reviewed corpus and editable decisions. |
| 138 MiB | `output/moonless_far_free_bspline_v2`, `output/moonless_far_free_bspline_v3`, `output/moonless_far_free_bspline_o2_anchor_v4` | Superseded by the native-grid v5 product. |
| 491 MiB | `experiments_oh_relative_weights_joint_adam_v1/outputs/corpus300_adam10000_v1` and `.../corpus300_adam_warm_restart_9900_v1` | Rejected/non-identifiable long run and failed warm restarts; the separately retained 100-step cache is the live dependency. |

Tier A is approximately 4.8 GiB including caches. The relative-weight cleanup
accounts for approximately 0.48 GiB of that total.

### Tier B -- archive compact evidence first, then remove bulk data

| Approx. size | Path | Required retained evidence |
|---:|---|---|
| 456 MiB | `experiments_moon_zodi_oh_palace_phase_lsf_adam_v1/outputs` | Final diagnostic summary, physical audit, trajectory fingerprint, and short result note. |
| 877 MiB | `experiments_moon_zodi_oh_lsf_joint_adam_v1/outputs` | Final run/channel/ablation/extreme-line summaries and executed notebook. |
| 608 MiB | `experiments_moon_zodi_oh_lsf_fluxnorm_adam_v2/outputs` | 10,000-step summary, checkpoint hash, channel/[O I] gate evidence, provenance, and notebook. |
| 570 MiB | `experiments_moon_zodi_oh_lsf_fluxnorm_adam_v3/outputs` | Final diagnostic, physical OH audit, provenance, and notebook. |
| 1.88 GiB | `experiments_moon_zodi_full_sky_joint_adam_300_v1/outputs` | Input selection manifest, final summaries for each objective family, selected 30,000-step endpoint summary, and representative diagnostics. |
| 158 MiB | `experiments_joint_diffuse_einstein_jax_v1/outputs/full_10spec` | Final endpoint/summary and parameter/degeneracy diagnostics. |
| 1.04 GiB | `experiments_joint_diffuse_einstein_jax_v2/outputs/full_50spec` | Final summary, comparison metrics, export manifest, and notebook. |

Tier B is approximately 5.5 GiB. The full Tier A plus Tier B proposal is about
10.3 GiB before filesystem block and compression differences.

### Explicitly protected from deletion for now

- `experiments_oh_upper_state_v1/outputs/all_sky_far_h11_v1`: canonical
  independent-line reference.
- `experiments_positive_sky_decomposition_v1` accepted selection, continuum,
  and merged line products.
- `experiments_oh_decoder_corpus_v1/outputs/corpus300_canonical_aij_gi_v1`:
  current corrected-decoder input to PALACE B/R/Z experiments.
- `experiments_oh_relative_weights_joint_adam_v1/outputs/corpus300_adam100_v1`:
  imported by later OH/LSF experiments despite its parent experiment's rejected
  long-run result.
- `experiments_moon_zodi_anchor_same_night_jax_v4/outputs/anchor_sequence_20260807_001`:
  read directly by v5.
- `experiments_moon_zodi_reviewed_clean_jax_v8/outputs/decomposition_lsf_surface_n5_joint_inputs_v3`
  and `.../spherical_model_inputs_q56_drp_nine_parameter_v5`: current reviewed
  inputs/dependency contracts.
- `experiments_moon_zodi_lsf_spline_grid_v1/outputs/run_v1`: accepted R6/Z6
  selection evidence.
- `experiments_moon_zodi_oh_hierarchical_lsf_chisq_adam_v1` corrected v2 run,
  `experiments_moon_zodi_oh_palace_brz_lsf_adam_v2`, and
  `experiments_moon_zodi_oh_palace_brz_lsf_transfer_adam_v3`: current diagnostic
  candidates still under interpretation.
- `output/moonless_far_free_bspline_o2_anchor_native_v5`: native-grid continuum
  and basis provenance, including the rejected pilot audit.
- `experiments_oh_data_driven_v1/.venv_jax`: still imported by later joint JAX
  experiments; remove only after those runtimes are decoupled.

## Required archive contract before deleting a scientific run

For every Tier B run, preserve at minimum:

1. experiment source and tests;
2. README and executed canonical notebook;
3. input selection and all source/table hashes;
4. exact wavelength shape, dtype, units, target, preprocessing, masks, and
   normalization provenance;
5. optimizer configuration, final parameter/checkpoint identity, and runtime;
6. final scientific status and every acceptance-gate metric;
7. a compact final parameter table or vector when it is needed to reproduce
   the reported model;
8. representative observed/model/residual/component and LSF diagnostics;
9. a machine-readable list of removed paths and their sizes.

The native scientific artifacts must continue to state element-for-element
wavelength-grid equality, 12,401 unchanged wavelength columns where applicable,
and zero removed or replaced samples. A transformed or diagnostic product must
not inherit a `native` or `full-grid` identity after its source arrays are
removed.

## Recommended next maintenance actions

1. Add a short README/result decision to
   `experiments_diffuse_refinement_v1`,
   `experiments_oh_lsf_flexibility_closure_v1`,
   `experiments_moon_zodi_diffuse_physical_jax_v3`, and
   `sky_model_search_framework3` before any scientific cleanup there.
2. Make v5 anchor-pair inputs self-contained so the failed 2.31 GiB v4 run can
   later be reduced to compact evidence.
3. Add a small machine-readable experiment registry in a future change only
   after this human-reviewed changelog is accepted; do not infer scientific
   status automatically from optimizer output.
4. Execute Tier A cleanup first, record exact reclaimed bytes, then review the
   Tier B archive bundles before approving deletion.
