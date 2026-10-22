# Sky-only block LSF subtraction

This method shares the median decomposition, ML predictor, sky-arm residual transfer and Science-median OH correction with `sky_subtraction`. Its LSF stage uses selected SkyE and SkyW fibers only, separately for each spectrograph in the R and Z channels, where sky lines constrain the profile. B is not fitted: every Science fiber keeps the survey LSF of the preparation there (it still enters the model, the flux and the residuals). Science spectra never enter this LSF fit.

## Command line

Use the existing scientific Python environment and install the small CLI dependencies:

```sh
conda activate lvmdrp311
python -m pip install -r skysub/sky_subtraction_block/requirements-cli.txt
python -m skysub.sky_subtraction_block.cli --help
```

Worker bootstrap requires PyTorch and `threadpoolctl`. Install PyTorch in the
same Python environment used to launch the CLI, then install all packages from
`requirements-cli.txt` above. Before executing a batch, the CLI checks that
`torch` and `threadpoolctl` can be imported and reports failures once, before
workers or run records are created. Help and dry runs do not import them.

Run the module from the repository root using the active environment's Python.
Root, `single`, `batch` and legacy help put examples immediately after the usage
and short description, before detailed option descriptions.

```sh
# The explicit one-exposure command; --ensemble remains a compatible alias.
python -m skysub.sky_subtraction_block.cli single /data/lvmCFrame-00011017.fits \
  --ml-model-file /models/trained-sky.pt --nproc 1 --cpus 1 \
  --output /results/lvmSFrame-skyml-block-00011017.fits

# Discover exposures through the standard SAS mirror and its DRP catalogue.
python -m skysub.sky_subtraction_block.cli batch \
  --redux-dir /data/sas/sdsswork/lvm/spectro/redux/1.3.2 \
  --expnum 11017,11018 --ml-model-file /models/trained-sky.pt \
  --gaia-dir /cache/gaia/fibers --nproc 1 --cpus 2 \
  --output-dir /results --preserve-drp-tree

# Inspect the plan before any fitting, network queries or output creation.
python -m skysub.sky_subtraction_block.cli batch \
  --mjd 60325:60330 --max-exposures 3 --ml-model-file /models/trained-sky.pt --dry-run
```

`--cpus C` sets the total logical-CPU budget. `--nproc N` sets workers per
exposure (default: 1 in batch, 8 in single/legacy). Batch automatically runs up to
`floor(C / N)` exposures simultaneously, with disjoint N-CPU allocations.
For example, `--cpus 16 --nproc 1` runs up to 16 exposures; `--cpus 16 --nproc 2`
runs up to 8 exposures with two workers each. Require N <= C. When C is omitted,
batch uses N CPUs and runs one exposure at a time. Serial stages use one worker;
native numerical threads remain one per process. Memory grows with concurrent
exposures. `--jobs` remains a compatible alias for `--nproc`.
Use `--cpu-list 8-11 --cpus 4` to select exact scheduler-allowed logical CPUs.

Use `--redux-dir DIR` for the input directory containing the catalogue and
`TILEGROUP/TILE/MJD` tree. The directory name is unrestricted. The CLI discovers
exactly one `drpall-*.fits` there, excluding `*-incomplete.fits`, and derives
the version from that filename. Zero or multiple catalogues are errors; use
`--drpall FILE` to select a complete catalogue explicitly. Catalogue rows and
CFrame headers must agree with the discovered version. There is no
`--drp-version` option. `--cframe-root` remains an alias for `--redux-dir`.

Without an explicit input directory, the CLI uses `LVM_SPECTRO_REDUX` or the
standard redux path below `--sas-base-dir` / `SAS_BASE_DIR`. It searches that
directory, then its immediate version directories if no catalogue is present
directly; ambiguous selections require an explicit directory/catalogue.

Preserving the tree is the default and does not prepend the DRP version.
Saved-state reuse follows the same `TILEGROUP/TILE/MJD` layout under
`--reuse-state-dir`. Existing outputs/snapshots from the older version-prefixed
layout are not moved automatically.

For a file of exposure IDs, use the same input-list option directly:

~~~sh
python -m skysub.sky_subtraction_block.cli batch \
  --redux-dir /data/sas/sdsswork/lvm/spectro/redux/1.3.2 \
  --input-list exposures.txt --ml-model-file /models/trained-sky.pt \
  --gaia-dir /cache/gaia/fibers --output-dir /results --cpus 4 --nproc 1
~~~

An ID list accepts one integer per line or whitespace/comma-separated integers.
A path list accepts one CFrame path per line, with relative paths resolved
against the list file. Both formats support blank lines and whole-line
`#` comments; do not mix IDs and paths in one file. ID lists use catalogue
selection and may combine with MJD/tile/exposure-time/calibration filters,
but not `--expnum` or `--all-exposures`. Path lists are explicit inputs and
cannot combine with catalogue filters.

Interactive batches show one live panel with processed/total, successful,
failed, running and queued exposures, elapsed wall time, successful exposures
per minute, approximate ETA and active exposure IDs. Reused/skipped inputs are
counted separately from new successes. ETA uses successful completions per
wall-clock time; failures do not increase that rate, and ETA is unknown until
the first success. This is exposure-level progress, not progress within a fit.

Completion/failure events print above the panel without overwriting it.
Routine start events appear only at DEBUG level. Full paths, PIDs and CPU
allocations are written to `OUTPUT_DIR/.batch/<run ID>/batch.log`; individual
worker logs stay beside it. The run-record directory is printed once at
startup. `--no-progress` and noninteractive runs retain compact event lines
without a live panel. Dry runs and DEBUG retain the full per-exposure plan.

The trained model is selected by `--ml-model-file FILE`; `--ensemble FILE`
remains accepted by all three entry points. Ordinary operation computes
intermediate state in memory and writes the final SFrame. `--keep-preparation`
explicitly enables an additional snapshot for later `--reuse-state FILE`
(compatible alias `--seed FILE`) or batch `--reuse-state-dir DIR` reruns.

`--gaia-dir` is the direct parent of eight-digit files such as
`/cache/gaia/fibers/lvmGAIA-fibers-00011017.fits`. Missing/identity-mismatched
cache triggers the existing ARI TAP query; malformed FITS/schema is an error.

Preserved output placement is
`OUTPUT_DIR/TILEGROUP/TILE/MJD/lvmSFrame-skyml-block-EXPNUM.fits`.
`--flat-output` writes directly in OUTPUT_DIR (default: `./skyml-block-output`).
For unavailable explicit inputs, `--missing-input skip` records the path; unknown
header identities use `-1` and the output path stays unset. `--on-existing skip` verifies
scientific fingerprints, complete FITS layout and checksums before reuse.
Run plans, results, status, exact invocations and individual logs are stored
under `OUTPUT_DIR/.batch/<run ID>`, or the explicit `--run-dir`.

Full flag definitions, directory trees, saved-state rules, selector examples
and CPU combinations are in [the CLI reference](CLI_PROPOSAL.md).
Local real-data commands and product checks are stored under
`analysis/cli_validation/` in the validation workspace.

## Method

1. Select hardware-good, Gaia-clean Sky fibers (Gaia / observed G < 0.01). Prepare masked telescope medians, decompose them, predict the Science component coefficients, and determine the OH corrections.
2. Decompose each selected Sky fiber with its telescope-median LSF and O2 temperature fixed. Each fiber gets its own nonnegative component coefficients. Freeze these coefficients for the block fit.
3. Solve one native-pixel weighted least-squares QP for each of six spectrograph × (R, Z) blocks. The LSF density is `sum H[k,q,s] M[k](offset) B[q](wavelength) S[s](fiberid)`. The initial slit model uses four cubic B-spline controls over the full spectrograph fiber-ID range. The wavelength basis is the existing four-control cubic basis. SkyE and SkyW use their own component coefficients and telescope transmission arrays.
4. Enforce `H >= 0` and `sum_k H[k,q,s] = 1` for every wavelength/slit control pair. The centroid is unconstrained by default; `--centered` (`centered=True`) also enforces `sum_k mu[k] H[k,q,s] = 0`, with `mu[k]` the exact first moment of the normalized cubic M-spline. Partition of unity propagates the equalities to every wavelength and slit position. There is no wavelength-solution shift or added shape/width constraint. Native IVAR supplies weights; input/protected-line masks apply. Five Angstroms at channel edges are trimmed from the fit because the affine designs hold neighboring-channel leakage fixed.
5. Evaluate the fitted R and Z surfaces at each Science fiber, keep that fiber's survey B LSF, and render the common ML coefficients with this LSF. Apply the sky-arm residual transfer and Science-median OH correction to every fiber, through the shared renderer. The transfer follows the current `residual_*` settings whatever the preparation recorded; the block product records them and the science-line velocities in `SKYML_META`. Render each usable SkyE/SkyW fiber as its telescope median model with the block R/Z surface at that fiber (B: telescope median LSF). Store the actual final SKY and FLUX = observed − SKY for Science and Sky fibers.

```sh
python -m skysub.sky_subtraction_block /path/lvmCFrame-00059750.fits \
  --ml-model-file /path/trusted_ensemble.pt --gaia-dir /path/gaia/fibers \
  --output /separate/products/lvmSFrame-skyml-block-00059750.fits --nproc 8
```

For a controlled LSF comparison, add `--seed /path/lvmSFrame-skyml-00059750.fits` to reuse that exact exposure's saved median/ML/correction preparation. The input SHA256 and row/wavelength identities must match. The seed's individual Science LSF fits are replaced completely. A run without a seed computes the same median/ML/correction stage (survey-only, no per-fiber LSF fits) in memory and writes no preparation file; `--keep-preparation` also writes it beside the output for reuse with `--seed`. A `--seed` file, or a preparation that already exists at that default path, is used as is and never modified or deleted. Existing destinations are protected; working filenames have no version suffixes.

## Basis dimensions and regularization

`--offset-count` (default 11, odd and ≥5), `--wave-count` (default 4, ≥4),
and `--slit-count` (default 4, ≥4) set the R/Z cubic tensor dimensions.
The offset support stays ±2.5 Å; B retains its exact original survey profile.
Sky-fiber component coefficients are fitted with the original median profiles;
the R/Z profiles are then projected onto the trial basis to construct the block
designs. Rendering and reconstruction use each channel's actual dimensions.

`--offset-roughness`, `--wave-roughness`, and `--slit-roughness` add squared
second-difference penalties on the corresponding coefficient axis.
`--reference-prior` adds squared distance to the Science-median profile projected
onto the requested basis, repeated across slit controls. Each strength multiplies
the mean diagonal of the native-IVAR data Hessian. All four defaults are zero.
The reported block χ² remains the data-only objective, excluding these penalties.
The older `roughness_fraction`, `offset_roughness_fraction`,
`fallback_prior_fraction`, and numerical ridge belong to other LSF fitters;
they were not regularizers of the block QP.

For example, run an independent single-core trial with:

```sh
taskset -c 0 python -m skysub.sky_subtraction_block /path/lvmCFrame-00014350.fits \
  --ml-model-file /path/trusted_ensemble.pt --gaia-dir /path/gaia/fibers \
  --nproc 1 --offset-count 11 --wave-count 4 --slit-count 4 \
  --offset-roughness 0.01 --wave-roughness 0 --slit-roughness 0 --reference-prior 0 \
  --output /separate/trial/lvmSFrame-skyml-block-00014350.fits
```

## Compact storage and diagnostics

The shared [compact product contract](../sky_subtraction/FITS_STRUCTURE.md) applies. `BLOCK_LSF_COEF` stores the six blocks × offset × wavelength × slit controls; `BLOCK_SLIT_KNOTS`, `BLOCK_WAVE_KNOTS` and `BLOCK_LSF_META` store exact domains, degrees, solver status, rank, input counts, chi-square, and moment checks. `SKY_FIBER_COEF` records every selected Sky fiber's identity, component vector and decomposition status. The usual `LSF_COEF/KNOTS/META` contain evaluations of the R/Z global solution at Science fibers together with their survey B LSF, so existing `SkyMLReconstructor.sky()` works directly. No intermediate model images are persisted. `SKYML_META` records preparation/input hashes, settings and implementation hashes. `SkyBlockReconstructor.global_state()` and `.sky_fiber()` reconstruct block profiles and selected Sky-fiber fits in memory.

## Sky coverage and explicit reference fallback

Channel-valid Sky fibers determine each R/Z block independently. An unavailable
individual Sky-component fit is recorded in `SKY_FIBER_COEF` with `usable=False`,
its reason and NaN coefficients. It does not enter the surface fit.
The coverage checks require at least twice the slit-basis count in distinct
positions, full slit rank, condition number ≤1000, a selected-fiber convex hull
covering the full physical spectrograph domain, and a maximum gap no larger than
twice the maximum physical Sky-fiber spacing in the spectrograph. The gap scale
comes from all Sky positions in SLITMAP, before photometric/hardware selection.

An unsupported or numerically unavailable block adopts the projected
Science-median reference LSF, constant across pseudo-slit. Fallback applies to
the whole block, avoiding a seam between fitted and reference intervals. This
is a reference profile, not a surface measured from the unavailable data.
The reference retains its centroid, including in a requested centered run.
The compact coefficient shapes remain the requested dimensions and B is unchanged.
`BLOCK_LSF_META` records `adopted`, `fallback_reason`, requested centering,
channel input/rejection counts, slit condition and gap diagnostics.
`LSF_QA` marks fitted-block acceptance separately from reference adoption.
Every stored fallback is reconstructed by the same block evaluator.

These are conservative support checks, not physical-profile validation.
Input IVAR describes observed noise only; coefficient/model covariance is not
propagated. A numerically feasible block and a lower subtraction χ² do not
establish physical LSF recovery. The 15 × 16 × 6 χ² minimum found on exposures
14350, 14771 and 14414 is withdrawn because of excessive flexibility and severe
channel-edge artifacts. No replacement is recommended from that ranking alone.
The fallback addresses the block stage: median/ML preparation still needs its
own usable input data.

## Conditional likelihood and tensor diagnostics

`BLOCK_LSF_META.data_chi2/data_n_pixels` evaluate the auxiliary Sky-component
model on channel-fit pixels for either a fitted block or the fixed reference
fallback. They exclude regularization penalties. `chi2` remains NaN for an
unfitted reference block. `nominal_df` is `(n_offset-1-centered)*n_wave*n_slit`
for a fitted block and zero for a fixed reference. `effective_df` is the local
active-face smoother trace on the equality tangent, using the data and penalized
Hessians. Active-face selection and shared median/component fitting uncertainty
are not included; this is not the effective dimension of the complete pipeline.
The positivity threshold and calculation time are also recorded.

The fixed-eleven analysis reports conditional Gaussian AIC/BIC with nominal and
effective dimensions as diagnostics, plus delivered Science/Sky residual scores.
The data likelihood is conditional on the preparation and individual component
coefficients. The analysis records full CLI wall/CPU time independently of its
later profile checks. `hessian_bytes`, `penalized_bytes`, `coefficient_bytes` and
`sample_design_bytes` describe named fitter arrays, not process peak memory.

Validation: `python -m pytest skysub/sky_subtraction_block/tests -q`. The injected test recovers a varying global surface from multiple fibers and checks continuous positivity, normalization and centering between control nodes.
