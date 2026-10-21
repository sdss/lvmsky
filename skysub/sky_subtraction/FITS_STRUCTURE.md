# Compact SkyML SFrame contract

Current contract, measured on 2026-10-15. Working filenames remain `lvmSFrame-skyml-<exposure>.fits`; `SKYMLVER=3` is internal schema provenance only.

## Size

| Exposure | Current MiB | Previous full-diagnostic MiB | Reduction |
|---|---:|---:|---:|
| 3469 | 445.02 | 1659.13 | 73.2% |
| 12630 | 434.82 | 1657.33 | 73.8% |
| 16819 | 454.17 | 1660.73 | 72.7% |
| 44808 | 438.51 | 1658.00 | 73.6% |
| 59750 | 442.46 | 1658.68 | 73.3% |

The standard frame content occupies about 304 MiB; final SKY adds 92 MiB. The rest stores coefficient/state/selection/context tables and observed median data. The standard LSF FWHM image is retained; full LSF surface coefficients are separate and compact.

## Layout shared by all methods

Every method (in-fiber, block) writes the same extensions in the same order, with the same column names: the lvmSFrame block `FLUX, IVAR, MASK, WAVE, LSF, SKY, FLUXCAL_STD, FLUXCAL_SCI, FLUXCAL_MOD, SLITMAP`, then the SkyML extensions in the order of the table below (`product.SHARED_ORDER`), then the method-specific ones, then `SKYML_META`. Headers, WAVE, FLUXCAL_*, and SLITMAP (with units, plus `skyml_applied`) are inherited from the lvmCFrame. `FILENAME` names the output file. Method-specific primary keywords are only `SKYMETH` (and `LSFCENT`, `BLOCKCH` for block). Deliberately absent against the DRP lvmSFrame: `SKY_IVAR` (see Storage boundary); against the lvmCFrame: `SKY_EAST/WEST` and their IVAR.

What differs by method is content, not structure:

| Extension | In-fiber | Block |
|---|---|---|
| `LSF_COEF/KNOTS/META` | Adopted per-fiber R/Z profile (in situ or survey) | R/Z block surface at the fiber; B survey |
| `CAND_LSF_*` | Last fitted candidate, adopted or not | Repeats `LSF_*` (every block profile is adopted) |
| `LSF_QA` | R/Z record per Science fiber from the fit | R/Z record per Science fiber with the same columns; fit-only columns are NaN/0, `reason` is `sky_block` or `sky_block_centered`, `converged` is the QP status |
| `FIT_LSF_HISTORY` | Present when per-fiber fits ran | Absent (no per-fiber fit); also absent for `--survey-only` |
| `BLOCK_LSF_*`, `SKY_FIBER_COEF` | Absent | Block surfaces and Sky-fiber coefficients |

### Sky fibers

`SKY` and `FLUX` cover the Science fibers and every usable (`fibstatus == 0`) SkyE and SkyW fiber; `SLITMAP.skyml_applied` marks these rows and `telescope` tells them apart. A Sky fiber's SKY is its telescope's median decomposition (`MED_COEF`) rendered with that fiber's LSF: the telescope median LSF (`MED_LSF_*`) for the in-fiber method; the block R/Z surface at the fiber plus the median B LSF for the block method. There are no ML, OH or Blue residual-transfer terms. `FLUX` is the observed minus that SKY. It is a residual of the fibers that build the median (and, for block, the LSF), so it is not an independent test. `LSF` holds the FWHM of the LSF used. Sky fibers are not in `LSF_FIBERS`, `LSF_COEF` or `LSF_QA`; `SkyMLReconstructor.sky(fiberid)` renders them in memory. Standard-star (`Spec`) and flagged fibers keep the observed CFrame flux.

## Stored extensions

Example: exposure 59750, 1944 input fibers × 12401 wavelengths; 1754 processed Science fibers. Shapes are NumPy dimensions. Sizes include padded FITS data blocks, excluding headers.

| Extension | Type / shape | MiB | Meaning |
|---|---|---:|---|
| `PRIMARY` | header; — | 0.000 | Exposure header; SKYMFMT=COEFFICIENTS; internal format SKYMLVER=3. |
| `FLUX` | f4; 1944 × 12401 | 91.964 | Delivered observed minus applied SKY on processed Science and Sky rows. |
| `IVAR` | f4; 1944 × 12401 | 91.964 | Observed-noise inverse variance; zero on invalid adopted-sky pixels. |
| `MASK` | uint8; 1944 × 12401 | 22.992 | Input flags plus invalid adopted-sky flags. |
| `WAVE` | f4; 12401 | 0.049 | Shared air wavelength (Angstrom). |
| `LSF` | f4; 1944 × 12401 | 91.964 | Standard per-pixel FWHM image (Angstrom), updated from adopted profiles. |
| `FLUXCAL_STD` | table; 12401 | 1.609 | Input standard-star calibration table. |
| `FLUXCAL_SCI` | table; 12401 | 1.609 | Input science calibration table. |
| `FLUXCAL_MOD` | table; 12401 | 1.609 | Input model calibration table. |
| `SLITMAP` | table; 1944 | 0.283 | Fiber identity/geometry/telescope/status and skyml_applied. |
| `SKY` | f4; 1944 × 12401 | 91.964 | Final sky actually subtracted (Science and Sky rows); the only stored full-resolution sky-model image. |
| `GAIA` | table; 1944 | 0.368 | Fiber photometry, Gaia contamination ratios, eligibility/selection/reasons. |
| `GAIA_SOURCES` | table; optional | — | Gaia source matches; present only when photometry came from ARI TAP rather than an aggregate cache. |
| `MED_LABELS` | table; 3 | 0.003 | Median spectrum/state order: Sci, SkyE, SkyW. |
| `MED_META` | table; 1 | 0.003 | Observing geometry, conditions, selection thresholds and fiber counts. |
| `MED_QA` | table; 3 | 0.003 | Median fit status, diagnostic objective and O2 temperature. |
| `MED_FLUX` | f8; 3 × 12401 | 0.286 | Actual observed masked median spectra. |
| `MED_IVAR` | f8; 3 × 12401 | 0.286 | Input median-error inverse variance, independent Gaussian approximation. |
| `MED_COUNT` | i8; 3 × 12401 | 0.286 | Contributors per median wavelength. |
| `LSF_FIBERS` | table; 1754 | 0.027 | Compact state index → input_row and physical fiberid. |
| `LSF_COEF` | f8; 1754 × 3 × 11 × 4 | 1.769 | Adopted/applied full profile coefficients: fiber × B/R/Z × offset basis × wavelength basis. |
| `LSF_KNOTS` | f8; 1754 × 3 × 8 | 0.321 | Adopted wavelength-knot vectors. |
| `LSF_META` | table; 5262 | 2.991 | Adopted LSF representation, dimensions, bounds, offset geometry and metadata. |
| `CAND_LSF_COEF` | f8; 1754 × 3 × 11 × 4 | 1.769 | Last fitted candidate coefficients, including candidates not adopted. |
| `CAND_LSF_KNOTS` | f8; 1754 × 3 × 8 | 0.321 | Candidate wavelength knots. |
| `CAND_LSF_META` | table; 5262 | 2.958 | Candidate representation metadata. |
| `MED_LSF_COEF` | f8; 3 × 3 × 11 × 4 | 0.005 | LSF coefficients from each median decomposition. |
| `MED_LSF_KNOTS` | f8; 3 × 3 × 8 | 0.003 | Median-fit wavelength knots. |
| `MED_LSF_META` | table; 9 | 0.008 | Median-fit representation metadata. |
| `LSF_QA` | table; 3508 | 0.711 | R/Z (Science fibers) convergence/adoption/profile guards, observed line-information pass, counts/median contrast/median SNR in four quarters. |
| `MED_FIT_MASK` | uint8; 12401 | 0.014 | Shared protected-line mask used in the median decomposition. |
| `MED_COEF` | f8; 3 × 388 | 0.011 | Median decomposition coefficients, using internal 1e14 flux scaling. |
| `COMPONENT_NAMES` | table; 7 | 0.003 | oh, moon, zodi, diffuse, atom, orc, o2. |
| `ML_COEF` | table; 388 | 0.014 | Named common IFU coef, ensemble coef_std and OH-corrected oh_coef; internal scaling. |
| `FIT_LSF_HISTORY` | table; 9735 | 33.434 | Per-iteration LSF, amplitudes, source/variance-continuum spline parameters and Huber scale; reconstructs spectral fits and robust weights. |
| `RENDER_META` | uint8; 1136 | 0.003 | UTF-8 JSON constructor settings, LSF configuration, per-median PWV/airmass/O2 temperature and telescope order. |
| `RENDER_TRANSMISSION` | f8; 3 × 12401 | 0.286 | Exact three telluric transmission vectors used by the median model constructors. |
| `LSF_FIT_MASK` | uint8; 12401 | 0.014 | Shared 1D protected-line exclusion mask used in individual R/Z fits. |
| `INIT_LSF_COEF` | f8; 3 × 3 × 11 × 4 | 0.005 | Initial selected-fiber mean LSFs, needed for exact auxiliary tensor/edge contributions. |
| `INIT_LSF_KNOTS` | f8; 3 × 3 × 8 | 0.003 | Initial mean-state wavelength knots. |
| `INIT_LSF_META` | table; 9 | 0.005 | Initial mean-state metadata. |
| `SKYML_META` | uint8; 301961 | 0.288 | UTF-8 JSON provenance: input/ensemble paths and SHA256, LSF `calibration` path and manifest, Gaia mode (and cache-miss reason), code SHA256 and git revision, all `Settings` operating points, corrections and timing. |

## Reconstruction API

Block basis searches can use distinct offset dimensions in B and R/Z.
`LSF_META.n_offset` records the actual offset count for each channel;
`LSF_COEF` is NaN-padded to the largest channel count. The reader accepts older
products without this column using the full stored offset axis. The existing
`n_basis` and `n_knots` columns likewise delimit wavelength coefficients and knots.
`SKYML_META.block_settings` records all three basis counts and the dimensionless
offset/wavelength/slit roughness and median-profile prior strengths, with their
Hessian scaling. Median and initial states retain their original bases.

```python
from skysub.sky_subtraction.reconstruction import SkyMLReconstructor

with SkyMLReconstructor("/path/to/lvmSFrame-skyml-00059750.fits") as result:
    sky = result.sky(680)
    diagnostic = result.fiber(680)
    iteration_5 = result.fiber(680, iteration=5)
```

- `sky()` renders the full final vector from recorded coefficients, LSF states, observation context and matching model/telluric assets. It needs neither the CFrame nor the trained ML archive (the `mlp_predictor` package is still imported for the residual transfer). The adopted rendering is checked against persisted SKY; incompatible rendering assets raise an error. Fit weights and the residual transfer use the operating points recorded in `SKYML_META.settings`.
- `fiber()` additionally reads the original CFrame (path and SHA256 recorded in SKYML_META) for observed spectra, native errors and masks. Supply `cframe=...` to the constructor after moving that input. It performs no fitting or ML inference.
- `diagnostic_hdul()` supplies report generators with an in-memory view of reconstructed arrays. Those arrays are never written back into the compact file.

### Reconstructed quantities

```text
SKY_ML   = common IFU ML coefficients rendered with adopted per-fiber LSF
SKY_OH   = OH-corrected rendering minus uncorrected OH rendering
SKY_RES  = sky-arm residual transfer rebuilt from observed sky medians and their decompositions
SKY      = SKY_ML + SKY_OH + SKY_RES
FLUX     = original observed spectrum - SKY

FIT_LINES = saved local line amplitude × candidate/iteration LSF line tensor
FIT_SOLAR = saved local solar amplitude × solar tensor + fixed diffuse/edge term
FIT_CONT  = saved nonnegative spline coefficients evaluated on saved knots
FIT_SKY   = FIT_LINES + FIT_SOLAR
FIT_MODEL = FIT_SKY + FIT_CONT
```

SKY_ML is the model **before corrections**. SKY_RAW is not a stored extension or a renderer API name. SKY_EPI is reconstructed from diagonal ML coefficient ensemble variance and is not total sky uncertainty. SKY_SURVEY can be rendered for comparison via `sky(fiberid, survey=True)`; it loads the survey LSF from the calibration directory recorded in SKYML_META and refuses a directory whose manifest differs.

The auxiliary fit uses saved iteration LSF/amplitudes/spline; the delivered sky uses adopted LSF and common IFU coefficients. Source spline and local amplitudes are never added to the final sky subtraction.

Robust fit weights are reconstructed from native CFrame IVAR, the stored pre-LSF variance-continuum spline, the iteration residual and the saved Huber scale. This preserves the original 3% continuum-dependent variance floor. These weights are not input IVAR; plotted normalized residuals use native flux errors.

## Sky-only block LSF variant

`skysub.sky_subtraction_block` writes `lvmSFrame-skyml-block-<exposure>.fits` with the shared layout above; individual Science fits are replaced by block surfaces. `LSF_COEF/KNOTS/META` hold the B survey LSF and the R/Z block surfaces evaluated at each Science fiber, so `SkyMLReconstructor.sky()` applies unchanged. `LSF_QA` is the in-fiber table (R/Z rows only). The following extensions are added before `SKYML_META`. Measured on exposure 59750:

| Extension | Type / shape | MiB | Meaning |
|---|---|---:|---|
| `BLOCK_LSF_META` | table; 6 | 0.008 | Spectrograph, channel, wavelength/slit degree, offset half-width, solver status, fiber/pixel counts, chi-square, slit rank, selected fiberid range, normalization/centroid errors, minimum coefficient, timing. |
| `BLOCK_LSF_COEF` | f8; 6 × 11 × 4 × 4 | 0.011 | Block × offset basis × wavelength basis × slit basis; rows follow `BLOCK_LSF_META`. |
| `BLOCK_SLIT_KNOTS` | f8; 6 × 8 | 0.005 | Cubic slit knots over the spectrograph's full fiberid range. |
| `BLOCK_WAVE_KNOTS` | f8; 6 × 8 | 0.005 | Wavelength knots the block was fitted on; R/Z `LSF_KNOTS` of every Science fiber equal these. |
| `SKY_FIBER_COEF` | table; 106 | 0.321 | Selected Sky fibers: input_row, fiberid, spectrograph, telescope, frozen component vector (internal 1e14 scaling) and QP status. |

The block run computes its preparation in memory and writes it only with `--keep-preparation` (reusable with `--seed`; a supplied `--seed` is never modified). The product does not read it. `SKYML_META` adds `method`, `preparation` path/SHA256 (`in-memory` and null when not written), `block_settings`, `block_timing` and `block_code_sha256`. `SkyBlockReconstructor.global_state(fiberid)` evaluates the block surfaces at any fiber; `.sky_fiber(fiberid)` renders a selected Sky fiber's fit in memory.

## Storage boundary

The compact file contains no OBSFLUX copy, zero-filled SKY_IVAR, per-pixel SKYML_VALID, intermediate full sky images, individual full fit images, or median model/component/residual images. Validity is represented by SLITMAP.skyml_applied, MASK, IVAR and finite final SKY. Coefficients and observed median data reconstruct the model curves and additive corrections. Invalid final pixels require the original CFrame to recover actual observed values.

Model arrays use physical flux units (erg s^-1 cm^-2 Angstrom^-1). Stored decomposition/ML coefficients and nuisance spline parameters retain the fitter's 1e14 internal scaling. Rendering converts them to physical flux; report display uses 1e-14 units.

## History columns

`input_row`, `fiberid`, `channel`, `iteration`, `chi2_red`, `line_amplitude`, `solar_amplitude`, `relative_width_step`, `centroid_step`, `coefficients`, `probe_wave`, `w50`, `w90`, `continuum_knots`, `continuum_coef`, `variance_continuum_coef`, `huber_scale`, `continuum_knots_n`, `continuum_coef_n`, `variance_continuum_coef_n`.

Each variable-length spline vector is padded with NaNs and has an explicit `_n` length. R/Z knots and coefficient counts can differ. The standard history row identifies input row, physical fiber, channel and iteration. The last row for that fiber/channel reconstructs its final auxiliary fit; intermediate rows reconstruct spectral evolution. A valid survey-only/no-iteration product may have no FIT_LSF_HISTORY; diagnostics then show unavailable fitted models rather than inventing them.

## Verification and limits

All five compact products were checked against their previous full-diagnostic products for complete median models and representative individual sky/fit components/weights. Additional user-selected 59750 fibers 680, 1447 and 971 reproduce stored spectral/weight arrays to <2e-6 of each array's peak; measured errors are at float32 rounding scale. Report generation checks reconstructed final SKY for every processed fiber. These are serialization/reconstruction checks, not independent LSF-recovery accuracy measurements.

Native observed-noise IVAR is retained. Full sky-model uncertainty is uncalibrated; there is no all-zero SKY_IVAR image. Exact reconstruction uses the matching shipped/configured numerical model and telluric assets.

## Block-LSF support decisions

Updated block products retain all six entries in `BLOCK_LSF_META`, including
reference fallbacks. `adopted=True` means the fitted surface passed the block
support/solver checks; `False` means its stored coefficients are the projected
Science-median reference, repeated across pseudo-slit. `fallback_reason` records
the decision. `input_fibers`, `rejected_channel_fibers`, `slit_condition`,
`expected_max_gap` and `selected_max_gap` accompany the existing counts/rank/range.
Fallback χ² fields are NaN. These fields assess support, not physical recovery.
`centered` describes the stored block; `requested_centered` records the CLI request.
The primary `LSFCENT` remains the requested run mode. A reference fallback is
never silently recentered. Science `LSF_QA.accepted/converged=False` and its reason
identify reference fallback. `SKY_FIBER_COEF.usable/reason` preserve unavailable
individual fits explicitly; their coefficients are NaN. Old products without
the new columns remain readable as historical fitted blocks.

Block likelihood diagnostics additionally store `data_chi2`, `data_n_pixels`,
`nominal_df`, `effective_df`, `active_coefficients`, `active_threshold` and
`df_seconds`. Data chi-square is unpenalized and includes evaluation of fixed
reference fallbacks; their nominal/effective dimension is zero. The older `chi2`
fit-status field stays NaN for a reference fallback. Local effective degrees of
freedom condition on the current positivity face and fixed preparation/component
coefficients; they do not describe the complete decomposition pipeline.
`hessian_bytes`, `penalized_bytes`, `coefficient_bytes` and `sample_design_bytes`
give array payload sizes in bytes. These are distinct from FITS size and peak RSS.
