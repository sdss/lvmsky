# LSF-block sky-subtraction CLI reference

Implemented interface: use the module invocation below or the checkout launcher ./skyml-block in the scientific environment. The original argparse invocation remains supported.

## 1. Current command and compatible single command

The current entry point uses standard-library argparse, not Typer:

~~~sh
python -m skysub.sky_subtraction_block /data/lvmCFrame-00011017.fits \
  --ml-model-file /models/ensemble.pt --gaia-dir /cache/gaia/fibers --nproc 8
~~~

It reads observed object-plus-sky flux from one calibrated lvmCFrame and calls pipeline.run(). It writes sky-subtracted flux and the applied sky, adopted LSF and reconstruction state to a compact lvmSFrame-skyml-block product. The current reader requires a CFrame; a previously sky-subtracted DRP SFrame is not an interchangeable input.

Keep this entry point, accepted arguments, defaults, validation, adjacent-state lookup and final-path stdout behavior compatible. Improve its English help descriptions while retaining its parser. Keep the Python pipeline.run() calling contract compatible.

The Typer application provides:

~~~text
python -m skysub.sky_subtraction_block.cli single INPUT [OPTIONS]
python -m skysub.sky_subtraction_block.cli batch [INPUTS...] [OPTIONS]
~~~

The checkout includes ./skyml-block. Add the checkout to PATH to use skyml-block directly. Examples below use this short name; the module invocation is equivalent.

The single command is the explicit name for the existing one-exposure operation. It accepts the existing arguments with their existing defaults and calls the same scientific runner. Batch selects exposures and schedules independent calls to that runner.

Provide -h/--help on the legacy command, new root, single and batch. CLI help, logs and documentation are English; examples precede the option groups.

## 2. Ordinary operation and optional saved-state reuse

### Ordinary operation

--ml-model-file is the canonical model-path option. The original --ensemble spelling remains accepted by legacy, single and batch commands with the same meaning.

~~~sh
skyml-block single /data/lvmCFrame-00011017.fits --ml-model-file /models/ensemble.pt
~~~

The ML model file option identifies an already trained ML model archive used to predict the Science IFU sky-component coefficients. The command performs inference; training is a separate operation. Use a trusted archive because loading restores Python/PyTorch objects.

An ordinary run computes Gaia-based selection, telescope medians and their decompositions, the ML prediction and correction state in memory. It then infers the block LSF, renders the applied sky and writes the final SFrame. **It writes no separate preliminary/preparation FITS by default.** The final compact product itself contains the normal reconstruction state.

### Special mode: reuse a previous same-exposure calculation

Saved-state reuse is useful when rerunning the block LSF with different basis dimensions, penalties or centering while holding the earlier median/ML calculation fixed. It avoids recomputing that stage and supports controlled comparisons.

Compatible inputs include:

1. An optional snapshot named lvmSFrame-skyml-preparation-00011017.fits, written by an earlier run that explicitly enabled --keep-preparation.
2. A compatible compact SkyML product from the same CFrame, such as lvmSFrame-skyml-00011017.fits or a previous block result, containing the required median/ML/correction state. An ordinary DRP SFrame lacks that SkyML state.

"Preparation" in the existing code means this internal state or its optional snapshot. It does not mean a required preliminary input or a default extra output.

Retain --seed FILE for backward compatibility. In the new single command, make --reuse-state FILE its descriptive canonical spelling, with --seed as an alias. Both select an existing same-exposure FITS state file. The value is not a random-number seed.

~~~sh
# Write a final product and an optional reusable snapshot.
skyml-block single /data/lvmCFrame-00011017.fits \
  --ml-model-file /models/ensemble.pt --keep-preparation \
  --output /results/reference/lvmSFrame-skyml-block-00011017.fits

# Reuse that snapshot for a controlled block-LSF rerun.
skyml-block single /data/lvmCFrame-00011017.fits \
  --reuse-state /results/reference/lvmSFrame-skyml-preparation-00011017.fits \
  --offset-roughness 0.01 \
  --output /results/trial/lvmSFrame-skyml-block-00011017.fits
~~~

Reuse keeps the saved selection, median decompositions, ML prediction and correction state; it recomputes block profiles and delivered spectra. The residual-transfer implementation uses the current residual-transfer settings. Input SHA256, physical fiber identities and wavelength grid must match. Supplied files remain read-only inputs. When reuse occurs, a simultaneously supplied ensemble does not recompute the saved ML prediction, and Gaia selection is reused too.

The existing automatic lookup remains compatible: when no explicit state file is supplied, the runner checks for lvmSFrame-skyml-preparation-<input suffix> beside its output. If found, it reuses that snapshot subject to identity checks. Otherwise it computes state in memory and requires --ml-model-file. Thus the default value of --seed is unset, while an existing adjacent snapshot can still select the special reuse path.

The log and dry-run plan must explicitly identify "compute state in memory", "reuse explicitly selected state" or "reuse adjacent saved state", including the path when applicable. The ordinary workflow needs only the ensemble and input; --keep-preparation is an explicit request for the extra snapshot.

## 3. Existing scientific options: defaults and detailed effects

Existing spellings remain supported by the legacy parser. Shared scientific options apply to batch as well; new aliases and CPU-budget flags belong to the Typer interface.

| Option | Default | Meaning and example |
| --- | --- | --- |
| INPUT | Required | One observed-flux CFrame path, e.g. /data/lvmCFrame-00011017.fits. Its wavelength grid and fiber identities must satisfy the existing scientific input contract. |
| --output FILE | Beside INPUT | Single-exposure final FITS path. The default replaces lvmCFrame- with lvmSFrame-skyml-block- and removes a trailing .gz. Explicit names must retain the lvmSFrame-skyml-block- prefix. Existing files are protected. |
| --ml-model-file FILE (alias --ensemble) | Unset | Trained ML archive for ordinary inference, e.g. /models/ensemble.pt. Required unless this exposure reuses compatible saved state. This does not enable writing an extra preliminary FITS. |
| --gaia-dir DIR | Unset | Direct parent of lvmGAIA-fibers-<exposure:08d>.fits, e.g. /cache/gaia/fibers. Exact schema and network behavior are in section 4. |
| --reuse-state FILE, alias --seed FILE | Unset | Special same-exposure reuse mode in section 2. Example: --seed /results/lvmSFrame-skyml-preparation-00011017.fits. The descriptive spelling is new; --seed retains its old meaning. |
| --keep-preparation | Off | Also save newly computed median/ML/correction state beside the final result as lvmSFrame-skyml-preparation-<input suffix>. Useful for subsequent controlled block reruns. Reused input state remains unchanged. |
| --nproc INTEGER | 8 in legacy/single; 1 in batch | Internal fiber-worker process limit per exposure on stages supporting it. At 1 those operations stay in the exposure process. Native numerical threads remain 1 per process. Section 8 defines the total CPU budget separately. |
| --slit-count INTEGER | 4 | Cubic spline controls across physical fiber position within a spectrograph; integer >=4. More controls allow finer slit variation and require stronger Sky-fiber coverage. Example: --slit-count 4. |
| --offset-count INTEGER | 11 | Cubic offset M-spline controls for R/Z profile shape; odd integer >=5. Support remains +/-2.5 Angstrom. This changes shape flexibility, not support width. Example: --offset-count 11. |
| --wave-count INTEGER | 4 | Cubic spline controls for R/Z profile variation with wavelength; integer >=4. More controls allow faster wavelength variation. Example: --wave-count 4. |
| --offset-roughness FLOAT | 0.0 | Dimensionless nonnegative penalty on squared second differences across offset coefficients. Larger values favor smoother coefficient patterns across profile shape. Zero disables this penalty. Example: --offset-roughness 0.01. |
| --wave-roughness FLOAT | 0.0 | Dimensionless nonnegative penalty on squared second differences along wavelength controls. Larger values favor smoother wavelength variation. Zero disables it. Example: --wave-roughness 0.01. |
| --slit-roughness FLOAT | 0.0 | Dimensionless nonnegative penalty on squared second differences along slit controls. Larger values favor smoother fiber-position variation. Zero disables it. Example: --slit-roughness 0.01. |
| --reference-prior FLOAT | 0.0 | Dimensionless nonnegative penalty on squared coefficient distance from the Science-median LSF projected onto the requested basis and repeated across slit position. Larger values pull profiles toward that reference. Zero disables this soft prior. Example: --reference-prior 0.01. |
| --centered | Off | Additionally constrain successfully fitted R/Z profiles to zero first moment. By default the centroid may vary. This does not introduce a wavelength-shift fit. A reference adopted for an unsupported block retains its own centroid and is identified in QA. |

All four strengths must be finite and >=0. They scale with the mean diagonal of the native-IVAR data Hessian, with a numerical floor, in the solver's quadratic-objective convention. The roughness terms penalize coefficient differences; they do not directly specify a physical smoothing length or FWHM. The reference prior is distinct from automatic reference adoption when Sky support is insufficient.

These options affect R/Z block fitting; B retains its adopted survey profile. Reported block data chi-square excludes penalties. The example strength 0.01 illustrates syntax rather than a calibrated recommendation. Increasing a basis count adds flexibility and computational cost; positivity and normalization constraints still apply. Definitions match regularized_normal() and the block-method README.

## 4. Gaia directory, files and actual fallback behavior

Pass the fibers directory itself:

~~~text
/cache/gaia/
  fibers/                          <-- --gaia-dir /cache/gaia/fibers
    lvmGAIA-fibers-00011017.fits
    lvmGAIA-fibers-00011018.fits
  sources/                         <-- not required for aggregate-cache input
    lvmGAIA-sources-00011017.fits
~~~

For exposure 11017, the lookup is exactly DIR/lvmGAIA-fibers-00011017.fits. The reader does not search subdirectories. These aggregate tables are produced by the existing lvm-medians fetch-gaia workflow described in the median-computation README.

Required content consumed by this pipeline:

- Primary-header EXPNUM equal to the CFrame EXPOSURE.
- A FIBERS table with unique fiberid values and the same physical fiber-ID set as SLITMAP.
- fiber_ra and fiber_dec matching SLITMAP.ra/dec in degrees, with absolute tolerance 1e-7 degree and matching NaN positions allowed.
- n_gaia_sources, gaia_g_n_valid_sources and gaia_g_flux measurements. Standard producer tables contain additional columns too.

Rows may be reordered; the reader aligns them by fiber ID. Observed G-band flux and Gaia-to-observed-flux ratios are recomputed from the current CFrame.

| Current condition | Action |
| --- | --- |
| --gaia-dir omitted | Query Gaia DR3 through ARI TAP. |
| DIR or the exact exposure file absent | Query ARI TAP. Passing /cache/gaia instead of /cache/gaia/fibers therefore causes a cache miss. |
| Readable table with mismatched exposure number, fiber-ID set or fiber coordinates | Record the mismatch and query ARI TAP. |
| Corrupt/unreadable FITS, absent FIBERS extension or absent/invalid consumed columns | Report a file/schema error. The current reader does not turn every malformed cache into a TAP fallback. |
| TAP succeeds | Use the source matches to derive fiber measurements; write the downloaded source table to OUTPUT_PARENT/.skyml-gaia/lvmGAIA-sources-00011017.fits. |
| All three TAP attempts fail | Report an exposure execution error; batch applies --on-error. |

Service: https://gaia.ari.uni-heidelberg.de/tap. A downloaded SOURCES table differs from the aggregate FIBERS cache supplied through --gaia-dir. Help must distinguish missing/identity-mismatched cache from malformed FITS/schema errors and must state the exact lookup path. This preserves the implemented behavior rather than promising a broader fallback.

## 5. Default paths from the SDSS SAS mirror

Use the local SAS mirror as the ordinary starting point, following the [LVM DRP directory convention](https://github.com/sdss/lvmdrp#running-the-drp):

~~~sh
export SAS_BASE_DIR=/data/sas
~~~

For the explicitly pinned project default version 1.3.2:

~~~text
$SAS_BASE_DIR/sdsswork/lvm/spectro/redux/1.3.2/
  drpall-1.3.2.fits
  1029XX/1029016/60325/lvmCFrame-00011017.fits
~~~

The version default reflects the current scientific input/grid contract. Selecting another path/version still requires compatible input and assets; the option alone does not establish support for a new reduction format.

| Batch path option | Default | Meaning and example |
| --- | --- | --- |
| --sas-base-dir DIR | Environment SAS_BASE_DIR | Local SAS mirror root containing sdsswork/, e.g. --sas-base-dir /data/sas. |
| --drp-version TEXT | 1.3.2 | Chosen reduction version for catalogue/data discovery. Example: --drp-version 1.3.2. No automatic newest-version selection. |
| --cframe-root DIR | Resolved version directory below | Advanced override for a nonstandard mirror: directory directly containing <tilegrp>/<tileid>/<mjd>, e.g. /scratch/lvm/1.3.2. |
| --drpall FILE | CFRAME_ROOT/drpall-<drp-version>.fits in catalogue mode | Exposure summary FITS. Example override: --drpall /catalogues/drpall-1.3.2.fits. A relocated catalogue does not move the data root. |

CFrame-root precedence: explicit --cframe-root; root constructed from explicit --sas-base-dir; environment LVM_SPECTRO_REDUX plus the chosen version; environment SAS_BASE_DIR plus sdsswork/lvm/spectro/redux/<version>. Explicit command-line paths are authoritative. An explicit --sas-base-dir also sets SAS_BASE_DIR for scientific assets; use --cframe-root when only the input-data mirror differs. The plan prints both resolved catalogue and CFrame-root paths.

For example, LVM_SPECTRO_REDUX=/data/sas/sdsswork/lvm/spectro/redux yields the same version root as SAS_BASE_DIR=/data/sas. Standard use needs no root option. Nonstandard storage can supply --cframe-root and --drpall explicitly. Known CFrame paths or an input-list do not require a SAS environment or catalogue.

A catalogue location may name a DRP SFrame and include an archive prefix. Extract its validated <tilegrp>/<tileid>/<mjd> suffix, place it below the selected version root and use the CFrame basename. Substitute the product basename prefix, rather than arbitrary directory substrings. Validate row/header exposure identity and version/path consistency.

## 6. Selecting a batch and defining its order

Choose one source: explicit INPUTS, --input-list, or catalogue selectors/--all-exposures. Catalogue filters apply only to catalogue mode. A catalogue path alone requires an explicit selection; it never schedules the entire catalogue silently.

| Batch selection option | Default | Meaning with concrete examples |
| --- | --- | --- |
| INPUTS... | Empty | Explicit files: batch /data/lvmCFrame-00011017.fits /data/lvmCFrame-00011018.fits. Shell-expanded globs become multiple paths. |
| --input-list FILE | Unset | One CFrame path per line. Example: batch --input-list cframes.txt. Relative paths are relative to that list file. Blank lines and lines starting with # are allowed. |
| --expnum VALUES | No exposure-number cut | One: --expnum 11017. Several: --expnum 11017,11018,11025. Inclusive range: --expnum 11017:11025. Combined: --expnum 11017,11020:11025. Repeatable. |
| --mjd VALUES | No night cut | Examples: --mjd 60325; --mjd 60325,60326; --mjd 60356:60401. Uses the catalogue observing-night identifier. Repeatable. |
| --tileid VALUES | No tile cut | Examples: --tileid 1029016; --tileid 1029016,1029017; --tileid 11112: for IDs >=11112. Repeatable. |
| --exptime SECONDS | No duration cut | Match the recorded exposure time, e.g. --exptime 900. Selects rows; does not rescale flux. |
| --catalog-fluxcal TEXT | No calibration-mode cut | Match drpall.fluxcal, e.g. --catalog-fluxcal MOD. Catalogue filtering only; does not calibrate or rescale flux. |
| --all-exposures | Off | Every exposure row in the selected catalogue. Example: batch --all-exposures --drp-version 1.3.2. Use instead of catalogue filters. |
| --max-exposures N | Unset: whole selected set | First N entries of the ordered selected list; N must be an integer >=1. Example: --mjd 60325 --max-exposures 3 selects that night's first three entries by the ordering below. |
| --missing-input error or skip | error | error stops planning if a selected local CFrame is unavailable; skip records it and continues with available members of the same selected set. |

Commas and repeated occurrences mean OR within a selector. Different selectors combine with AND. Ranges are inclusive, and open ranges are supported. There are no implicit exposure-time, MJD, tile or calibration-mode cuts. The old driver's strict 60355 < mjd < 60402 becomes --mjd 60356:60401.

Plan in this order:

1. Read the chosen source and apply requested filters.
2. Deduplicate identical exposure records and report conflicting records.
3. Sort by DRP version, numeric MJD, numeric exposure number, numeric tile ID, then normalized path. This is reproducible ordering, not a quality ranking or random sample.
4. Apply --max-exposures if supplied. Otherwise retain the whole selected list.
5. Resolve and validate input/output/state paths, missing inputs, existing products and collisions.
6. Schedule the resulting execution list.

For [11017, 11018, 11025, 11030], --max-exposures 2 retains [11017, 11018]. Later skips are recorded without replacing them with later entries. A requested exposure absent from the catalogue and an empty initial selection are planning errors. Initial catalogue support targets the current DRP 1.3.2 FITS workflow; HDF5 needs its own reader/schema contract.

The catalogue-calibration filter matches the existing driver's t['fluxcal'] == 'MOD'. It describes the DRP calibration mode that produced the input, not a sky-subtraction mode or a new calibration step. The new name --catalog-fluxcal makes that scope explicit. Leave it unset for ordinary processing across modes; request MOD only to select that metadata subset. Match actual catalogue values rather than inventing a universal mode list across versions. All selected inputs must still satisfy the physical-flux/grid contract.

## 7. Final output placement and optional batch reuse

| Batch option | Default | Meaning and example |
| --- | --- | --- |
| --output-dir DIR | ./skyml-block-output | Root for generated products, relative to the launch directory. Example: --output-dir /results/skyml-block. |
| --preserve-drp-tree / --flat-output | Preserve DRP tree | Mutually exclusive choices: recreate VERSION/TILEGROUP/TILE/MJD below output-dir, or write final FITS directly there. |
| --on-existing error or skip | error | error reports a destination conflict. skip reuses a readable complete product only when input identity and scientific fingerprint match; mismatches remain conflicts. |
| --reuse-state-dir DIR | Unset | Special reuse mode: existing per-exposure snapshots from --keep-preparation, with exact placement below. |
| --run-dir DIR | OUTPUT_DIR/.batch/<UTC run ID> | Plan, status, results and per-exposure logs, e.g. --run-dir /logs/skyml-night60325. |

~~~text
Input:
  $SAS_BASE_DIR/sdsswork/lvm/spectro/redux/1.3.2/1029XX/1029016/60325/lvmCFrame-00011017.fits

--output-dir /results --preserve-drp-tree:
  /results/1.3.2/1029XX/1029016/60325/lvmSFrame-skyml-block-00011017.fits

--output-dir /results --flat-output:
  /results/lvmSFrame-skyml-block-00011017.fits
~~~

The tree contains generated products and sidecars; input files are not copied. Product filenames have no automatic version suffixes. Whole-plan collision checks include flat-output collisions. Preserve-tree output needs validated layout metadata; explicit paths without it can use flat output.

With --reuse-state-dir /previous, preserve-tree mode expects:

~~~text
/previous/
  1.3.2/1029XX/1029016/60325/
    lvmSFrame-skyml-preparation-00011017.fits
~~~

Flat mode expects:

~~~text
/previous/
  lvmSFrame-skyml-preparation-00011017.fits
  lvmSFrame-skyml-preparation-00011018.fits
~~~

These snapshots were explicitly written beside prior final results using --keep-preparation. Batch basenames use eight-digit exposure numbers. For exposures needing execution, an explicit reuse directory must provide each matching snapshot; planning reports the exact unavailable or incompatible path. Existing SHA256/fiber/wavelength checks apply. This directory is an input for the special reuse mode, not a default output directory.

Without this option, exposures compute state in memory or follow the compatible adjacent-snapshot lookup. A mixed batch needs an ensemble for exposures computing state anew. --keep-preparation can save those newly computed snapshots beside their final results.

## 8. CPU budget and workers per exposure

Use two controls: `--cpus C` is the whole-run logical-CPU budget; `--nproc N`
is the worker-process count per exposure. Batch automatically runs up to
`floor(C / N)` exposures at once on disjoint N-CPU subsets. N must be <= C.
Batch defaults to N=1 and, when C is omitted, C=N (one concurrent exposure).
Single/legacy retain N=8; single without C inherits the allowed CPU set.
`--jobs` remains a compatible alias for `--nproc`. There is no `--processes` option.

| Options | Concurrent exposures | Workers per exposure |
| --- | ---: | ---: |
| --cpus 16 --nproc 1 | 16 | 1 |
| --cpus 16 --nproc 2 | 8 | 2 |
| --cpus 16 --nproc 4 | 4 | 4 |
| --cpus 7 --nproc 2 | 3 | 2 |

Serial stages use one worker. Each native numerical pool uses one thread.
The coordinator and descendants share the total CPU affinity. CPU allocation
bounds usage; it does not guarantee continuous utilization or linear speedup.
Memory scales with concurrent exposures. Fewer queued exposures use fewer slots.
Optional `--cpu-list IDS` selects exactly C scheduler-allowed CPU IDs.

## 9. Run controls and help requirements

| Option | Default | Detailed meaning |
| --- | --- | --- |
| --on-error continue or stop | continue in batch | continue records an exposure error and schedules independent remaining work. stop stops new submissions after an error while active exposures finish. |
| --dry-run | Off; batch | Resolve selection, exact paths, reuse decisions and CPU allocation; print the plan and exit before fitting, TAP requests or output creation. |
| --log-level LEVEL | INFO; new commands | DEBUG: detailed diagnostics; INFO: ordinary stages/progress; WARNING: conditions needing attention; ERROR: execution errors; CRITICAL: fatal application conditions. Sets console/file logging threshold. |
| --no-progress | Off; batch | Hide terminal progress while retaining logs and status. Plain output is automatic for unsupported/noninteractive terminals. |
| -h, --help | Show help and exit | Root and command documentation with defaults, valid values, units, examples and combined-option rules. |
| --version | Print and exit | New root: package version and implementation revision when available. |
| --install-completion | Unset | New root: install Typer shell completion. |
| --show-completion | Unset | New root: print completion code. |

Help/completion work without scientific-pipeline imports, a model archive or an available fitting runtime. Root help explicitly identifies single as the compatible name for the original one-exposure operation. Single help includes all option meanings and the cache/reuse details above, its J=8 default and its added optional CPU budget. Legacy help improves descriptions without changing accepted invocations.

For every option show its effect, datatype/units, accepted choices or bounds, default and a concrete syntax example where needed. All advanced options remain visible. Defaults shown in help must come from the same definitions as parser defaults.

### Batch help structure

~~~text
Usage: skyml-block batch [OPTIONS] [INPUTS...]

Subtract sky from calibrated CFrame exposures; write one final SkyML SFrame each.
Default batch execution: one active exposure, one computational thread.
Ordinary operation writes final products; saved-state reuse is optional.

EXAMPLES
  skyml-block batch --expnum 11017,11018 --ml-model-file /models/ensemble.pt
  skyml-block batch --mjd 60325:60330 --ml-model-file /models/ensemble.pt \
    --cpus 8 --output-dir /results --preserve-drp-tree
  skyml-block batch --input-list cframes.txt --ml-model-file /models/ensemble.pt \
    --flat-output --cpus 4
  skyml-block batch --all-exposures --ml-model-file /models/ensemble.pt --dry-run

INPUT AND SDSS PATHS
  INPUTS...                 Explicit CFrame paths: a.fits b.fits.
  --input-list FILE         One path per line; relative to that list file.
  --sas-base-dir DIR        SAS mirror root containing sdsswork/.
                            Example: /data/sas. [env: SAS_BASE_DIR]
  --drp-version TEXT        Catalogue/data version. [default: 1.3.2]
  --cframe-root DIR         Override version directory containing
                            <tilegrp>/<tileid>/<mjd>.
  --drpall FILE             Summary FITS. [catalogue-mode default:
                            CFRAME_ROOT/drpall-<drp-version>.fits]

SELECT CATALOGUE EXPOSURES
  --expnum VALUES           11017 or 11017,11018 or 11017:11025.
                            Inclusive ranges; repeatable. [no cut]
  --mjd VALUES              60325 or 60325:60330; repeatable. [no cut]
  --tileid VALUES           1029016,1029017 or 11112:; repeatable. [no cut]
  --exptime SECONDS         Match duration, e.g. 900. [no cut]
  --catalog-fluxcal TEXT    Match drpall.fluxcal, e.g. MOD. [no cut]
                            Catalogue selection only; no flux recalibration.
  --all-exposures           Every row in the selected catalogue. [off]
                            Use instead of filters; not needed for INPUTS.
  --max-exposures N         First N entries after selection, deduplication and
                            version/MJD/expnum/tile/path sorting. [all selected]
                            N must be an integer >=1.
                            A reproducible subset, not a quality ranking.
  --missing-input POLICY    error: stop planning; skip: record unavailable input.
                            [choices: error, skip; default: error]

FINAL PRODUCTS
  --output-dir DIR          Generated product root. [default: ./skyml-block-output]
  --preserve-drp-tree       VERSION/TILEGROUP/TILE/MJD below output-dir. [default]
  --flat-output             Final FITS directly in output-dir.
  --on-existing POLICY      error: conflict; skip: verify compatible product.
                            [choices: error, skip; default: error]
                            No automatic filename version suffixes.

MODEL, GAIA AND OPTIONAL REUSE
  --ml-model-file FILE           Trained ML archive for ordinary inference.
                            Required for exposures computing state anew.
                            Separate preliminary FITS is not written by default.
  --gaia-dir DIR            Direct parent of lvmGAIA-fibers-EXPNUM.fits.
                            EXPNUM has eight digits, e.g. 00011017.
                            Example: /cache/gaia/fibers. [unset: ARI TAP query]
                            Requires PRIMARY.EXPNUM and a FIBERS table.
                            Fiber IDs/coordinates must match the input SLITMAP.
                            Missing file or identity mismatch also queries TAP.
                            Service: https://gaia.ari.uni-heidelberg.de/tap
                            Corrupt FITS/invalid schema is a cache error.
  --reuse-state-dir DIR     Advanced: read optional snapshots from earlier
                            --keep-preparation runs, in the output tree style.
                            Basename: lvmSFrame-skyml-preparation-00011017.fits
                            [unset: ordinary/compatible adjacent-state lookup]
  --keep-preparation        Also save a new median/ML/correction snapshot beside
                            the final product for future controlled reruns. [off]

R/Z BLOCK-LSF FLEXIBILITY AND REGULARIZATION
  --offset-count N          Cubic offset controls over +/-2.5 Angstrom.
                            [odd integer >=5; default: 11]
  --wave-count N            Cubic wavelength controls. [integer >=4; default: 4]
  --slit-count N            Cubic fiber-position controls. [integer >=4; default: 4]
  --offset-roughness X      Squared second differences across offset coefficients;
                            larger X favors smoother profile coefficient patterns.
  --wave-roughness X        Squared second differences along wavelength controls;
                            larger X favors smoother wavelength variation.
  --slit-roughness X        Squared second differences along slit controls;
                            larger X favors smoother fiber-position variation.
  --reference-prior X       Squared coefficient distance to the projected
                            Science-median LSF, constant across slit position.
                            Larger X pulls toward that reference.
                            All four strengths: finite dimensionless X>=0.
                            Default 0 (off); scaled by mean data-Hessian diagonal.
                            R/Z only; reported data chi-square excludes penalties.
  --centered               Zero first moment for fitted profiles. [off]
                            Adopted reference fallback retains its own centroid.

PARALLELISM AND TOTAL CPU BUDGET
  --cpus C                  Total logical CPUs; default: N.
  --nproc N                 Workers per exposure; default: 1 in batch.
                            Concurrent exposures: floor(C / N). Require N <= C.
  --cpu-list IDS            Optional exact scheduler-allowed CPU IDs.

RUN CONTROL AND LOGGING
  --on-error POLICY         continue: schedule remaining independent exposures;
                            stop: finish active work and stop new submissions.
                            [choices: continue, stop; default: continue]
  --run-dir DIR             Plan, results, status and per-exposure logs.
                            [default: OUTPUT_DIR/.batch/<UTC run ID>]
  --dry-run                 Print resolved plan; exit before fitting, network
                            queries and output creation. [off]
  --log-level LEVEL         DEBUG, INFO, WARNING, ERROR, CRITICAL. [default: INFO]
                            DEBUG: diagnostics; INFO: stages; WARNING: attention;
                            ERROR: execution errors; CRITICAL: fatal conditions.
  --no-progress             Hide progress; retain logs and status. [off]
  -h, --help                Show this documentation and exit.

PATH AND SELECTION RULES
  CFrame-root precedence: --cframe-root; explicit --sas-base-dir plus standard
  suffix; LVM_SPECTRO_REDUX/<version>; SAS_BASE_DIR plus standard suffix.
  Standard suffix: sdsswork/lvm/spectro/redux/<version>.
  For a full catalogue request --all-exposures explicitly.
  Commas/repeated values mean OR; different selectors mean AND.
  No max-exposures means all SELECTED entries, not the whole catalogue.

CPU EXAMPLES
  --cpus 16 --nproc 1        Up to 16 exposures, one worker each.
  --cpus 16 --nproc 2        Up to eight exposures, two workers each.

~~~

A real run records plan.ecsv, status.json, results.ecsv and individual logs. Preserve scientific QA states separately from execution status: an adopted reference LSF is not itself an execution error. Reuse fingerprints include input identity, state/ensemble identity, scientific settings and relevant implementation provenance; scheduling parameters are recorded separately.

Exit codes: 0 successful completion/verified reuse; 1 exposure execution errors; 2 invalid arguments/planning conflicts; 130 user interruption. Stop new submissions and terminate owned exposure processes on interruption, finalize status and retain atomic final-product publication. A subsequent --on-existing skip can reuse compatible completed products.

## 10. Complete examples and implementation validation

Standard SAS mirror, preserved output tree and one CPU per exposure:

~~~sh
export SAS_BASE_DIR=/data/sas
skyml-block batch \
  --expnum 11017,11018,11020:11025 \
  --ml-model-file /models/mlp_ensemble_split_zodi_cont-1.3.2.pt \
  --gaia-dir /cache/gaia/fibers \
  --output-dir /results/skyml-block --preserve-drp-tree \
  --nproc 1 --cpus 8
~~~

Nonstandard data/catalogue locations with deliberate catalogue cuts:

~~~sh
skyml-block batch \
  --cframe-root /scratch/lvm/1.3.2 --drpall /catalogues/drpall-1.3.2.fits \
  --mjd 60356:60401 --tileid 11112: --exptime 900 --catalog-fluxcal MOD \
  --ml-model-file /models/ensemble.pt \
  --output-dir /results/selected --flat-output \
  --nproc 1 --cpus 10 --cpu-list 8-17
~~~

Release checks cover legacy argv/defaults and improved help, single/legacy scientific equivalence, input-source rules and selectors, deterministic maximum selection, SAS/environment/explicit-path precedence, actual Gaia cache/error cases, same-exposure state checks, mirrored/flat collision handling, compatible-product reuse, atomic publication and interruptions. Test P/J/C arithmetic and observe the actual process-tree affinity and native threads during both J=1 and explicitly parallel exposure runs.

The revision was checked against the argparse entry point, runner, Gaia producer/consumer, regularization code, shared threading bootstrap and SDSS directory documentation. Local real-data commands, measured process resources and product checks are stored under `analysis/cli_validation/` in the validation workspace.
