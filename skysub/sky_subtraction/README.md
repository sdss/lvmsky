# CFrame sky subtraction (SkyML)

Subtracts sky from one calibrated `lvmCFrame` and writes `lvmSFrame-skyml-<exposure>.fits` beside it. The product contract (HDUs, shapes, history columns, reconstruction API) is in [FITS_STRUCTURE.md](FITS_STRUCTURE.md).

## Running

From the repository root, in the lvmdrp Python 3.11 environment (`lvmdrp311` locally, `lvmdrp_dev_311` on the server; see the root README). It provides `lvmdrp`, `TelluricCalculator`, `SAS_BASE_DIR`/`LVMCORE_DIR` assets and the predictor dependencies:

```sh
python -m skysub.sky_subtraction /path/to/lvmCFrame-00012630.fits \
  --ensemble /path/to/mlp_ensemble_split_zodi_cont-1.3.2.pt \
  --gaia-dir /path/to/gaia/fibers --jobs 6
```

- Existing destinations and standard DRP filenames are never overwritten; `--output` selects another path.
- `--lsf-calibration` selects another survey ASCII LSF directory on the same wavelength grid (bundled by default).
- `--survey-only` keeps survey R/Z LSFs (validation baseline).
- Without a matching `--gaia-dir` cache the ARI TAP Gaia DR3 query runs (three attempts); its source table is cached in `.skyml-gaia/` beside the output and the cache-miss reason is logged and recorded.
- The ensemble archive must be trusted: its loader restores Python/PyTorch objects.
- The CLI limits native threads to one per worker. When calling `pipeline.run()` from Python, set `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 before importing numpy.

Supported inputs are calibrated CFrames on the bundled native DRP 1.3.2 wavelength grid with physical fiber IDs.

## Processing

1. Join Gaia to SLITMAP by fiber ID, recompute the observed G-band flux, select Gaia-clean fibers (Science falls back to its 50 lowest-ratio fibers) and form masked medians with contributor counts.
2. Fit the Science, SkyE and SkyW medians with the corpus decomposition: Blue fixed, R/Z seeded with the selected fibers' mean survey LSF.
3. Predict one coefficient vector at the Science IFU centre. Estimate the Blue sky-arm residual transfer (tapered off at 5000–5100 Å) and common OH band scales/upper-level-energy tilts on the Science median. OH proposals with negative populations revert to the raw prediction.
4. Fit the R/Z LSF of every usable Science fiber against the full sky model (emission plus lunar/zodiacal continuum) with a nonnegative source spline; adopt it only if it passes the checks below, else keep the fiber's survey LSF.
5. Render each fiber's sky with its adopted LSF, add both corrections, subtract and publish atomically.

The source spline and the local line/solar amplitudes are nuisance parameters; they never enter the subtracted sky.

## LSF acceptance

All thresholds are fields of [`Settings`](settings.py) and are recorded in `SKYML_META.settings`. Checks are evaluated at nine wavelength probes per channel, in this order; `LSF_QA.reason` names the first failure:

| Reason | Condition |
|---|---|
| `insufficient_lines` | fewer than 5 line pixels per LSF coefficient, or a wavelength basis function with ≤5 supporting pixels |
| `continuum_solver`, `lsf_solver_*` | a solver failed; the survey LSF is kept and the candidate holds the last successful iteration |
| `invalid_profile` | coefficients violate nonnegativity, unit sum or central-component shape |
| `line_information` | each of four equal wavelength quarters needs ≥2 observed sky lines ≥10 Å apart with fitted peak/continuum ≥1 and peak S/N ≥5 (local intercept and slope removed; contrast against the larger of observed-sideband and model-continuum 75th percentiles, with a noise floor; input masks apply) |
| `not_converged` | W50/W90 step ≥1% or centroid step ≥0.015 Å after 12 iterations |
| `profile_guard` | FWHM outside 0.9–3.0 Å, W50/W90 change >35% from the survey reference or >50% from the Science median, or centroid change >0.25 Å |
| `insitu` | accepted |
| `fiber_error:<type>` | unexpected error in this fiber's fit; logged, survey LSF kept |
| `unusable_survey_fallback` | rejected and the survey reference lacks measured quality: the channel's sky is NaN, MASK set and IVAR zero |

These are conservative operating points, not a calibrated recovery probability. Per-quarter counts, median contrast and median S/N are stored in LSF_QA.

## Approximations

- The equal mean of the selected fibers' LSF coefficients stands in for the LSF of their masked median; contributor sets can differ by wavelength.
- Median IVAR is the independent heteroscedastic Gaussian asymptote, without spatial covariance or intrinsic fiber-to-fiber scatter.
- IVAR describes observed noise only. Reconstructed SKY_EPI propagates the diagonal ML coefficient variance; full sky-model uncertainty is uncalibrated.
- Bundled survey references are quality-screened per-fiber/channel averages; epoch variation remains possible.

## Code layout

| Module | Role |
|---|---|
| `__main__.py` | CLI |
| `settings.py` | operating points (frozen dataclass) |
| `io.py` | CFrame reading, Gaia photometry, fiber selection, medians |
| `decomp_adapter.py` | the only calls into private `decompose_parallel`, `sky_decomp` and `lvm_medians` helpers |
| `ml.py` | the only import point of `mlp_predictor` (needs `skysub/` on `sys.path`) |
| `rendering.py` | sky rendering, channel bounds and robust weights shared by fit and reconstruction |
| `fiber_lsf.py` | per-fiber R/Z LSF fit and its design tensors |
| `lsf_information.py` | observed line-information gate |
| `pipeline.py` | exposure stages and the fiber process pool |
| `product.py` | FITS assembly and atomic publication |
| `reconstruction.py` | `SkyMLReconstructor` and `diagnostic_hdul` |

## Validation

```sh
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 python -m pytest -q skysub/sky_subtraction/tests
python -m ruff check skysub/sky_subtraction && python -m ruff format --check skysub/sky_subtraction
```

The controlled LSF tests use known profiles under broad/narrow stellar structure and a documented STScI Pickles M5V fixture. The report `analysis/sky_subtraction_validation/00_MAIN_sky_subtraction_report.html` shows direct SFrame FLUX comparisons, individual fiber fits, line-scale windows, LSF surfaces and iteration histories; regenerate it with `python analysis/sky_subtraction_validation/report_fibers.py` (`run_manifest_v2.json` lists the five outputs) and publish with `publish_report.py`. Presentation rules: [PLOTTING_AND_REPORT_RULES.md](../../PLOTTING_AND_REPORT_RULES.md). The sample demonstrates execution and quality controls; it does not establish survey-wide generalization or independence from ML training. The development history is in `analysis/sky_subtraction_validation/history/IMPLEMENTATION_PLAN.md`.
