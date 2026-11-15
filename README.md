# lvmsky
This code is for LVM sky subtraction

## Plotting and report rules

Figure styling, spectral residuals, interactive plots, and report publication layout are defined in [PLOTTING_AND_REPORT_RULES.md](PLOTTING_AND_REPORT_RULES.md). Read it before creating or updating a scientific report.



## Solar-constrained blue LSF

Use `palacecorr-aijc-vnf-split-zodi-lsf-spline2d-solarblue` for bright solar
continuum. The batch default fixes Blue to the bundled per-fiber CSV reference, matching Sky Subtraction.
The model renders sky lines and high-resolution solar absorption on the original
input wavelength grid. `--blue-lsf solar` fits Blue from those solar lines;
reference mode holds Blue fixed while still rendering through that profile.

The LSF representation and refinement settings are:

- B/R/Z: 11 offset M-splines × 4 wavelength B-splines, cubic, Δλ ±2.5 Å.
- Offset coefficient second-difference penalty: α=0.01, scaled by the mean
  data Hessian diagonal. Wavelength roughness, reference prior and numerical
  ridge are zero; continuum/OH priors and science-line masks are unchanged.
- Convergence at every native wavelength in every arm: barycenter change
  <0.015 Å, relative W50 and W90 changes <1%, relative data χ² change <1%.
  Maximum 30 cycles. `fit_summary` and channel metrics report convergence;
  a solved QP alone does not imply convergence.
- Nonnegative, normalized profiles with coefficient monotonicity constraints;
  these are not exact continuous derivative constraints at ±0.5 Å.

The server needs the existing `lvmdrp` / `TelluricCalculator` runtime and its
transmission assets under the configured `SAS_BASE_DIR`. Solar and PALACE
source tables are bundled in `skysub/sky_decomp/data`.

From the repository root in a configured Python 3.11 runtime. For an existing
uv environment, use `.venv/bin/python` in place of `python` in these commands:

```bash
python -m skysub.decompose_parallel /path/to/median-stack.fits \
  --fit-model palacecorr-aijc-vnf-split-zodi-lsf-spline2d-solarblue \
  --compact-only --n-workers 4 --chunk-size 1 \
  --output-dir /path/to/solarblue-results
```

Input is the existing median-stack format (WAVE, FLUX/LSF SCI/SKY_NEAR/SKY_FAR,
META), not a CFrame. `--limit 1` is useful for a server smoke run. Compact
outputs retain per-row coefficients, all three continuous LSF surfaces and
configuration, with resumable row caches. Source hashes include `solar_lsf.py`
to invalidate stale cache entries. Do not reuse an output directory belonging
to a different method. Each worker builds and caches its own solar tensor;
start with a small worker count and measure server memory before increasing it.
The default refinement limit is 30; convergence may stop earlier.

Python entry point: `skysub.sky_decomp.solar_lsf.SkyDecompPalaceCorrSolarBlueLSF`.
An explicit `offset_roughness_fraction=0.0` disables the offset penalty.

## SkyML SFrame

The compact output contract and measured extension sizes are in [skysub/sky_subtraction/FITS_STRUCTURE.md](skysub/sky_subtraction/FITS_STRUCTURE.md). Model/fit diagnostics are reconstructed from recorded coefficients by [reconstruction.py](skysub/sky_subtraction/reconstruction.py); full diagnostic spectral images are transient report data. The working output filename is `lvmSFrame-skyml-<exposure>.fits`.

## Shared reference LSF for batch decomposition and Sky Subtraction

Both paths use `sky_decomp.lsf_reference.load_fiber_references` and
`mean_reference` with the same bundled `skysub/sky_decomp/data/lsf` CSV bundle.
For each median, physical IDs from `FIBER_ROLE` select the contributor profiles.
Only measured profiles enter each channel's equal-weight coefficient mean.
Blue stays fixed; R/Z are seeded with those means and refined in situ. An empty
usable channel is recorded as a failed input, without a telescope-mean fallback.
CSV median stacks without valid `FIBER_ROLE` are rejected before fitting.

From any directory, with the repository's runtime available:

```bash
python /path/to/lvmsky/skysub/decompose_parallel.py /path/to/median-stack.fits \
  --blue-lsf reference \
  --lsf-calibration /path/to/lvmsky/skysub/sky_decomp/data/lsf \
  --n-workers 20 \
  --output-dir /path/to/new-decomposition
```

The first two options are defaults. `--blue-lsf solar` explicitly fits Blue;
an explicit legacy calibration directory containing `survey/lsf_survey.fits`
retains the previous FITS-reference behavior. CSV checksums and grid identity
are validated; the manifest identity and channel contributor counts are saved.
Use a new output directory when changing the reference. This change does not
make existing trained MLP checkpoints automatically compatible with a new corpus.

See [CHANGELOG.md](CHANGELOG.md) for pipeline changes.

The [sky subtraction procedure](SKY_SUBTRACTION_PROCEDURE.md) records the scan and median-build
commands reported by the user and the next batch-decomposition command.
