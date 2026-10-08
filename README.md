# lvmsky
This code is for LVM sky subtraction


## Solar-constrained blue LSF

Use `palacecorr-aijc-vnf-split-zodi-lsf-spline2d-solarblue` for bright solar
continuum. The ordinary production method remains the CLI default.
The new method fits sky lines in all arms and high-resolution solar absorption
in B with one fitted convolution, on the original input wavelength grid.

Defaults in both the Python class and decomposition worker:

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

From the repository root in the `lvmdrp_dev_311` environment:

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
The default CLI cycle value 5 maps to 30 for solarblue; any other positive
`--n-refinement-cycles` value sets its safety limit explicitly.

Python entry point: `skysub.sky_decomp.solar_lsf.SkyDecompPalaceCorrSolarBlueLSF`.
An explicit `offset_roughness_fraction=0.0` disables the offset penalty.
