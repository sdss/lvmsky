# Baseline-decomposed Moon and Zodi continuum v1

This experiment keeps the existing 356-spectrum strict lunar cohort and the
existing analytic Moon and zodiacal-light model. It changes only the training
target and pixel weights:

1. Run the baseline wavelength-dependent LSF decomposition independently for
   every selected spectrum.
2. Subtract the fitted OH, atomic, Orc, and O2 emission-line components from
   the observed spectrum on the unchanged 12,401-pixel LVM wavelength grid.
3. Fit one shared JAX Moon+Zodi parameter vector to all line-subtracted spectra.
4. Use the continuum weights returned by the baseline decomposition. Pixels
   covered by skylines receive the baseline relative multiplier `0.0005` and
   therefore contribute very little to the global objective.

The decomposition is frozen preprocessing. The deployable Moon+Zodi forward
model and global fitting objective are JAX-compatible. No wavelength binning,
resampling, smoothing, or per-exposure fitted normalization is used.

Run from `lvmsky/` in the `lvmdrp_dev_311` environment:

```bash
python -m skysub.experiments_moon_zodi_baseline_decomp_jax_v1.model all
```

Regenerate the native-grid and 20 Angstrom OOF galleries without refitting:

```bash
python -m skysub.experiments_moon_zodi_baseline_decomp_jax_v1.model render
```

Outputs are written only below this experiment's `outputs/` directory.
