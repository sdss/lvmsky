# Niv-integrated decomposition v1

This experiment keeps Niv's frozen split Moon/Zodiacal and diffuse-continuum
priors while retaining the telluric-aware line models and continuous B/R/Z LSF
surface from the target branch.

Production method IDs:

- `adam25k-telluric-niv-continuum`: Adam25k OH, no PCA, 388 coefficients.
- `palace-aijc-vnf-pca30-niv-continuum`: PALACE-Aijc VNF plus the retrained
  signed individual-line PCA30 basis, 419 coefficients.

Both methods require the exact 12,401-pixel native wavelength grid. The Niv
continuum contract is fixed in `sky_decomp/niv_continuum.py`; incompatible
overrides raise instead of changing the method silently.

Rebuild the PCA30 asset from 1,000 deterministic far-sky spectra:

```bash
conda run -n lvmdrp_dev_311 python -m \
  skysub.experiments_niv_integrated_decomposition_v1.train_pca30 \
  /Users/ik52/obs/sas/sdsswork/users/u6058164/lvmsframe_median_stack_1.2.1_gaia1over100.fits \
  --workers 8
```

Run the resumable full-stack coefficient and continuous-LSF decomposition:

```bash
conda run -n lvmdrp_dev_311 python -m skysub.decompose_parallel \
  /Users/ik52/obs/sas/sdsswork/users/u6058164/lvmsframe_median_stack_1.2.1_gaia1over100.fits \
  --fit-model adam25k-telluric-niv-continuum \
  --compact-only --n-workers 8 --chunk-size 1 --max-in-flight 8 \
  --output-dir skysub/experiments_niv_integrated_decomposition_v1/outputs/full_adam
```

The hidden per-row NPZ cache makes an identical rerun resumable. A run
fingerprint covers the input, data manifest, code, and fit parameters. Invalid
telluric rows stay aligned in the final FITS with NaN coefficients and an
explicit error reason.

Train the merged ten-member Niv MLP variant on the new 388-column corpus:

```bash
conda run -n lvmdrp_dev_311 python -m \
  skysub.experiments_niv_integrated_decomposition_v1.train_mlp \
  /Users/ik52/obs/sas/sdsswork/users/u6058164/lvmsframe_median_stack_1.2.1_gaia1over100.fits \
  skysub/experiments_niv_integrated_decomposition_v1/outputs/full_adam \
  --reference-checkpoint /Users/ik52/Downloads/mlp_ensemble_split_zodi_cont.pt
```

The first training run also builds the frozen physical Moon/Zodiacal context
cache beside the input stack. It can be prepared explicitly:

```bash
PYTHONPATH=skysub:. conda run -n lvmdrp_dev_311 python -m \
  mlp_predictor.moon_model_cache \
  /Users/ik52/obs/sas/sdsswork/users/u6058164/lvmsframe_median_stack_1.2.1_gaia1over100 \
  --n-workers 8
```

The canonical executed report is
`notebook_niv_integrated_decomposition_output.ipynb`.
