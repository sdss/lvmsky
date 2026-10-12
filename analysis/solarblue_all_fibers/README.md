# SolarBlue all-fiber run

This run applies the validated SolarBlue Sky Decomposition and continuous B/R/Z LSF recovery to every Science, SkyE, and SkyW fiber in each available DRP 1.3.2 exposure satisfying:

- Moon elevation > 10 degrees
- Moon FLI > 0.5

The selection snapshot, progress, and completed-exposure ledger are stored in `exposure_selection_audit.csv`, `status.json`, and `completed_exposures.csv` under the output directory. Compact parameter FITS files contain `COEF`, `COEF_ERR`, `LSF_COEF`, `LSF_KNOTS`, `LSF_META`, and the identifying `INPUT_META` table. Full-vector reference FITS files contain the same fit parameters plus the native wavelength grid, fitted spectra, residuals, uncertainties, covariance blocks, and all model-component vectors for the established nine Science fibers (three slit positions in each spectrograph).

The runner groups 32 exposures per resumable shard. An interrupted active shard resumes from its per-row compact cache; completed FITS products are retained and skipped.

```bash
python -m analysis.solarblue_all_fibers.run_all_fibers --workers 112 --exposures-per-shard 32
```
