# Frozen Moon/Zodi physical sub-bundle

This directory contains the four physical assets and frozen parameter manifest
used by `SkyDecompMoonZodiLSFSurfaceIterative`. The complete decomposition data
root is the parent directory; see `../README.md` for the PALACE inputs, original
Meftah source, citations, units, transformation details, redistribution notes,
and external-bundle usage.

Runtime code verifies `model_manifest.json` and all four asset SHA-256 digests
before evaluating the model. It never loads experiment checkpoints, imports
JAX, downloads data, or falls back to the historical solar-only Moon model.

The bundle retains only the 12 global parameters from the frozen 30,000-step
checkpoint. Training-local arrays with shapes `100 x 401`, `100 x 8`, and
`100 x 3` are deliberately excluded. The source model failed its original
decomposition non-regression gate, so its scientific status remains
`diagnostic_only_decomposition_nonregression_failed` until an independent
holdout assessment says otherwise.

The file
`meftah_solar_hrs_disk_integrated_v1_1_vacuum_velocity_step_2kms.npz`
contains 157,365 float64 samples from 3500 to 9999.9748 Angstrom on a 2 km/s
logarithmic vacuum grid. The explicit `velocity_step_2kms` replaces the former
ambiguous abbreviation `dv2`; it is a sampling interval, not a version number.
The source text digest, model checkpoint digest, native LVM-grid digest,
training selection, model formula, conventions, and per-asset attribution are
recorded in the manifest.
