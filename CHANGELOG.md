# Changelog

## Unreleased

### Changed

- Set a three-second default socket timeout before batch decomposition imports,
  including spawned workers, to limit optional DRP/Astroquery network waits.
- Batch decomposition and Sky Subtraction share CSV LSF loading, physical-fiber
  selection averaging, and fixed-Blue validation. Blue is held fixed; R/Z start
  from the same selected-fiber CSV means and are refined on the target spectrum.
- Reference mode is the batch default and uses the bundled `sky_decomp/data/lsf`
  calibration. CSV median inputs require `FIBER_ROLE`; explicit solar fitting
  and explicitly selected legacy FITS calibration directories remain supported.
- Reference provenance records the CSV manifest SHA256, averaging convention,
  membership use, and per-channel contributor counts in LSF metadata. Previous
  compact caches are invalidated by source and reference fingerprints.

### Validation

- 35 targeted tests passed in the Python 3.11 uv runtime: shared reference and
  worker/API/compact contracts (27), plus numerical fixed-Blue/R/Z refinement,
  profile-mean, strong-continuum and SkyML compact checks (8).
- No new full-corpus calculation or neural-network training was performed.
