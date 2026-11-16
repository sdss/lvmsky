# Changelog

## Unreleased

### Changed

- Queue block-CLI path lists without scanning CFrame headers or contents first.
  Workers resolve each input and output path independently, with batch progress
  active during preparation. Missing/invalid inputs follow per-worker policies.
- Remove CFrame content hashing from ordinary block calculations and use
  size/mtime metadata for CLI reuse fingerprints. Legacy preparation reuse keeps
  its hash check when no metadata identity exists. Old CLI fingerprints require
  explicit recomputation; metadata identity does not detect same-size/mtime edits.
- Set a three-second default socket timeout when importing `sky_decomp`, covering
  decomposition and Sky Subtraction workers and existing-product validation
  before DRP/Astroquery imports.
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

- 128 block-CLI, scheduling, block-LSF and compact-product tests passed in the
  Python 3.11 uv runtime, including deferred subprocess workers with successful,
  missing and malformed FITS inputs. Scientific fitting was stubbed in the new
  scheduling integration test; no full Utah exposure was recomputed.
- 35 targeted tests passed in the Python 3.11 uv runtime: shared reference and
  worker/API/compact contracts (27), plus numerical fixed-Blue/R/Z refinement,
  profile-mean, strong-continuum and SkyML compact checks (8).
- No new full-corpus calculation or neural-network training was performed.
