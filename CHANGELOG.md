# Changelog

## Unreleased

### Changed

- `mlp_predictor`: optional split of the `ionospheric` group into one group per
  emission mechanism (`data.IONOSPHERIC_SPLIT`, `data.split_ionospheric_group()`),
  a per-group `head_dim` mapping in the model, and a `group_weight_overrides`
  trainer key. Off by default; the deployed ensemble loads and predicts as before.
  Used by `experiments_ionospheric_split_heads_v1` (not adopted).

- `mlp_predictor`: optional isolated airglow branches (`isolated_branches`,
  `isolated_branch_inputs` model/trainer keys): small per-group MLPs that bypass the
  trunk, read the full sci/near/far context plus chosen groups' near/far scores, and
  add to the group's alpha blend. Also `data.split_atomic_group()` moving [O I] 5577
  into its own `oi_green` group (95 km). Off by default; the deployed ensemble loads
  and predicts as before. Used by `experiments_ionospheric_branches_v1`.

- `mlp_predictor`: optional `fixed_blend_alpha` model/trainer key holding a group's
  near/far blend weight fixed instead of learned. Off by default.

- Extend the least-contaminated-fiber median fallback to SkyE/SkyW, targeting
  10 eligible fibers per Sky telescope (Science keeps 50). Persist fallback headers, selection counts
  and a prediction reliability warning. Retain finite `InsufficientProgress`
  median coefficients with their original status and a propagated fit-failure
  bit; invalid inputs and nonfinite coefficients remain fatal.
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
