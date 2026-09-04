# Moon/Zodi Split — Identifiability Priors for `lsf-surface-iterative-split-zodi`

**Date:** 2026-09-03
**Branch:** `moon-zodi-shape-priors`
**Files changed:** `skysub/sky_decomp/fit.py`, `skysub/sky_decomp/moon_zodi_model.py`,
`skysub/decompose_parallel.py`, `notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi_module.ipynb`

---

## 1. Summary

The `split_zodi=True` decomposition was degenerate in a way that fit quality
could not detect. On 200 lunation-stratified sky spectra the deployed
configuration put the moon and zodi **colours in the wrong order on 164 of 168
moon-up spectra**, gave the moon family 45% of the continuum with the moon 37°
*below* the horizon, and produced a "zodi" whose amplitude tracked lunar
illumination at ρ = +0.95 while retaining almost none of its Leinert B500
dependence. The ML pipeline was being handed moon-labelled zodi spectra and
vice versa.

Four constraints fix this, now the defaults for that fit model:

| knob | value | fixes |
|---|---|---|
| `moon_ratio_bound` | 0.7 | colour ordering |
| `zodi_ratio_bound` | 0.7 | colour ordering (the decisive one) |
| `amp_prior_tol` | 3.0 | dark-time spurious moon |
| `zodi_amp_bound` | 2.0 | zodi amplitude geometry |

They are **not interchangeable** — see §4.

### Before → after (200 spectra, 168 moon-up, 28 dark)

| metric | before | after | target |
|---|---|---|---|
| colour reversals (moon-up) | 164 | **2** | 0 |
| separation margin, p10 | −7.71 | **+0.98** | > 0 |
| moon log-log slope | +0.21 | **−4.10** | −3.7 |
| zodi log-log slope | −4.04 | **−1.64** | −0.3 |
| moon share, moon >5° below horizon | 0.45 | **0.02** | ~0 |
| ρ(fitted zodi, Leinert B500) | 0.09 | **0.90** | ~1 |
| ρ(fitted zodi, FLI) | 0.94 | **0.12** | 0.10 |
| ρ(moon fraction, predicted) | 0.09 | **0.94** | 1 |
| fitted/predicted zodi | 4.34 | **2.00** | 1 |
| rms ratio, median / p90 | 1.000 | 1.010 / 1.457 | 1 |

Cost is concentrated entirely at bright moon (×1.49 median for FLI > 0.8,
≈1% of the continuum, blue-weighted) and is partly the *old* fit overfitting —
see §5.

---

## 2. Implementation

### `fit.py` — three constraint families, all hard linear inequalities

Appended to the same Clarabel non-negative cone as `c ≥ 0`. Being homogeneous
in `c` (except the absolute zodi bracket) they are scale-invariant and never
bid against the likelihood — unlike a quadratic penalty, which carries a hidden
`data_scale²` (documented on the surviving D² blocks: effective strength scales
as flux², so the same λ acts ~170× more strongly over a 13× brightness range).

1. **Adjacent-knot ratio bounds** — `β·c_k ≤ c_{k+1} ≤ c_k/β` on native
   coefficients, per family. Caps how fast a family's multiplier may vary with
   wavelength, so neither family can reproduce the other's colour. Also
   prevents a spline zeroing out mid-band (§5).
2. **Moon-fraction bracket** (`amp_prior_tol`, `amp_prior_floor=0.02`) —
   `f_lo ≤ ∫moon/∫(moon+zodi) ≤ f_hi` around the geometry prediction.
   Calibration-free: both carriers leave the physical model through the same
   solid-angle/flux conversion, so the ratio is independent of throughput and
   of `physical_to_fit_flux_scale`.
3. **Absolute Leinert bracket** (`zodi_amp_bound`) —
   `Z_pred/κ ≤ ∫zodi ≤ κ·Z_pred`. The only rows with a non-zero RHS, hence
   stated in native flux units; `zodi_total` must carry the same FACTOR as the
   flux being fitted.

Setters mirror `set_target_airmass`: `set_amplitude_prior(moon_fraction,
zodi_total)` installs scalars only (no design-matrix rebuild).

### `moon_zodi_model.py` — `geometry_amplitude_prior()`

Returns `(moon_fraction, zodi_total, target_airmass)` from a **single**
`predict` call so the shared conversion cancels in the ratio. Uses a cached
predictor with the learned scale factors zeroed **except**
`moon_below_horizon_suppression` (a geometric taper, not a calibration).

Why the other ten are dropped: they were fitted with
`CORRECTION_SCOPE = 'moon_plus_zodi'`, so only the *sum* was ever constrained
and the split they imply is uncalibrated. Worse,
`zodi_target_airmass_log = 1.28` (zodi growing as X^1.28, where extinction
should make it *fall*) looks like it absorbed the same moon-into-zodi leakage
the prior exists to remove. The physics-only and learned variants agree on the
moon fraction to a few percent, so this costs nothing.

### `decompose_parallel.py`

`SPLIT_ZODI_*` constants + `_install_split_zodi_amplitude_prior()` per
spectrum; split-zodi now loads LSF and META like the model-based mode.
Knot defaults set to **11 moon / 1 zodi** to match the validated configuration.

---

## 3. Why the prior is trustworthy

Three independent checks before any code was written:

* **Absolute scale.** Fitted moon+zodi total is 1.49× the prediction moon-up,
  1.06× in dark time — the model's absolute calibration is right to ~50%, so
  the continuum these splines fit really is mostly moon+zodi.
* **Dark-time validation.** With the shape bounds on and no moon,
  fitted/predicted `∫zodi` = **0.94 with 0.16 dex scatter**. That is an
  independent confirmation of the Leinert absolute scale on this data, and it
  is the licence the removed-2026-08-26 attempt lacked (it self-calibrated
  α on moon-down rows that carry the 43% spurious moon share).
* **Predicted zodi does not track lunation** (ρ = +0.10, p = 0.39) and tracks
  B500 at ρ = +0.996, exactly as physics requires.

**Do not tighten `zodi_amp_bound` below ~1.5** — that is the measured spread of
the anchor itself.

---

## 4. What each knob does, separately

Established by ablation, not assumption:

* **Reversals are fixed by the zodi shape bound, not the anchor.** With
  `zodi_ratio_bound` 0.5 → 0.7 reversals collapse 36 → 2. A run with the
  tightened bound and the anchor *off* also gets 2 reversals — but
  ρ(zodi, B500) = 0.125 and ρ(zodi, FLI) = 0.932.
* **Amplitude geometry is fixed by the anchor, not the bounds.** Only the
  anchor moves ρ(zodi, B500) 0.13 → 0.90 and ρ(zodi, FLI) 0.93 → 0.12.
* **The fraction bracket alone fixes dark time** (0.45 → 0.02) at
  rms ×1.001, and is free.
* **κ is front-loaded.** Turning the anchor on at all costs +0.24 in
  bright-moon rms for +0.42 in ρ(B500); κ=5 → κ=2 costs only +0.08 more for
  another +0.37. If used at all, use κ=2.

---

## 5. The rms comparison is biased toward the old fit

With β = 0 nothing stops a spline coefficient going to zero mid-band. The old
fit exploits this: on row 1442_near the baseline moon component **collapses to
~0.02 across 5100–6900 Å and re-emerges in the red**, while "zodi" carries the
mid-band and rises toward 9800 Å. In the 5400–5700 Å band the moon changes by
×64–×179 between configurations while the *total* continuum changes by ×1.10 —
pure relabelling. The ratio bounds forbid this, so part of the measured "+49%
rms at bright moon" is the baseline overfitting via a degree of freedom it
should not have.

---

## 6. Things tried that did NOT work

* **Two-channel scattered-moonlight envelope** (aerosol channel + multiple
  scattering + path terms, `moon_scatter_envelope`). Implemented, physically
  correct, and **neutral**: rms vs champion 1.000, bright-moon cost unchanged
  (1.488 → 1.486). Reason: 15 moon knots at β=0.7 can already manufacture the
  separation-dependent colour, so the envelope was never the binding
  constraint. At the model's own coefficient (3.5) it makes reversals *worse*
  (2 → 10). Code retained, **off by default**.
* **Blue relaxation of the moon spline** (`moon_ratio_relax_window`, plus a
  moon-dominance `moon_relax_gate`). The gate works correctly, but the
  relaxation buys ~1.5%: `b=0.3` and fully-free give bit-identical results
  because the moon spline never reaches its blue bound. Shape is not the
  binding constraint — the anchor is. Code retained, off by default.
* **A third smooth continuum component.** Withdrawn: the residual is
  oscillatory and sign-changing at the 1–3% level, not a smooth pedestal, so
  such a component would fit only the part degenerate with the diffuse family.
* **D¹ slope penalties.** Removed entirely (87 lines) — could not reach the
  colour targets without 27% flux drift.

---

## 7. Corrections made during this work (read before trusting older notes)

* **`res.resid` / `res.bestfit` are the PRE-refinement seed fit.** The refined
  model is `res.bestfit_lsf`; refined per-family spectra are in
  `res.components`. `rms_resid`, `reduced_chi2` and `r2` describe the *refined*
  fit. Reading `.resid` faked an "LSF bug" (apparent full-spectrum RMS 1.46 vs
  0.67 by seed) and a "systematic red over-prediction of −8…−22%". Both are
  artifacts: `rms(flux − bestfit_lsf)` is 0.5449 vs 0.5450 between seeds.
  **`lsf_sigma` is only a starting value — the iterative class converges to the
  same LSF surface regardless.** No LSF fix, no re-decomposition needed on that
  account.
* Config comparisons in this campaign used `rms_resid` and `res.coef`, both
  refined, so they are unaffected.
* Residual localisation, recomputed on the refined model: **96.5%** of the
  extra residual is on continuum pixels (not 89%), and 3900–5400 Å carries
  **66%** (not the 95% first reported).

---

## 8. Deployment notes

* **The design matrix is unchanged.** `set_target_airmass` is deliberately NOT
  applied: it rebuilds the zodi envelope per row, and eight reconstruction
  sites (`mlp_predictor/data.py`, `trainer.py`, `wavelengths.py`, four
  diagnostics cells) build their decomposer with no airmass — enabling it would
  silently desynchronise the basis the ML coefficients are defined against.
  Measured cost of omitting it: reversals 2 either way, ρ(zodi,B500) 0.903
  both, dark share 0.020 both, per-spectrum rms change 0.01% median.
* **Knot counts matter.** Ratio bounds are per adjacent knot *pair*, so the
  same β is looser with more knots. Everything here is validated at
  11 moon / 1 zodi. Re-validate the bounds if these change.
* **Corpus regeneration + ML retrain required**: both the split and the
  moon+zodi *sum* move (sum ratio q10 0.82, q90 1.48 per spectrum).

---

## 9. Open items

* **Bright-moon residual, unexplained.** Scales with FLI (ρ = +0.91), 96% on
  continuum pixels, blue-weighted, ~1% of the continuum. Left free, the fit
  puts zodi at 12–36× Leinert, so there is blue continuum under a bright moon
  that is neither Leinert-amplitude zodi nor solar-shaped moon. Not
  instrumental: the residual shape does **not** collapse across geometry
  (fractional spread 0.84, median pairwise shape correlation 0.38).
* **Telluric O₂ is the largest single defect** and is unrelated to this work:
  A-band (7590–7690 Å) sits at −6.5% of the continuum, B-band (6860–6950 Å) at
  +3.1%, systematically across spectra.
* **Step at the r|z arm boundary**: −4.7% of the continuum at 7453.5 Å
  (p10 −9.5%, p90 +0.6%) with O₂ excluded from the windows; b|r shows nothing
  (+0.22%). Looks like a reduction-side issue.
* **Row 12** (moon separation 18–25°) fits catastrophically, rms ~20× the
  sample median. Inside the model's own `near_moon < 20°` flag. Not chased.
* **ML side**: the Phase-F additive moon–zodi coupling and possibly the blend α
  may have been compensating for the role reversal. Worth re-testing whether
  they still earn their keep once the corpus is regenerated.

---

## 10. Reproduction

Scripts and CSVs live in the session scratchpad (not preserved):
`shape200.csv` (κ_zodi sweep), `shape200c.csv` (zodi bound sweep),
`shape200d.csv` (scatter envelope), `airmass_test.csv`,
`missing_shape.npy`, plus `bright_resid.png`, `bright_resid2.png`,
`resid_shape.png`, `missing_shape.png`.

Sample construction: 100 rows valid on **both** sky telescopes, stratified
14/18/18/18/18 across FLI bins with the moon >15°, plus 14 dark rows
(moon < −5°), spread over airmass within each bin; fitted on both SKY_NEAR and
SKY_FAR = 200 spectra. Metrics: colour = log-log slope of each fitted
component; reversal = `moon_slope > zodi_slope`; separation margin =
`zodi_slope − moon_slope`.
