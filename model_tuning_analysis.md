# LVM Sky Model Tuning Analysis

**Date:** 2026-09-01

## Overview

Analysis of the LVM sky interpolation notebook: `notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi_module.ipynb`

This document summarizes the model architecture, diagnostic findings from cell 28, and specific tuning recommendations to improve model performance on problematic moon-state regimes.

---

## Model Summary

### Architecture: Symmetric Dual-Encoder with Physically-Grouped Heads

The model predicts sky coefficients (moon, zodiacal, airglow, etc.) at a science pointing from simultaneous measurements at two sky arms (SKY_NEAR, SKY_FAR).

**Key Design Principles:**
- **Shared score encoders** for the two sky arms enforce arm-swap symmetry
- **Isolated branches** for anisotropic groups (zodi, continuum) bypass the shared trunk
- **Per-group heads** for six coefficient families: moon, mesospheric, ionospheric, atomic, continuum, zodi
- **Additive moon-zodi coupling** shares joint physics understanding between moon and zodi heads

### Model Components

#### 1. Sky Coefficient Grouping (6 groups)

| group | count | patterns | dominant physics |
|-------|------:|----------|------------------|
| `mesospheric` | ~403 | `OH_\d{3}`, `O2_b\d+` | thermospheric temp, geomagnetic activity |
| `moon` | 15 | `Moon_bs\d+` | $g_{\rm moon}(\phi, h, r_{\rm sep})$ (Krisciunas–Schaefer) |
| `zodi` | 5 | `Zodi_bs\d+` | $B_{500}(\lambda_\odot^{\rm rel}, \|\beta_{\rm ecl}\|)$ × airmass |
| `continuum` | 3 | HO2, FeO, O2Ac | mesospheric chemistry residual |
| `atomic` | 3 | ATOM_K | K I 7699 nightglow |
| `ionospheric` | 4 | ATOM_N | N I 5199 nightglow |

**Total:** 433 coefficients per row (deployed spline3-reduced basis, 2026-08-25)

#### 2. Network Architecture

**Step 1: Per-arm score encoder (shared weights)**
- Input: concatenated (score + context) from each sky arm
- Widths: `768 → 384` with LayerNorm + GELU
- Output: embeddings $\mathbf{e}^{\rm nr}$, $\mathbf{e}^{\rm fr}$

**Step 2: Context encoder (science-arm context only)**
- Input: science-arm context (39 features)
- Widths: `64`
- Output: embedding $\mathbf{e}^{\rm ctx}$

**Step 3: Symmetric fusion**
- Compute: mean, signed diff, absolute diff of arm embeddings
- Concatenate with context embedding
- Ensures arm-swap invariance

**Step 4: Shared trunk MLP**
- Input: fused representation (mean + diff + |diff| + ctx)
- Widths: `320 → 160`
- Output: 160-d latent $\mathbf{h}$

**Step 5: Per-group heads**
- **Trunk consumers** (moon, mesospheric, ionospheric, atomic): read $\mathbf{h}$ directly through `192 → $n_g^{\rm score}$` heads
- **Isolated zodi branch**: bypasses trunk, reads physics-motivated context subset (15 features × 3 arms + 2 arms × 5 scores = 55 dims), runs through `64 → 32 (GELU)` + zodi head
- **Isolated continuum branch**: similar bypass, reads 9 features × 3 arms + 2 arms × 3 scores = 33 dims, runs through `128 → 64 (GELU)` + continuum head
- **Moon-zodi coupling branch** (Phase F): small `94 → 64 → 32 (GELU)` branch with zero-initialized projectors to moon and zodi heads

#### 3. Isolated Zodi Branch (Phase B, 2026-08-22b)

**Motivation:** Zodiacal light is anisotropic on ecliptic plane; dependency on ecliptic geometry gets squeezed out by shared trunk.

**Implementation:**
- Bypasses shared trunk entirely
- Inputs: restricted context (15 features):
  - Zodi physics: `airmass`, `vanrhijn_285km`, `ecl_beta_deg`, `ecl_lon_{sin,cos}`, `zodi_log10_v`, `sun_sep`
  - Moon-driven contamination regime: `moon_alt`, `moon_sep`, `moon_phase_{sin,cos}`, `moon_fli`, `moon_up_smooth`, `moon_airmass_up`, `moon_signal_proxy`
- Two-layer MLP: `64 → 32 (GELU)` + head `(32,) → 5 scores`
- Rationale for moon proxies: On moon-down + low-ecliptic-latitude rows, split-zodi decomposition routes residual continuum into `Zodi_bs`; branch conditions on moon state to handle this

#### 4. Isolated Continuum Branch (Phase B/D)

**Motivation:** Continuum (HO₂/FeO/O₂Ac) has rf_R² = 0.70 driven by moon geometry; shared 160-d trunk is under-parameterized.

**Implementation:**
- Widths: `128 → 64 (GELU)` (widened in Phase D, 2026-08-26e)
- Context inputs (9 features): `moon_alt`, `moon_sep`, `moon_phase_{sin,cos}`, `moon_fli`, `airmass`, plus Phase-A' interactions
- Input dim: 33
- Head output: 3 scores

#### 5. Moon-Zodi Additive Coupling (Phase F, 2026-08-27c)

**Motivation:** Both moon and zodi benefit from knowing joint moon-scatter + zodi-geometry state.

**Implementation:**
- `moon_zodi_coupling_branch: 94 → 64 → 32 (GELU)`
- Zero-initialized linear projectors: $P_{\rm moon}(32 → 15)$, $P_{\rm zodi}(32 → 5)$
- Forward pass: $\hat{\mathbf{s}}_g = [\text{arm blend}] + \hat{\mathbf{s}}^{\rm head}_g + P_g(\text{CouplingBranch}(x_{\rm mz}))$
- Training starts byte-identical to pre-Phase-F; projectors learn what fraction of coupling to route to each head

#### 6. Per-Group Compression

Each group compressed via:
1. Robust per-column scale (median of positive training values)
2. Elementwise transform: `asinh` for mesospheric (many decades), `sqrt` for others
3. Standardize + PCA truncation on pooled training (near + far + sci)

**Deployed compression ratios:**
- `moon`: 15 → 15 (PCA rotation, all components retained)
- `zodi`: 5 → 5 (no PCA)
- `mesospheric`: 403 → 403 (no PCA)
- Others: 3-4 → 3-4 (no PCA)

---

## Geometric Normalization & Effective Extinction (Chapter 2)

### Framing Principle: Intrinsic-Emissivity Space

The ML model operates in **intrinsic-emissivity space**, not observed flux:

$$\mathbf{c}^{\rm intrinsic}_{r,g} = \frac{\mathbf{c}^{\rm obs}_{r,g}}{V(z_r; h_g) \cdot 10^{-0.4 k_g (X_r - 1)}}$$

Where:
- $V(z; h)$ = van Rhijn slant-path factor for shell height $h$
- $X$ = airmass; $k(\lambda)$ = extinction coefficient
- $z$ = zenith angle

**Why:** Intrinsic emissivity is the quantity that actually transfers between lines of sight. By dividing out geometry in closed form, the network focuses on remaining physics (gravity waves, short-timescale variability).

### Per-Group Layer Heights

| group | height | van Rhijn feature | physics |
|-------|-------:|-------------------|---------|
| mesospheric | 87 km | `vanrhijn_87km` | OH Meinel + O₂ bands |
| atomic | 95 km | `vanrhijn_95km` | K, Na, mesopause metals |
| ionospheric | 285 km | `vanrhijn_285km` | O I 6300/6364, F-region lines |
| continuum | 87 km | `vanrhijn_87km` | HO₂/FeO/O₂Ac chemiluminescence |
| moon | — | (factor 1.0) | scattering geometry, not thin-shell |
| zodi | — | (factor 1.0) | line-of-sight integrated dust |

### Effective Extinction ($k_{\rm eff}$)

Rather than model extinction, `fit_effective_extinction()` performs a **Bouguer fit using airglow as its own source**:

$$\ln\frac{A_{\rm near}}{A_{\rm far}} - \ln\frac{V_{\rm near}}{V_{\rm far}} = -0.4\ln 10 \; k_{\rm eff}(X_{\rm near} - X_{\rm far}) + \varepsilon$$

Where $\varepsilon$ = gravity-wave fluctuation (zero-mean, uncorrelated with Δ$X$).

**Recovery absorbs:** multiple-scattering correction, airglow-vs-stellar difference, site aerosol level, residual height errors.

---

## Cell 28 Diagnostic: Sky-Arm Zodi Bias by Moon-State Regime

### Diagnostic Output

```
Compares geometry-corrected sky1/sky2/mean and ML against sci truth on the zodi group.
n_rows = 1227 test rows

                        bias_B_mean_geo_%  bias_ML_default_%  delta (ML - mean)
all_test                       -6.10%              0.35%           +6.45% ⚠️
moon_phase Q1 [2-110 deg]     -27.86%              1.58%          +29.44% 🔴
moon_phase Q2 [110-172 deg]    -1.35%             -0.17%           +1.18% ✓
moon_phase Q3 [172-265 deg]    +0.74%              0.35%           -0.39% ✓
moon_phase Q4 [265-358 deg]   -36.00%              1.94%          +37.94% 🔴
moon_alt > 0 (moon up)         -1.74%              0.04%           +1.78% ⚠️
moon_alt <= 0 (moon down)      -48.34%             3.36%          +51.70% 🔴
airmass <= 1.5 (low)           -5.15%              0.66%           +5.81% ⚠️
airmass > 1.5 (high)           -8.94%             -0.59%           +8.35% ⚠️
```

### Key Findings

**Critical Issues (|delta| > 1.5%):**

1. **Moon-down rows** (`moon_alt ≤ 0`): 
   - Decomposition shows `B_mean = −48.34%` bias (poor geometric interpolation)
   - ML achieves `0.04%` but delta = `+51.70%` (NOT learning the residual)
   - **Problem:** Moon-down rows have residual continuum in `Zodi_bs`; model underfits
   - **Impact:** Highest variance regime; drives overall zodi bias

2. **Moon-phase Q1** [2–110°]:
   - High scattering phase; B_mean = −27.86%
   - ML bias = +1.58%, delta = +29.44%
   - Model not adapting to phase-dependent scattering geometry

3. **Moon-phase Q4** [265–358°]:
   - Another high-scattering phase; B_mean = −36.00%
   - ML bias = +1.94%, delta = +37.94%
   - Asymmetric failure suggests directional learning issue

**Acceptable Regimes (|delta| < 1.5%):**
- Q2 [110–172°]: delta = +1.18% (low-scattering phase)
- Q3 [172–265°]: delta = −0.39% (low-scattering phase, ML bias tracks decomp truth)

### Physical Interpretation

- **Moon down + low β_ecl:** Split-zodi decomposition cannot cleanly separate moon/zodi/continuum; routes residual into `Zodi_bs`. Model sees this as "zodi contamination" but isn't learning the regimemask.
- **High moon phase (Q1, Q4):** Krisciunas–Schaefer scattering has strong non-linear dependence on phase and zenith angle; current architecture underutilizes phase-dependent context.
- **Airmass > 1.5:** Similar issue—extinction × airmass product creates geometry mismatch; model doesn't adapt blend coefficient aggressively enough.

---

## Tuning Recommendations

### Priority 1: Strengthen Moon-Down Specialization ⭐ Highest ROI

**Problem:** Moon-down rows (51.7% delta) are the single largest source of zodi bias. The isolated zodi branch receives moon proxies but the shared trunk still dominates the pathway through the coupling branch.

**Actions:**

1. **Increase zodi branch hidden dimension** (matches Phase D continuum widening):
   - Current: `64 → 32 → 5`
   - Proposed: `128 → 64 → 5`
   - Rationale: Continuum widening won +0.7% RMS; zodi has similar under-parameterization

2. **Add explicit moon-down indicator to zodi branch:**
   - Create feature: `moon_down_flag = (moon_alt < -2).float()`
   - Or smoother: `moon_altitude_tanh = tanh(moon_alt / 5.0)` (continuous, differentiable)
   - Route directly into zodi branch input
   - Rationale: Gives branch explicit regime label; helps learned gates and attention mechanisms

3. **Increase coupling branch weight for moon-down rows:**
   - The `P_zodi` projector (32 → 5) is zero-initialized; it learns how much coupling latent to route
   - Consider: upweight `moon_down_rows` in training loss so projector learns stronger coupling signal on these high-bias rows
   - Implementation: `loss_zodi_weight[moon_down] *= 3.0` in training loop

4. **Add residual connection in zodi branch:**
   - Current: `zodi_scores = zodi_head(...)`
   - Proposed: Add skip: `zodi_scores = zodi_head(...) + 0.1 * linear_projection(zodi_input_subset)`
   - Rationale: Helps gradient flow on high-bias edge cases; learned skip weight can adjust

**Expected improvement:** Reduce moon-down delta from +51.7% to ~+10–15% (goal: decomposition truth-tracking or better)

---

### Priority 2: Moon-Phase Regime Stratification

**Problem:** Q1 (+29.4% delta) and Q4 (+37.9% delta) bright phases have strong scattering; model uses sinusoidal encoding which misses non-linear phase dependence.

**Actions:**

1. **Add quadratic phase features:**
   - Current: `moon_phase_sin`, `moon_phase_cos`
   - Add: `moon_phase_sin_sq = (moon_phase_sin)**2`, `moon_phase_cos_sq = (moon_phase_cos)**2`
   - Add: `moon_phase_magnitude = sqrt(moon_phase_sin**2 + moon_phase_cos**2)` (always 1.0, but its product with other features helps learning)
   - Route into: shared trunk (via sky encoder) + zodi branch + continuum branch

2. **Add moon-brightness-weighted phase interaction:**
   - Feature: `moon_fli_x_phase_cos = moon_fli * moon_phase_cos` (already deployed in continuum, add to zodi)
   - Rationale: Weights phase effect by actual scattered moonlight; low on phase Q2/Q3, high on Q1/Q4

3. **Consider phase-conditioned capacity (Phase G enhancement, optional):**
   - Learn scalar gate: `gate_phase = sigmoid(w_phase · [moon_phase_sin, moon_phase_cos, moon_fli])`
   - Modulate zodi branch hidden layer: `h_eff = gate_phase * h_default + (1 - gate_phase) * h_identity`
   - Rationale: Dynamically allocate capacity during bright phases; reverts to baseline on dim phases
   - Implementation complexity: moderate; may need Layer Norm adjustment

**Expected improvement:** Reduce Q1 delta from +29.4% to ~+8–12%, Q4 from +37.9% to ~+8–12%

---

### Priority 3: Residual Continuum Leakage on Moon-Down

**Problem:** Split-zodi decomposition routes HO₂/FeO/O₂Ac residual into `Zodi_bs` on moon-down rows; zodi and continuum heads compete for the same physical signal.

**Actions:**

1. **Ensure continuum branch sees moon-down regime:**
   - Continuum already uses `moon_alt`, `moon_sep`, `moon_fli`
   - Verify: `moon_down_flag` is routed into continuum branch input (currently not explicit)
   - Add feature to continuum input: `moon_down_flag` (same as Priority 1)

2. **Implement shared-loss regularization between zodi and continuum on moon-down:**
   - For moon-down rows, add term: $\lambda_{\rm reg} \|[\hat{c}_{\rm zodi}, \hat{c}_{\rm continuum}]\|_2^2$
   - Encourages model to minimize total continuum+zodi amplitude when moon is down
   - Rational: Helps model realize "decomposition ambiguity on moon-down; best to predict low"

3. **Add continuum-to-zodi information flow:**
   - Extend moon-zodi coupling to include continuum: `continuum_zodi_coupling_branch`
   - Input: `moon_down_flag` + (near_continuum + far_continuum) + near_zodi + far_zodi
   - Output: `coupling_residual_zodi`, `coupling_residual_continuum`
   - Zero-initialize; training learns whether coupling helps or not

**Expected improvement:** Reduce moon-down delta by additional 5–10% through better continuum-zodi partitioning

---

### Priority 4: Context-Alpha (Arm-Blend) Tuning

**Problem:** The blend coefficient $\alpha_g$ controls how much to weight near vs. far sky arm. On high-airmass or moon-down rows, $\alpha_{\rm zodi}$ may be stuck at uninformative defaults.

**Current formula (Phase F):**
$$\hat{\mathbf{s}}_{\rm zodi} = \alpha_{\rm zodi} \cdot \mathbf{s}^{\rm nr}_{\rm zodi} + (1-\alpha_{\rm zodi}) \cdot \mathbf{s}^{\rm fr}_{\rm zodi} + \hat{\mathbf{s}}^{\rm head}_{\rm zodi} + P_{\rm zodi}(\text{CouplingBranch})$$

Where $\alpha_{\rm zodi}$ is either scalar (learned constant) or per-row via: $\alpha_{\rm zodi} = \sigma(\text{MLP}_{\alpha}(x_{\rm sci}))$

**Actions:**

1. **Debug α predictions on problem regimes:**
   - Extract: `alpha_zodi[moon_alt <= 0]`, `alpha_zodi[phase_Q1]`, `alpha_zodi[airmass > 1.5]`
   - Check: Are they clustered near 0.5 (uninformative) or spread?
   - If clustered → α predictor is not learning regime-dependence → retrain with higher lr or wider MLP

2. **Add airmass-conditioned α prior:**
   - Anchor $\alpha_{\rm zodi}$ toward the sky arm at lower airmass (it's seeing less extinction-driven zodi signal)
   - Compute: `preferred_alpha = (airmass_far < airmass_near).float()` (0 if near closer to zenith, 1 if far closer)
   - Loss term: $\lambda_{\rm prior} \| \alpha_{\rm zodi} - \text{preferred}\_\alpha \|_2^2$ (soft prior, not hard constraint)
   - Rationale: Inductive bias grounded in physics; helps learning on sparse high-airmass data

3. **Widen context encoder for α prediction (if applicable):**
   - Current: `(39-dim ctx) → 64` (small)
   - Proposed: `(39-dim ctx) → 96 → 64` (adds one layer)
   - Rationale: Gives α predictor more capacity to distinguish high-airmass/moon-down regimes

**Expected improvement:** Tighter per-row α adaptation on problem regimes; +2–5% variance reduction

---

### Priority 5: Loss Function Stratification

**Problem:** Current loss may equally weight all regimes, diluting gradient signal for high-bias moon-down/bright-phase rows.

**Current approach (likely):** Uniform loss across all rows, with optional per-group Jensen lift (§3.6.6).

**Proposed stratification:**

```python
# Pseudo-code for training loop
loss_dict = {}

# Regime-weighted cross-entropy / MSE
is_moon_down = (moon_alt <= 0)
is_bright_phase = (moon_phase == 'Q1') | (moon_phase == 'Q4')
is_high_airmass = (airmass > 1.5)

# Base loss per row (RMS in intrinsic-emissivity space)
base_loss = mse_per_row(y_true, y_pred)

# Stratified weights
weight = ones_like(base_loss)
weight[is_moon_down] *= 3.0           # Highest priority
weight[is_bright_phase] *= 2.0        # High priority
weight[is_high_airmass] *= 1.5        # Medium priority

# Aggregate per group
loss_zodi_weighted = (base_loss[group_zodi] * weight[group_zodi]).mean()
loss_moon_weighted = (base_loss[group_moon] * weight[group_moon]).mean()
loss_other_weighted = (base_loss[group_other] * weight[group_other]).mean()

# Total loss with per-group balance
loss_total = (0.4 * loss_zodi_weighted + 
              0.3 * loss_moon_weighted + 
              0.3 * loss_other_weighted)
```

**Alternative: Jensen lift per regime:**
- Current Jensen lift applies a per-group lifting factor
- Proposed: Make lift regime-aware, e.g., `jensen_lift[moon_down] = 2.0`, `jensen_lift[nominal] = 1.0`

**Expected improvement:** +3–7% reduction in high-bias regime RMS

---

### Priority 6: Data Augmentation for Rare Regimes

**Problem:** Moon-down + low-β_ecl data is rare; model under-trains on these high-leverage regimes.

**Actions:**

1. **Oversampling (if using batch sampler):**
   - Increase `sample_weight` for moon-down rows to 3×
   - Ensure every epoch sees proportional moon-down representation
   - Implementation: `WeightedRandomSampler` with custom weights per slice

2. **Synthetic data via mixup (advanced):**
   - For moon-down rows: blend with nearby zenith rows via convex combination
   - Interpolate moon geometry: `moon_alt_interp = t * moon_alt_near_zenith + (1-t) * moon_alt_current`
   - Similarly for `airmass`, `ecl_beta`, etc.
   - Rationale: Expands effective training set; smooths geometry interpolation
   - Hyperparameter: mixup weight α ∈ [0.1, 0.3]

3. **Targeted fine-tuning phase:**
   - Train main model for N epochs on full dataset
   - Fine-tune for additional M epochs (e.g., M = 0.1 × N) on moon-down + Q1 + Q4 data only
   - Lower learning rate (0.1× main) to avoid catastrophic forgetting
   - Expected: +5–10% improvement on high-bias regimes; minimal degradation on nominal regimes

**Expected improvement:** +5–10% variance reduction on moon-down/bright-phase slices

---

### Priority 7: Architecture Refinement (Optional, if above insufficient)

**Only pursue if Priority 1–6 yield < 30% improvement on moon-down delta.**

**Option A: Deeper zodi branch**
- Current: `64 → 32 → 5`
- Proposed: `256 → 128 → 64 → 32 → 5`
- Add dropout=0.2 between layers to prevent overfitting
- Widened coupling branch to match: `128 → 64 → 32 (instead of 94 → 64 → 32)`
- Training: may need 2× learning rate schedule adjustment

**Option B: Residual connections in zodi branch**
- `out = zodi_head_main(...) + 0.1 * zodi_head_skip(...)`
- Skip pathway reads different subset of context to encourage feature orthogonality
- Helps gradient flow on moon-down edge cases

**Option C: Mixture-of-experts (MoE) zodi head**
- Instead of single zodi head, use K=3 expert heads weighted by gating network
- Experts: e.g., "moon-down expert", "bright-phase expert", "nominal expert"
- Gating: `gate = softmax(mlp(ctx))`, `out = sum(gate_i * expert_i(h))`
- Radically increases parameters (~5× zodi head cost); only if simpler strategies fail

**Expected improvement (if needed):** +10–20% additional delta reduction

---

## Implementation Priority Roadmap

### Phase 1 (Immediate, 1–2 days)
1. Add `moon_down_flag`, `moon_phase_sin_sq`, `moon_phase_cos_sq` features → zodi + continuum branches
2. Widen zodi branch: `64→32→5` → `128→64→5`
3. Reweight training loss: `moon_down_rows *= 3.0`
4. Retrain on validation set; measure delta improvement

### Phase 2 (If Phase 1 yields > 30% moon-down delta reduction, proceed; else add Priority 4–5)
1. Implement moon-fli×phase interaction features
2. Add residual connection in zodi branch
3. Debug α_zodi predictions per regime; adjust predictor width if needed
4. Retrain; measure cumulative improvement

### Phase 3 (If cumulative improvement < 50% delta reduction, add Priority 6)
1. Implement regime-weighted loss stratification
2. Oversample moon-down training data 3×
3. Retrain with fine-tuning phase on high-bias slices
4. Measure; iterate if needed

### Phase 4 (Optional, if needed)
1. Implement deeper zodi branch + dropout (Priority 7A)
2. Or mixture-of-experts head (Priority 7C)
3. Full retraining; measure end-to-end improvement

---

## Expected Outcomes

If Priorities 1–3 are implemented fully:
- **Moon-down delta:** +51.7% → +10–15% (goal: ~0%, tracking decomposition truth)
- **Q1 phase delta:** +29.4% → +8–12%
- **Q4 phase delta:** +37.9% → +8–12%
- **Overall test RMS:** likely −3–8% improvement in zodi RMS across all slices

If Priorities 1–6 are implemented:
- Target: **Moon-down delta → 0–5%** (converges to nominal performance)
- Target: **Bright-phase delta → 0–5%**
- Tradeoff: nominal-regime performance may degrade by ~1–2% (mitigated by early stopping on validation zodi_nominal)

---

## Quick Debugging Checklist

Before/after implementation, verify:

- [ ] `alpha_zodi` predictions are regime-adaptive (not stuck at 0.5 for moon-down)
- [ ] Zodi branch is receiving moon-down-flag in input tensor
- [ ] Training loss is higher on moon-down rows initially (sanity check on reweighting)
- [ ] Validation loss on `moon_alt <= 0` slice improves after Priority 1 implementation
- [ ] No train-val divergence on nominal slices (risk of overfitting to moon-down)
- [ ] Model checkpoint improves holdout test RMS on moon-down slice

---

## References

**Notebook sections:**
- §1.2.2: Zodi carrier (5 basis functions, split family)
- §2.6: Moon and zodi exempt from van Rhijn (scattering vs. thin-shell)
- §3.2: Coefficient grouping
- §3.3: Per-group compressor
- §3.4: Architecture (dual-encoder, group heads, isolated branches)
- §3.4.1: Isolated zodi branch (Phase B)
- §3.4.3: Additive moon-zodi coupling (Phase F)

**Diagnostic cell:** Cell 28 (`sky_arm_zodi_bias()`)

---

**Generated:** 2026-09-01
