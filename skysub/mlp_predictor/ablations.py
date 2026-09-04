"""Ablation configs for the post-reversal-fix simplification pass (2026-09-04).

The decomposition corpus was regenerated with the moon/zodi identifiability
priors (see ``moon_zodi_split_priors_2026-09-03.md``).  Several ML components
were added while the decomposition was reversing the moon and zodi roles on
~98% of moon-up rows, and may have been compensating for that.  Each entry
below is a config override to A/B against the deployed default.

Usage in the training-config cell::

    from mlp_predictor import ablations
    train_cfg.update(ablations.ABLATIONS["A1_no_coupling"])   # or omit for the default
    print(ablations.describe("A1_no_coupling"))

Run one ablation per training run, 10 seeds, same split.

READING THE RESULT -- the aggregate alone is not enough
-------------------------------------------------------
On this corpus the seed-to-seed ``mean_eRMSE`` std is **0.331** and the
10-seed ensemble stderr is **0.105** (~1.3% of the 8.289 mean).  So a change
in the aggregate below ~0.33 is seed noise, not signal.

Judge each ablation on the group-level diagnostics, which is where a
reversal-compensating component would show up:

* ``naive_baseline`` (cell 22) per-group ML gain -- the cross-corpus-safe
  metric.  Deployed reference: moon +23.0%, zodi +41.0%, continuum +10.6%,
  mesospheric +3.9%, ionospheric +5.6%, atomic +20.1%; group-equal +5.5%.
* ``sky_arm_zodi_bias`` (cell 27) -- deployed reference: 0 slices carry an
  ML-specific zodi bias.  Any ablation that puts a slice back over threshold
  is a real regression.
* ``wavelength_residual_atlas`` (cell 39) per-band RMS|frac| -- deployed
  reference: blue 0.59%, mid 1.11%, NIR 0.90%.
* ``sky_arm_disagreement_floor`` (cell 38) err/delta -- compare WITHIN this
  corpus only.  Do not compare against the pre-fix numbers in the changelog:
  the near-far moon disagreement fell 12x when the reversal was removed
  (p50 delta 0.0554 -> 0.0045), so the denominator changed and the ratio is
  not commensurable across corpora.
"""

from __future__ import annotations

ABLATIONS: dict[str, dict] = {
    # ---------------------------------------------------------------- A1
    "A1_no_coupling": {
        "moon_zodi_coupling_enabled": False,
    },
    # ---------------------------------------------------------------- A2
    "A2_moon_weight_down": {
        "moon_group_weight": 1.5,
    },
    # ---------------------------------------------------------------- A3
    "A3_no_flux_mse": {
        "flux_mse_groups": (),
    },
    # ---------------------------------------------------------------- A4
    "A4_thin_zodi_head": {
        "zodi_head_extra_dims": (),
    },
    # ---------------------------------------------------------------- B1/B2
    "B1_flux_moon_only": {
        "flux_mse_groups": ("moon",),
    },
    "B2_flux_zodi_only": {
        "flux_mse_groups": ("zodi",),
    },
    # ---------------------------------------------------------------- A5
    "A5_no_coupling_thin_zodi": {
        "moon_zodi_coupling_enabled": False,
        "zodi_head_extra_dims": (),
    },
}

#: Config keys the deployed notebook still sets but the trainer no longer
#: reads -- it warns about these on every run.  Listed here because two of them
#: are active hazards for this simplification pass: ``moon_zodi_mode`` looks
#: like the coupling's off-switch and is not (hence A1 needed a real one), and
#: ``alpha_ctx_groups`` looks like it restricts ctx-alpha and does not.  An A/B
#: driven through any of these measures nothing at all.
DEAD_CONFIG_KEYS: tuple[str, ...] = (
    "alpha_ctx_groups", "block_cov_loss_groups", "flux_mse_eps_frac",
    "head_extra_dims", "high_airmass_boost", "moon_down_ecliptic_beta_deg",
    "moon_down_ecliptic_boost", "moon_zodi_branch_dims", "moon_zodi_mode",
    "moon_zodi_moon_head_extra_dims", "moon_zodi_zodi_head_extra_dims",
    "relative_mse_eps_frac", "relative_mse_groups", "use_coef_err_weights",
)

_NOTES: dict[str, str] = {
    "A1_no_coupling": """
A1 -- remove the additive moon-zodi coupling (~9k params).  THE PRIME SUSPECT.

  A shared 94->64->32 latent is projected additively into BOTH the moon and the
  zodi head via zero-init linear projectors.  Architecturally that is exactly a
  mechanism for shifting flux between the two families in a correlated way,
  i.e. for undoing a role swap.

  What the trained weights already say: the projectors did NOT stay at zero --
  |W|_F is 0.94 (moon) and 0.25 (zodi), about 10.6% of the moon head's own
  weight norm.  So the coupling is doing something.  But the loading is
  ASYMMETRIC, 3.7x more into moon than zodi, whereas pure reversal
  compensation should be roughly equal-and-opposite.  That is more consistent
  with the stated Phase-F intent -- "give the moon head the restricted
  moon-scatter context it otherwise never sees" -- which would survive the
  corpus fix.  Weight norms cannot settle it; this run can.

  Note the off-switch had to be restored: the no-coupling path was deleted in
  the 2026-08-27 trim, so the changelog's "set the moon_zodi_* entries to None"
  raised a ValueError.  With the same seed the disabled model is now
  BIT-IDENTICAL at initialisation to the enabled one (verified), so the A/B
  isolates the coupling and not a different init.

  Expect if it was reversal compensation: no material change, or an
  improvement.  Expect if it was the context fix: moon gain drops from +23%.
""",
    "A2_moon_weight_down": """
A2 -- lower the moon group weight 3.0 -> 1.5.  THE NEW PROBLEM, not an old one.

  The priors made the moon target strongly heavy-tailed.  Median moon_tot is
  146 against p99 = 1.2e5 (a factor of 830), and the population median moon
  SHARE of the moon+zodi continuum is 0.020 -- the prior floor.  Over half the
  corpus now has essentially no moon, and all the moon signal is concentrated
  in a minority of rows.  Before the fix the moon was spread across every row
  (median coefficient sum 2.14 vs 0.153 now).

  An absolute MSE on such a target is dominated by the few bright-moon rows,
  which is the likely reason new large outliers appeared in the moon panel
  (rows 1183 and 1222 are both moon-dominated, moon share 0.69-0.96).  The
  weight m_moon = 3.0 was tuned when the moon was spread over every row; on a
  target this concentrated the same weight puts proportionally more of the
  gradient on the tail and on the shared trunk.

  Per-row relative MSE would be the better instrument, but ``relative_mse_groups``
  was REMOVED from the trainer in the 2026-08-27 trim -- it is one of the 14
  dead config keys the Trainer now warns about.  Reinstating it is a code
  change, not an ablation; this weight sweep is the config-only probe.

  If lowering the weight improves the moon panel tail without costing the
  moon's +23% baseline gain, the heavy tail is the issue and reinstating
  relative MSE is worth the code change.
""",
    "A3_no_flux_mse": """
A3 -- drop the flux-space MSE.  MEASURED: catastrophic.  Do not revisit.

  RESULT (10 seeds): moon gain over B0_copy_near +24.0% -> +5.5%, zodi
  +41.3% -> +23.2%, and on moon-up rows zodi went NEGATIVE (-11.2%, i.e. worse
  than copying the near arm).  Aggregate mean_eRMSE +0.537 and seed std tripled
  (0.276 -> 0.806).  Continuum -7.1pp and mesospheric -4.5pp too, groups that
  have no flux term at all -- so the effect propagates through the shared
  trunk.  The term is a primary driver of accuracy, not overhead.

  The hypothesis this ablation was built on was WRONG, and the correction
  matters for anyone planning further work here: the term does NOT constrain
  the moon+zodi SUM, so it was never "insurance that survives a role swap".
  Reading ``compressed_loss``: it is strictly per-group and it REPLACES that
  group's coefficient-space smooth_l1 loss (if/else, not an added term).  Each
  group's own coefficients are inverted through its own compressor and basis to
  per-pixel flux and scored there.  A swap would not be invisible to it.

  Consequence for the planned full-spectrum flux loss with its own lambda:
  that is a genuinely different design, not a tuning knob.  Today the flux term
  SUPPLANTS the coefficient loss for moon and zodi; an additive lambda-weighted
  term would keep both signals, and would need the per-group weights re-tuned
  because the current group_loss_weight values were fitted against a loss where
  moon and zodi are scored in flux units and everything else in coefficient
  units.
""",
    "A4_thin_zodi_head": """
A4 -- remove the extra 32-d hidden layer from the isolated zodi head.

  The zodi head is now comfortably the best-behaved group: err/delta = 0.354
  all-regime (the diagnostic's own oracle floor is ~1.41, so it is far inside
  the irreducible band), a +41% gain over the best naive baseline, and zero
  ML-specific bias on all nine slices of the zodi-bias diagnostic.

  A group sitting that far below its floor does not need more capacity -- extra
  capacity there chases irreducible noise and drags the shared trunk, which is
  the caution the disagreement-floor cell prints.  This is the cheap partial
  test.  Removing the isolated zodi BRANCH entirely would need a code switch
  like A1's; do that only if A4 shows the capacity is genuinely unused.
""",
    "B1_flux_moon_only": """
B1/B2 -- give the flux-space loss to ONE group at a time.

  Follow-up to A3, which showed the flux-space term is load-bearing: removing
  it entirely cost moon 18.5pp and zodi 18.1pp of their gain over
  B0_copy_near, took moon-up zodi NEGATIVE (-11.2%, worse than copying the
  near arm), and tripled the seed std (0.28 -> 0.81).

  Note what the term actually is, because it is easy to misread: it is
  strictly PER-GROUP and it REPLACES that group's coefficient-space
  smooth_l1 loss (an if/else in ``compressed_loss``, not an added term).  Each
  group's own coefficients are pushed through its own basis matrix to per-pixel
  flux and scored there.  No moon+zodi SUM is constrained anywhere -- so the
  term is not "insurance on the total", and asking whether the gain comes from
  the sum or from per-coefficient weighting is not a well-posed question.

  What the split does test is SEPARABILITY.  If the mechanism is purely
  per-group, B1 should hold the moon near the default +24% while zodi falls
  back toward A3's +23%, and B2 should mirror that.  If instead both groups
  degrade in both arms, the effect is routed through the shared trunk -- which
  A3 already hints at, since it also cost mesospheric 4.5pp and continuum
  7.1pp, groups that have no flux term at all.

  Only ``moon`` and ``zodi`` have basis matrices
  (``_precompute_flux_basis_and_geometry``), so adding other groups to
  ``flux_mse_groups`` silently skips them -- supplying more bases is part of
  the full-spectrum flux-loss work, not a config change.
""",
    "A5_no_coupling_thin_zodi": """
A5 -- A1 + A4 together.

  Run only after A1 and A4 have each been measured alone.  Both touch the
  moon/zodi pair, so a combined run cannot attribute a change to either, but if
  both are individually neutral this confirms they are neutral together rather
  than compensating for one another.
""",
}


def describe(name: str) -> str:
    """Rationale, prior evidence, and what to watch for one ablation."""
    if name not in ABLATIONS:
        raise KeyError(f"unknown ablation {name!r}; have {sorted(ABLATIONS)}")
    return _NOTES.get(name, "(no notes)").strip("\n")


def apply(train_cfg: dict, name: str, verbose: bool = True) -> dict:
    """Apply one ablation's overrides in place and report what changed."""
    if name not in ABLATIONS:
        raise KeyError(f"unknown ablation {name!r}; have {sorted(ABLATIONS)}")
    overrides = ABLATIONS[name]
    # Validate against the trainer's own consumed-key set, not just against
    # train_cfg: a key can be present in the config and still be dead (14 of
    # them currently are).  Setting a dead knob would silently produce a null
    # A/B, which is exactly the failure this pass exists to avoid.
    from . import trainer as _trainer
    live = getattr(_trainer.Trainer, "_CONSUMED_CFG_KEYS", None)
    for key in overrides:
        if key not in train_cfg:
            raise KeyError(
                f"ablation {name!r} sets {key!r}, which is not in train_cfg; "
                "the config schema has changed and this ablation needs updating."
            )
        if live is not None and key not in live:
            raise KeyError(
                f"ablation {name!r} sets {key!r}, which the trainer does NOT read "
                f"(see ablations.DEAD_CONFIG_KEYS). The run would be identical to "
                f"the default -- fix the ablation or plumb the knob through."
            )
    if verbose:
        print(f"=== ablation {name} ===")
        for key, new in overrides.items():
            print(f"  {key}: {train_cfg[key]!r} -> {new!r}")
    train_cfg.update(overrides)
    return train_cfg


__all__ = ["ABLATIONS", "describe", "apply"]
