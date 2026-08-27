"""Ensemble trainer + inference for the compressed dual-encoder group-head MLP.

Extracted from cell ``18659296`` of the split-zodi notebook:
- ``DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP`` : per-group sigma floors
  (fraction of per-column median finite sigma) used by the coef-err
  weighted training loss.
- ``train_compressed_group_mlp`` : trains one seed of the compressed
  dual-encoder + group-head model.  Big function (~1150 lines) kept in
  one piece so its many implicit invariants stay adjacent to the code
  that maintains them.
- ``predict_sci_coefficients_default`` : ensemble-mean inference in
  compressed-score space with clip-to-upper-bound.
- ``default_dual_group_config`` : the deployed 10-seed config.
- ``Trainer`` : thin orchestration class wrapping the notebook's
  ensemble loop (config resolution, flux-MSE basis matrix precompute,
  seed loop, per-seed test metric report).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from astropy.io import fits

from sky_decomp.lsf_surface_iterative import SkyDecompLSFSurfaceIterative

from .compressor import (
    compress_coef_err_to_score_sigma,
    compress_coefs_to_scores,
    expand_scores_to_coefs,
)
from .data import _infer_base_dir_for_reconstruction, airglow_geometry_scale
from .metrics import metric_row
from .ml_utils import (
    RobustScaler,
    moon_phase_deg_from_ctx,
    set_reproducibility,
    split_indices_by_moon_phase,
)
from .model import DualEncoderGroupHeadMLPCompressed

# Aliases so the extracted function body works verbatim without dropping the leading underscore.
_moon_phase_deg_from_ctx = moon_phase_deg_from_ctx
_set_reproducibility = set_reproducibility

DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP = {
    'zodi': 0.05,
    # Per-group floors on the relative sigma (fraction of the per-column
    # median finite sigma).  These were set on 2026-08-16 from the sigma /
    # weight diagnostic (§7): moon / continuum / ionospheric have
    # well-behaved tails and take the historical 5% floor; mesospheric
    # (403 lines) and atomic have p99 / p50 ratios of 10^3-10^6 in
    # compressed-score space, so they need a 20% floor to keep
    # w_p99 / w_median below ~30 per column.
    'moon': 0.05,
    'continuum': 0.05,
    'mesospheric': 0.20,
    'ionospheric': 0.05,
    'atomic': 0.20,
}


def train_compressed_group_mlp(
    filtered, compressors, group_indices, geom_kwargs,
    split_indices=None,
    n_epochs=50, batch_size=256, lr=7e-4,
    encoder_dims=(384, 192), ctx_dims=(64,),
    trunk_dims=(320, 160), head_dim=192,
    head_extra_dims=(),
    zodi_head_extra_dims=(),
    weight_decay=1e-4, grad_clip=1.0, patience=4,
    seed=42, drop_vanrhijn_from_context=False,
    moon_group_weight=1.0,
    zodi_group_weight=1.0,
    continuum_group_weight=1.0,
    mesospheric_group_weight=1.0,
    ionospheric_group_weight=1.0,
    blend_init_alpha=0.7,
    blend_optim='direct',
    moon_alt_conditional_alpha=False,
    zodi_ctx_restriction=None,
    continuum_ctx_restriction=None,
    continuum_head_extra_dims=(),
    continuum_branch_dims=(64, 32),
    alpha_ctx_features=None,
    alpha_ctx_groups=None,
    # 2026-08-25: opt-in covariance-aware loss for moon and/or zodi.  Empty
    # tuple keeps the existing per-element weighted Huber loss (backwards
    # compat).  Setting e.g. ('moon',) precomputes Sigma_scaled^-1 per row
    # from filtered_triplet['coef_cov_moon_sci'] (persisted by fit.py) and
    # replaces the diagonal-weighted loss on that group with the block
    # Mahalanobis form L_g = r^T Sigma^-1 r / n_g, scaled to match the
    # diagonal loss on train rows at init so downstream group weights
    # (moon_group_weight, ...) keep their meaning.  Requires the group's
    # compressor to be per-element (no PCA); moon and zodi both are.
    # 2026-08-25b: default flipped to () after relative_mse_groups
    # became the preferred loss-magnitude fix for moon/zodi (see below).
    block_cov_loss_groups=(),
    # 2026-08-25b: opt-in relative-error MSE per group. Per group g in this
    # tuple the loss becomes L_g(row) = ||pred_raw - true_raw||^2 /
    # (||true_raw||^2 + eps_g^2) on the RobustScaler-undone score vector,
    # where eps_g^2 = (relative_mse_eps_frac)^2 * median_train||true_raw||^2.
    # Takes precedence over block_cov_loss_groups for the same group.
    # Rationale: COEF_ERR sigma is ~30x too large so 1/sigma^2 weighting
    # misweights low- vs high-flux rows; relative MSE couples the loss to
    # target amplitude directly.
    relative_mse_groups=(),
    relative_mse_eps_frac=0.05,
    # 2026-08-25c: opt-in per-pixel flux MSE for moon and/or zodi.  Takes
    # precedence over relative_mse_groups and block_cov_loss_groups.
    # ``flux_basis_matrices`` is a dict {group: (n_coef_g, n_wave_ds)}
    # numpy array carrying the pre-LSF spline basis matrix (from
    # ``SkyDecompLSFSurfaceIterative.matrix_moon`` / ``matrix_zodi``).
    # ``flux_geom_sc_sci`` is (n_row, n_coef) with the per-row per-coef
    # geometry scale from ``airglow_geometry_scale(ctx_sci, ...)``.  The
    # torch inverse compressor runs sqrt-inverse + per-group PCA-inverse
    # + destandardise + geometry, mirroring ``inverse_group_compressor``.
    flux_mse_groups=(),
    flux_mse_eps_frac=0.0,
    flux_basis_matrices=None,
    flux_geom_sc_sci=None,
    high_airmass_boost=1.0,
    # 2026-08-21: per-regime row-weight boosts for the two rare-but-important
    # sky regimes identified by the row-741-vs-830 diagnostic (see §12).  All
    # default to 1.0 (disabled); flip in default_dual_group_config to enable.
    moon_down_ecliptic_boost=1.0,
    moon_down_ecliptic_beta_deg=15.0,
    bright_moon_close_boost=1.0,
    bright_moon_close_fli_min=0.85,
    bright_moon_close_sep_max_deg=45.0,
    # 2026-08-21b: physics-informed continuous zodi-asymmetry boost.
    # Multiplier = 1 + K * moon_down_gate(moon_alt) * max(0, log10(zodi_v_sci/zodi_v_far)).
    # 0.0 = disabled (additive convention).
    # New (2026-08-16): heteroscedastic-Gaussian loss weighting from the
    # decomposition COEF_ERR (§7 of the top methods cell).  The floor is
    # now per group; the mesospheric and atomic groups need a larger floor
    # because their p99/p50 sigma ratio is 10^3-10^6 (§7 diagnostic).
    use_coef_err_weights=True,
    coef_err_sigma_floor_rel=None,
    # 2026-08-24: per-group row-ceiling on sigma. sigma_eff = min(sigma,
    # beta_g * col_median_sigma[j]) per column j. beta_g dict per group
    # (or scalar for all).  None = disabled (default). Breaks the
    # sigma ~ |coef| scaling which under-weights bright-moon rows.
    coef_err_row_ceiling_factor=None,
):
    _set_reproducibility(seed)

    ctx_near = np.asarray(filtered['ctx_near'], dtype=np.float32)
    ctx_far = np.asarray(filtered['ctx_far'], dtype=np.float32)
    ctx_sci = np.asarray(filtered['ctx_sci'], dtype=np.float32)

    if split_indices is not None:
        train_idx, val_idx, test_idx = split_indices
    else:
        _moon_phase_col = _moon_phase_deg_from_ctx(filtered)
        train_idx, val_idx, test_idx = split_indices_by_moon_phase(
            filtered['obstime_mjd'], _moon_phase_col, seed=seed)
    train_idx = np.asarray(train_idx, dtype=int)
    val_idx = np.asarray(val_idx, dtype=int)
    test_idx = np.asarray(test_idx, dtype=int)

    coef_near = np.asarray(filtered['coef_near'], dtype=np.float64)
    coef_far = np.asarray(filtered['coef_far'], dtype=np.float64)
    coef_sci = np.asarray(filtered['coef_sci'], dtype=np.float64)

    scores_near, slices_near = compress_coefs_to_scores(
        coef_near, ctx_near, compressors, geom_kwargs, group_indices)
    scores_far, slices_far = compress_coefs_to_scores(
        coef_far, ctx_far, compressors, geom_kwargs, group_indices)
    scores_sci, slices_sci = compress_coefs_to_scores(
        coef_sci, ctx_sci, compressors, geom_kwargs, group_indices)
    assert slices_near == slices_far == slices_sci, 'score slices disagree'
    score_slices = slices_near
    n_input_score = scores_near.shape[1]
    group_score_dims = {g: (hi - lo) for g, (lo, hi) in score_slices.items()}

    print(f'Compressed input dimension: {n_input_score} '
          f'(uncompressed was {coef_near.shape[1]}); per-group scores: '
          + ', '.join(f'{g}={n}' for g, n in group_score_dims.items()))

    # Per-group multiplier m_g on top of 1/sqrt(n_g_score); defaults to 1.0.
    _group_multipliers = {'moon': float(moon_group_weight),
                          'zodi': float(zodi_group_weight),
                          'continuum': float(continuum_group_weight),
                          'mesospheric': float(mesospheric_group_weight),
                          'ionospheric': float(ionospheric_group_weight)}
    group_loss_weight = {
        g: (1.0 / max(float(np.sqrt(max(int(n), 1))), 1.0))
           * _group_multipliers.get(g, 1.0)
        for g, n in group_score_dims.items()
    }
    print(f'Loss weights (moon_group_weight={float(moon_group_weight):.2f}, '
          f'zodi_group_weight={float(zodi_group_weight):.2f}, '
          f'continuum_group_weight={float(continuum_group_weight):.2f}, '
          f'mesospheric_group_weight={float(mesospheric_group_weight):.2f}, '
          f'ionospheric_group_weight={float(ionospheric_group_weight):.2f}):')
    print(f"  {'group':<14s} {'n_score':>7s} {'m_g':>6s} {'w_g':>8s} {'w_g/n':>10s}")
    for _g, _n in group_score_dims.items():
        _m = _group_multipliers.get(_g, 1.0)
        _w = group_loss_weight[_g]
        print(f'  {_g:<14s} {int(_n):>7d} {_m:>6.2f} {_w:>8.4f} {(_w / max(int(_n), 1)):>10.5f}')

    score_scaler = RobustScaler().fit(np.vstack([
        scores_near[train_idx], scores_far[train_idx], scores_sci[train_idx],
    ]).astype(np.float32))
    ctx_scaler = RobustScaler().fit(np.vstack([
        ctx_near[train_idx], ctx_far[train_idx], ctx_sci[train_idx],
    ]))

    near_s = np.clip(score_scaler.transform(scores_near.astype(np.float32)),
                     -25.0, 25.0).astype(np.float32)
    far_s = np.clip(score_scaler.transform(scores_far.astype(np.float32)),
                    -25.0, 25.0).astype(np.float32)
    sci_s = np.clip(score_scaler.transform(scores_sci.astype(np.float32)),
                    -25.0, 25.0).astype(np.float32)

    # -----------------------------------------------------------------
    # Per-element inverse-variance weights from decomposition COEF_ERR.
    # Propagate through the compressor with a first-order Jacobian, then
    # divide out the score_scaler scale so the weights live in the same
    # scaled-score space as `sci_s`; finally clip to a per-column floor
    # (frac of median finite sigma) and normalise so E_train[w_pe] = 1
    # per column.  Missing/boundary sigmas fall through to the floor.
    # See §7 of the top methods cell for the derivation.
    coef_err_sci_local = np.asarray(
        filtered.get('coef_err_sci', np.full_like(coef_sci, np.nan)),
        dtype=np.float64,
    )
    _n_finite_err = int(np.sum(np.isfinite(coef_err_sci_local)))
    _n_total_err = int(coef_err_sci_local.size)
    _weights_enabled = bool(use_coef_err_weights) and _n_finite_err > 0
    if _weights_enabled:
        _sigma_scores_sci, _sigma_slices_sci = compress_coef_err_to_score_sigma(
            coef_sci, coef_err_sci_local, ctx_sci,
            compressors, geom_kwargs, group_indices,
        )
        assert _sigma_slices_sci == score_slices, 'sigma slices disagree with score slices'
        _scale_ = np.asarray(getattr(score_scaler, 'scale_',
                                    np.ones(n_input_score)),
                             dtype=np.float64)
        _scale_ = np.where(_scale_ > 0.0, _scale_, 1.0)
        _sigma_scaled = _sigma_scores_sci / _scale_[None, :]
        _finite = np.where(np.isfinite(_sigma_scaled) & (_sigma_scaled > 0.0),
                           _sigma_scaled, np.nan)
        _floor_col = np.nanmedian(_finite, axis=0)
        _floor_col = np.where(np.isfinite(_floor_col) & (_floor_col > 0.0),
                              _floor_col, 1.0)
        # Resolve per-group floors: accept scalar (broadcast), dict
        # (per-group), or None (built-in DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP).
        _floor_arg = coef_err_sigma_floor_rel
        if _floor_arg is None:
            _floor_arg = DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP
        if isinstance(_floor_arg, dict):
            _floor_by_group = {str(k): float(v) for k, v in _floor_arg.items()}
        else:
            _floor_by_group = {g: float(_floor_arg) for g in score_slices}
        _floor_rel_col = np.ones(_floor_col.shape[0], dtype=np.float64)
        _missing_groups = []
        for _g, (_lo, _hi) in score_slices.items():
            if _g not in _floor_by_group:
                _missing_groups.append(_g)
                _fr = 0.05
            else:
                _fr = _floor_by_group[_g]
            _floor_rel_col[_lo:_hi] = _fr
        if _missing_groups:
            print(f'  warning: coef_err_sigma_floor_rel missing groups '
                  f'{_missing_groups}; using fallback 0.05 for each.')
        _floor_col = _floor_rel_col * _floor_col
        _sigma_scaled = np.where(_sigma_scaled > 0.0, _sigma_scaled,
                                 _floor_col[None, :])
        _sigma_scaled = np.maximum(_sigma_scaled, _floor_col[None, :])
        if coef_err_row_ceiling_factor is not None:
            if isinstance(coef_err_row_ceiling_factor, dict):
                _ceil_by_group = {str(k): float(v)
                                  for k, v in coef_err_row_ceiling_factor.items()}
            else:
                _ceil_by_group = {g: float(coef_err_row_ceiling_factor)
                                  for g in score_slices}
            _col_med_sigma = np.nanmedian(_sigma_scaled[train_idx], axis=0)
            _col_med_sigma = np.where(np.isfinite(_col_med_sigma)
                                      & (_col_med_sigma > 0.0),
                                      _col_med_sigma, 1.0)
            _col_ceil = np.full(_sigma_scaled.shape[1], np.inf, dtype=np.float64)
            for _g, (_lo, _hi) in score_slices.items():
                _bg = _ceil_by_group.get(_g)
                if _bg is not None and np.isfinite(_bg) and _bg > 0.0:
                    _col_ceil[_lo:_hi] = float(_bg) * _col_med_sigma[_lo:_hi]
            _n_capped = int(np.sum(_sigma_scaled[train_idx]
                                    > _col_ceil[None, :]))
            _sigma_scaled = np.minimum(_sigma_scaled, _col_ceil[None, :])
            _ceil_report = ', '.join(
                f'{g}=x{_ceil_by_group.get(g, np.inf):.2g}'
                for g in score_slices)
            print(f'Per-group sigma row-ceiling: {_ceil_report}')
            print(f'  fraction of train coef-cells capped: '
                  f'{_n_capped / max(_sigma_scaled[train_idx].size, 1):.3%}')
        w_pe_np = 1.0 / (_sigma_scaled ** 2)
        _wcm = np.mean(w_pe_np[train_idx], axis=0)
        _wcm = np.where(_wcm > 0.0, _wcm, 1.0)
        w_pe_np = (w_pe_np / _wcm[None, :]).astype(np.float32)
        _floor_report = ', '.join(f'{g}={_floor_by_group.get(g, 0.05):.3g}'
                                  for g in score_slices)
        print(
            f'Per-element loss weights from coef_err_sci: '
            f'finite fraction={_n_finite_err / max(_n_total_err, 1):.3f} in native space, '
            f'compressed w median={float(np.median(w_pe_np[train_idx])):.3f}, '
            f'1-99% = [{float(np.percentile(w_pe_np[train_idx], 1)):.2f}, '
            f'{float(np.percentile(w_pe_np[train_idx], 99)):.2f}]'
        )
        print(f'  per-group floor_rel: {_floor_report}')
        # Save the resolved dict so the artifact/config can round-trip it.
        _resolved_floor_by_group = {g: float(_floor_by_group.get(g, 0.05))
                                    for g in score_slices}
    else:
        # Uniform per-element weights: reduces to the previous unweighted loss.
        w_pe_np = np.ones_like(sci_s, dtype=np.float32)
        print('Per-element loss weights: disabled '
              '(use_coef_err_weights=False or no finite COEF_ERR); '
              'compressed_loss falls back to the unweighted SmoothL1.')

    ctx_near_n = np.clip(ctx_scaler.transform(ctx_near), -25.0, 25.0).astype(np.float32)
    ctx_far_n = np.clip(ctx_scaler.transform(ctx_far), -25.0, 25.0).astype(np.float32)
    ctx_sci_n = np.clip(ctx_scaler.transform(ctx_sci), -25.0, 25.0).astype(np.float32)

    # Row weights: boost rare-but-important regimes so training sees them
    # proportional to deployment relevance rather than raw frequency.  Boosts
    # multiply and the final vector is normalized so the mean training weight
    # = 1.  See §12 (2026-08-12) for the earlier ablation that rejected
    # high-F10.7 boosts; 2026-08-21 added moon_down_ecliptic + bright_moon_close
    # boosts for the two failure modes identified in the row-741-vs-830
    # diagnostic (moon-down on-ecliptic zodi red slope + bright-moon
    # extreme-blue).  All boosts default to 1.0 (disabled).
    _ctx_names_local = list(filtered['ctx_names'])
    _row_weights_np = np.ones(ctx_sci.shape[0], dtype=np.float32)

    def _apply_row_boost(_mask, _factor, _label):
        nonlocal _row_weights_np
        _mask = np.asarray(_mask, dtype=bool)
        _factor = float(_factor)
        _n_all = int(_mask.sum())
        _n_tr = int(_mask[train_idx].sum())
        _n_va = int(_mask[val_idx].sum())
        if _factor <= 1.0 or _n_all == 0:
            if _factor > 1.0:
                print(f'Row weights: {_label} requested (x{_factor:.2f}) but '
                      f'no rows match; skipping.')
            return
        _row_weights_np = np.where(_mask, _row_weights_np * _factor,
                                    _row_weights_np).astype(np.float32)
        print(f'Row weights: {_label} x{_factor:.2f}  '
              f'n_all={_n_all}, n_train={_n_tr}, n_val={_n_va}')

    if 'airmass' in _ctx_names_local:
        _ai = _ctx_names_local.index('airmass')
        _apply_row_boost(ctx_sci[:, _ai] > 1.5, high_airmass_boost, 'airmass>1.5')

    _has_ma = 'moon_alt' in _ctx_names_local
    _has_eb = 'ecl_beta_deg' in _ctx_names_local
    if _has_ma and _has_eb:
        _ma = ctx_sci[:, _ctx_names_local.index('moon_alt')]
        _eb = ctx_sci[:, _ctx_names_local.index('ecl_beta_deg')]
        _mask_mde = (_ma < 0.0) & (np.abs(_eb) < float(moon_down_ecliptic_beta_deg))
        _apply_row_boost(_mask_mde, moon_down_ecliptic_boost,
                         f'moon_down + |ecl_beta|<{float(moon_down_ecliptic_beta_deg):.0f}')
    elif float(moon_down_ecliptic_boost) > 1.0:
        print('Row weights: moon_down_ecliptic boost requested but '
              'moon_alt or ecl_beta_deg missing from ctx; skipping.')

    _has_ms = 'moon_sep' in _ctx_names_local
    _has_pc = ('moon_phase_sin' in _ctx_names_local
               and 'moon_phase_cos' in _ctx_names_local)
    if _has_ma and _has_ms and _has_pc:
        _ma = ctx_sci[:, _ctx_names_local.index('moon_alt')]
        _ms = ctx_sci[:, _ctx_names_local.index('moon_sep')]
        _pc = ctx_sci[:, _ctx_names_local.index('moon_phase_cos')]
        # moon_fli = (1 - cos(phase))/2 (fill-illumination convention).
        _fli = np.clip(0.5 * (1.0 - _pc), 0.0, 1.0)
        _mask_bmc = ((_ma > 0.0)
                     & (_fli >= float(bright_moon_close_fli_min))
                     & (_ms <= float(bright_moon_close_sep_max_deg)))
        _apply_row_boost(_mask_bmc, bright_moon_close_boost,
                         f'bright_moon (fli>={float(bright_moon_close_fli_min):.2f}, '
                         f'sep<={float(bright_moon_close_sep_max_deg):.0f}, alt>0)')
    elif float(bright_moon_close_boost) > 1.0:
        print('Row weights: bright_moon_close boost requested but '
              'moon_alt/moon_sep/moon_phase_{sin,cos} missing from ctx; skipping.')

    _train_mean_w = float(np.mean(_row_weights_np[train_idx]))
    if _train_mean_w > 1e-9 and abs(_train_mean_w - 1.0) > 1e-6:
        _row_weights_np = (_row_weights_np / _train_mean_w).astype(np.float32)
        print(f'Row weights: normalized to train-mean=1 '
              f'(pre-norm train-mean={_train_mean_w:.3f}, '
              f'post-norm max={float(np.max(_row_weights_np[train_idx])):.3f}).')

    # ------------------------------------------------------------------
    # Block-Mahalanobis loss precomputation.
    # For each configured block group g, build Sigma_scaled_inv[r] per row
    # by propagating the persisted native COEF_COV block through the
    # diagonal per-element Jacobian implied by _sigma_scaled / sigma_native.
    # ------------------------------------------------------------------
    _block_cov_inv_by_group = {}
    _block_cov_scale_by_group = {}
    if block_cov_loss_groups:
        _sigma_native_sci_bc = np.abs(
            np.asarray(filtered['coef_err_sci'], dtype=np.float64))
        for _g_bc in block_cov_loss_groups:
            _key_bc = f'coef_cov_{_g_bc}_sci'
            _cov_native_bc = filtered.get(_key_bc)
            if _cov_native_bc is None:
                print(f'  [block-cov-loss] {_g_bc}: no {_key_bc} in filtered; '
                      f'skipping (diagonal Huber will apply).')
                continue
            _gidx_bc = np.asarray(group_indices[_g_bc], dtype=int)
            _lo_bc, _hi_bc = score_slices[_g_bc]
            _n_g_bc = int(_hi_bc - _lo_bc)
            if _n_g_bc != _gidx_bc.size:
                print(f'  [block-cov-loss] {_g_bc}: PCA-compressed group not '
                      f'supported (n_score={_n_g_bc} != n_native='
                      f'{_gidx_bc.size}); skipping.')
                continue
            _n_row_bc = int(coef_sci.shape[0])
            _cov_inv_g = np.zeros((_n_row_bc, _n_g_bc, _n_g_bc), dtype=np.float32)
            _sigma_scaled_g = _sigma_scaled[:, _lo_bc:_hi_bc]
            _sigma_native_g = _sigma_native_sci_bc[:, _gidx_bc]
            _J_g = _sigma_scaled_g / np.maximum(_sigma_native_g, 1e-30)
            _n_active_row = np.zeros(_n_row_bc, dtype=np.int32)
            _n_success = 0
            _quad_train = []
            for _r_bc in range(_n_row_bc):
                _cov_native_r = np.asarray(_cov_native_bc[_r_bc], dtype=np.float64)
                _diag_native = np.diagonal(_cov_native_r)
                _active_bc = np.isfinite(_diag_native) & (_diag_native > 0.0)
                if _active_bc.sum() < 2:
                    _n_active_row[_r_bc] = 0
                    continue
                _J_r = _J_g[_r_bc]
                _cov_scaled_full = (_J_r[:, None] * _J_r[None, :]) * _cov_native_r
                _idx_active = np.flatnonzero(_active_bc)
                _sub = _cov_scaled_full[np.ix_(_idx_active, _idx_active)]
                _sub = 0.5 * (_sub + _sub.T) + 1e-12 * np.eye(_sub.shape[0])
                try:
                    _inv_sub = np.linalg.pinv(_sub)
                except np.linalg.LinAlgError:
                    continue
                _cov_inv_g[_r_bc][np.ix_(_idx_active, _idx_active)] = (
                    _inv_sub.astype(np.float32))
                _n_active_row[_r_bc] = int(_active_bc.sum())
                _n_success += 1
                # Baseline quadratic form on truth vs zero (for scale match):
                # for a well-calibrated Sigma, E[c^T Sigma^-1 c] = n_active
                # under H0.  On real data we use the median TRUE quadratic
                # over train rows to match the diagonal loss magnitude.
                if _r_bc < 8 or (train_idx is not None and _r_bc in train_idx):
                    _true_r = _sigma_scaled_g[_r_bc]
                    # Placeholder; actual scale matching done after loop.
            # Scale-matching: on train rows, compute
            #   median( c^T Sigma^-1 c / n_active )   where c = truth in scaled score.
            # We rescale Sigma^-1 so this median equals the median of the
            # per-row diagonal loss  mean_k( w_pe_g * (c_k)^2 )  on train rows,
            # so switching to Mahalanobis leaves the overall scale of the
            # group's contribution roughly unchanged.
            _true_scaled = sci_s[:, _lo_bc:_hi_bc]  # score-space TRUE (scaled)
            _train_mask_bc = np.zeros(_n_row_bc, dtype=bool)
            _train_mask_bc[np.asarray(train_idx, dtype=int)] = True
            _train_mask_bc &= (_n_active_row > 1)
            _idx_tr_bc = np.flatnonzero(_train_mask_bc)
            _tr_maha = np.empty(_idx_tr_bc.size, dtype=np.float64)
            for _i_bc, _r_bc in enumerate(_idx_tr_bc):
                _t = _true_scaled[_r_bc]
                _M = _cov_inv_g[_r_bc]
                _tr_maha[_i_bc] = float(_t @ (_M @ _t)) / max(int(_n_active_row[_r_bc]), 1)
            _tr_diag = np.mean(
                w_pe_np[_idx_tr_bc, _lo_bc:_hi_bc]
                * (_true_scaled[_idx_tr_bc] ** 2), axis=1)
            _med_maha = float(np.nanmedian(_tr_maha)) if _tr_maha.size else 1.0
            _med_diag = float(np.nanmedian(_tr_diag)) if _tr_diag.size else 1.0
            if _med_maha > 0 and _med_diag > 0:
                _scale = _med_diag / _med_maha
            else:
                _scale = 1.0
            _cov_inv_g *= np.float32(_scale)
            _block_cov_inv_by_group[_g_bc] = _cov_inv_g
            _block_cov_scale_by_group[_g_bc] = float(_scale)
            print(f'  [block-cov-loss] {_g_bc}: precomputed Sigma_scaled^-1 for '
                  f'{_n_success}/{_n_row_bc} rows ({_n_g_bc}x{_n_g_bc} per row); '
                  f'median n_active={int(np.median(_n_active_row[_n_active_row>0]))}, '
                  f'scale-match factor={_scale:.4g}  '
                  f'(median diag={_med_diag:.4g}, median Mahalanobis={_med_maha:.4g}).')

    # ------------------------------------------------------------------
    # 2026-08-25b: relative-MSE precomputation for moon / zodi.
    # Undo the RobustScaler per group so the loss lives in unscaled score
    # space:
    #   L_g(row) = ||pred_raw - true_raw||^2 / (||true_raw||^2 + eps_g^2)
    # where eps_g^2 = (relative_mse_eps_frac)^2 * median_train||true_raw||^2.
    # A scale-match factor brings the per-row relative loss onto the same
    # magnitude as the diagonal Huber term it replaces (median truth-vs-zero
    # convention, matching block_cov_loss_groups above).
    # ------------------------------------------------------------------
    _rel_mse_state_by_group = {}
    if relative_mse_groups:
        _center_ = np.asarray(getattr(score_scaler, 'center_',
                                      getattr(score_scaler, 'mean_',
                                              np.zeros(n_input_score))),
                              dtype=np.float64)
        for _g_rm in relative_mse_groups:
            if _g_rm not in score_slices:
                print(f'  [relative-mse] {_g_rm}: not in score_slices; skipping.')
                continue
            _lo_rm, _hi_rm = score_slices[_g_rm]
            _sc_rm = _scale_[_lo_rm:_hi_rm]
            _mu_rm = _center_[_lo_rm:_hi_rm]
            _t_scaled_tr = sci_s[train_idx, _lo_rm:_hi_rm].astype(np.float64)
            _t_raw_tr = _t_scaled_tr * _sc_rm[None, :] + _mu_rm[None, :]
            _norm2_tr = np.nansum(_t_raw_tr ** 2, axis=1)
            _n2_ok = _norm2_tr[np.isfinite(_norm2_tr) & (_norm2_tr > 0.0)]
            _median_norm2 = float(np.median(_n2_ok)) if _n2_ok.size else 1.0
            _eps_sq_rm = ((float(relative_mse_eps_frac) ** 2)
                          * max(_median_norm2, 1e-30))
            # Scale-match on truth-vs-zero (matches block-cov convention).
            _tr_rel = _norm2_tr / (_norm2_tr + _eps_sq_rm)
            _tr_diag_rm = np.nanmean(
                w_pe_np[train_idx, _lo_rm:_hi_rm]
                * (sci_s[train_idx, _lo_rm:_hi_rm].astype(np.float64) ** 2),
                axis=1)
            _med_rel = float(np.nanmedian(_tr_rel)) if _tr_rel.size else 1.0
            _med_diag_rm = (float(np.nanmedian(_tr_diag_rm))
                            if _tr_diag_rm.size else 1.0)
            if _med_rel > 0 and _med_diag_rm > 0:
                _scale_rm = _med_diag_rm / _med_rel
            else:
                _scale_rm = 1.0
            _rel_mse_state_by_group[_g_rm] = {
                'scale': _sc_rm.astype(np.float32),
                'center': _mu_rm.astype(np.float32),
                'eps_sq': float(_eps_sq_rm),
                'scale_match': float(_scale_rm),
            }
            print(f'  [relative-mse] {_g_rm}: n_score={_hi_rm - _lo_rm}, '
                  f'median ||target_raw||^2={_median_norm2:.4g}, '
                  f'eps^2={_eps_sq_rm:.4g} '
                  f'(eps_frac={float(relative_mse_eps_frac):.3f}), '
                  f'scale-match factor={_scale_rm:.4g}  '
                  f'(median diag={_med_diag_rm:.4g}, median rel={_med_rel:.4g}).')

    # ------------------------------------------------------------------
    # 2026-08-25c: per-pixel flux MSE precomputation for moon / zodi.
    # For each group g, undo the compressor + geometry to recover native
    # coefs, then A_g @ c gives flux.  Loss = mean_lambda((f_pred - f_true)^2)
    # per row, scale-matched to the diagonal Huber magnitude on train rows.
    # Assumes n_kept == n_full for the group (holds for moon+zodi in this
    # notebook — see compressor cell); asserts otherwise.
    # ------------------------------------------------------------------
    _flux_mse_state_by_group = {}
    if flux_mse_groups and flux_basis_matrices is not None and flux_geom_sc_sci is not None:
        _center_ = np.asarray(getattr(score_scaler, 'center_',
                                      getattr(score_scaler, 'mean_',
                                              np.zeros(n_input_score))),
                              dtype=np.float64)
        for _g_fm in flux_mse_groups:
            if _g_fm not in score_slices:
                print(f'  [flux-mse] {_g_fm}: not in score_slices; skipping.')
                continue
            _A_np = flux_basis_matrices.get(_g_fm)
            if _A_np is None:
                print(f'  [flux-mse] {_g_fm}: no basis matrix supplied; skipping.')
                continue
            _comp_fm = compressors[_g_fm]
            _gidx_fm = np.asarray(group_indices[_g_fm], dtype=int)
            _lo_fm, _hi_fm = score_slices[_g_fm]
            _n_kept_fm = int(_hi_fm - _lo_fm)
            _n_full_fm = int(_comp_fm['basis'].shape[0])
            if _n_kept_fm != _n_full_fm:
                print(f'  [flux-mse] {_g_fm}: kept={_n_kept_fm} != full={_n_full_fm}; '
                      f'PCA-subset flux MSE not implemented, skipping.')
                continue
            if str(_comp_fm.get('kind')) != 'sqrt':
                print(f'  [flux-mse] {_g_fm}: only sqrt compressor is supported '
                      f'(got {_comp_fm.get("kind")}); skipping.')
                continue
            _A_g = np.asarray(_A_np, dtype=np.float32)
            if _A_g.shape[0] != _gidx_fm.size:
                print(f'  [flux-mse] {_g_fm}: basis rows {_A_g.shape[0]} != '
                      f'n_coef {_gidx_fm.size}; skipping.')
                continue
            _sc_fm = _scale_[_lo_fm:_hi_fm]
            _mu_fm = _center_[_lo_fm:_hi_fm]
            _sc_sci_g_np = np.asarray(flux_geom_sc_sci[:, _gidx_fm], dtype=np.float32)
            # Scale-match on train rows: median flux^2 vs median diagonal Huber.
            _c_true_train = coef_sci[train_idx][:, _gidx_fm].astype(np.float64)
            _flux_true_train = _c_true_train @ _A_g.astype(np.float64)
            _row_flux_norm_train = np.mean(_flux_true_train ** 2, axis=1)
            _fin_row = np.isfinite(_row_flux_norm_train) & (_row_flux_norm_train > 0.0)
            _median_flux_norm = (float(np.median(_row_flux_norm_train[_fin_row]))
                                 if _fin_row.any() else 1.0)
            _tr_diag_fm = np.nanmean(
                w_pe_np[train_idx, _lo_fm:_hi_fm]
                * (sci_s[train_idx, _lo_fm:_hi_fm].astype(np.float64) ** 2),
                axis=1)
            _med_diag_fm = (float(np.nanmedian(_tr_diag_fm))
                            if _tr_diag_fm.size else 1.0)
            _scale_match_fm = (_med_diag_fm / max(_median_flux_norm, 1e-30)
                               if _med_diag_fm > 0 else 1.0)
            _flux_mse_state_by_group[_g_fm] = {
                'A': _A_g,
                'basis_T': _comp_fm['basis'].T.astype(np.float32),
                'sd_vec': _comp_fm['sd_vec'].astype(np.float32),
                'mean_vec': _comp_fm['mean_vec'].astype(np.float32),
                'score_scale': _sc_fm.astype(np.float32),
                'score_center': _mu_fm.astype(np.float32),
                'g_scale_sci': _sc_sci_g_np,
                'scale_match': float(_scale_match_fm),
                'median_flux_norm': float(_median_flux_norm),
            }
            print(f'  [flux-mse] {_g_fm}: n_coef={_gidx_fm.size}, '
                  f'n_wave_ds={_A_g.shape[1]}, '
                  f'median mean(flux^2)={_median_flux_norm:.3g}, '
                  f'scale-match={_scale_match_fm:.4g} '
                  f'(median diag={_med_diag_fm:.4g}).')

    # Moon_alt sign per row (physical ctx_sci). 0 = moon down (alt<=0), 1 = up.
    _moon_alt_idx_local = _ctx_names_local.index('moon_alt')
    _moon_alt_sign_all = (ctx_sci[:, _moon_alt_idx_local] > 0.0).astype(np.int64)

    if torch.cuda.is_available():
        device = 'cuda'
    elif getattr(torch.backends, 'mps', None) is not None and torch.backends.mps.is_available():
        device = 'mps'
    else:
        device = 'cpu'

    # 2026-08-24: skip DataLoader for the ~O(10MB) in-memory training set.
    # Move all tensors to device once and iterate with torch.randperm --
    # avoids ~180k per-batch CPU->device transfers per full 10-seed run
    # (~2x wall-time speedup on MPS).
    # Deterministic order for the extra block-cov-inv tensors staged
    # in the tuple below (kept as a module-level for the batch unpack).
    _block_cov_group_order = tuple(
        g for g in block_cov_loss_groups if g in _block_cov_inv_by_group)
    # 2026-08-25b: relative-MSE state on device (float32 buffers per group).
    _rel_mse_torch_by_group = {}
    for _g_rm, _st_rm in _rel_mse_state_by_group.items():
        _rel_mse_torch_by_group[_g_rm] = {
            'scale': torch.from_numpy(_st_rm['scale']).to(device),
            'center': torch.from_numpy(_st_rm['center']).to(device),
            'eps_sq': float(_st_rm['eps_sq']),
            'scale_match': float(_st_rm['scale_match']),
        }
    # 2026-08-25c: per-pixel flux MSE state on device (float32 buffers per group).
    _flux_mse_torch_by_group = {}
    for _g_fm, _st_fm in _flux_mse_state_by_group.items():
        _flux_mse_torch_by_group[_g_fm] = {
            'A':           torch.from_numpy(_st_fm['A']).to(device),
            'basis_T':     torch.from_numpy(_st_fm['basis_T']).to(device),
            'sd_vec':      torch.from_numpy(_st_fm['sd_vec']).to(device),
            'mean_vec':    torch.from_numpy(_st_fm['mean_vec']).to(device),
            'score_scale':  torch.from_numpy(_st_fm['score_scale']).to(device),
            'score_center': torch.from_numpy(_st_fm['score_center']).to(device),
            'g_scale_sci':  torch.from_numpy(_st_fm['g_scale_sci']).to(device),
            'scale_match':  float(_st_fm['scale_match']),
        }
    _row_idx_all_np = np.arange(int(coef_sci.shape[0]), dtype=np.int64)
    def _stage_on_device(_idx):
        _base = (
            torch.from_numpy(near_s[_idx]).to(device),
            torch.from_numpy(far_s[_idx]).to(device),
            torch.from_numpy(ctx_near_n[_idx]).to(device),
            torch.from_numpy(ctx_far_n[_idx]).to(device),
            torch.from_numpy(ctx_sci_n[_idx]).to(device),
            torch.from_numpy(_moon_alt_sign_all[_idx]).to(device),
            torch.from_numpy(_row_weights_np[_idx]).to(device),
            torch.from_numpy(w_pe_np[_idx]).to(device),
            torch.from_numpy(sci_s[_idx]).to(device),
            torch.from_numpy(_row_idx_all_np[_idx]).to(device),
        )
        _cov = tuple(
            torch.from_numpy(_block_cov_inv_by_group[_g][_idx]).to(device)
            for _g in _block_cov_group_order)
        return _base + _cov
    _tr_tensors = _stage_on_device(train_idx)
    _va_tensors = _stage_on_device(val_idx)
    _tr_bs = int(batch_size)
    _va_bs = 512
    _n_train = _tr_tensors[0].shape[0]
    _n_val = _va_tensors[0].shape[0]

    model = DualEncoderGroupHeadMLPCompressed(
        n_score=n_input_score, n_ctx=ctx_near_n.shape[1],
        group_score_dims=group_score_dims,
        ctx_names=[str(x) for x in filtered['ctx_names']],
        encoder_dims=tuple(int(v) for v in encoder_dims),
        ctx_dims=tuple(int(v) for v in ctx_dims),
        trunk_dims=tuple(int(v) for v in trunk_dims),
        head_dim=int(head_dim),
        head_extra_dims=tuple(int(v) for v in head_extra_dims),
        zodi_head_extra_dims=tuple(int(v) for v in zodi_head_extra_dims),
        continuum_head_extra_dims=tuple(int(v) for v in continuum_head_extra_dims),
        continuum_branch_dims=tuple(int(v) for v in continuum_branch_dims),
        drop_vanrhijn_from_context=bool(drop_vanrhijn_from_context),
        blend_init_alpha=float(blend_init_alpha),
        blend_use_direct=(str(blend_optim).lower()
                          in ('direct', 'prefit_freeze', 'prefit_warmstart')),
        moon_alt_conditional_alpha=bool(moon_alt_conditional_alpha),
        alpha_ctx_features=alpha_ctx_features,
        alpha_ctx_groups=alpha_ctx_groups,
        zodi_ctx_restriction=zodi_ctx_restriction,
        continuum_ctx_restriction=continuum_ctx_restriction,
    ).to(device)

    # Optimizer construction depends on blend_optim mode:
    #   'default' -> single AdamW param group.
    #   'no_wd'   -> blend params in a group with weight_decay=0.
    #   'lr10x'   -> blend params in a group with lr=10*base_lr.
    #   'direct'  -> alpha_g stored directly (no sigmoid); blend group has weight_decay=0
    #                and is clamped to [eps, 1-eps] after each step.
    _blend_mode = str(blend_optim).lower()
    _ALLOWED_BLEND_MODES = ('default', 'no_wd', 'lr10x', 'direct',
                            'prefit_freeze', 'prefit_warmstart')
    if _blend_mode not in _ALLOWED_BLEND_MODES:
        raise ValueError(f"Unknown blend_optim={blend_optim!r}; "
                         f"expected one of {_ALLOWED_BLEND_MODES}.")
    # Pre-fit alpha_g in scaled score space on training rows. Closed-form OLS:
    #   alpha_g = <s_near - s_far, s_sci - s_far> / ||s_near - s_far||^2
    # summed across all coefficients in group g. 'prefit_freeze' locks the
    # value; 'prefit_warmstart' uses it as an init and lets AdamW adjust.
    _prefit_alphas = None
    if _blend_mode in ('prefit_freeze', 'prefit_warmstart'):
        _prefit_alphas = {}
        for _g_pf, (_lo_pf, _hi_pf) in score_slices.items():
            if model.moon_alt_conditional_alpha and _g_pf == 'moon':
                continue  # conditional moon alpha does not use the closed-form prefit.
            _sn = near_s[train_idx, _lo_pf:_hi_pf].astype(np.float64)
            _sf = far_s[train_idx, _lo_pf:_hi_pf].astype(np.float64)
            _ss = sci_s[train_idx, _lo_pf:_hi_pf].astype(np.float64)
            _diff = _sn - _sf
            _target = _ss - _sf
            _num = float(np.sum(_diff * _target))
            _den = float(np.sum(_diff * _diff))
            if _den <= 1e-12:
                # near ≈ far for this group; blend is ill-defined -> midpoint.
                _alpha_star = 0.5
            else:
                _alpha_star = float(np.clip(_num / _den,
                                            float(model.blend_alpha_eps),
                                            1.0 - float(model.blend_alpha_eps)))
            _prefit_alphas[_g_pf] = _alpha_star
            with torch.no_grad():
                model.blend_alpha_direct[str(_g_pf)].data.fill_(float(_alpha_star))
        if _blend_mode == 'prefit_freeze':
            for _g_pf in group_score_dims:
                model.blend_alpha_direct[str(_g_pf)].requires_grad_(False)
        print(f'Pre-fit alpha per group ({_blend_mode}): '
              + '  '.join(f'{_g_pf}={_prefit_alphas[_g_pf]:.3f}'
                          for _g_pf in group_score_dims))

    blend_pnames = set()
    blend_params = []
    if model.blend_use_direct:
        for _k, _p in model.blend_alpha_direct.items():
            blend_params.append(_p)
            blend_pnames.add(f'blend_alpha_direct.{_k}')
        if model.moon_alt_conditional_alpha and hasattr(model, 'blend_alpha_moon_by_alt'):
            blend_params.append(model.blend_alpha_moon_by_alt)
            blend_pnames.add('blend_alpha_moon_by_alt')
    else:
        for _k, _p in model.blend_logit.items():
            blend_params.append(_p)
            blend_pnames.add(f'blend_logit.{_k}')
        if model.moon_alt_conditional_alpha and hasattr(model, 'blend_moon_by_alt_logit'):
            blend_params.append(model.blend_moon_by_alt_logit)
            blend_pnames.add('blend_moon_by_alt_logit')
    other_params = [p for n, p in model.named_parameters() if n not in blend_pnames]
    if _blend_mode == 'default':
        opt = torch.optim.AdamW(model.parameters(), lr=float(lr),
                                weight_decay=float(weight_decay))
    elif _blend_mode == 'no_wd':
        opt = torch.optim.AdamW(
            [{'params': other_params, 'weight_decay': float(weight_decay)},
             {'params': blend_params, 'weight_decay': 0.0}],
            lr=float(lr))
    elif _blend_mode == 'lr10x':
        opt = torch.optim.AdamW(
            [{'params': other_params, 'lr': float(lr)},
             {'params': blend_params, 'lr': float(lr) * 10.0}],
            weight_decay=float(weight_decay))
    elif _blend_mode == 'prefit_freeze':
        # blend params carry the OLS alpha and are frozen (requires_grad=False).
        opt = torch.optim.AdamW(other_params, lr=float(lr),
                                weight_decay=float(weight_decay))
    else:  # 'direct' or 'prefit_warmstart'
        opt = torch.optim.AdamW(
            [{'params': other_params, 'weight_decay': float(weight_decay)},
             {'params': blend_params, 'weight_decay': 0.0}],
            lr=float(lr))
    print(f"Blend optim mode: {_blend_mode} "
          f"(use_direct={bool(model.blend_use_direct)}, "
          f"n_blend_params={sum(p.numel() for p in blend_params)})")

    def compressed_loss(pred_dict, yb, w_row, w_pe, row_idx_b=None, block_cov_inv_batch=None):
        loss = torch.tensor(0.0, device=yb.device)
        for g, y_head in pred_dict.items():
            lo, hi = score_slices[g]
            target = yb[:, lo:hi].contiguous()
            res = y_head - target
            if g in _flux_mse_torch_by_group and row_idx_b is not None:
                # 2026-08-25c: per-pixel flux MSE for moon and/or zodi.
                # Inverse compressor in torch: undo RobustScaler, PCA, sqrt,
                # geometry.  Reconstruct flux via c @ A and MSE per pixel.
                _st_fm = _flux_mse_torch_by_group[g]
                _raw_target = target * _st_fm['score_scale'] + _st_fm['score_center']
                _raw_pred   = y_head * _st_fm['score_scale'] + _st_fm['score_center']
                _cent_true = _raw_target @ _st_fm['basis_T']
                _cent_pred = _raw_pred   @ _st_fm['basis_T']
                _z_true = _cent_true * _st_fm['sd_vec'] + _st_fm['mean_vec']
                _z_pred = _cent_pred * _st_fm['sd_vec'] + _st_fm['mean_vec']
                _em_true = torch.clamp(_z_true, min=0.0) ** 2
                _em_pred = torch.clamp(_z_pred, min=0.0) ** 2
                _g_row = _st_fm['g_scale_sci'][row_idx_b, :]
                _c_true = _em_true * _g_row
                _c_pred = _em_pred * _g_row
                _flux_true = _c_true @ _st_fm['A']
                _flux_pred = _c_pred @ _st_fm['A']
                _per_row = ((_flux_pred - _flux_true) ** 2).mean(dim=1) * float(_st_fm['scale_match'])
            elif g in _rel_mse_torch_by_group:
                # 2026-08-25b: relative-error MSE on the RobustScaler-undone
                # score vector.  L_g(row) = ||pred_raw - true_raw||^2 /
                # (||true_raw||^2 + eps_g^2), then scaled to match the
                # diagonal Huber magnitude on train rows so
                # group_loss_weight[g] stays calibrated across loss types.
                # Takes precedence over block_cov_inv_batch.
                _st_rm = _rel_mse_torch_by_group[g]
                _sc_rm = _st_rm['scale']
                _mu_rm = _st_rm['center']
                target_raw = target * _sc_rm + _mu_rm
                res_raw = res * _sc_rm
                _num = (res_raw * res_raw).sum(dim=1)
                _den = (target_raw * target_raw).sum(dim=1) + _st_rm['eps_sq']
                _per_row = (_num / _den) * float(_st_rm['scale_match'])
            elif block_cov_inv_batch is not None and g in block_cov_inv_batch:
                # Block Mahalanobis in scaled score space.  The precomputed
                # Sigma^-1 already includes the scale-match factor so this
                # term is on the same magnitude as the diagonal Huber term
                # it replaces.
                _Sinv = block_cov_inv_batch[g]
                _per_row = torch.einsum('bi,bij,bj->b', res, _Sinv, res)
                _per_row = _per_row / float(max(int(_Sinv.shape[-1]), 1))
            else:
                w_pe_g = w_pe[:, lo:hi].contiguous()
                _per_elem = F.smooth_l1_loss(y_head, target, reduction='none') * w_pe_g
                _per_row = _per_elem.mean(dim=1)
            _weighted = (w_row * _per_row).mean()
            _w_g = group_loss_weight[g]
            loss = loss + float(_w_g) * _weighted
        return loss / max(len(pred_dict), 1)

    def _snapshot_alpha():
        result = {}
        for g in group_score_dims:
            if model.moon_alt_conditional_alpha and g == 'moon':
                if model.blend_use_direct:
                    _a = model.blend_alpha_moon_by_alt.detach().cpu().numpy()
                else:
                    _a = torch.sigmoid(model.blend_moon_by_alt_logit).detach().cpu().numpy()
                result['moon_dn'] = float(_a[0])
                result['moon_up'] = float(_a[1])
            elif model.blend_use_direct:
                result[g] = float(model.blend_alpha_direct[str(g)].item())
            else:
                result[g] = float(torch.sigmoid(model.blend_logit[str(g)]).item())
        return result
    history = []
    blend_history = []  # per-epoch snapshot of learned alpha_g
    _init_blend = _snapshot_alpha()
    blend_history.append({'epoch': 0, **_init_blend})
    best_val = np.inf
    best_epoch = -1
    best_state = None
    stale = 0
    for ep in range(1, int(n_epochs) + 1):
        model.train()
        tr_loss = 0.0
        tr_n = 0
        _perm = torch.randperm(_n_train, device=device)
        for _b0 in range(0, _n_train, _tr_bs):
            _sel = _perm[_b0:_b0 + _tr_bs]
            _batch_tr = tuple(t[_sel] for t in _tr_tensors)
            near_b, far_b, near_ctx_b, far_ctx_b, ctx_b, ma_sign_b, w_row_b, w_pe_b, yb, row_idx_b = _batch_tr[:10]
            _block_cov_batch_tr = {
                _g: _batch_tr[10 + _i]
                for _i, _g in enumerate(_block_cov_group_order)
            }
            pred = model(near_b, far_b, near_ctx_b, far_ctx_b, ctx_b,
                         sci_moon_alt_sign=(ma_sign_b
                                            if model.moon_alt_conditional_alpha else None))
            loss = compressed_loss(pred, yb, w_row_b, w_pe_b,
                                    row_idx_b=row_idx_b,
                                    block_cov_inv_batch=(_block_cov_batch_tr
                                                         if _block_cov_batch_tr else None))
            if not torch.isfinite(loss):
                continue
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), float(grad_clip))
            opt.step()
            if model.blend_use_direct:
                with torch.no_grad():
                    _eps = float(model.blend_alpha_eps)
                    for _g in group_score_dims:
                        if model.moon_alt_conditional_alpha and _g == 'moon':
                            continue
                        model.blend_alpha_direct[str(_g)].clamp_(_eps, 1.0 - _eps)
                    if model.moon_alt_conditional_alpha and hasattr(model, 'blend_alpha_moon_by_alt'):
                        model.blend_alpha_moon_by_alt.clamp_(_eps, 1.0 - _eps)
            tr_loss += float(loss.item())
            tr_n += 1

        model.eval()
        va_loss = 0.0
        va_n = 0
        with torch.no_grad():
            for _b0 in range(0, _n_val, _va_bs):
                _sel = slice(_b0, _b0 + _va_bs)
                _batch_va = tuple(t[_sel] for t in _va_tensors)
                near_b, far_b, near_ctx_b, far_ctx_b, ctx_b, ma_sign_b, w_row_b, w_pe_b, yb, row_idx_b = _batch_va[:10]
                _block_cov_batch_va = {
                    _g: _batch_va[10 + _i]
                    for _i, _g in enumerate(_block_cov_group_order)
                }
                pred = model(near_b, far_b, near_ctx_b, far_ctx_b, ctx_b,
                             sci_moon_alt_sign=(ma_sign_b
                                                if model.moon_alt_conditional_alpha else None))
                loss = compressed_loss(pred, yb, w_row_b, w_pe_b,
                                        row_idx_b=row_idx_b,
                                        block_cov_inv_batch=(_block_cov_batch_va
                                                             if _block_cov_batch_va else None))
                if torch.isfinite(loss):
                    va_loss += float(loss.item())
                    va_n += 1

        tr_mean = tr_loss / max(tr_n, 1)
        va_mean = va_loss / max(va_n, 1)
        history.append({'epoch': ep, 'train_loss': tr_mean, 'val_loss': va_mean})
        blend_history.append({'epoch': ep, **_snapshot_alpha()})
        if va_mean < best_val:
            best_val = va_mean
            best_epoch = ep
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if ep == 1 or ep % 5 == 0 or ep == int(n_epochs):
            print(f'[compressed] epoch={ep:03d} train={tr_mean:.6f} val={va_mean:.6f}')
        if stale >= int(patience):
            print(f'Early stopping at epoch {ep} (patience={patience})')
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    _final_blend = _snapshot_alpha()
    _parts = []
    for _g in group_score_dims:
        if model.moon_alt_conditional_alpha and _g == 'moon':
            _parts.append(f'moon_dn={_final_blend["moon_dn"]:.3f}')
            _parts.append(f'moon_up={_final_blend["moon_up"]:.3f}')
        else:
            _parts.append(f'{_g}={_final_blend[_g]:.3f}')
    print(f"Learned per-group near-arm blend alpha at best epoch [{_blend_mode}]: "
          + '  '.join(_parts))
    if _init_blend is not None:
        init_s = float(blend_init_alpha)
        _deltas = []
        for _g in group_score_dims:
            if model.moon_alt_conditional_alpha and _g == 'moon':
                _deltas.append(f'moon_dn={_final_blend["moon_dn"] - _init_blend["moon_dn"]:+.3f}')
                _deltas.append(f'moon_up={_final_blend["moon_up"] - _init_blend["moon_up"]:+.3f}')
            else:
                _deltas.append(f'{_g}={_final_blend[_g] - _init_blend[_g]:+.3f}')
        print(f'  init alpha={init_s:.3f} (uniform)  |  delta alpha at best epoch: '
              + '  '.join(_deltas))

    # --- Fit per-group empirical mean-bias calibration on train + val rows ---
    # Since 2026-08-16: fit on train+val (test excluded), applied to ALL groups.
    # Previously fit on val only and gated to asinh/log groups -- that left
    # linear-compressor groups (continuum, ionospheric, atomic) with the raw
    # SmoothL1 median-bias, producing a -1.77% systematic mean bias on continuum
    # in the group-bias diagnostic.  The lift is a scalar mean_true/mean_pred
    # ratio: valid for any group with a coherent non-zero mean, so the asinh/log
    # gate was unnecessarily restrictive.  Train+val gives ~85% of rows for
    # groups as small as n=3 (continuum, atomic) where val-only sampling noise
    # dominated the estimator.  Stored as a uniform per-coefficient array so
    # `inverse_group_compressor` continues to accept the same shape.
    jensen_corrections = {}
    _calib_idx = np.concatenate([train_idx, val_idx]).astype(int)
    model.eval()
    with torch.no_grad():
        _calib_pred_dict = model(
            torch.from_numpy(near_s[_calib_idx]).to(device),
            torch.from_numpy(far_s[_calib_idx]).to(device),
            torch.from_numpy(ctx_near_n[_calib_idx]).to(device),
            torch.from_numpy(ctx_far_n[_calib_idx]).to(device),
            torch.from_numpy(ctx_sci_n[_calib_idx]).to(device),
            sci_moon_alt_sign=(torch.from_numpy(_moon_alt_sign_all[_calib_idx]).to(device)
                               if model.moon_alt_conditional_alpha else None),
        )
    _calib_scores_scaled = np.zeros((_calib_idx.size, n_input_score), dtype=np.float64)
    for _g, (_lo, _hi) in score_slices.items():
        _calib_scores_scaled[:, _lo:_hi] = _calib_pred_dict[_g].detach().cpu().numpy().astype(np.float64)
    _calib_pred_scores = score_scaler.inverse_transform(
        _calib_scores_scaled.astype(np.float32)).astype(np.float64)

    _calib_ctx_sci_phys = np.asarray(ctx_sci, dtype=np.float64)[_calib_idx]
    _calib_pred_phys_naive = expand_scores_to_coefs(
        _calib_pred_scores, _calib_ctx_sci_phys,
        compressors, group_indices, geom_kwargs,
        int(coef_sci.shape[1]), score_slices,
        jensen_corrections=None,
    )
    _calib_true_phys = np.asarray(coef_sci, dtype=np.float64)[_calib_idx]

    _CALIB_LIFT_CLIP = (0.5, 2.0)  # sanity bounds; outside indicates a broken group
    # 2026-08-19: moon uses a per-coefficient lift; every other group uses the historical
    # scalar lift.  Rationale: cell 27 shows the moon residual has a spectral tilt bias
    # (+0.054 blue / -0.018 NIR) that a single per-group scalar cannot correct.  Per-coef
    # lift on the 29 Moon_bs knots nulls the tilt by construction.  Tight [0.7, 1.4] clip
    # so a broken knot doesn't cascade.
    # 2026-08-24: zodi promoted to per-coef too. Residual coherent bias is
    # ~4% moon-down / ~5.7% Q4 phase after the ctx expansion; it is shape-
    # dependent (7 Zodi_bs knots have different spectral response) so a
    # single scalar cannot correct it.
    _MOON_LIFT_CLIP = (0.7, 1.4)
    _ZODI_LIFT_CLIP = (0.7, 1.4)
    _PER_COEF_LIFT_GROUPS = ('moon', 'zodi')  # 2026-08-24e: zodi re-enabled with regime-aware per-coef lift (moon_up / moon_horizon / moon_down buckets).
    _ZODI_REGIME_THR_UP = 10.0
    _ZODI_REGIME_THR_DN = -10.0
    _ZODI_REGIME_BOUNDARY_DEG = 5.0
    _ZODI_REGIME_MIN_ROWS = 30
    _ctx_names_calib = list(filtered['ctx_names'])
    _moon_alt_idx_calib = _ctx_names_calib.index('moon_alt') if 'moon_alt' in _ctx_names_calib else None
    if _moon_alt_idx_calib is not None:
        _calib_moon_alt = np.asarray(_calib_ctx_sci_phys[:, _moon_alt_idx_calib], dtype=np.float64)
    else:
        _calib_moon_alt = np.zeros(_calib_ctx_sci_phys.shape[0], dtype=np.float64)
    for _gname, _comp in compressors.items():
        _gidx = np.asarray(_comp['coef_indices'], dtype=int)
        if _gname in _PER_COEF_LIFT_GROUPS:
            _pc_clip = _MOON_LIFT_CLIP if _gname == 'moon' else _ZODI_LIFT_CLIP
            def _per_coef_lift(_mask=None):
                # Per-coef lift on the specified row mask; falls back to global if mask is None.
                if _mask is None:
                    _mt = np.mean(_calib_true_phys[:, _gidx], axis=0).astype(np.float64)
                    _mp = np.mean(_calib_pred_phys_naive[:, _gidx], axis=0).astype(np.float64)
                else:
                    _mt = np.mean(_calib_true_phys[_mask][:, _gidx], axis=0).astype(np.float64)
                    _mp = np.mean(_calib_pred_phys_naive[_mask][:, _gidx], axis=0).astype(np.float64)
                _rel = np.abs(_mp) / np.maximum(np.abs(_mt), 1e-30)
                _lift = np.where(
                    (_mt > 0) & (_mp > 0) & (_rel >= 0.05),
                    _mt / np.where(_mp != 0.0, _mp, 1.0),
                    1.0,
                )
                return np.clip(_lift, _pc_clip[0], _pc_clip[1]).astype(np.float64)
            if _gname == 'zodi':
                # 1D regime lift: 3 moon_alt buckets (moon_up/moon_horizon/moon_down).
                _regime_masks = {
                    'moon_up':      _calib_moon_alt > _ZODI_REGIME_THR_UP,
                    'moon_horizon': (_calib_moon_alt >= _ZODI_REGIME_THR_DN)
                                    & (_calib_moon_alt <= _ZODI_REGIME_THR_UP),
                    'moon_down':    _calib_moon_alt < _ZODI_REGIME_THR_DN,
                }
                _lift_global = _per_coef_lift(None)
                _regime_lifts = {}
                for _rname, _mask in _regime_masks.items():
                    _n_regime = int(_mask.sum())
                    if _n_regime < _ZODI_REGIME_MIN_ROWS:
                        _regime_lifts[_rname] = _lift_global.copy()
                        print(f'Calibration: zodi/{_rname} n={_n_regime} '
                              f'(< {_ZODI_REGIME_MIN_ROWS}); using global lift.')
                        continue
                    _regime_lifts[_rname] = _per_coef_lift(_mask)
                    _l = _regime_lifts[_rname]
                    print(f'Calibration: zodi/{_rname} n={_n_regime} per-coef '
                          f'lift range=[{_l.min():.3f}, {_l.max():.3f}] '
                          f'median={float(np.median(_l)):.3f}')
                jensen_corrections[_gname] = {
                    'moon_up': _regime_lifts['moon_up'],
                    'moon_horizon': _regime_lifts['moon_horizon'],
                    'moon_down': _regime_lifts['moon_down'],
                    'moon_up_threshold_deg': _ZODI_REGIME_THR_UP,
                    'moon_down_threshold_deg': _ZODI_REGIME_THR_DN,
                    'boundary_scale_deg': _ZODI_REGIME_BOUNDARY_DEG,
                }
                continue
            _lift_pc_clipped = _per_coef_lift(None)
            _lift_pc_unclipped_check = _lift_pc_clipped.copy()  # already clipped inside helper
            _n_clipped = 0  # helper always clips; count separately if needed
            jensen_corrections[_gname] = _lift_pc_clipped.astype(np.float64)
            print(f'Calibration: {_gname} per-coef lift range=[{_lift_pc_clipped.min():.3f}, '
                  f'{_lift_pc_clipped.max():.3f}] median={float(np.median(_lift_pc_clipped)):.3f} '
                  f'(clip range {_pc_clip})')
            continue
        _mean_true = float(np.mean(_calib_true_phys[:, _gidx]))
        _mean_pred = float(np.mean(_calib_pred_phys_naive[:, _gidx]))
        _rel_mag = abs(_mean_pred) / max(abs(_mean_true), 1e-30)
        if (not np.isfinite(_mean_true) or not np.isfinite(_mean_pred)
                or _mean_true * _mean_pred <= 0.0
                or _rel_mag < 0.05):
            print(f'Calibration: {_gname} skipped (degenerate or near-zero mean; '
                  f'true={_mean_true:.4g}, pred_naive={_mean_pred:.4g}).')
            continue
        _raw_lift = _mean_true / _mean_pred
        _lift = float(np.clip(_raw_lift, _CALIB_LIFT_CLIP[0], _CALIB_LIFT_CLIP[1]))
        if _raw_lift != _lift:
            print(f'Calibration: {_gname} lift {_raw_lift:.4f} clipped to {_lift:.4f}.')
        jensen_corrections[_gname] = np.full(_gidx.size, _lift, dtype=np.float64)

    if jensen_corrections:
        print('Empirical per-group mean-bias calibration (train+val rows, uniform per-group scalar):')
        print(f"  {'group':<14s} {'n_g':>4s} {'mean_true':>10s} {'mean_pred_naive':>16s} "
              f"{'lift':>7s} {'delta_%':>8s}")
        for _gname, _corr in jensen_corrections.items():
            _gidx = np.asarray(compressors[_gname]['coef_indices'], dtype=int)
            _mean_true = float(np.mean(_calib_true_phys[:, _gidx]))
            _mean_pred = float(np.mean(_calib_pred_phys_naive[:, _gidx]))
            # 2026-08-24e: regime lift dict summary uses the moon_horizon per-coef vector.
            if isinstance(_corr, dict) and ('moon_horizon' in _corr or 'phase_q2_moon_horizon' in _corr):
                _lift_vec = np.asarray((_corr['moon_horizon'] if 'moon_horizon' in _corr else _corr['phase_q2_moon_horizon']), dtype=np.float64)
                _lift = float(np.median(_lift_vec))
                _n_show = int(_lift_vec.size)
                _delta_pct = 100.0 * (_lift - 1.0)
                print(f'  {_gname:<14s} {_n_show:>4d} {_mean_true:>10.4g} '
                      f'{_mean_pred:>16.4g} {_lift:>7.4f} {_delta_pct:>+7.2f}% (regime, median)')
                continue
            _lift = float(_corr[0])
            _delta_pct = 100.0 * (_lift - 1.0)
            print(f'  {_gname:<14s} {len(_corr):>4d} {_mean_true:>10.4g} '
                  f'{_mean_pred:>16.4g} {_lift:>7.4f} {_delta_pct:>+7.2f}%')

    # Per-group upper cap = 3.0 x max(coef_sci_train, axis=0) per coefficient.
    # Defensive guard applied at inference in expand_scores_to_coefs (§11 item 11).
    coef_upper_bound = {}
    _coef_sci_train_arr = np.asarray(coef_sci, dtype=np.float64)[train_idx]
    for _gname, _gidx_local in group_indices.items():
        _gidx_local = np.asarray(_gidx_local, dtype=int)
        if _gidx_local.size == 0:
            continue
        _max_train = np.max(_coef_sci_train_arr[:, _gidx_local], axis=0)
        coef_upper_bound[_gname] = (
            3.0 * np.maximum(_max_train, 0.0)).astype(np.float32)

    return {
        'model': model,
        'device': device,
        'score_scaler': score_scaler,
        'ctx_scaler': ctx_scaler,
        'compressors': compressors,
        'jensen_corrections': jensen_corrections,
        'coef_upper_bound': coef_upper_bound,
        'geom_kwargs': geom_kwargs,
        'group_indices': group_indices,
        'score_slices': score_slices,
        'group_score_dims': group_score_dims,
        'n_input_score': n_input_score,
        'history': history,
        'blend_history': blend_history,
        'blend_init_alpha': float(blend_init_alpha),
        'blend_optim': str(_blend_mode),
        'blend_use_direct': bool(model.blend_use_direct),
        'blend_prefit_alphas': (dict(_prefit_alphas)
                                if _prefit_alphas is not None else None),
        'best_val_loss': float(best_val),
        'best_epoch': int(best_epoch),
        'train_idx': train_idx, 'val_idx': val_idx, 'test_idx': test_idx,
        'coef_names': [str(x) for x in filtered['coef_names']],
        'ctx_names': [str(x) for x in filtered['ctx_names']],
        'config': {
            'n_epochs': int(n_epochs), 'batch_size': int(batch_size),
            'lr': float(lr), 'encoder_dims': tuple(int(v) for v in encoder_dims),
            'ctx_dims': tuple(int(v) for v in ctx_dims),
            'trunk_dims': tuple(int(v) for v in trunk_dims),
            'head_dim': int(head_dim), 'weight_decay': float(weight_decay),
            'head_extra_dims': tuple(int(v) for v in head_extra_dims),
            'zodi_head_extra_dims': tuple(int(v) for v in zodi_head_extra_dims),
            'continuum_head_extra_dims': tuple(int(v) for v in continuum_head_extra_dims),
            'continuum_branch_dims': tuple(int(v) for v in continuum_branch_dims),
            'alpha_ctx_features': (tuple(str(x) for x in alpha_ctx_features)
                                    if alpha_ctx_features else None),
            'alpha_ctx_groups': (tuple(str(x) for x in alpha_ctx_groups)
                                    if alpha_ctx_groups else None),
            'patience': int(patience), 'seed': int(seed),
            'drop_vanrhijn_from_context': bool(drop_vanrhijn_from_context),
            'moon_group_weight': float(moon_group_weight),
            'zodi_group_weight': float(zodi_group_weight),
            'continuum_group_weight': float(continuum_group_weight),
            'mesospheric_group_weight': float(mesospheric_group_weight),
            'ionospheric_group_weight': float(ionospheric_group_weight),
            'blend_init_alpha': float(blend_init_alpha),
            'blend_optim': str(_blend_mode),
            'moon_alt_conditional_alpha': bool(moon_alt_conditional_alpha),
            'high_airmass_boost': float(high_airmass_boost),
            'moon_down_ecliptic_boost': float(moon_down_ecliptic_boost),
            'moon_down_ecliptic_beta_deg': float(moon_down_ecliptic_beta_deg),
            'bright_moon_close_boost': float(bright_moon_close_boost),
            'bright_moon_close_fli_min': float(bright_moon_close_fli_min),
            'bright_moon_close_sep_max_deg': float(bright_moon_close_sep_max_deg),
            'use_coef_err_weights': bool(_weights_enabled),
            # 2026-08-25 opt-in covariance-aware loss; empty tuple = the
            # baseline diagonal Huber loss (backwards compat).
            'block_cov_loss_groups': tuple(block_cov_loss_groups),
            # 2026-08-25b opt-in relative-error MSE per group; takes
            # precedence over block_cov_loss_groups.
            'relative_mse_groups': tuple(relative_mse_groups),
            'relative_mse_eps_frac': float(relative_mse_eps_frac),
            # 2026-08-25c opt-in per-pixel flux MSE per group; takes
            # precedence over relative_mse and block_cov.
            'flux_mse_groups': tuple(flux_mse_groups),
            'flux_mse_eps_frac': float(flux_mse_eps_frac),
            'coef_err_sigma_floor_rel': (
                _resolved_floor_by_group
                if _weights_enabled
                else dict(DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP)
            ),
        },
    }


def predict_sci_coefficients_default(artifacts, coef_near_phys, coef_far_phys,
                                     ctx_near_phys, ctx_far_phys, ctx_sci_phys):
    """Default sky-to-science coefficient predictor (compressed model, see methods §5.5).

    Pipeline: physical coefficients -> divide by geometry factor (§4) ->
    per-group forward compressor (asinh / linear / sqrt + PCA rotation +
    xarm-selected retained subspace) -> encoder + fusion + trunk + heads ->
    inverse compressor -> multiply by science-pointing geometry factor ->
    non-negativity clip.

    When ``artifacts`` is an N-seed ensemble (``is_ensemble=True``), predicts
    with each member and returns the arithmetic mean in physical space
    (the deployed default since 2026-08-11; see §12).
    """
    if artifacts.get('is_ensemble', False):
        _preds = [predict_sci_coefficients_default(
                      m, coef_near_phys, coef_far_phys,
                      ctx_near_phys, ctx_far_phys, ctx_sci_phys)
                  for m in artifacts['members']]
        return np.mean(np.stack(_preds, axis=0), axis=0).astype(np.float32)
    model = artifacts['model']
    device = artifacts['device']
    score_scaler = artifacts['score_scaler']
    ctx_scaler = artifacts['ctx_scaler']
    compressors = artifacts['compressors']
    geom_kwargs = artifacts['geom_kwargs']
    group_indices = artifacts['group_indices']
    score_slices = artifacts['score_slices']
    n_coef = len(artifacts['coef_names'])

    scores_near, _ = compress_coefs_to_scores(
        np.asarray(coef_near_phys, dtype=np.float64),
        np.asarray(ctx_near_phys, dtype=np.float64),
        compressors, geom_kwargs, group_indices)
    scores_far, _ = compress_coefs_to_scores(
        np.asarray(coef_far_phys, dtype=np.float64),
        np.asarray(ctx_far_phys, dtype=np.float64),
        compressors, geom_kwargs, group_indices)

    near_s = np.clip(score_scaler.transform(scores_near.astype(np.float32)),
                     -25.0, 25.0).astype(np.float32)
    far_s = np.clip(score_scaler.transform(scores_far.astype(np.float32)),
                    -25.0, 25.0).astype(np.float32)
    ctx_near_n = np.clip(ctx_scaler.transform(np.asarray(ctx_near_phys, dtype=np.float32)),
                         -25.0, 25.0).astype(np.float32)
    ctx_far_n = np.clip(ctx_scaler.transform(np.asarray(ctx_far_phys, dtype=np.float32)),
                        -25.0, 25.0).astype(np.float32)
    ctx_sci_n = np.clip(ctx_scaler.transform(np.asarray(ctx_sci_phys, dtype=np.float32)),
                        -25.0, 25.0).astype(np.float32)

    _pred_ma_sign = None
    if getattr(model, 'moon_alt_conditional_alpha', False):
        _ctx_names_p = list(artifacts.get('ctx_names', []))
        if 'moon_alt' in _ctx_names_p:
            _ma_idx = _ctx_names_p.index('moon_alt')
            _ma_sign_np = (np.asarray(ctx_sci_phys, dtype=np.float32)[:, _ma_idx] > 0.0
                           ).astype(np.int64)
            _pred_ma_sign = torch.from_numpy(_ma_sign_np).to(device)

    with torch.no_grad():
        pred_dict = model(
            torch.from_numpy(near_s).to(device),
            torch.from_numpy(far_s).to(device),
            torch.from_numpy(ctx_near_n).to(device),
            torch.from_numpy(ctx_far_n).to(device),
            torch.from_numpy(ctx_sci_n).to(device),
            sci_moon_alt_sign=_pred_ma_sign,
        )

    n_rows = near_s.shape[0]
    n_score_total = artifacts['n_input_score']
    pred_scaled = np.zeros((n_rows, n_score_total), dtype=np.float64)
    for g, (lo, hi) in score_slices.items():
        pred_scaled[:, lo:hi] = pred_dict[g].detach().cpu().numpy().astype(np.float64)

    pred_scores = score_scaler.inverse_transform(
        pred_scaled.astype(np.float32)).astype(np.float64)

    coef_predicted = expand_scores_to_coefs(
        pred_scores, np.asarray(ctx_sci_phys, dtype=np.float64),
        compressors, group_indices, geom_kwargs, n_coef, score_slices,
        jensen_corrections=artifacts.get('jensen_corrections'),
        coef_upper_bound=artifacts.get('coef_upper_bound'))
    return coef_predicted.astype(np.float32)


# --- Train the compressed model as an N-seed ensemble, register as the default predictor --
default_dual_group_config: dict[str, Any] = {
    "name": "dual_group_mlp_compressed",
    "n_epochs": 50,
    "batch_size": 512,
    "lr": 1.0e-3,
    "encoder_dims": (768, 384),
    "ctx_dims": (96,),
    "trunk_dims": (320, 160),
    "head_dim": 192,
    # Phase C (2026-08-26): extra hidden layer inside per-group heads.  Empirically
    # neutral on shared-trunk heads (moon, continuum, etc) but a small win on the
    # isolated zodi head; default is therefore no extras on shared heads and one
    # extra layer on the zodi isolated head.  Set to (96,) to reactivate globally.
    "head_extra_dims": (),
    "zodi_head_extra_dims": (32,),
    "continuum_head_extra_dims": (),
    "continuum_branch_dims": (64, 32),
    "weight_decay": 1.0e-4,
    "patience": 12,
    "moon_group_weight": 4.0,
    "zodi_group_weight": 1.0,
    "continuum_group_weight": 1.0,
    "mesospheric_group_weight": 1.0,
    "ionospheric_group_weight": 1.0,
    "block_cov_loss_groups": (),
    "relative_mse_groups": (),
    "relative_mse_eps_frac": 0.05,
    "flux_mse_groups": ("moon", "zodi"),
    "flux_mse_eps_frac": 0.0,
    "high_airmass_boost": 1.0,
    "moon_down_ecliptic_boost": 1.0,
    "moon_down_ecliptic_beta_deg": 15.0,
    "bright_moon_close_boost": 1.0,
    "bright_moon_close_fli_min": 0.85,
    "bright_moon_close_sep_max_deg": 45.0,
    "use_coef_err_weights": True,
    "coef_err_sigma_floor_rel": dict(DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP),
    "ensemble_seeds": (42, 43, 44, 45, 46, 47, 48, 49, 50, 51),
    "zodi_ctx_restriction": (
        "airmass", "vanrhijn_285km",
        "ecl_beta_deg", "ecl_lon_sin", "ecl_lon_cos",
        "zodi_log10_v", "sun_sep",
        "moon_alt", "moon_sep",
        "moon_phase_sin", "moon_phase_cos",
        "moon_fli", "moon_up_smooth",
        "moon_airmass_up", "moon_signal_proxy",
    ),
    # Phase B (2026-08-26): route the continuum head through an isolated branch
    # keyed on moon geometry.  residual_ctx_attribution flagged continuum with
    # RF R^2 = 0.70 driven by moon_fli + moon_phase_cos + moon_sep.  Empty tuple
    # / None reverts to the shared-trunk head.
    "continuum_ctx_restriction": (
        "moon_alt", "moon_sep",
        "moon_phase_sin", "moon_phase_cos",
        "moon_fli", "airmass",
    ),
    # Phase A (2026-08-26): context-dependent per-group blend alpha.
    # Empty tuple / None reverts to a single learnable scalar per group.
    "alpha_ctx_features": (
        "moon_up_smooth", "ecl_beta_deg", "airmass",
    ),
    # Phase A restriction: which groups get a ctx-alpha predictor.  Only the
    # groups whose sky-to-sci transfer is anisotropic (moon, zodi, continuum).
    # Mesospheric / ionospheric / atomic are LOS-integrated thin shells that
    # want alpha near 0.5 with no ctx dependence.  None -> all groups.
    "alpha_ctx_groups": ("moon", "zodi", "continuum"),
}


_WAVE_STRIDE_FLUX_LOSS = 5


def _precompute_flux_basis_and_geometry(
    *,
    filtered_triplet,
    compress_geom_kwargs,
    input_fits_for_basis,
    n_moon_knots,
    split_zodi,
    n_zodi_knots,
    palace_oh_suffix="_joint_v2_updated",
    palace_diffuse_suffix="_joint_native_adam_invsky_p2_10000iter",
    verbose=True,
):
    """Materialise moon + zodi flux basis matrices (stride 5) and per-row geometry."""

    with fits.open(str(input_fits_for_basis)) as hdul:
        wave_ref = np.asarray(hdul["WAVE"].data, dtype=np.float64)
    model = SkyDecompLSFSurfaceIterative(
        wave_ref, lsf_sigma=1.0, n_spline_knots=n_moon_knots,
        base_dir=_infer_base_dir_for_reconstruction(),
        palace_oh_suffix=palace_oh_suffix,
        palace_diffuse_suffix=palace_diffuse_suffix,
        split_zodi=split_zodi, n_zodi_spline_knots=n_zodi_knots,
    )
    stride = slice(None, None, int(_WAVE_STRIDE_FLUX_LOSS))
    flux_basis_matrices = {
        "moon": np.asarray(model.matrix_moon[:, stride], dtype=np.float32),
    }
    if split_zodi and model.matrix_zodi.shape[0] > 0:
        flux_basis_matrices["zodi"] = np.asarray(
            model.matrix_zodi[:, stride], dtype=np.float32
        )
    flux_geom_sc_sci = airglow_geometry_scale(
        filtered_triplet["ctx_sci"], **compress_geom_kwargs
    ).astype(np.float32)
    if verbose:
        print(
            f"[flux-mse prep] wave grid: n_full={wave_ref.size}, "
            f"stride={_WAVE_STRIDE_FLUX_LOSS}, "
            f"n_ds={flux_basis_matrices['moon'].shape[1]}"
        )
        for k, m in flux_basis_matrices.items():
            print(
                f"[flux-mse prep] {k}: basis shape={m.shape}, "
                f"||A[k]|| range=[{float(np.linalg.norm(m, axis=1).min()):.3g}, "
                f"{float(np.linalg.norm(m, axis=1).max()):.3g}]"
            )
        print(
            f"[flux-mse prep] geom sc_sci: shape={flux_geom_sc_sci.shape}, "
            f"median={float(np.median(flux_geom_sc_sci)):.3g}"
        )
    return flux_basis_matrices, flux_geom_sc_sci


@dataclass
class EnsembleArtifacts:
    """Combined artifacts of a full ensemble train run."""

    seeds: list
    members: list
    mlp_artifacts: dict
    per_seed_test_metrics: pd.DataFrame
    seed_std_rmse: float
    ensemble_stderr: float


class Trainer:
    """Orchestrates the ensemble training loop.

    Notebook workflow condenses to::

        trainer = Trainer(cfg=default_dual_group_config)
        artifacts = trainer.run_ensemble(
            filtered_triplet, compressors, group_indices, geom_kwargs,
            input_fits_for_basis=..., n_moon_knots=..., split_zodi=..., n_zodi_knots=...,
        )
    """

    def __init__(self, cfg=None):
        self.cfg = dict(cfg) if cfg is not None else dict(default_dual_group_config)

    def _shared_train_kwargs(self, *, flux_basis_matrices, flux_geom_sc_sci):
        c = self.cfg
        return dict(
            n_epochs=int(c["n_epochs"]),
            batch_size=int(c["batch_size"]),
            lr=float(c["lr"]),
            encoder_dims=tuple(int(v) for v in c["encoder_dims"]),
            ctx_dims=tuple(int(v) for v in c["ctx_dims"]),
            trunk_dims=tuple(int(v) for v in c["trunk_dims"]),
            head_dim=int(c["head_dim"]),
            head_extra_dims=tuple(int(v) for v in c.get("head_extra_dims", ())),
            zodi_head_extra_dims=tuple(int(v) for v in c.get("zodi_head_extra_dims", ())),
            continuum_head_extra_dims=tuple(int(v) for v in c.get("continuum_head_extra_dims", ())),
            continuum_branch_dims=tuple(int(v) for v in c.get("continuum_branch_dims", (64, 32))),
            weight_decay=float(c["weight_decay"]),
            patience=int(c["patience"]),
            moon_group_weight=float(c["moon_group_weight"]),
            zodi_group_weight=float(c["zodi_group_weight"]),
            continuum_group_weight=float(c["continuum_group_weight"]),
            mesospheric_group_weight=float(c["mesospheric_group_weight"]),
            ionospheric_group_weight=float(c["ionospheric_group_weight"]),
            high_airmass_boost=float(c["high_airmass_boost"]),
            moon_down_ecliptic_boost=float(c["moon_down_ecliptic_boost"]),
            moon_down_ecliptic_beta_deg=float(c["moon_down_ecliptic_beta_deg"]),
            bright_moon_close_boost=float(c["bright_moon_close_boost"]),
            bright_moon_close_fli_min=float(c["bright_moon_close_fli_min"]),
            bright_moon_close_sep_max_deg=float(c["bright_moon_close_sep_max_deg"]),
            use_coef_err_weights=bool(c["use_coef_err_weights"]),
            coef_err_sigma_floor_rel=c["coef_err_sigma_floor_rel"],
            zodi_ctx_restriction=c.get("zodi_ctx_restriction"),
            continuum_ctx_restriction=c.get("continuum_ctx_restriction"),
            alpha_ctx_features=c.get("alpha_ctx_features"),
            alpha_ctx_groups=c.get("alpha_ctx_groups"),
            block_cov_loss_groups=tuple(c.get("block_cov_loss_groups", ())),
            relative_mse_groups=tuple(c.get("relative_mse_groups", ())),
            relative_mse_eps_frac=float(c.get("relative_mse_eps_frac", 0.05)),
            flux_mse_groups=tuple(c.get("flux_mse_groups", ())),
            flux_mse_eps_frac=float(c.get("flux_mse_eps_frac", 0.0)),
            flux_basis_matrices=flux_basis_matrices,
            flux_geom_sc_sci=flux_geom_sc_sci,
        )

    def run_ensemble(
        self,
        filtered_triplet,
        compressors,
        group_indices,
        geom_kwargs,
        *,
        input_fits_for_basis,
        n_moon_knots,
        split_zodi,
        n_zodi_knots,
        palace_oh_suffix="_joint_v2_updated",
        palace_diffuse_suffix="_joint_native_adam_invsky_p2_10000iter",
        verbose=True,
    ):
        """Fit the full seed ensemble, assemble artifacts, report per-seed metrics."""

        seeds = tuple(int(s) for s in self.cfg["ensemble_seeds"])
        if verbose:
            print(
                f"=== Training compressed dual-encoder group-head MLP "
                f"({len(seeds)}-seed ensemble default) ==="
            )
            print(self.cfg)

        flux_basis_matrices = None
        flux_geom_sc_sci = None
        if self.cfg.get("flux_mse_groups"):
            flux_basis_matrices, flux_geom_sc_sci = _precompute_flux_basis_and_geometry(
                filtered_triplet=filtered_triplet,
                compress_geom_kwargs=geom_kwargs,
                input_fits_for_basis=input_fits_for_basis,
                n_moon_knots=n_moon_knots,
                split_zodi=split_zodi,
                n_zodi_knots=n_zodi_knots,
                palace_oh_suffix=palace_oh_suffix,
                palace_diffuse_suffix=palace_diffuse_suffix,
                verbose=verbose,
            )

        shared = self._shared_train_kwargs(
            flux_basis_matrices=flux_basis_matrices,
            flux_geom_sc_sci=flux_geom_sc_sci,
        )
        split_for_members = (
            filtered_triplet["compress_train_idx"],
            filtered_triplet["compress_val_idx"],
            filtered_triplet["compress_test_idx"],
        )
        members = []
        for seed in seeds:
            if verbose:
                print(f"\n--- Ensemble member seed={seed} ---")
            member = train_compressed_group_mlp(
                filtered_triplet, compressors, group_indices, geom_kwargs,
                split_indices=split_for_members,
                seed=int(seed),
                **shared,
            )
            members.append(member)

        first = members[0]
        mlp_artifacts = {
            "is_ensemble": True,
            "seeds": list(seeds),
            "members": members,
            "compressors": first["compressors"],
            "coef_upper_bound": first["coef_upper_bound"],
            "geom_kwargs": first["geom_kwargs"],
            "group_indices": first["group_indices"],
            "score_slices": first["score_slices"],
            "group_score_dims": first["group_score_dims"],
            "n_input_score": first["n_input_score"],
            "coef_names": first["coef_names"],
            "ctx_names": first["ctx_names"],
            "train_idx": first["train_idx"],
            "val_idx": first["val_idx"],
            "test_idx": first["test_idx"],
            "config": first["config"],
            "best_epochs": [int(m["best_epoch"]) for m in members],
            "best_val_losses": [float(m["best_val_loss"]) for m in members],
        }
        if verbose:
            print(
                f"\nEnsemble assembled: {len(members)} members, "
                f"best_epochs={mlp_artifacts['best_epochs']}, "
                f"best_val_losses={[f'{v:.5f}' for v in mlp_artifacts['best_val_losses']]}"
            )

        per_seed_df, seed_std, stderr = self._report_test_metrics(
            filtered_triplet, mlp_artifacts, group_indices, verbose=verbose,
        )
        return EnsembleArtifacts(
            seeds=list(seeds),
            members=members,
            mlp_artifacts=mlp_artifacts,
            per_seed_test_metrics=per_seed_df,
            seed_std_rmse=seed_std,
            ensemble_stderr=stderr,
        )

    def _report_test_metrics(self, filtered_triplet, mlp_artifacts, group_indices, verbose=True):
        test_idx = np.asarray(mlp_artifacts["test_idx"], dtype=int)
        coef_near = np.asarray(filtered_triplet["coef_near"], dtype=np.float32)
        coef_far = np.asarray(filtered_triplet["coef_far"], dtype=np.float32)
        coef_sci = np.asarray(filtered_triplet["coef_sci"], dtype=np.float32)
        ctx_near = np.asarray(filtered_triplet["ctx_near"], dtype=np.float32)
        ctx_far = np.asarray(filtered_triplet["ctx_far"], dtype=np.float32)
        ctx_sci = np.asarray(filtered_triplet["ctx_sci"], dtype=np.float32)
        coef_err_sci = np.asarray(
            filtered_triplet.get("coef_err_sci", np.full_like(coef_sci, np.nan)),
            dtype=np.float32,
        )
        y_te = coef_sci[test_idx]
        sig_te = coef_err_sci[test_idx]
        floor_te = dict(DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP)

        coef_pred_det = predict_sci_coefficients_default(
            mlp_artifacts,
            coef_near_phys=coef_near[test_idx],
            coef_far_phys=coef_far[test_idx],
            ctx_near_phys=ctx_near[test_idx],
            ctx_far_phys=ctx_far[test_idx],
            ctx_sci_phys=ctx_sci[test_idx],
        ).astype(np.float32)

        rows = []
        for seed, member in zip(mlp_artifacts["seeds"], mlp_artifacts["members"]):
            pred_seed = predict_sci_coefficients_default(
                member,
                coef_near_phys=coef_near[test_idx],
                coef_far_phys=coef_far[test_idx],
                ctx_near_phys=ctx_near[test_idx],
                ctx_far_phys=ctx_far[test_idx],
                ctx_sci_phys=ctx_sci[test_idx],
            ).astype(np.float32)
            rows.append({
                "variant": f"seed={seed}",
                **metric_row(y_te, pred_seed, f"seed={seed}",
                             sigma=sig_te, group_indices=group_indices,
                             floor_by_group=floor_te),
            })
        rows.append({
            "variant": f"{len(mlp_artifacts['members'])}-seed ensemble (default)",
            **metric_row(y_te, coef_pred_det, "ensemble",
                         sigma=sig_te, group_indices=group_indices,
                         floor_by_group=floor_te),
        })
        df = pd.DataFrame(rows)
        if verbose:
            print(
                f"\nPer-seed and ensemble test metrics on the night-held-out split "
                f"(n_test = {test_idx.size} rows):"
            )
            print(df.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
        per_seed_rmse = np.array([r["mean_eRMSE"] for r in rows[:len(mlp_artifacts["members"])]])
        seed_std = float(np.std(per_seed_rmse, ddof=1))
        ensemble_rmse = float(df.iloc[-1]["mean_eRMSE"])
        stderr = seed_std / (len(mlp_artifacts["members"]) ** 0.5)
        if verbose:
            print(
                f"\nSeed-to-seed test mean_eRMSE std: {seed_std:.3f}  "
                f"(single-seed noise scale on this dataset)"
            )
            print(
                f"Ensemble-mean stderr = std/sqrt(N={len(mlp_artifacts['members'])}) = "
                f"{stderr:.3f}  (~{100.0 * stderr / max(ensemble_rmse, 1e-30):.1f}% of "
                f"ensemble mean_eRMSE {ensemble_rmse:.3f})."
            )
        return df, seed_std, stderr


__all__ = [
    "DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP",
    "EnsembleArtifacts",
    "Trainer",
    "default_dual_group_config",
    "predict_sci_coefficients_default",
    "train_compressed_group_mlp",
]
