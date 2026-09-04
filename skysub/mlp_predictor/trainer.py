"""Ensemble trainer and inference for the dual-encoder group-head MLP.

Every knob in ``default_dual_group_config`` is read by ``Trainer``; the class
asserts this on construction and warns about any config key it does not
consume.

Training
--------
* Targets are per-group compressed *scores*, not raw coefficients: each group
  is transformed (``sqrt`` for the continuum families, ``asinh`` for the 358
  mesospheric coefficients), centred and scaled.  ``compress_coefs_to_scores``
  and ``expand_scores_to_coefs`` are the two directions.
* The loss is a per-group mean, weighted by ``group_loss_weight`` = the config
  weight over ``sqrt(n_group)`` so a 358-coefficient group does not swamp a
  3-coefficient one.  Row weights are uniform.
* ``moon`` and ``zodi`` are scored in FLUX space rather than coefficient
  space: their predicted scores are inverted through the compressor and their
  own basis matrices to per-pixel flux, and the MSE is taken there.  This
  replaces (not augments) the coefficient-space ``smooth_l1`` for those two
  groups, and it is load-bearing -- removing it costs the moon 18 percentage
  points of its gain over copying the near arm, and also degrades the
  continuum and mesospheric groups, which have no flux term of their own,
  through the shared trunk.
* Remaining groups use ``smooth_l1`` weighted per element by
  ``1/sigma^2`` from the decomposition's own ``COEF_ERR``, with a per-group
  relative floor (``coef_err_sigma_floor_rel``) so near-zero uncertainties
  cannot dominate.
* The blend alphas train with their own learning-rate multiplier
  (``alpha_lr_mult``); at the shared rate they do not move measurably.

Ensembling and calibration
--------------------------
``run_ensemble`` trains one member per seed on a night-disjoint,
moon-phase-stratified split, then fits a per-group scalar mean-bias
correction (a Jensen-style lift, since the compressor transforms are convex)
on the training and validation rows.
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

DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP = {
    # Per-group floors on the relative sigma (fraction of the per-column
    # median finite sigma).  Set on 2026-08-16 from the sigma / weight
    # diagnostic (§7): moon / continuum / ionospheric take the historical
    # 5% floor; mesospheric (403 lines) and atomic have p99 / p50 sigma
    # ratios of 10^3-10^6 in compressed-score space, so they need a 20%
    # floor to keep w_p99 / w_median below ~30 per column.
    'moon': 0.05,
    'zodi': 0.05,
    'continuum': 0.05,
    'mesospheric': 0.20,
    'ionospheric': 0.05,
    'atomic': 0.20,
}


def train_compressed_group_mlp(
    filtered, compressors, group_indices, geom_kwargs,
    split_indices=None,
    n_epochs=50, batch_size=256, lr=7e-4,
    encoder_dims=(768, 384), ctx_dims=(96,),
    trunk_dims=(320, 160), head_dim=192,
    zodi_head_extra_dims=(32,),
    continuum_head_extra_dims=(64,),
    continuum_branch_dims=(128, 64),
    moon_zodi_coupling_dims=(64, 32),
    weight_decay=1e-4, grad_clip=1.0, patience=4,
    seed=42,
    moon_group_weight=1.0,
    zodi_group_weight=1.0,
    continuum_group_weight=1.0,
    mesospheric_group_weight=1.0,
    ionospheric_group_weight=1.0,
    blend_init_alpha=0.7,
    # 2026-09-01: learning-rate multiplier for the blend-alpha parameters
    # (``blend_alpha_direct`` + the ctx-alpha predictors).  At the default 1.0
    # alpha barely moves from its init within the early-stopping budget -- the
    # loss surface in alpha is nearly flat because the additive head absorbs any
    # systematic part of a mis-set blend, and the ctx path is damped further by
    # the sigmoid (sigma' = 0.25 at alpha=0.5, 0.128 at 0.85).  Measured travel
    # over a run: scalar groups +0.08..+0.24, ctx groups +0.01..+0.05, always
    # toward the near arm and still climbing at early stop.  Raising this lets
    # alpha reach its own optimum so the init stops being a hyperparameter.
    alpha_lr_mult=1.0,
    alpha_ctx_features=("moon_up_smooth", "ecl_beta_deg", "airmass"),
    zodi_ctx_restriction=(),
    continuum_ctx_restriction=(),
    moon_zodi_ctx_restriction=(),
    # 2026-08-25c per-pixel flux MSE for moon and/or zodi (deployed default).
    flux_mse_groups=(),
    flux_basis_matrices=None,
    flux_geom_sc_sci=None,
    # 2026-08-21 tight-mask bright-moon row boost (Phase A'').
    # Heteroscedastic-Gaussian loss weighting from decomposition COEF_ERR.
    coef_err_sigma_floor_rel=None,
):
    """Train one seed of the compressed dual-encoder group-head MLP."""
    set_reproducibility(seed)

    ctx_near = np.asarray(filtered['ctx_near'], dtype=np.float32)
    ctx_far = np.asarray(filtered['ctx_far'], dtype=np.float32)
    ctx_sci = np.asarray(filtered['ctx_sci'], dtype=np.float32)

    if split_indices is not None:
        train_idx, val_idx, test_idx = split_indices
    else:
        _moon_phase_col = moon_phase_deg_from_ctx(filtered)
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

    # Per-group multiplier m_g on top of 1/sqrt(n_g_score).
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
    # divide out the score_scaler scale so weights live in the same
    # scaled-score space as sci_s; clip to per-column floor and normalise
    # so E_train[w_pe] = 1 per column.  Missing/boundary sigmas fall
    # through to the floor.  See §7 of the methods cell for derivation.
    coef_err_sci_local = np.asarray(
        filtered.get('coef_err_sci', np.full_like(coef_sci, np.nan)),
        dtype=np.float64,
    )
    _n_finite_err = int(np.sum(np.isfinite(coef_err_sci_local)))
    _n_total_err = int(coef_err_sci_local.size)
    if _n_finite_err == 0:
        raise RuntimeError(
            'coef_err_sci is entirely non-finite; the deployed pipeline requires '
            'the pipeline COEF_ERR HDU to be present.')
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
    _floor_arg = coef_err_sigma_floor_rel or DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP
    _floor_by_group = {str(k): float(v) for k, v in _floor_arg.items()}
    _floor_rel_col = np.ones(_floor_col.shape[0], dtype=np.float64)
    _missing_groups = []
    for _g, (_lo, _hi) in score_slices.items():
        _fr = _floor_by_group.get(_g)
        if _fr is None:
            _missing_groups.append(_g)
            _fr = 0.05
        _floor_rel_col[_lo:_hi] = _fr
    if _missing_groups:
        print(f'  warning: coef_err_sigma_floor_rel missing groups '
              f'{_missing_groups}; using fallback 0.05 for each.')
    _floor_col = _floor_rel_col * _floor_col
    _sigma_scaled = np.where(_sigma_scaled > 0.0, _sigma_scaled,
                             _floor_col[None, :])
    _sigma_scaled = np.maximum(_sigma_scaled, _floor_col[None, :])
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
    _resolved_floor_by_group = {g: float(_floor_by_group.get(g, 0.05))
                                for g in score_slices}

    ctx_near_n = np.clip(ctx_scaler.transform(ctx_near), -25.0, 25.0).astype(np.float32)
    ctx_far_n = np.clip(ctx_scaler.transform(ctx_far), -25.0, 25.0).astype(np.float32)
    ctx_sci_n = np.clip(ctx_scaler.transform(ctx_sci), -25.0, 25.0).astype(np.float32)

    # Row weights are uniform.  Every per-row reweighting scheme tried was
    # rejected on measurement: high_airmass and moon_down_ecliptic boosts never
    # helped, and the bright-moon-close boost (fli>=0.90, sep<=30 deg, alt>0,
    # x1.5) was carried at 1.0 (= off) because it leaked into the blue atlas
    # band (+23% RMS|frac|) for a mid-band gain the flux-space loss now
    # provides directly.  Kept as a named constant rather than a knob so the
    # loss has one less silent degree of freedom.
    _row_weights_np = np.ones(ctx_sci.shape[0], dtype=np.float32)

    # ------------------------------------------------------------------
    # Per-pixel flux MSE precomputation for moon / zodi (deployed default).
    # For each group g, undo the compressor + geometry to recover native
    # coefs; A_g @ c gives flux.  Loss = mean_lambda((f_pred - f_true)^2)
    # per row, scale-matched to the diagonal Huber magnitude on train rows.
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

    # Device selection.
    if torch.cuda.is_available():
        device = 'cuda'
    elif getattr(torch.backends, 'mps', None) is not None and torch.backends.mps.is_available():
        device = 'mps'
    else:
        device = 'cpu'

    # 2026-08-24: skip DataLoader for the ~O(10MB) in-memory training set.
    # Move all tensors to device once and iterate with torch.randperm --
    # avoids ~180k per-batch CPU->device transfers per full 10-seed run.
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
        return (
            torch.from_numpy(near_s[_idx]).to(device),
            torch.from_numpy(far_s[_idx]).to(device),
            torch.from_numpy(ctx_near_n[_idx]).to(device),
            torch.from_numpy(ctx_far_n[_idx]).to(device),
            torch.from_numpy(ctx_sci_n[_idx]).to(device),
            torch.from_numpy(_row_weights_np[_idx]).to(device),
            torch.from_numpy(w_pe_np[_idx]).to(device),
            torch.from_numpy(sci_s[_idx]).to(device),
            torch.from_numpy(_row_idx_all_np[_idx]).to(device),
        )
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
        zodi_head_extra_dims=tuple(int(v) for v in zodi_head_extra_dims),
        continuum_head_extra_dims=tuple(int(v) for v in continuum_head_extra_dims),
        continuum_branch_dims=tuple(int(v) for v in continuum_branch_dims),
        moon_zodi_coupling_dims=tuple(int(v) for v in moon_zodi_coupling_dims),
        blend_init_alpha=float(blend_init_alpha),
        alpha_ctx_features=alpha_ctx_features,
        zodi_ctx_restriction=zodi_ctx_restriction,
        continuum_ctx_restriction=continuum_ctx_restriction,
        moon_zodi_ctx_restriction=moon_zodi_ctx_restriction,
    ).to(device)

    # Blend params (per-group alpha) get weight_decay=0, an optional LR boost,
    # and -- for the direct-parametrised ones -- clamping to [eps, 1-eps] after
    # each optimizer step.  The ctx-alpha predictors join this group: they are
    # blend parameters too, and leaving them under weight decay was never the
    # intent.  (At alpha_lr_mult=1.0 that regrouping is numerically a no-op:
    # decoupled decay would have shrunk them by lr*wd per step, i.e. <0.03% over
    # a full run.)
    blend_pnames = {f'blend_alpha_direct.{_k}' for _k in model.blend_alpha_direct}
    blend_params, other_params = [], []
    for _n, _p in model.named_parameters():
        if _n in blend_pnames or _n.startswith('alpha_predictors.'):
            blend_params.append(_p)
        else:
            other_params.append(_p)
    _alpha_lr = float(lr) * float(alpha_lr_mult)
    opt = torch.optim.AdamW(
        [{'params': other_params, 'weight_decay': float(weight_decay), 'lr': float(lr)},
         {'params': blend_params, 'weight_decay': 0.0, 'lr': _alpha_lr}],
        lr=float(lr))
    print(f"Blend optim: direct (n_blend_params={sum(p.numel() for p in blend_params)}, "
          f"alpha_lr={_alpha_lr:.2e} = {float(alpha_lr_mult):g} x lr)")

    def compressed_loss(pred_dict, yb, w_row, w_pe, row_idx_b):
        loss = torch.tensor(0.0, device=yb.device)
        for g, y_head in pred_dict.items():
            lo, hi = score_slices[g]
            target = yb[:, lo:hi].contiguous()
            if g in _flux_mse_torch_by_group:
                # Per-pixel flux MSE for moon and/or zodi.  Inverse compressor
                # in torch: undo RobustScaler, PCA, sqrt, geometry.
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
                _flux_true = (_em_true * _g_row) @ _st_fm['A']
                _flux_pred = (_em_pred * _g_row) @ _st_fm['A']
                _per_row = ((_flux_pred - _flux_true) ** 2).mean(dim=1) * float(_st_fm['scale_match'])
            else:
                w_pe_g = w_pe[:, lo:hi].contiguous()
                _per_elem = F.smooth_l1_loss(y_head, target, reduction='none') * w_pe_g
                _per_row = _per_elem.mean(dim=1)
            loss = loss + float(group_loss_weight[g]) * (w_row * _per_row).mean()
        return loss / max(len(pred_dict), 1)

    def _snapshot_alpha():
        return {str(g): float(model.blend_alpha_direct[str(g)].item())
                for g in group_score_dims}

    history = []
    blend_history = []
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
            near_b, far_b, near_ctx_b, far_ctx_b, ctx_b, w_row_b, w_pe_b, yb, row_idx_b = tuple(
                t[_sel] for t in _tr_tensors)
            pred = model(near_b, far_b, near_ctx_b, far_ctx_b, ctx_b)
            loss = compressed_loss(pred, yb, w_row_b, w_pe_b, row_idx_b)
            if not torch.isfinite(loss):
                continue
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), float(grad_clip))
            opt.step()
            with torch.no_grad():
                _eps = float(model.blend_alpha_eps)
                for _g in group_score_dims:
                    model.blend_alpha_direct[str(_g)].clamp_(_eps, 1.0 - _eps)
            tr_loss += float(loss.item())
            tr_n += 1

        model.eval()
        va_loss = 0.0
        va_n = 0
        with torch.no_grad():
            for _b0 in range(0, _n_val, _va_bs):
                _sel = slice(_b0, _b0 + _va_bs)
                near_b, far_b, near_ctx_b, far_ctx_b, ctx_b, w_row_b, w_pe_b, yb, row_idx_b = tuple(
                    t[_sel] for t in _va_tensors)
                pred = model(near_b, far_b, near_ctx_b, far_ctx_b, ctx_b)
                loss = compressed_loss(pred, yb, w_row_b, w_pe_b, row_idx_b)
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
    print(f"Learned per-group near-arm blend alpha at best epoch: "
          + '  '.join(f'{_g}={_final_blend[_g]:.3f}' for _g in group_score_dims))
    _deltas = [f'{_g}={_final_blend[_g] - _init_blend[_g]:+.3f}' for _g in group_score_dims]
    print(f'  init alpha={float(blend_init_alpha):.3f} (uniform)  |  delta alpha at best epoch: '
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
        'alpha_lr_mult': float(alpha_lr_mult),
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
            'zodi_head_extra_dims': tuple(int(v) for v in zodi_head_extra_dims),
            'continuum_head_extra_dims': tuple(int(v) for v in continuum_head_extra_dims),
            'continuum_branch_dims': tuple(int(v) for v in continuum_branch_dims),
            'zodi_ctx_restriction': tuple(str(x) for x in zodi_ctx_restriction),
            'continuum_ctx_restriction': tuple(str(x) for x in continuum_ctx_restriction),
            'moon_zodi_ctx_restriction': tuple(str(x) for x in moon_zodi_ctx_restriction),
            'moon_zodi_coupling_dims': tuple(int(v) for v in moon_zodi_coupling_dims),
            'alpha_ctx_features': tuple(str(x) for x in alpha_ctx_features),
            'patience': int(patience), 'seed': int(seed),
            'moon_group_weight': float(moon_group_weight),
            'zodi_group_weight': float(zodi_group_weight),
            'continuum_group_weight': float(continuum_group_weight),
            'mesospheric_group_weight': float(mesospheric_group_weight),
            'ionospheric_group_weight': float(ionospheric_group_weight),
            'blend_init_alpha': float(blend_init_alpha),
            'alpha_lr_mult': float(alpha_lr_mult),
            'flux_mse_groups': tuple(flux_mse_groups),
            'coef_err_sigma_floor_rel': _resolved_floor_by_group,
        },
    }


def predict_sci_coefficients_default(artifacts, coef_near_phys, coef_far_phys,
                                     ctx_near_phys, ctx_far_phys, ctx_sci_phys):
    """Default sky-to-science coefficient predictor (compressed model, §5.5).

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

    with torch.no_grad():
        pred_dict = model(
            torch.from_numpy(near_s).to(device),
            torch.from_numpy(far_s).to(device),
            torch.from_numpy(ctx_near_n).to(device),
            torch.from_numpy(ctx_far_n).to(device),
            torch.from_numpy(ctx_sci_n).to(device),
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


# --- Deployed ensemble config (matches the shipped mlp_ensemble_split_zodi_current.pt) ---
default_dual_group_config: dict[str, Any] = {
    "name": "dual_group_mlp_compressed",
    "n_epochs": 50,
    "batch_size": 512,
    "lr": 1.0e-3,
    "encoder_dims": (768, 384),
    "ctx_dims": (96,),
    "trunk_dims": (320, 160),
    "head_dim": 192,
    "zodi_head_extra_dims": (32,),
    "continuum_head_extra_dims": (64,),
    "continuum_branch_dims": (128, 64),
    "moon_zodi_coupling_dims": (64, 32),
    "blend_init_alpha": 0.7,
    "alpha_lr_mult": 1.0,
    "weight_decay": 1.0e-4,
    "patience": 12,
    "moon_group_weight": 2.0,
    "zodi_group_weight": 2.0,
    "continuum_group_weight": 1.0,
    "mesospheric_group_weight": 1.0,
    "ionospheric_group_weight": 1.0,
    "flux_mse_groups": ("moon", "zodi"),
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
    "continuum_ctx_restriction": (
        "moon_alt", "moon_sep",
        "moon_phase_sin", "moon_phase_cos",
        "moon_fli", "airmass",
        "moon_fli_x_phase_cos",
        "moon_sig_x_lon_cos", "moon_sig_x_lon_sin",
    ),
    "moon_zodi_ctx_restriction": (
        "airmass", "vanrhijn_285km",
        "ecl_beta_deg", "ecl_lon_sin", "ecl_lon_cos",
        "zodi_log10_v", "sun_sep",
        "moon_alt", "moon_sep",
        "moon_phase_sin", "moon_phase_cos",
        "moon_fli", "moon_up_smooth",
        "moon_airmass_up", "moon_signal_proxy",
        "moon_fli_x_phase_cos",
        "moon_sig_x_lon_cos", "moon_sig_x_lon_sin",
    ),
    "alpha_ctx_features": ("moon_up_smooth", "ecl_beta_deg", "airmass"),
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
    palace_oh_suffix=None,
    palace_diffuse_suffix=None,
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

    #: Every config key the Trainer actually reads.  ``_shared_train_kwargs``
    #: reads all of these except ``ensemble_seeds`` / ``name``, which
    #: ``run_ensemble`` handles.  Keep in sync when adding a knob -- the
    #: ``test_consumed_cfg_keys_match_source`` check greps the class body for
    #: ``c["..."]`` / ``c.get("...")`` and compares against this set.
    _CONSUMED_CFG_KEYS = frozenset({
        "name", "ensemble_seeds",
        "n_epochs", "batch_size", "lr", "weight_decay", "patience",
        "encoder_dims", "ctx_dims", "trunk_dims", "head_dim",
        "zodi_head_extra_dims", "continuum_head_extra_dims",
        "continuum_branch_dims", "moon_zodi_coupling_dims",
        "moon_zodi_ctx_restriction", "zodi_ctx_restriction",
        "continuum_ctx_restriction", "alpha_ctx_features",
        "blend_init_alpha", "alpha_lr_mult",
        "moon_group_weight", "zodi_group_weight", "continuum_group_weight",
        "mesospheric_group_weight", "ionospheric_group_weight",
        "coef_err_sigma_floor_rel", "flux_mse_groups",
    })

    def __init__(self, cfg=None):
        self.cfg = dict(cfg) if cfg is not None else dict(default_dual_group_config)
        # Warn loudly about knobs the Trainer will not read.  The 2026-08-27d
        # dead-code sweep removed a batch of them from the trainer signature
        # but callers (notebook cell 9) still carry the pre-sweep list, so an
        # edit to one of those entries silently does nothing -- exactly how a
        # blend_init_alpha=0.5 A/B ran at the 0.7 default on 2026-09-01.
        _ignored = sorted(set(self.cfg) - self._CONSUMED_CFG_KEYS)
        if _ignored:
            print(
                f"[Trainer] WARNING: {len(_ignored)} config key(s) are NOT read by "
                f"the trainer and will have no effect on this run:\n"
                f"           {', '.join(_ignored)}\n"
                f"           Remove them, or check the spelling if you meant to "
                f"change behaviour."
            )
        _missing = sorted(self._CONSUMED_CFG_KEYS
                          - {"name", "ensemble_seeds"} - set(self.cfg))
        if _missing:
            raise KeyError(
                f"Trainer cfg is missing required key(s): {', '.join(_missing)}. "
                f"Start from mlp_predictor.trainer.default_dual_group_config.")

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
            zodi_head_extra_dims=tuple(int(v) for v in c["zodi_head_extra_dims"]),
            continuum_head_extra_dims=tuple(int(v) for v in c["continuum_head_extra_dims"]),
            continuum_branch_dims=tuple(int(v) for v in c["continuum_branch_dims"]),
            moon_zodi_ctx_restriction=c["moon_zodi_ctx_restriction"],
            moon_zodi_coupling_dims=tuple(int(v) for v in c["moon_zodi_coupling_dims"]),
            blend_init_alpha=float(c["blend_init_alpha"]),
            alpha_lr_mult=float(c["alpha_lr_mult"]),
            weight_decay=float(c["weight_decay"]),
            patience=int(c["patience"]),
            moon_group_weight=float(c["moon_group_weight"]),
            zodi_group_weight=float(c["zodi_group_weight"]),
            continuum_group_weight=float(c["continuum_group_weight"]),
            mesospheric_group_weight=float(c["mesospheric_group_weight"]),
            ionospheric_group_weight=float(c["ionospheric_group_weight"]),
            coef_err_sigma_floor_rel=c["coef_err_sigma_floor_rel"],
            zodi_ctx_restriction=c["zodi_ctx_restriction"],
            continuum_ctx_restriction=c["continuum_ctx_restriction"],
            alpha_ctx_features=c["alpha_ctx_features"],
            flux_mse_groups=tuple(c.get("flux_mse_groups", ())),
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
        palace_oh_suffix=None,
        palace_diffuse_suffix=None,
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
