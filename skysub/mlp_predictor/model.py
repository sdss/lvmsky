"""Symmetric dual-encoder group-head MLP.

Verbatim extraction of the model definition from cell ``5dbf36fe`` of the
split-zodi notebook.  The ``VAN_RHIJN_FEATURES`` constant lives here because
this is the module that decides whether to drop those columns from the
context vector; downstream ``data`` can import from this module.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F

VAN_RHIJN_FEATURES: frozenset[str] = frozenset({
    "vanrhijn_87km",
    "vanrhijn_95km",
    "vanrhijn_285km",
})


class DualEncoderGroupHeadMLP(nn.Module):
    """Same-time SCI coefficient predictor.

    Architecture (unchanged from the deployed split-zodi variant):
    - a shared MLP encodes each sky arm's ``(coef, ctx)`` concatenation;
    - the two arm embeddings are combined as ``(mean, diff, |diff|)``;
    - a ctx-only MLP encodes the science pointing's context vector;
    - a trunk MLP fuses the four embeddings and feeds per-group heads,
      each producing softplus-positive intrinsic emissivities.
    """

    def __init__(
        self,
        n_coef: int,
        n_ctx: int,
        group_indices: Mapping[str, Sequence[int]],
        ctx_names: Sequence[str] | None = None,
        encoder_dims: tuple[int, ...] = (384, 192),
        ctx_dims: tuple[int, ...] = (64,),
        trunk_dims: tuple[int, ...] = (256, 128),
        head_dim: int = 128,
        head_extra_dims: tuple[int, ...] = (),
        drop_vanrhijn_from_context: bool = False,
    ) -> None:
        super().__init__()
        self.group_indices = group_indices
        self.ctx_names = (
            [str(n).strip().lower() for n in ctx_names] if ctx_names is not None else None
        )

        # Airglow geometry is divided out analytically before anything reaches
        # this model, so the vanrhijn_* columns are redundant FOR THE AIRGLOW
        # GROUPS.  They are not redundant for the moon and continuum groups,
        # which receive no analytic correction and must learn their airmass
        # dependence from context — and the encoders are shared across groups,
        # so dropping the columns takes those smooth non-linear functions of
        # altitude away from every head.  Default is therefore to keep them.
        if self.ctx_names is not None and bool(drop_vanrhijn_from_context):
            self.ctx_keep_idx = [
                i for i, name in enumerate(self.ctx_names)
                if name not in VAN_RHIJN_FEATURES
            ]
        elif self.ctx_names is not None:
            self.ctx_keep_idx = list(range(len(self.ctx_names)))
        else:
            self.ctx_keep_idx = list(range(int(n_ctx)))

        self.ctx_dim = len(self.ctx_keep_idx)
        self.sky_encoder = self._make_mlp(n_coef + self.ctx_dim, encoder_dims)
        sky_embed_dim = int(encoder_dims[-1])

        self.ctx_encoder = self._make_mlp(self.ctx_dim, ctx_dims)
        ctx_embed_dim = int(ctx_dims[-1]) if len(ctx_dims) > 0 else self.ctx_dim

        fusion_in = 3 * sky_embed_dim + ctx_embed_dim
        self.trunk = self._make_mlp(fusion_in, trunk_dims)
        trunk_out = int(trunk_dims[-1]) if len(trunk_dims) > 0 else fusion_in

        self.heads = nn.ModuleDict({
            g: self._make_group_head(trunk_out, int(head_dim),
                                     tuple(int(v) for v in head_extra_dims),
                                     len(idx))
            for g, idx in group_indices.items()
        })

    @staticmethod
    def _make_group_head(in_dim: int, head_dim: int,
                         extra_dims: Sequence[int], n_out: int) -> nn.Module:
        """trunk_out -> head_dim -> [extra_dims...] -> n_out with GELU between."""
        layers: list[nn.Module] = [nn.Linear(int(in_dim), int(head_dim)), nn.GELU()]
        d0 = int(head_dim)
        for d1 in extra_dims:
            d1 = int(d1)
            layers += [nn.Linear(d0, d1), nn.GELU()]
            d0 = d1
        layers += [nn.Linear(d0, int(n_out))]
        return nn.Sequential(*layers)

    @staticmethod
    def _make_mlp(in_dim: int, dims: Sequence[int]) -> nn.Module:
        layers: list[nn.Module] = []
        d0 = int(in_dim)
        for d1 in dims:
            d1 = int(d1)
            layers += [nn.Linear(d0, d1), nn.LayerNorm(d1), nn.GELU()]
            d0 = d1
        return nn.Sequential(*layers) if layers else nn.Identity()

    def _select_ctx(self, ctx: torch.Tensor) -> torch.Tensor:
        if self.ctx_names is None or len(self.ctx_keep_idx) == ctx.shape[1]:
            return ctx
        if len(self.ctx_keep_idx) == 0:
            return ctx[:, :0]
        return ctx[:, self.ctx_keep_idx]

    def forward(
        self,
        near_coef: torch.Tensor,
        far_coef: torch.Tensor,
        near_ctx: torch.Tensor,
        far_ctx: torch.Tensor,
        sci_ctx: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        near_ctx_model = self._select_ctx(near_ctx)
        far_ctx_model = self._select_ctx(far_ctx)
        sci_ctx_model = self._select_ctx(sci_ctx)

        e_near = self.sky_encoder(torch.cat([near_coef, near_ctx_model], dim=1))
        e_far = self.sky_encoder(torch.cat([far_coef, far_ctx_model], dim=1))
        e_mean = 0.5 * (e_near + e_far)
        e_diff = e_near - e_far
        e_absdiff = torch.abs(e_diff)
        e_ctx = self.ctx_encoder(sci_ctx_model)

        z = torch.cat([e_mean, e_diff, e_absdiff, e_ctx], dim=1)
        h = self.trunk(z)
        # Intrinsic emissivities in model space; the caller multiplies by
        # ``airglow_geometry_scale(ctx_sci)`` to get observed amplitudes.
        return {g: F.softplus(head(h)) for g, head in self.heads.items()}


import numpy as np

# Compressed dual-encoder group-head MLP.  THIS is the notebook's default model.
#
# See §5.5 of the methods cell for the compression rationale and the empirical
# argument for switching from the uncompressed baseline of §6 to this one.
# The uncompressed model is kept only as the parent class of
# `DualEncoderGroupHeadMLPCompressed`; its architecture, geometry handling,
# splitting and optimiser are inherited unchanged.
#
# The compressed model differs from the parent baseline in three ways:
#   * input to the encoder is per-arm compressed scores (n_input_score) rather
#     than sqrt-scaled coefficients (n_coef).
#   * group heads emit SIGNED score-space predictions (linear, no Softplus);
#     non-negativity is recovered by clipping the reconstructed physical
#     emissivity after inverse projection through the per-group compressor.
#   * loss is SmoothL1 on scaled score-space vectors, per group, with the
#     1/sqrt(n_g_score) balancing described in §7.  Because compressed groups
#     have comparable n_g_score, the per-coefficient weight imbalance between
#     continuum and OH of the uncompressed loss collapses by an order of
#     magnitude. Both `moon_group_weight` (default 3.0) and
#     `continuum_group_weight` (default 1.0) are exposed as per-group
#     multipliers $m_g$ on top of $1/\sqrt{n_g}$. `moon_smooth_lambda`
#     is retired (subsumed by the PCA truncation on the moon spline).
# ============================================================================


class DualEncoderGroupHeadMLPCompressed(DualEncoderGroupHeadMLP):
    """Subclass of the baseline arch with per-group heads emitting signed scores."""

    def __init__(self, *, n_score, n_ctx, group_score_dims, ctx_names=None,
                 encoder_dims=(384, 192), ctx_dims=(64,),
                 trunk_dims=(320, 160), head_dim=192,
                 head_extra_dims=(),
                 zodi_head_extra_dims=(),
                 continuum_head_extra_dims=(),
                 continuum_branch_dims=(64, 32),
                 drop_vanrhijn_from_context=False,
                 blend_init_alpha=0.7,
                 blend_use_direct=False,
                 moon_alt_conditional_alpha=False,
                 alpha_ctx_features=None,
                 alpha_ctx_groups=None,
                 zodi_ctx_restriction=None,
                 continuum_ctx_restriction=None,
                 moon_zodi_ctx_restriction=None,
                 moon_zodi_branch_dims=(128, 64),
                 moon_zodi_moon_head_extra_dims=(64,),
                 moon_zodi_zodi_head_extra_dims=(32,),
                 moon_zodi_mode='additive',
                 moon_zodi_coupling_dims=(64, 32)):
        fake_group_indices = {g: np.arange(int(n))
                              for g, n in group_score_dims.items()}
        super().__init__(
            n_coef=int(n_score),
            n_ctx=int(n_ctx),
            group_indices=fake_group_indices,
            ctx_names=ctx_names,
            encoder_dims=encoder_dims,
            ctx_dims=ctx_dims,
            trunk_dims=trunk_dims,
            head_dim=head_dim,
            head_extra_dims=tuple(int(v) for v in head_extra_dims),
            drop_vanrhijn_from_context=drop_vanrhijn_from_context,
        )
        self.group_score_dims = dict(group_score_dims)
        # Cumulative offsets track how compress_coefs_to_scores concatenates the per-group score blocks.
        _offsets = np.cumsum([0] + [int(group_score_dims[g]) for g in group_score_dims])
        self._score_offsets = {g: (int(_offsets[i]), int(_offsets[i + 1]))
                               for i, g in enumerate(group_score_dims)}
        # Per-group learnable blend weight. Two parametrizations:
        #   blend_use_direct=False (default) -> alpha_g = sigmoid(logit_g).
        #   blend_use_direct=True            -> alpha_g stored directly, clamped in [eps, 1-eps] after each opt step.
        self.blend_use_direct = bool(blend_use_direct)
        self.blend_alpha_eps = 1e-3
        # Moon-alt-conditional alpha: 'moon' blend uses a per-row 2-value lookup
        # selected by the sign of moon_alt at the science pointing (2026-08-12).
        self.moon_alt_conditional_alpha = bool(moon_alt_conditional_alpha)
        _init_alpha = float(np.clip(blend_init_alpha, self.blend_alpha_eps, 1 - self.blend_alpha_eps))
        _init_logit = float(np.log(_init_alpha / (1.0 - _init_alpha)))
        _scalar_alpha_groups = [str(g) for g in group_score_dims
                                if not (self.moon_alt_conditional_alpha and g == 'moon')]
        if self.blend_use_direct:
            self.blend_alpha_direct = nn.ParameterDict({
                str(g): nn.Parameter(torch.tensor(_init_alpha))
                for g in _scalar_alpha_groups
            })
            if self.moon_alt_conditional_alpha and 'moon' in group_score_dims:
                # Index 0 = moon_alt <= 0 (moon down); 1 = moon_alt > 0 (moon up).
                self.blend_alpha_moon_by_alt = nn.Parameter(
                    torch.tensor([_init_alpha, _init_alpha]))
        else:
            self.blend_logit = nn.ParameterDict({
                str(g): nn.Parameter(torch.tensor(_init_logit))
                for g in _scalar_alpha_groups
            })
            if self.moon_alt_conditional_alpha and 'moon' in group_score_dims:
                self.blend_moon_by_alt_logit = nn.Parameter(
                    torch.tensor([_init_logit, _init_logit]))

        # -----------------------------------------------------------
        # Phase A (2026-08-26): context-dependent per-group alpha.
        # When `alpha_ctx_features` is a non-empty iterable of ctx feature
        # names, each group's blend is computed as
        #     alpha_g(x_sci) = sigmoid(w_g . x_sci_alpha + b_g)
        # per row, where x_sci_alpha slices the requested features from the
        # science-arm ctx.  The linear predictors are initialised so that at
        # t=0 alpha_g == blend_init_alpha for every group, so training starts
        # byte-identical to the scalar-alpha path.  Overrides both the
        # scalar and moon_alt_conditional_alpha paths when active.
        # -----------------------------------------------------------
        self.alpha_ctx_features = None
        self.alpha_ctx_groups = None
        self.alpha_predictors = None
        if alpha_ctx_features and ctx_names is not None:
            _ctx_names_lower = [str(n).strip().lower() for n in ctx_names]
            _want = [str(x).strip().lower() for x in alpha_ctx_features]
            _idx = [i for i, n in enumerate(_ctx_names_lower) if n in _want]
            _missing = [w for w in _want if w not in _ctx_names_lower]
            if len(_idx) == 0:
                print('DualEncoderGroupHeadMLPCompressed: alpha_ctx_features '
                      f'requested {_want!r} but none of those features are in '
                      f'ctx_names; falling back to scalar alpha.')
            else:
                if _missing:
                    print('DualEncoderGroupHeadMLPCompressed: alpha_ctx_features '
                          f'missing features {_missing!r} (they will be skipped).')
                self.alpha_ctx_features = tuple(_ctx_names_lower[i] for i in _idx)
                self.register_buffer('alpha_ctx_idx',
                                     torch.tensor(_idx, dtype=torch.long))
                _n_alpha = len(_idx)
                # Which groups get a ctx-alpha predictor.  Default = all.
                if alpha_ctx_groups:
                    _groups_ctx = tuple(str(g) for g in alpha_ctx_groups
                                         if str(g) in group_score_dims)
                else:
                    _groups_ctx = tuple(str(g) for g in group_score_dims)
                self.alpha_ctx_groups = _groups_ctx
                self.alpha_predictors = nn.ModuleDict({
                    str(g): nn.Linear(_n_alpha, 1)
                    for g in _groups_ctx
                })
                # Zero the weights and set the bias so that sigmoid(bias)=init_alpha.
                with torch.no_grad():
                    for _p in self.alpha_predictors.values():
                        _p.weight.zero_()
                        _p.bias.fill_(_init_logit)
                print(f'DualEncoderGroupHeadMLPCompressed: context-dependent alpha active, '
                      f'features = {self.alpha_ctx_features}, groups = {self.alpha_ctx_groups}.')

        # -----------------------------------------------------------
        # Optional isolated zodi branch (physical-prior architecture).
        # When `zodi_ctx_restriction` is a non-empty iterable of ctx
        # feature names, the zodi head is routed through a dedicated
        # sub-network that sees ONLY those ctx features plus the two
        # arm zodi score blocks. This enforces the physical prior that
        # zodi brightness depends on ecliptic geometry and airmass,
        # not on moon phase / geomagnetic activity / time-of-night.
        # See §12 (zodi context restriction) in the intro cell.
        # -----------------------------------------------------------
        self.zodi_ctx_restriction = None
        self.zodi_branch = None
        self.zodi_head_isolated = None
        if zodi_ctx_restriction and 'zodi' in group_score_dims and ctx_names is not None:
            _ctx_names_lower = [str(n).strip().lower() for n in ctx_names]
            _want = [str(x).strip().lower() for x in zodi_ctx_restriction]
            _idx = [i for i, n in enumerate(_ctx_names_lower) if n in _want]
            _missing = [w for w in _want if w not in _ctx_names_lower]
            if len(_idx) == 0:
                print('DualEncoderGroupHeadMLPCompressed: zodi_ctx_restriction '
                      f'requested {_want!r} but none of those features are '
                      f'in ctx_names; falling back to shared trunk for zodi.')
            else:
                if _missing:
                    print('DualEncoderGroupHeadMLPCompressed: zodi_ctx_restriction '
                          f'missing features {_missing!r} (they will be skipped).')
                self.zodi_ctx_restriction = tuple(_ctx_names_lower[i] for i in _idx)
                self.register_buffer('zodi_ctx_idx',
                                     torch.tensor(_idx, dtype=torch.long))
                _n_zodi_score = int(group_score_dims['zodi'])
                # sci + near + far each contribute `len(_idx)` restricted ctx features;
                # the two arm zodi score blocks provide the per-arm coefficient signal.
                _zodi_in = 3 * len(_idx) + 2 * _n_zodi_score
                _hidden_a, _hidden_b = 64, 32
                self.zodi_branch = nn.Sequential(
                    nn.Linear(_zodi_in, _hidden_a),
                    nn.LayerNorm(_hidden_a), nn.GELU(),
                    nn.Linear(_hidden_a, _hidden_b),
                    nn.LayerNorm(_hidden_b), nn.GELU(),
                )
                _zh_extras = tuple(int(v) for v in zodi_head_extra_dims)
                if _zh_extras:
                    _layers: list[nn.Module] = []
                    _din = int(_hidden_b)
                    for _dh in _zh_extras:
                        _layers += [nn.Linear(_din, _dh), nn.GELU()]
                        _din = _dh
                    _layers += [nn.Linear(_din, _n_zodi_score)]
                    self.zodi_head_isolated = nn.Sequential(*_layers)
                else:
                    self.zodi_head_isolated = nn.Linear(_hidden_b, _n_zodi_score)
                print(f'DualEncoderGroupHeadMLPCompressed: isolated zodi branch active, '
                      f'ctx features = {self.zodi_ctx_restriction} (fed per-arm sci+near+far), '
                      f'n_score={_n_zodi_score}, head input dim={_zodi_in}, '
                      f'zodi_head_extra_dims={_zh_extras}.')

        # -----------------------------------------------------------
        # Optional isolated continuum branch (moon-regime conditioning).
        # Same architecture as the zodi branch above.  Phase B (2026-08-26)
        # attacks the ~0.70 residual RF R^2 on continuum by giving the head
        # direct access to moon geometry, which the shared trunk squeezes.
        # -----------------------------------------------------------
        self.continuum_ctx_restriction = None
        self.continuum_branch = None
        self.continuum_head_isolated = None
        if continuum_ctx_restriction and 'continuum' in group_score_dims and ctx_names is not None:
            _ctx_names_lower = [str(n).strip().lower() for n in ctx_names]
            _want = [str(x).strip().lower() for x in continuum_ctx_restriction]
            _idx = [i for i, n in enumerate(_ctx_names_lower) if n in _want]
            _missing = [w for w in _want if w not in _ctx_names_lower]
            if len(_idx) == 0:
                print('DualEncoderGroupHeadMLPCompressed: continuum_ctx_restriction '
                      f'requested {_want!r} but none of those features are '
                      f'in ctx_names; falling back to shared trunk for continuum.')
            else:
                if _missing:
                    print('DualEncoderGroupHeadMLPCompressed: continuum_ctx_restriction '
                          f'missing features {_missing!r} (they will be skipped).')
                self.continuum_ctx_restriction = tuple(_ctx_names_lower[i] for i in _idx)
                self.register_buffer('continuum_ctx_idx',
                                     torch.tensor(_idx, dtype=torch.long))
                _n_c_score = int(group_score_dims['continuum'])
                _c_in = 3 * len(_idx) + 2 * _n_c_score
                _cb_dims = tuple(int(v) for v in continuum_branch_dims)
                if len(_cb_dims) != 2:
                    raise ValueError(
                        f'continuum_branch_dims must be a 2-tuple (hidden_a, hidden_b); got {_cb_dims!r}.')
                _hidden_a, _hidden_b = _cb_dims
                self.continuum_branch = nn.Sequential(
                    nn.Linear(_c_in, _hidden_a),
                    nn.LayerNorm(_hidden_a), nn.GELU(),
                    nn.Linear(_hidden_a, _hidden_b),
                    nn.LayerNorm(_hidden_b), nn.GELU(),
                )
                _ch_extras = tuple(int(v) for v in continuum_head_extra_dims)
                if _ch_extras:
                    _layers = []
                    _din = int(_hidden_b)
                    for _dh in _ch_extras:
                        _layers += [nn.Linear(_din, _dh), nn.GELU()]
                        _din = _dh
                    _layers += [nn.Linear(_din, _n_c_score)]
                    self.continuum_head_isolated = nn.Sequential(*_layers)
                else:
                    self.continuum_head_isolated = nn.Linear(_hidden_b, _n_c_score)
                print(f'DualEncoderGroupHeadMLPCompressed: isolated continuum branch active, '
                      f'ctx features = {self.continuum_ctx_restriction} (fed per-arm sci+near+far), '
                      f'n_score={_n_c_score}, head input dim={_c_in}, '
                      f'continuum_branch_dims={_cb_dims}, '
                      f'continuum_head_extra_dims={_ch_extras}.')

        # -----------------------------------------------------------
        # 2026-08-27b (Phase E) / 2026-08-27c (Phase F): moon-zodi coupling.
        # When `moon_zodi_ctx_restriction` is a non-empty iterable of ctx feature
        # names, `moon_zodi_mode` selects between:
        #   'shared_branch' -- Phase E: retire the shared-trunk moon head AND
        #     the isolated zodi branch; moon and zodi both consume the shared
        #     `moon_zodi_branch` latent via new dedicated heads.  A/B on this
        #     path (2026-08-27b) landed +moon wins but atlas mid regressed +35%.
        #   'additive' (default) -- Phase F: keep both existing paths (moon on
        #     shared trunk, zodi on isolated branch).  Add a small coupling
        #     branch producing a shared latent, projected additively into the
        #     moon and zodi head outputs via zero-init linear projectors.  At
        #     init the projector outputs are zero, so training starts byte-
        #     identical to Phase A''; the projectors learn to route the
        #     moon-zodi coupling signal without breaking zodi's isolation.
        # -----------------------------------------------------------
        self.moon_zodi_ctx_restriction = None
        self.moon_zodi_mode = None
        self.moon_zodi_branch = None
        self.moon_zodi_moon_head = None
        self.moon_zodi_zodi_head = None
        self.moon_zodi_coupling_branch = None
        self.moon_zodi_moon_projector = None
        self.moon_zodi_zodi_projector = None
        if (moon_zodi_ctx_restriction and 'moon' in group_score_dims
                and 'zodi' in group_score_dims and ctx_names is not None):
            _mz_mode = str(moon_zodi_mode).lower() if moon_zodi_mode else 'additive'
            if _mz_mode not in ('shared_branch', 'additive'):
                raise ValueError(
                    f'moon_zodi_mode must be "additive" or "shared_branch"; got {_mz_mode!r}.')
            _ctx_names_lower = [str(n).strip().lower() for n in ctx_names]
            _want = [str(x).strip().lower() for x in moon_zodi_ctx_restriction]
            _idx = [i for i, n in enumerate(_ctx_names_lower) if n in _want]
            _missing = [w for w in _want if w not in _ctx_names_lower]
            if len(_idx) == 0:
                print('DualEncoderGroupHeadMLPCompressed: moon_zodi_ctx_restriction '
                      f'requested {_want!r} but none of those features are '
                      f'in ctx_names; falling back to per-group paths for moon+zodi.')
            else:
                if _missing:
                    print('DualEncoderGroupHeadMLPCompressed: moon_zodi_ctx_restriction '
                          f'missing features {_missing!r} (they will be skipped).')
                self.moon_zodi_ctx_restriction = tuple(_ctx_names_lower[i] for i in _idx)
                self.moon_zodi_mode = _mz_mode
                self.register_buffer('moon_zodi_ctx_idx',
                                     torch.tensor(_idx, dtype=torch.long))
                _n_m = int(group_score_dims['moon'])
                _n_z = int(group_score_dims['zodi'])
                _mz_in = 3 * len(_idx) + 2 * (_n_m + _n_z)

                def _build_head(_n_out, _extras, _in_dim):
                    _extras = tuple(int(v) for v in _extras)
                    if not _extras:
                        return nn.Linear(_in_dim, _n_out)
                    _layers = []
                    _din = int(_in_dim)
                    for _dh in _extras:
                        _layers += [nn.Linear(_din, _dh), nn.GELU()]
                        _din = _dh
                    _layers += [nn.Linear(_din, _n_out)]
                    return nn.Sequential(*_layers)

                if _mz_mode == 'shared_branch':
                    _mzb_dims = tuple(int(v) for v in moon_zodi_branch_dims)
                    if len(_mzb_dims) != 2:
                        raise ValueError(
                            f'moon_zodi_branch_dims must be a 2-tuple; got {_mzb_dims!r}.')
                    _hidden_a, _hidden_b = _mzb_dims
                    self.moon_zodi_branch = nn.Sequential(
                        nn.Linear(_mz_in, _hidden_a),
                        nn.LayerNorm(_hidden_a), nn.GELU(),
                        nn.Linear(_hidden_a, _hidden_b),
                        nn.LayerNorm(_hidden_b), nn.GELU(),
                    )
                    _mh_extras = tuple(int(v) for v in moon_zodi_moon_head_extra_dims)
                    _zh_extras = tuple(int(v) for v in moon_zodi_zodi_head_extra_dims)
                    self.moon_zodi_moon_head = _build_head(_n_m, _mh_extras, _hidden_b)
                    self.moon_zodi_zodi_head = _build_head(_n_z, _zh_extras, _hidden_b)
                    # Precedence: shared branch overrides isolated zodi.
                    if self.zodi_branch is not None or self.zodi_head_isolated is not None:
                        print('DualEncoderGroupHeadMLPCompressed: moon_zodi shared branch '
                              'takes precedence; disabling the pre-existing isolated zodi branch.')
                        self.zodi_branch = None
                        self.zodi_head_isolated = None
                        self.zodi_ctx_restriction = None
                    print(f'DualEncoderGroupHeadMLPCompressed: shared moon-zodi branch active, '
                          f'ctx features = {self.moon_zodi_ctx_restriction} (fed per-arm sci+near+far), '
                          f'n_moon_score={_n_m}, n_zodi_score={_n_z}, branch input dim={_mz_in}, '
                          f'moon_zodi_branch_dims={_mzb_dims}, '
                          f'moon_head_extra_dims={_mh_extras}, zodi_head_extra_dims={_zh_extras}.')
                else:
                    # additive mode: small coupling branch + zero-init projectors.
                    _cdims = tuple(int(v) for v in moon_zodi_coupling_dims)
                    if len(_cdims) != 2:
                        raise ValueError(
                            f'moon_zodi_coupling_dims must be a 2-tuple; got {_cdims!r}.')
                    _c_hidden_a, _c_hidden_b = _cdims
                    self.moon_zodi_coupling_branch = nn.Sequential(
                        nn.Linear(_mz_in, _c_hidden_a),
                        nn.LayerNorm(_c_hidden_a), nn.GELU(),
                        nn.Linear(_c_hidden_a, _c_hidden_b),
                        nn.LayerNorm(_c_hidden_b), nn.GELU(),
                    )
                    self.moon_zodi_moon_projector = nn.Linear(_c_hidden_b, _n_m)
                    self.moon_zodi_zodi_projector = nn.Linear(_c_hidden_b, _n_z)
                    # Zero the projector weights and biases so at init the
                    # additive contribution is exactly zero and training starts
                    # byte-identical to Phase A''.
                    with torch.no_grad():
                        for _p in (self.moon_zodi_moon_projector,
                                   self.moon_zodi_zodi_projector):
                            _p.weight.zero_()
                            _p.bias.zero_()
                    print(f'DualEncoderGroupHeadMLPCompressed: additive moon-zodi coupling active, '
                          f'ctx features = {self.moon_zodi_ctx_restriction} (fed per-arm sci+near+far), '
                          f'n_moon_score={_n_m}, n_zodi_score={_n_z}, branch input dim={_mz_in}, '
                          f'coupling_dims={_cdims}, projectors zero-init.')

    def forward(self, near_score, far_score, near_ctx, far_ctx, sci_ctx,
                sci_moon_alt_sign=None):
        near_ctx_m = self._select_ctx(near_ctx)
        far_ctx_m = self._select_ctx(far_ctx)
        sci_ctx_m = self._select_ctx(sci_ctx)
        e_near = self.sky_encoder(torch.cat([near_score, near_ctx_m], dim=1))
        e_far = self.sky_encoder(torch.cat([far_score, far_ctx_m], dim=1))
        e_mean = 0.5 * (e_near + e_far)
        e_diff = e_near - e_far
        e_absdiff = torch.abs(e_diff)
        e_ctx = self.ctx_encoder(sci_ctx_m)
        z = torch.cat([e_mean, e_diff, e_absdiff, e_ctx], dim=1)
        h = self.trunk(z)
        # Predict a residual on top of an explicit per-group blend of the two sky arms
        # so delta=0 is the interpolation baseline rather than the training mean.
        # moon-zodi latent computed once per batch when either the shared-branch
        # (Phase E) or additive coupling (Phase F) is active.
        _mz_latent = None
        _mz_branch_active = self.moon_zodi_branch is not None
        _mz_coupling_active = self.moon_zodi_coupling_branch is not None
        if _mz_branch_active or _mz_coupling_active:
            _mzctx_sci  = sci_ctx[:,  self.moon_zodi_ctx_idx]
            _mzctx_near = near_ctx[:, self.moon_zodi_ctx_idx]
            _mzctx_far  = far_ctx[:,  self.moon_zodi_ctx_idx]
            _lo_m, _hi_m = self._score_offsets['moon']
            _lo_z, _hi_z = self._score_offsets['zodi']
            _mzin = torch.cat([
                _mzctx_sci, _mzctx_near, _mzctx_far,
                near_score[:, _lo_m:_hi_m], far_score[:, _lo_m:_hi_m],
                near_score[:, _lo_z:_hi_z], far_score[:, _lo_z:_hi_z],
            ], dim=1)
            _mz_latent = (self.moon_zodi_branch(_mzin) if _mz_branch_active
                          else self.moon_zodi_coupling_branch(_mzin))
        out = {}
        for g, head in self.heads.items():
            lo, hi = self._score_offsets[g]
            if (self.alpha_predictors is not None
                    and str(g) in self.alpha_predictors):
                # Phase A: alpha_g(sci_ctx) as (batch, 1), broadcast over n_g.
                _x_alpha = sci_ctx[:, self.alpha_ctx_idx]
                alpha = torch.sigmoid(self.alpha_predictors[str(g)](_x_alpha))
                blend = (alpha * near_score[:, lo:hi]
                         + (1.0 - alpha) * far_score[:, lo:hi])
            elif (self.moon_alt_conditional_alpha and g == 'moon'
                    and sci_moon_alt_sign is not None):
                if self.blend_use_direct:
                    _alpha_row = self.blend_alpha_moon_by_alt[sci_moon_alt_sign]
                else:
                    _alpha_row = torch.sigmoid(
                        self.blend_moon_by_alt_logit[sci_moon_alt_sign])
                blend = (_alpha_row.unsqueeze(1) * near_score[:, lo:hi]
                         + (1.0 - _alpha_row.unsqueeze(1)) * far_score[:, lo:hi])
            else:
                if self.blend_use_direct:
                    alpha = self.blend_alpha_direct[str(g)]
                else:
                    alpha = torch.sigmoid(self.blend_logit[str(g)])
                blend = (alpha * near_score[:, lo:hi]
                         + (1.0 - alpha) * far_score[:, lo:hi])
            if (_mz_latent is not None and g == 'moon'
                    and self.moon_zodi_moon_head is not None):
                # Phase E shared_branch: moon head is replaced.
                out[g] = blend + self.moon_zodi_moon_head(_mz_latent)
            elif (_mz_latent is not None and g == 'zodi'
                    and self.moon_zodi_zodi_head is not None):
                # Phase E shared_branch: zodi head is replaced.
                out[g] = blend + self.moon_zodi_zodi_head(_mz_latent)
            elif (g == 'zodi' and self.zodi_branch is not None
                    and self.zodi_head_isolated is not None):
                _zctx_sci  = sci_ctx[:,  self.zodi_ctx_idx]
                _zctx_near = near_ctx[:, self.zodi_ctx_idx]
                _zctx_far  = far_ctx[:,  self.zodi_ctx_idx]
                _zin = torch.cat([_zctx_sci, _zctx_near, _zctx_far,
                                  near_score[:, lo:hi],
                                  far_score[:, lo:hi]], dim=1)
                _z_out = self.zodi_head_isolated(self.zodi_branch(_zin))
                if (_mz_latent is not None
                        and self.moon_zodi_zodi_projector is not None):
                    _z_out = _z_out + self.moon_zodi_zodi_projector(_mz_latent)
                out[g] = blend + _z_out
            elif (g == 'continuum' and self.continuum_branch is not None
                    and self.continuum_head_isolated is not None):
                _cctx_sci  = sci_ctx[:,  self.continuum_ctx_idx]
                _cctx_near = near_ctx[:, self.continuum_ctx_idx]
                _cctx_far  = far_ctx[:,  self.continuum_ctx_idx]
                _cin = torch.cat([_cctx_sci, _cctx_near, _cctx_far,
                                  near_score[:, lo:hi],
                                  far_score[:, lo:hi]], dim=1)
                out[g] = blend + self.continuum_head_isolated(self.continuum_branch(_cin))
            else:
                _head_out = head(h)
                if (g == 'moon' and _mz_latent is not None
                        and self.moon_zodi_moon_projector is not None):
                    _head_out = _head_out + self.moon_zodi_moon_projector(_mz_latent)
                out[g] = blend + _head_out
        return out


__all__ = [
    "DualEncoderGroupHeadMLP",
    "DualEncoderGroupHeadMLPCompressed",
    "VAN_RHIJN_FEATURES",
]
