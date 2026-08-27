"""Serialize / deserialize trained ensembles for inference.

The training pipeline in ``mlp_predictor.trainer`` produces an
``EnsembleArtifacts`` dataclass whose ``mlp_artifacts`` dict carries the
per-seed torch models, the RobustScaler instances for scores + ctx, the
per-group compressors, the Jensen post-training lift, and every downstream
knob the inference path needs.  This module writes those objects to a single
``.pt`` archive on disk and reads them back into a shape that is drop-in
compatible with :func:`mlp_predictor.trainer.predict_sci_coefficients_default`
and :func:`mlp_predictor.inference.predict_sky_from_minimal_inputs`.

Only inference-necessary state is persisted; training-time bookkeeping
(``train_idx``/``val_idx``/``test_idx``, history, blend_history, best epoch)
is intentionally dropped so the saved file stays small and doesn't leak the
training split.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Any

import torch

from .model import DualEncoderGroupHeadMLPCompressed

FORMAT_VERSION = 1


def _cpu_state_dict(model: torch.nn.Module) -> dict:
    """Detach + move a model state_dict to CPU tensors for portable serialisation."""
    return {k: v.detach().cpu() for k, v in model.state_dict().items()}


def save_ensemble(ensemble_artifacts: Mapping[str, Any], out_path: str | Path) -> Path:
    """Persist an ensemble artifact to a single ``.pt`` archive.

    Parameters
    ----------
    ensemble_artifacts:
        The ``mlp_artifacts`` dict produced by :class:`~mlp_predictor.trainer.Trainer`
        (equivalently ``EnsembleArtifacts.mlp_artifacts``).  Must have
        ``is_ensemble=True`` and a ``members`` list.
    out_path:
        Destination file path.  Parent directories are created automatically.

    Returns
    -------
    Path
        The absolute path to the written archive.
    """
    if not ensemble_artifacts.get("is_ensemble", False):
        raise ValueError(
            "save_ensemble expects an ensemble artifact "
            "(is_ensemble=True); got a single-seed artifact.")
    members = ensemble_artifacts["members"]
    if not members:
        raise ValueError("Ensemble has no members to save.")

    first = members[0]
    payload = {
        "format_version": FORMAT_VERSION,
        "config": dict(first["config"]),
        "seeds": list(ensemble_artifacts["seeds"]),
        "member_state_dicts": [_cpu_state_dict(m["model"]) for m in members],
        "score_scaler": first["score_scaler"],
        "ctx_scaler":   first["ctx_scaler"],
        "compressors":  ensemble_artifacts["compressors"],
        "jensen_corrections": first["jensen_corrections"],
        "coef_upper_bound":   ensemble_artifacts["coef_upper_bound"],
        "geom_kwargs":  ensemble_artifacts["geom_kwargs"],
        "group_indices":    ensemble_artifacts["group_indices"],
        "score_slices":     ensemble_artifacts["score_slices"],
        "group_score_dims": ensemble_artifacts["group_score_dims"],
        "n_input_score":    int(ensemble_artifacts["n_input_score"]),
        "coef_names": list(ensemble_artifacts["coef_names"]),
        "ctx_names":  list(ensemble_artifacts["ctx_names"]),
    }
    out_path = Path(out_path).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(payload, str(out_path))
    print(f"[save_ensemble] wrote {len(members)}-member ensemble to {out_path} "
          f"(format v{FORMAT_VERSION}).")
    return out_path


def _pick_device(explicit: str | None = None) -> str:
    if explicit:
        return str(explicit)
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def _build_model_from_config(cfg: Mapping[str, Any], *, n_input_score: int,
                             n_ctx: int, group_score_dims: Mapping[str, int],
                             ctx_names: list[str]) -> DualEncoderGroupHeadMLPCompressed:
    """Reconstruct a model instance with the architecture recorded at training time."""
    return DualEncoderGroupHeadMLPCompressed(
        n_score=int(n_input_score), n_ctx=int(n_ctx),
        group_score_dims=dict(group_score_dims),
        ctx_names=list(ctx_names),
        encoder_dims=tuple(int(v) for v in cfg.get("encoder_dims", (384, 192))),
        ctx_dims=tuple(int(v) for v in cfg.get("ctx_dims", (64,))),
        trunk_dims=tuple(int(v) for v in cfg.get("trunk_dims", (320, 160))),
        head_dim=int(cfg.get("head_dim", 192)),
        head_extra_dims=tuple(int(v) for v in cfg.get("head_extra_dims", ())),
        zodi_head_extra_dims=tuple(int(v) for v in cfg.get("zodi_head_extra_dims", ())),
        continuum_head_extra_dims=tuple(int(v) for v in cfg.get("continuum_head_extra_dims", ())),
        continuum_branch_dims=tuple(int(v) for v in cfg.get("continuum_branch_dims", (64, 32))),
        moon_zodi_ctx_restriction=cfg.get("moon_zodi_ctx_restriction"),
        moon_zodi_branch_dims=tuple(int(v) for v in cfg.get("moon_zodi_branch_dims", (128, 64))),
        moon_zodi_moon_head_extra_dims=tuple(int(v) for v in cfg.get("moon_zodi_moon_head_extra_dims", (64,))),
        moon_zodi_zodi_head_extra_dims=tuple(int(v) for v in cfg.get("moon_zodi_zodi_head_extra_dims", (32,))),
        moon_zodi_mode=str(cfg.get("moon_zodi_mode", "additive")),
        moon_zodi_coupling_dims=tuple(int(v) for v in cfg.get("moon_zodi_coupling_dims", (64, 32))),
        drop_vanrhijn_from_context=bool(cfg.get("drop_vanrhijn_from_context", False)),
        blend_init_alpha=float(cfg.get("blend_init_alpha", 0.7)),
        blend_use_direct=(str(cfg.get("blend_optim", "direct")).lower()
                          in ("direct", "prefit_freeze", "prefit_warmstart")),
        moon_alt_conditional_alpha=bool(cfg.get("moon_alt_conditional_alpha", False)),
        alpha_ctx_features=cfg.get("alpha_ctx_features"),
        alpha_ctx_groups=cfg.get("alpha_ctx_groups"),
        zodi_ctx_restriction=cfg.get("zodi_ctx_restriction"),
        continuum_ctx_restriction=cfg.get("continuum_ctx_restriction"),
    )


def load_ensemble(path: str | Path, *, device: str | None = None,
                  weights_only: bool = False) -> dict:
    """Load an ensemble saved with :func:`save_ensemble`.

    Parameters
    ----------
    path:
        Path to the ``.pt`` archive.
    device:
        Torch device string ('cpu', 'cuda', 'mps').  Defaults to the best
        available device.
    weights_only:
        Passed through to ``torch.load``.  Default False because the payload
        contains custom Python objects (``RobustScaler`` instances, dicts of
        numpy arrays) that ``torch.load(weights_only=True)`` would reject.

    Returns
    -------
    dict
        An ``mlp_artifacts`` dict shaped exactly like the one that
        :class:`~mlp_predictor.trainer.Trainer` returns.  Ready to hand to
        :func:`~mlp_predictor.trainer.predict_sci_coefficients_default`
        or :func:`~mlp_predictor.inference.predict_sky_from_minimal_inputs`.
    """
    dev = _pick_device(device)
    payload = torch.load(str(Path(path).expanduser()), map_location=dev,
                          weights_only=bool(weights_only))
    fv = int(payload.get("format_version", -1))
    if fv != FORMAT_VERSION:
        raise ValueError(f"Unsupported format_version {fv!r}; expected {FORMAT_VERSION}.")

    cfg = dict(payload["config"])
    ctx_names = list(payload["ctx_names"])
    group_score_dims = dict(payload["group_score_dims"])
    n_input_score = int(payload["n_input_score"])

    # Backwards-compat: ensembles saved before 2026-08-27 did not persist
    # `zodi_ctx_restriction` / `continuum_ctx_restriction` in their config
    # dict.  Recover them from the state_dict's `zodi_ctx_idx` /
    # `continuum_ctx_idx` buffers when present so the rebuild picks up the
    # same isolated branches the training run used.
    if payload["member_state_dicts"]:
        _sd0 = payload["member_state_dicts"][0]
        for _key, _cfg_key in (("zodi_ctx_idx", "zodi_ctx_restriction"),
                                ("continuum_ctx_idx", "continuum_ctx_restriction")):
            if cfg.get(_cfg_key) is None and _key in _sd0:
                _idx = _sd0[_key].detach().cpu().numpy().astype(int).tolist()
                cfg[_cfg_key] = tuple(str(ctx_names[i]).strip().lower() for i in _idx)

    members: list[dict] = []
    for sd in payload["member_state_dicts"]:
        model = _build_model_from_config(
            cfg,
            n_input_score=n_input_score,
            n_ctx=len(ctx_names),
            group_score_dims=group_score_dims,
            ctx_names=ctx_names,
        ).to(dev)
        model.load_state_dict(sd)
        model.eval()
        members.append({
            "model": model,
            "device": dev,
            "score_scaler": payload["score_scaler"],
            "ctx_scaler":   payload["ctx_scaler"],
            "compressors":  payload["compressors"],
            "jensen_corrections": payload["jensen_corrections"],
            "coef_upper_bound":   payload["coef_upper_bound"],
            "geom_kwargs":  payload["geom_kwargs"],
            "group_indices":    payload["group_indices"],
            "score_slices":     payload["score_slices"],
            "group_score_dims": group_score_dims,
            "n_input_score":    n_input_score,
            "coef_names": list(payload["coef_names"]),
            "ctx_names":  ctx_names,
            "config": cfg,
        })

    print(f"[load_ensemble] restored {len(members)}-member ensemble from {path} "
          f"onto device={dev}.")
    return {
        "is_ensemble": True,
        "seeds": list(payload["seeds"]),
        "members": members,
        "compressors":  payload["compressors"],
        "coef_upper_bound":   payload["coef_upper_bound"],
        "geom_kwargs":  payload["geom_kwargs"],
        "group_indices":    payload["group_indices"],
        "score_slices":     payload["score_slices"],
        "group_score_dims": group_score_dims,
        "n_input_score":    n_input_score,
        "coef_names": list(payload["coef_names"]),
        "ctx_names":  ctx_names,
        "config": cfg,
    }


__all__ = ["save_ensemble", "load_ensemble", "FORMAT_VERSION"]
