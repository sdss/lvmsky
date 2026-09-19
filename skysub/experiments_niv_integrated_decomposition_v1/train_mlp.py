"""Train the existing coefficient-transfer MLP on the integrated Adam corpus."""

# ruff: noqa: E402 -- expose the historical top-level ``sky_decomp`` package.

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

# The extracted MLP package still imports ``sky_decomp`` as a top-level package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from skysub.mlp_predictor import (
    compressor,
    data,
    moon_model_cache,
    serialization,
    trainer,
    wavelengths,
)
from skysub.mlp_predictor.config import DEFAULT_CONTEXT_COLUMNS
from skysub.mlp_predictor.ml_utils import (
    moon_phase_deg_from_ctx,
    split_indices_by_moon_phase,
)
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_PALACE_DIFFUSE_SUFFIX,
    DEFAULT_PALACE_OH_SUFFIX,
    file_sha256,
)


SUFFIX = "_adam25k_telluric_niv_continuum"


def _paths(stack: Path, output_dir: Path) -> dict[str, Path]:
    stem = stack.stem
    return {
        role: output_dir / f"{stem}_{role}_meta_coef{SUFFIX}.fits"
        for role in ("sky1", "sky2", "sci")
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stack", type=Path)
    parser.add_argument("decomposition_dir", type=Path)
    parser.add_argument("--seeds", default="42,43,44,45,46,47,48,49,50,51")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference-checkpoint", type=Path)
    args = parser.parse_args()
    seeds = tuple(int(value) for value in args.seeds.split(","))
    if not seeds or args.epochs < 1:
        raise ValueError("At least one seed and one epoch are required")
    paths = _paths(args.stack, args.decomposition_dir)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing compact decomposition: " + ", ".join(missing))

    triplet = data.build_triplet_coef_dataset(
        input_fits_path=args.stack,
        sky_near_decomp_fits_path=paths["sky1"],
        sky_far_decomp_fits_path=paths["sky2"],
        sci_decomp_fits_path=paths["sci"],
        context_columns=list(DEFAULT_CONTEXT_COLUMNS),
        return_chi2=True,
    )
    filtered = data.apply_triplet_filters(
        triplet,
        thin_every_n=1,
        chi2_qmax=90.0,
        chi2_min=0.0,
        chi2_max=10.0,
        hard_coef_bounds={"feo": (0.0, 36.01), "atom_k": (0.0, 10.01)},
        kappa=8.0,
        kappa_iter=3,
        oh_kappa=6.0,
        oh_kappa_iter=3,
        exclude_field_regions=[data.LMC_EXCLUSION, data.SMC_EXCLUSION],
        colour_excess_input_fits=args.stack,
        colour_excess_max=data.SCI_COLOUR_EXCESS_MAX,
        diffuse_zeroed_frac=data.DIFFUSE_ZEROED_FRAC,
        show_plots=False,
    )
    if args.max_rows is not None and len(filtered["row_index"]) > args.max_rows:
        n_before = len(filtered["row_index"])
        keep = np.rint(
            np.linspace(0, n_before - 1, args.max_rows)
        ).astype(int)
        for key, value in list(filtered.items()):
            if isinstance(value, np.ndarray) and value.shape[:1] == (n_before,):
                filtered[key] = value[keep]
        filtered["n_rows"] = int(args.max_rows)

    data._augment_triplet_with_ecliptic(
        filtered, force=True, meta_fits_path=args.stack
    )
    data._augment_triplet_with_physics_priors(filtered, force=True)
    data._augment_triplet_with_moon_model(
        filtered, args.stack.with_suffix(""), force=True
    )
    extinction = wavelengths.resolve_wavelengths_and_extinction(
        filtered,
        input_fits_for_basis=args.stack,
        use_fitted_extinction=True,
        palace_oh_suffix=DEFAULT_PALACE_OH_SUFFIX,
        palace_diffuse_suffix=DEFAULT_PALACE_DIFFUSE_SUFFIX,
        verbose=True,
    )
    group_indices = extinction.group_indices
    n_moon, split_zodi, n_zodi = wavelengths.infer_spline_knots(
        filtered["coef_names"]
    )
    split_train, split_val, split_test = split_indices_by_moon_phase(
        filtered["obstime_mjd"], moon_phase_deg_from_ctx(filtered), seed=42
    )
    group_compressors, geom_kwargs = compressor.fit_all_group_compressors(
        filtered,
        group_indices,
        train_idx=split_train,
        held_idx=split_val,
        xarm_threshold=compressor.COMPRESSION_XARM_THRESHOLD,
        verbose=True,
    )
    filtered["compress_train_idx"] = split_train
    filtered["compress_val_idx"] = split_val
    filtered["compress_test_idx"] = split_test

    train_config = dict(trainer.default_dual_group_config)
    reference = None
    if args.reference_checkpoint is not None:
        old = serialization.load_ensemble(args.reference_checkpoint, device="cpu")
        old_config = dict(old["config"])
        for key in train_config:
            if key in old_config:
                train_config[key] = old_config[key]
        if "flux_pixel_weighted" in old_config:
            train_config["flux_pixel_weighting"] = bool(
                old_config["flux_pixel_weighted"]
            )
        reference = {
            "path": str(args.reference_checkpoint.resolve()),
            "sha256": file_sha256(args.reference_checkpoint),
            "seeds": list(old["seeds"]),
            "coefficient_schema_matches": list(old["coef_names"])
            == list(filtered["coef_names"]),
            "context_schema_matches": list(old["ctx_names"])
            == list(filtered["ctx_names"]),
            "training_config_reused": True,
        }
        if not (
            reference["coefficient_schema_matches"]
            and reference["context_schema_matches"]
        ):
            raise ValueError(
                "Reference checkpoint coefficient/context schema does not "
                "match the integrated Adam training corpus"
            )
        del old
    train_config["ensemble_seeds"] = seeds
    train_config["n_epochs"] = args.epochs
    artifacts = trainer.Trainer(cfg=train_config).run_ensemble(
        filtered,
        group_compressors,
        group_indices,
        geom_kwargs,
        input_fits_for_basis=args.stack,
        input_fits_flux=args.stack,
        n_moon_knots=n_moon,
        split_zodi=split_zodi,
        n_zodi_knots=n_zodi,
        palace_oh_suffix=DEFAULT_PALACE_OH_SUFFIX,
        palace_diffuse_suffix=DEFAULT_PALACE_DIFFUSE_SUFFIX,
        verbose=True,
    )
    output = args.output or args.decomposition_dir / "mlp_ensemble_niv_adam_v1.pt"
    serialization.save_ensemble(artifacts.mlp_artifacts, output)
    artifacts.per_seed_test_metrics.to_csv(
        args.decomposition_dir / "mlp_test_metrics_niv_adam_v1.csv", index=False
    )

    test_idx = np.asarray(artifacts.mlp_artifacts["test_idx"], dtype=int)
    prediction = trainer.predict_sci_coefficients_default(
        artifacts.mlp_artifacts,
        coef_near_phys=filtered["coef_near"][test_idx],
        coef_far_phys=filtered["coef_far"][test_idx],
        ctx_near_phys=filtered["ctx_near"][test_idx],
        ctx_far_phys=filtered["ctx_far"][test_idx],
        ctx_sci_phys=filtered["ctx_sci"][test_idx],
    ).astype(np.float32)
    example = int(test_idx[0])
    np.savez_compressed(
        args.decomposition_dir / "mlp_example_niv_adam_v1.npz",
        source_row=np.asarray(filtered["row_index"][example], dtype=np.int64),
        coef_names=np.asarray(filtered["coef_names"]),
        coef_near=np.asarray(filtered["coef_near"][example]),
        coef_far=np.asarray(filtered["coef_far"][example]),
        coef_true=np.asarray(filtered["coef_sci"][example]),
        coef_pred=np.asarray(prediction[0]),
    )
    summary = {
        "checkpoint": str(output.resolve()),
        "checkpoint_sha256": file_sha256(output),
        "input_stack": str(args.stack.resolve()),
        "input_stack_sha256": file_sha256(args.stack),
        "decomposition_sha256": {
            role: file_sha256(path) for role, path in paths.items()
        },
        "moon_model_cache": {
            "path": str(
                moon_model_cache.cache_path(args.stack.with_suffix("")).resolve()
            ),
            "sha256": file_sha256(
                moon_model_cache.cache_path(args.stack.with_suffix(""))
            ),
        },
        "solar_activity_source": {
            "path": str(data.SOLAR_ACTIVITY_RAW_PATH.resolve()),
            "sha256": file_sha256(data.SOLAR_ACTIVITY_RAW_PATH),
        },
        "input_rows": int(triplet["n_rows"]),
        "filtered_rows": int(len(filtered["row_index"])),
        "train_rows": int(split_train.size),
        "validation_rows": int(split_val.size),
        "test_rows": int(split_test.size),
        "seeds": list(seeds),
        "epochs": args.epochs,
        "native_wave_pixels": 12_401,
        "wavelength_stride": 1,
        "component_method": "adam25k-telluric-niv-continuum",
        "pca30_used": False,
        "compact_component_cube_filters_skipped": ["moon_zodi_reversal"],
        "reference_niv_checkpoint": reference,
        "training_config": train_config,
        "test_metrics": artifacts.per_seed_test_metrics.to_dict("records"),
    }
    (args.decomposition_dir / "mlp_training_summary_niv_adam_v1.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
