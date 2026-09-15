"""Train PCA30 above PALACE VNF with the frozen Niv continuum contract."""

# ruff: noqa: E402 -- thread limits must be set before importing NumPy.

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import time
import traceback

for name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
):
    os.environ[name] = "1"

import numpy as np
import pandas as pd
import scipy.sparse as sp
from astropy.io import fits
from scipy.sparse.linalg import lsqr

import skysub.decompose_parallel as decompose
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_DIFFUSE_SUFFIX,
    DEFAULT_PALACE_OH_SUFFIX,
    file_sha256,
    wave_sha256,
)
from skysub.sky_decomp.niv_continuum import NIV_CONTINUUM_KWARGS
from skysub.sky_decomp.residual_pca import (
    NIV_VNF_LINE_AMPLITUDE_PCA_ASSET,
    SkyDecompPalaceAijcVNFNivContinuumLSFSpline2D,
    _individual_line_design,
    _individual_line_names,
)


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "outputs" / "pca30"
SELECTION_PATH = OUTPUT_DIR / "selection_1000.csv"
CACHE_DIR = OUTPUT_DIR / "line_amplitudes"
PCA_ASSET = DEFAULT_DATA_ROOT / NIV_VNF_LINE_AMPLITUDE_PCA_ASSET
RIDGE_LAMBDA = 1.0e-4
OUTLIER_SIGMA = 7.0
PCA_COMPONENTS = 30
SOLVER_ID = "scipy.sparse.linalg.lsqr-unit-integral-line-ridge-v3"
ROW_ESTIMATOR_ID = "native-finite-mask-lsqr-v2"
_RUN_FINGERPRINT = None


def build_selection(stack_path: Path, count: int = 1000) -> pd.DataFrame:
    """Select evenly across all rows with valid native telluric metadata."""
    with fits.open(stack_path, memmap=True, lazy_load_hdus=True) as hdul:
        meta = hdul["META"].data
        labels = np.char.lower(np.char.strip(meta["sky_far_label"].astype(str)))
        source_airmass = np.where(
            labels == "skye", meta["skye_airmass"], meta["skyw_airmass"]
        ).astype(np.float64)
        valid = (
            np.isin(labels, ("skye", "skyw"))
            & np.isfinite(meta["pwv_med"])
            & (meta["pwv_med"] > 0.0)
            & np.isfinite(meta["sci_airmass"])
            & (meta["sci_airmass"] > 0.0)
            & np.isfinite(source_airmass)
            & (source_airmass > 0.0)
        )
        candidates = np.flatnonzero(valid)
        if candidates.size < count:
            raise ValueError(f"Only {candidates.size} valid rows are available")
        chosen = candidates[
            np.rint(np.linspace(0, candidates.size - 1, count)).astype(int)
        ]
        frame = pd.DataFrame(
            {
                "source_row": chosen,
                "expnum": np.asarray(meta["expnum"][chosen], dtype=np.int64),
                "mjd": np.asarray(meta["mjd"][chosen], dtype=np.int64),
                "pwv_mm": np.asarray(meta["pwv_med"][chosen], dtype=np.float64),
                "sci_airmass": np.asarray(
                    meta["sci_airmass"][chosen], dtype=np.float64
                ),
                "source_airmass": source_airmass[chosen],
                "sky_far_label": labels[chosen],
            }
        )
    if frame["source_row"].nunique() != count or frame["expnum"].nunique() != count:
        raise ValueError("The PCA selection must contain unique rows and exposures")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SELECTION_PATH, index=False)
    return frame


def _source_assets() -> dict[str, str]:
    relative = (
        f"palace/PMD/pmd_popmodel_OH{DEFAULT_PALACE_OH_SUFFIX}.dat",
        f"palace/PMD/pmd_refcont{DEFAULT_PALACE_DIFFUSE_SUFFIX}.dat",
        "palace/PMD/pmd_intdata_atom.dat",
        "palace/PMD/pmd_intmodel_Orc.dat",
        "palace/PMD/pmd_popmodel_O2.dat",
        "Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt",
        "moon_zodi/eso_skycalc_rolo_moon_albedo.dat",
    )
    return {name: file_sha256(DEFAULT_DATA_ROOT / name) for name in relative}


def _run_provenance(stack_path: Path, input_sha256: str) -> dict[str, object]:
    source_root = HERE.parent / "sky_decomp"
    payload = {
        "schema_version": 1,
        "model": "SkyDecompPalaceAijcVNFNivContinuumLSFSpline2D",
        "input_path": str(stack_path.resolve()),
        "input_sha256": input_sha256,
        "selection_sha256": file_sha256(SELECTION_PATH),
        "continuum_contract": NIV_CONTINUUM_KWARGS,
        "source_assets_sha256": _source_assets(),
        "source_code_sha256": {
            name: file_sha256(source_root / name)
            for name in (
                "fit.py",
                "lsf_spline2d.py",
                "lsf_surface_iterative.py",
                "niv_continuum.py",
                "residual_pca.py",
                "telluric_corrected_lines.py",
            )
        },
        "native_grid_only": True,
        "removed_or_replaced_wavelength_pixels": 0,
        "science_line_mask": "production IVAR=0 mask centred on measured Halpha",
        "ridge_lambda": RIDGE_LAMBDA,
        "amplitude_solver": SOLVER_ID,
        "row_estimator": ROW_ESTIMATOR_ID,
    }
    fingerprint = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()
    ).hexdigest()
    return payload | {"run_fingerprint": fingerprint}


def _init_worker(stack_path: str, wave: np.ndarray, run_fingerprint: str) -> None:
    global _RUN_FINGERPRINT
    _RUN_FINGERPRINT = run_fingerprint
    decompose.init_worker(
        wave,
        0.5,
        str(DEFAULT_DATA_ROOT),
        1.0e14,
        stack_path,
        fit_model=decompose.ADAM25K_NIV_CONTINUUM_FIT_MODEL,
        n_refinement_cycles=5,
        n_spline_knots=NIV_CONTINUUM_KWARGS["n_spline_knots"],
        n_zodi_spline_knots=NIV_CONTINUUM_KWARGS["n_zodi_spline_knots"],
        zodi_smooth_lambda=NIV_CONTINUUM_KWARGS["zodi_smooth_lambda"],
    )
    # Reuse the exact production input, mask, prior, and LSF path; only the OH
    # source class changes from Adam25k to PALACE Aijc VNF for PCA training.
    decompose._WORKER_DECOMPOSER = SkyDecompPalaceAijcVNFNivContinuumLSFSpline2D


def _cache_path(source_row: int) -> Path:
    return CACHE_DIR / f"row-{source_row:05d}.npz"


def _row_fingerprint(source_row: int) -> str:
    return hashlib.sha256(f"{_RUN_FINGERPRINT}:{source_row}".encode()).hexdigest()


def _fit_one(row: dict[str, object]) -> dict[str, object]:
    source_row = int(row["source_row"])
    output = _cache_path(source_row)
    fingerprint = _row_fingerprint(source_row)
    if output.is_file():
        try:
            with np.load(output, allow_pickle=False) as cached:
                if (
                    str(cached["fingerprint"].item()) == fingerprint
                    and np.isfinite(cached["amplitude"]).all()
                ):
                    return json.loads(str(cached["summary_json"].item())) | {
                        "from_cache": True
                    }
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass

    started = time.perf_counter()
    flux = (
        np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=np.float64)
        * decompose._WORKER_FACTOR
    )
    ivar = np.ones_like(flux)
    if decompose._WORKER_SCIENCE_LINE_MASK is not None:
        ivar[decompose._science_line_mask_for_row(source_row)] = 0.0
    model = decompose._telluric_decomposer("sky2", source_row)
    decompose._install_split_zodi_amplitude_prior(model, "sky2", source_row)
    result = model.fit(flux, ivar, verbose=False)
    if result.fit_status not in {"Solved", "AlmostSolved"}:
        raise RuntimeError(result.fit_summary)
    residual = flux - np.asarray(result.bestfit_lsf, dtype=np.float64)
    design = _individual_line_design(model)
    names = _individual_line_names(model)
    use = np.isfinite(residual) & np.isfinite(ivar) & (ivar > 0.0)
    fitted_design = design[use]
    column_norm = np.sqrt(np.asarray(fitted_design.power(2).sum(axis=0)).ravel())
    active = column_norm > 0.0
    normalized = fitted_design[:, active] @ sp.diags(1.0 / column_norm[active])
    solution = lsqr(
        normalized,
        residual[use],
        damp=np.sqrt(RIDGE_LAMBDA),
        atol=1.0e-6,
        btol=1.0e-6,
        iter_lim=2000,
    )
    amplitude = np.zeros(names.size, dtype=np.float64)
    amplitude[active] = solution[0] / column_norm[active]
    summary = {
        "source_row": source_row,
        "expnum": int(row["expnum"]),
        "fit_status": result.fit_status,
        "baseline_rms": float(np.sqrt(np.mean(residual[use] ** 2))),
        "solver_iterations": int(solution[2]),
        "solver_condition": float(solution[6]),
        "elapsed_sec": time.perf_counter() - started,
        "from_cache": False,
    }
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        fingerprint=np.asarray(fingerprint),
        summary_json=np.asarray(json.dumps(summary, sort_keys=True)),
        line_names=names,
        line_wave=np.asarray(model._line_wave, dtype=np.float64),
        line_group=np.asarray(model._line_group, dtype=np.int64),
        active_line=active,
        amplitude=amplitude,
    )
    os.replace(temporary, output)
    return summary


def _fit_one_safe(row: dict[str, object]) -> dict[str, object]:
    try:
        return _fit_one(row)
    except Exception as error:
        return {
            "source_row": int(row["source_row"]),
            "expnum": int(row["expnum"]),
            "status": "error",
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
            "from_cache": False,
        }


def fit_corpus(
    stack_path: Path,
    wave: np.ndarray,
    selection: pd.DataFrame,
    workers: int,
    run_fingerprint: str,
) -> pd.DataFrame:
    started = time.perf_counter()
    results = []
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=mp.get_context("spawn"),
        initializer=_init_worker,
        initargs=(str(stack_path), wave, run_fingerprint),
    ) as executor:
        futures = [
            executor.submit(_fit_one_safe, row)
            for row in selection.to_dict("records")
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            results.append(future.result())
            if completed == 1 or completed % 25 == 0 or completed == len(futures):
                failures = sum(item.get("status") == "error" for item in results)
                print(
                    f"completed={completed}/{len(futures)} "
                    f"elapsed={(time.perf_counter() - started) / 60.0:.1f} min "
                    f"failures={failures}",
                    flush=True,
                )
    frame = pd.DataFrame(results).sort_values("source_row").reset_index(drop=True)
    frame.to_csv(OUTPUT_DIR / "fit_manifest.csv", index=False)
    return frame


def _canonicalize(components: np.ndarray) -> np.ndarray:
    components = components.copy()
    for index in range(components.shape[0]):
        pivot = int(np.argmax(np.abs(components[index])))
        if components[index, pivot] < 0.0:
            components[index] *= -1.0
    return components


def build_asset(
    stack_path: Path,
    selection: pd.DataFrame,
    fit_manifest: pd.DataFrame,
    provenance: dict[str, object],
) -> dict[str, object]:
    failures = (
        set(
            fit_manifest.loc[
                fit_manifest["status"] == "error", "source_row"
            ].astype(int)
        )
        if "status" in fit_manifest
        else set()
    )
    rows = []
    summaries = []
    active_rows = []
    names_ref = wave_ref = group_ref = None
    for row in selection.to_dict("records"):
        source_row = int(row["source_row"])
        if source_row in failures:
            continue
        with np.load(_cache_path(source_row), allow_pickle=False) as fit:
            names = np.asarray(fit["line_names"])
            line_wave = np.asarray(fit["line_wave"], dtype=np.float64)
            line_group = np.asarray(fit["line_group"], dtype=np.int64)
            active = np.asarray(fit["active_line"], dtype=bool)
            if names_ref is None:
                names_ref, wave_ref, group_ref = names, line_wave, line_group
            else:
                np.testing.assert_array_equal(names, names_ref)
                np.testing.assert_array_equal(line_wave, wave_ref)
                np.testing.assert_array_equal(line_group, group_ref)
            rows.append(np.asarray(fit["amplitude"], dtype=np.float64))
            active_rows.append(active)
            summaries.append(json.loads(str(fit["summary_json"].item())))
    amplitude = np.stack(rows)
    active_matrix = np.stack(active_rows)
    metrics = pd.DataFrame(summaries)
    rms = metrics["baseline_rms"].to_numpy(dtype=np.float64)
    center = float(np.median(rms))
    robust_sigma = float(1.4826 * np.median(np.abs(rms - center)))
    keep = np.ones(rms.size, dtype=bool)
    if robust_sigma > 0.0:
        keep = rms <= center + OUTLIER_SIGMA * robust_sigma
    outlier_rows = metrics.loc[~keep, "source_row"].astype(int).tolist()
    amplitude = amplitude[keep]
    active_count = active_matrix[keep].sum(axis=0, dtype=np.int64)
    active_global = active_count > 0
    metrics["keep_for_pca"] = keep
    metrics.to_csv(OUTPUT_DIR / "outlier_analysis.csv", index=False)
    if amplitude.shape[0] <= PCA_COMPONENTS or amplitude.shape[1] != 11_552:
        raise ValueError(f"Invalid PCA training matrix: {amplitude.shape}")
    if np.any(~np.isfinite(amplitude)):
        raise ValueError("The PCA training matrix contains non-finite values")

    mean = amplitude.mean(axis=0, dtype=np.float64)
    centered = amplitude - mean
    gram = centered @ centered.T
    eigenvalue, eigenvector = np.linalg.eigh(gram)
    order = np.argsort(eigenvalue)[::-1]
    eigenvalue = np.maximum(eigenvalue[order], 0.0)
    positive = eigenvalue > max(float(eigenvalue[0]) * 1.0e-14, 0.0)
    eigenvalue = eigenvalue[positive]
    eigenvector = eigenvector[:, order][:, positive]
    components = _canonicalize(
        (eigenvector.T @ centered) / np.sqrt(eigenvalue)[:, None]
    )
    stored = components[:PCA_COMPONENTS]
    np.testing.assert_allclose(
        stored @ stored.T, np.eye(PCA_COMPONENTS), atol=5.0e-11
    )
    explained_ratio = eigenvalue / np.sum(centered**2)
    with fits.open(stack_path, memmap=True, lazy_load_hdus=True) as hdul:
        native_wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()

    metadata = {
        **provenance,
        "basis_id": "palace-aijc-vnf-niv-continuum-line-amplitude-pca30-v1",
        "source_model": "SkyDecompPalaceAijcVNFNivContinuumLSFSpline2D",
        "source_residual": "FLUX_SKY_FAR * 1e14 - bestfit_lsf",
        "source_oh_strength": "PALACE Aijc * gi",
        "source_oh_group_keys": ["v_upper", "N_upper", "F_upper"],
        "input_spectra": int(len(selection)),
        "successful_input_fits": int(len(rows)),
        "training_spectra": int(amplitude.shape[0]),
        "failed_source_rows": sorted(failures),
        "rms_outlier_source_rows": outlier_rows,
        "rms_outlier_rule": "baseline RMS > median + 7 * 1.4826 * MAD",
        "selection_method": "1000 positions evenly spaced over source-row-sorted valid telluric metadata",
        "native_wave_pixels": int(native_wave.size),
        "wave_sha256": wave_sha256(native_wave),
        "line_amplitude_columns": int(amplitude.shape[1]),
        "observable_line_transitions": int(np.count_nonzero(active_global)),
        "zero_support_edge_transitions": int(np.count_nonzero(~active_global)),
        "variable_support_line_transitions": int(
            np.count_nonzero(
                (active_count > 0) & (active_count < amplitude.shape[0])
            )
        ),
        "inactive_amplitude_estimator": "zero minimum-norm ridge solution",
        "line_names_sha256": hashlib.sha256(
            "\n".join(names_ref.tolist()).encode()
        ).hexdigest(),
        "amplitude_solver": SOLVER_ID,
        "amplitude_weights": "production IVAR mask on all 12401 native pixels",
        "amplitude_constraints": "none; all 11552 individual amplitudes signed",
        "ridge_lambda": RIDGE_LAMBDA,
        "dtype": "float64",
        "preprocessing": "per-amplitude column mean subtraction only",
        "amplitude_normalization": "unit native-grid integral of each exact telluric-plus-LSF profile",
        "removed_or_replaced_wavelength_pixels": 0,
        "pca_algorithm": "exact eigendecomposition of the sample Gram matrix",
        "stored_components": PCA_COMPONENTS,
        "selected_components": PCA_COMPONENTS,
        "available_component_counts": [PCA_COMPONENTS],
        "cumulative_explained_variance_30": float(
            explained_ratio[:PCA_COMPONENTS].sum()
        ),
    }
    PCA_ASSET.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        PCA_ASSET,
        wave=native_wave,
        line_names=names_ref,
        line_wave=wave_ref,
        line_group=group_ref,
        active_line=active_global,
        active_line_training_rows=active_count,
        amplitude_mean=mean,
        components=stored,
        explained_variance=eigenvalue[:PCA_COMPONENTS] / (amplitude.shape[0] - 1),
        explained_variance_ratio=explained_ratio[:PCA_COMPONENTS],
        all_explained_variance_ratio=explained_ratio,
        parameter_scale=np.std(amplitude, axis=0, ddof=1),
        training_source_row=metrics.loc[keep, "source_row"].to_numpy(dtype=np.int64),
        training_expnum=metrics.loc[keep, "expnum"].to_numpy(dtype=np.int64),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    summary = metadata | {
        "asset": str(PCA_ASSET),
        "asset_sha256": file_sha256(PCA_ASSET),
    }
    (OUTPUT_DIR / "pca_asset_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stack", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--asset-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or (args.limit is not None and args.limit < 1):
        raise ValueError("--workers and --limit must be positive")

    selection = build_selection(args.stack)
    input_sha256 = file_sha256(args.stack)
    provenance = _run_provenance(args.stack, input_sha256)
    (OUTPUT_DIR / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    selected = selection if args.limit is None else selection.iloc[: args.limit]
    with fits.open(args.stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    if not args.asset_only:
        manifest = fit_corpus(
            args.stack,
            wave,
            selected,
            args.workers,
            provenance["run_fingerprint"],
        )
    else:
        manifest = pd.read_csv(OUTPUT_DIR / "fit_manifest.csv")
    if args.limit is not None:
        return
    print(json.dumps(build_asset(args.stack, selection, manifest, provenance), indent=2))


if __name__ == "__main__":
    main()
