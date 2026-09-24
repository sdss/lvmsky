"""Train PCA30 from quality-filtered Far-Sky decomposition residuals.

The expensive physical fit is not repeated. This program reads the observed
Far-Sky flux and continuous per-spectrum LSF from a completed ``sky2``
decomposition, recomputes each residual against ``BESTFIT_LSF``, projects it
onto the exact 11,552-line design, then trains PCA30 in signed line-amplitude
space.
"""

# ruff: noqa: E402 -- native thread limits must be set before NumPy imports.

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

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
from skysub.sky_decomp.fit import LSF_CHANNELS
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT, file_sha256, wave_sha256
from skysub.sky_decomp.residual_pca import _individual_line_design, _individual_line_names
from skysub.sky_decomp.result_io import load_lsf_surface_state


PCA_COMPONENTS = 30
RIDGE_LAMBDA = 1.0e-4
QC_SIGMA = 7.0
FINITE_FRACTION_MIN = 0.995
DISPLAY_BIN_PIXELS = 8
MIN_LINE_SUPPORT = 0.5
SOLVER_ID = "scipy.sparse.linalg.lsqr-unit-integral-line-ridge-v5-support"
PCA_ID = "deterministic-randomized-svd-column-centred-v1"

OUTPUT_DIR: Path | None = None
CACHE_DIR: Path | None = None
DECOMP_PATH: Path | None = None
_DECOMP_HDU = None
_RUN_FINGERPRINT: str | None = None
_CATALOG_SHA256: str | None = None


def _text_array(values) -> np.ndarray:
    return np.char.strip(np.asarray(values).astype(str))


def _robust_location_scale(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan, np.nan
    centre = float(np.median(finite))
    scale = float(1.4826 * np.median(np.abs(finite - centre)))
    if not np.isfinite(scale) or scale == 0.0:
        scale = float(np.std(finite))
    return centre, scale


def _robust_keep(
    values: np.ndarray, sigma: float, *, log_positive: bool = False
) -> tuple[np.ndarray, dict[str, object]]:
    values = np.asarray(values, dtype=np.float64)
    transformed = np.full(values.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(values) & ((values > 0.0) if log_positive else True)
    transformed[valid] = np.log10(values[valid]) if log_positive else values[valid]
    centre, scale = _robust_location_scale(transformed)
    keep = np.isfinite(transformed)
    if np.isfinite(scale) and scale > 0.0:
        keep &= np.abs(transformed - centre) <= float(sigma) * scale
    convert = (lambda value: 10.0**value) if log_positive else (lambda value: value)
    return keep, {
        "centre": convert(centre),
        "robust_sigma": scale,
        "lower": convert(centre - float(sigma) * scale),
        "upper": convert(centre + float(sigma) * scale),
        "transform": "log10" if log_positive else "identity",
    }


def _binned_mean(values: np.ndarray, width: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return np.asarray(
        [
            np.nanmean(values[start : min(start + width, values.size)])
            for start in range(0, values.size, width)
        ],
        dtype=np.float32,
    )


def _metric_windows(wave: np.ndarray) -> dict[str, np.ndarray]:
    windows = {
        channel: (wave >= (-np.inf if lower is None else lower))
        & (wave < (np.inf if upper is None else upper))
        for channel, lower, upper in LSF_CHANNELS
    }
    for boundary in (5787.0, 7454.0):
        key = str(int(boundary))
        windows[f"join_{key}_left"] = (wave >= boundary - 25.0) & (wave < boundary - 5.0)
        windows[f"join_{key}_right"] = (wave > boundary + 5.0) & (wave <= boundary + 25.0)
    return windows


def build_quality_selection(
    stack_path: Path,
    decomp_path: Path,
    output_dir: Path,
    qc_sigma: float = QC_SIGMA,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, object]]:
    """Measure reduction-artifact metrics on the original Far-Sky spectra."""
    with fits.open(stack_path, memmap=True, lazy_load_hdus=True) as source, fits.open(
        decomp_path, memmap=True, lazy_load_hdus=True
    ) as fitted:
        wave = np.asarray(source["WAVE"].data, dtype=np.float64).copy()
        flux = source["FLUX_SKY_FAR"].data
        meta = source["META"].data
        fit_meta = fitted["META"].data
        if len(meta) != len(fit_meta) or fitted["BESTFIT_LSF"].data.shape != flux.shape:
            raise ValueError("The stack and sky2 decomposition are not row/grid aligned")
        n_rows = len(meta)
        windows = _metric_windows(wave)
        values = {
            name: np.full(n_rows, np.nan, dtype=np.float64)
            for name in ("median_B", "median_R", "median_Z", "join_BR", "join_RZ")
        }
        finite_fraction = np.zeros(n_rows, dtype=np.float64)
        for start in range(0, n_rows, 256):
            stop = min(start + 256, n_rows)
            block = np.asarray(flux[start:stop], dtype=np.float64)
            finite_fraction[start:stop] = np.isfinite(block).mean(axis=1)
            for channel in ("B", "R", "Z"):
                values[f"median_{channel}"][start:stop] = np.nanmedian(
                    block[:, windows[channel]], axis=1
                )
            for label, boundary in (("BR", "5787"), ("RZ", "7454")):
                left = np.nanmedian(block[:, windows[f"join_{boundary}_left"]], axis=1)
                right = np.nanmedian(block[:, windows[f"join_{boundary}_right"]], axis=1)
                denominator = np.abs(left) + np.abs(right)
                values[f"join_{label}"][start:stop] = np.divide(
                    2.0 * (right - left),
                    denominator,
                    out=np.full_like(left, np.nan),
                    where=denominator > 0.0,
                )

        status = _text_array(fit_meta["fit_status"])
        labels = np.char.lower(_text_array(meta["sky_far_label"]))
        lsf_meta = fitted["LSF_META"].data
        lsf_available_count = np.bincount(
            np.asarray(lsf_meta["spectrum_index"], dtype=int),
            weights=np.asarray(lsf_meta["available"], dtype=bool),
            minlength=n_rows,
        )
        base_valid = (
            np.isin(status, ("Solved", "AlmostSolved"))
            & np.isin(labels, ("skye", "skyw"))
            & (finite_fraction >= FINITE_FRACTION_MIN)
            & (lsf_available_count == 3)
        )
        frame = pd.DataFrame(
            {
                "source_row": np.arange(n_rows, dtype=np.int64),
                "expnum": np.asarray(meta["expnum"], dtype=np.int64),
                "mjd": np.asarray(meta["mjd"], dtype=np.int64),
                "sky_far_label": labels,
                "fit_status": status,
                "finite_fraction": finite_fraction,
                **values,
                "pass_base_contract": base_valid,
            }
        )

    thresholds: dict[str, dict[str, object]] = {}
    metric_passes = []
    for name in ("median_B", "median_R", "median_Z", "join_BR", "join_RZ"):
        keep = np.zeros(len(frame), dtype=bool)
        keep_valid, thresholds[name] = _robust_keep(
            frame.loc[base_valid, name].to_numpy(dtype=np.float64),
            qc_sigma,
            log_positive=name.startswith("median_"),
        )
        keep[np.flatnonzero(base_valid)] = keep_valid
        frame[f"pass_{name}"] = keep
        metric_passes.append(keep)
    frame["keep_for_pca"] = base_valid & np.logical_and.reduce(metric_passes)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "far_sky_quality_selection.csv", index=False)
    summary = {
        "input_rows": int(len(frame)),
        "base_contract_rows": int(base_valid.sum()),
        "retained_rows": int(frame["keep_for_pca"].sum()),
        "rejected_rows": int((~frame["keep_for_pca"]).sum()),
        "qc_sigma": float(qc_sigma),
        "finite_fraction_min": FINITE_FRACTION_MIN,
        "metrics": thresholds,
        "join_definition": "2*(right-left)/(abs(left)+abs(right)); 20 A median windows separated by 10 A",
        "channel_boundaries_angstrom": [5787.0, 7454.0],
    }
    (output_dir / "far_sky_quality_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if summary["retained_rows"] <= PCA_COMPONENTS:
        raise ValueError(f"Only {summary['retained_rows']} spectra survive Far-Sky QC")
    return frame.loc[frame["keep_for_pca"]].reset_index(drop=True), wave, summary


def _catalog_hash(
    line_names: np.ndarray, line_wave: np.ndarray, line_group: np.ndarray
) -> str:
    digest = hashlib.sha256()
    digest.update("\n".join(np.asarray(line_names).astype(str)).encode())
    digest.update(np.asarray(line_wave, dtype="<f8").tobytes())
    digest.update(np.asarray(line_group, dtype="<i8").tobytes())
    return digest.hexdigest()


def _init_worker(
    stack_path: str,
    decomp_path: str,
    wave: np.ndarray,
    output_dir: str,
    run_fingerprint: str,
    catalog_sha256: str | None,
    worker_counter=None,
    pin_cpu: bool = False,
) -> None:
    global OUTPUT_DIR, CACHE_DIR, DECOMP_PATH, _DECOMP_HDU, _RUN_FINGERPRINT, _CATALOG_SHA256
    OUTPUT_DIR = Path(output_dir)
    CACHE_DIR = OUTPUT_DIR / "line_amplitudes"
    DECOMP_PATH = Path(decomp_path)
    _RUN_FINGERPRINT = run_fingerprint
    _CATALOG_SHA256 = catalog_sha256
    decompose.init_worker(
        wave,
        0.5,
        str(DEFAULT_DATA_ROOT),
        1.0e14,
        stack_path,
        fit_model=decompose.PALACE_VNF_SPLIT_ZODI_FIT_MODEL,
        n_refinement_cycles=5,
        n_spline_knots=11,
        n_zodi_spline_knots=1,
        zodi_smooth_lambda=0.1,
        worker_counter=worker_counter,
        pin_cpu=pin_cpu,
    )
    _DECOMP_HDU = fits.open(DECOMP_PATH, memmap=True, lazy_load_hdus=True)


def _close_worker_files() -> None:
    global _DECOMP_HDU
    if _DECOMP_HDU is not None:
        _DECOMP_HDU.close()
        _DECOMP_HDU = None
    if decompose._WORKER_HDU is not None:
        decompose._WORKER_HDU.close()


def _cache_path(source_row: int) -> Path:
    assert CACHE_DIR is not None
    return CACHE_DIR / f"row-{source_row:05d}.npz"


def _row_fingerprint(source_row: int) -> str:
    return hashlib.sha256(f"{_RUN_FINGERPRINT}:{source_row}".encode()).hexdigest()


def _model_and_design(source_row: int):
    assert DECOMP_PATH is not None
    model = decompose._telluric_decomposer("sky2", source_row)
    model._set_lsf_state(load_lsf_surface_state(DECOMP_PATH, source_row))
    return model, _individual_line_design(model)


def _observable_line_columns(
    design: sp.csr_matrix, retained_design: sp.csr_matrix, min_support: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    full_norm = np.sqrt(np.asarray(design.power(2).sum(axis=0)).ravel())
    retained_norm = np.sqrt(
        np.asarray(retained_design.power(2).sum(axis=0)).ravel()
    )
    support = np.divide(
        retained_norm,
        full_norm,
        out=np.zeros_like(retained_norm),
        where=full_norm > 0.0,
    )
    return retained_norm, support, (full_norm > 0.0) & (support >= min_support)


def _fit_one(
    row: dict[str, object],
    ridge_lambda: float,
    display_bin: int,
    min_line_support: float,
) -> dict[str, object]:
    source_row = int(row["source_row"])
    output = _cache_path(source_row)
    fingerprint = _row_fingerprint(source_row)
    if output.is_file():
        try:
            with np.load(output, allow_pickle=False) as cached:
                if str(cached["fingerprint"].item()) == fingerprint:
                    return json.loads(str(cached["summary_json"].item())) | {"from_cache": True}
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass

    started = time.perf_counter()
    observed = (
        np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=np.float64)
        * decompose._WORKER_FACTOR
    )
    residual = observed - np.asarray(
        _DECOMP_HDU["BESTFIT_LSF"].data[source_row], dtype=np.float64
    )
    model, design = _model_and_design(source_row)
    names = _individual_line_names(model)
    catalog_sha256 = _catalog_hash(names, model._line_wave, model._line_group)
    if _CATALOG_SHA256 is not None and catalog_sha256 != _CATALOG_SHA256:
        raise ValueError("The individual-line catalog changed between workers")
    use = np.isfinite(residual)
    if decompose._WORKER_SCIENCE_LINE_MASK is not None:
        use &= ~decompose._science_line_mask_for_row(source_row)
    fitted_design = design[use]
    column_norm, _, active = _observable_line_columns(
        design, fitted_design, min_line_support
    )
    normalized = fitted_design[:, active] @ sp.diags(1.0 / column_norm[active])
    solution = lsqr(
        normalized,
        residual[use],
        damp=np.sqrt(ridge_lambda),
        atol=1.0e-6,
        btol=1.0e-6,
        iter_lim=2000,
    )
    amplitude = np.zeros(names.size, dtype=np.float64)
    amplitude[active] = solution[0] / column_norm[active]
    line_fit = np.asarray(design @ amplitude).ravel()
    post_residual = residual - line_fit
    summary = {
        "source_row": source_row,
        "expnum": int(row["expnum"]),
        "sky_far_label": str(row["sky_far_label"]),
        "baseline_rms": float(np.sqrt(np.mean(residual[use] ** 2))),
        "post_line_rms": float(np.sqrt(np.mean(post_residual[use] ** 2))),
        "solver_istop": int(solution[1]),
        "solver_iterations": int(solution[2]),
        "solver_condition": float(solution[6]),
        "active_lines": int(active.sum()),
        "elapsed_sec": time.perf_counter() - started,
        "from_cache": False,
    }
    assert CACHE_DIR is not None
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        fingerprint=np.asarray(fingerprint),
        catalog_sha256=np.asarray(catalog_sha256),
        summary_json=np.asarray(json.dumps(summary, sort_keys=True)),
        active_line=active,
        amplitude=amplitude,
        residual_display=_binned_mean(residual, display_bin),
        post_residual_display=_binned_mean(post_residual, display_bin),
    )
    os.replace(temporary, output)
    return summary


def _fit_one_safe(
    row: dict[str, object],
    ridge_lambda: float,
    display_bin: int,
    min_line_support: float,
):
    try:
        return _fit_one(row, ridge_lambda, display_bin, min_line_support)
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
    decomp_path: Path,
    wave: np.ndarray,
    selection: pd.DataFrame,
    output_dir: Path,
    workers: int,
    run_fingerprint: str,
    catalog_sha256: str,
    ridge_lambda: float,
    display_bin: int,
    min_line_support: float,
) -> pd.DataFrame:
    started = time.perf_counter()
    results = []
    context = mp.get_context("spawn")
    worker_counter = context.Value("i", 0)
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=_init_worker,
        initargs=(
            str(stack_path), str(decomp_path), wave, str(output_dir), run_fingerprint,
            catalog_sha256, worker_counter, True,
        ),
    ) as executor:
        futures = [
            executor.submit(
                _fit_one_safe, row, ridge_lambda, display_bin, min_line_support
            )
            for row in selection.to_dict("records")
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            results.append(future.result())
            if completed == 1 or completed % 100 == 0 or completed == len(futures):
                failures = sum(item.get("status") == "error" for item in results)
                print(
                    f"completed={completed}/{len(futures)} "
                    f"elapsed={(time.perf_counter() - started) / 60.0:.1f} min "
                    f"failures={failures}", flush=True,
                )
    frame = pd.DataFrame(results).sort_values("source_row").reset_index(drop=True)
    frame.to_csv(output_dir / "line_amplitude_fit_manifest.csv", index=False)
    return frame


def _canonicalize(components: np.ndarray) -> np.ndarray:
    components = components.copy()
    for index in range(components.shape[0]):
        pivot = int(np.argmax(np.abs(components[index])))
        if components[index, pivot] < 0.0:
            components[index] *= -1.0
    return components


def randomized_pca(
    matrix: np.ndarray,
    components: int,
    *,
    seed: int = 20260918,
    oversample: int = 15,
    power_iterations: int = 2,
    block_rows: int = 256,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Column-centred randomized SVD without a second full matrix copy."""
    n_rows, n_columns = matrix.shape
    if not 0 < components < min(n_rows, n_columns):
        raise ValueError("PCA component count must be below both matrix dimensions")
    mean = np.zeros(n_columns, dtype=np.float64)
    for start in range(0, n_rows, block_rows):
        mean += np.asarray(matrix[start : start + block_rows], dtype=np.float64).sum(axis=0)
    mean /= n_rows

    rng = np.random.default_rng(seed)
    width = min(components + oversample, min(n_rows, n_columns))
    omega = rng.standard_normal((n_columns, width))

    def right_multiply(right: np.ndarray) -> np.ndarray:
        output = np.empty((n_rows, right.shape[1]), dtype=np.float64)
        for start in range(0, n_rows, block_rows):
            stop = min(start + block_rows, n_rows)
            output[start:stop] = (
                np.asarray(matrix[start:stop], dtype=np.float64) - mean
            ) @ right
        return output

    def transpose_multiply(right: np.ndarray) -> np.ndarray:
        output = np.zeros((n_columns, right.shape[1]), dtype=np.float64)
        for start in range(0, n_rows, block_rows):
            stop = min(start + block_rows, n_rows)
            output += (
                np.asarray(matrix[start:stop], dtype=np.float64) - mean
            ).T @ right[start:stop]
        return output

    sample = right_multiply(omega)
    for _ in range(power_iterations):
        sample, _ = np.linalg.qr(sample, mode="reduced")
        feature_sample, _ = np.linalg.qr(transpose_multiply(sample), mode="reduced")
        sample = right_multiply(feature_sample)
    basis, _ = np.linalg.qr(sample, mode="reduced")
    projected = np.zeros((basis.shape[1], n_columns), dtype=np.float64)
    feature_ss = np.zeros(n_columns, dtype=np.float64)
    total_ss = 0.0
    for start in range(0, n_rows, block_rows):
        stop = min(start + block_rows, n_rows)
        centred = np.asarray(matrix[start:stop], dtype=np.float64) - mean
        projected += basis[start:stop].T @ centred
        feature_ss += np.sum(centred**2, axis=0)
        total_ss += float(np.sum(centred**2))
    _, singular_value, component = np.linalg.svd(projected, full_matrices=False)
    component = _canonicalize(component[:components])
    singular_value = singular_value[:components]
    explained_variance = singular_value**2 / (n_rows - 1)
    explained_ratio = singular_value**2 / total_ss
    parameter_scale = np.sqrt(feature_ss / (n_rows - 1))
    np.testing.assert_allclose(component @ component.T, np.eye(components), atol=2.0e-10)
    return mean, component, explained_variance, explained_ratio, parameter_scale


def _matrix_scores(
    matrix: np.ndarray, mean: np.ndarray, components: np.ndarray, block_rows: int = 256
) -> np.ndarray:
    scores = np.empty((matrix.shape[0], components.shape[0]), dtype=np.float32)
    for start in range(0, matrix.shape[0], block_rows):
        stop = min(start + block_rows, matrix.shape[0])
        scores[start:stop] = (
            (np.asarray(matrix[start:stop], dtype=np.float64) - mean) @ components.T
        ).astype(np.float32)
    return scores


def build_asset_and_analysis(
    selection: pd.DataFrame,
    fit_manifest: pd.DataFrame,
    wave: np.ndarray,
    output_dir: Path,
    provenance: dict[str, object],
    catalog: dict[str, np.ndarray],
    display_bin: int,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    failures = (
        set(fit_manifest.loc[fit_manifest["status"] == "error", "source_row"].astype(int))
        if "status" in fit_manifest else set()
    )
    retained = selection.loc[~selection["source_row"].isin(failures)].reset_index(drop=True)
    n_rows = len(retained)
    n_lines = len(catalog["line_names"])
    if n_rows <= PCA_COMPONENTS or n_lines != 11_552:
        raise ValueError(f"Invalid PCA training matrix shape: {(n_rows, n_lines)}")
    amplitude = np.lib.format.open_memmap(
        output_dir / "line_amplitude_matrix.npy", mode="w+", dtype=np.float64,
        shape=(n_rows, n_lines),
    )
    display_size = (wave.size + display_bin - 1) // display_bin
    residual_display = np.empty((n_rows, display_size), dtype=np.float32)
    post_display = np.empty_like(residual_display)
    active_count = np.zeros(n_lines, dtype=np.int64)
    summaries = []
    for index, source_row in enumerate(retained["source_row"].astype(int)):
        with np.load(
            output_dir / "line_amplitudes" / f"row-{source_row:05d}.npz",
            allow_pickle=False,
        ) as fit:
            if str(fit["catalog_sha256"].item()) != provenance["line_catalog_sha256"]:
                raise ValueError(f"Line catalog mismatch at source row {source_row}")
            amplitude[index] = np.asarray(fit["amplitude"], dtype=np.float64)
            residual_display[index] = np.asarray(fit["residual_display"], dtype=np.float32)
            post_display[index] = np.asarray(fit["post_residual_display"], dtype=np.float32)
            active_count += np.asarray(fit["active_line"], dtype=bool)
            summaries.append(json.loads(str(fit["summary_json"].item())))
    globally_active = active_count == n_rows
    if np.count_nonzero(globally_active) <= PCA_COMPONENTS:
        raise ValueError("Too few globally observable line-amplitude columns")
    amplitude[:, ~globally_active] = 0.0
    amplitude.flush()
    metrics = pd.DataFrame(summaries)
    metrics.to_csv(output_dir / "line_amplitude_metrics.csv", index=False)
    mean, components, variance, ratio, parameter_scale = randomized_pca(
        amplitude, PCA_COMPONENTS
    )
    scores = _matrix_scores(amplitude, mean, components)
    order = np.argsort(scores[:, 0], kind="stable")
    display_wave = np.asarray(
        [
            np.mean(wave[start : min(start + display_bin, wave.size)])
            for start in range(0, wave.size, display_bin)
        ], dtype=np.float64,
    )
    analysis_path = output_dir / "far_sky_pca30_analysis.npz"
    np.savez_compressed(
        analysis_path,
        wave_display=display_wave,
        source_row=retained["source_row"].to_numpy(dtype=np.int64),
        sky_far_label=retained["sky_far_label"].to_numpy(dtype="U8"),
        scores=scores,
        sort_order=order,
        residual_display=residual_display,
        post_residual_display=post_display,
        baseline_rms=metrics["baseline_rms"].to_numpy(dtype=np.float64),
        post_line_rms=metrics["post_line_rms"].to_numpy(dtype=np.float64),
    )
    metadata = {
        **provenance,
        "basis_id": "far-sky-palace-aijc-vnf-split-zodi-line-amplitude-pca30-v1",
        "source_residual": "recomputed FLUX_SKY_FAR * 1e14 - BESTFIT_LSF (the stored legacy RESID is not used)",
        "quality_filter": "original full-spectrum B/R/Z medians and B-R/R-Z joins before residual analysis",
        "training_spectra": n_rows,
        "native_wave_pixels": int(wave.size),
        "wave_sha256": wave_sha256(wave),
        "line_amplitude_columns": n_lines,
        "observable_line_transitions": int(np.count_nonzero(globally_active)),
        "line_observability": (
            "retained/full profile L2 norm >= min_line_support in every training row"
        ),
        "min_line_support": float(provenance["min_line_support"]),
        "amplitude_solver": SOLVER_ID,
        "amplitude_constraints": "none; all individual-line residual amplitudes are signed",
        "ridge_lambda": float(provenance["ridge_lambda"]),
        "preprocessing": "per-line amplitude column mean subtraction only",
        "pca_algorithm": PCA_ID,
        "pca_random_seed": 20260918,
        "pca_oversample": 15,
        "pca_power_iterations": 2,
        "stored_components": PCA_COMPONENTS,
        "cumulative_explained_variance_30": float(ratio.sum()),
    }
    asset_path = output_dir / "far_sky_line_amplitude_pca30.npz"
    np.savez_compressed(
        asset_path,
        wave=wave,
        line_names=catalog["line_names"],
        line_wave=catalog["line_wave"],
        line_group=catalog["line_group"],
        active_line=globally_active,
        active_line_training_rows=active_count,
        amplitude_mean=mean,
        components=components,
        explained_variance=variance,
        explained_variance_ratio=ratio,
        parameter_scale=parameter_scale,
        training_source_row=retained["source_row"].to_numpy(dtype=np.int64),
        training_expnum=retained["expnum"].to_numpy(dtype=np.int64),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    summary = metadata | {
        "asset": str(asset_path),
        "asset_sha256": file_sha256(asset_path),
        "analysis": str(analysis_path),
        "successful_line_fits": n_rows,
        "failed_line_fits": len(failures),
        "median_baseline_rms": float(np.median(metrics["baseline_rms"])),
        "median_post_line_rms": float(np.median(metrics["post_line_rms"])),
    }
    (output_dir / "pca30_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary, retained["source_row"].to_numpy(dtype=np.int64), scores


def build_example_spectra(
    stack_path: Path,
    decomp_path: Path,
    wave: np.ndarray,
    output_dir: Path,
    run_fingerprint: str,
    catalog_sha256: str,
    source_rows: np.ndarray,
    scores: np.ndarray,
) -> None:
    order = np.argsort(scores[:, 0])
    chosen = order[np.rint(np.linspace(0.05, 0.95, 9) * (len(order) - 1)).astype(int)]
    _init_worker(
        str(stack_path), str(decomp_path), wave, str(output_dir), run_fingerprint,
        catalog_sha256,
    )
    residuals, models, posts = [], [], []
    try:
        for index in chosen:
            source_row = int(source_rows[index])
            observed = (
                np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=np.float64)
                * decompose._WORKER_FACTOR
            )
            residual = observed - np.asarray(
                _DECOMP_HDU["BESTFIT_LSF"].data[source_row], dtype=np.float64
            )
            with np.load(_cache_path(source_row), allow_pickle=False) as fit:
                amplitude = np.asarray(fit["amplitude"], dtype=np.float64)
            _, design = _model_and_design(source_row)
            model = np.asarray(design @ amplitude).ravel()
            residuals.append(residual)
            models.append(model)
            posts.append(residual - model)
    finally:
        _close_worker_files()
    np.savez_compressed(
        output_dir / "typical_residual_examples.npz",
        wave=wave,
        source_row=source_rows[chosen],
        pc1_score=scores[chosen, 0],
        residual=np.asarray(residuals),
        line_fit=np.asarray(models),
        post_residual=np.asarray(posts),
    )


def _provenance(
    stack_path: Path,
    decomp_path: Path,
    selection_path: Path,
    catalog_sha256: str,
    ridge_lambda: float,
    qc_sigma: float,
    display_bin: int,
    min_line_support: float,
) -> dict[str, object]:
    source_paths = [Path(__file__).resolve(), REPO_ROOT / "skysub/decompose_parallel.py"]
    source_paths.extend(sorted((REPO_ROOT / "skysub/sky_decomp").glob("*.py")))
    source_sha256 = {
        str(path.relative_to(REPO_ROOT)): file_sha256(path) for path in source_paths
    }
    source_rel = list(source_sha256)
    try:
        repository_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        relevant_status = subprocess.run(
            ["git", "status", "--porcelain=v1", "--", *source_rel],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        repository_commit = "unavailable"
        relevant_status = "unavailable"
    science_line_mask = {
        "enabled": bool(decompose.SCIENCE_LINE_MASK_ENABLED),
        "lines_air_angstrom": [
            [str(name), float(wavelength)]
            for name, wavelength in decompose.SCIENCE_EMISSION_LINES
        ],
        "fwhm_multiple": float(decompose.SCIENCE_LINE_MASK_FWHM_MULTIPLE),
        "velocity_km_s": float(decompose.SCIENCE_LINE_MASK_VELOCITY_KM_S),
        "min_half_width_angstrom": float(
            decompose.SCIENCE_LINE_MASK_MIN_HALF_WIDTH_A
        ),
        "centre_on_halpha": bool(decompose.SCIENCE_LINE_MASK_CENTRE_ON_HALPHA),
        "max_shift_km_s": float(decompose.SCIENCE_LINE_MASK_MAX_SHIFT_KM_S),
        "centre_min_snr": float(decompose.SCIENCE_LINE_MASK_CENTRE_MIN_SNR),
    }
    payload = {
        "schema_version": 2,
        "stack_path": str(stack_path.resolve()),
        "stack_sha256": file_sha256(stack_path),
        "decomposition_path": str(decomp_path.resolve()),
        "decomposition_sha256": file_sha256(decomp_path),
        "selection_sha256": file_sha256(selection_path),
        "source_code_sha256": file_sha256(Path(__file__)),
        "source_sha256": source_sha256,
        "repository_commit": repository_commit,
        "repository_relevant_dirty": bool(relevant_status),
        "repository_relevant_status_sha256": hashlib.sha256(
            relevant_status.encode()
        ).hexdigest(),
        "bundle_manifest_sha256": file_sha256(
            DEFAULT_DATA_ROOT / "bundle_manifest.json"
        ),
        "science_line_mask": science_line_mask,
        "line_catalog_sha256": catalog_sha256,
        "far_sky_role": "sky2 / FLUX_SKY_FAR",
        "ridge_lambda": ridge_lambda,
        "qc_sigma": qc_sigma,
        "display_bin_pixels": display_bin,
        "min_line_support": min_line_support,
    }
    return payload | {
        "run_fingerprint": hashlib.sha256(
            json.dumps(payload, sort_keys=True).encode()
        ).hexdigest()
    }


def self_test() -> None:
    values = np.array([0.0, 0.1, -0.1, 0.05, 50.0])
    keep, _ = _robust_keep(values, 7.0)
    assert keep.tolist() == [True, True, True, True, False]
    np.testing.assert_allclose(_binned_mean(np.arange(5.0), 2), [0.5, 2.5, 4.0])
    design = sp.csr_matrix([[1.0, 0.0], [1.0, 1.0]])
    norm, support, active = _observable_line_columns(
        design, design[[True, False]], 0.5
    )
    np.testing.assert_allclose(norm, [1.0, 0.0])
    np.testing.assert_allclose(support, [1.0 / np.sqrt(2.0), 0.0])
    assert active.tolist() == [True, False]
    rng = np.random.default_rng(7)
    left = rng.normal(size=(120, 4))
    right, _ = np.linalg.qr(rng.normal(size=(30, 4)))
    matrix = left @ right.T + 1.0e-5 * rng.normal(size=(120, 30))
    _, components, _, ratio, _ = randomized_pca(
        matrix, 4, seed=3, oversample=4, power_iterations=2, block_rows=31
    )
    np.testing.assert_allclose(components @ components.T, np.eye(4), atol=1.0e-10)
    assert ratio.sum() > 0.999
    print("self-test: ok")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stack", nargs="?", type=Path)
    parser.add_argument("decomposition", nargs="?", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, default=64)
    parser.add_argument("--limit", type=int, help="Development-only cap after quality filtering")
    parser.add_argument("--ridge-lambda", type=float, default=RIDGE_LAMBDA)
    parser.add_argument("--qc-sigma", type=float, default=QC_SIGMA)
    parser.add_argument("--display-bin", type=int, default=DISPLAY_BIN_PIXELS)
    parser.add_argument("--min-line-support", type=float, default=MIN_LINE_SUPPORT)
    parser.add_argument("--asset-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.stack is None or args.decomposition is None:
        parser.error("stack and decomposition are required unless --self-test is used")
    if args.workers < 1 or args.display_bin < 1 or args.ridge_lambda < 0.0:
        raise ValueError("workers/display-bin must be positive and ridge-lambda non-negative")
    if not 0.0 < args.min_line_support <= 1.0:
        raise ValueError("min-line-support must be in (0, 1]")
    if args.limit is not None and args.limit <= PCA_COMPONENTS:
        raise ValueError(f"--limit must exceed {PCA_COMPONENTS}")
    output_dir = args.output_dir or args.decomposition.parent / "far_sky_pca30"
    output_dir.mkdir(parents=True, exist_ok=True)
    selection, wave, quality_summary = build_quality_selection(
        args.stack, args.decomposition, output_dir, args.qc_sigma
    )
    if args.limit is not None:
        selection = selection.iloc[: args.limit].copy()

    provisional = "catalog-bootstrap"
    _init_worker(
        str(args.stack), str(args.decomposition), wave, str(output_dir), provisional, None
    )
    try:
        first_row = int(selection.iloc[0]["source_row"])
        model, _ = _model_and_design(first_row)
        catalog = {
            "line_names": _individual_line_names(model),
            "line_wave": np.asarray(model._line_wave, dtype=np.float64),
            "line_group": np.asarray(model._line_group, dtype=np.int64),
        }
    finally:
        _close_worker_files()
    catalog_sha256 = _catalog_hash(**catalog)
    np.savez_compressed(output_dir / "line_catalog.npz", **catalog)
    provenance = _provenance(
        args.stack, args.decomposition,
        output_dir / "far_sky_quality_selection.csv",
        catalog_sha256, args.ridge_lambda, args.qc_sigma, args.display_bin,
        args.min_line_support,
    )
    provenance["quality_summary"] = quality_summary
    (output_dir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    if args.asset_only:
        manifest = pd.read_csv(output_dir / "line_amplitude_fit_manifest.csv")
    else:
        manifest = fit_corpus(
            args.stack, args.decomposition, wave, selection, output_dir, args.workers,
            str(provenance["run_fingerprint"]), catalog_sha256, args.ridge_lambda,
            args.display_bin, args.min_line_support,
        )
    summary, source_rows, scores = build_asset_and_analysis(
        selection, manifest, wave, output_dir, provenance, catalog, args.display_bin
    )
    build_example_spectra(
        args.stack, args.decomposition, wave, output_dir,
        str(provenance["run_fingerprint"]), catalog_sha256, source_rows, scores,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
