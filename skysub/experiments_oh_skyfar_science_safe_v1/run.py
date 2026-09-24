"""Skyfar-only, science-line-masked linear calibration of PALACE OH Aijc ratios."""

# ruff: noqa: E402 -- cap native thread pools before importing NumPy.

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import sys
import time
import traceback

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

for _name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
import pandas as pd
import scipy.sparse as sp
from astropy.io import fits

import skysub.decompose_parallel as decompose
from skysub.experiments_niv_integrated_decomposition_v1 import train_pca30 as base
from skysub.experiments_oh_branch_corrections_v1 import (
    fit_global_oh_branch_corrections as global_fit,
)
from skysub.sky_decomp.lsf_spline2d import (
    evaluate_lsf_density,
    evaluate_lsf_diagnostics,
)
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_OH_SUFFIX,
    SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX,
    file_sha256,
    wave_sha256,
)


SCHEMA = "oh-skyfar-fwhm-v3"
FACTOR = 1.0e14
R_LOWER_A = 5787.0
DISPLAY_BIN_PIXELS = 4
RIDGE_GRID = (1.0e-2, 1.0e-1, 1.0, 10.0, 100.0)
INFORMATION_QUANTILE = 0.8
MAX_CORRECTION_FACTOR = 4.0
LSF_REFERENCE_SAMPLE_ROWS = 64
PALACE_OUTPUT_SUFFIX = SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX

# Air wavelengths explicitly requested for protection.
SCIENCE_LINES = (
    ("[S III] 6312", 6312.06),
    ("[N II] 6548", 6548.05),
    ("H alpha", 6562.80),
    ("[N II] 6583", 6583.45),
    ("[S II] 6716", 6716.44),
    ("[S II] 6731", 6730.82),
    ("[S III] 9069", 9068.60),
    ("[S III] 9531", 9530.60),
)

_WAVE: np.ndarray | None = None
_RZ_PIXEL: np.ndarray | None = None
_SELECTED_LINE: np.ndarray | None = None
_SELECTED_Q: np.ndarray | None = None
_SELECTED_PARENT: np.ndarray | None = None
_LOADING: sp.csr_matrix | None = None
_OH_COEF_NAMES: tuple[str, ...] = ()


def _atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _file_identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def deterministic_mjd_split(
    mjd: np.ndarray,
    seed: int = 17,
    fractions: tuple[float, float, float] = (0.70, 0.15, 0.15),
) -> np.ndarray:
    """Return a deterministic, MJD-disjoint train/validation/holdout split."""
    values = np.asarray(mjd, dtype=np.int64)
    nights = np.unique(values)
    if nights.size < 3:
        raise ValueError("At least three MJD values are required")
    if not np.isclose(sum(fractions), 1.0):
        raise ValueError("Split fractions must sum to one")
    order = sorted(
        nights.tolist(),
        key=lambda value: hashlib.sha256(f"{seed}:{value}".encode()).digest(),
    )
    n_train = max(1, int(round(fractions[0] * len(order))))
    n_validation = max(1, int(round(fractions[1] * len(order))))
    if n_train + n_validation >= len(order):
        n_train, n_validation = len(order) - 2, 1
    lookup = {
        night: label
        for label, subset in (
            ("train", order[:n_train]),
            ("validation", order[n_train : n_train + n_validation]),
            ("holdout", order[n_train + n_validation :]),
        )
        for night in subset
    }
    result = np.asarray([lookup[int(value)] for value in values], dtype="U10")
    for label in ("train", "validation", "holdout"):
        if not np.any(result == label):
            raise AssertionError(f"Empty split: {label}")
    return result


def line_protection(
    wavelength: np.ndarray, reference_fwhm: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    protected = np.zeros(len(wavelength), dtype=bool)
    reason = np.full(len(wavelength), "", dtype=object)
    for (name, centre), fwhm in zip(SCIENCE_LINES, reference_fwhm, strict=True):
        hit = np.abs(wavelength - centre) <= float(fwhm)
        protected |= hit
        for index in np.flatnonzero(hit):
            reason[index] = "; ".join(filter(None, (str(reason[index]), name)))
    return protected, reason.astype(str)


def _individual_loading(
    catalog: pd.DataFrame,
) -> tuple[sp.csr_matrix, int]:
    """Full-rank q-centred contrasts, with protected lines fixed at zero."""
    eligible = ~catalog["protected_science_line"].to_numpy(bool)
    parent = catalog["parent_index"].to_numpy(int)
    q = catalog["q_aijc_gi"].to_numpy(float)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    parameter = 0
    for parent_id in np.unique(parent):
        line = np.flatnonzero((parent == parent_id) & eligible)
        if line.size < 2:
            continue
        reference = int(line[np.argmax(q[line])])
        for index in line:
            if index == reference:
                continue
            rows.extend((int(index), reference))
            columns.extend((parameter, parameter))
            values.extend((1.0, -float(q[index] / q[reference])))
            parameter += 1
    loading = sp.coo_matrix(
        (values, (rows, columns)),
        shape=(len(catalog), parameter),
    ).tocsr()
    return loading, int(eligible.sum())


def build_corpus(
    stack: Path,
    decomposition: Path,
    output: Path,
    qc_sigma: float,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, object]]:
    quality, wave, quality_summary = base.build_quality_selection(
        stack, decomposition, output, qc_sigma
    )
    row = quality["source_row"].to_numpy(dtype=int)
    with fits.open(stack, memmap=True, lazy_load_hdus=True) as hdul:
        meta = hdul["META"].data
        label = quality["sky_far_label"].astype(str).str.lower().to_numpy()
        pwv = np.asarray(meta["pwv_med"][row], dtype=float)
        fallback = np.asarray(meta["pwv_fallback"][row], dtype=bool)
        sun_alt = np.asarray(meta["sun_alt"][row], dtype=float)
        fraction = np.asarray(meta["fiberfrac_sky_far"][row], dtype=float)
        east = np.asarray(meta["skye_airmass"][row], dtype=float)
        west = np.asarray(meta["skyw_airmass"][row], dtype=float)
    far_airmass = np.where(label == "skye", east, west)
    strict = (
        np.isfinite(pwv)
        & (pwv > 0.0)
        & ~fallback
        & np.isfinite(sun_alt)
        & (sun_alt <= -18.0)
        & np.isfinite(far_airmass)
        & (far_airmass <= 1.8)
        & np.isfinite(fraction)
        & (fraction >= 0.8)
    )
    corpus = quality.loc[strict].copy().reset_index(drop=True)
    corpus["pwv_mm"] = pwv[strict]
    corpus["sun_alt_deg"] = sun_alt[strict]
    corpus["skyfar_airmass"] = far_airmass[strict]
    corpus["skyfar_fibre_fraction"] = fraction[strict]
    corpus["split"] = deterministic_mjd_split(corpus["mjd"].to_numpy(), seed=17)
    if corpus.empty or set(corpus["split"]) != {"train", "validation", "holdout"}:
        raise ValueError("Strict Skyfar corpus does not contain all three MJD splits")
    corpus.to_csv(output / "skyfar_corpus.csv", index=False)
    summary = {
        "upstream_quality": quality_summary,
        "strict_rows": int(len(corpus)),
        "strict_nights": int(corpus["mjd"].nunique()),
        "split_rows": corpus["split"].value_counts().sort_index().to_dict(),
        "split_nights": corpus.groupby("split")["mjd"].nunique().sort_index().to_dict(),
        "strict_contract": {
            "source": "FLUX_SKY_FAR only",
            "pwv": "finite positive measured PWV; fallback rows rejected",
            "sun_altitude_deg": "<= -18",
            "skyfar_airmass": "<= 1.8",
            "skyfar_fibre_fraction": ">= 0.8",
            "continuum_and_join_outliers": f"existing robust Far-Sky QC at {qc_sigma:g} sigma",
        },
    }
    _atomic_json(output / "corpus_summary.json", summary)
    return corpus, wave, summary


def build_catalog(
    stack: Path,
    decomposition: Path,
    output: Path,
    wave: np.ndarray,
    training_source_rows: np.ndarray,
) -> tuple[pd.DataFrame, sp.csr_matrix, np.ndarray, np.ndarray, tuple[str, ...]]:
    base._init_worker(
        str(stack), str(decomposition), wave, str(output), "catalog", None
    )
    try:
        sample_index = np.unique(
            np.rint(
                np.linspace(
                    0, len(training_source_rows) - 1, LSF_REFERENCE_SAMPLE_ROWS
                )
            ).astype(int)
        )
        science_wave = np.asarray([item[1] for item in SCIENCE_LINES], dtype=float)
        sampled_fwhm = np.asarray(
            [
                evaluate_lsf_diagnostics(
                    global_fit._model(
                        int(training_source_rows[index])
                    ).lsf_surface_state,
                    science_wave,
                )["fwhm_angstrom"]
                for index in sample_index
            ]
        )
        reference_fwhm = np.nanmedian(sampled_fwhm, axis=0)
        if not np.all(np.isfinite(reference_fwhm) & (reference_fwhm > 0.0)):
            raise ValueError("Cannot determine finite positive science-line LSF FWHM")
        first_row = int(training_source_rows[0])
        model = global_fit._model(first_row)
        oh_slice = model._group_slices["oh"]
        oh_line = np.flatnonzero(
            (model._line_group >= oh_slice.start) & (model._line_group < oh_slice.stop)
        )
        full = global_fit._production_oh_catalog(wave).reset_index(drop=True)
        np.testing.assert_allclose(model._line_wave[oh_line], full["wave"])
        np.testing.assert_allclose(model._base_line_weight[oh_line], full["q_aijc_gi"])
        keep = full["wave"].to_numpy(float) >= R_LOWER_A
        catalog = full.loc[keep].copy().reset_index().rename(columns={"index": "oh_index"})
        catalog["model_line_index"] = oh_line[keep]
        catalog["parent_model_index"] = (
            model._line_group[oh_line[keep]] - oh_slice.start
        ).astype(int)
        catalog["q"] = model._base_line_weight[oh_line[keep]]
        protected, reason = line_protection(
            catalog["wave"].to_numpy(float), reference_fwhm
        )
        catalog["protected_science_line"] = protected
        catalog["protection_reason"] = reason
        loading, eligible_count = _individual_loading(catalog)
        selected_line = catalog["model_line_index"].to_numpy(dtype=int)
        selected_parent = catalog["parent_model_index"].to_numpy(dtype=int)
        names = tuple(
            f"OH_{index:03d}" for index in range(oh_slice.stop - oh_slice.start)
        )
        catalog.to_csv(output / "oh_rz_catalog.csv", index=False)
        _atomic_json(
            output / "science_line_protection.json",
            {
                "science_lines": [
                    {
                        "name": name,
                        "air_wavelength_angstrom": wavelength,
                        "reference_fwhm_angstrom": float(reference_fwhm[index]),
                        "sample_min_fwhm_angstrom": float(
                            np.nanmin(sampled_fwhm[:, index])
                        ),
                        "sample_max_fwhm_angstrom": float(
                            np.nanmax(sampled_fwhm[:, index])
                        ),
                    }
                    for index, (name, wavelength) in enumerate(SCIENCE_LINES)
                ],
                "rule": (
                    "Each row masks exactly +/-1 literal FWHM from its fitted 2-D LSF, "
                    "centred on each protected rest wavelength. OH transitions "
                    "within +/-1 median training-sample FWHM of a protected rest wavelength "
                    "are fixed at factor 1."
                ),
                "reference_lsf_rows": int(len(sample_index)),
                "reference_lsf_split": "train",
                "protected_oh_transitions": int(protected.sum()),
                "eligible_oh_transitions": eligible_count,
            },
        )
        return catalog, loading, selected_line, selected_parent, names
    finally:
        base._close_worker_files()


def _init_worker(
    stack: str,
    decomposition: str,
    wave: np.ndarray,
    output: str,
    fingerprint: str,
    selected_line: np.ndarray,
    selected_q: np.ndarray,
    selected_parent: np.ndarray,
    loading: sp.csr_matrix,
    oh_coef_names: tuple[str, ...],
    worker_counter=None,
) -> None:
    global _WAVE, _RZ_PIXEL, _SELECTED_LINE
    global _SELECTED_Q, _SELECTED_PARENT, _LOADING, _OH_COEF_NAMES
    _WAVE = np.asarray(wave, dtype=float)
    _RZ_PIXEL = _WAVE >= R_LOWER_A
    _SELECTED_LINE = selected_line
    _SELECTED_Q = selected_q
    _SELECTED_PARENT = selected_parent
    _LOADING = loading
    _OH_COEF_NAMES = oh_coef_names
    base._init_worker(
        stack,
        decomposition,
        wave,
        output,
        fingerprint,
        None,
        worker_counter=worker_counter,
        pin_cpu=True,
    )


def _row_science_mask(
    source_row: int, model
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the rest-centred row mask and literal fitted-LSF FWHM."""
    assert _WAVE is not None
    centres = np.asarray([item[1] for item in SCIENCE_LINES], dtype=float)
    fwhm = evaluate_lsf_diagnostics(
        model.lsf_surface_state, centres
    )["fwhm_angstrom"]
    if not np.all(np.isfinite(fwhm) & (fwhm > 0.0)):
        raise ValueError(f"Invalid science-line LSF FWHM for source row {source_row}")
    mask = np.zeros(_WAVE.size, dtype=bool)
    for centre, width in zip(centres, fwhm, strict=True):
        mask |= np.abs(_WAVE - centre) <= width
    return mask, centres, fwhm


def _row_design(
    source_row: int,
) -> tuple[np.ndarray, sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    assert _RZ_PIXEL is not None
    assert _SELECTED_LINE is not None and _SELECTED_Q is not None
    assert _SELECTED_PARENT is not None and _LOADING is not None
    observed = (
        np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=float) * FACTOR
    )
    residual = observed - np.asarray(
        base._DECOMP_HDU["BESTFIT_LSF"].data[source_row], dtype=float
    )
    model = global_fit._model(source_row)
    science_mask, _, _ = _row_science_mask(source_row, model)
    raw = global_fit._raw_line_design(model)[:, _SELECTED_LINE][_RZ_PIXEL].tocsr()
    coefficient_row = base._DECOMP_HDU["COEF"].data[source_row]
    oh_coefficient = np.asarray(
        [coefficient_row[name] for name in _OH_COEF_NAMES], dtype=float
    )
    amplitude = _SELECTED_Q * oh_coefficient[_SELECTED_PARENT]
    design = (raw @ sp.diags(amplitude) @ _LOADING).tocsr()
    ivar = decompose._fit_ivar_row("sky2", source_row, observed)[_RZ_PIXEL]
    y = residual[_RZ_PIXEL]
    valid = (
        np.isfinite(y)
        & np.isfinite(ivar)
        & (ivar > 0.0)
        & ~science_mask[_RZ_PIXEL]
        & (np.asarray(raw.getnnz(axis=1)).ravel() > 0)
    )
    return y, design, ivar, valid, science_mask[_RZ_PIXEL]


def _empty_equations(n_parameters: int) -> dict[str, object]:
    return {
        "gram": sp.csr_matrix((n_parameters, n_parameters), dtype=float),
        "rhs": np.zeros(n_parameters, dtype=float),
        "weighted_sum_sq": 0.0,
        "unweighted_sum_sq": 0.0,
        "weight_sum": 0.0,
        "pixels": 0,
        "full_rz_weighted_sum_sq": 0.0,
        "full_rz_weight_sum": 0.0,
        "full_rz_pixels": 0,
        "rows": 0,
    }


def _accumulate_chunk(task: tuple[str, int, np.ndarray]) -> dict[str, object]:
    split, chunk_id, source_rows = task
    assert _LOADING is not None
    result = _empty_equations(_LOADING.shape[1])
    errors = []
    started = time.perf_counter()
    for source_row in source_rows:
        try:
            residual, design, ivar, use, science_mask = _row_design(int(source_row))
            root_weight = np.sqrt(ivar[use])
            x = design[use].multiply(root_weight[:, None]).tocsr()
            y = residual[use]
            yw = y * root_weight
            result["gram"] += x.T @ x
            result["rhs"] += np.asarray(x.T @ yw).ravel()
            result["weighted_sum_sq"] += float(yw @ yw)
            result["unweighted_sum_sq"] += float(y @ y)
            result["weight_sum"] += float(ivar[use].sum())
            result["pixels"] += int(use.sum())
            full_use = (
                np.isfinite(residual)
                & np.isfinite(ivar)
                & (ivar > 0.0)
                & ~science_mask
            )
            result["full_rz_weighted_sum_sq"] += float(
                np.sum(ivar[full_use] * residual[full_use] ** 2)
            )
            result["full_rz_weight_sum"] += float(ivar[full_use].sum())
            result["full_rz_pixels"] += int(full_use.sum())
            result["rows"] += 1
        except Exception as error:
            errors.append(
                {
                    "source_row": int(source_row),
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc(),
                }
            )
    result["gram"].eliminate_zeros()
    return result | {
        "split": split,
        "chunk_id": chunk_id,
        "errors": errors,
        "elapsed_sec": time.perf_counter() - started,
    }


def _pool(tasks, worker, workers: int, initializer: tuple, label: str):
    context = mp.get_context("spawn")
    counter = context.Value("i", 0)
    started = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=_init_worker,
        initargs=initializer + (counter,),
    ) as executor:
        futures = [executor.submit(worker, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), 1):
            yield future.result()
            if completed == 1 or completed % 10 == 0 or completed == len(futures):
                print(
                    f"{label}: chunks={completed}/{len(futures)} "
                    f"elapsed={(time.perf_counter() - started) / 60.0:.1f} min",
                    flush=True,
                )


def accumulate(
    corpus: pd.DataFrame,
    chunk_size: int,
    workers: int,
    initializer: tuple,
    n_parameters: int,
) -> dict[str, dict[str, object]]:
    equations = {
        split: _empty_equations(n_parameters)
        for split in ("train", "validation", "holdout")
    }
    tasks = []
    chunk_id = 0
    for split in equations:
        rows = corpus.loc[corpus["split"] == split, "source_row"].to_numpy(np.int64)
        for start in range(0, len(rows), chunk_size):
            tasks.append((split, chunk_id, rows[start : start + chunk_size]))
            chunk_id += 1
    errors = []
    for result in _pool(tasks, _accumulate_chunk, workers, initializer, "normal"):
        errors.extend(result["errors"])
        target = equations[result["split"]]
        for name in (
            "weighted_sum_sq",
            "unweighted_sum_sq",
            "weight_sum",
            "full_rz_weighted_sum_sq",
            "full_rz_weight_sum",
        ):
            target[name] += float(result[name])
        for name in ("pixels", "full_rz_pixels", "rows"):
            target[name] += int(result[name])
        target["gram"] += result["gram"]
        target["rhs"] += result["rhs"]
    if errors:
        raise RuntimeError(json.dumps(errors[:10], indent=2))
    return equations


def _objective(equation: dict[str, object], theta: np.ndarray) -> dict[str, float]:
    post = (
        float(equation["weighted_sum_sq"])
        - 2.0 * float(theta @ equation["rhs"])
        + float(theta @ (equation["gram"] @ theta))
    )
    baseline = float(equation["weighted_sum_sq"])
    improvement = baseline - post
    full_baseline = float(equation["full_rz_weighted_sum_sq"])
    full_post = full_baseline - improvement
    return {
        "rows": int(equation["rows"]),
        "pixels": int(equation["pixels"]),
        "baseline_weighted_sum_sq": baseline,
        "corrected_weighted_sum_sq": post,
        "chi2_ratio": post / baseline,
        "chi2_reduction_percent": 100.0 * (1.0 - post / baseline),
        "baseline_weighted_rms": float(
            np.sqrt(baseline / float(equation["weight_sum"]))
        ),
        "corrected_weighted_rms": float(
            np.sqrt(post / float(equation["weight_sum"]))
        ),
        "full_rz_pixels": int(equation["full_rz_pixels"]),
        "full_rz_baseline_weighted_sum_sq": full_baseline,
        "full_rz_corrected_weighted_sum_sq": full_post,
        "full_rz_chi2_ratio": full_post / full_baseline,
        "full_rz_chi2_reduction_percent": 100.0 * improvement / full_baseline,
        "full_rz_baseline_weighted_rms": float(
            np.sqrt(full_baseline / float(equation["full_rz_weight_sum"]))
        ),
        "full_rz_corrected_weighted_rms": float(
            np.sqrt(full_post / float(equation["full_rz_weight_sum"]))
        ),
    }


def _solve_projected_ridge(
    gram: sp.csr_matrix,
    rhs: np.ndarray,
    ridge: float,
    loading: sp.csr_matrix,
    catalog: pd.DataFrame,
) -> tuple[np.ndarray, dict[str, object]]:
    """Solve scaled ridge equations, then enforce physical VNF line factors."""
    diagonal = np.asarray(gram.diagonal(), dtype=float)
    usable = np.isfinite(diagonal) & (diagonal > 0.0)
    information_floor = float(np.quantile(diagonal[usable], INFORMATION_QUANTILE))
    active = usable & (diagonal >= information_floor)
    ridge_scale = float(np.median(diagonal[active]))
    normal = gram[active][:, active]
    normal = (
        (normal + normal.T) * 0.5
        + ridge * ridge_scale * sp.eye(active.sum())
    ).tocsc()
    z = sp.linalg.spsolve(normal, rhs[active])
    theta = np.zeros_like(rhs)
    theta[active] = np.asarray(z, dtype=float)

    raw_factor = 1.0 + np.asarray(loading @ theta).ravel()
    protected = catalog["protected_science_line"].to_numpy(bool)
    parent = catalog["parent_index"].to_numpy(int)
    delta = raw_factor - 1.0
    for parent_id in np.unique(parent):
        use = (parent == parent_id) & ~protected
        value = delta[use]
        if value.size < 2 or not np.any(value):
            continue
        bounds = [1.0]
        if np.any(value < -1.0):
            bounds.append(float(np.min(-1.0 / value[value < -1.0])))
        if np.any(value > MAX_CORRECTION_FACTOR - 1.0):
            bounds.append(
                float(
                    np.min(
                        (MAX_CORRECTION_FACTOR - 1.0)
                        / value[value > MAX_CORRECTION_FACTOR - 1.0]
                    )
                )
            )
        delta[use] *= max(0.0, min(bounds))

    columns = loading.tocsc()
    for column in range(columns.shape[1]):
        start, stop = columns.indptr[column : column + 2]
        rows = columns.indices[start:stop]
        positive = rows[columns.data[start:stop] > 0.0]
        if positive.size != 1:
            raise AssertionError("Each contrast must have one non-reference line")
        theta[column] = delta[int(positive[0])]
    np.testing.assert_allclose(loading @ theta, delta, rtol=1.0e-9, atol=1.0e-10)
    curvature = float(theta @ (gram @ theta))
    slope = float(theta @ rhs)
    data_scale = min(1.0, max(0.0, slope / curvature)) if curvature > 0.0 else 0.0
    theta *= data_scale
    factor = 1.0 + np.asarray(loading @ theta).ravel()
    return theta, {
        "solver": "scipy.sparse.linalg.spsolve + per-VNF nonnegative q-closure projection",
        "active_parameters": int(active.sum()),
        "information_quantile": INFORMATION_QUANTILE,
        "information_floor": information_floor,
        "physical_ridge_scale": ridge_scale,
        "raw_minimum_factor": float(np.nanmin(raw_factor)),
        "minimum_factor": float(factor.min()),
        "maximum_factor": float(factor.max()),
        "data_ray_scale": data_scale,
        "projected_factors": int(
            np.count_nonzero(np.abs(factor - raw_factor) > 1.0e-10)
        ),
    }


def fit_ridge_grid(
    equations: dict[str, dict[str, object]],
    loading: sp.csr_matrix,
    catalog: pd.DataFrame,
) -> tuple[np.ndarray, float, pd.DataFrame, dict[str, object]]:
    rows = []
    candidates: dict[float, np.ndarray] = {}
    solver = {}
    for ridge in RIDGE_GRID:
        theta, info = _solve_projected_ridge(
            equations["train"]["gram"],
            equations["train"]["rhs"],
            ridge,
            loading,
            catalog,
        )
        candidates[ridge] = theta
        solver[ridge] = info
        for split in ("train", "validation"):
            rows.append({"ridge_lambda": ridge, "split": split, **_objective(equations[split], theta)})
    screen = pd.DataFrame(rows)
    validation = screen.loc[screen["split"] == "validation"]
    winner = float(validation.loc[validation["chi2_ratio"].idxmin(), "ridge_lambda"])
    development = {
        name: equations["train"][name] + equations["validation"][name]
        for name in (
            "weighted_sum_sq",
            "unweighted_sum_sq",
            "weight_sum",
            "pixels",
            "full_rz_weighted_sum_sq",
            "full_rz_weight_sum",
            "full_rz_pixels",
            "rows",
        )
    }
    development["gram"] = equations["train"]["gram"] + equations["validation"]["gram"]
    development["rhs"] = equations["train"]["rhs"] + equations["validation"]["rhs"]
    theta, final_solver = _solve_projected_ridge(
        development["gram"], development["rhs"], winner, loading, catalog
    )
    details = {
        "selection": "minimum validation row-normalized photon-weighted SSE; MJD-disjoint split",
        "winner_ridge_lambda": winner,
        "screen_solver": {str(key): value for key, value in solver.items()},
        "final_solver": final_solver,
        "active_parameters": int(final_solver["active_parameters"]),
        "parameters": int(theta.size),
    }
    return theta, winner, screen, details


def _binned(values: np.ndarray, width: int) -> np.ndarray:
    return global_fit._binned_nanmean(values, width)


def artifact_proximity_order(baseline: np.ndarray) -> np.ndarray:
    """Cluster normalized residual morphologies and return adjacent artifact families."""
    from scipy.cluster.hierarchy import leaves_list, linkage

    matrix = np.asarray(baseline, dtype=float).copy()
    column_median = np.nanmedian(matrix, axis=0)
    column_median[~np.isfinite(column_median)] = 0.0
    bad = ~np.isfinite(matrix)
    matrix[bad] = np.broadcast_to(column_median, matrix.shape)[bad]
    matrix -= np.median(matrix, axis=0)
    scale = np.sqrt(np.mean(matrix**2, axis=1))
    matrix /= np.where(scale > 0.0, scale, 1.0)[:, None]
    components = min(12, min(matrix.shape) - 1)
    left, singular, _ = sp.linalg.svds(
        sp.csr_matrix(matrix), k=components, rng=np.random.default_rng(0)
    )
    features = left * singular
    return leaves_list(linkage(features, method="ward", optimal_ordering=True))


def _render_chunk(task: tuple[int, np.ndarray, np.ndarray]) -> dict[str, object]:
    chunk_id, source_rows, theta = task
    assert _RZ_PIXEL is not None and _WAVE is not None
    r_band = _WAVE[_RZ_PIXEL] < 7454.0
    result = {
        "chunk_id": chunk_id,
        "source_row": [],
        "baseline": [],
        "corrected": [],
        "science_mask": [],
        "baseline_rms": [],
        "corrected_rms": [],
        "baseline_weighted_rms": [],
        "corrected_weighted_rms": [],
        "mask_max_abs_correction": [],
        "mask_rms_correction": [],
        "mask_nonzero_correction_pixels": [],
        "r_baseline_weighted_sse": [],
        "r_corrected_weighted_sse": [],
        "z_baseline_weighted_sse": [],
        "z_corrected_weighted_sse": [],
        "pixels": [],
        "errors": [],
    }
    for source_row in source_rows:
        try:
            residual, design, ivar, use, science_mask = _row_design(int(source_row))
            correction = np.asarray(design @ theta).ravel()
            corrected = residual - correction
            finite_mask = science_mask & np.isfinite(correction)
            full_use = (
                np.isfinite(residual)
                & np.isfinite(ivar)
                & (ivar > 0.0)
                & ~science_mask
            )
            result["source_row"].append(int(source_row))
            result["baseline"].append(_binned(residual, DISPLAY_BIN_PIXELS))
            result["corrected"].append(_binned(corrected, DISPLAY_BIN_PIXELS))
            result["science_mask"].append(
                _binned(science_mask.astype(float), DISPLAY_BIN_PIXELS) > 0.0
            )
            result["baseline_rms"].append(float(np.sqrt(np.mean(residual[use] ** 2))))
            result["corrected_rms"].append(float(np.sqrt(np.mean(corrected[use] ** 2))))
            weight = ivar[use]
            result["baseline_weighted_rms"].append(
                float(np.sqrt(np.sum(weight * residual[use] ** 2) / np.sum(weight)))
            )
            result["corrected_weighted_rms"].append(
                float(np.sqrt(np.sum(weight * corrected[use] ** 2) / np.sum(weight)))
            )
            result["mask_max_abs_correction"].append(
                float(np.max(np.abs(correction[finite_mask])))
            )
            result["mask_rms_correction"].append(
                float(np.sqrt(np.mean(correction[finite_mask] ** 2)))
            )
            result["mask_nonzero_correction_pixels"].append(
                int(np.count_nonzero(np.abs(correction[finite_mask]) > 1.0e-12))
            )
            for prefix, band in (("r", r_band), ("z", ~r_band)):
                score = full_use & band
                result[f"{prefix}_baseline_weighted_sse"].append(
                    float(np.sum(ivar[score] * residual[score] ** 2))
                )
                result[f"{prefix}_corrected_weighted_sse"].append(
                    float(np.sum(ivar[score] * corrected[score] ** 2))
                )
            result["pixels"].append(int(use.sum()))
        except Exception as error:
            result["errors"].append(
                {
                    "source_row": int(source_row),
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc(),
                }
            )
    for name in ("baseline", "corrected", "science_mask"):
        result[name] = np.asarray(result[name], dtype=np.float32)
    return result


def render_holdout(
    corpus: pd.DataFrame,
    theta: np.ndarray,
    wave: np.ndarray,
    chunk_size: int,
    workers: int,
    initializer: tuple,
    output: Path,
) -> pd.DataFrame:
    holdout = corpus.loc[corpus["split"] == "holdout"].copy()
    rows = holdout["source_row"].to_numpy(np.int64)
    tasks = [
        (index, rows[start : start + chunk_size], theta)
        for index, start in enumerate(range(0, len(rows), chunk_size))
    ]
    parts = []
    errors = []
    for result in _pool(tasks, _render_chunk, workers, initializer, "render"):
        errors.extend(result["errors"])
        parts.append(result)
    if errors:
        raise RuntimeError(json.dumps(errors[:10], indent=2))
    source_row = np.concatenate([np.asarray(part["source_row"], int) for part in parts])
    order = np.argsort(source_row)
    baseline = np.concatenate([part["baseline"] for part in parts], axis=0)[order]
    corrected = np.concatenate([part["corrected"] for part in parts], axis=0)[order]
    science_mask = np.concatenate(
        [part["science_mask"] for part in parts], axis=0
    )[order].astype(bool)
    artifact_order = artifact_proximity_order(baseline)
    artifact_rank = np.empty(len(artifact_order), dtype=int)
    artifact_rank[artifact_order] = np.arange(len(artifact_order))
    metrics = pd.DataFrame(
        {
            "source_row": source_row[order],
            **{
                name: np.concatenate([np.asarray(part[name]) for part in parts])[order]
                for name in (
                    "baseline_rms",
                    "corrected_rms",
                    "baseline_weighted_rms",
                    "corrected_weighted_rms",
                    "mask_max_abs_correction",
                    "mask_rms_correction",
                    "mask_nonzero_correction_pixels",
                    "r_baseline_weighted_sse",
                    "r_corrected_weighted_sse",
                    "z_baseline_weighted_sse",
                    "z_corrected_weighted_sse",
                    "pixels",
                )
            },
        }
    ).merge(
        holdout[["source_row", "expnum", "mjd", "sky_far_label"]],
        on="source_row",
        how="left",
        validate="one_to_one",
    )
    metrics["artifact_rank"] = artifact_rank
    metrics.to_csv(output / "holdout_spectrum_metrics.csv", index=False)
    rz = wave >= R_LOWER_A
    np.savez_compressed(
        output / "holdout_residual_maps.npz",
        source_row=metrics["source_row"].to_numpy(int),
        mjd=metrics["mjd"].to_numpy(int),
        sky_far_label=metrics["sky_far_label"].astype(str).to_numpy(str),
        wave=_binned(wave[rz], DISPLAY_BIN_PIXELS),
        baseline=baseline,
        corrected=corrected,
        difference=baseline - corrected,
        artifact_order=artifact_order,
        protected_pixel=science_mask,
    )
    return metrics


def export_decomposition_examples(
    metrics: pd.DataFrame,
    theta: np.ndarray,
    wave: np.ndarray,
    initializer: tuple,
    output: Path,
) -> Path:
    ranked = metrics.sort_values("artifact_rank").reset_index(drop=True)
    quantiles = np.array([0.05, 0.20, 0.40, 0.60, 0.80, 0.95])
    selected = ranked.iloc[
        np.rint(quantiles * (len(ranked) - 1)).astype(int)
    ].copy()
    arrays: dict[str, list[np.ndarray]] = {
        name: []
        for name in (
            "observed",
            "baseline_bestfit",
            "corrected_bestfit",
            "baseline_residual",
            "corrected_residual",
            "oh_baseline",
            "oh_corrected",
            "moon",
            "zodi",
            "o2",
            "atom",
            "orc",
            "diffuse",
            "minor_molecular",
            "science_mask",
        )
    }
    science_centres: list[np.ndarray] = []
    science_fwhm: list[np.ndarray] = []
    science_lsf_density: list[np.ndarray] = []
    lsf_delta = np.linspace(-5.0, 5.0, 401)
    _init_worker(*initializer)
    try:
        rz = wave >= R_LOWER_A
        for source_row in selected["source_row"].to_numpy(int):
            observed = (
                np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=float)
                * FACTOR
            )
            bestfit = np.asarray(
                base._DECOMP_HDU["BESTFIT_LSF"].data[source_row], dtype=float
            )
            residual, design, _, _, science_mask = _row_design(source_row)
            correction = np.zeros_like(wave)
            correction[rz] = np.asarray(design @ theta).ravel()
            arrays["observed"].append(observed)
            arrays["baseline_bestfit"].append(bestfit)
            arrays["corrected_bestfit"].append(bestfit + correction)
            arrays["baseline_residual"].append(observed - bestfit)
            arrays["corrected_residual"].append(observed - bestfit - correction)
            full_mask = np.zeros_like(wave, dtype=bool)
            full_mask[rz] = science_mask
            arrays["science_mask"].append(full_mask)
            oh = np.asarray(base._DECOMP_HDU["COMP_OH"].data[source_row], dtype=float)
            arrays["oh_baseline"].append(oh)
            arrays["oh_corrected"].append(oh + correction)
            for key, extension in (
                ("moon", "COMP_MOON"),
                ("zodi", "COMP_ZODI"),
                ("o2", "COMP_O2"),
                ("atom", "COMP_ATOM"),
                ("orc", "COMP_ORC"),
                ("diffuse", "COMP_DIFFUSE"),
            ):
                arrays[key].append(
                    np.asarray(base._DECOMP_HDU[extension].data[source_row], dtype=float)
                )
            arrays["minor_molecular"].append(
                sum(
                    np.asarray(base._DECOMP_HDU[name].data[source_row], dtype=float)
                    for name in ("COMP_HO2", "COMP_FEO", "COMP_O2AC")
                )
            )
            model = global_fit._model(source_row)
            _, centres, fwhm = _row_science_mask(source_row, model)
            science_centres.append(centres)
            science_fwhm.append(fwhm)
            science_lsf_density.append(
                evaluate_lsf_density(model.lsf_surface_state, centres, lsf_delta)
            )
    finally:
        base._close_worker_files()
    selected["artifact_quantile"] = quantiles
    selected.to_csv(output / "decomposition_examples.csv", index=False)
    path = output / "decomposition_examples.npz"
    np.savez_compressed(
        path,
        wave=wave,
        lsf_delta=lsf_delta,
        science_centres=np.asarray(science_centres, dtype=np.float64),
        science_fwhm=np.asarray(science_fwhm, dtype=np.float64),
        science_lsf_density=np.asarray(science_lsf_density, dtype=np.float32),
        source_row=selected["source_row"].to_numpy(int),
        **{name: np.asarray(value, dtype=np.float32) for name, value in arrays.items()},
    )
    return path


def export_palace_table(
    output: Path,
    catalog: pd.DataFrame,
    factor: np.ndarray,
    provenance: dict[str, object],
    publish_canonical: bool = True,
) -> tuple[Path, Path, dict[str, object]]:
    source = (
        DEFAULT_DATA_ROOT
        / "palace"
        / "PMD"
        / f"pmd_popmodel_OH{DEFAULT_PALACE_OH_SUFFIX}.dat"
    )
    canonical = source.with_name(f"pmd_popmodel_OH{PALACE_OUTPUT_SUFFIX}.dat")
    artifact = output / canonical.name
    mapping = {
        str(identifier).strip(): float(value)
        for identifier, value in zip(catalog["ID"], factor, strict=True)
    }
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    header_replacements = {
        "# lvmsky derived OH coefficients": (
            "# lvmsky empirical OH coefficients for Skyfar-response calibration."
        ),
        "# Source checkpoint:": f"# Parent PALACE table: {source.name}",
        "# Source export:": f"# Run fingerprint: {provenance['run_fingerprint']}",
        "# Optimization:": (
            "# Calibration: fixed-baseline linear weighted SSE with VNF q-closure; "
            f"ridge lambda={provenance['winner_ridge_lambda']:g}."
        ),
        "# Corpus:": (
            f"# Corpus: {provenance['fit_corpus_rows']} FLUX_SKY_FAR spectra; "
            "MJD-disjoint train/validation/holdout."
        ),
        "# Mapping:": (
            "# Mapping: 22058 rows; only Aijc and I changed; "
            "sum(Aijc * gi) preserved per (vi, Ni, Fi) VNF group."
        ),
    }
    identifiers = []
    rendered = []
    started = False
    changed = 0
    for raw in lines:
        line = raw.rstrip("\r\n")
        ending = raw[len(line) :]
        replacement = next(
            (
                value
                for prefix, value in header_replacements.items()
                if line.startswith(prefix)
            ),
            None,
        )
        if replacement is not None:
            rendered.append(replacement + ending)
            continue
        if line.startswith("lam ID "):
            started = True
            rendered.append(raw)
            continue
        if not line or line.startswith("#") or not started:
            rendered.append(raw)
            continue
        fields = line.split()
        identifier = fields[1].strip()
        identifiers.append(identifier)
        value = mapping.get(identifier, 1.0)
        if identifier not in mapping or value == 1.0:
            rendered.append(raw)
            continue
        replacements = {
            8: f"{float(fields[8]) * value:.12g}",
            11: f"{float(fields[11]) * value:.12g}",
        }
        spans = list(re.finditer(r"\S+", line))
        for index in sorted(replacements, reverse=True):
            span = spans[index].span()
            line = line[: span[0]] + replacements[index] + line[span[1] :]
        rendered.append(line + ending)
        changed += 1
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("PALACE transition IDs are not unique")
    missing = sorted(set(mapping) - set(identifiers))
    if missing:
        raise ValueError(f"Corrected transition IDs missing from PALACE table: {missing[:5]}")
    text = "".join(rendered)
    if len(rendered) != len(lines):
        raise AssertionError("PALACE row structure changed")
    for old, new in zip(lines, rendered, strict=True):
        old_fields, new_fields = old.split(), new.split()
        try:
            float(old_fields[0])
        except (IndexError, ValueError):
            continue
        if len(old_fields) != len(new_fields) or any(
            left != right
            for index, (left, right) in enumerate(zip(old_fields, new_fields, strict=True))
            if index not in (8, 11)
        ):
            raise AssertionError("PALACE fields outside Aijc and I changed")
    destinations = (canonical, artifact) if publish_canonical else (artifact,)
    for destination in destinations:
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, destination)
    validation = {
        "source": str(source.resolve()),
        "canonical_output": str(canonical.resolve()) if publish_canonical else None,
        "published_copy": str(artifact.resolve()),
        "same_line_count": True,
        "rows": len(identifiers),
        "changed_rows": changed,
        "changed_fields": ["Aijc", "I"],
        "all_other_fields_preserved": True,
        "provenance_header_updated": True,
    }
    return canonical if publish_canonical else artifact, artifact, validation


def build_report(output: Path) -> tuple[Path, Path]:
    import nbformat as nbf
    from nbconvert import HTMLExporter
    from nbconvert.preprocessors import ExecutePreprocessor

    notebook_path = output / "oh_skyfar_fwhm_v3.ipynb"
    html_path = output / "oh_skyfar_fwhm_v3.html"
    out = repr(str(output.resolve()))
    cells = [
        nbf.v4.new_markdown_cell(
            "# Skyfar-only PALACE OH calibration with ±1-FWHM protection\n\n"
            "This report estimates empirical corrections to the relative OH line "
            "strengths inside each production `(v_upper, N_upper, F_upper)` VNF group. "
            "It starts from the canonical PALACE `Aijc × gi` ratios used by "
            "`palace-aijc-vnf-split-zodi-lsf-spline2d`. No alternative theoretical "
            "coefficient set is considered.\n\n"
            "**Scope.** This is a fixed-baseline linear replay: VNF amplitudes, telluric "
            "solution, continuum, and 2-D LSF are held at their baseline decomposition "
            "values. The reported improvement is therefore a matched weighted-SSE result, "
            "not an end-to-end production claim. A fresh held-out redecomposition with the "
            "exported table is required before deployment.\n\n"
            "The model is linear in fixed decomposition amplitudes. Training, validation, "
            "and holdout data are exclusively `FLUX_SKY_FAR`, split by MJD. Each VNF "
            "group preserves its total `Aijc × gi`; line factors are non-negative. "
            "Each protected window is centred on its rest wavelength and extends exactly ±1 "
            "literal FWHM of that row's fitted 2-D LSF. OH transitions within ±1 training-reference "
            "FWHM of Hα, [N II], [S II], or [S III] are fixed at factor 1. Ridge strength "
            "is selected on validation nights, then refitted on train+validation and scored "
            "once on held-out nights. In symbols, `r_s = D_s L theta + epsilon`, with "
            "`q_g^T L_g theta = 0` for every VNF group. Here `D_s` is the fixed line-profile "
            "design, `L` maps contrasts to transitions, `theta` contains fitted contrasts, "
            "and `q_g` contains baseline line strengths. Ridge screening activates the top "
            "20% of contrasts by training information; the final train+validation refit "
            "recomputes that threshold without using holdout data. A group-preserving ray "
            "projection limits factors to [0, 4].\n\n"
            "**Terminology.** PALACE is the production airglow-line catalogue; VNF is the "
            "upper-state grouping used by the decomposition; `Aijc` is its corrected "
            "transition coefficient and `gi` is statistical weight. Skyfar is the far-sky "
            "fibre spectrum, and PWV is precipitable water vapour. Residual *fit units* are "
            "the input flux-density units multiplied by `1e14`. The objective uses "
            "row-normalized photon-weight shapes, so it is a weighted sum of squared errors "
            "(SSE), not a calibrated statistical chi-square. For each spectrum, photon weights "
            "are built from absolute Skyfar flux with the production variance floor, then divided "
            "by their mean over the finite full-band row before R/Z and science masks are applied. "
            "The data objective is "
            "`sum(w * (r - D L theta)^2)`, plus `lambda` times the squared active contrasts. "
            "The information score is the diagonal of `X.T W X`, and contrasts at or above its "
            "80th percentile are active. Machine-readable `chi2_*` fields are retained as legacy "
            "aliases for these weighted-SSE ratios."
        ),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import warnings\n"
            "warnings.simplefilter('ignore')\n"
            "import json, numpy as np, pandas as pd, matplotlib.pyplot as plt\n"
            "import plotly.graph_objects as go, plotly.io as pio\n"
            "from plotly.subplots import make_subplots\n"
            "pio.renderers.default = 'notebook_connected'\n"
            f"OUT = Path({out})\n"
            "summary = json.loads((OUT/'summary.json').read_text())\n"
            "corpus = pd.read_csv(OUT/'skyfar_corpus.csv')\n"
            "lines = pd.read_csv(OUT/'oh_line_corrections.csv')\n"
            "screen = pd.read_csv(OUT/'ridge_screen.csv')\n"
            "metrics = pd.read_csv(OUT/'holdout_spectrum_metrics.csv')\n"
            "maps = np.load(OUT/'holdout_residual_maps.npz')\n"
            "example_meta = pd.read_csv(OUT/'decomposition_examples.csv')\n"
            "examples = np.load(OUT/'decomposition_examples.npz')\n"
            "protection = json.loads((OUT/'science_line_protection.json').read_text())"
        ),
        nbf.v4.new_markdown_cell(
            "## Wall-clock computation time\n\n"
            "Times include worker initialization and I/O. The total covers corpus selection, "
            "LSF-reference construction, fitting, holdout rendering, and example export; HTML "
            "rendering itself is outside the scientific-compute total. The run used 64 workers "
            "on the host that served this report."
        ),
        nbf.v4.new_code_cell(
            "timing=pd.Series(summary['timing_seconds'],name='seconds').drop('workers')\n"
            "display(pd.DataFrame({'seconds':timing,'minutes':timing/60}))\n"
            "print(f\"Workers: {int(summary['timing_seconds']['workers'])}\")"
        ),
        nbf.v4.new_markdown_cell(
            "## Corpus and safeguards\n\n"
            "The upstream robust continuum/join filter removes level and arm-join outliers. "
            "The stricter corpus also requires measured positive PWV, astronomical night, "
            "Skyfar airmass ≤ 1.8, and Skyfar fibre retention ≥ 0.8. Unique MJDs are assigned "
            "70%/15%/15% to train/validation/holdout by sorting `SHA256('17:' + MJD)`; therefore "
            "the split is deterministic and MJD-disjoint. Protection is intentionally limited "
            "to the eight lines in the table below. Every other astrophysical line, including "
            "all Ar lines and other weak species, remains unprotected and contributes normally."
        ),
        nbf.v4.new_code_cell(
            "display(pd.DataFrame({\n"
            " 'rows': corpus.groupby('split').size(),\n"
            " 'nights': corpus.groupby('split').mjd.nunique(),\n"
            "}))\n"
            "display(pd.DataFrame({\n"
            " 'OH transitions': [len(lines)],\n"
            " 'eligible': [(~lines.protected_science_line).sum()],\n"
            " 'protected at factor 1': [lines.protected_science_line.sum()],\n"
            " 'VNF groups': [lines.parent_index.nunique()],\n"
            "}))\n"
            "guard=pd.DataFrame(protection['science_lines']).rename(columns={'name':'protected line','air_wavelength_angstrom':'air wavelength (Å)','reference_fwhm_angstrom':'reference ±mask half-width (Å)','sample_min_fwhm_angstrom':'sample minimum FWHM (Å)','sample_max_fwhm_angstrom':'sample maximum FWHM (Å)'})\n"
            "display(guard[['protected line','air wavelength (Å)','reference ±mask half-width (Å)','sample minimum FWHM (Å)','sample maximum FWHM (Å)']])"
        ),
        nbf.v4.new_markdown_cell("## Validation-selected regularization"),
        nbf.v4.new_code_cell(
            "fig, ax = plt.subplots(figsize=(7,4))\n"
            "for split, part in screen.groupby('split'):\n"
            "    ax.plot(part.ridge_lambda, part.chi2_ratio, marker='o', label=split)\n"
            "ax.axvline(summary['winner_ridge_lambda'], color='k', ls='--', lw=1, label='selected')\n"
            "ax.set(xscale='log', xlabel='Physical line-factor ridge strength', ylabel='OH-support weighted-SSE ratio')\n"
            "ax.grid(alpha=.2); ax.legend(); plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## Learned relative intensities\n\n"
            "Protected transitions are shown at exactly one. The export changes only `Aijc` "
            "and the derived reference intensity `I`; `Aij` and all other PALACE fields are unchanged."
        ),
        nbf.v4.new_code_cell(
            "fig, ax = plt.subplots(figsize=(12,4))\n"
            "protected = lines.protected_science_line\n"
            "changed = ~protected & ~np.isclose(lines.correction_factor,1)\n"
            "fixed = ~protected & ~changed\n"
            "ax.scatter(lines.loc[fixed,'wave'], lines.loc[fixed,'correction_factor'], s=4, alpha=.18, label=f'eligible, retained at 1 ({fixed.sum():,})')\n"
            "ax.scatter(lines.loc[changed,'wave'], lines.loc[changed,'correction_factor'], s=7, alpha=.45, label=f'estimated, non-unit ({changed.sum():,})')\n"
            "ax.scatter(lines.loc[protected,'wave'], lines.loc[protected,'correction_factor'], s=16, marker='x', color='crimson', label=f'science-line protected ({protected.sum():,})')\n"
            "ax.axhline(1, color='k', lw=1); ax.set(xlabel='Air wavelength (Å)', ylabel='Aijc correction factor', ylim=(-.05,4.05))\n"
            "ax.grid(alpha=.2); ax.legend(); plt.show()"
        ),
        nbf.v4.new_code_cell(
            "display(pd.DataFrame({\n"
            " 'active contrasts':[summary['active_parameters']],\n"
            " 'projected factors':[summary['final_solver']['projected_factors']],\n"
            " 'factors at 0':[np.isclose(lines.correction_factor,0).sum()],\n"
            " 'factors at 4':[np.isclose(lines.correction_factor,4).sum()],\n"
            " 'changed by >10%':[summary['factor_statistics']['changed_over_10_percent']],\n"
            " 'maximum relative VNF closure error':[summary['factor_statistics']['maximum_relative_vnf_closure_error']],\n"
            "}))"
        ),
        nbf.v4.new_markdown_cell(
            "## Mask width versus the fitted LSF\n\n"
            "The curves are the exact fitted 2-D LSF densities, median-normalized across the "
            "six held-out examples. Gold shading is the mask. Its half-width is one full "
            "literal FWHM, so its total width is 2 × FWHM; the dotted horizontal line marks "
            "half maximum. The displayed FWHM is the six-example median; the fit uses each "
            "row's own value. Masked pixels are excluded from the objective, not guaranteed "
            "to receive zero correction."
        ),
        nbf.v4.new_code_cell(
            "delta=examples['lsf_delta']; density=examples['science_lsf_density']; widths=examples['science_fwhm']\n"
            "fig=make_subplots(rows=4,cols=2,subplot_titles=[item['name'] for item in protection['science_lines']],vertical_spacing=.10)\n"
            "for j,item in enumerate(protection['science_lines']):\n"
            "    row,col=divmod(j,2); profile=density[:,j,:]/np.nanmax(density[:,j,:],axis=1)[:,None]\n"
            "    median_profile=np.nanmedian(profile,axis=0); width=float(np.nanmedian(widths[:,j]))\n"
            "    fig.add_vrect(x0=-width,x1=width,fillcolor='gold',opacity=.18,line_width=0,row=row+1,col=col+1)\n"
            "    fig.add_trace(go.Scattergl(x=delta,y=median_profile,mode='lines',line=dict(color='black',width=1.5),name='Exact fitted LSF',showlegend=(j==0),hovertemplate='Δλ=%{x:.3f} Å<br>normalized LSF=%{y:.3f}<extra></extra>'),row=row+1,col=col+1)\n"
            "    fig.add_hline(y=.5,line_color='gray',line_dash='dot',line_width=.8,row=row+1,col=col+1)\n"
            "    fig.add_vline(x=-width,line_color='darkgoldenrod',line_width=1,row=row+1,col=col+1); fig.add_vline(x=width,line_color='darkgoldenrod',line_width=1,row=row+1,col=col+1)\n"
            "    fig.update_xaxes(title_text=f'Δλ (Å); mask half-width = 1 FWHM = {width:.2f} Å',showspikes=False,row=row+1,col=col+1)\n"
            "    fig.update_yaxes(range=[0,1.05],title_text='Normalized LSF',showspikes=False,row=row+1,col=col+1)\n"
            "fig.add_trace(go.Scatter(x=[None],y=[None],mode='markers',marker=dict(symbol='square',size=12,color='rgba(255,180,0,.35)'),name='Mask: ±1 FWHM (total 2 FWHM)',hoverinfo='skip'),row=1,col=1)\n"
            "fig.update_layout(template='plotly_white',height=1050,hovermode='closest',title='Exact fitted 2-D LSF profiles and science-line masks',legend=dict(orientation='h',y=1.04))\n"
            "fig.show(config={'responsive':True,'displaylogo':False})"
        ),
        nbf.v4.new_markdown_cell(
            "## Interactive standard decomposition examples\n\n"
            "Six held-out Skyfar spectra span the cluster leaf ordering from the 5th "
            "to the 95th percentile. Each interactive Plotly panel shows the observed spectrum, "
            "the original and corrected total models, the corrected OH component, the standard "
            "decomposition components, and vertically shifted residuals. Use the Full/R/Z buttons "
            "or drag to zoom; double-click resets the view."
        ),
        nbf.v4.new_code_cell(
            "plot_wave=examples['wave']\n"
            "styles={\n"
            " 'observed':('Observed Skyfar','black',1.0,None),\n"
            " 'baseline_bestfit':('Baseline total','darkorange',1.2,'dash'),\n"
            " 'corrected_bestfit':('Corrected total','crimson',1.2,None),\n"
            " 'oh_corrected':('OH corrected','#1f77b4',1.0,None),\n"
            " 'moon':('Moon','#9467bd',.8,None), 'zodi':('Zodiacal','#8c564b',.8,None),\n"
            " 'o2':('O₂','#2ca02c',.8,None), 'atom':('Atomic','#17becf',.8,None),\n"
            " 'orc':('Oxygen recombination (ORC)','#bcbd22',.8,None), 'diffuse':('Diffuse','#7f7f7f',.8,None),\n"
            " 'minor_molecular':('HO₂+FeO+O₂AC','#e377c2',.8,None),\n"
            "}\n"
            "for i,row in example_meta.reset_index(drop=True).iterrows():\n"
            "    br=examples['baseline_residual'][i]; cr=examples['corrected_residual'][i]\n"
            "    offset=-1.25*np.nanpercentile(np.abs(np.r_[br,cr]),99)\n"
            "    fig=go.Figure()\n"
            "    for key,(label,color,width,dash) in styles.items():\n"
            "        line=dict(color=color,width=width); line.update(dash=dash or 'solid')\n"
            "        fig.add_trace(go.Scattergl(x=plot_wave,y=examples[key][i],mode='lines',name=label,line=line,hovertemplate='%{x:.2f} Å<br>%{y:.4g}<extra>'+label+'</extra>'))\n"
            "    fig.add_trace(go.Scattergl(x=plot_wave,y=br+offset,customdata=br,mode='lines',name='Baseline residual + offset',line=dict(color='royalblue',width=.8),hovertemplate='%{x:.2f} Å<br>residual=%{customdata:.4g}<br>display offset='+f'{offset:.4g}'+'<extra>baseline residual</extra>'))\n"
            "    fig.add_trace(go.Scattergl(x=plot_wave,y=cr+offset,customdata=cr,mode='lines',name='Corrected residual + offset',line=dict(color='seagreen',width=.8),hovertemplate='%{x:.2f} Å<br>residual=%{customdata:.4g}<br>display offset='+f'{offset:.4g}'+'<extra>corrected residual</extra>'))\n"
            "    fig.add_hline(y=offset,line_width=.6,line_color='gray')\n"
            "    for centre,width in zip(examples['science_centres'][i],examples['science_fwhm'][i]):\n"
            "        fig.add_vrect(x0=centre-width,x1=centre+width,fillcolor='gold',opacity=.16,line_width=.7,line_color='darkgoldenrod')\n"
            "    fig.add_trace(go.Scatter(x=[None],y=[None],mode='markers',marker=dict(symbol='square',size=11,color='rgba(255,180,0,.35)'),name='Mask: ±1 row FWHM',hoverinfo='skip'))\n"
            "    buttons=[dict(label='Full',method='relayout',args=[{'xaxis.range':[float(plot_wave[0]),float(plot_wave[-1])]}]),dict(label='R',method='relayout',args=[{'xaxis.range':[5787,7454]}]),dict(label='Z',method='relayout',args=[{'xaxis.range':[7454,float(plot_wave[-1])]}]),dict(label='SIII 6312',method='relayout',args=[{'xaxis.range':[6303,6321]}]),dict(label='Hα/NII/SII',method='relayout',args=[{'xaxis.range':[6538,6741]}]),dict(label='SIII 9069',method='relayout',args=[{'xaxis.range':[9059,9078]}]),dict(label='SIII 9531',method='relayout',args=[{'xaxis.range':[9521,9540]}])]\n"
            "    fig.update_layout(template='plotly_white',height=620,hovermode='closest',margin=dict(l=60,r=210,t=105,b=50),title=f\"EXP {int(row.expnum)} · MJD {int(row.mjd)} · cluster leaf rank {int(row.artifact_rank)} ({row.artifact_quantile:.0%})\",xaxis_title='Air wavelength (Å)',yaxis_title='Flux density × 10¹⁴',legend=dict(orientation='v',x=1.01,y=1,font=dict(size=9)),updatemenus=[dict(type='buttons',direction='right',x=0,xanchor='left',y=1.13,buttons=buttons)])\n"
            "    fig.update_xaxes(showspikes=False)\n"
            "    fig.update_yaxes(showspikes=False)\n"
            "    fig.show(config={'responsive':True,'displaylogo':False})"
        ),
        nbf.v4.new_markdown_cell(
            "## Held-out Skyfar residuals\n\n"
            "All scores below use nights that were unavailable to coefficient estimation and "
            "ridge selection. They measure the fixed-baseline replay described above. The "
            "objective excludes the exact row-specific mask extending one literal fitted FWHM "
            "to either side of each protected rest wavelength. Protected OH transitions within "
            "a training-reference ±1-FWHM interval retain factor 1 exactly. This prevents science-line "
            "pixels from teaching the coefficients; it is not a zero-change constraint inside the "
            "window, because corrected neighbouring OH profiles can still contribute there. "
            "OH-support pixels are pixels reached by at least one rendered OH line profile."
        ),
        nbf.v4.new_code_cell(
            "h=summary['split_metrics']['holdout']\n"
            "display(pd.DataFrame({\n"
            " 'scope':['OH-support pixels','All unprotected R/Z pixels'],\n"
            " 'baseline weighted SSE':[h['baseline_weighted_sum_sq'],h['full_rz_baseline_weighted_sum_sq']],\n"
            " 'corrected weighted SSE':[h['corrected_weighted_sum_sq'],h['full_rz_corrected_weighted_sum_sq']],\n"
            " 'SSE ratio':[h['chi2_ratio'],h['full_rz_chi2_ratio']],\n"
            " 'reduction (%)':[h['chi2_reduction_percent'],h['full_rz_chi2_reduction_percent']],\n"
            "}))\n"
            "arm=pd.DataFrame({\n"
            " 'arm':['R: 5787–7454 Å','Z: 7454–9800 Å'],\n"
            " 'baseline weighted SSE':[metrics.r_baseline_weighted_sse.sum(),metrics.z_baseline_weighted_sse.sum()],\n"
            " 'corrected weighted SSE':[metrics.r_corrected_weighted_sse.sum(),metrics.z_corrected_weighted_sse.sum()],\n"
            "})\n"
            "arm['SSE ratio']=arm['corrected weighted SSE']/arm['baseline weighted SSE']\n"
            "arm['reduction (%)']=100*(1-arm['SSE ratio'])\n"
            "display(arm)\n"
            "ratio=metrics.corrected_weighted_rms/metrics.baseline_weighted_rms\n"
            "night=pd.DataFrame({'mjd':metrics.mjd,'ratio':ratio}).groupby('mjd').ratio.median()\n"
            "display(pd.DataFrame({\n"
            " 'held-out nights':[night.size],\n"
            " 'nights with median weighted-RMS ratio < 1':[(night<1).sum()],\n"
            " 'night-median ratio, 10th/50th/90th percentile':[' / '.join(f'{x:.3f}' for x in night.quantile([.1,.5,.9]))],\n"
            " 'maximum correction inside masked pixels':[metrics.mask_max_abs_correction.max()],\n"
            " 'masked pixels with |correction| > 1e-12':[metrics.mask_nonzero_correction_pixels.sum()],\n"
            "}))"
        ),
        nbf.v4.new_code_cell(
            "wave=maps['wave']; b=maps['baseline']; c=maps['corrected']; p=maps['protected_pixel']\n"
            "fig, ax = plt.subplots(figsize=(13,4))\n"
            "ax.plot(wave, np.nanmedian(b,axis=0), lw=.8, label='baseline')\n"
            "ax.plot(wave, np.nanmedian(c,axis=0), lw=.8, label='corrected')\n"
            "for index,item in enumerate(protection['science_lines']):\n"
            "    x=item['air_wavelength_angstrom']; half=item['reference_fwhm_angstrom']\n"
            "    if x+half>=wave[0]: ax.axvspan(x-half,x+half,color='gold',alpha=.16,label='Reference ±1-FWHM window' if index==0 else None)\n"
            "ax.set(xlabel='Air wavelength (Å)', ylabel='Median residual (input flux-density × 10¹⁴)')\n"
            "ax.grid(alpha=.2); ax.legend(); plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "### Residual maps ordered by residual-pattern clustering\n\n"
            "Rows are ordered from their residual fingerprints, not by MJD: wavelength-centered "
            "baseline residuals are normalized by row RMS, compressed to 12 singular-vector "
            "scores, and ordered with Ward linkage plus optimal leaf ordering. Adjacent rows "
            "therefore have similar artifact morphology; rank direction and percentile are "
            "arbitrary leaf coordinates, not artifact severity."
        ),
        nbf.v4.new_code_cell(
            "order=maps['artifact_order']; d=b-c\n"
            "v=np.nanquantile(np.abs(np.r_[b.ravel(),c.ravel()]), .995)\n"
            "vd=np.nanquantile(np.abs(d), .995)\n"
            "fig, axes=plt.subplots(3,1,figsize=(14,10),sharex=True,constrained_layout=True)\n"
            "for ax,z,title,lim in zip(axes,[b[order],c[order],d[order]],['Baseline residual','Corrected residual','Removed OH structure'],[v,v,vd]):\n"
            "    im=ax.imshow(z,aspect='auto',origin='lower',extent=[wave[0],wave[-1],0,len(order)],cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True)\n"
            "    mask_image=np.where(p[order],1.0,np.nan)\n"
            "    ax.imshow(mask_image,aspect='auto',origin='lower',extent=[wave[0],wave[-1],0,len(order)],cmap='autumn',vmin=0,vmax=1,alpha=.18,interpolation='nearest',rasterized=True)\n"
            "    ax.set(ylabel='Spectrum index in cluster leaf order',title=title); fig.colorbar(im,ax=ax,pad=.01,label='Residual (input flux-density × 10¹⁴)')\n"
            "axes[-1].set_xlabel('Air wavelength (Å)'); plt.show()\n"
            "print('Gold overlay is the exact row-specific ±1-FWHM mask after four-pixel display binning.')"
        ),
        nbf.v4.new_markdown_cell(
            "## Candidate table for end-to-end validation\n\n"
            "`pmd_popmodel_OH_skyfar_linear_ridge_0p1_v1.dat` changes "
            "only the relative PALACE table; use suffix "
            "`_skyfar_linear_ridge_0p1_v1`. The decomposition algorithm and VNF "
            "grouping remain unchanged. A fresh held-out production decomposition is the "
            "required deployment check because VNF amplitudes and other nuisance terms will be "
            "refitted there. The mask-and-transition-freeze contract applies only to Hα, [N II] 6548/6583, "
            "[S II] 6716/6731, and [S III] 6312/9069/9531; all other astrophysical lines are "
            "deliberately unprotected. Run with "
            "`--fit-model palacecorr-aijc-vnf-split-zodi-lsf-spline2d` "
            "and verify that the decomposition metadata records this suffix before comparison."
        ),
        nbf.v4.new_code_cell(
            "guard_counts=[]\n"
            "for item in protection['science_lines']:\n"
            "    use=lines.protection_reason.str.contains(item['name'],regex=False,na=False)\n"
            "    guard_counts.append({'protected line':item['name'],'air wavelength (Å)':item['air_wavelength_angstrom'],'protected OH transitions':int(use.sum()),'minimum factor':lines.loc[use,'correction_factor'].min(),'maximum factor':lines.loc[use,'correction_factor'].max()})\n"
            "display(pd.DataFrame(guard_counts))\n"
            "print('Canonical table:', summary['exported_palace_table'])\n"
            "print('Published copy:', summary['published_palace_table'])"
        ),
    ]
    notebook = nbf.v4.new_notebook(cells=cells)
    notebook.metadata.kernelspec = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    executor = ExecutePreprocessor(timeout=900, kernel_name="python3")
    executor.preprocess(notebook, {"metadata": {"path": str(output)}})
    nbf.write(notebook, notebook_path)
    exporter = HTMLExporter(template_name="lab")
    exporter.exclude_input = True
    body, _ = exporter.from_notebook_node(notebook)
    body = body.replace("<title>Notebook</title>", "<title>Skyfar-only PALACE OH calibration</title>", 1)
    temporary = html_path.with_suffix(".tmp")
    temporary.write_text(body, encoding="utf-8")
    os.replace(temporary, html_path)
    return notebook_path, html_path


def self_test() -> None:
    wave = np.array([6550.0, 6562.8, 7135.79, 8950.0, 9068.6, 9500.0, 9530.6])
    protected, _ = line_protection(wave, np.full(len(SCIENCE_LINES), 1.5))
    np.testing.assert_array_equal(
        protected, [False, True, False, False, True, False, True]
    )
    order = artifact_proximity_order(np.array([[0, 1, 0], [2, 0, -1], [0, 1.1, 0], [-2, 0, 1.1]]))
    np.testing.assert_array_equal(np.sort(order), np.arange(4))
    catalog = pd.DataFrame(
        {
            "parent_index": [0, 0, 0, 1, 1],
            "q_aijc_gi": [1.0, 2.0, 3.0, 4.0, 1.0],
            "protected_science_line": [False, True, False, False, False],
        }
    )
    loading, _ = _individual_loading(catalog)
    theta = np.arange(loading.shape[1], dtype=float) / 10.0
    delta = np.asarray(loading @ theta).ravel()
    assert delta[1] == 0.0
    for parent in (0, 1):
        use = catalog["parent_index"].to_numpy() == parent
        np.testing.assert_allclose(
            catalog.loc[use, "q_aijc_gi"].to_numpy() @ delta[use], 0.0, atol=1e-14
        )
    solved, _ = _solve_projected_ridge(
        sp.eye(loading.shape[1], format="csr"),
        np.full(loading.shape[1], 0.2),
        0.1,
        loading,
        catalog,
    )
    solved_factor = 1.0 + np.asarray(loading @ solved).ravel()
    assert solved_factor.min() >= -1.0e-12
    assert solved_factor.max() <= MAX_CORRECTION_FACTOR + 1.0e-12
    equation = _empty_equations(loading.shape[1])
    equation.update(
        gram=sp.eye(loading.shape[1], format="csr"),
        rhs=np.full(loading.shape[1], 0.2),
        weighted_sum_sq=10.0,
        weight_sum=5.0,
        pixels=5,
        full_rz_weighted_sum_sq=20.0,
        full_rz_weight_sum=10.0,
        full_rz_pixels=10,
        rows=1,
    )
    metric = _objective(equation, solved)
    assert metric["full_rz_corrected_weighted_sum_sq"] <= 20.0
    print("self-test: ok")


def run(args: argparse.Namespace) -> None:
    run_started = time.perf_counter()
    timing: dict[str, float] = {}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stage_started = time.perf_counter()
    corpus, wave, corpus_summary = build_corpus(
        args.stack, args.decomposition, args.output_dir, args.qc_sigma
    )
    timing["corpus_selection_seconds"] = time.perf_counter() - stage_started
    if args.max_rows is not None and len(corpus) > args.max_rows:
        pieces = []
        for split in ("train", "validation", "holdout"):
            part = corpus.loc[corpus["split"] == split]
            take = max(1, args.max_rows * len(part) // len(corpus))
            pieces.append(part.iloc[:take])
        corpus = pd.concat(pieces).sort_values("source_row").reset_index(drop=True)
        corpus.to_csv(args.output_dir / "skyfar_corpus.csv", index=False)
    stage_started = time.perf_counter()
    catalog, loading, selected_line, selected_parent, oh_names = build_catalog(
        args.stack,
        args.decomposition,
        args.output_dir,
        wave,
        corpus.loc[corpus["split"] == "train", "source_row"].to_numpy(int),
    )
    timing["catalog_and_lsf_reference_seconds"] = time.perf_counter() - stage_started
    selected_q = catalog["q"].to_numpy(float)
    provenance = {
        "schema": SCHEMA,
        "source_sha256": file_sha256(Path(__file__)),
        "stack": _file_identity(args.stack),
        "decomposition": _file_identity(args.decomposition),
        "wave_sha256": wave_sha256(wave),
        "baseline_palace_oh_suffix": DEFAULT_PALACE_OH_SUFFIX,
        "external_theoretical_coefficients": "not considered",
        "data_roles": ["sky2 / FLUX_SKY_FAR"],
        "fit_corpus_rows": int(len(corpus)),
        "corpus": corpus_summary,
        "ridge_grid": list(RIDGE_GRID),
        "active_information_quantile": INFORMATION_QUANTILE,
        "correction_factor_bounds": [0.0, MAX_CORRECTION_FACTOR],
        "science_line_catalog_size": len(SCIENCE_LINES),
        "science_mask": "row-specific +/-1 literal fitted 2-D-LSF FWHM, rest-centred",
        "protected_oh_transitions": int(catalog["protected_science_line"].sum()),
        "eligible_oh_transitions": int((~catalog["protected_science_line"]).sum()),
        "linear_model": "fixed baseline VNF amplitudes and exact telluric+2-D-LSF line profiles",
    }
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    provenance["run_fingerprint"] = fingerprint
    _atomic_json(args.output_dir / "run_provenance.json", provenance)
    initializer = (
        str(args.stack),
        str(args.decomposition),
        wave,
        str(args.output_dir),
        fingerprint,
        selected_line,
        selected_q,
        selected_parent,
        loading,
        oh_names,
    )
    stage_started = time.perf_counter()
    equations = accumulate(
        corpus, args.chunk_size, args.workers, initializer, loading.shape[1]
    )
    timing["normal_equation_accumulation_seconds"] = time.perf_counter() - stage_started
    stage_started = time.perf_counter()
    theta, ridge, screen, fit_details = fit_ridge_grid(equations, loading, catalog)
    timing["ridge_selection_and_refit_seconds"] = time.perf_counter() - stage_started
    screen.to_csv(args.output_dir / "ridge_screen.csv", index=False)
    np.savez_compressed(
        args.output_dir / "solution.npz",
        theta=theta,
        ridge_lambda=np.asarray(ridge),
    )
    delta = np.asarray(loading @ theta).ravel()
    factor = 1.0 + delta
    protected = catalog["protected_science_line"].to_numpy(bool)
    factor[protected] = 1.0
    if np.any(factor < -1.0e-7) or not np.all(np.isfinite(factor)):
        raise ValueError("Final line factors violate non-negativity or finiteness")
    factor = np.maximum(factor, 0.0)
    line_output = catalog.copy()
    line_output["correction_delta"] = factor - 1.0
    line_output["correction_factor"] = factor
    line_output.to_csv(args.output_dir / "oh_line_corrections.csv", index=False)
    closure = (
        line_output.assign(weighted_delta=line_output["q"] * (factor - 1.0))
        .groupby("parent_index")
        .agg(weighted_delta=("weighted_delta", "sum"), q_sum=("q", "sum"))
    )
    closure["relative_error"] = closure["weighted_delta"].abs() / closure["q_sum"]
    canonical_table, published_table, table_validation = export_palace_table(
        args.output_dir,
        catalog,
        factor,
        provenance | fit_details,
        publish_canonical=args.max_rows is None,
    )
    stage_started = time.perf_counter()
    metrics = render_holdout(
        corpus,
        theta,
        wave,
        args.chunk_size,
        args.workers,
        initializer,
        args.output_dir,
    )
    timing["holdout_render_seconds"] = time.perf_counter() - stage_started
    stage_started = time.perf_counter()
    examples = export_decomposition_examples(
        metrics, theta, wave, initializer, args.output_dir
    )
    timing["example_export_seconds"] = time.perf_counter() - stage_started
    timing["science_compute_total_seconds"] = time.perf_counter() - run_started
    timing["workers"] = float(args.workers)
    _atomic_json(args.output_dir / "timing.json", timing)
    split_metrics = {
        split: _objective(equations[split], theta)
        for split in ("train", "validation", "holdout")
    }
    summary = provenance | fit_details | {
        "winner_ridge_lambda": ridge,
        "split_metrics": split_metrics,
        "factor_statistics": {
            "minimum": float(factor.min()),
            "median": float(np.median(factor)),
            "maximum": float(factor.max()),
            "at_lower_bound": int(np.count_nonzero(np.isclose(factor, 0.0))),
            "at_upper_bound": int(np.count_nonzero(np.isclose(factor, MAX_CORRECTION_FACTOR))),
            "changed_over_10_percent": int(np.count_nonzero(np.abs(factor - 1.0) > 0.1)),
            "protected_exactly_one": bool(np.all(factor[protected] == 1.0)),
            "maximum_relative_vnf_closure_error": float(closure["relative_error"].max()),
        },
        "heldout_spectra_improved": int(
            np.count_nonzero(metrics["corrected_weighted_rms"] < metrics["baseline_weighted_rms"])
        ),
        "heldout_spectra": int(len(metrics)),
        "exported_palace_table": str(canonical_table.resolve()),
        "published_palace_table": str(published_table.resolve()),
        "palace_structure_validation": table_validation,
        "decomposition_examples": str(examples.resolve()),
        "deployment_suffix": PALACE_OUTPUT_SUFFIX,
        "timing_seconds": timing,
    }
    _atomic_json(args.output_dir / "summary.json", summary)
    notebook, html = build_report(args.output_dir)
    print(json.dumps({"summary": summary, "notebook": str(notebook), "html": str(html)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stack", nargs="?", type=Path)
    parser.add_argument("decomposition", nargs="?", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, default=64)
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--qc-sigma", type=float, default=5.0)
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.report_only:
        if args.output_dir is None:
            parser.error("--report-only requires --output-dir")
        build_report(args.output_dir)
        return
    if args.stack is None or args.decomposition is None or args.output_dir is None:
        parser.error("stack, decomposition, and --output-dir are required")
    if args.workers < 1 or args.chunk_size < 1:
        raise ValueError("workers and chunk-size must be positive")
    run(args)


if __name__ == "__main__":
    main()
