"""Fit global linear OH branching-ratio corrections to Far-Sky residuals."""

# ruff: noqa: E402 -- native thread limits must precede NumPy imports.

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
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
):
    os.environ[name] = "1"

import numpy as np
import pandas as pd
import scipy.sparse as sp
import clarabel
from astropy.io import fits

import skysub.decompose_parallel as decompose
from skysub.experiments_niv_integrated_decomposition_v1 import train_pca30 as base
from skysub.sky_decomp.fit import (
    CAP_WAVE,
    LSF_CHANNELS,
    decode_hitran_id,
    read_static_table,
    vac_to_air,
)
from skysub.sky_decomp.lsf_spline2d import (
    LSF_OFFSET_BASIS_COUNT,
    _component_masses,
)
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_OH_SUFFIX,
    file_sha256,
    wave_sha256,
)
from skysub.sky_decomp.result_io import load_lsf_surface_state


R_LOWER_A = 5787.0
RIDGE_LAMBDA = 1.0e-4
DISPLAY_BIN_PIXELS = 4
SCHEMA = "global-oh-branch-corrections-cd-v1"
VARIANTS = ("C", "D")

LOADINGS: dict[str, sp.csr_matrix] = {}
SELECTED_LINE: np.ndarray | None = None
SELECTED_Q: np.ndarray | None = None
SELECTED_PARENT: np.ndarray | None = None
RZ_PIXEL: np.ndarray | None = None
OH_COEF_NAMES: tuple[str, ...] = ()
WAVE: np.ndarray | None = None


def _fingerprint(payload: dict[str, object]) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _load_quality(path: Path, max_rows: int | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame = frame.loc[frame["keep_for_pca"].astype(bool)].sort_values("source_row")
    if max_rows is not None:
        frame = frame.iloc[:max_rows]
    return frame.reset_index(drop=True)


def _raw_line_design(model) -> sp.csc_matrix:
    """Return exact non-normalized profiles for the model's ordered line catalog."""
    if model.lsf_surface_state is None:
        raise ValueError("An LSF surface is required for line profiles")
    transmission = model._line_transmission(model._line_wave)
    n_lines = model._line_wave.size
    blocks = []
    for channel, _, _ in LSF_CHANNELS:
        indices = model._line_indices[channel]
        masses = _component_masses(
            model.lsf_surface_state, channel, model._line_wave[indices]
        )
        transform = sp.coo_matrix(
            (
                (masses * transmission[indices, None]).ravel(),
                (
                    np.arange(indices.size * LSF_OFFSET_BASIS_COUNT),
                    np.repeat(indices, LSF_OFFSET_BASIS_COUNT),
                ),
            ),
            shape=(indices.size * LSF_OFFSET_BASIS_COUNT, n_lines),
        ).tocsc()
        blocks.append(model._line_components[channel] @ transform)
    return sum(blocks[1:], start=blocks[0]).tocsc()


def _model(source_row: int):
    model = decompose._telluric_decomposer("sky2", source_row)
    model._set_lsf_state(load_lsf_surface_state(base.DECOMP_PATH, source_row))
    return model


def _production_oh_catalog(wave: np.ndarray) -> pd.DataFrame:
    path = (
        DEFAULT_DATA_ROOT
        / "palace"
        / "PMD"
        / f"pmd_popmodel_OH{DEFAULT_PALACE_OH_SUFFIX}.dat"
    )
    table = read_static_table(path)
    table["wave"] = vac_to_air(np.asarray(table["lam"], float) * 1.0e4)
    keep = (table["wave"] >= wave.min() - CAP_WAVE) & (
        table["wave"] <= wave.max() + CAP_WAVE
    )
    grouped = decode_hitran_id(table[keep]).group_by(
        ["v_upper", "N_upper", "F_upper"]
    )
    frame = grouped.to_pandas()
    for column in ("ID", "branch_N", "branch_J", "parity"):
        frame[column] = frame[column].astype(str).str.strip()
    frame["q_aijc_gi"] = frame["Aijc"].astype(float) * frame["gi"].astype(float)
    frame["parent_index"] = (
        frame.groupby(["v_upper", "N_upper", "F_upper"], sort=False)
        .ngroup()
        .astype(int)
    )
    frame["branch_family"] = (
        frame["branch_N"]
        + frame["branch_J"]
        + frame["F_upper"].astype(str)
        + frame["F_lower"].astype(str)
    )
    frame["N_bin"] = np.select(
        [frame["N_upper"] <= 4, frame["N_upper"] <= 8],
        ["01-04", "05-08"],
        default="09+",
    )
    frame["band"] = frame["v_upper"].astype(str) + "-" + frame["v_lower"].astype(str)
    return frame


def _group_loading(
    catalog: pd.DataFrame, key_columns: list[str]
) -> tuple[sp.csr_matrix, pd.DataFrame, np.ndarray]:
    """Build q-centred line-to-group contrasts within every VNF parent."""
    keys = pd.MultiIndex.from_frame(catalog[key_columns])
    group_id, unique = pd.factorize(keys, sort=True)
    n_lines, n_groups = len(catalog), len(unique)
    parent = catalog["parent_index"].to_numpy(dtype=int)
    q = catalog["q_aijc_gi"].to_numpy(dtype=float)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for parent_id in np.unique(parent):
        line_index = np.flatnonzero(parent == parent_id)
        local_group = group_id[line_index]
        total = float(q[line_index].sum())
        if not total > 0.0:
            raise ValueError(f"Non-positive Aijc*gi total for parent {parent_id}")
        groups, inverse = np.unique(local_group, return_inverse=True)
        shares = np.bincount(inverse, weights=q[line_index], minlength=groups.size) / total
        for line in line_index:
            rows.append(int(line))
            columns.append(int(group_id[line]))
            values.append(1.0)
        for line in line_index:
            rows.extend([int(line)] * groups.size)
            columns.extend(groups.tolist())
            values.extend((-shares).tolist())
    loading = sp.coo_matrix(
        (values, (rows, columns)), shape=(n_lines, n_groups)
    ).tocsr()
    loading.sum_duplicates()
    loading.eliminate_zeros()
    centred = np.asarray(q @ loading).ravel()
    if not np.allclose(centred, 0.0, atol=1.0e-12 * max(q.sum(), 1.0)):
        raise AssertionError("Global q-weighted loading is not centred")
    for parent_id in np.unique(parent):
        use = parent == parent_id
        if not np.allclose(
            np.asarray(q[use] @ loading[use]).ravel(),
            0.0,
            atol=1.0e-12 * max(q[use].sum(), 1.0),
        ):
            raise AssertionError(f"Parent {parent_id} loading is not centred")
    groups = unique.to_frame(index=False)
    groups.columns = key_columns
    groups.insert(0, "group_id", np.arange(n_groups, dtype=int))
    return loading, groups, group_id.astype(int)


def build_catalog_and_loadings(
    wave: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, sp.csr_matrix], dict[str, pd.DataFrame]]:
    full = _production_oh_catalog(wave)
    selected = (full["wave"] >= R_LOWER_A) & (full["wave"] <= wave.max())
    catalog = full.loc[selected].copy().reset_index().rename(columns={"index": "oh_index"})
    definitions = {
        "C": ["v_upper", "v_lower", "branch_family", "N_bin"],
        "D": ["v_upper", "N_upper", "F_upper", "branch_family"],
    }
    loadings, groups = {}, {}
    for variant, columns in definitions.items():
        loadings[variant], groups[variant], ids = _group_loading(catalog, columns)
        catalog[f"group_{variant}"] = ids
    return catalog, loadings, groups


def _init_worker(
    stack_path: str,
    decomposition_path: str,
    wave: np.ndarray,
    output_dir: str,
    run_fingerprint: str,
    loadings: dict[str, sp.csr_matrix],
    selected_line: np.ndarray,
    selected_q: np.ndarray,
    selected_parent: np.ndarray,
    rz_pixel: np.ndarray,
    oh_coef_names: tuple[str, ...],
    worker_counter=None,
) -> None:
    global LOADINGS, SELECTED_LINE, SELECTED_Q, SELECTED_PARENT, RZ_PIXEL
    global OH_COEF_NAMES, WAVE
    LOADINGS = loadings
    SELECTED_LINE = selected_line
    SELECTED_Q = selected_q
    SELECTED_PARENT = selected_parent
    RZ_PIXEL = rz_pixel
    OH_COEF_NAMES = oh_coef_names
    WAVE = np.asarray(wave, dtype=float)
    base._init_worker(
        stack_path,
        decomposition_path,
        wave,
        output_dir,
        run_fingerprint,
        None,
        worker_counter=worker_counter,
        pin_cpu=True,
    )


def _row_design(source_row: int) -> tuple[np.ndarray, dict[str, sp.csr_matrix], np.ndarray]:
    assert SELECTED_LINE is not None
    assert SELECTED_Q is not None
    assert SELECTED_PARENT is not None
    assert RZ_PIXEL is not None
    observed = (
        np.asarray(decompose._WORKER_FLUX["sky2"][source_row], dtype=float)
        * decompose._WORKER_FACTOR
    )
    residual = observed - np.asarray(
        base._DECOMP_HDU["BESTFIT_LSF"].data[source_row], dtype=float
    )
    model = _model(source_row)
    raw = _raw_line_design(model)[:, SELECTED_LINE][RZ_PIXEL].tocsr()
    coefficient_row = base._DECOMP_HDU["COEF"].data[source_row]
    oh_coefficient = np.asarray(
        [coefficient_row[name] for name in OH_COEF_NAMES], dtype=float
    )
    amplitude = SELECTED_Q * oh_coefficient[SELECTED_PARENT]
    weighted = raw @ sp.diags(amplitude)
    designs = {
        variant: (weighted @ LOADINGS[variant]).tocsr() for variant in VARIANTS
    }
    science = decompose._science_line_mask_for_row(source_row)
    valid = np.isfinite(residual[RZ_PIXEL])
    if science is not None:
        valid &= ~science[RZ_PIXEL]
    support = np.asarray(raw.getnnz(axis=1)).ravel() > 0
    return residual[RZ_PIXEL], designs, valid & support


def _accumulate_chunk(task: tuple[int, np.ndarray]) -> dict[str, object]:
    chunk_id, source_rows = task
    grams = {
        variant: sp.csr_matrix((LOADINGS[variant].shape[1],) * 2, dtype=float)
        for variant in VARIANTS
    }
    rhs = {
        variant: np.zeros(LOADINGS[variant].shape[1], dtype=float)
        for variant in VARIANTS
    }
    residual_sum_sq = 0.0
    n_pixels = 0
    errors = []
    started = time.perf_counter()
    for source_row in source_rows:
        try:
            residual, designs, use = _row_design(int(source_row))
            y = residual[use]
            residual_sum_sq += float(y @ y)
            n_pixels += int(use.sum())
            for variant in VARIANTS:
                x = designs[variant][use]
                grams[variant] += x.T @ x
                rhs[variant] += np.asarray(x.T @ y).ravel()
        except Exception as error:
            errors.append(
                {
                    "source_row": int(source_row),
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc(),
                }
            )
    for gram in grams.values():
        gram.eliminate_zeros()
    return {
        "chunk_id": chunk_id,
        "rows": len(source_rows),
        "grams": grams,
        "rhs": rhs,
        "residual_sum_sq": residual_sum_sq,
        "n_pixels": n_pixels,
        "errors": errors,
        "elapsed_sec": time.perf_counter() - started,
    }


def _run_pool(
    tasks,
    worker,
    workers: int,
    initializer_args: tuple,
    label: str,
):
    context = mp.get_context("spawn")
    counter = context.Value("i", 0)
    started = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=_init_worker,
        initargs=initializer_args + (counter,),
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


def accumulate_normal_equations(
    tasks,
    workers: int,
    initializer_args: tuple,
) -> tuple[dict[str, sp.csr_matrix], dict[str, np.ndarray], float, int]:
    grams = {
        variant: sp.csr_matrix((LOADINGS[variant].shape[1],) * 2, dtype=float)
        for variant in VARIANTS
    }
    rhs = {
        variant: np.zeros(LOADINGS[variant].shape[1], dtype=float)
        for variant in VARIANTS
    }
    residual_sum_sq = 0.0
    n_pixels = 0
    errors = []
    for result in _run_pool(tasks, _accumulate_chunk, workers, initializer_args, "normal"):
        residual_sum_sq += float(result["residual_sum_sq"])
        n_pixels += int(result["n_pixels"])
        errors.extend(result["errors"])
        for variant in VARIANTS:
            grams[variant] += result["grams"][variant]
            rhs[variant] += result["rhs"][variant]
    if errors:
        raise RuntimeError(json.dumps(errors[:10], indent=2))
    return grams, rhs, residual_sum_sq, n_pixels


def _solve(
    gram: sp.csr_matrix,
    rhs: np.ndarray,
    ridge_lambda: float,
    loading: sp.csr_matrix,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    diagonal = np.asarray(gram.diagonal(), dtype=float)
    active = np.isfinite(diagonal) & (diagonal > 0.0)
    scale = np.zeros_like(diagonal)
    scale[active] = np.sqrt(diagonal[active])
    inverse_scale = sp.diags(1.0 / scale[active])
    normalized = inverse_scale @ gram[active][:, active] @ inverse_scale
    normalized = normalized + ridge_lambda * sp.eye(active.sum(), format="csr")
    p = sp.triu((normalized + normalized.T) * 0.5).tocsc()
    q = -rhs[active] / scale[active]
    # Clarabel uses A z + s = b with s >= 0.  With theta=D^-1 z,
    # A=-L D^-1 and b=1 impose 1 + L theta >= 0 in the global solve.
    constraint = -(loading[:, active] @ inverse_scale).tocsc()
    bound = np.ones(loading.shape[0], dtype=float)
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    solver = clarabel.DefaultSolver(
        p,
        np.asarray(q, dtype=float),
        constraint,
        bound,
        [clarabel.NonnegativeConeT(loading.shape[0])],
        settings,
    )
    result = solver.solve()
    status = str(result.status)
    if status not in ("Solved", "AlmostSolved"):
        raise RuntimeError(f"Clarabel failed: {status}")
    z = np.asarray(result.x, dtype=float)
    theta = np.zeros_like(rhs)
    theta[active] = z / scale[active]
    factor = 1.0 + np.asarray(loading @ theta).ravel()
    if not np.all(np.isfinite(factor)) or float(factor.min()) < -1.0e-7:
        raise RuntimeError(
            f"Clarabel returned an infeasible line factor: min={factor.min():.9g}"
        )
    info = {
        "solver": "clarabel",
        "solver_version": clarabel.__version__,
        "status": status,
        "iterations": int(getattr(result, "iterations", -1)),
        "solve_time_sec": float(getattr(result, "solve_time", np.nan)),
        "primal_residual": float(getattr(result, "r_prim", np.nan)),
        "dual_residual": float(getattr(result, "r_dual", np.nan)),
        "absolute_gap": float(
            abs(getattr(result, "obj_val", np.nan) - getattr(result, "obj_val_dual", np.nan))
        ),
        "relative_gap": float(
            abs(getattr(result, "obj_val", np.nan) - getattr(result, "obj_val_dual", np.nan))
            / max(abs(getattr(result, "obj_val_dual", np.nan)), 1.0)
        ),
        "minimum_line_factor": float(factor.min()),
        "active_factor_constraints": int(np.count_nonzero(factor <= 1.0e-7)),
    }
    return theta, scale, active, info


def _binned_nanmean(values: np.ndarray, width: int) -> np.ndarray:
    size = (values.size + width - 1) // width
    padded = np.full(size * width, np.nan, dtype=float)
    padded[: values.size] = values
    reshaped = padded.reshape(size, width)
    count = np.isfinite(reshaped).sum(axis=1)
    return np.divide(
        np.nansum(reshaped, axis=1),
        count,
        out=np.full(size, np.nan),
        where=count > 0,
    )


def _render_chunk(task: tuple[int, np.ndarray, dict[str, np.ndarray], int]):
    chunk_id, source_rows, theta, display_bin = task
    assert RZ_PIXEL is not None
    result = {
        "chunk_id": chunk_id,
        "source_rows": source_rows,
        "baseline": [],
        "oh_flux": [],
        "support_pixels": [],
        "baseline_rms": [],
        "post_C": [],
        "post_D": [],
        "post_C_rms": [],
        "post_D_rms": [],
        "errors": [],
    }
    for source_row in source_rows:
        try:
            residual, designs, use = _row_design(int(source_row))
            display = residual.copy()
            display[~np.isfinite(display)] = np.nan
            science = decompose._science_line_mask_for_row(int(source_row))
            if science is not None:
                display[science[RZ_PIXEL]] = np.nan
            result["baseline"].append(_binned_nanmean(display, display_bin))
            result["baseline_rms"].append(float(np.sqrt(np.mean(residual[use] ** 2))))
            result["support_pixels"].append(int(use.sum()))
            component = np.asarray(
                base._DECOMP_HDU["COMP_OH"].data[int(source_row)], dtype=float
            )
            result["oh_flux"].append(float(np.nansum(component[RZ_PIXEL])))
            for variant in VARIANTS:
                post = residual - np.asarray(designs[variant] @ theta[variant]).ravel()
                post_display = post.copy()
                post_display[~np.isfinite(display)] = np.nan
                result[f"post_{variant}"].append(
                    _binned_nanmean(post_display, display_bin)
                )
                result[f"post_{variant}_rms"].append(
                    float(np.sqrt(np.mean(post[use] ** 2)))
                )
        except Exception as error:
            result["errors"].append(
                {
                    "source_row": int(source_row),
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc(),
                }
            )
    for key in ("baseline", "post_C", "post_D"):
        result[key] = np.asarray(result[key], dtype=np.float32)
    for key in (
        "oh_flux",
        "baseline_rms",
        "post_C_rms",
        "post_D_rms",
    ):
        result[key] = np.asarray(result[key], dtype=float)
    result["support_pixels"] = np.asarray(result["support_pixels"], dtype=np.int32)
    return result


def render_diagnostics(
    tasks,
    theta: dict[str, np.ndarray],
    workers: int,
    initializer_args: tuple,
    n_rows: int,
    n_bins: int,
    display_bin: int,
) -> dict[str, np.ndarray]:
    arrays = {
        "baseline": np.empty((n_rows, n_bins), dtype=np.float32),
        "post_C": np.empty((n_rows, n_bins), dtype=np.float32),
        "post_D": np.empty((n_rows, n_bins), dtype=np.float32),
        "oh_flux": np.empty(n_rows),
        "support_pixels": np.empty(n_rows, dtype=np.int32),
        "baseline_rms": np.empty(n_rows),
        "post_C_rms": np.empty(n_rows),
        "post_D_rms": np.empty(n_rows),
    }
    errors = []
    render_tasks = [(i, rows, theta, display_bin) for i, rows in tasks]
    offsets = {chunk_id: sum(len(rows) for _, rows in tasks[:chunk_id]) for chunk_id, _ in tasks}
    for result in _run_pool(
        render_tasks, _render_chunk, workers, initializer_args, "render"
    ):
        errors.extend(result["errors"])
        start = offsets[int(result["chunk_id"])]
        stop = start + len(result["source_rows"])
        for key in arrays:
            arrays[key][start:stop] = result[key]
    if errors:
        raise RuntimeError(json.dumps(errors[:10], indent=2))
    return arrays


def _group_results(
    catalog: pd.DataFrame,
    groups: pd.DataFrame,
    loading: sp.csr_matrix,
    theta: np.ndarray,
    variant: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    delta = np.asarray(loading @ theta).ravel()
    frame = catalog.copy()
    frame["correction_delta"] = delta
    frame["correction_factor"] = 1.0 + delta
    rows = []
    for group_id, part in frame.groupby(f"group_{variant}", sort=True):
        q = part["q_aijc_gi"].to_numpy(dtype=float)
        value = part["correction_delta"].to_numpy(dtype=float)
        mean = float(np.average(value, weights=q))
        rows.append(
            {
                "group_id": int(group_id),
                "theta": float(theta[int(group_id)]),
                "effective_delta": mean,
                "correction_factor": 1.0 + mean,
                "within_group_delta_std": float(
                    np.sqrt(np.average((value - mean) ** 2, weights=q))
                ),
                "line_count": len(part),
                "wave_min_air_angstrom": float(part["wave"].min()),
                "wave_max_air_angstrom": float(part["wave"].max()),
                "q_sum": float(q.sum()),
            }
        )
    result = groups.merge(pd.DataFrame(rows), on="group_id", validate="one_to_one")
    return result, delta


def _closure_check(
    stack_path: Path,
    decomposition_path: Path,
    wave: np.ndarray,
    output_dir: Path,
    source_row: int,
) -> dict[str, float]:
    base._init_worker(
        str(stack_path),
        str(decomposition_path),
        wave,
        str(output_dir),
        "closure",
        None,
    )
    try:
        model = _model(source_row)
        oh_slice = model._group_slices["oh"]
        line = np.flatnonzero(
            (model._line_group >= oh_slice.start) & (model._line_group < oh_slice.stop)
        )
        coefficient_row = base._DECOMP_HDU["COEF"].data[source_row]
        names = tuple(f"OH_{index:03d}" for index in range(oh_slice.stop - oh_slice.start))
        coefficient = np.asarray([coefficient_row[name] for name in names], dtype=float)
        amplitude = model._base_line_weight[line] * coefficient[
            model._line_group[line] - oh_slice.start
        ]
        reconstructed = np.asarray(_raw_line_design(model)[:, line] @ amplitude).ravel()
        reference = np.asarray(base._DECOMP_HDU["COMP_OH"].data[source_row], dtype=float)
        difference = reconstructed - reference
        return {
            "source_row": source_row,
            "max_abs": float(np.nanmax(np.abs(difference))),
            "relative_l2": float(np.linalg.norm(difference) / np.linalg.norm(reference)),
            "sum_ratio": float(reconstructed.sum() / reference.sum()),
        }
    finally:
        base._close_worker_files()


def self_test() -> None:
    catalog = pd.DataFrame(
        {
            "parent_index": [0, 0, 0, 1, 1],
            "q_aijc_gi": [1.0, 2.0, 3.0, 4.0, 1.0],
            "family": ["a", "a", "b", "a", "b"],
        }
    )
    loading, _, _ = _group_loading(catalog, ["family"])
    for parent in (0, 1):
        use = catalog["parent_index"].to_numpy() == parent
        np.testing.assert_allclose(
            catalog.loc[use, "q_aijc_gi"].to_numpy() @ loading[use],
            0.0,
            atol=1.0e-14,
        )
    rng = np.random.default_rng(4)
    x1 = rng.normal(size=(7, 3))
    x2 = rng.normal(size=(5, 3))
    y1 = rng.normal(size=7)
    y2 = rng.normal(size=5)
    x = np.vstack([x1, x2])
    y = np.r_[y1, y2]
    np.testing.assert_allclose(x.T @ x, x1.T @ x1 + x2.T @ x2)
    np.testing.assert_allclose(x.T @ y, x1.T @ y1 + x2.T @ y2)
    values = np.array([1.0, np.nan, 3.0, 5.0, 7.0])
    qp_theta, _, _, qp_info = _solve(
        sp.eye(2, format="csr"),
        np.array([-2.0, 0.5]),
        1.0e-8,
        sp.eye(2, format="csr"),
    )
    np.testing.assert_allclose(qp_theta, [-1.0, 0.5], atol=1.0e-6)
    assert qp_info["minimum_line_factor"] >= -1.0e-7
    np.testing.assert_allclose(_binned_nanmean(values, 2), [1.0, 4.0, 7.0])
    print("self-test: ok")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stack", nargs="?", type=Path)
    parser.add_argument("decomposition", nargs="?", type=Path)
    parser.add_argument("--quality", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--ridge-lambda", type=float, default=RIDGE_LAMBDA)
    parser.add_argument("--display-bin", type=int, default=DISPLAY_BIN_PIXELS)
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.stack is None or args.decomposition is None:
        parser.error("stack and decomposition are required")
    if args.workers < 1 or args.chunk_size < 1 or args.display_bin < 1:
        raise ValueError("workers, chunk-size, and display-bin must be positive")
    if args.ridge_lambda < 0.0:
        raise ValueError("ridge-lambda must be non-negative")

    output_dir = args.output_dir or args.decomposition.parent / "far_sky_oh_branch_corrections_cd_v1"
    quality_path = args.quality or args.decomposition.parent / "far_sky_pca30" / "far_sky_quality_selection.csv"
    output_dir.mkdir(parents=True, exist_ok=True)
    selection = _load_quality(quality_path, args.max_rows)
    with fits.open(args.stack, memmap=True) as stack_hdu:
        wave = np.asarray(stack_hdu["WAVE"].data, dtype=float).copy()
    catalog, loadings, groups = build_catalog_and_loadings(wave)

    closure = _closure_check(
        args.stack, args.decomposition, wave, output_dir, int(selection.iloc[0]["source_row"])
    )
    if closure["relative_l2"] > 1.0e-10:
        raise ValueError(f"Exact OH closure check failed: {closure}")

    base._init_worker(
        str(args.stack), str(args.decomposition), wave, str(output_dir), "catalog", None
    )
    try:
        model = _model(int(selection.iloc[0]["source_row"]))
        oh_slice = model._group_slices["oh"]
        oh_line = np.flatnonzero(
            (model._line_group >= oh_slice.start) & (model._line_group < oh_slice.stop)
        )
        np.testing.assert_allclose(model._line_wave[oh_line], _production_oh_catalog(wave)["wave"])
        np.testing.assert_allclose(model._base_line_weight[oh_line], _production_oh_catalog(wave)["q_aijc_gi"])
        full_to_model = {int(full): int(model_index) for full, model_index in enumerate(oh_line)}
        selected_line = np.asarray([full_to_model[int(i)] for i in catalog["oh_index"]], dtype=int)
        selected_parent = model._line_group[selected_line] - oh_slice.start
        if not np.array_equal(selected_parent, catalog["parent_index"].to_numpy(dtype=int)):
            raise ValueError("Catalog parent ordering does not match the production model")
        selected_q = model._base_line_weight[selected_line].copy()
        oh_coef_names = tuple(
            f"OH_{index:03d}" for index in range(oh_slice.stop - oh_slice.start)
        )
    finally:
        base._close_worker_files()

    rz_pixel = wave >= R_LOWER_A
    rows = selection["source_row"].to_numpy(dtype=np.int64)
    tasks = [
        (index, rows[start : start + args.chunk_size])
        for index, start in enumerate(range(0, rows.size, args.chunk_size))
    ]
    payload = {
        "schema": SCHEMA,
        "source_code_sha256": file_sha256(Path(__file__)),
        "stack_path": str(args.stack.resolve()),
        "stack_sha256": file_sha256(args.stack),
        "decomposition_path": str(args.decomposition.resolve()),
        "decomposition_sha256": file_sha256(args.decomposition),
        "quality_path": str(quality_path.resolve()),
        "quality_sha256": file_sha256(quality_path),
        "wave_sha256": wave_sha256(wave),
        "ridge_lambda": args.ridge_lambda,
        "display_bin_pixels": args.display_bin,
        "retained_spectra": len(selection),
        "variant_C": "(v_upper,v_lower,branch_N+branch_J+F_upper+F_lower,N_upper bin)",
        "variant_D": "(v_upper,N_upper,F_upper,branch_N+branch_J+F_upper+F_lower)",
        "solver": "Clarabel sparse convex quadratic programming",
        "line_factor_constraint": "1 + delta_line >= 0 with delta_line = L theta",
        "contrast": "Aijc*gi weighted zero-sum within each production (v_upper,N_upper,F_upper) parent",
    }
    provenance = payload | {"run_fingerprint": _fingerprint(payload), "closure": closure}
    (output_dir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    global LOADINGS
    LOADINGS = loadings
    initializer_args = (
        str(args.stack),
        str(args.decomposition),
        wave,
        str(output_dir),
        provenance["run_fingerprint"],
        loadings,
        selected_line,
        selected_q,
        selected_parent,
        rz_pixel,
        oh_coef_names,
    )
    grams, rhs, residual_sum_sq, n_pixels = accumulate_normal_equations(
        tasks, args.workers, initializer_args
    )
    theta, scales, active, solver_info = {}, {}, {}, {}
    objective = {}
    for variant in VARIANTS:
        theta[variant], scales[variant], active[variant], solver_info[variant] = _solve(
            grams[variant],
            rhs[variant],
            args.ridge_lambda,
            loadings[variant],
        )
        post_sum_sq = (
            residual_sum_sq
            - 2.0 * float(theta[variant] @ rhs[variant])
            + float(theta[variant] @ (grams[variant] @ theta[variant]))
        )
        objective[variant] = {
            **solver_info[variant],
            "parameters": int(theta[variant].size),
            "active_parameters": int(active[variant].sum()),
            "baseline_rms": float(np.sqrt(residual_sum_sq / n_pixels)),
            "post_rms": float(np.sqrt(max(post_sum_sq, 0.0) / n_pixels)),
            "sum_squared_residual": float(post_sum_sq),
        }
        sp.save_npz(output_dir / f"normal_matrix_{variant}.npz", grams[variant])
        np.savez_compressed(
            output_dir / f"solution_{variant}.npz",
            theta=theta[variant],
            parameter_scale=scales[variant],
            active_parameter=active[variant],
            rhs=rhs[variant],
            solver_info_json=np.asarray(json.dumps(solver_info[variant], sort_keys=True)),
        )

    line_output = catalog.copy()
    parity_spread = {}
    factor_diagnostics = {}
    for variant in VARIANTS:
        group_result, delta = _group_results(
            catalog, groups[variant], loadings[variant], theta[variant], variant
        )
        group_result.to_csv(output_dir / f"groups_{variant}.csv", index=False)
        line_output[f"delta_{variant}"] = delta
        line_output[f"factor_{variant}"] = 1.0 + delta
        parity_spread[variant] = float(
            line_output.groupby(line_output["ID"].str[:-1])[f"delta_{variant}"].agg(
                lambda values: values.max() - values.min()
            ).max()
        )
        factor = 1.0 + delta
        factor_diagnostics[variant] = {
            "minimum": float(factor.min()),
            "at_zero_boundary": int(np.count_nonzero(factor <= 1.0e-7)),
            "below_feasibility_tolerance": int(np.count_nonzero(factor < -1.0e-7)),
            "negative_from_solver_tolerance": int(np.count_nonzero(factor < 0.0)),
        }
    line_output.to_csv(output_dir / "oh_line_corrections.csv", index=False)

    n_bins = int(np.ceil(rz_pixel.sum() / args.display_bin))
    diagnostics = render_diagnostics(
        tasks,
        theta,
        args.workers,
        initializer_args,
        len(selection),
        n_bins,
        args.display_bin,
    )
    binned_wave = _binned_nanmean(wave[rz_pixel], args.display_bin)
    np.savez_compressed(
        output_dir / "residual_maps.npz",
        source_row=rows,
        expnum=selection["expnum"].to_numpy(dtype=np.int64),
        sky_far_label=selection["sky_far_label"].astype(str).to_numpy(dtype=str),
        wave=binned_wave,
        **diagnostics,
    )
    metrics = selection[["source_row", "expnum", "sky_far_label"]].copy()
    for key in (
        "oh_flux",
        "support_pixels",
        "baseline_rms",
        "post_C_rms",
        "post_D_rms",
    ):
        metrics[key] = diagnostics[key]
    metrics.to_csv(output_dir / "spectrum_metrics.csv", index=False)

    summary = provenance | {
        "selected_oh_lines": len(catalog),
        "production_oh_parents": int(catalog["parent_index"].nunique()),
        "fit_pixels_total": n_pixels,
        "normal_equation_representation": "sum_s X_s.T X_s and sum_s X_s.T r_s",
        "group_counts": {variant: len(groups[variant]) for variant in VARIANTS},
        "parity_pair_max_delta_spread": parity_spread,
        "line_factor_diagnostics": factor_diagnostics,
        "objective": objective,
        "median_per_spectrum_rms": {
            "baseline": float(np.median(diagnostics["baseline_rms"])),
            "C": float(np.median(diagnostics["post_C_rms"])),
            "D": float(np.median(diagnostics["post_D_rms"])),
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
