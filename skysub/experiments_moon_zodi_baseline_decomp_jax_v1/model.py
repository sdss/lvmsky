"""Fit Moon+Zodi to baseline-decomposed full-resolution LVM spectra."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from hashlib import sha256
import json
import multiprocessing as mp
from pathlib import Path
import time

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

jax.config.update("jax_enable_x64", True)

from skysub.experiments_moon_zodi_corpus_jax_v1 import model as corpus_model  # noqa: E402
from skysub.experiments_moon_zodi_corpus_jax_v2 import model as previous  # noqa: E402
from skysub.sky_decomp.lsf_surface_iterative import (  # noqa: E402
    LSFSurfaceIterativeConfig,
    SkyDecompLSFSurfaceIterative,
)


EXPERIMENT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
DECOMPOSITION_CACHE = OUTPUT_DIR / "baseline_decomposition.npz"
DECOMPOSITION_METRICS = OUTPUT_DIR / "decomposition_metrics.csv"
FROZEN_PARAMETERS = OUTPUT_DIR / "frozen_parameters.json"
FAMILY = "A3"
FLUX_FACTOR = 1.0e14
FLUX_SCALE = previous.FLUX_SCALE
N_FOLDS = previous.N_FOLDS
N_STARTS = previous.N_STARTS
RANDOM_SEED = previous.RANDOM_SEED
DISPLAY_PERCENTILES = (1.0, 99.0)
DISPLAY_PADDING_FRACTION = 0.08
RESIDUAL_RMS_MULTIPLIER = 7.0
DISPLAY_BIN_WIDTH_ANGSTROM = 20.0
DISPLAY_BIN_MIN_PIXELS = 5
EXAMPLE_DPI = 260
CONTACT_SHEET_DPI = 200
BANDS = (
    ("B", 3600.0, 5800.0),
    ("R", 5800.0, 7550.0),
    ("Z", 7550.0, 9800.0),
)
LINE_COMPONENTS = ("oh", "atom", "orc", "o2")
HISTORICAL_LSF_CONFIG = LSFSurfaceIterativeConfig(n_refinement_cycles=4)


@dataclass(slots=True)
class PreparedData:
    """Physical bases plus baseline-decomposed observed targets and weights."""

    base: previous.TrainingData
    target: np.ndarray
    line_model: np.ndarray
    weights: np.ndarray
    skyline_mask: np.ndarray


@dataclass(frozen=True, slots=True)
class FitResult:
    """One deterministic bounded fit of the shared analytic parameters."""

    parameters: previous.FrozenParameters
    loss: float
    success: bool
    message: str
    iterations: int
    start_index: int


_WORKER_STACK = None
_WORKER_DECOMPOSER = None


def _json_safe(value: object) -> object:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _strict_corpus_rows() -> tuple[corpus_model.CorpusData, np.ndarray, np.ndarray]:
    corpus = corpus_model.load_corpus()
    corpus_rows = np.flatnonzero(previous.strict_lunar_mask(corpus))
    if corpus_rows.size != 356:
        raise RuntimeError(
            f"Expected the frozen 356-spectrum cohort, found {corpus_rows.size}"
        )
    mjd = previous._mjd_by_stack_row(corpus)[corpus_rows]
    return corpus, corpus_rows, mjd


def _init_decomposition_worker(source_stack: str, wave: np.ndarray) -> None:
    global _WORKER_STACK, _WORKER_DECOMPOSER
    _WORKER_STACK = fits.open(source_stack, memmap=True, lazy_load_hdus=True)
    _WORKER_DECOMPOSER = SkyDecompLSFSurfaceIterative(
        wave,
        lsf_sigma=0.5,
        base_dir=EXPERIMENT_DIR.parent,
        moon_smooth_lambda=0.1,
        config=HISTORICAL_LSF_CONFIG,
    )


def _line_subtracted_target(
    observed: np.ndarray,
    line_model: np.ndarray,
) -> np.ndarray:
    observed = np.asarray(observed, dtype=float)
    line_model = np.asarray(line_model, dtype=float)
    if observed.shape != line_model.shape:
        raise ValueError("Observed and line-model arrays must match")
    return observed - line_model


def _decompose_one(task: tuple[int, int, int]) -> dict[str, object]:
    output_row, stack_row, exposure = task
    if _WORKER_STACK is None or _WORKER_DECOMPOSER is None:
        raise RuntimeError("The decomposition worker is not initialized")
    observed = np.asarray(
        _WORKER_STACK["FLUX_SKY_FAR"].data[stack_row], dtype=np.float64
    ).copy()
    valid = np.isfinite(observed)
    observed_for_fit = np.where(valid, observed, 0.0) * FLUX_FACTOR
    result = _WORKER_DECOMPOSER.fit(
        observed_for_fit,
        valid.astype(float),
        verbose=False,
    )
    line_scaled = np.sum(
        np.vstack([result.components[name] for name in LINE_COMPONENTS]), axis=0
    )
    line_model = line_scaled / FLUX_FACTOR
    target = _line_subtracted_target(observed, line_model)
    weights = np.asarray(_WORKER_DECOMPOSER.continuum_weights, dtype=float)
    skyline_mask = np.asarray(_WORKER_DECOMPOSER.skyline_mask, dtype=bool)
    reconstructed = target + line_model
    closure = float(np.nanmax(np.abs(reconstructed - observed)))
    return {
        "output_row": output_row,
        "stack_row": stack_row,
        "exposure": exposure,
        "status": result.fit_status,
        "elapsed_seconds": float(result.fit_elapsed_sec),
        "r2": float(result.r2),
        "rms_residual_scaled": float(result.rms_resid),
        "closure_max_abs": closure,
        "line_model": line_model.astype(np.float32),
        "target": target.astype(np.float32),
        "weights": weights.astype(np.float32),
        "skyline_mask": skyline_mask,
        "skyline_fraction": float(np.mean(skyline_mask)),
        "weight_min": float(np.min(weights)),
        "weight_median": float(np.median(weights)),
        "weight_max": float(np.max(weights)),
    }


def prepare_decomposition(
    output: str | Path = DECOMPOSITION_CACHE,
    *,
    n_workers: int = 4,
) -> Path:
    """Run baseline LSF-surface decomposition for the frozen lunar cohort."""
    output = Path(output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    corpus, corpus_rows, mjd = _strict_corpus_rows()
    source_stack = Path(str(corpus.provenance["source_stack"])).resolve()
    with fits.open(source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    stack_rows = np.asarray(corpus.stack_row[corpus_rows], dtype=int)
    exposures = np.asarray(corpus.exposure[corpus_rows], dtype=int)
    tasks = [
        (index, int(stack_row), int(exposure))
        for index, (stack_row, exposure) in enumerate(
            zip(stack_rows, exposures, strict=True)
        )
    ]
    shape = (len(tasks), wave.size)
    targets = np.empty(shape, dtype=np.float32)
    line_models = np.empty(shape, dtype=np.float32)
    weights = np.empty(shape, dtype=np.float32)
    skyline_masks = np.empty(shape, dtype=bool)
    metrics: list[dict[str, object] | None] = [None] * len(tasks)
    parts_dir = output.parent / "decomposition_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    def part_path(index: int) -> Path:
        return parts_dir / f"part_{index:04d}.npz"

    pending_tasks = []
    for task in tasks:
        index, stack_row, exposure = task
        path = part_path(index)
        if not path.is_file():
            pending_tasks.append(task)
            continue
        with np.load(path, allow_pickle=False) as part:
            if int(part["stack_row"]) != stack_row or int(part["exposure"]) != exposure:
                raise RuntimeError(f"Checkpoint row mismatch: {path}")
            targets[index] = part["target"]
            line_models[index] = part["line_model"]
            weights[index] = part["weights"]
            skyline_masks[index] = part["skyline_mask"]
            metrics[index] = {
                "output_row": index,
                "stack_row": stack_row,
                "exposure": exposure,
                "status": str(part["status"]),
                "elapsed_seconds": float(part["elapsed_seconds"]),
                "r2": float(part["r2"]),
                "rms_residual_scaled": float(part["rms_residual_scaled"]),
                "closure_max_abs": float(part["closure_max_abs"]),
                "skyline_fraction": float(part["skyline_fraction"]),
                "weight_min": float(part["weight_min"]),
                "weight_median": float(part["weight_median"]),
                "weight_max": float(part["weight_max"]),
            }
    started = time.perf_counter()
    if pending_tasks:
        print(
            f"resuming with {len(tasks) - len(pending_tasks)} cached and "
            f"{len(pending_tasks)} pending spectra",
            flush=True,
        )
        context = mp.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=context,
            initializer=_init_decomposition_worker,
            initargs=(str(source_stack), wave),
        ) as executor:
            futures = [executor.submit(_decompose_one, task) for task in pending_tasks]
            for completed, future in enumerate(as_completed(futures), start=1):
                item = future.result()
                index = int(item["output_row"])
                target = np.asarray(item.pop("target"))
                line_model = np.asarray(item.pop("line_model"))
                continuum_weight = np.asarray(item.pop("weights"))
                skyline_mask = np.asarray(item.pop("skyline_mask"), dtype=bool)
                np.savez_compressed(
                    part_path(index),
                    target=target,
                    line_model=line_model,
                    weights=continuum_weight,
                    skyline_mask=skyline_mask,
                    **item,
                )
                targets[index] = target
                line_models[index] = line_model
                weights[index] = continuum_weight
                skyline_masks[index] = skyline_mask
                metrics[index] = item
                total_completed = len(tasks) - len(pending_tasks) + completed
                if completed % 10 == 0 or completed == len(pending_tasks):
                    elapsed = time.perf_counter() - started
                    print(
                        f"decomposed {total_completed}/{len(tasks)} spectra "
                        f"in {elapsed:.1f} s",
                        flush=True,
                    )
    if any(item is None for item in metrics):
        raise RuntimeError("The decomposition did not return every spectrum")
    metric_rows = [dict(item) for item in metrics if item is not None]
    failed = [
        row for row in metric_rows if row["status"] not in {"Solved", "AlmostSolved"}
    ]
    if failed:
        raise RuntimeError(f"Baseline decomposition failed for {len(failed)} spectra")
    np.savez_compressed(
        output,
        wave=wave,
        corpus_rows=corpus_rows,
        stack_rows=stack_rows,
        exposures=exposures,
        mjd=mjd,
        target=targets,
        line_model=line_models,
        weights=weights,
        skyline_mask=skyline_masks,
    )
    metrics_path = output.parent / DECOMPOSITION_METRICS.name
    _write_csv(metrics_path, metric_rows)
    provenance = {
        "experiment": EXPERIMENT_DIR.name,
        "source_stack": str(source_stack),
        "source_stack_size": source_stack.stat().st_size,
        "source_stack_mtime_ns": source_stack.stat().st_mtime_ns,
        "selection": corpus.provenance["selection"],
        "n_spectra": len(tasks),
        "n_nights": int(np.unique(mjd).size),
        "grid": {
            "n_pixels": wave.size,
            "minimum_angstrom": float(wave[0]),
            "maximum_angstrom": float(wave[-1]),
            "sampling_angstrom": float(np.median(np.diff(wave))),
            "dtype": str(wave.dtype),
        },
        "target": "FLUX_SKY_FAR minus baseline OH+atomic+Orc+O2 line model",
        "decomposition": {
            "implementation": "SkyDecompLSFSurfaceIterative",
            "flux_factor": FLUX_FACTOR,
            "line_weight": HISTORICAL_LSF_CONFIG.line_weight,
            "skyline_cumulative_fraction": HISTORICAL_LSF_CONFIG.skyline_cumulative_fraction,
            "skyline_half_width_angstrom": HISTORICAL_LSF_CONFIG.skyline_half_width_angstrom,
            "huber_transition_sigma": HISTORICAL_LSF_CONFIG.huber_transition_sigma,
            "n_refinement_cycles": HISTORICAL_LSF_CONFIG.n_refinement_cycles,
        },
        "cache": str(output),
        "cache_sha256": _sha256(output),
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output.parent / "decomposition_provenance.json").write_text(
        json.dumps(_json_safe(provenance), indent=2), encoding="utf-8"
    )
    return output


def load_prepared_data(
    cache: str | Path = DECOMPOSITION_CACHE,
) -> PreparedData:
    """Load frozen decomposition products and rebuild the physical JAX bases."""
    cache = Path(cache).expanduser().resolve()
    if not cache.is_file():
        raise FileNotFoundError(f"Run the decomposition first: {cache}")
    base = previous.build_training_data()
    with np.load(cache, allow_pickle=False) as payload:
        wave = np.asarray(payload["wave"], dtype=np.float64)
        stack_rows = np.asarray(payload["stack_rows"], dtype=int)
        target = np.asarray(payload["target"], dtype=np.float64)
        line_model = np.asarray(payload["line_model"], dtype=np.float64)
        weights = np.asarray(payload["weights"], dtype=np.float64)
        skyline_mask = np.asarray(payload["skyline_mask"], dtype=bool)
    with fits.open(base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        source_wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    if not np.array_equal(wave, source_wave):
        raise RuntimeError("The decomposition grid differs from the source grid")
    if not np.array_equal(stack_rows, base.stack_row):
        raise RuntimeError("The decomposition rows differ from the frozen cohort")
    expected = (len(base.exposure), wave.size)
    if (
        target.shape != expected
        or line_model.shape != expected
        or weights.shape != expected
    ):
        raise RuntimeError("The decomposition arrays have unexpected shapes")
    if skyline_mask.shape != expected:
        raise RuntimeError("The exposure-specific skyline masks have unexpected shape")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
        raise RuntimeError("Continuum weights must be finite and nonnegative")
    observed_native = np.asarray(base.observed_native, dtype=float)
    closure = np.nanmax(np.abs(target + line_model - observed_native))
    closure_tolerance = (
        16.0 * np.finfo(np.float32).eps * np.nanmax(np.abs(observed_native))
    )
    if not np.isfinite(closure) or closure > closure_tolerance:
        raise RuntimeError(
            f"Line-subtraction closure failed: {closure} > {closure_tolerance}"
        )
    return PreparedData(base, target, line_model, weights, skyline_mask)


def _geometry_inputs(data: PreparedData, rows: np.ndarray) -> previous.GeometryInputs:
    return previous._geometry_inputs(data.base, np.asarray(rows, dtype=int))


def _night_weights(mjd: np.ndarray) -> np.ndarray:
    return previous._night_weights(np.asarray(mjd, dtype=int))


def _weighted_loss(
    residual: jax.Array,
    weights: jax.Array,
    delta: float = previous.PSEUDO_HUBER_DELTA,
) -> jax.Array:
    residual = jnp.where(weights > 0.0, residual, 0.0)
    robust = delta**2 * (jnp.sqrt(1.0 + (residual / delta) ** 2) - 1.0)
    return jnp.sum(weights * robust, axis=1) / jnp.maximum(
        jnp.sum(weights, axis=1), 1.0e-30
    )


def _masked_observed(values: jax.Array, weights: jax.Array) -> jax.Array:
    """Remove invalid zero-weight values before they enter the JAX graph."""
    return jnp.where(weights > 0.0, values, 0.0)


def fit_parameters(
    data: PreparedData,
    rows: np.ndarray,
    *,
    n_starts: int = N_STARTS,
    max_iterations: int = 500,
) -> FitResult:
    """Fit one global A3 parameter vector using baseline continuum weights."""
    rows = np.asarray(rows, dtype=int)
    weights = jnp.asarray(data.weights[rows])
    observed = _masked_observed(jnp.asarray(data.target[rows] / FLUX_SCALE), weights)
    rayleigh = jnp.asarray(data.base.rayleigh_native[rows] / FLUX_SCALE)
    aerosol = jnp.asarray(data.base.aerosol_native[rows] / FLUX_SCALE)
    zodi = jnp.asarray(data.base.zodi_native[rows] / FLUX_SCALE)
    with fits.open(data.base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave_air = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    wave_ratio = jnp.asarray(
        previous.physical.air_to_vacuum(wave_air)[None, :] / 5000.0
    )
    geometry = _geometry_inputs(data, rows)
    night_weights = jnp.asarray(_night_weights(data.base.mjd[rows]))
    line_free = ~data.skyline_mask[rows] & np.isfinite(data.target[rows])
    level = np.nanmedian(
        np.where(line_free, data.target[rows] / FLUX_SCALE, np.nan), axis=1
    )
    if np.any(~np.isfinite(level) | (level <= 0.0)):
        raise RuntimeError("Every training spectrum needs a positive line-free level")
    epsilon = jnp.asarray(0.02 * level[:, None])

    def loss_function(values: jax.Array) -> jax.Array:
        full = previous._expand_family(values, FAMILY)
        prediction = previous.predict_from_basis_jax(
            full,
            rayleigh,
            aerosol,
            zodi,
            wave_ratio,
            geometry,
        )[2]
        denominator = jnp.maximum(prediction + observed + epsilon, epsilon)
        residual = 2.0 * (prediction - observed) / denominator
        per_exposure = _weighted_loss(residual, weights)
        return jnp.sum(night_weights * per_exposure)

    objective_jax = jax.jit(jax.value_and_grad(loss_function))

    def objective(values: np.ndarray) -> tuple[float, np.ndarray]:
        value, gradient = objective_jax(jnp.asarray(values, dtype=jnp.float64))
        value.block_until_ready()
        return float(value), np.asarray(gradient, dtype=np.float64)

    bounds = previous.FAMILY_BOUNDS[FAMILY]
    old = previous.load_frozen_parameters().values
    starts = [np.asarray(old, dtype=float), np.zeros(len(bounds), dtype=float)]
    rng = np.random.default_rng(RANDOM_SEED)
    while len(starts) < n_starts:
        starts.append(
            np.asarray(
                [rng.uniform(0.25 * lower, 0.25 * upper) for lower, upper in bounds],
                dtype=float,
            )
        )
    candidates = []
    for start_index, initial in enumerate(starts[:n_starts]):
        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": max_iterations, "ftol": 1.0e-12, "gtol": 1.0e-8},
        )
        candidates.append((float(result.fun), start_index, result))
    finite_candidates = [item for item in candidates if np.isfinite(item[0])]
    if not finite_candidates:
        raise RuntimeError("Every optimization start returned a non-finite loss")
    _, start_index, best = min(finite_candidates, key=lambda item: (item[0], item[1]))
    return FitResult(
        parameters=previous.FrozenParameters(
            FAMILY, tuple(float(value) for value in best.x)
        ),
        loss=float(best.fun),
        success=bool(best.success),
        message=str(best.message),
        iterations=int(best.nit),
        start_index=int(start_index),
    )


def predict_components(
    data: PreparedData,
    rows: np.ndarray,
    parameters: previous.FrozenParameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.asarray(rows, dtype=int)
    with fits.open(data.base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave_air = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    result = previous.predict_from_basis_jit(
        jnp.asarray(parameters.full_values()),
        jnp.asarray(data.base.rayleigh_native[rows]),
        jnp.asarray(data.base.aerosol_native[rows]),
        jnp.asarray(data.base.zodi_native[rows]),
        jnp.asarray(previous.physical.air_to_vacuum(wave_air)[None, :] / 5000.0),
        _geometry_inputs(data, rows),
    )
    result[2].block_until_ready()
    return tuple(np.asarray(item) for item in result)  # type: ignore[return-value]


def _fractional_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
    valid: np.ndarray,
) -> dict[str, np.ndarray | float]:
    level = np.nanmedian(np.where(valid, observed, np.nan), axis=1)
    denominator = np.maximum(observed, 0.02 * level[:, None])
    signed = np.where(valid, (predicted - observed) / denominator, np.nan)
    absolute = np.abs(signed)
    e50 = np.nanmedian(absolute, axis=1)
    e90 = np.nanpercentile(absolute, 90.0, axis=1)
    bias = np.nanmedian(signed, axis=1)
    return {
        "e50": e50,
        "e90": e90,
        "bias": bias,
        "median_e50": float(np.nanmedian(e50)),
        "p90_e50": float(np.nanpercentile(e50, 90.0)),
        "median_e90": float(np.nanmedian(e90)),
        "median_bias": float(np.nanmedian(bias)),
        "fraction_e50_le_5pct": float(np.mean(e50 <= 0.05)),
    }


def band_metrics(
    data: PreparedData,
    predicted: np.ndarray,
) -> dict[str, dict[str, np.ndarray | float | bool]]:
    with fits.open(data.base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    output = {}
    for band, lower, upper in BANDS:
        use = (wave >= lower) & (wave < upper)
        valid = (
            np.isfinite(data.target[:, use])
            & np.isfinite(predicted[:, use])
            & (data.target[:, use] > 0.0)
            & ~data.skyline_mask[:, use]
        )
        metrics = _fractional_metrics(data.target[:, use], predicted[:, use], valid)
        metrics["passed"] = bool(
            float(metrics["p90_e50"]) <= 0.05
            and abs(float(metrics["median_bias"])) <= 0.02
        )
        metrics["median_line_free_pixels_per_spectrum"] = int(
            np.median(np.sum(~data.skyline_mask[:, use], axis=1))
        )
        output[band] = metrics
    return output


def _select_examples(
    data: PreparedData, e50: np.ndarray, count: int = 20
) -> np.ndarray:
    features = np.column_stack(
        (
            data.base.target_airmass,
            data.base.moon_fli,
            data.base.moon_altitude_deg,
            data.base.moon_separation_deg,
            data.base.signed_phase_deg,
            e50,
        )
    )
    selected: list[int] = []
    for column in range(features.shape[1]):
        for index in (
            int(np.nanargmin(features[:, column])),
            int(np.nanargmax(features[:, column])),
        ):
            if index not in selected:
                selected.append(index)
    for quantile in np.linspace(0.0, 1.0, count):
        value = float(np.nanquantile(e50, quantile))
        index = int(np.nanargmin(np.abs(e50 - value)))
        if index not in selected:
            selected.append(index)
        if len(selected) == count:
            break
    for index in np.argsort(data.base.exposure):
        if int(index) not in selected:
            selected.append(int(index))
        if len(selected) == count:
            break
    return np.asarray(selected[:count], dtype=int)


def _continuum_display_limits(
    target: np.ndarray,
    skyline_mask: np.ndarray,
) -> tuple[float, float]:
    """Return a shared robust y-range determined only by line-free continuum."""
    target = np.asarray(target, dtype=float)
    skyline_mask = np.asarray(skyline_mask, dtype=bool)
    valid = np.isfinite(target) & ~skyline_mask
    if not np.any(valid):
        raise ValueError("A continuum display needs finite line-free pixels")
    lower, upper = np.nanpercentile(target[valid], DISPLAY_PERCENTILES)
    span = float(upper - lower)
    if not np.isfinite(span) or span <= 0.0:
        span = float(np.nanmax(np.abs(target[valid])))
    if not np.isfinite(span) or span <= 0.0:
        span = 1.0
    padding = DISPLAY_PADDING_FRACTION * span
    return min(0.0, float(lower - padding)), max(0.0, float(upper + padding))


def _line_free_residual_rms(
    residual: np.ndarray,
    skyline_mask: np.ndarray,
) -> float:
    """Compute ordinary RMS after excluding persisted skyline-covered pixels."""
    residual = np.asarray(residual, dtype=float)
    skyline_mask = np.asarray(skyline_mask, dtype=bool)
    valid = np.isfinite(residual) & ~skyline_mask
    if not np.any(valid):
        raise ValueError("A residual display needs finite line-free pixels")
    rms = float(np.sqrt(np.mean(residual[valid] ** 2)))
    if not np.isfinite(rms) or rms <= 0.0:
        raise ValueError("Line-free residual RMS must be positive and finite")
    return rms


def _display_nodes(
    wave: np.ndarray,
    values: np.ndarray,
    valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return display-only 20 Angstrom medians without changing model inputs."""
    wave = np.asarray(wave, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = np.asarray(valid, dtype=bool) & np.isfinite(values)
    edges = np.arange(
        BANDS[0][1],
        BANDS[-1][2] + DISPLAY_BIN_WIDTH_ANGSTROM,
        DISPLAY_BIN_WIDTH_ANGSTROM,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    nodes = np.full(centers.size, np.nan, dtype=float)
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        use = valid & (wave >= lower) & (wave < upper)
        if np.count_nonzero(use) >= DISPLAY_BIN_MIN_PIXELS:
            nodes[index] = float(np.median(values[use]))
    return centers, nodes


def save_example_png(
    data: PreparedData,
    index: int,
    moon: np.ndarray,
    zodi: np.ndarray,
    total: np.ndarray,
    output: str | Path,
) -> Path:
    """Save raw and line-subtracted observed/model/residual diagnostics."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with fits.open(data.base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    observed = np.asarray(data.base.observed_native[index], dtype=float)
    target = data.target[index]
    skyline_mask = data.skyline_mask[index]
    residual = target - total
    upper_limits = _continuum_display_limits(target, skyline_mask)
    residual_rms = _line_free_residual_rms(residual, skyline_mask)
    residual_limits = tuple(
        value * RESIDUAL_RMS_MULTIPLIER for value in (-residual_rms, residual_rms)
    )
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(15, 7),
        sharex="col",
        gridspec_kw={"height_ratios": (3.0, 1.0)},
        constrained_layout=True,
    )
    axes[0, 0].plot(wave, observed, color="black", lw=0.65, label="Observed", zorder=2)
    axes[0, 0].plot(
        wave,
        total,
        color="crimson",
        lw=1.15,
        label="Moon + Zodi continuum",
        zorder=10,
    )
    axes[0, 0].set_title("Raw spectrum with continuum prediction")
    axes[1, 0].plot(wave, residual, color="royalblue", lw=0.65)
    axes[0, 1].plot(
        wave,
        target,
        color="black",
        lw=0.65,
        label="Observed - skylines",
        zorder=2,
    )
    axes[0, 1].plot(
        wave, moon, color="darkorange", lw=0.75, ls="--", label="Moon", zorder=4
    )
    axes[0, 1].plot(
        wave, zodi, color="purple", lw=0.75, ls="--", label="Zodi", zorder=4
    )
    axes[0, 1].plot(
        wave, total, color="crimson", lw=1.15, label="Moon + Zodi", zorder=10
    )
    axes[0, 1].set_title("Continuum training target")
    axes[1, 1].plot(wave, residual, color="royalblue", lw=0.65)
    for axis in axes[1]:
        axis.axhline(0.0, color="black", lw=0.7, alpha=0.7)
        axis.set_xlabel("Air wavelength (Angstrom)")
        axis.set_ylabel("Residual")
    for axis in axes[0]:
        axis.set_ylabel("Flux per Far-Sky fibre")
        axis.set_ylim(upper_limits)
        axis.legend(fontsize=7, ncol=2, loc="upper right")
        axis.grid(alpha=0.15)
    for axis in axes[1]:
        axis.set_ylim(residual_limits)
        axis.grid(alpha=0.15)
        axis.text(
            0.01,
            0.94,
            f"range = +/-7 line-free RMS; RMS={residual_rms:.2e}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )
    figure.suptitle(
        f"exp={int(data.base.exposure[index])}  MJD={int(data.base.mjd[index])}  "
        f"FLI={data.base.moon_fli[index]:.3f}  "
        f"airmass={data.base.target_airmass[index]:.2f}  "
        f"Moon alt={data.base.moon_altitude_deg[index]:.1f} deg  "
        f"separation={data.base.moon_separation_deg[index]:.1f} deg"
    )
    figure.savefig(output, dpi=EXAMPLE_DPI)
    plt.close(figure)
    return output


def save_20a_example_png(
    data: PreparedData,
    index: int,
    moon: np.ndarray,
    zodi: np.ndarray,
    total: np.ndarray,
    output: str | Path,
) -> Path:
    """Save a display-only 20 Angstrom view for one held-out example."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with fits.open(data.base.source_stack, memmap=True, lazy_load_hdus=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
    observed = np.asarray(data.base.observed_native[index], dtype=float)
    target = np.asarray(data.target[index], dtype=float)
    skyline_mask = np.asarray(data.skyline_mask[index], dtype=bool)
    all_valid = np.isfinite(wave)
    clean_valid = all_valid & ~skyline_mask
    centers, observed_nodes = _display_nodes(wave, observed, all_valid)
    _, target_nodes = _display_nodes(wave, target, clean_valid)
    _, moon_nodes = _display_nodes(wave, moon, clean_valid)
    _, zodi_nodes = _display_nodes(wave, zodi, clean_valid)
    _, total_nodes = _display_nodes(wave, total, clean_valid)
    continuum_residual = target_nodes - total_nodes
    node_valid = np.isfinite(target_nodes)
    upper_limits = _continuum_display_limits(target_nodes, ~node_valid)
    residual_rms = _line_free_residual_rms(continuum_residual, ~node_valid)
    residual_limits = tuple(
        value * RESIDUAL_RMS_MULTIPLIER for value in (-residual_rms, residual_rms)
    )
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(15, 7),
        sharex="col",
        gridspec_kw={"height_ratios": (3.0, 1.0)},
        constrained_layout=True,
    )
    axes[0, 0].plot(
        centers,
        observed_nodes,
        "o-",
        color="black",
        ms=2.0,
        lw=0.7,
        label="Observed",
        zorder=2,
    )
    axes[0, 0].plot(
        centers,
        total_nodes,
        "o-",
        color="crimson",
        ms=2.0,
        lw=1.1,
        label="Moon + Zodi continuum",
        zorder=10,
    )
    axes[0, 0].set_title("Raw spectrum: continuum model at 20 Angstrom nodes")
    axes[1, 0].plot(
        centers, continuum_residual, "o-", color="royalblue", ms=1.8, lw=0.7
    )
    axes[0, 1].plot(
        centers,
        target_nodes,
        "o-",
        color="black",
        ms=2.0,
        lw=0.7,
        label="Observed - skylines",
        zorder=2,
    )
    axes[0, 1].plot(
        centers,
        moon_nodes,
        "o--",
        color="darkorange",
        ms=1.8,
        lw=0.8,
        label="Moon",
        zorder=4,
    )
    axes[0, 1].plot(
        centers,
        zodi_nodes,
        "o--",
        color="purple",
        ms=1.8,
        lw=0.8,
        label="Zodi",
        zorder=4,
    )
    axes[0, 1].plot(
        centers,
        total_nodes,
        "o-",
        color="crimson",
        ms=2.0,
        lw=1.1,
        label="Moon + Zodi",
        zorder=10,
    )
    axes[0, 1].set_title("Continuum target: line-free 20 Angstrom medians")
    axes[1, 1].plot(
        centers, continuum_residual, "o-", color="royalblue", ms=1.8, lw=0.7
    )
    for axis in axes[0]:
        axis.set_ylim(upper_limits)
        axis.set_ylabel("Flux per Far-Sky fibre")
        axis.legend(fontsize=7, ncol=2, loc="upper right")
        axis.grid(alpha=0.15)
    for axis in axes[1]:
        axis.axhline(0.0, color="black", lw=0.7, alpha=0.7)
        axis.set_ylim(residual_limits)
        axis.set_xlabel("Air wavelength (Angstrom)")
        axis.set_ylabel("Residual")
        axis.grid(alpha=0.15)
    figure.suptitle(
        f"exp={int(data.base.exposure[index])}  MJD={int(data.base.mjd[index])}  "
        f"FLI={data.base.moon_fli[index]:.3f}  "
        f"airmass={data.base.target_airmass[index]:.2f}  "
        f"Moon alt={data.base.moon_altitude_deg[index]:.1f} deg  "
        f"separation={data.base.moon_separation_deg[index]:.1f} deg"
    )
    figure.savefig(output, dpi=EXAMPLE_DPI)
    plt.close(figure)
    return output


def _save_contact_sheet(paths: list[Path], output: Path) -> None:
    images = [plt.imread(path) for path in paths]
    figure, axes = plt.subplots(5, 4, figsize=(20, 17), constrained_layout=True)
    for axis, image, path in zip(axes.ravel(), images, paths, strict=True):
        axis.imshow(image)
        axis.set_title(path.stem, fontsize=7)
        axis.axis("off")
    figure.savefig(output, dpi=CONTACT_SHEET_DPI)
    plt.close(figure)


def save_ingredients_png(
    result: previous.ModelResult,
    output: str | Path,
) -> Path:
    """Save the complete full-LVM analytic ingredient curves."""
    output = Path(output)
    prepared = corpus_model.MoonZodiModel().prepare(result.spectrum)
    base_prepared = prepared.base_prepared
    wave = base_prepared.wave_hr_vacuum_angstrom
    tau_aerosol = previous.AOD500 * (wave / 5000.0) ** (-previous.AEROSOL_ALPHA)
    correction = np.asarray(
        previous._geometry_correction(
            jnp.asarray(result.parameters.full_values()),
            jnp.asarray(wave / 5000.0),
            previous.GeometryInputs(
                target_airmass=jnp.asarray(result.geometry.target_airmass),
                moon_airmass=jnp.asarray(result.geometry.moon_airmass),
                signed_phase_deg=jnp.asarray(result.geometry.signed_phase_deg),
                cos_separation=jnp.asarray(
                    np.cos(np.deg2rad(result.geometry.moon_separation_deg))
                ),
            ),
        )
    )

    def normalized(values: np.ndarray) -> np.ndarray:
        reference = (wave >= 4950.0) & (wave <= 5050.0)
        return values / np.nanmedian(values[reference])

    ingredients = (
        ("Meftah solar carrier", normalized(base_prepared.moon_solar_hr), "black"),
        ("ROLO lunar albedo", normalized(base_prepared.lunar_albedo_hr), "darkorange"),
        (
            "Rayleigh optical depth",
            normalized(base_prepared.tau_rayleigh_hr),
            "steelblue",
        ),
        ("Aerosol optical depth", normalized(tau_aerosol), "seagreen"),
        ("Fitted analytic correction", normalized(correction), "crimson"),
    )
    figure, axis = plt.subplots(figsize=(13, 5), constrained_layout=True)
    for label, values, color in ingredients:
        axis.plot(wave, values, lw=1.0, color=color, label=label)
    axis.set(
        xlim=(3600.0, 9800.0),
        xlabel="Vacuum wavelength (Angstrom)",
        ylabel="Factor normalized near 5000 Angstrom",
        title=f"Full-LVM analytic ingredients for exposure {result.spectrum.exposure}",
    )
    axis.grid(alpha=0.2)
    axis.legend(fontsize=8, ncol=2)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return output


def render_diagnostics(
    output_dir: str | Path = OUTPUT_DIR,
    *,
    cache: str | Path = DECOMPOSITION_CACHE,
) -> list[dict[str, object]]:
    """Regenerate full-resolution and 20 Angstrom OOF galleries without refitting."""
    output_dir = Path(output_dir).expanduser().resolve()
    data = load_prepared_data(cache)
    manifest_path = output_dir / "example_manifest.csv"
    predictions_path = output_dir / "oof_predictions.npz"
    if not manifest_path.is_file() or not predictions_path.is_file():
        raise FileNotFoundError("Run the fit before rendering held-out diagnostics")
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle))
    if len(manifest) != 20:
        raise RuntimeError(f"Expected 20 manifest examples, found {len(manifest)}")
    with np.load(predictions_path, allow_pickle=False) as predictions:
        oof_moon = np.asarray(predictions["moon"], dtype=float)
        oof_zodi = np.asarray(predictions["zodi"], dtype=float)
        oof_total = np.asarray(predictions["total"], dtype=float)
    example_dir = output_dir / "examples"
    example_20a_dir = output_dir / "examples_20a"
    full_paths = []
    node_paths = []
    updated_rows = []
    for record in manifest:
        number = int(record["display_order"])
        index = int(record["cohort_index"])
        exposure = int(data.base.exposure[index])
        if exposure != int(record["exposure"]):
            raise RuntimeError("Example manifest no longer matches the frozen cohort")
        full_path = example_dir / f"{number:02d}_exp{exposure}_full_lvm.png"
        node_path = example_20a_dir / f"{number:02d}_exp{exposure}_20a.png"
        full_paths.append(
            save_example_png(
                data,
                index,
                oof_moon[index],
                oof_zodi[index],
                oof_total[index],
                full_path,
            )
        )
        node_paths.append(
            save_20a_example_png(
                data,
                index,
                oof_moon[index],
                oof_zodi[index],
                oof_total[index],
                node_path,
            )
        )
        updated_rows.append(
            {
                "display_order": number,
                "cohort_index": index,
                "exposure": exposure,
                "mjd": int(data.base.mjd[index]),
                "filename": full_path.name,
                "filename_20a": node_path.name,
            }
        )
    _write_csv(manifest_path, updated_rows)
    _save_contact_sheet(full_paths, output_dir / "examples_contact_sheet.png")
    _save_contact_sheet(node_paths, output_dir / "examples_20a_contact_sheet.png")
    return updated_rows


def run_fit(
    output_dir: str | Path = OUTPUT_DIR,
    *,
    cache: str | Path = DECOMPOSITION_CACHE,
) -> dict[str, object]:
    """Run MJD-held-out fitting, freeze the final fit, and write diagnostics."""
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    data = load_prepared_data(cache)
    n_spectra = len(data.base.exposure)
    oof_moon = np.full_like(data.target, np.nan, dtype=np.float32)
    oof_zodi = np.full_like(data.target, np.nan, dtype=np.float32)
    oof_total = np.full_like(data.target, np.nan, dtype=np.float32)
    fold_payload = []
    for fold in range(N_FOLDS):
        train_rows = np.flatnonzero(data.base.fold != fold)
        valid_rows = np.flatnonzero(data.base.fold == fold)
        fit = fit_parameters(data, train_rows)
        if not fit.success or not np.isfinite(fit.loss):
            raise RuntimeError(
                f"Fold {fold} optimization failed: {fit.message}; loss={fit.loss}"
            )
        moon, zodi, total = predict_components(data, valid_rows, fit.parameters)
        oof_moon[valid_rows] = moon
        oof_zodi[valid_rows] = zodi
        oof_total[valid_rows] = total
        fold_payload.append(
            {
                "fold": fold,
                "n_train": train_rows.size,
                "n_validation": valid_rows.size,
                "loss": fit.loss,
                "success": fit.success,
                "message": fit.message,
                "iterations": fit.iterations,
                "start_index": fit.start_index,
                "parameters": dict(
                    zip(fit.parameters.names, fit.parameters.values, strict=True)
                ),
            }
        )
        print(f"completed fold {fold + 1}/{N_FOLDS}", flush=True)
    metrics = band_metrics(data, oof_total)
    final_fit = fit_parameters(data, np.arange(n_spectra))
    if not final_fit.success or not np.isfinite(final_fit.loss):
        raise RuntimeError(
            f"Final optimization failed: {final_fit.message}; loss={final_fit.loss}"
        )
    parameter_rows = [
        {
            "name": name,
            "value": value,
            "lower_bound": bounds[0],
            "upper_bound": bounds[1],
            "role": "shared fitted parameter",
        }
        for name, value, bounds in zip(
            final_fit.parameters.names,
            final_fit.parameters.values,
            previous.FAMILY_BOUNDS[FAMILY],
            strict=True,
        )
    ]
    _write_csv(output_dir / "parameter_table.csv", parameter_rows)
    frozen = {
        "experiment": EXPERIMENT_DIR.name,
        "family": FAMILY,
        "parameter_names": list(final_fit.parameters.names),
        "values": list(final_fit.parameters.values),
        "bounds": [list(item) for item in previous.FAMILY_BOUNDS[FAMILY]],
        "fit": {
            "loss": final_fit.loss,
            "success": final_fit.success,
            "message": final_fit.message,
            "iterations": final_fit.iterations,
            "start_index": final_fit.start_index,
        },
        "fixed_physical_parameters": {
            "aod500": previous.AOD500,
            "aerosol_alpha": previous.AEROSOL_ALPHA,
            "aerosol_single_scattering_albedo": previous.AEROSOL_SINGLE_SCATTERING_ALBEDO,
            "aerosol_hg_g": previous.AEROSOL_HG_G,
            "zodi_scale": previous.ZODI_SCALE,
            "zodi_gamma": previous.ZODI_GAMMA,
            "kappa_blue": previous.KAPPA_BLUE,
        },
    }
    (output_dir / FROZEN_PARAMETERS.name).write_text(
        json.dumps(_json_safe(frozen), indent=2), encoding="utf-8"
    )
    (output_dir / "fold_parameters.json").write_text(
        json.dumps(_json_safe(fold_payload), indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        output_dir / "oof_predictions.npz",
        moon=oof_moon,
        zodi=oof_zodi,
        total=oof_total,
    )
    exposure_rows = []
    per_band = {}
    for band, lower, upper in BANDS:
        with fits.open(
            data.base.source_stack, memmap=True, lazy_load_hdus=True
        ) as hdul:
            wave = np.asarray(hdul["WAVE"].data, dtype=np.float64).copy()
        use = (wave >= lower) & (wave < upper)
        valid = (
            np.isfinite(data.target[:, use])
            & np.isfinite(oof_total[:, use])
            & (data.target[:, use] > 0.0)
            & ~data.skyline_mask[:, use]
        )
        per_band[band] = _fractional_metrics(
            data.target[:, use], oof_total[:, use], valid
        )
    for index in range(n_spectra):
        row: dict[str, object] = {
            "exposure": int(data.base.exposure[index]),
            "stack_row": int(data.base.stack_row[index]),
            "mjd": int(data.base.mjd[index]),
            "fold": int(data.base.fold[index]),
            "sky_far_label": str(data.base.sky_far_label[index]),
            "target_airmass": float(data.base.target_airmass[index]),
            "moon_fli": float(data.base.moon_fli[index]),
            "moon_altitude_deg": float(data.base.moon_altitude_deg[index]),
            "moon_separation_deg": float(data.base.moon_separation_deg[index]),
            "signed_phase_deg": float(data.base.signed_phase_deg[index]),
        }
        for band in ("B", "R", "Z"):
            row[f"e50_{band.lower()}"] = float(np.asarray(per_band[band]["e50"])[index])
            row[f"e90_{band.lower()}"] = float(np.asarray(per_band[band]["e90"])[index])
            row[f"bias_{band.lower()}"] = float(
                np.asarray(per_band[band]["bias"])[index]
            )
        exposure_rows.append(row)
    _write_csv(output_dir / "per_exposure_oof_metrics.csv", exposure_rows)
    e50_mean = np.nanmean(
        np.column_stack([per_band[band]["e50"] for band in ("B", "R", "Z")]),
        axis=1,
    )
    examples = _select_examples(data, e50_mean, 20)
    example_dir = output_dir / "examples"
    example_20a_dir = output_dir / "examples_20a"
    example_paths = []
    example_20a_paths = []
    example_rows = []
    for number, index in enumerate(examples, start=1):
        path = (
            example_dir
            / f"{number:02d}_exp{int(data.base.exposure[index])}_full_lvm.png"
        )
        example_paths.append(
            save_example_png(
                data,
                int(index),
                oof_moon[index],
                oof_zodi[index],
                oof_total[index],
                path,
            )
        )
        path_20a = (
            example_20a_dir
            / f"{number:02d}_exp{int(data.base.exposure[index])}_20a.png"
        )
        example_20a_paths.append(
            save_20a_example_png(
                data,
                int(index),
                oof_moon[index],
                oof_zodi[index],
                oof_total[index],
                path_20a,
            )
        )
        example_rows.append(
            {
                "display_order": number,
                "cohort_index": int(index),
                "exposure": int(data.base.exposure[index]),
                "mjd": int(data.base.mjd[index]),
                "filename": path.name,
                "filename_20a": path_20a.name,
            }
        )
    _write_csv(output_dir / "example_manifest.csv", example_rows)
    _save_contact_sheet(example_paths, output_dir / "examples_contact_sheet.png")
    _save_contact_sheet(
        example_20a_paths,
        output_dir / "examples_20a_contact_sheet.png",
    )
    representative = int(np.nanargmin(np.abs(e50_mean - np.nanmedian(e50_mean))))
    runtime_model = previous.AnalyticMoonZodiModel(final_fit.parameters)
    spectrum = corpus_model.load_spectrum(
        data.base.source_stack,
        row=int(data.base.stack_row[representative]),
        load_flux=True,
    )
    save_ingredients_png(
        runtime_model.evaluate(spectrum), output_dir / "analytic_ingredients.png"
    )
    summary = {
        "experiment": EXPERIMENT_DIR.name,
        "status": "accepted"
        if all(bool(metrics[band]["passed"]) for band in metrics)
        else "failed_5pct_gate",
        "approved_contract": {
            "cohort": "unchanged strict 356-spectrum lunar cohort",
            "target": "observed FLUX_SKY_FAR minus baseline OH+atomic+Orc+O2",
            "grid": "unchanged 12401-pixel source wavelength grid",
            "fit_range_angstrom": [3600.0, 9800.0],
            "weights": "baseline LSF-surface iterative continuum weights",
            "line_weight": HISTORICAL_LSF_CONFIG.line_weight,
            "per_exposure_fitted_normalization": False,
            "split": "five folds grouped by MJD",
            "model": "unchanged analytic A3 Moon+fixed-Zodi family",
        },
        "n_spectra": n_spectra,
        "n_nights": int(np.unique(data.base.mjd).size),
        "n_native_pixels": data.target.shape[1],
        "skyline_fraction": float(np.mean(data.skyline_mask)),
        "oof_band_metrics": {
            band: {
                key: value for key, value in metrics[band].items() if np.isscalar(value)
            }
            for band in metrics
        },
        "folds": fold_payload,
        "final_fit": frozen["fit"],
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output_dir / "fit_results.json").write_text(
        json.dumps(_json_safe(summary), indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "fit", "render", "all"))
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--n-workers", type=int, default=4)
    args = parser.parse_args()
    cache = args.output_dir / DECOMPOSITION_CACHE.name
    if args.command in {"prepare", "all"}:
        prepare_decomposition(cache, n_workers=args.n_workers)
    if args.command in {"fit", "all"}:
        summary = run_fit(args.output_dir, cache=cache)
        print(json.dumps(_json_safe(summary), indent=2))
    if args.command == "render":
        rows = render_diagnostics(args.output_dir, cache=cache)
        print(f"rendered {len(rows)} full-resolution and 20 Angstrom examples")


if __name__ == "__main__":
    main()
