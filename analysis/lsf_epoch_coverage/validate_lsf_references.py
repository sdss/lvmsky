"""Validate published means, weight accounting, covariance factors and HTTP assets."""
from pathlib import Path
import json
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.lsf_epoch_coverage.build_lsf_references import DEST
from skysub.sky_decomp.lsf_reference import validate_coefficients, load_reference_lsf, reference_covariance
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_density, evaluate_lsf_diagnostics


def validate(destination=DEST):
    destination = Path(destination)
    cache = np.load(destination.parent / "lsf_reference_audit/cache/input_coefficients.npy", mmap_mode="r")
    catalog = pd.read_parquet(destination.parent / "lsf_reference_audit/publication_records/reference_catalog.parquet")
    paths = [destination / name for name in catalog.file.unique()]
    largest_mean_error = 0.0
    max_weight_error = 0.0
    inspected = 0
    for path in paths:
        with fits.open(path, memmap=True) as compact:
            archive = path.parent / compact[0].header["AUDIT"]
            with fits.open(archive) as full:
                for name in ("LSF_COEF", "LSF_KNOTS", "WAVE", "COEF_SCATTER", "COEF_CONDERR"):
                    np.testing.assert_array_equal(compact[name].data, full[name].data)
        with fits.open(archive, memmap=True) as hdus:
            refs = hdus["REFERENCES"].data
            coef = np.asarray(hdus["LSF_COEF"].data, dtype=float)
            validate_coefficients(coef, tolerance=2e-12)
            contributors = hdus["CONTRIBUTORS"].data
            refids = np.asarray(contributors["reference_id"], dtype=int)
            weights = np.asarray(contributors["weight"], dtype=float)
            if not np.isfinite(weights).all() or np.any(weights <= 0):
                raise ValueError("Contributing weights must be finite and positive")
            sums = np.bincount(refids, weights=weights, minlength=len(refs))
            np.testing.assert_allclose(sums, 1, rtol=0, atol=2e-12)
            max_weight_error = max(max_weight_error, float(np.max(np.abs(sums-1))))
            np.testing.assert_array_equal(np.bincount(refids, minlength=len(refs)), refs["n_exposures"])
            np.testing.assert_allclose(1/np.bincount(refids, weights=weights**2), refs["n_eff"], rtol=2e-12)
            expected = np.zeros_like(coef)
            variance = np.zeros((len(refs), 44))
            for start in range(0, len(contributors), 20000):
                end = start + 20000
                x = np.asarray(cache[contributors["sample_id"][start:end]], dtype=float)
                x /= x.sum(axis=1, keepdims=True)
                ids = refids[start:end]
                alpha = weights[start:end]
                np.add.at(expected, ids, alpha[:, None, None]*x)
                factors = np.asarray(hdus["SCATTER_FACTORS"].data[start:end], dtype=float)
                np.testing.assert_allclose(factors, np.sqrt(alpha)[:, None]*(x-coef[ids]).reshape(len(x), 44), rtol=2e-6, atol=2e-8)
                np.add.at(variance, ids, factors**2)
            np.testing.assert_allclose(expected, coef, rtol=0, atol=2e-12)
            largest_mean_error = max(largest_mean_error, float(np.max(np.abs(expected-coef))))
            denominator = refs["scatter_denominator"].astype(float)
            good = denominator > 0
            actual_scatter = np.asarray(hdus["COEF_SCATTER"].data).reshape(len(refs), 44)
            np.testing.assert_allclose(variance[good]/denominator[good, None], actual_scatter[good]**2, rtol=2e-6, atol=1e-15)
            if not np.isnan(actual_scatter[~good]).all():
                raise ValueError("Singleton scatter must be undefined")
            nightids = np.asarray(hdus["NIGHTS"].data["reference_id"], dtype=int)
            nightfactors = np.asarray(hdus["NIGHT_FACTORS"].data, dtype=float)
            np.testing.assert_array_equal(np.bincount(nightids, minlength=len(refs)), refs["n_nights"])
            nightvariance = np.zeros((len(refs), 44))
            np.add.at(nightvariance, nightids, nightfactors**2)
            count = refs["n_nights"].astype(float)
            good = count > 1
            errors = np.asarray(hdus["COEF_CONDERR"].data).reshape(len(refs), 44)
            np.testing.assert_allclose(nightvariance[good]*count[good, None]/(count[good, None]-1), errors[good]**2, rtol=2e-6, atol=1e-15)
            if not np.isnan(errors[~good]).all():
                raise ValueError("Single-night conditional errors must be undefined")
            candidates = pd.DataFrame(dict(fiberid=refs["fiberid"].astype(np.int64), channel=np.asarray(refs["channel"], dtype=str)))
            coverage = candidates.groupby("fiberid").channel.nunique()
            complete = coverage[coverage == 3]
            if len(complete):
                fiber = int(complete.index[0])
                state = load_reference_lsf(path, fiber, allow_epoch_pooling=True)
                wavelengths = np.array([4700., 6600., 8500.])
                diag = evaluate_lsf_diagnostics(state, wavelengths)
                density = evaluate_lsf_density(state, wavelengths, np.linspace(-2.5, 2.5, 501))
                if np.any(density < 0) or not np.isfinite(density).all():
                    raise ValueError("Loaded density must be finite and nonnegative")
                for ci, channel in enumerate(("B", "R", "Z")):
                    i = state.metrics[channel]["reference_id"]
                    np.testing.assert_array_equal(state.coefficients[channel], coef[i])
                    for metric in ("fwhm_angstrom", "w50_angstrom", "w90_angstrom", "centroid_angstrom"):
                        np.testing.assert_allclose(diag[metric][ci], refs[i][metric+"_reference"], rtol=1e-9, atol=1e-9)
                    cov = reference_covariance(path, i)
                    if refs[i]["n_exposures"] > 1:
                        direction = np.zeros(44)
                        direction[::4] = 1
                        np.testing.assert_allclose(cov @ direction, 0, atol=2e-8)
                inspected += 1
        print(f"Validated {path.relative_to(destination)}", flush=True)
    assets = [p for p in destination.rglob("*") if p.is_file() and p.suffix in (".fits", ".md", ".csv", ".html", ".png", ".json", ".parquet", ".py", ".js")]
    base = "http://lumen.abudhabi.nyu.edu/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132/lsf_calibration/"
    def check(path):
        relative = str(path.relative_to(destination))
        with urllib.request.urlopen(urllib.request.Request(base+relative, method="HEAD"), timeout=30) as response:
            if response.status != 200 or int(response.headers["Content-Length"]) != path.stat().st_size:
                raise ValueError(f"HTTP asset differs: {relative}")
        return relative
    with ThreadPoolExecutor(max_workers=8) as pool:
        checked = list(pool.map(check, assets))
    with urllib.request.urlopen(base, timeout=30) as response:
        listing = response.read().decode()
        if "lsf_reference_surfaces.html" not in listing or "README.md" not in listing:
            raise ValueError("Directory listing must remain available")
    result = dict(fits_files=len(paths), reference_surfaces=len(catalog),
                  compact_arrays_identical_to_archive=True,
                  checked_mean_and_weight_rows=len(catalog), max_mean_abs_error=largest_mean_error,
                  max_weight_sum_abs_error=max_weight_error, files_with_loaded_B_R_Z_state=inspected,
                  http_assets=len(checked), directory_listing_verified=True,
                  covariance_diagonals_checked="all references", native_evaluator_agreement="one complete fiber per file",
                  independent_prediction_validation="not performed")
    (destination.parent / "lsf_reference_audit/publication_records/validation.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == "__main__":
    validate()
