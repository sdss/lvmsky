"""Publish constrained, dispersion-weighted per-fiber LSF reference FITS products."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table, vstack
import matplotlib
matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.lsf_epoch_coverage.build_coverage import OUT, DATA, fit_files, frame
from analysis.lsf_epoch_coverage.analyze_stability import diagnostics
from analysis.solarblue_lsf_serious_corpus.analyze_lsf import ANCHORS, CENTRAL_ANCHOR
from skysub.sky_decomp.lsf_spline2d import evaluate_bspline_basis, LSF_SPLINE2D_REPRESENTATION
from skysub.sky_decomp.lsf_surface_iterative import _wave_fingerprint
from skysub.sky_decomp.lsf_reference import validate_coefficients

DEST = OUT / "lsf_calibration"
METRICS = ["fwhm_angstrom", "w50_angstrom", "w90_angstrom", "centroid_angstrom"]
TRANSITION = 2.5
SCALE_FLOOR = 1e-6


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def dispersion_average(coefficients, nights):
    """A single robust weight per surface; equal night mass before downweighting."""
    x = np.asarray(coefficients, dtype=float)
    validate_coefficients(x)
    flat = x.reshape(len(x), -1)
    night_ids, inverse, counts = np.unique(nights, return_inverse=True, return_counts=True)
    center = np.median(flat, axis=0)
    scale = np.maximum(1.4826 * np.median(np.abs(flat - center), axis=0), SCALE_FLOOR)
    distance = np.sqrt(np.mean(((flat - center) / scale) ** 2, axis=1))
    robust = np.minimum(1, TRANSITION / np.maximum(distance, np.finfo(float).tiny)) if len(x) >= 3 else np.ones(len(x))
    base = 1 / counts[inverse]
    weights = base * robust
    weights /= weights.sum()
    average = np.einsum("n,nkb->kb", weights, x)
    validate_coefficients(average)
    deviation = flat - average.ravel()
    factors = np.sqrt(weights)[:, None] * deviation
    denominator = max(0.0, 1 - weights @ weights)
    scatter = np.sqrt(np.sum(factors ** 2, axis=0) / denominator) if len(x) > 1 and denominator > 0 else np.full(flat.shape[1], np.nan)
    night_factors = np.zeros((len(night_ids), flat.shape[1]))
    np.add.at(night_factors, inverse, weights[:, None] * deviation)
    conditional_error = (np.sqrt(np.sum(night_factors ** 2, axis=0) * len(night_ids) / (len(night_ids) - 1))
                         if len(night_ids) > 1 else np.full(flat.shape[1], np.nan))
    return dict(coefficient=average, weights=weights, base=base, robust=robust, distance=distance,
                scale=scale.reshape(x.shape[1:]), scatter=scatter.reshape(x.shape[1:]),
                conditional_error=conditional_error.reshape(x.shape[1:]), factors=factors,
                night_ids=night_ids, night_factors=night_factors, denominator=denominator,
                n_eff=1 / (weights @ weights))


def verified_selection():
    """Check the exact Gaia source, checksums and current SFrame cache provenance."""
    provenance_path = OUT / "gaia_photometry_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if Path(provenance["source_directory"]).resolve() != (DATA / "gaia/fibers").resolve():
        raise ValueError("Gaia provenance is not the user-authorized source")
    for record in provenance["files"]:
        gaia = Path(record["gaia_path"])
        if not gaia.resolve().is_relative_to((DATA / "gaia/fibers").resolve()) or sha256(gaia) != record["gaia_sha256"]:
            raise ValueError(f"Gaia provenance changed: {gaia}")
        stat = Path(record["sframe_path"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (record["sframe_size"], record["sframe_mtime_ns"]):
            raise ValueError("Current FLUX+SKY source changed; rebuild the photometry cache")
    rows = pd.read_parquet(OUT / "fiber_exposure_channel.parquet")
    if set(rows.expnum.unique()) != {r["expnum"] for r in provenance["files"]}:
        raise ValueError("Selection exposures disagree with verified provenance")
    for threshold, column in ((.01, "gaia_clean_1pct"), (.1, "gaia_clean_10pct")):
        expected = rows.gaia_verified & np.isfinite(rows.gaia_ratio) & rows.gaia_ratio.ge(0) & rows.gaia_ratio.lt(threshold)
        if not np.array_equal(expected, rows[column]):
            raise ValueError(f"Cached {column} disagrees with the strict flux-ratio criterion")
    selection = rows[rows.gaia_clean_10pct & rows.quality_ok & rows.calibration_epoch.ge(0)].copy()
    if selection.duplicated(["expnum", "fiberid", "channel"]).any():
        raise ValueError("Duplicate exposure-fiber-channel measurements")
    paths = fit_files()
    manifest = dict(gaia_source=str(DATA.resolve()), gaia_provenance_sha256=sha256(provenance_path),
                    selection_sha256=sha256(OUT / "fiber_exposure_channel.parquet"),
                    source_fits=[dict(path=str(p), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns,
                                     sha256=sha256(p)) for p in paths],
                    threshold=.1, strict_threshold=.01, transition=TRANSITION, scale_floor=SCALE_FLOOR,
                    builder_sha256=sha256(__file__), loader_sha256=sha256(ROOT / "skysub/sky_decomp/lsf_reference.py"))
    manifest["extraction_fingerprint"] = hashlib.sha256(json.dumps(
        {k: manifest[k] for k in ("gaia_provenance_sha256", "selection_sha256", "source_fits", "threshold")}, sort_keys=True).encode()).hexdigest()
    manifest["fingerprint"] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    print(f"Verified {len(provenance['files'])} Gaia/SFrame sources; selected {len(selection):,} channels", flush=True)
    return selection, paths, manifest


def extract(selection, paths, destination):
    wanted = selection.set_index(["expnum", "fiberid", "channel"])
    arrays, tables, source_configs, configuration, canonical_knots = [], [], [], None, {}
    for source_id, path in enumerate(paths):
        with fits.open(path, memmap=True) as hdus:
            inputs = frame(hdus["INPUT_META"].data)
            config = Table(hdus["LSF_CONFIG"].data)
            source_config = config.copy()
            source_config["source_id"] = np.full(len(config), source_id)
            source_configs.append(source_config)
            if configuration is None:
                configuration = config.copy()
            for ci, channel in enumerate(("B", "R", "Z")):
                keys = pd.MultiIndex.from_arrays([inputs.expnum, inputs.fiberid, np.full(len(inputs), channel)],
                                                names=["expnum", "fiberid", "channel"])
                indices = np.flatnonzero(keys.isin(wanted.index))
                if not len(indices):
                    continue
                cf = config[ci]
                canonical = configuration[ci]
                for key in ("lower", "upper", "degree", "n_basis", "n_knots", "tap_offsets",
                            "schema_version", "legacy_kernel_representation", "wave_n", "wave_min",
                            "wave_max", "wave_sha256", "offset_half_width_angstrom"):
                    a, b = np.asarray(cf[key]), np.asarray(canonical[key])
                    if a.dtype.kind in "fiu":
                        equal = np.array_equal(a, b, equal_nan=True)
                    else:
                        equal = np.array_equal(a, b)
                    if not equal:
                        raise ValueError(f"Incompatible source configuration: {path.name}, {channel}, {key}")
                if str(cf["legacy_kernel_representation"]).strip() != LSF_SPLINE2D_REPRESENTATION:
                    raise ValueError("Expected continuous M-spline LSF")
                coefficient = np.asarray(hdus["LSF_COEF"].data[indices, ci, :, :int(cf["n_basis"])], dtype=float)
                validate_coefficients(coefficient)
                knot = np.asarray(hdus["LSF_KNOTS"].data[indices, ci, :int(cf["n_knots"])], dtype=float)
                if not np.all(knot == knot[0]):
                    raise ValueError("Source knot vectors differ between exposures")
                if channel in canonical_knots and not np.array_equal(canonical_knots[channel], knot[0]):
                    raise ValueError("Source knot vectors differ between shards")
                canonical_knots[channel] = knot[0].copy()
                # Float32 archival roundoff is repaired once per unit-sum column.
                coefficient = np.maximum(coefficient, 0)
                coefficient /= coefficient.sum(axis=1, keepdims=True)
                validate_coefficients(coefficient, tolerance=2e-12)
                meta = wanted.loc[keys[indices]].reset_index()
                meta["source_id"], meta["fit_row"] = source_id, indices
                basis = evaluate_bspline_basis(np.array([CENTRAL_ANCHOR[channel]]), knot[0], int(cf["degree"]))
                values = {k: [] for k in METRICS}
                for block in np.array_split(coefficient, max(1, (len(coefficient) + 511) // 512)):
                    diag = diagnostics(np.einsum("wb,ntb->nwt", basis, block)[:, 0], float(cf["offset_half_width_angstrom"]))
                    for k in METRICS:
                        values[k].append(diag[k])
                for k in METRICS:
                    meta[k] = np.concatenate(values[k])
                meta["wavelength_angstrom"] = CENTRAL_ANCHOR[channel]
                arrays.append(coefficient.astype(np.float32))
                tables.append(meta)
        print(f"Extracted {path.name}", flush=True)
    meta = pd.concat(tables, ignore_index=True)
    coefficient = np.concatenate(arrays)
    if len(meta) != len(selection) or meta.duplicated(["expnum", "fiberid", "channel"]).any():
        raise ValueError("Extracted measurement coverage differs from the selection")
    if not np.isfinite(meta[METRICS]).all().all():
        raise ValueError("Nonfinite LSF diagnostics in selected measurements")
    meta["sample_id"] = np.arange(len(meta))
    meta.to_parquet(destination / "input_measurements.parquet", index=False)
    np.save(destination / "input_coefficients.npy", coefficient)
    configuration.write(destination / "source_configuration.fits", overwrite=True)
    vstack(source_configs).write(destination / "source_configurations_all.fits", overwrite=True)
    with fits.open(paths[0]) as hdus:
        first_source = str(hdus["INPUT_META"].data[0]["source_path"]).strip()
    with fits.open(first_source) as hdus:
        wave = np.asarray(hdus["WAVE"].data, dtype=float)
    cf = configuration[0]
    if _wave_fingerprint(wave) != (int(cf["wave_n"]), float(cf["wave_min"]), float(cf["wave_max"]), str(cf["wave_sha256"]).strip()):
        raise ValueError("Native wavelength grid disagrees with source LSF fingerprint")
    np.save(destination / "native_wave.npy", wave)
    (destination / "canonical_knots.json").write_text(json.dumps({c: k.tolist() for c, k in canonical_knots.items()}))
    return meta, coefficient, configuration, wave, canonical_knots


def write_product(scope, partition, meta, coefficients, configuration, wave, knots, manifest, destination):
    """One file per epoch/lunation/phase, with individual fiber-channel rows."""
    label = "survey" if not partition else "_".join(str(v).replace(":", "-").replace("%", "pct") for v in partition)
    path = destination / scope / f"lsf_{label}.fits"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        with fits.open(path) as previous:
            current = previous[0].header.get("BUILDFP") == manifest["fingerprint"]
        if current:
            with fits.open(path) as saved:
                reference = frame(saved["REFERENCES"].data)
            reference["scope"], reference["file"] = scope, str(path.relative_to(destination))
            return reference
    group_columns = ["fiberid", "channel"]
    order = np.lexsort((meta.channel.to_numpy(dtype=str), meta.fiberid.to_numpy()))
    data = meta.iloc[order].reset_index(drop=True)
    groups = data.groupby(group_columns, sort=False).indices
    means, scales, scatter, errors, knot_vectors, refs = [], [], [], [], [], []
    all_factors = np.empty((len(data), 44), dtype=np.float32)
    all_weights, all_base, all_robust, all_distance = (np.empty(len(data)) for _ in range(4))
    ref_ids = np.empty(len(data), dtype=np.int32)
    night_rows, night_factors = [], []
    for ref_id, ((fiber, channel), indices) in enumerate(groups.items()):
        rows = data.iloc[indices]
        x = np.asarray(coefficients[rows.sample_id.to_numpy()], dtype=float)
        x /= x.sum(axis=1, keepdims=True)
        result = dispersion_average(x, rows.mjd.to_numpy())
        means.append(result["coefficient"])
        scales.append(result["scale"])
        scatter.append(result["scatter"])
        errors.append(result["conditional_error"])
        knot_vectors.append(knots[channel])
        alpha = result["weights"]
        all_weights[indices], all_base[indices], all_robust[indices], all_distance[indices] = alpha, result["base"], result["robust"], result["distance"]
        all_factors[indices], ref_ids[indices] = result["factors"], ref_id
        center = rows.iloc[0]
        epoch = int(center.calibration_epoch) if scope != "survey" else -1
        n_epochs = rows.calibration_epoch.nunique()
        record = dict(reference_id=ref_id, fiberid=int(fiber), channel=channel,
                      telescope=str(center.telescope), spectrographid=int(center.spectrographid), calibration_epoch=epoch,
                      n_epochs=int(n_epochs), full_moon_utc=str(center.full_moon_utc) if scope == "lunation" else "ALL",
                      phase_bin="ALL", n_exposures=len(rows),
                      n_nights=len(result["night_ids"]), n_strict_1pct=int(rows.gaia_clean_1pct.sum()),
                      first_mjd=float(rows.observation_mjd_utc.min()), last_mjd=float(rows.observation_mjd_utc.max()),
                      weighted_mjd=float(alpha @ rows.observation_mjd_utc.to_numpy()), n_eff=float(result["n_eff"]),
                      scatter_denominator=float(result["denominator"]), downweighted=int(np.sum(result["robust"] < 1)),
                      robust_weighting=len(rows) >= 3, scatter_available=len(rows) > 1,
                      conditional_error_available=len(result["night_ids"]) > 1,
                      wavelength_angstrom=float(CENTRAL_ANCHOR[channel]))
        for metric in METRICS:
            values = rows[metric].to_numpy()
            average = alpha @ values
            record[metric + "_input_mean"] = float(average)
            record[metric + "_input_scatter"] = float(np.sqrt(alpha @ (values - average) ** 2 / result["denominator"])) if result["denominator"] > 0 else np.nan
        refs.append(record)
        for night, factor in zip(result["night_ids"], result["night_factors"]):
            night_rows.append(dict(reference_id=ref_id, mjd=int(night)))
            night_factors.append(factor)
    reference = pd.DataFrame(refs)
    means = np.asarray(means)
    validate_coefficients(means, tolerance=2e-12)
    for channel in ("B", "R", "Z"):
        indices = np.flatnonzero(reference.channel.eq(channel))
        if not len(indices):
            continue
        ci = ("B", "R", "Z").index(channel)
        basis = evaluate_bspline_basis(np.array([CENTRAL_ANCHOR[channel]]), knots[channel], int(configuration[ci]["degree"]))
        diag = diagnostics(np.einsum("wb,ntb->nwt", basis, means[indices])[:, 0], float(configuration[ci]["offset_half_width_angstrom"]))
        for metric in METRICS:
            reference.loc[indices, metric + "_reference"] = diag[metric]
    primary = fits.PrimaryHDU()
    primary.header.update(REFSCHEM=1, SCOPE=scope, APPLYOK=scope != "survey", GAIALIM=.1, STRICTG=.01,
                          WGTMETH="night-balanced robust RMS", HUBER=TRANSITION, SCALEMIN=SCALE_FLOOR,
                          FLATORD="C: offset then wavelength", BUILDFP=manifest["fingerprint"],
                          COVTYPE="empirical; conditional night mean", NSOURCE=len(data), NREF=len(reference))
    primary.header["COMMENT"] = "Survey pooling is descriptive; epoch/lunation products preserve calibration epochs."
    contributor_columns = ["sample_id", "expnum", "fiberid", "channel", "mjd", "observation_mjd_utc", "calibration_epoch",
                           "source_id", "fit_row", "gaia_ratio", "gaia_clean_1pct", "gaia_clean_10pct"]
    contributors = data[contributor_columns].copy()
    contributors["reference_id"], contributors["weight"] = ref_ids, all_weights
    contributors["base_weight"], contributors["robust_weight"], contributors["dispersion_distance"] = all_base, all_robust, all_distance
    hdus = [primary, fits.BinTableHDU(Table.from_pandas(reference), name="REFERENCES"),
            fits.ImageHDU(means, name="LSF_COEF"), fits.ImageHDU(np.asarray(knot_vectors), name="LSF_KNOTS"),
            fits.BinTableHDU(configuration, name="LSF_CONFIG"), fits.ImageHDU(wave, name="WAVE"),
            fits.ImageHDU(np.asarray(scales), name="COEF_SCALE"), fits.ImageHDU(np.asarray(scatter), name="COEF_SCATTER"),
            fits.ImageHDU(np.asarray(errors), name="COEF_CONDERR"),
            fits.BinTableHDU(Table.from_pandas(contributors), name="CONTRIBUTORS"),
            fits.ImageHDU(all_factors, name="SCATTER_FACTORS"),
            fits.BinTableHDU(Table(rows=night_rows), name="NIGHTS"),
            fits.ImageHDU(np.asarray(night_factors, dtype=np.float32), name="NIGHT_FACTORS"),
            fits.BinTableHDU(Table(rows=manifest["source_fits"]), name="SOURCES")]
    all_config = destination.parent / "lsf_reference_audit/cache/source_configurations_all.fits"
    if all_config.exists():
        hdus.append(fits.BinTableHDU(Table.read(all_config), name="SOURCE_CONFIG"))
    hdus[2].header["BUNIT"] = "dimensionless"
    hdus[2].header["AXES"] = "reference,offset_component,wavelength_basis"
    hdus[3].header["BUNIT"] = "Angstrom"
    hdus[5].header["BUNIT"] = "Angstrom"
    from analysis.lsf_epoch_coverage.compact_lsf_references import write_reference_hdus
    write_reference_hdus(hdus, path, destination)
    reference["scope"], reference["file"] = scope, str(path.relative_to(destination))
    print(f"Wrote {path.relative_to(destination)}: {len(reference):,} references", flush=True)
    return reference


def publish(catalog, manifest, selection, destination):
    from analysis.lsf_epoch_coverage.reference_views import publish_reference_views
    return publish_reference_views(catalog, manifest, selection, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEST)
    args = parser.parse_args()
    destination = args.output.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    selection, paths, manifest = verified_selection()
    cache = destination.parent / "lsf_reference_audit/cache"
    cache.mkdir(parents=True, exist_ok=True)
    (destination.parent / "lsf_reference_audit/publication_records").mkdir(parents=True, exist_ok=True)
    provenance_path = destination.parent / "lsf_reference_audit/publication_records/provenance.json"
    old = json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    if old.get("extraction_fingerprint") == manifest["extraction_fingerprint"] and (cache / "canonical_knots.json").exists():
        meta = pd.read_parquet(cache / "input_measurements.parquet")
        coefficients = np.load(cache / "input_coefficients.npy", mmap_mode="r")
        configuration = Table.read(cache / "source_configuration.fits")
        wave = np.load(cache / "native_wave.npy")
        knots = {c: np.asarray(k) for c, k in json.loads((cache / "canonical_knots.json").read_text()).items()}
    else:
        meta, coefficients, configuration, wave, knots = extract(selection, paths, cache)
    provenance_path.write_text(json.dumps(manifest, indent=2) + "\n")
    catalogs = []
    columns = dict(survey=[], epoch=["calibration_epoch"], lunation=["calibration_epoch", "full_moon_utc"])
    for scope in columns:
        parts = [((), meta)] if not columns[scope] else meta.groupby(columns[scope], sort=True)
        for partition, rows in parts:
            if not isinstance(partition, tuple):
                partition = (partition,)
            catalogs.append(write_product(scope, partition, rows, coefficients, configuration, wave, knots, manifest, destination))
    catalog = pd.concat(catalogs, ignore_index=True)
    statistics = publish(catalog, manifest, selection, destination)
    print(json.dumps(statistics, indent=2), flush=True)


if __name__ == "__main__":
    main()
