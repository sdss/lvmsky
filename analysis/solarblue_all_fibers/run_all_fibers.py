#!/usr/bin/env python3
"""Resume SolarBlue decomposition for every Science/Sky fiber in bright-Moon exposures."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Column, Table, vstack
from astropy.time import Time

from analysis.solarblue_lsf_serious_corpus.build_corpus import (
    STACK,
    _combined_ivar,
    _flux_scale,
    gaia_columns,
    local_path,
    repeated_meta,
    slit_indices,
)


OUT = Path("/data/sas/sdsswork/users/u6058164/lvm_solarblue_all_fibers_drp132")
FIT_MODEL = "palacecorr-aijc-vnf-split-zodi-lsf-spline2d-solarblue"
FIT_SUFFIX = "_palacecorr_aijc_vnf_split_zodi_lsf_spline2d_solarblue"
MOON_ELEVATION_MIN_DEG = 10.0
MOON_FLI_MIN = 0.5
TELESCOPES = ("Sci", "SkyE", "SkyW")
REFERENCE_COUNT = 9
DROP_PARAMETER_HDUS = {"COEF_ERR", "COEF_COV_MOON", "COEF_COV_ZODI"}
DROP_REFERENCE_HDUS = DROP_PARAMETER_HDUS | {"BESTFIT", "FLUX_SIGMA_TOTAL", "VECTOR_O2"}
FIT_META_COLUMNS = (
    "source_row", "role", "input_valid", "fit_elapsed_sec", "fit_status",
    "o2_fit_elapsed_sec", "o2_fit_status", "o2_prefit_amp", "o2_valid_frac",
    "peak_memory_mb", "r2", "reduced_chi2", "resid_level", "rms_resid",
    "t_o2", "t_o2_err",
)
LSF_DIAGNOSTIC_COLUMNS = (
    "spectrum_index", "channel", "available", "completed_cycles", "converged",
    "n_valid_pixels", "chi2_red", "rms_resid", "center_pix_min",
    "center_pix_median", "center_pix_max", "sigma_pix_min", "sigma_pix_median",
    "sigma_pix_max", "runtime_sec",
)
LSF_CONFIG_COLUMNS = (
    "channel", "lower", "upper", "degree", "n_basis", "n_knots", "tap_offsets",
    "schema_version", "requested_cycles", "legacy_kernel_representation", "wave_n",
    "wave_min", "wave_max", "wave_sha256", "offset_half_width_angstrom",
    "offset_roughness_fraction", "line_weight", "skyline_cumulative_fraction",
    "skyline_half_width_angstrom", "huber_transition_sigma", "config_n_basis",
    "config_degree", "roughness_fraction", "fallback_prior_fraction",
    "information_prior_max_boost", "background_degree", "blue_fit_lower",
    "convergence_lsf_rtol", "convergence_chi2_rtol",
    "convergence_centroid_atol_angstrom", "convergence_lsf_metric", "n_pixels",
    "fit_lower", "fit_upper", "knots_fixed", "model",
)
INPUT_META_COLUMNS = (
    "source_path", "expnum", "mjd", "obstime", "date_obs", "pwv_med",
    "sci_ra", "sci_dec", "skye_ra", "skye_dec", "skyw_ra", "skyw_dec",
    "sky_near_label", "sky_far_label", "sci_airmass", "skye_airmass", "skyw_airmass",
    "sci_moon_sep", "skye_moon_sep", "skyw_moon_sep", "moon_alt", "sun_alt",
    "moon_phase", "moon_fli", "source_fiber_index", "reference_fiber",
    "fit_valid_fraction", "fiberid", "spectrographid", "targettype", "telescope",
    "ypix_b", "ypix_r", "ypix_z", "fibstatus", "ra", "dec", "gaia_status",
    "gaia_clean_1pct", "gaia_ratio", "gaia_n_sources", "corpus_row",
)


def exposure_selection(
    stack: Path = STACK,
    moon_alt_min: float = MOON_ELEVATION_MIN_DEG,
    moon_fli_min: float = MOON_FLI_MIN,
) -> pd.DataFrame:
    """Select exposures above the requested Moon elevation and FLI limits."""
    with fits.open(stack, memmap=True) as hdul:
        frame = Table(hdul["META"].data).to_pandas()
    for name in frame.select_dtypes(include=["object", "str"]):
        frame[name] = frame[name].map(
            lambda value: value.decode().strip() if isinstance(value, bytes) else str(value).strip())
    frame["source_path"] = frame.path.map(lambda value: str(local_path(value)))
    frame["source_available"] = frame.source_path.map(lambda value: Path(value).is_file())
    frame["selected"] = (
        (frame.moon_alt > float(moon_alt_min))
        & (frame.moon_fli > float(moon_fli_min))
        & frame.source_available
    )
    return frame


def write_immutable(path: Path, text: str) -> None:
    """Create a selection artifact once; later resumes must reproduce it exactly."""
    if path.is_file():
        if path.read_text(encoding="utf-8") != text:
            raise RuntimeError(f"Immutable selection artifact differs: {path}")
        return
    path.write_text(text, encoding="utf-8")
    path.chmod(0o444)


def read_selected_manifest(path: Path, moon_alt_min: float, moon_fli_min: float) -> pd.DataFrame:
    """Load and validate an explicitly frozen exposure selection."""
    selected = pd.read_csv(path).sort_values(["mjd", "expnum"])
    required = {"expnum", "mjd", "moon_alt", "moon_fli", "source_path", "source_available"}
    if missing := required.difference(selected.columns):
        raise ValueError(f"Selected manifest lacks columns: {sorted(missing)}")
    valid = (
        (selected.moon_alt > float(moon_alt_min))
        & (selected.moon_fli > float(moon_fli_min))
        & selected.source_available
    )
    if selected.expnum.duplicated().any() or not valid.all():
        raise ValueError("Selected manifest contains duplicate or ineligible exposures")
    return selected


def selected_fibers(slitmap: Table) -> tuple[np.ndarray, np.ndarray]:
    """Return all Science/Sky rows, with the established nine references first."""
    telescope = np.asarray(slitmap["telescope"]).astype(str)
    all_rows = np.flatnonzero(np.isin(telescope, TELESCOPES))
    reference = np.asarray([item[0] for item in slit_indices(slitmap)], dtype=int)
    remaining = all_rows[~np.isin(all_rows, reference)]
    return np.r_[reference, remaining], reference


def build_exposure_input(record: dict, destination: Path) -> dict:
    """Write one temporary mixed-fiber fit input; the first nine rows are references."""
    source = Path(record["source_path"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp.fits")
    with fits.open(source, memmap=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
        slitmap = Table(hdul["SLITMAP"].data)
        indices, reference = selected_fibers(slitmap)
        scale = _flux_scale(hdul[0].header, hdul["FLUX"])
        flux = np.asarray(hdul["FLUX"].data[indices], dtype=np.float32) * scale
        sky = np.asarray(hdul["SKY"].data[indices], dtype=np.float32) * scale
        ivar = np.asarray(hdul["IVAR"].data[indices], dtype=np.float32) / scale**2
        sky_ivar = np.asarray(hdul["SKY_IVAR"].data[indices], dtype=np.float32) / scale**2
        mask = np.asarray(hdul["MASK"].data[indices], dtype=np.uint8)
        lsf = np.asarray(hdul["LSF"].data[indices], dtype=np.float32)
    observed = flux + sky
    combined_ivar = _combined_ivar(ivar, sky_ivar)
    good = (np.isfinite(observed) & np.isfinite(combined_ivar) & (combined_ivar > 0)
            & np.isfinite(lsf) & (lsf > 0) & (mask == 0))
    combined_ivar[~good] = 0.0
    fiber_meta = slitmap[indices].copy()
    fiber_meta.add_column(Column(indices, name="source_fiber_index"), index=0)
    fiber_meta.add_column(Column(np.arange(len(indices)) < REFERENCE_COUNT,
                                 name="reference_fiber"), index=1)
    fiber_meta.add_column(Column(good.mean(axis=1), name="fit_valid_fraction"), index=2)
    for name, values in gaia_columns(int(record["expnum"]),
                                     np.asarray(fiber_meta["fiberid"])) .items():
        fiber_meta[name] = values
    meta = repeated_meta(pd.Series(record), fiber_meta)
    meta["corpus_row"] = np.arange(len(meta), dtype=np.int64)
    primary = fits.PrimaryHDU()
    primary.header["DRPVER"] = "1.3.2"
    primary.header["EXPNUM"] = int(record["expnum"])
    primary.header["MOONALT"] = float(record["moon_alt"])
    primary.header["MOONFLI"] = float(record["moon_fli"])
    primary.header["NFIBER"] = len(indices)
    primary.header["NREF"] = REFERENCE_COUNT
    primary.header["CREATED"] = Time.now().isot
    fits.HDUList([
        primary,
        fits.ImageHDU(wave, name="WAVE"),
        fits.ImageHDU(observed, name="FLUX_SCI"),
        fits.ImageHDU(combined_ivar, name="IVAR_SCI"),
        fits.ImageHDU(mask, name="MASK_SCI"),
        fits.ImageHDU(lsf, name="LSF_SCI"),
        fits.BinTableHDU(meta, name="META"),
    ]).writeto(temporary, overwrite=True, checksum=True)
    os.replace(temporary, destination)
    return {
        "fibers": len(indices),
        "science": int(np.sum(np.asarray(fiber_meta["telescope"]).astype(str) == "Sci")),
        "sky": int(np.sum(np.isin(np.asarray(fiber_meta["telescope"]).astype(str), ("SkyE", "SkyW")))),
        "reference_fiberids": np.asarray(slitmap[reference]["fiberid"], dtype=int).tolist(),
    }


def output_paths(output: Path, expnum: int) -> tuple[Path, Path, Path, Path]:
    stem = f"solarblue-all-fibers-{expnum:08d}"
    input_path = output / "work" / f"{stem}.fits"
    parameters = output / "parameters" / f"{stem}_sci_meta_coef{FIT_SUFFIX}.fits"
    reference = output / "reference_models" / f"{stem}_decomp_sci{FIT_SUFFIX}.fits"
    cache = output / "parameters" / f".{stem}{FIT_SUFFIX}_compact_cache"
    return input_path, parameters, reference, cache


def build_shard_input(records: list[dict], destination: Path) -> list[dict]:
    """Combine exposure inputs, placing every nine-row reference block first."""
    part_dir = destination.parent / f".{destination.stem}-parts"
    part_dir.mkdir(parents=True, exist_ok=True)
    parts, info = [], []
    for record in records:
        part = part_dir / f"{int(record['expnum']):08d}.fits"
        detail = build_exposure_input(record, part) if not part.is_file() else {
            "fibers": int(fits.getheader(part, 0)["NFIBER"])}
        parts.append(part)
        info.append(dict(expnum=int(record["expnum"]), **detail))
    opened = [fits.open(path, memmap=True) for path in parts]
    try:
        wave = np.asarray(opened[0]["WAVE"].data, dtype=np.float32)
        for hdul in opened[1:]:
            if not np.array_equal(wave, hdul["WAVE"].data):
                raise ValueError("WAVE differs between exposures in one shard")
        order = [(hdul, slice(0, REFERENCE_COUNT)) for hdul in opened]
        order += [(hdul, slice(REFERENCE_COUNT, None)) for hdul in opened]
        arrays = {
            name: np.concatenate([np.asarray(hdul[name].data[rows]) for hdul, rows in order])
            for name in ("FLUX_SCI", "IVAR_SCI", "MASK_SCI", "LSF_SCI")
        }
        meta = vstack([Table(hdul["META"].data[rows]) for hdul, rows in order],
                      metadata_conflicts="silent")
        meta["corpus_row"] = np.arange(len(meta), dtype=np.int64)
        primary = fits.PrimaryHDU()
        primary.header["DRPVER"] = "1.3.2"
        primary.header["NEXP"] = len(records)
        primary.header["NFIBER"] = len(meta)
        primary.header["NREF"] = REFERENCE_COUNT * len(records)
        primary.header["CREATED"] = Time.now().isot
        temporary = destination.with_suffix(".tmp.fits")
        fits.HDUList([
            primary, fits.ImageHDU(wave, name="WAVE"),
            *(fits.ImageHDU(arrays[name], name=name) for name in arrays),
            fits.BinTableHDU(meta, name="META"),
        ]).writeto(temporary, overwrite=True, checksum=True)
        os.replace(temporary, destination)
    finally:
        for hdul in opened:
            hdul.close()
        shutil.rmtree(part_dir, ignore_errors=True)
    return info


def shard_output_paths(output: Path, shard_index: int, records: list[dict]):
    first, last = int(records[0]["expnum"]), int(records[-1]["expnum"])
    stem = f"solarblue-shard-{shard_index:04d}-{first:08d}-{last:08d}"
    input_path = output / "work" / f"{stem}.fits"
    parameters = output / "parameters" / f"{stem}_sci_meta_coef{FIT_SUFFIX}.fits"
    reference = output / "reference_models" / f"{stem}_decomp_sci{FIT_SUFFIX}.fits"
    cache = output / "parameters" / f".{stem}{FIT_SUFFIX}_compact_cache"
    return input_path, parameters, reference, cache


def attach_input_context(input_path: Path, output_path: Path, rows=None, include_wave=False):
    """Atomically attach the fiber identity table and optional wavelength grid."""
    with fits.open(input_path, memmap=True) as source, fits.open(output_path, memmap=False) as product:
        vector_names = ("INPUT_IVAR", "INPUT_MASK") if include_wave else ()
        if ("INPUT_META" in product and (not include_wave or "WAVE" in product)
                and all(name in product for name in vector_names)):
            return
        hdus = [hdu.copy() for hdu in product]
        if include_wave and "WAVE" not in product:
            hdus.append(fits.ImageHDU(np.asarray(source["WAVE"].data), name="WAVE"))
        if include_wave and "INPUT_IVAR" not in product:
            hdus.append(fits.ImageHDU(
                np.asarray(source["IVAR_SCI"].data[:rows], dtype=np.float32), name="INPUT_IVAR"))
        if include_wave and "INPUT_MASK" not in product:
            hdus.append(fits.ImageHDU(
                np.asarray(source["MASK_SCI"].data[:rows], dtype=np.uint8), name="INPUT_MASK"))
        if "INPUT_META" not in product:
            data = source["META"].data if rows is None else source["META"].data[:rows]
            hdus.append(fits.BinTableHDU(Table(data), name="INPUT_META"))
    temporary = output_path.with_suffix(".context.tmp.fits")
    fits.HDUList(hdus).writeto(temporary, overwrite=True, checksum=True)
    os.replace(temporary, output_path)


def compact_product(path: Path, reference=False):
    """Keep requested science products and store floating-point values as float32."""
    drop = DROP_REFERENCE_HDUS if reference else DROP_PARAMETER_HDUS
    with fits.open(path, memmap=False) as source:
        if source[0].header.get("STORFMT") == "minimal-f32-v2":
            return
        hdus = []
        extras = []
        for hdu in source:
            if hdu.name in drop:
                continue
            if isinstance(hdu, fits.PrimaryHDU):
                copied = hdu.copy()
                copied.header["STORFMT"] = "minimal-f32-v2"
            elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                data = np.asarray(hdu.data)
                if data.dtype.kind == "f" and data.dtype.itemsize > 4:
                    data = data.astype(np.float32)
                copied = fits.ImageHDU(data, header=hdu.header, name=hdu.name)
            elif isinstance(hdu, fits.BinTableHDU):
                table = Table(hdu.data)
                if hdu.name == "META":
                    failed = np.char.strip(np.asarray(table["fit_status"]).astype(str)) != "Solved"
                    if np.any(failed):
                        extras.append(fits.BinTableHDU(table[failed], name="FIT_FAILURES"))
                    table = table[[name for name in FIT_META_COLUMNS if name in table.colnames]]
                elif hdu.name == "LSF_META":
                    failed = np.char.strip(np.asarray(table["status"]).astype(str)) != "Solved"
                    if np.any(failed):
                        extras.append(fits.BinTableHDU(table[failed], name="LSF_FAILURES"))
                    channels = np.asarray(table["channel"]).astype(str)
                    first = [np.flatnonzero(channels == channel)[0]
                             for channel in np.unique(channels)]
                    config = table[first][
                        [name for name in LSF_CONFIG_COLUMNS if name in table.colnames]]
                    extras.append(fits.BinTableHDU(config, name="LSF_CONFIG"))
                    table = table[
                        [name for name in LSF_DIAGNOSTIC_COLUMNS if name in table.colnames]]
                elif hdu.name == "INPUT_META":
                    table = table[[name for name in INPUT_META_COLUMNS if name in table.colnames]]
                for name in table.colnames:
                    dtype = np.asarray(table[name]).dtype
                    if dtype.kind == "f" and dtype.itemsize > 4:
                        table[name] = np.asarray(table[name], dtype=np.float32)
                copied = fits.BinTableHDU(table, name=hdu.name)
            else:
                copied = hdu.copy()
            hdus.append(copied)
        hdus.extend(extras)
    temporary = path.with_suffix(".compact.tmp.fits")
    fits.HDUList(hdus).writeto(temporary, overwrite=True, checksum=True)
    os.replace(temporary, path)


def fit_shard(records: list[dict], shard_index: int, output: Path,
              workers: int, chunk_size: int) -> list[dict]:
    """Fit one resumable multi-exposure shard and return per-exposure status rows."""
    from skysub import decompose_parallel as dp

    input_path, parameters, reference, cache = shard_output_paths(output, shard_index, records)
    started = time.time()
    if not parameters.is_file() or not reference.is_file():
        input_path.parent.mkdir(parents=True, exist_ok=True)
        info = build_shard_input(records, input_path) if not input_path.is_file() else [
            {"expnum": int(record["expnum"]), "fibers": 1920} for record in records]
    else:
        info = [{"expnum": int(record["expnum"]), "fibers": 1920} for record in records]
    common = dict(data_file=input_path, palace_dir=None, lsf_sigma=0.5, factor=1e14,
                  chunk_size=chunk_size, max_in_flight=workers, fit_model=FIT_MODEL,
                  kinds=("sci",), start_method="fork", pin_workers=True)
    n_reference = REFERENCE_COUNT * len(records)
    if not reference.is_file():
        reference.parent.mkdir(parents=True, exist_ok=True)
        dp.run(n_workers=min(workers, n_reference), output_dir=reference.parent,
               limit=n_reference, compact_only=False, **common)
    attach_input_context(input_path, reference, rows=n_reference, include_wave=True)
    compact_product(reference, reference=True)
    if not parameters.is_file():
        parameters.parent.mkdir(parents=True, exist_ok=True)
        dp.run(n_workers=workers, output_dir=parameters.parent, limit=None,
               compact_only=True, **common)
    attach_input_context(input_path, parameters)
    compact_product(parameters)
    if parameters.is_file() and reference.is_file():
        input_path.unlink(missing_ok=True)
        shutil.rmtree(cache, ignore_errors=True)
    elapsed = time.time() - started
    return [dict(row, status="complete", shard=shard_index,
                 elapsed_seconds=elapsed, parameter_file=str(parameters),
                 reference_file=str(reference)) for row in info]


def fit_exposure(record: dict, output: Path, workers: int, chunk_size: int) -> dict:
    """Fit or resume one exposure, preserving compact all-fiber and full reference outputs."""
    from skysub import decompose_parallel as dp

    expnum = int(record["expnum"])
    input_path, parameters, reference, cache = output_paths(output, expnum)
    started = time.time()
    info = None
    if not parameters.is_file() or not reference.is_file():
        info = build_exposure_input(record, input_path) if not input_path.is_file() else {
            "fibers": fits.getheader(input_path, 0)["NFIBER"]}
    common = dict(
        data_file=input_path, palace_dir=None, lsf_sigma=0.5, factor=1e14,
        chunk_size=chunk_size, max_in_flight=workers, fit_model=FIT_MODEL,
        kinds=("sci",), start_method="fork", pin_workers=True,
    )
    if not reference.is_file():
        reference.parent.mkdir(parents=True, exist_ok=True)
        dp.run(n_workers=min(workers, REFERENCE_COUNT), output_dir=reference.parent,
               limit=REFERENCE_COUNT, compact_only=False, **common)
    attach_input_context(input_path, reference, rows=REFERENCE_COUNT, include_wave=True)
    compact_product(reference, reference=True)
    if not parameters.is_file():
        parameters.parent.mkdir(parents=True, exist_ok=True)
        dp.run(n_workers=workers, output_dir=parameters.parent, limit=None,
               compact_only=True, **common)
    attach_input_context(input_path, parameters)
    compact_product(parameters)
    if parameters.is_file() and reference.is_file():
        input_path.unlink(missing_ok=True)
        shutil.rmtree(cache, ignore_errors=True)
    return dict(expnum=expnum, status="complete", elapsed_seconds=time.time() - started,
                parameter_file=str(parameters), reference_file=str(reference), **(info or {}))


def write_status(output: Path, selection: pd.DataFrame, completed: list[dict], active=None):
    done = {int(row["expnum"]): row for row in completed}
    payload = {
        "updated_utc": Time.now().isot,
        "selected_exposures": int(len(selection)),
        "completed_exposures": len(done),
        "completed_fits": int(sum(int(row.get("fibers", 0)) for row in done.values())),
        "active_expnum": active,
        "fraction_exposures": len(done) / len(selection) if len(selection) else 0.0,
    }
    temporary = output / "status.tmp.json"
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, output / "status.json")
    if done:
        pd.DataFrame(done.values()).sort_values("expnum").to_csv(
            output / "completed_exposures.csv", index=False)
    print(json.dumps(payload), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--workers", type=int, default=96)
    parser.add_argument("--chunk-size", type=int, default=1)
    parser.add_argument("--limit-exposures", type=int)
    parser.add_argument("--expnum", type=int)
    parser.add_argument("--exposures-per-shard", type=int, default=32)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--moon-alt-min", type=float, default=MOON_ELEVATION_MIN_DEG)
    parser.add_argument("--moon-fli-min", type=float, default=MOON_FLI_MIN)
    parser.add_argument("--selected-manifest", type=Path)
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.selected_manifest:
        manifest_path = args.selected_manifest.expanduser().resolve()
        selected = read_selected_manifest(
            manifest_path,
            args.moon_alt_min,
            args.moon_fli_min,
        )
        selected_csv = manifest_path.read_text(encoding="utf-8")
    else:
        audit = exposure_selection(
            moon_alt_min=args.moon_alt_min,
            moon_fli_min=args.moon_fli_min,
        )
        write_immutable(output / "exposure_selection_audit.csv", audit.to_csv(index=False))
        selected = audit[audit.selected].sort_values(["mjd", "expnum"])
        selected_csv = selected.to_csv(index=False)
    write_immutable(output / "selected_exposures.csv", selected_csv)
    full_selected = selected
    if args.expnum is not None:
        selected = selected[selected.expnum == args.expnum]
    if args.limit_exposures is not None:
        selected = selected.head(args.limit_exposures)
    summary = {
        "selection": (
            f"moon_alt > {args.moon_alt_min:g} deg and "
            f"moon_fli > {args.moon_fli_min:g}"
        ),
        "selected_exposures": int(len(full_selected)),
        "selected_nights": int(full_selected.mjd.nunique()),
        "fit_model": FIT_MODEL,
        "fiber_telescopes": TELESCOPES,
        "full_vector_reference_fibers_per_exposure": REFERENCE_COUNT,
        "exposures_per_shard": args.exposures_per_shard,
    }
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    provenance = {
        **summary,
        "source_stack": str(Path(STACK).resolve()),
        "selection_inventory": (
            str(args.selected_manifest.expanduser().resolve())
            if args.selected_manifest else str(Path(STACK).resolve())
        ),
        "selected_manifest": "selected_exposures.csv",
        "selected_manifest_sha256": hashlib.sha256(
            selected_csv.encode("utf-8")
        ).hexdigest(),
        "comparison_operators": {"moon_alt": ">", "moon_fli": ">"},
    }
    write_immutable(
        output / "selection_provenance.json",
        json.dumps(provenance, indent=2) + "\n",
    )
    print(json.dumps(summary, indent=2), flush=True)
    if args.selection_only:
        return
    completed_path = output / "completed_exposures.csv"
    completed = (pd.read_csv(completed_path).to_dict("records")
                 if completed_path.is_file() else [])
    completed_ids = {int(row["expnum"]) for row in completed}
    remaining = [record for record in selected.to_dict("records")
                 if int(record["expnum"]) not in completed_ids]
    if args.expnum is not None:
        for record in remaining:
            write_status(output, full_selected, completed, active=[int(record["expnum"])])
            completed.append(fit_exposure(record, output, args.workers, args.chunk_size))
            write_status(output, full_selected, completed)
        return
    for shard_index, start in enumerate(range(0, len(remaining), args.exposures_per_shard)):
        records = remaining[start:start + args.exposures_per_shard]
        active = [int(record["expnum"]) for record in records]
        write_status(output, full_selected, completed, active=active)
        completed.extend(fit_shard(records, shard_index, output, args.workers, args.chunk_size))
        write_status(output, full_selected, completed)


if __name__ == "__main__":
    main()
