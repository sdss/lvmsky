#!/usr/bin/env python3
"""Build the DRP 1.3.2 bright-Moon individual-Science-fiber corpus."""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import json
import mmap
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Column, Table, vstack
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).parents[2] / "skysub/medians_computation/src"))
from lvm_medians.stack import (  # noqa: E402
    _combined_ivar,
    _flux_scale,
    gaia_ratio_clean_mask,
)


SAS = Path(os.environ.get("SAS_BASE_DIR", "/data/sas"))
DATA = SAS / "sdsswork/users/u6058164/gaia_fiber_matching_1.3.2"
STACK = DATA / "lvmsframe_median_stack_1.3.2_gaia1over100.fits"
OUT = SAS / "sdsswork/users/u6058164/lvm_solarblue_lsf_serious_corpus_drp132"
GAIA_DIRS = (OUT / "gaia/fibers", DATA / "gaia/fibers",
             Path(__file__).parents[1] / "bright_moon_triplet/gaia/fibers")
QUANTILES = (0.1, 0.5, 0.9)
CONTINUUM_WINDOW = (5400.0, 5500.0)
MAX_DAYS_FROM_FULL = 5.0
MIN_RELATIVE_CONTINUUM = 0.3
MIN_VALID_FRACTION = 0.8


def local_path(path: str) -> Path:
    return SAS / ("sdsswork/" + str(path).split("/sdsswork/", 1)[1])


def native(values):
    values = np.asarray(values)
    return values.astype(values.dtype.newbyteorder("="), copy=True)


def select_exposures(stack: Path = STACK) -> pd.DataFrame:
    from analysis.bright_moon.analyze import full_moons

    """Apply the agreed geometry, lunation and measured-continuum selection."""
    with fits.open(stack, memmap=True) as hdul:
        meta = hdul["META"].data
        frame = Table(meta).to_pandas()
        for name in frame.select_dtypes("object"):
            frame[name] = frame[name].map(
                lambda value: value.decode().strip() if isinstance(value, bytes) else str(value).strip())
        midpoint = Time(frame.obstime.to_numpy(dtype=str), format="isot", scale="utc").mjd
        moons = full_moons(midpoint.min(), midpoint.max())
        nearest = np.abs(midpoint[:, None] - moons).argmin(axis=1)
        frame["midpoint_mjd_utc"] = midpoint
        frame["full_moon_mjd_utc"] = moons[nearest]
        frame["full_moon_utc"] = Time(moons[nearest], format="mjd").isot
        frame["days_from_full_moon"] = midpoint - moons[nearest]
        frame["geometry_selected"] = (
            (frame.moon_alt > 0.0)
            & (frame.sun_alt < -18.0)
            & (frame.days_from_full_moon.abs() <= MAX_DAYS_FROM_FULL)
        )
        wave = np.asarray(hdul["WAVE"].data)
        band = np.flatnonzero((wave >= CONTINUUM_WINDOW[0]) & (wave <= CONTINUUM_WINDOW[1]))
        window = slice(band[0], band[-1] + 1)
        rows = np.flatnonzero(frame.geometry_selected)
        continuum = np.full(len(frame), np.nan)
        hdul._file._mmap.madvise(mmap.MADV_RANDOM)
        for chunk in np.array_split(rows, max(1, (len(rows) + 255) // 256)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                continuum[chunk] = np.nanmedian(
                    np.asarray(hdul["FLUX_SCI"].data[chunk, window], dtype=float), axis=1)
    frame["continuum_5400_5500"] = continuum
    log_continuum = np.log10(frame.continuum_5400_5500.where(frame.continuum_5400_5500 > 0))
    pool = log_continuum[frame.geometry_selected]
    centre = pool.median()
    mad_sigma = 1.4826 * (pool - centre).abs().median()
    frame["continuum_gross_outlier"] = log_continuum - centre > 8.0 * mad_sigma
    usable = frame.geometry_selected & ~frame.continuum_gross_outlier
    frame["lunation_reference_continuum"] = np.nan
    frame.loc[usable, "lunation_reference_continuum"] = frame.loc[usable].groupby(
        "full_moon_mjd_utc").continuum_5400_5500.transform(
            lambda values: float(np.nanmedian(np.sort(values.dropna())[-3:]))
            if values.notna().any() else np.nan)
    frame["relative_lunar_continuum"] = (
        frame.continuum_5400_5500 / frame.lunation_reference_continuum)
    paths = frame.path.map(local_path)
    frame["source_available"] = paths.map(Path.is_file)
    frame["source_path"] = paths.map(str)
    frame["selected"] = (
        usable & frame.source_available
        & (frame.relative_lunar_continuum >= MIN_RELATIVE_CONTINUUM)
    )
    return frame


def slit_indices(slitmap: Table) -> list[tuple[int, int, float, int]]:
    """Return three good Science fibers at 10/50/90% of each spectrograph slit."""
    telescope = np.asarray(slitmap["telescope"]).astype(str)
    selected = []
    for spectrograph in (1, 2, 3):
        candidates = np.flatnonzero(
            (telescope == "Sci")
            & (np.asarray(slitmap["spectrographid"]) == spectrograph)
            & (np.asarray(slitmap["fibstatus"]) == 0)
        )
        if len(candidates) < 3:
            raise ValueError(f"spectrograph {spectrograph} has fewer than three good Science fibers")
        candidates = candidates[np.argsort(np.asarray(slitmap["ypix_b"])[candidates], kind="stable")]
        positions = np.rint(np.asarray(QUANTILES) * (len(candidates) - 1)).astype(int)
        selected.extend(
            (int(candidates[position]), spectrograph, quantile, rank)
            for rank, (quantile, position) in enumerate(zip(QUANTILES, positions), 1)
        )
    return selected


def gaia_path(expnum: int) -> Path | None:
    name = f"lvmGAIA-fibers-{expnum:08d}.fits"
    return next((directory / name for directory in GAIA_DIRS if (directory / name).is_file()), None)


def gaia_columns(expnum: int, fiberids: np.ndarray) -> dict[str, np.ndarray]:
    path = gaia_path(expnum)
    count = len(fiberids)
    result = {
        "gaia_status": np.full(count, "unavailable", dtype="U16"),
        "gaia_clean_1pct": np.zeros(count, dtype=bool),
        "gaia_ratio": np.full(count, np.nan),
        "gaia_n_sources": np.full(count, -1, dtype=np.int32),
    }
    if path is None:
        return result
    try:
        table = Table.read(path, hdu="FIBERS")
    except KeyError:
        table = Table.read(path, hdu=1)
    by_id = {int(value): index for index, value in enumerate(table["fiberid"])}
    clean = gaia_ratio_clean_mask(table, np.nextafter(0.01, 0.0))
    for out_index, fiberid in enumerate(fiberids):
        source_index = by_id.get(int(fiberid))
        if source_index is None:
            result["gaia_status"][out_index] = "missing_fiber"
            continue
        n_sources = int(table["n_gaia_sources"][source_index])
        ratio = float(table["ratio_gaia_to_lvm_flux_plus_sky_g"][source_index])
        result["gaia_n_sources"][out_index] = n_sources
        result["gaia_ratio"][out_index] = ratio
        result["gaia_clean_1pct"][out_index] = bool(clean[source_index])
        result["gaia_status"][out_index] = (str(table['gaia_status'][source_index])
            if 'gaia_status' in table.colnames else "clean" if clean[source_index] else "contaminated")
    return result


def extract_exposure(item: tuple[int, dict]) -> dict:
    slot, record = item
    path = Path(record["source_path"])
    try:
        with fits.open(path, memmap=True) as hdul:
            comments = hdul[0].header.get("COMMENT", [])
            comments = [comments] if isinstance(comments, str) else comments
            epochs = re.findall(r"/calib/(\d+)/", "\n".join(comments))
            calibration_epoch = max(map(int, epochs)) if epochs else -1
            wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
            scale = _flux_scale(hdul[0].header, hdul["FLUX"])
            slitmap = Table(hdul["SLITMAP"].data)
            chosen = slit_indices(slitmap)
            indices = np.asarray([row[0] for row in chosen])
            flux = np.asarray(hdul["FLUX"].data[indices], dtype=np.float32) * scale
            sky = np.asarray(hdul["SKY"].data[indices], dtype=np.float32) * scale
            ivar = np.asarray(hdul["IVAR"].data[indices], dtype=np.float32) / scale**2
            sky_ivar = np.asarray(hdul["SKY_IVAR"].data[indices], dtype=np.float32) / scale**2
            mask = np.asarray(hdul["MASK"].data[indices], dtype=np.uint8)
            lsf = np.asarray(hdul["LSF"].data[indices], dtype=np.float32)
        combined = flux + sky
        combined_ivar = _combined_ivar(ivar, sky_ivar)
        valid = (
            np.isfinite(combined) & np.isfinite(combined_ivar) & (combined_ivar > 0)
            & np.isfinite(lsf) & (lsf > 0) & (mask == 0)
        )
        combined_ivar[~valid] = 0.0
        valid_fraction = valid.mean(axis=1)
        if np.any(valid_fraction < MIN_VALID_FRACTION):
            raise ValueError(
                f"selected fiber valid fraction below {MIN_VALID_FRACTION}: {valid_fraction.min():.3f}")
        fiber_meta = slitmap[indices].copy()
        fiber_meta.add_column(Column(indices, name="source_fiber_index"), index=0)
        fiber_meta.add_column(Column([row[3] for row in chosen], name="slit_rank"), index=1)
        fiber_meta.add_column(Column([row[2] for row in chosen], name="slit_quantile"), index=2)
        fiber_meta.add_column(Column(valid_fraction, name="fit_valid_fraction"), index=3)
        for name, values in gaia_columns(int(record["expnum"]), np.asarray(fiber_meta["fiberid"])) .items():
            fiber_meta[name] = values
        return dict(slot=slot, status="OK", wave=wave, flux=combined,
                    ivar=combined_ivar, mask=mask, lsf=lsf, fiber_meta=fiber_meta,
                    calibration_epoch=calibration_epoch)
    except Exception as exc:  # one bad exposure must not stop the corpus
        return dict(slot=slot, status="SKIP", error=f"{type(exc).__name__}: {exc}")


def bounded_extract(items, workers: int):
    iterator = iter(items)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = set()
        for _ in range(min(len(items), workers * 2)):
            try:
                pending.add(pool.submit(extract_exposure, next(iterator)))
            except StopIteration:
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                try:
                    pending.add(pool.submit(extract_exposure, next(iterator)))
                except StopIteration:
                    pass


def repeated_meta(row, fiber_meta: Table) -> Table:
    exposure = Table(rows=[tuple(row[name] for name in row.index)], names=list(row.index))
    exposure = vstack([exposure] * len(fiber_meta), metadata_conflicts="silent")
    for column in fiber_meta.colnames:
        if column not in exposure.colnames:
            exposure[column] = fiber_meta[column]
    return exposure


def build(output: Path, workers: int, pilot: bool = False) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    audit = select_exposures()
    audit.to_csv(output / "exposure_selection_audit.csv", index=False)
    selected = audit[audit.selected].copy().sort_values(["full_moon_mjd_utc", "midpoint_mjd_utc"])
    if pilot:
        cycles = selected.full_moon_mjd_utc.drop_duplicates().to_numpy()
        chosen_cycles = cycles[[0, -1]] if len(cycles) > 1 else cycles
        selected = (selected[selected.full_moon_mjd_utc.isin(chosen_cycles)]
                    .sort_values("relative_lunar_continuum", ascending=False)
                    .groupby("full_moon_mjd_utc", as_index=False).head(1)
                    .sort_values("midpoint_mjd_utc"))
    selected.to_csv(output / "selected_exposures.csv", index=False)
    max_rows = 9 * len(selected)
    if not max_rows:
        raise RuntimeError("No selected local exposures")
    with fits.open(STACK, memmap=True) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
    required = max_rows * len(wave) * (4 + 4 + 1 + 4) * 2
    if shutil.disk_usage(output).free < required:
        raise OSError(f"Need approximately {required / 2**30:.1f} GiB free")
    temporary = Path(tempfile.mkdtemp(prefix=".fiber-corpus-", dir=output))
    arrays = {
        name: np.lib.format.open_memmap(temporary / f"{name}.npy", mode="w+", dtype=dtype,
                                        shape=(max_rows, len(wave)))
        for name, dtype in (("FLUX_SCI", "f4"), ("IVAR_SCI", "f4"),
                            ("MASK_SCI", "u1"), ("LSF_SCI", "f4"))
    }
    records = selected.to_dict("records")
    results = [None] * len(records)
    for count, result in enumerate(bounded_extract(list(enumerate(records)), workers), 1):
        results[result["slot"]] = result
        if count % 50 == 0:
            print(f"Extracted {count}/{len(records)} exposures", flush=True)
    status_rows, meta_tables, cursor = [], [], 0
    for slot, record in enumerate(records):
        result = results[slot]
        if result["status"] != "OK":
            status_rows.append({"expnum": record["expnum"], "status": "SKIP", "error": result["error"]})
            continue
        if not np.allclose(result["wave"], wave, rtol=0.0, atol=1e-5):
            status_rows.append({"expnum": record["expnum"], "status": "SKIP", "error": "WAVE mismatch"})
            continue
        stop = cursor + 9
        for name in arrays:
            arrays[name][cursor:stop] = result[name.lower().replace("_sci", "")]
        record = dict(record)
        record["calibration_epoch"] = result.get("calibration_epoch", -1)
        meta_tables.append(repeated_meta(pd.Series(record), result["fiber_meta"]))
        status_rows.append({"expnum": record["expnum"], "status": "OK", "error": ""})
        cursor = stop
    pd.DataFrame(status_rows).to_csv(output / "extraction_status.csv", index=False)
    if not cursor:
        raise RuntimeError("Every selected exposure failed extraction")
    meta = vstack(meta_tables, metadata_conflicts="silent")
    meta["corpus_row"] = np.arange(cursor, dtype=np.int64)
    primary = fits.PrimaryHDU()
    primary.header["DRPVER"] = "1.3.2"
    primary.header["NFIBEXP"] = 9
    primary.header["RELTHR"] = MIN_RELATIVE_CONTINUUM
    primary.header["DAYWIN"] = MAX_DAYS_FROM_FULL
    primary.header["CREATED"] = Time.now().isot
    hdus = [primary, fits.ImageHDU(wave, name="WAVE")]
    hdus.extend(fits.ImageHDU(np.asarray(arrays[name][:cursor]), name=name) for name in arrays)
    hdus.append(fits.BinTableHDU(meta, name="META"))
    destination = output / ("individual_science_fibers_pilot.fits" if pilot else "individual_science_fibers.fits")
    tmp_fits = destination.with_suffix(".tmp.fits")
    fits.HDUList(hdus).writeto(tmp_fits, overwrite=True, checksum=True)
    os.replace(tmp_fits, destination)
    shutil.rmtree(temporary)
    summary = dict(input_stack=str(STACK), output=str(destination), pilot=pilot,
                   selected_exposures=len(selected), extracted_exposures=cursor // 9,
                   rows=cursor, fibers_per_exposure=9, slit_quantiles=QUANTILES,
                   continuum_window_angstrom=CONTINUUM_WINDOW,
                   max_days_from_full_moon=MAX_DAYS_FROM_FULL,
                   minimum_relative_continuum=MIN_RELATIVE_CONTINUUM,
                   minimum_valid_fraction=MIN_VALID_FRACTION)
    (output / ("pilot_summary.json" if pilot else "corpus_summary.json")).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    build(args.output.expanduser().resolve(), args.workers, args.pilot)


if __name__ == "__main__":
    main()
