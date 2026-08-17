#!/usr/bin/env python
"""Build a median spectral stack from lvmSFrame FITS files."""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table


DRP_VERSION = "1.2.1"
DEFAULT_RELATIVE_LIST = f"sdsswork/users/u6058164/lvmsframe_{DRP_VERSION}.dat"
DEFAULT_RELATIVE_OUTPUT = (
    f"sdsswork/users/u6058164/lvmsframe_median_stack_{DRP_VERSION}.fits"
)
OUTPUT_ARRAYS = (
    "FLUX_SCI",
    "FLUX_SKY_NEAR",
    "FLUX_SKY_FAR",
    "FLUX_SCI_NOSKY",
    "LSF_SCI",
    "LSF_SKY_NEAR",
    "LSF_SKY_FAR",
)


REFERENCE_WAVE: np.ndarray | None = None


def sas_path(relative_path: str) -> Path:
    sas_base = os.environ.get("SAS_BASE_DIR")
    if not sas_base:
        raise RuntimeError("SAS_BASE_DIR is not set")
    return Path(sas_base) / relative_path


def default_file_list() -> Path:
    return sas_path(DEFAULT_RELATIVE_LIST)


def percentile_label(value: float) -> str:
    if float(value).is_integer():
        return str(int(value))
    return f"{value:g}".replace(".", "p")


def default_output_path(
    sci_percentile: float,
    sky_percentile: float,
    every_nth: int,
) -> Path:
    base = sas_path(DEFAULT_RELATIVE_OUTPUT)
    every_suffix = f"_every{every_nth}" if every_nth > 1 else ""
    if sci_percentile == sky_percentile:
        suffix = f"p{percentile_label(sci_percentile)}{every_suffix}"
    else:
        suffix = (
            f"psci{percentile_label(sci_percentile)}"
            f"_psky{percentile_label(sky_percentile)}{every_suffix}"
        )
    return base.with_name(f"{base.stem}_{suffix}{base.suffix}")


def init_worker(reference_wave: np.ndarray) -> None:
    global REFERENCE_WAVE
    REFERENCE_WAVE = np.asarray(reference_wave, dtype=np.float32)


def read_file_list(
    path: Path,
    every_nth: int = 1,
    limit: int | None = None,
) -> list[tuple[int, Path]]:
    with path.open("r", encoding="utf-8") as handle:
        paths = [Path(line.strip()) for line in handle if line.strip()]
    indexed_paths = list(enumerate(paths))[::every_nth]
    if limit is not None:
        return indexed_paths[:limit]
    return indexed_paths


def first_valid_wave(paths: list[Path]) -> tuple[np.ndarray, Path]:
    errors = []
    for path in paths:
        try:
            with fits.open(path, memmap=True) as hdul:
                wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
                if wave.ndim != 1:
                    raise ValueError(f"WAVE is not one-dimensional: {wave.shape}")
                return wave.copy(), path
        except Exception as exc:  # pragma: no cover - depends on external data
            errors.append(f"{path}: {exc}")
    joined = "\n".join(errors[:10])
    raise RuntimeError(f"could not read WAVE from any input file:\n{joined}")


def finite_float(value: Any, default: float = np.nan) -> float:
    if value is None:
        return default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    if math.isfinite(result):
        return result
    return default


def finite_int(value: Any, default: int = -1) -> int:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return default


def expnum_from_path(path: Path) -> int:
    match = re.search(r"lvmSFrame-(\d+)\.fits$", path.name)
    if not match:
        return -1
    return finite_int(match.group(1))


def expnum_value(path: Path, header: fits.Header) -> int:
    for key in ("EXPNUM", "EXPOSURE"):
        value = finite_int(header.get(key))
        if value >= 0:
            return value
    return expnum_from_path(path)


def text_value(value: Any, default: str = "") -> str:
    if value is None:
        return default
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    return str(value).strip()


def header_float(header: fits.Header, keys: tuple[str, ...]) -> float:
    for key in keys:
        value = finite_float(header.get(key))
        if np.isfinite(value):
            return value
    return np.nan


def header_text(header: fits.Header, keys: tuple[str, ...]) -> str:
    for key in keys:
        value = text_value(header.get(key))
        if value:
            return value
    return ""


def coord_pair(
    header: fits.Header,
    primary_keys: tuple[str, str],
    fallback_keys: tuple[str, str],
) -> tuple[float, float]:
    ra = header_float(header, (primary_keys[0],))
    dec = header_float(header, (primary_keys[1],))
    if np.isfinite(ra) and np.isfinite(dec):
        return ra, dec
    return (
        header_float(header, (fallback_keys[0],)),
        header_float(header, (fallback_keys[1],)),
    )


def angular_sep_deg(
    ra1: float,
    dec1: float,
    ra2: float,
    dec2: float,
) -> float:
    values = (ra1, dec1, ra2, dec2)
    if not all(np.isfinite(values)):
        return np.nan
    ra1_rad, dec1_rad, ra2_rad, dec2_rad = np.deg2rad(values)
    sin_d_dec = np.sin((dec2_rad - dec1_rad) / 2.0)
    sin_d_ra = np.sin((ra2_rad - ra1_rad) / 2.0)
    a = sin_d_dec**2 + np.cos(dec1_rad) * np.cos(dec2_rad) * sin_d_ra**2
    return float(np.rad2deg(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))))


def string_column(values: Any) -> np.ndarray:
    return np.char.strip(np.asarray(values).astype(str))


def selected_indices(
    flux_plus_sky: np.ndarray,
    group_mask: np.ndarray,
    faint_fiber_percentile: float,
) -> tuple[np.ndarray, int, int]:
    good_indices = np.flatnonzero(group_mask)
    if good_indices.size == 0:
        return good_indices, 0, 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        brightness = np.nanmedian(flux_plus_sky[good_indices], axis=1)

    finite = np.isfinite(brightness)
    finite_indices = good_indices[finite]
    finite_brightness = brightness[finite]
    if finite_indices.size == 0:
        return finite_indices, int(good_indices.size), 0

    keep_count = max(
        1,
        int(np.ceil(faint_fiber_percentile / 100.0 * finite_indices.size)),
    )
    order = np.argsort(finite_brightness, kind="stable")
    selected = finite_indices[order[:keep_count]]
    return selected, int(good_indices.size), int(selected.size)


def median_or_nan(data: np.ndarray, indices: np.ndarray, n_wave: int) -> np.ndarray:
    if indices.size == 0:
        return np.full(n_wave, np.nan, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = np.nanmedian(data[indices], axis=0)
    return np.asarray(result, dtype=np.float32)


def lsf_or_nan(lsf: np.ndarray, indices: np.ndarray, n_wave: int) -> np.ndarray:
    if indices.size == 0:
        return np.full(n_wave, np.nan, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = sigma_clipped_stats(
            lsf[indices],
            axis=0,
            sigma=3.0,
            maxiters=5,
        )[0]
    return np.asarray(result, dtype=np.float32)


def empty_meta(
    path: Path,
    status: str,
    error: str = "",
    *,
    input_index: int = -1,
    sci_percentile: float = np.nan,
    sky_percentile: float = np.nan,
    every_nth: int = 1,
) -> dict[str, Any]:
    return {
        "path": str(path),
        "input_index": input_index,
        "sci_faint_fiber_percentile": sci_percentile,
        "sky_faint_fiber_percentile": sky_percentile,
        "every_nth": every_nth,
        "exposure": -1,
        "expnum": -1,
        "mjd": -1,
        "tile_id": -1,
        "obstime": "",
        "date_obs": "",
        "sci_ra": np.nan,
        "sci_dec": np.nan,
        "skye_ra": np.nan,
        "skye_dec": np.nan,
        "skyw_ra": np.nan,
        "skyw_dec": np.nan,
        "sky_near_label": "",
        "sky_far_label": "",
        "sky_near_ra": np.nan,
        "sky_near_dec": np.nan,
        "sky_far_ra": np.nan,
        "sky_far_dec": np.nan,
        "skye_sep_deg": np.nan,
        "skyw_sep_deg": np.nan,
        "sky_near_sep_deg": np.nan,
        "sky_far_sep_deg": np.nan,
        "sci_alt": np.nan,
        "skye_alt": np.nan,
        "skyw_alt": np.nan,
        "sci_airmass": np.nan,
        "skye_airmass": np.nan,
        "skyw_airmass": np.nan,
        "sci_moon_sep": np.nan,
        "skye_moon_sep": np.nan,
        "skyw_moon_sep": np.nan,
        "moon_alt": np.nan,
        "sun_alt": np.nan,
        "moon_ra": np.nan,
        "moon_dec": np.nan,
        "moon_phase": np.nan,
        "moon_fli": np.nan,
        "moon_illum": np.nan,
        "fibers_sci_good": 0,
        "fibers_sci_used": 0,
        "fibers_skye_good": 0,
        "fibers_skye_used": 0,
        "fibers_skyw_good": 0,
        "fibers_skyw_used": 0,
        "fibers_sky_near_good": 0,
        "fibers_sky_near_used": 0,
        "fibers_sky_far_good": 0,
        "fibers_sky_far_used": 0,
        "fibers_sci": 0,
        "fiberfrac_sci": np.nan,
        "fibers_sky_near": 0,
        "fiberfrac_sky_near": np.nan,
        "fibers_sky_far": 0,
        "fiberfrac_sky_far": np.nan,
        "worker_pid": os.getpid(),
        "status": status,
        "error": error[:512],
    }


def extract_meta(
    path: Path,
    input_index: int,
    sci_percentile: float,
    sky_percentile: float,
    every_nth: int,
    header: fits.Header,
    sky_labels: dict[str, str],
    coords: dict[str, tuple[float, float]],
    separations: dict[str, float],
    counts: dict[str, int],
    status: str = "OK",
    error: str = "",
) -> dict[str, Any]:
    meta = empty_meta(
        path,
        status=status,
        error=error,
        input_index=input_index,
        sci_percentile=sci_percentile,
        sky_percentile=sky_percentile,
        every_nth=every_nth,
    )
    near_label = sky_labels["near"]
    far_label = sky_labels["far"]

    def used_count(label: str) -> int:
        return int(counts[f"{label}_used"])

    def used_fraction(label: str) -> float:
        total = int(counts[f"{label}_good"])
        if total <= 0:
            return np.nan
        return float(counts[f"{label}_used"] / total)

    meta.update(
        {
            "exposure": finite_int(header.get("EXPOSURE")),
            "expnum": expnum_value(path, header),
            "mjd": finite_int(header.get("MJD")),
            "tile_id": finite_int(header.get("TILE_ID")),
            "obstime": header_text(header, ("OBSTIME", "DATE-OBS")),
            "date_obs": header_text(header, ("DATE-OBS", "OBSTIME")),
            "sci_ra": coords["Sci"][0],
            "sci_dec": coords["Sci"][1],
            "skye_ra": coords["SkyE"][0],
            "skye_dec": coords["SkyE"][1],
            "skyw_ra": coords["SkyW"][0],
            "skyw_dec": coords["SkyW"][1],
            "sky_near_label": sky_labels["near"],
            "sky_far_label": sky_labels["far"],
            "sky_near_ra": coords[sky_labels["near"]][0],
            "sky_near_dec": coords[sky_labels["near"]][1],
            "sky_far_ra": coords[sky_labels["far"]][0],
            "sky_far_dec": coords[sky_labels["far"]][1],
            "skye_sep_deg": separations["SkyE"],
            "skyw_sep_deg": separations["SkyW"],
            "sky_near_sep_deg": separations[sky_labels["near"]],
            "sky_far_sep_deg": separations[sky_labels["far"]],
            "sci_alt": header_float(header, ("SKY SCI_ALT", "SCIALT", "ALT")),
            "skye_alt": header_float(header, ("SKY SKYE_ALT", "SKYEALT")),
            "skyw_alt": header_float(header, ("SKY SKYW_ALT", "SKYWALT")),
            "sci_airmass": header_float(header, ("SCIAM", "TESCIAM", "AIRMASS")),
            "skye_airmass": header_float(header, ("SKYEAM", "TESKYEAM")),
            "skyw_airmass": header_float(header, ("SKYWAM", "TESKYWAM")),
            "sci_moon_sep": header_float(header, ("SKY SCI_MOON_SEP",)),
            "skye_moon_sep": header_float(header, ("SKY SKYE_MOON_SEP",)),
            "skyw_moon_sep": header_float(header, ("SKY SKYW_MOON_SEP",)),
            "moon_alt": header_float(header, ("SKY MOON_ALT", "MOONALT")),
            "sun_alt": header_float(header, ("SKY SUN_ALT", "SUNALT")),
            "moon_ra": header_float(header, ("SKY MOON_RA", "MOONRA")),
            "moon_dec": header_float(header, ("SKY MOON_DEC", "MOONDEC")),
            "moon_phase": header_float(header, ("SKY MOON_PHASE", "MOONPHA")),
            "moon_fli": header_float(header, ("SKY MOON_FLI",)),
            "moon_illum": header_float(header, ("SKY MOON_FLI", "MOONILL")),
            "fibers_sci_good": counts["Sci_good"],
            "fibers_sci_used": counts["Sci_used"],
            "fibers_skye_good": counts["SkyE_good"],
            "fibers_skye_used": counts["SkyE_used"],
            "fibers_skyw_good": counts["SkyW_good"],
            "fibers_skyw_used": counts["SkyW_used"],
            "fibers_sky_near_good": counts[f"{near_label}_good"],
            "fibers_sky_near_used": counts[f"{near_label}_used"],
            "fibers_sky_far_good": counts[f"{far_label}_good"],
            "fibers_sky_far_used": counts[f"{far_label}_used"],
            "fibers_sci": used_count("Sci"),
            "fiberfrac_sci": used_fraction("Sci"),
            "fibers_sky_near": used_count(near_label),
            "fiberfrac_sky_near": used_fraction(near_label),
            "fibers_sky_far": used_count(far_label),
            "fiberfrac_sky_far": used_fraction(far_label),
        }
    )
    return meta


def sky_geometry(header: fits.Header) -> tuple[
    dict[str, tuple[float, float]],
    dict[str, float],
    dict[str, str],
]:
    coords = {
        "Sci": coord_pair(header, ("SCIRA", "SCIDEC"), ("POSCIRA", "POSCIDE")),
        "SkyE": coord_pair(header, ("SKYERA", "SKYEDEC"), ("POSKYERA", "POSKYEDE")),
        "SkyW": coord_pair(header, ("SKYWRA", "SKYWDEC"), ("POSKYWRA", "POSKYWDE")),
    }
    separations = {
        "SkyE": header_float(header, ("SKY SCI_SKYE_SEP", "SCI_SKYE_SEP")),
        "SkyW": header_float(header, ("SKY SCI_SKYW_SEP", "SCI_SKYW_SEP")),
    }
    if not np.isfinite(separations["SkyE"]):
        separations["SkyE"] = angular_sep_deg(*coords["Sci"], *coords["SkyE"])
    if not np.isfinite(separations["SkyW"]):
        separations["SkyW"] = angular_sep_deg(*coords["Sci"], *coords["SkyW"])

    if np.isfinite(separations["SkyE"]) and np.isfinite(separations["SkyW"]):
        near = "SkyE" if separations["SkyE"] <= separations["SkyW"] else "SkyW"
    else:
        near = "SkyE"
    far = "SkyW" if near == "SkyE" else "SkyE"
    return coords, separations, {"near": near, "far": far}


def process_one(
    index: int,
    input_index: int,
    path: Path,
    sci_percentile: float,
    sky_percentile: float,
    every_nth: int,
) -> dict[str, Any]:
    reference_wave = REFERENCE_WAVE
    if reference_wave is None:
        raise RuntimeError("worker reference wavelength grid was not initialized")
    n_wave = int(reference_wave.size)

    try:
        with fits.open(path, memmap=True) as hdul:
            header = hdul[0].header
            wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
            if wave.shape != reference_wave.shape or not np.allclose(
                wave,
                reference_wave,
                rtol=0.0,
                atol=1.0e-5,
                equal_nan=True,
            ):
                raise ValueError("WAVE grid does not match the reference file")

            flux = np.asarray(hdul["FLUX"].data, dtype=np.float32)
            sky = np.asarray(hdul["SKY"].data, dtype=np.float32)
            lsf = np.asarray(hdul["LSF"].data, dtype=np.float32)
            if flux.shape != sky.shape or flux.shape != lsf.shape:
                raise ValueError(
                    f"FLUX/SKY/LSF shape mismatch: {flux.shape}, {sky.shape}, {lsf.shape}"
                )
            if flux.ndim != 2 or flux.shape[1] != n_wave:
                raise ValueError(f"unexpected FLUX shape: {flux.shape}")

            slitmap = Table.read(hdul, hdu="SLITMAP")
            if len(slitmap) != flux.shape[0]:
                raise ValueError(
                    f"SLITMAP rows ({len(slitmap)}) do not match fibers ({flux.shape[0]})"
                )

            telescopes = string_column(slitmap["telescope"])
            fibstatus = np.asarray(slitmap["fibstatus"])
            good = fibstatus == 0
            flux_plus_sky = flux + sky

            selected: dict[str, np.ndarray] = {}
            counts: dict[str, int] = {}
            for telescope in ("Sci", "SkyE", "SkyW"):
                pct = sci_percentile if telescope == "Sci" else sky_percentile
                indices, n_good, n_used = selected_indices(
                    flux_plus_sky,
                    good & (telescopes == telescope),
                    pct,
                )
                selected[telescope] = indices
                counts[f"{telescope}_good"] = n_good
                counts[f"{telescope}_used"] = n_used

            coords, separations, sky_labels = sky_geometry(header)
            near_label = sky_labels["near"]
            far_label = sky_labels["far"]

            arrays = {
                "FLUX_SCI": median_or_nan(flux_plus_sky, selected["Sci"], n_wave),
                "FLUX_SKY_NEAR": median_or_nan(
                    flux_plus_sky,
                    selected[near_label],
                    n_wave,
                ),
                "FLUX_SKY_FAR": median_or_nan(
                    flux_plus_sky,
                    selected[far_label],
                    n_wave,
                ),
                "FLUX_SCI_NOSKY": median_or_nan(flux, selected["Sci"], n_wave),
                "LSF_SCI": lsf_or_nan(lsf, selected["Sci"], n_wave),
                "LSF_SKY_NEAR": lsf_or_nan(lsf, selected[near_label], n_wave),
                "LSF_SKY_FAR": lsf_or_nan(lsf, selected[far_label], n_wave),
            }
            meta = extract_meta(
                path,
                input_index,
                sci_percentile,
                sky_percentile,
                every_nth,
                header,
                sky_labels,
                coords,
                separations,
                counts,
            )
            return {"index": index, "path": str(path), "arrays": arrays, "meta": meta}
    except Exception as exc:
        arrays = {
            name: np.full(n_wave, np.nan, dtype=np.float32) for name in OUTPUT_ARRAYS
        }
        meta = empty_meta(
            path,
            status="ERROR",
            error=repr(exc),
            input_index=input_index,
            sci_percentile=sci_percentile,
            sky_percentile=sky_percentile,
            every_nth=every_nth,
        )
        return {"index": index, "path": str(path), "arrays": arrays, "meta": meta}


def table_from_meta(rows: list[dict[str, Any]]) -> Table:
    names = list(rows[0].keys())
    columns = {name: [row[name] for row in rows] for name in names}
    table = Table(columns)
    for name in ("path", "error"):
        table[name] = np.asarray(table[name], dtype="U512")
    for name in ("status",):
        table[name] = np.asarray(table[name], dtype="U16")
    for name in ("obstime", "date_obs"):
        table[name] = np.asarray(table[name], dtype="U32")
    for name in ("sky_near_label", "sky_far_label"):
        table[name] = np.asarray(table[name], dtype="U8")
    return table


def create_memmaps(
    tmpdir: Path,
    n_files: int,
    n_wave: int,
) -> dict[str, np.memmap]:
    tmpdir.mkdir(parents=True, exist_ok=True)
    return {
        name: np.lib.format.open_memmap(
            tmpdir / f"{name.lower()}.npy",
            mode="w+",
            dtype=np.float32,
            shape=(n_files, n_wave),
        )
        for name in OUTPUT_ARRAYS
    }


def write_output(
    output_path: Path,
    wave: np.ndarray,
    arrays: dict[str, np.memmap],
    meta_rows: list[dict[str, Any]],
    input_list: Path,
    reference_path: Path,
    sci_percentile: float,
    sky_percentile: float,
    every_nth: int,
    overwrite: bool,
) -> None:
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"output exists, pass --overwrite to replace: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    primary_header = fits.Header()
    primary_header["CREATED"] = datetime.now(timezone.utc).isoformat()
    primary_header["DRPVER"] = DRP_VERSION
    primary_header["NFILES"] = len(meta_rows)
    primary_header["FAINTSCI"] = sci_percentile
    primary_header["FAINTSKY"] = sky_percentile
    primary_header["EVERYNTH"] = every_nth
    primary_header["INLIST"] = str(input_list)
    primary_header["REFWAVE"] = str(reference_path)
    primary_header["COMMENT"] = "Median SFrame stack built from FLUX, SKY, LSF, and SLITMAP."

    hdus: list[fits.hdu.base.ExtensionHDU | fits.PrimaryHDU] = [
        fits.PrimaryHDU(header=primary_header),
        fits.ImageHDU(data=np.asarray(wave, dtype=np.float32), name="WAVE"),
    ]
    for name in OUTPUT_ARRAYS:
        hdus.append(fits.ImageHDU(data=arrays[name], name=name))
    hdus.append(fits.BinTableHDU(data=table_from_meta(meta_rows), name="META"))
    fits.HDUList(hdus).writeto(output_path, overwrite=overwrite, checksum=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a median spectra stack from lvmSFrame FITS files."
    )
    parser.add_argument(
        "--input-list",
        type=Path,
        default=None,
        help=f"Input .dat list. Default: $SAS_BASE_DIR/{DEFAULT_RELATIVE_LIST}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output FITS. Default name encodes the percentile(s) and appends "
            "_every{N} only when N > 1."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of worker processes.",
    )
    parser.add_argument(
        "--faint-fiber-percentile",
        type=float,
        default=70.0,
        help=(
            "Fallback faintest-fiber percentile applied to any telescope group "
            "without a per-group override. Fibers are ranked by median "
            "FLUX+SKY and the faintest fraction is kept. Default: 70."
        ),
    )
    parser.add_argument(
        "--sci-faint-fiber-percentile",
        type=float,
        default=None,
        help=(
            "Overrides --faint-fiber-percentile for the science telescope. "
            "Default: fall back to --faint-fiber-percentile."
        ),
    )
    parser.add_argument(
        "--sky-faint-fiber-percentile",
        type=float,
        default=None,
        help=(
            "Overrides --faint-fiber-percentile for both sky telescopes. "
            "Default: fall back to --faint-fiber-percentile."
        ),
    )
    parser.add_argument(
        "--every-nth",
        type=int,
        default=1,
        help="Process every Nth input-list entry, starting with the first. Default: 1.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N files from the list.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output FITS file.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary memmap arrays after a successful write.",
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        default=None,
        help="Directory for temporary memmap arrays. Default: output directory.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print progress after this many completed files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_list = args.input_list if args.input_list is not None else default_file_list()
    sci_percentile = (
        args.sci_faint_fiber_percentile
        if args.sci_faint_fiber_percentile is not None
        else args.faint_fiber_percentile
    )
    sky_percentile = (
        args.sky_faint_fiber_percentile
        if args.sky_faint_fiber_percentile is not None
        else args.faint_fiber_percentile
    )
    output_path = (
        args.output
        if args.output is not None
        else default_output_path(sci_percentile, sky_percentile, args.every_nth)
    )
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive")
    if args.workers <= 0:
        raise ValueError("--workers must be positive")
    for label, value in (("science", sci_percentile), ("sky", sky_percentile)):
        if not 0.0 < value <= 100.0:
            raise ValueError(
                f"{label}-arm faint-fiber percentile must be in (0, 100]; got {value}"
            )
    if args.every_nth <= 0:
        raise ValueError("--every-nth must be positive")

    indexed_paths = read_file_list(
        input_list, every_nth=args.every_nth, limit=args.limit
    )
    if not indexed_paths:
        raise RuntimeError(f"no input paths found in {input_list}")
    paths = [path for _, path in indexed_paths]

    wave, reference_path = first_valid_wave(paths)
    n_files = len(paths)
    n_wave = int(wave.size)
    base_temp_dir = args.temp_dir if args.temp_dir is not None else output_path.parent
    tmpdir = Path(
        tempfile.mkdtemp(
            prefix=f"{output_path.stem}_tmp_",
            dir=base_temp_dir,
        )
    )

    arrays = create_memmaps(tmpdir, n_files, n_wave)
    meta_rows: list[dict[str, Any] | None] = [None] * n_files
    workers = min(args.workers, n_files)
    completed = 0
    print(
        f"processing {n_files} files with {workers} workers; "
        f"output rows have {n_wave} wavelength pixels"
    )
    print(f"temporary arrays: {tmpdir}")

    try:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=init_worker,
            initargs=(wave,),
        ) as executor:
            futures = [
                executor.submit(
                    process_one,
                    index,
                    input_index,
                    path,
                    sci_percentile,
                    sky_percentile,
                    args.every_nth,
                )
                for index, (input_index, path) in enumerate(indexed_paths)
            ]
            for future in as_completed(futures):
                result = future.result()
                row_index = int(result["index"])
                for name in OUTPUT_ARRAYS:
                    arrays[name][row_index, :] = result["arrays"][name]
                meta_rows[row_index] = result["meta"]
                completed += 1
                if (
                    args.progress_every > 0
                    and (completed % args.progress_every == 0 or completed == n_files)
                ):
                    print(f"completed {completed}/{n_files}")

        missing = [index for index, row in enumerate(meta_rows) if row is None]
        if missing:
            raise RuntimeError(f"missing metadata for rows: {missing[:10]}")
        for memmap in arrays.values():
            memmap.flush()
        write_output(
            output_path=output_path,
            wave=wave,
            arrays=arrays,
            meta_rows=[row for row in meta_rows if row is not None],
            input_list=input_list,
            reference_path=reference_path,
            sci_percentile=sci_percentile,
            sky_percentile=sky_percentile,
            every_nth=args.every_nth,
            overwrite=args.overwrite,
        )
        print(f"wrote {output_path}")
    finally:
        if args.keep_temp:
            print(f"kept temporary arrays in {tmpdir}")
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
