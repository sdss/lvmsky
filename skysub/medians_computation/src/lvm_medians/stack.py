"""SFrame discovery and offline stack builders."""

from __future__ import annotations

import logging
import math
import os
import re
import shutil
import tempfile
import warnings
from collections.abc import Callable, Iterator
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Column, Table, vstack

from .metadata import median_meta

TELESCOPES = ("Sci", "SkyE", "SkyW")
OUTPUT_ARRAYS = (
    "FLUX_SCI",
    "FLUX_SKY_NEAR",
    "FLUX_SKY_FAR",
    "FLUX_SCI_NOSKY",
    "LSF_SCI",
    "LSF_SKY_NEAR",
    "LSF_SKY_FAR",
)
FLUX_UNIT = u.erg / (u.Angstrom * u.s * u.cm**2)
GAIA_VALID_FRACTION_MIN = 0.95
REFERENCE_WAVE: np.ndarray | None = None


class SkipExposure(Exception):
    """An input is valid but cannot contribute to the requested stack."""


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as fh:
        fh.write(text)
        temporary = Path(fh.name)
    os.replace(temporary, path)


def scan_sframes(root: Path, progress: Callable[[int], None] | None = None) -> list[Path]:
    if not root.is_dir():
        raise NotADirectoryError(f"SFrame root is not a directory: {root}")
    paths: list[Path] = []
    for directory, directories, filenames in os.walk(root):
        directories.sort()
        for filename in sorted(filenames):
            if filename.startswith("lvmSFrame-") and filename.endswith((".fits", ".fits.gz")):
                paths.append(Path(directory) / filename)
                if progress:
                    progress(len(paths))
    return paths


def write_manifest(paths: list[Path], destination: Path) -> None:
    atomic_text(destination, "".join(f"{path.resolve()}\n" for path in paths))


def expnum_from_path(path: Path) -> int:
    match = re.search(r"lvmSFrame-(\d+)\.fits(?:\.gz)?$", path.name)
    return int(match.group(1)) if match else -1


def read_manifest(
    path: Path, every_nth: int = 1, limit: int | None = None
) -> list[tuple[int, Path]]:
    if every_nth < 1:
        raise ValueError("--every-nth must be positive")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be positive")
    base = path.resolve().parent
    entries: list[Path] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        item = Path(value).expanduser()
        entries.append(item if item.is_absolute() else base / item)
    selected = list(enumerate(entries))[::every_nth]
    if limit is not None:
        selected = selected[:limit]
    if not selected:
        raise ValueError(f"SFrame list contains no selected paths: {path}")
    missing = [item for _, item in selected if not item.is_file()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} SFrame-list inputs do not exist; first: {missing[0]}"
        )
    duplicates = len({item.resolve() for _, item in selected}) != len(selected)
    if duplicates:
        raise ValueError("SFrame list contains duplicate selected paths")
    expnums = [expnum_from_path(item) for _, item in selected]
    known = [value for value in expnums if value >= 0]
    if len(known) != len(set(known)):
        raise ValueError("SFrame list contains duplicate exposure numbers")
    return selected


def first_wave(paths: list[Path]) -> tuple[np.ndarray, Path, str]:
    errors: list[str] = []
    for path in paths:
        try:
            with fits.open(path, memmap=True) as hdul:
                wave = np.asarray(hdul["WAVE"].data, dtype=np.float32)
                if wave.ndim != 1:
                    raise ValueError(f"WAVE has shape {wave.shape}")
                return wave.copy(), path, str(hdul[0].header.get("DRPVER", "unknown"))
        except Exception as exc:  # noqa: BLE001 - continue to the next external FITS
            errors.append(f"{path}: {exc}")
    raise RuntimeError("no readable WAVE array; " + "; ".join(errors[:3]))


def _init_worker(wave: np.ndarray) -> None:
    global REFERENCE_WAVE
    REFERENCE_WAVE = np.asarray(wave, dtype=np.float32)


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode(errors="replace").strip()
    return "" if value is None else str(value).strip()


def _strings(values: Any) -> np.ndarray:
    return np.char.strip(np.asarray(values).astype(str))


def _finite_float(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return np.nan
    return number if math.isfinite(number) else np.nan


def _finite_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return -1


def _expnum(path: Path, header: fits.Header) -> int:
    for key in ("EXPNUM", "EXPOSURE"):
        value = _finite_int(header.get(key))
        if value >= 0:
            return value
    return expnum_from_path(path)


def _flux_scale(primary: fits.Header, image: fits.ImageHDU) -> float:
    raw = image.header.get("BUNIT", primary.get("BUNIT"))
    if raw is None:
        raise ValueError("FLUX BUNIT is missing")
    try:
        return float(u.Unit(raw).to(FLUX_UNIT))
    except (ValueError, TypeError) as exc:
        raise ValueError(f"FLUX BUNIT is not convertible to {FLUX_UNIT}: {raw!r}") from exc


def _validate_wave(wave: np.ndarray) -> None:
    if REFERENCE_WAVE is None:
        raise RuntimeError("worker wavelength grid is not initialized")
    if wave.shape != REFERENCE_WAVE.shape or not np.allclose(
        wave, REFERENCE_WAVE, rtol=0.0, atol=1e-5, equal_nan=True
    ):
        raise ValueError("WAVE grid differs from the reference input")


def selected_indices(
    flux_plus_sky: np.ndarray, group_mask: np.ndarray, percentile: float
) -> tuple[np.ndarray, int, int]:
    candidates = np.flatnonzero(group_mask)
    if candidates.size == 0:
        return candidates, 0, 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        brightness = np.nanmedian(flux_plus_sky[candidates], axis=1)
    finite = np.isfinite(brightness)
    candidates, brightness = candidates[finite], brightness[finite]
    count = max(1, int(np.ceil(percentile / 100.0 * candidates.size))) if candidates.size else 0
    selected = candidates[np.argsort(brightness, kind="stable")[:count]]
    return selected, int(group_mask.sum()), int(selected.size)


def gaia_clean_mask(data: Any, sigma: int) -> np.ndarray:
    n_sources = np.asarray(data["n_gaia_sources"], dtype=np.int64)
    n_valid = np.asarray(data["gaia_g_n_valid_sources"], dtype=np.int64)
    gaia_flux = np.asarray(data["gaia_g_flux"], dtype=float)
    error = np.asarray(data["lvm_flux_g_err"], dtype=float)
    fraction = np.asarray(data["lvm_flux_g_valid_fraction"], dtype=float)
    measured = np.isfinite(error) & (error > 0) & (fraction >= GAIA_VALID_FRACTION_MIN)
    complete = (n_valid == n_sources) & (n_valid > 0) & np.isfinite(gaia_flux) & (gaia_flux > 0)
    return measured & ((n_sources == 0) | (complete & (gaia_flux < sigma * error)))


def gaia_ratio_clean_mask(data: Any, threshold: float) -> np.ndarray:
    n_sources = np.asarray(data["n_gaia_sources"], dtype=np.int64)
    n_valid = np.asarray(data["gaia_g_n_valid_sources"], dtype=np.int64)
    observed = np.asarray(data["lvm_flux_plus_sky_g"], dtype=float)
    fraction = np.asarray(data["lvm_flux_plus_sky_g_valid_fraction"], dtype=float)
    ratio = np.asarray(data["ratio_gaia_to_lvm_flux_plus_sky_g"], dtype=float)
    measured = np.isfinite(observed) & (observed > 0) & (fraction >= GAIA_VALID_FRACTION_MIN)
    complete = (n_valid == n_sources) & (n_valid > 0) & np.isfinite(ratio)
    return measured & ((n_sources == 0) | (complete & (ratio <= threshold)))


def _read_gaia_cache(path: Path) -> Any:
    try:
        return fits.getdata(path, "FIBERS")
    except KeyError:
        return fits.getdata(path, 1)


def _median(data: np.ndarray, indices: np.ndarray) -> np.ndarray:
    if not indices.size:
        return np.full(data.shape[1], np.nan, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.asarray(np.nanmedian(data[indices], axis=0), dtype=np.float32)


def _lsf(data: np.ndarray, indices: np.ndarray) -> np.ndarray:
    if not indices.size:
        return np.full(data.shape[1], np.nan, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.asarray(
            sigma_clipped_stats(data[indices], axis=0, sigma=3.0, maxiters=5)[0],
            dtype=np.float32,
        )


def _coords(header: fits.Header, prefix: str) -> tuple[float, float]:
    choices = {
        "Sci": (("SCIRA", "SCIDEC"), ("POSCIRA", "POSCIDE")),
        "SkyE": (("SKYERA", "SKYEDEC"), ("POSKYERA", "POSKYEDE")),
        "SkyW": (("SKYWRA", "SKYWDEC"), ("POSKYWRA", "POSKYWDE")),
    }
    for ra_key, dec_key in choices[prefix]:
        ra, dec = _finite_float(header.get(ra_key)), _finite_float(header.get(dec_key))
        if np.isfinite(ra) and np.isfinite(dec):
            return ra, dec
    return np.nan, np.nan


def _separation(a: tuple[float, float], b: tuple[float, float]) -> float:
    if not np.all(np.isfinite((*a, *b))):
        return np.nan
    ra1, dec1, ra2, dec2 = np.deg2rad((*a, *b))
    value = np.sin((dec2 - dec1) / 2) ** 2
    value += np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2) ** 2
    return float(np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(value, 0, 1)))))


def _median_worker(
    slot: int,
    input_index: int,
    path: Path,
    sci_percentile: float,
    sky_percentile: float,
    gaia_dir: Path | None,
    gaia_sigma: int | None,
    gaia_ratio_threshold: float | None,
) -> dict[str, Any]:
    try:
        with fits.open(path, memmap=True) as hdul:
            header = hdul[0].header
            _validate_wave(np.asarray(hdul["WAVE"].data, dtype=np.float32))
            scale = _flux_scale(header, hdul["FLUX"])
            flux = np.asarray(hdul["FLUX"].data, dtype=np.float32) * scale
            sky = np.asarray(hdul["SKY"].data, dtype=np.float32) * scale
            lsf = np.asarray(hdul["LSF"].data, dtype=np.float32)
            if flux.ndim != 2 or flux.shape != sky.shape or flux.shape != lsf.shape:
                raise ValueError(
                    f"FLUX/SKY/LSF shape mismatch: {flux.shape}/{sky.shape}/{lsf.shape}"
                )
            slitmap = Table.read(hdul, hdu="SLITMAP")
            if len(slitmap) != flux.shape[0]:
                raise ValueError("SLITMAP row count differs from spectral arrays")
            telescopes = _strings(slitmap["telescope"])
            good = np.asarray(slitmap["fibstatus"]) == 0
            combined = flux + sky
            gaia_selection = gaia_sigma is not None or gaia_ratio_threshold is not None
            gaia_clean = np.ones(len(slitmap), dtype=bool)
            expnum = _expnum(path, header)
            if gaia_selection:
                if gaia_dir is None:
                    raise ValueError("Gaia fiber cache directory is required")
                cache = gaia_dir / f"lvmGAIA-fibers-{expnum:08d}.fits"
                table = _read_gaia_cache(cache)
                fiberids = np.asarray(table["fiberid"], dtype=np.int64)
                if len(np.unique(fiberids)) != len(fiberids):
                    raise ValueError(f"duplicate fiberid in {cache}")
                sframe_fiberids = np.asarray(slitmap["fiberid"], dtype=np.int64)
                if len(fiberids) != len(sframe_fiberids) or not np.array_equal(
                    np.sort(fiberids), np.sort(sframe_fiberids)
                ):
                    raise ValueError(f"incomplete or mismatched fiberid coverage in {cache}")
                if "expnum" not in (table.dtype.names or ()) or not np.all(
                    np.asarray(table["expnum"], dtype=np.int64) == expnum
                ):
                    raise ValueError(f"expnum does not match {expnum} in {cache}")
                clean = (
                    gaia_clean_mask(table, gaia_sigma)
                    if gaia_sigma is not None
                    else gaia_ratio_clean_mask(table, float(gaia_ratio_threshold))
                )
                gaia_clean = np.isin(
                    np.asarray(slitmap["fiberid"], dtype=np.int64), fiberids[clean]
                )

            selected: dict[str, np.ndarray] = {}
            counts: dict[str, int] = {}
            for telescope in TELESCOPES:
                group = good & (telescopes == telescope)
                if gaia_selection:
                    indices = np.flatnonzero(group & gaia_clean)
                    selected[telescope] = indices
                    counts[f"{telescope}_good"] = int(group.sum())
                    counts[f"{telescope}_used"] = int(indices.size)
                else:
                    percentile = sci_percentile if telescope == "Sci" else sky_percentile
                    indices, total, used = selected_indices(combined, group, percentile)
                    selected[telescope] = indices
                    counts[f"{telescope}_good"], counts[f"{telescope}_used"] = total, used
            if gaia_selection and any(not selected[name].size for name in TELESCOPES):
                raise SkipExposure("at least one telescope has no Gaia-clean fibers")

            positions = {name: _coords(header, name) for name in TELESCOPES}
            separations = {
                name: _separation(positions["Sci"], positions[name]) for name in ("SkyE", "SkyW")
            }
            near = min(
                separations,
                key=lambda key: separations[key] if np.isfinite(separations[key]) else np.inf,
            )
            far = "SkyW" if near == "SkyE" else "SkyE"
            arrays = {
                "FLUX_SCI": _median(combined, selected["Sci"]),
                "FLUX_SKY_NEAR": _median(combined, selected[near]),
                "FLUX_SKY_FAR": _median(combined, selected[far]),
                "FLUX_SCI_NOSKY": _median(flux, selected["Sci"]),
                "LSF_SCI": _lsf(lsf, selected["Sci"]),
                "LSF_SKY_NEAR": _lsf(lsf, selected[near]),
                "LSF_SKY_FAR": _lsf(lsf, selected[far]),
            }
            meta = median_meta(
                header=header,
                path=path,
                input_index=input_index,
                expnum=expnum,
                positions=positions,
                separations=separations,
                near=near,
                far=far,
                counts=counts,
                sci_percentile=sci_percentile,
                sky_percentile=sky_percentile,
                gaia_sigma=gaia_sigma,
                gaia_ratio_threshold=gaia_ratio_threshold,
            )
            return {"slot": slot, "status": "OK", "arrays": arrays, "meta": meta, "expnum": expnum}
    except Exception as exc:  # noqa: BLE001 - record one bad exposure, keep the run alive
        return {
            "slot": slot,
            "status": "SKIP" if isinstance(exc, SkipExposure) else "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
            "expnum": expnum_from_path(path),
        }


def _combined_ivar(ivar: np.ndarray, sky_ivar: np.ndarray) -> np.ndarray:
    result = np.zeros_like(ivar, dtype=np.float32)
    valid = np.isfinite(ivar) & (ivar > 0) & np.isfinite(sky_ivar) & (sky_ivar > 0)
    result[valid] = 1.0 / (1.0 / ivar[valid] + 1.0 / sky_ivar[valid])
    return result


def _prepend(table: Table, name: str, values: Any, dtype: str) -> None:
    table.add_column(Column(values, name=name, dtype=dtype), index=0)


def _faint_worker(
    slot: int, input_index: int, path: Path, fibers_per_telescope: int
) -> dict[str, Any]:
    try:
        with fits.open(path, memmap=True) as hdul:
            header = hdul[0].header
            _validate_wave(np.asarray(hdul["WAVE"].data, dtype=np.float32))
            scale = _flux_scale(header, hdul["FLUX"])
            flux = np.asarray(hdul["FLUX"].data, dtype=np.float32) * scale
            sky = np.asarray(hdul["SKY"].data, dtype=np.float32) * scale
            ivar = np.asarray(hdul["IVAR"].data, dtype=np.float32) / scale**2
            sky_ivar = np.asarray(hdul["SKY_IVAR"].data, dtype=np.float32) / scale**2
            if len({array.shape for array in (flux, sky, ivar, sky_ivar)}) != 1:
                raise ValueError("FLUX/SKY/IVAR/SKY_IVAR shapes differ")
            slitmap = Table.read(hdul, hdu="SLITMAP")
            if len(slitmap) != flux.shape[0]:
                raise ValueError("SLITMAP row count differs from spectral arrays")
            combined = flux + sky
            telescopes = _strings(slitmap["telescope"])
            good = np.asarray(slitmap["fibstatus"]) == 0
            selected: list[np.ndarray] = []
            brightnesses: list[np.ndarray] = []
            for telescope in TELESCOPES:
                candidates = np.flatnonzero(good & (telescopes == telescope))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    brightness = np.nanmedian(combined[candidates], axis=1)
                finite = np.isfinite(brightness)
                candidates, brightness = candidates[finite], brightness[finite]
                if len(candidates) < fibers_per_telescope:
                    raise SkipExposure(
                        f"{telescope} has {len(candidates)} usable fibers; {fibers_per_telescope} requested"
                    )
                order = np.argsort(brightness, kind="stable")[:fibers_per_telescope]
                selected.append(candidates[order])
                brightnesses.append(brightness[order])
            indices = np.concatenate(selected)
            meta = slitmap[indices].copy()
            rows = len(indices)
            expnum = _expnum(path, header)
            _prepend(meta, "selection_median_flux", np.concatenate(brightnesses), "f4")
            _prepend(
                meta, "selection_rank", np.tile(np.arange(1, fibers_per_telescope + 1), 3), "i2"
            )
            _prepend(meta, "source_fiber_index", indices, "i4")
            _prepend(
                meta,
                "date_obs",
                [_text(header.get("DATE-OBS", header.get("OBSTIME")))] * rows,
                "U32",
            )
            _prepend(meta, "expnum", [expnum] * rows, "i8")
            _prepend(meta, "input_index", [input_index] * rows, "i8")
            _prepend(meta, "source_path", [str(path)] * rows, "U512")
            return {
                "slot": slot,
                "status": "OK",
                "flux": np.asarray(combined[indices], dtype=np.float32),
                "ivar": _combined_ivar(ivar[indices], sky_ivar[indices]),
                "meta": meta,
                "expnum": expnum,
            }
    except Exception as exc:  # noqa: BLE001 - record one bad exposure, keep the run alive
        return {
            "slot": slot,
            "status": "SKIP" if isinstance(exc, SkipExposure) else "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
            "expnum": expnum_from_path(path),
        }


def _bounded_results(
    worker: Callable[..., dict[str, Any]],
    tasks: list[tuple[Any, ...]],
    workers: int,
    wave: np.ndarray,
) -> Iterator[dict[str, Any]]:
    iterator = iter(tasks)
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(wave,)
    ) as pool:
        pending = set()

        def submit() -> bool:
            try:
                args = next(iterator)
            except StopIteration:
                return False
            pending.add(pool.submit(worker, *args))
            return True

        for _ in range(min(len(tasks), max(1, workers * 2))):
            submit()
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                submit()


def _memmaps(
    directory: Path, names: tuple[str, ...], shape: tuple[int, int]
) -> dict[str, np.memmap]:
    return {
        name: np.lib.format.open_memmap(
            directory / f"{name.lower()}.npy", mode="w+", dtype=np.float32, shape=shape
        )
        for name in names
    }


def _check_space(path: Path, required: int) -> None:
    free = shutil.disk_usage(path).free
    if free < int(required * 1.1):
        raise OSError(f"insufficient free space in {path}: need about {required / 2**30:.1f} GiB")


def _atomic_fits(hdus: fits.HDUList, output: Path, overwrite: bool) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {output}")
    with tempfile.NamedTemporaryFile(suffix=".fits", dir=output.parent, delete=False) as fh:
        temporary = Path(fh.name)
    try:
        hdus.writeto(temporary, overwrite=True, checksum=True)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _status_table(indexed: list[tuple[int, Path]], results: list[dict[str, Any] | None]) -> Table:
    rows = []
    output_index = 0
    for slot, (input_index, path) in enumerate(indexed):
        result = results[slot]
        status = result["status"] if result else "ERROR"
        rows.append(
            {
                "input_index": input_index,
                "path": str(path),
                "expnum": result.get("expnum", expnum_from_path(path))
                if result
                else expnum_from_path(path),
                "status": status,
                "error": result.get("error", "worker returned no result")
                if result
                else "worker returned no result",
                "output_index": output_index if status == "OK" else -1,
            }
        )
        output_index += status == "OK"
    return Table(rows=rows)


def build_stack(
    manifest: Path,
    output: Path,
    *,
    mode: str = "median",
    workers: int = 4,
    every_nth: int = 1,
    limit: int | None = None,
    sci_percentile: float = 70.0,
    sky_percentile: float = 70.0,
    gaia_dir: Path | None = None,
    gaia_sigma: int | None = None,
    gaia_ratio_threshold: float | None = None,
    fibers_per_telescope: int = 3,
    temp_dir: Path | None = None,
    keep_temp: bool = False,
    overwrite: bool = False,
    progress: Callable[[dict[str, int]], None] | None = None,
) -> dict[str, Any]:
    if workers < 1:
        raise ValueError("--workers must be positive")
    if mode not in {"median", "faint-fibers"}:
        raise ValueError(f"unknown build mode: {mode}")
    if gaia_sigma is not None and gaia_sigma not in {1, 3}:
        raise ValueError("--gaia-sigma must be 1 or 3")
    if gaia_sigma is not None and gaia_ratio_threshold is not None:
        raise ValueError("Gaia sigma and ratio selections are mutually exclusive")
    if mode != "median" and (gaia_sigma is not None or gaia_ratio_threshold is not None):
        raise ValueError("Gaia selection is only available in median mode")
    if not (0 < sci_percentile <= 100 and 0 < sky_percentile <= 100):
        raise ValueError("faint-fiber percentiles must be in (0, 100]")
    if gaia_ratio_threshold is not None and not 0 < gaia_ratio_threshold < 1:
        raise ValueError("--gaia-ratio-threshold must be greater than 0 and less than 1")
    if fibers_per_telescope < 1:
        raise ValueError("--fibers-per-telescope must be positive")
    if output.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {output}")

    indexed = read_manifest(manifest, every_nth, limit)
    paths = [path for _, path in indexed]
    wave, reference, drpver = first_wave(paths)
    output.parent.mkdir(parents=True, exist_ok=True)
    base = (temp_dir or output.parent).resolve()
    base.mkdir(parents=True, exist_ok=True)
    rows_per_input = 1 if mode == "median" else 3 * fibers_per_telescope
    array_count = len(OUTPUT_ARRAYS) if mode == "median" else 2
    bytes_per_copy = len(indexed) * rows_per_input * len(wave) * array_count * 4
    if base == output.parent.resolve():
        _check_space(base, bytes_per_copy * 2)
    else:
        _check_space(base, bytes_per_copy)
        _check_space(output.parent.resolve(), bytes_per_copy)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output.stem}-", dir=base))
    results: list[dict[str, Any] | None] = [None] * len(indexed)
    counts = {"completed": 0, "ok": 0, "skip": 0, "error": 0, "total": len(indexed)}

    try:
        if mode == "median":
            arrays = _memmaps(temporary, OUTPUT_ARRAYS, (len(indexed), len(wave)))
            tasks = [
                (
                    slot,
                    input_index,
                    path,
                    sci_percentile,
                    sky_percentile,
                    gaia_dir,
                    gaia_sigma,
                    gaia_ratio_threshold,
                )
                for slot, (input_index, path) in enumerate(indexed)
            ]
            iterator = _bounded_results(_median_worker, tasks, min(workers, len(tasks)), wave)
        else:
            arrays = _memmaps(
                temporary,
                ("FLUX", "IVAR"),
                (len(indexed) * rows_per_input, len(wave)),
            )
            tasks = [
                (slot, input_index, path, fibers_per_telescope)
                for slot, (input_index, path) in enumerate(indexed)
            ]
            iterator = _bounded_results(_faint_worker, tasks, min(workers, len(tasks)), wave)

        meta_chunks: list[Table | dict[str, Any] | None] = [None] * len(indexed)
        for result in iterator:
            slot = int(result["slot"])
            results[slot] = result
            status = result["status"].lower()
            counts[status] += 1
            counts["completed"] += 1
            if result["status"] != "OK":
                logging.getLogger("lvm_medians").error(
                    "input %s expnum=%s status=%s error=%s",
                    indexed[slot][1],
                    result.get("expnum"),
                    result["status"],
                    result.get("error", ""),
                )
            if result["status"] == "OK":
                if mode == "median":
                    for name in OUTPUT_ARRAYS:
                        arrays[name][slot] = result["arrays"][name]
                else:
                    start = slot * rows_per_input
                    stop = start + rows_per_input
                    arrays["FLUX"][start:stop] = result["flux"]
                    arrays["IVAR"][start:stop] = result["ivar"]
                meta_chunks[slot] = result["meta"]
            if progress:
                progress(counts.copy())

        good_slots = [i for i, result in enumerate(results) if result and result["status"] == "OK"]
        if not good_slots:
            first = next((result.get("error", "") for result in results if result), "")
            raise RuntimeError(f"no input produced output rows; first error: {first}")
        for array in arrays.values():
            array.flush()
        status_table = _status_table(indexed, results)
        header = fits.Header()
        header["CREATED"] = datetime.now(UTC).isoformat()
        header["DRPVER"] = drpver
        header["AGGMODE"] = mode.upper()
        header["NINPUT"] = len(indexed)
        header["NOUTPUT"] = len(good_slots)
        header["NSKIP"] = counts["skip"]
        header["NERROR"] = counts["error"]
        header["COMPLETE"] = counts["error"] == 0
        header["INLIST"] = str(manifest.resolve())
        header["REFWAVE"] = str(reference)
        hdus: list[Any] = [fits.PrimaryHDU(header=header)]
        wave_hdu = fits.ImageHDU(data=wave, name="WAVE")
        wave_hdu.header["BUNIT"] = "Angstrom"
        hdus.append(wave_hdu)
        if mode == "median":
            if gaia_sigma is not None:
                hdus[0].header["GAIAMODE"] = "SIGMA"
                hdus[0].header["GAIASIG"] = gaia_sigma
            elif gaia_ratio_threshold is not None:
                hdus[0].header["GAIAMODE"] = "RATIO"
                hdus[0].header["GAIATHR"] = gaia_ratio_threshold
            else:
                hdus[0].header["FAINTSCI"] = sci_percentile
                hdus[0].header["FAINTSKY"] = sky_percentile
            for destination, source in enumerate(good_slots):
                if destination != source:
                    for array in arrays.values():
                        array[destination] = array[source]
            for name in OUTPUT_ARRAYS:
                hdu = fits.ImageHDU(data=arrays[name][: len(good_slots)], name=name)
                if name.startswith("FLUX"):
                    hdu.header["BUNIT"] = FLUX_UNIT.to_string()
                hdus.append(hdu)
            meta = Table(rows=[meta_chunks[i] for i in good_slots])
        else:
            hdus[0].header["NFIBTEL"] = fibers_per_telescope
            for destination, source in enumerate(good_slots):
                src = slice(source * rows_per_input, (source + 1) * rows_per_input)
                dst = slice(destination * rows_per_input, (destination + 1) * rows_per_input)
                if destination != source:
                    for array in arrays.values():
                        array[dst] = array[src]
            n_rows = len(good_slots) * rows_per_input
            hdus[0].header["NROWS"] = n_rows
            for name in ("FLUX", "IVAR"):
                hdu = fits.ImageHDU(data=arrays[name][:n_rows], name=name)
                hdu.header["BUNIT"] = (
                    FLUX_UNIT.to_string() if name == "FLUX" else (FLUX_UNIT**-2).to_string()
                )
                hdus.append(hdu)
            meta = vstack([meta_chunks[i] for i in good_slots], metadata_conflicts="silent")
        hdus.append(fits.BinTableHDU(data=meta, name="META"))
        hdus.append(fits.BinTableHDU(data=status_table, name="INPUT_STATUS"))
        _atomic_fits(fits.HDUList(hdus), output, overwrite)
    finally:
        if not keep_temp:
            shutil.rmtree(temporary, ignore_errors=True)

    return {**counts, "output": str(output), "complete": counts["error"] == 0}
