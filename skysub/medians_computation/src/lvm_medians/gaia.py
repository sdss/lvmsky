"""Resumable Gaia TAP download and per-fiber photometry cache."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import random
import tempfile
import threading
import time
from collections.abc import Callable
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np
import requests
from astropy import units as u
from astropy.io import fits
from astropy.io.votable import parse
from astropy.table import Table, vstack

from .stack import _expnum, _flux_scale, expnum_from_path, read_manifest

TAP_SERVICES = {
    "ari": "https://gaia.ari.uni-heidelberg.de/tap",
    "esa": "https://gea.esac.esa.int/tap-server/tap",
    "aip": "https://gaia.aip.de/tap",
}
APERTURE_ARCSEC = 35.3 / 2
QUERY = f"""
SELECT f.expnum, f.fiberid, g.source_id, g.ra, g.dec,
       g.phot_g_mean_mag, g.phot_g_mean_flux_over_error
FROM TAP_UPLOAD.fibers AS f
JOIN gaiadr3.gaia_source AS g
  ON 1=CONTAINS(
      POINT('ICRS', g.ra, g.dec),
      CIRCLE('ICRS', f.ra, f.dec, {APERTURE_ARCSEC / 3600.0:.12f}))
""".strip()
_LOCAL = threading.local()


class TimeoutSession(requests.Session):
    def __init__(self, timeout: float, token: str | None = None):
        super().__init__()
        self.timeout = timeout
        if token:
            self.headers["Authorization"] = f"Bearer {token}"

    def request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        kwargs.setdefault("timeout", self.timeout)
        return super().request(method, url, **kwargs)


def resolve_service(value: str) -> str:
    return TAP_SERVICES.get(value.lower(), value.rstrip("/"))


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _atomic_table(table: Table, path: Path, header: fits.Header, name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".fits", dir=path.parent, delete=False) as fh:
        temporary = Path(fh.name)
    try:
        fits.HDUList(
            [fits.PrimaryHDU(header=header), fits.BinTableHDU(data=table, name=name)]
        ).writeto(temporary, overwrite=True, checksum=True)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _load_table(path: Path, name: str) -> Table:
    with fits.open(path, memmap=True) as hdul:
        return Table(hdul[name if name in hdul else 1].data)


def _fiber_cache_valid(path: Path, expnum: int, fiberids: np.ndarray) -> bool:
    if not path.is_file():
        return False
    required = {
        "expnum",
        "fiberid",
        "n_gaia_sources",
        "gaia_g_flux",
        "gaia_g_n_valid_sources",
        "lvm_flux_g_err",
        "lvm_flux_g_valid_fraction",
        "lvm_flux_plus_sky_g",
        "lvm_flux_plus_sky_g_valid_fraction",
        "ratio_gaia_to_lvm_flux_plus_sky_g",
    }
    try:
        table = _load_table(path, "FIBERS")
        cached_ids = np.asarray(table["fiberid"], dtype=np.int64)
        return (
            required.issubset(table.colnames)
            and len(table) == len(fiberids)
            and np.array_equal(np.sort(cached_ids), np.sort(fiberids.astype(np.int64)))
            and np.all(np.asarray(table["expnum"], dtype=np.int64) == expnum)
        )
    except Exception:  # noqa: BLE001 - corrupt/incompatible cache is a miss
        return False


def _passband(path: Path | None) -> tuple[np.ndarray, np.ndarray, float, float]:
    source = path or Path(str(files("lvm_medians").joinpath("data/GAIA3_G.vot")))
    resource = parse(source).get_first_table()
    params = {parameter.name: parameter.value for parameter in resource.params}
    table = resource.to_table()
    wavelength_column = table["Wavelength"]
    if getattr(wavelength_column, "unit", None) is not None:
        wavelength = wavelength_column.quantity.to_value(u.Angstrom)
    else:
        wavelength = np.asarray(wavelength_column, dtype=float)
    transmission = np.asarray(table["Transmission"], dtype=float)
    finite = np.isfinite(wavelength) & np.isfinite(transmission) & (transmission >= 0)
    wavelength = wavelength[finite]
    transmission = transmission[finite]
    order = np.argsort(wavelength)
    wavelength = wavelength[order]
    transmission = transmission[order]
    if wavelength.size < 2 or np.any(np.diff(wavelength) <= 0):
        raise ValueError("passband wavelengths must be finite, unique, and increasing")
    if str(params.get("ZeroPointUnit", "")).lower() != "jy":
        raise ValueError("passband ZeroPointUnit must be Jy")
    if int(float(params.get("DetectorType", -1))) != 1:
        raise ValueError("passband DetectorType must identify a photon counter (1)")
    return (
        wavelength,
        transmission,
        float(params["ZeroPoint"]),
        float(params.get("WavelengthPivot") or params["WavelengthRef"]),
    )


def _trapezoid_weights(x: np.ndarray) -> np.ndarray:
    weights = np.empty_like(x, dtype=float)
    weights[0] = (x[1] - x[0]) / 2
    weights[-1] = (x[-1] - x[-2]) / 2
    weights[1:-1] = (x[2:] - x[:-2]) / 2
    return weights


def _synthetic(
    spectra: np.ndarray,
    variance: np.ndarray,
    weights: np.ndarray,
    pixel_valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid = pixel_valid & np.isfinite(spectra) & np.isfinite(variance) & (variance >= 0)
    effective = valid * weights[None, :]
    normalization = effective.sum(axis=1)
    values = np.full(len(spectra), np.nan)
    errors = np.full(len(spectra), np.nan)
    coverage = normalization / weights.sum()
    good = normalization > 0
    coverage[~good] = np.nan
    values[good] = (
        np.sum(np.where(valid, spectra, 0) * effective, axis=1)[good] / normalization[good]
    )
    normalized = np.divide(
        effective,
        normalization[:, None],
        out=np.zeros_like(effective),
        where=normalization[:, None] > 0,
    )
    errors[good] = np.sqrt(np.sum(np.where(valid, variance, 0) * normalized**2, axis=1)[good])
    return values, errors, coverage


def derive_fiber_table(
    sframe: Path,
    sources: Table,
    passband: Path | None = None,
) -> tuple[Table, fits.Header]:
    with fits.open(sframe, memmap=True) as hdul:
        primary = hdul[0].header
        expnum = _expnum(sframe, primary)
        wave = np.asarray(hdul["WAVE"].data, dtype=float)
        scale = _flux_scale(primary, hdul["FLUX"])
        flux = np.asarray(hdul["FLUX"].data, dtype=float) * scale
        sky = np.asarray(hdul["SKY"].data, dtype=float) * scale
        ivar = np.asarray(hdul["IVAR"].data, dtype=float) / scale**2
        sky_ivar = np.asarray(hdul["SKY_IVAR"].data, dtype=float) / scale**2
        mask = np.asarray(hdul["MASK"].data) if "MASK" in hdul else np.zeros_like(flux)
        slitmap = Table.read(hdul, hdu="SLITMAP")
    if len({array.shape for array in (flux, sky, ivar, sky_ivar, mask)}) != 1:
        raise ValueError("FLUX/SKY/IVAR/SKY_IVAR/MASK shapes differ")
    if flux.ndim != 2 or flux.shape[1] != wave.size or len(slitmap) != len(flux):
        raise ValueError("SLITMAP row count differs from spectral arrays")

    filter_wave, transmission, zero_point_jy, pivot = _passband(passband)
    response = np.interp(wave, filter_wave, transmission, left=0.0, right=0.0)
    weights = response * wave * _trapezoid_weights(wave)
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Gaia G passband does not overlap the SFrame wavelength grid")
    filter_widths = _trapezoid_weights(filter_wave)
    in_lvm = (filter_wave >= np.nanmin(wave)) & (filter_wave <= np.nanmax(wave))
    photon_native = filter_widths * filter_wave * transmission
    energy_native = filter_widths * transmission
    photon_coverage = photon_native[in_lvm].sum() / photon_native.sum()
    energy_coverage = energy_native[in_lvm].sum() / energy_native.sum()
    variance = np.divide(1.0, ivar, out=np.full_like(ivar, np.nan), where=ivar > 0)
    sky_variance = np.divide(1.0, sky_ivar, out=np.full_like(sky_ivar, np.nan), where=sky_ivar > 0)
    lvm_flux, lvm_error, lvm_coverage = _synthetic(flux, variance, weights, mask == 0)
    total_flux, total_error, total_coverage = _synthetic(
        flux + sky, variance + sky_variance, weights, mask == 0
    )

    fiberids = np.asarray(slitmap["fiberid"], dtype=np.int64)
    position = {int(fiberid): index for index, fiberid in enumerate(fiberids)}
    n_sources = np.zeros(len(slitmap), dtype=np.int32)
    n_valid = np.zeros(len(slitmap), dtype=np.int32)
    gaia_flux = np.zeros(len(slitmap), dtype=float)
    if len(sources):
        source_fibers = np.asarray(sources["fiberid"], dtype=np.int64)
        magnitudes = np.asarray(sources["phot_g_mean_mag"], dtype=float)
        snr = np.asarray(sources["phot_g_mean_flux_over_error"], dtype=float)
        density_zero = zero_point_jy * 1e-23 * 2.99792458e18 / pivot**2
        source_flux = density_zero * 10 ** (-0.4 * magnitudes)
        for fiberid, magnitude, source_snr, value in zip(
            source_fibers, magnitudes, snr, source_flux
        ):
            index = position.get(int(fiberid))
            if index is None:
                continue
            n_sources[index] += 1
            if (
                np.isfinite(magnitude)
                and np.isfinite(source_snr)
                and source_snr > 0
                and np.isfinite(value)
            ):
                n_valid[index] += 1
                gaia_flux[index] += value
    gaia_mag = np.full(len(slitmap), np.nan)
    valid_gaia = gaia_flux > 0
    density_zero = zero_point_jy * 1e-23 * 2.99792458e18 / pivot**2
    gaia_mag[valid_gaia] = -2.5 * np.log10(gaia_flux[valid_gaia] / density_zero)

    def column(name: str, default: Any) -> np.ndarray:
        return (
            np.asarray(slitmap[name])
            if name in slitmap.colnames
            else np.full(len(slitmap), default)
        )

    output = Table(
        {
            "expnum": np.full(len(slitmap), expnum, dtype=np.int64),
            "fiberid": fiberids,
            "fiber_ra": column("ra", np.nan),
            "fiber_dec": column("dec", np.nan),
            "fibstatus": column("fibstatus", -1),
            "telescope": column("telescope", ""),
            "targettype": column("targettype", ""),
            "n_gaia_sources": n_sources,
            "gaia_g_flux": gaia_flux,
            "gaia_g_mag": gaia_mag,
            "gaia_g_n_valid_sources": n_valid,
            "svo_g_zero_point_jy": np.full(len(slitmap), zero_point_jy),
            "svo_g_pivot_angstrom": np.full(len(slitmap), pivot),
            "svo_g_photon_coverage": np.full(len(slitmap), photon_coverage),
            "svo_g_energy_coverage": np.full(len(slitmap), energy_coverage),
            "lvm_flux_g": lvm_flux,
            "lvm_flux_g_err": lvm_error,
            "lvm_flux_g_valid_fraction": lvm_coverage,
            "ratio_gaia_to_lvm_flux_g": np.divide(
                gaia_flux, lvm_flux, out=np.full(len(slitmap), np.nan), where=lvm_flux > 0
            ),
            "lvm_flux_plus_sky_g": total_flux,
            "lvm_flux_plus_sky_g_err": total_error,
            "lvm_flux_plus_sky_g_valid_fraction": total_coverage,
            "ratio_gaia_to_lvm_flux_plus_sky_g": np.divide(
                gaia_flux, total_flux, out=np.full(len(slitmap), np.nan), where=total_flux > 0
            ),
        }
    )
    header = fits.Header()
    header["CACHEVER"] = 1
    header["EXPNUM"] = expnum
    header["CREATED"] = _now()
    header["APRAD"] = APERTURE_ARCSEC
    header["PASSBAND"] = "GAIA/GAIA3.G"
    header["SFRAME"] = str(sframe)
    return output, header


def _service(url: str, timeout: float, token: str | None) -> Any:
    key = (url, timeout, token)
    if getattr(_LOCAL, "key", None) != key:
        try:
            import pyvo
        except ImportError as exc:
            raise RuntimeError(
                "pyvo is required for fetch-gaia; install the project dependencies"
            ) from exc
        _LOCAL.key = key
        _LOCAL.service = pyvo.dal.TAPService(url, session=TimeoutSession(timeout, token))
    return _LOCAL.service


def _download(
    sframe: Path,
    source_path: Path,
    service_url: str,
    timeout: float,
    token: str | None,
    maxrec: int,
) -> Table:
    with fits.open(sframe, memmap=True) as hdul:
        header = hdul[0].header
        expnum = _expnum(sframe, header)
        slitmap = Table.read(hdul, hdu="SLITMAP")
    upload = Table(
        {
            "expnum": np.full(len(slitmap), expnum, dtype=np.int64),
            "fiberid": np.asarray(slitmap["fiberid"], dtype=np.int64),
            "ra": np.asarray(slitmap["ra"], dtype=float),
            "dec": np.asarray(slitmap["dec"], dtype=float),
        }
    )
    result = (
        _service(service_url, timeout, token)
        .run_async(QUERY, uploads={"fibers": upload}, maxrec=maxrec)
        .to_table()
    )
    query_hash = hashlib.sha256(QUERY.encode()).hexdigest()[:16]
    cache_header = fits.Header()
    cache_header["CACHEVER"] = 1
    cache_header["EXPNUM"] = expnum
    cache_header["CREATED"] = _now()
    cache_header["TAPURL"] = service_url
    cache_header["QUERYID"] = query_hash
    cache_header["APRAD"] = APERTURE_ARCSEC
    _atomic_table(result, source_path, cache_header, "SOURCES")
    return result


def _load_failures(path: Path) -> dict[int, dict[str, Any]]:
    if not path.exists():
        return {}
    failures: dict[int, dict[str, Any]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            failures[int(row["expnum"])] = row
    return failures


def _write_failures(path: Path, failures: dict[int, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(json.dumps(failures[key], sort_keys=True) + "\n" for key in sorted(failures))
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as fh:
        fh.write(text)
        temporary = Path(fh.name)
    os.replace(temporary, path)


def failed_sframes(
    sframe_list: Path,
    cache_root: Path,
    every_nth: int = 1,
    limit: int | None = None,
) -> list[tuple[int, Path]]:
    """Return only SFrames recorded in the persistent failure ledger."""
    failures = _load_failures(cache_root / "gaia-failures.jsonl")
    return [
        item
        for item in read_manifest(sframe_list, every_nth, limit)
        if expnum_from_path(item[1]) in failures
    ]


def fetch_gaia(
    manifest: Path,
    cache_root: Path,
    *,
    service: str = "ari",
    workers: int = 5,
    retries: int = 3,
    timeout: float = 120,
    maxrec: int = 1_000_000,
    token: str | None = None,
    passband: Path | None = None,
    every_nth: int = 1,
    limit: int | None = None,
    retry_failed: bool = False,
    progress: Callable[[dict[str, int]], None] | None = None,
) -> dict[str, int]:
    if workers < 1 or retries < 0 or timeout <= 0 or maxrec < 1:
        raise ValueError("workers/timeout/maxrec must be positive and retries non-negative")
    indexed = (
        failed_sframes(manifest, cache_root, every_nth, limit)
        if retry_failed
        else read_manifest(manifest, every_nth, limit)
    )
    sources_dir = cache_root / "sources"
    fibers_dir = cache_root / "fibers"
    sources_dir.mkdir(parents=True, exist_ok=True)
    fibers_dir.mkdir(parents=True, exist_ok=True)
    failures_path = cache_root / "gaia-failures.jsonl"
    failures = _load_failures(failures_path)
    service_url = resolve_service(service)
    counts = {"total": len(indexed), "completed": 0, "cached": 0, "downloaded": 0, "failed": 0}
    if not indexed:
        return counts

    def one(item: tuple[int, Path]) -> tuple[int, str, str]:
        _, sframe = item
        with fits.open(sframe, memmap=True) as hdul:
            expnum = _expnum(sframe, hdul[0].header)
            fiberids = np.asarray(Table.read(hdul, hdu="SLITMAP")["fiberid"], dtype=np.int64)
        fiber_path = fibers_dir / f"lvmGAIA-fibers-{expnum:08d}.fits"
        source_path = sources_dir / f"lvmGAIA-sources-{expnum:08d}.fits"
        if _fiber_cache_valid(fiber_path, expnum, fiberids):
            return expnum, "cached", ""
        last_error = ""
        for attempt in range(retries + 1):
            try:
                if source_path.exists():
                    try:
                        source_table = _load_table(source_path, "SOURCES")
                    except Exception:  # noqa: BLE001 - replace a corrupt raw cache
                        source_table = _download(
                            sframe, source_path, service_url, timeout, token, maxrec
                        )
                else:
                    source_table = _download(
                        sframe, source_path, service_url, timeout, token, maxrec
                    )
                table, header = derive_fiber_table(sframe, source_table, passband)
                header["TAPURL"] = service_url
                _atomic_table(table, fiber_path, header, "FIBERS")
                return expnum, "downloaded", ""
            except Exception as exc:  # noqa: BLE001 - retry the complete TAP/cache operation
                last_error = f"{type(exc).__name__}: {exc}"
                logging.getLogger("lvm_medians").warning(
                    "Gaia expnum=%s attempt=%s/%s error=%s",
                    expnum,
                    attempt + 1,
                    retries + 1,
                    last_error,
                )
                if attempt < retries:
                    time.sleep(min(30.0, 2**attempt + random.random()))
        return expnum, "failed", last_error

    iterator = iter(indexed)
    with ThreadPoolExecutor(max_workers=min(workers, len(indexed))) as pool:
        pending = set()

        def submit() -> bool:
            try:
                item = next(iterator)
            except StopIteration:
                return False
            pending.add(pool.submit(one, item))
            return True

        for _ in range(min(len(indexed), workers * 2)):
            submit()
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                expnum, status, error = future.result()
                counts[status] += 1
                counts["completed"] += 1
                ledger_changed = False
                if status == "failed":
                    previous = failures.get(expnum, {})
                    failures[expnum] = {
                        "expnum": expnum,
                        "attempts": int(previous.get("attempts", 0)) + retries + 1,
                        "error": error,
                        "updated": _now(),
                    }
                    ledger_changed = True
                elif failures.pop(expnum, None) is not None:
                    ledger_changed = True
                if ledger_changed or not failures_path.exists():
                    _write_failures(failures_path, failures)
                if progress:
                    progress(counts.copy())
                submit()
    return counts


def combine_gaia_tables(
    manifest: Path,
    cache_root: Path,
    output: Path,
    *,
    table_kind: str = "fibers",
    every_nth: int = 1,
    limit: int | None = None,
    overwrite: bool = False,
    progress: Callable[[dict[str, int]], None] | None = None,
) -> dict[str, Any]:
    output = output.expanduser().resolve()
    if table_kind not in {"fibers", "sources"}:
        raise ValueError("--table must be fibers or sources")
    if output.suffix.lower() not in {".fit", ".fits", ".parquet"}:
        raise ValueError("--output must end in .fit, .fits, or .parquet")
    if output.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {output}")

    indexed = read_manifest(manifest, every_nth, limit)
    expnums = [expnum_from_path(path) for _, path in indexed]
    if any(expnum < 0 for expnum in expnums):
        raise ValueError("every SFrame-list path must contain an exposure number")

    extension = table_kind.upper()
    prefix = f"lvmGAIA-{table_kind}-"
    counts = {
        "total": len(expnums),
        "completed": 0,
        "combined": 0,
        "failed": 0,
        "rows": 0,
    }
    tables: list[Table] = []
    logger = logging.getLogger("lvm_medians")
    for expnum in expnums:
        path = cache_root / table_kind / f"{prefix}{expnum:08d}.fits"
        try:
            table = _load_table(path, extension)
            if len(table) and (
                "expnum" not in table.colnames
                or not np.all(np.asarray(table["expnum"], dtype=np.int64) == expnum)
            ):
                raise ValueError(f"expnum column does not match {expnum}")
        except Exception as exc:  # noqa: BLE001 - report corrupt/missing cache and continue
            counts["failed"] += 1
            logger.error("cannot combine Gaia cache %s: %s: %s", path, type(exc).__name__, exc)
        else:
            tables.append(table)
            counts["combined"] += 1
            counts["rows"] += len(table)
        counts["completed"] += 1
        if progress:
            progress(counts.copy())

    if not tables:
        raise RuntimeError(f"no readable Gaia {table_kind} cache files")
    combined = vstack(tables, metadata_conflicts="silent")
    complete = counts["failed"] == 0
    created = _now()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() in {".fit", ".fits"}:
        header = fits.Header()
        header["CREATED"] = created
        header["TABKIND"] = table_kind
        header["NEXP"] = counts["combined"]
        header["NMISSING"] = counts["failed"]
        header["NROWS"] = len(combined)
        header["COMPLETE"] = complete
        header["INLIST"] = str(manifest.resolve())
        _atomic_table(combined, output, header, extension)
    else:
        try:
            import pyarrow  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "Parquet output requires: python -m pip install '.[parquet]'"
            ) from exc
        combined.meta.update(
            {
                "created": created,
                "table_kind": table_kind,
                "sframe_list": str(manifest.resolve()),
                "exposures": counts["combined"],
                "missing": counts["failed"],
                "complete": complete,
            }
        )
        with tempfile.NamedTemporaryFile(
            suffix=".parquet", dir=output.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
        try:
            combined.write(temporary, format="parquet", overwrite=True)
            os.replace(temporary, output)
        finally:
            temporary.unlink(missing_ok=True)

    return {
        **counts,
        "complete": complete,
        "table": table_kind,
        "output": str(output),
    }


def cache_status(manifest: Path, cache_root: Path) -> dict[str, int]:
    indexed = read_manifest(manifest)
    expected = {expnum_from_path(path) for _, path in indexed}
    cached = {
        expnum_from_name(path) for path in (cache_root / "fibers").glob("lvmGAIA-fibers-*.fits")
    }
    failures = _load_failures(cache_root / "gaia-failures.jsonl")
    return {
        "sframes": len(expected),
        "ready": len(expected & cached),
        "awaiting_download": len(expected - cached),
        "failures": len(expected & set(failures)),
    }


def expnum_from_name(path: Path) -> int:
    digits = path.stem.rsplit("-", 1)[-1]
    return int(digits) if digits.isdigit() else -1
