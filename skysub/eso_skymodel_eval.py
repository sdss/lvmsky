#!/usr/bin/env python3
"""
Short, importable evaluator for the ESO Sky Model (local calcskymodel backend).

- Importable API: evaluate_single(), evaluate_batch()
- Click CLI: single and batch modes
- Uses ESO_SKY_MODEL env var or --eso-sky-model path
- Runs in an isolated temporary working directory; no project tree mutations
- Avoids writing intermediates; optionally writes final FITS only

Docstrings are concise per user rule.
"""

import os
import sys
import tempfile
import subprocess
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple, Union

import numpy as np
import click
from astropy.io import fits
from astropy.table import Table, join
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import (
    get_body,
    solar_system_ephemeris,
    AltAz,
    EarthLocation,
    SkyCoord,
    GeocentricTrueEcliptic,
)
from tqdm.auto import tqdm

# Columns available from the generated table for per-array FITS output
_AVAILABLE_COLS = "FLUX,DFLUX1,DFLUX2,MOON,FLUX_SSL,ZODI,FLUX_TIE,FLUX_TME,LINES,DIFFUSE,TRANS,DTRANS1,DTRANS2,TRANS_MA,TRANS_O3,TRANS_RS,TRANS_MS,CONT"

# ---------- Time utilities ----------

def convert_time(time_input: Union[str, float, int, "datetime.datetime"], output_format: str = "iso_ms"):
    """Convert iso/mjd/jd/datetime to requested format (iso, iso_ms, mjd, jd, datetime)."""
    import datetime as _dt
    import re as _re

    if isinstance(time_input, _dt.datetime):
        t = Time(time_input, scale="utc")
    else:
        s = str(time_input)
        try:
            val = float(s)
            if val > 2400000:
                t = Time(val, format="jd", scale="utc")
            else:
                t = Time(val, format="mjd", scale="utc")
        except ValueError:
            # accept YYYY-MM-DDTHH:MM:SS[.sss]
            if _re.search(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?", s):
                t = Time(s, format="isot", scale="utc")
            else:
                t = Time(s, scale="utc")
    if output_format == "datetime":
        return t.to_datetime()
    if output_format == "iso":
        return t.iso
    if output_format == "iso_ms":
        return t.to_datetime().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
    if output_format == "mjd":
        return float(t.mjd)
    if output_format == "jd":
        return float(t.jd)
    raise ValueError(f"Unsupported output_format: {output_format}")

# ---------- Sky geometry at La Silla ----------

_LA_SILLA = EarthLocation(lat=-29.2583 * u.deg, lon=-70.7346 * u.deg, height=2400 * u.m)


def _info_la_silla(datetime_utc: str, ra: float, dec: float) -> dict:
    """Compute solar/lunar geometry and ecliptic coords for ESO inputs."""
    obs_time = Time(datetime_utc)
    # Use ICRS for the input source; we'll compute separations in AltAz to avoid frame-mixing warnings.
    src = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    with solar_system_ephemeris.set("builtin"):
        moon = get_body("moon", obs_time, location=_LA_SILLA)
        sun = get_body("sun", obs_time, location=_LA_SILLA)
    altaz = AltAz(obstime=obs_time, location=_LA_SILLA)
    moon_altaz = moon.transform_to(altaz)
    sun_altaz = sun.transform_to(altaz)
    src_altaz = src.transform_to(altaz)
    moon_ecl = moon.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    sun_ecl = sun.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    src_ecl = src.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    # illuminated fraction proxy (same frame: GCRS bodies are comparable)
    phase_angle = moon.separation(sun).radian
    illum = (1 - np.cos(phase_angle)) / 2.0
    mdiff = (moon.ra - sun.ra).wrap_at(360 * u.deg).value
    mphase = illum / 2.0 if mdiff > 0 else 1 - illum / 2.0
    mean_moon_dist = 384400 * u.km
    md = moon.distance.to(u.km)
    return {
        "SunEclipLon": float(sun_ecl.lon.deg),
        "SunEclipLat": float(sun_ecl.lat.deg),
        "MoonAlt": float(moon_altaz.alt.deg),
        "MoonAz": float(moon_altaz.az.deg),
        "MoonEclipLon": float(moon_ecl.lon.deg - (360 if moon_ecl.lon.deg > 180 else 0)),
        "MoonEclipLat": float(moon_ecl.lat.deg),
        "MoonIll": float(illum * 100.0),
        "MoonPhas": float(mphase),
        "MoonDistanceInMeanUnits": float((md / mean_moon_dist).value),
        "SourceAlt": float(src_altaz.alt.deg),
        "SourceAz": float(src_altaz.az.deg),
        "SourceEclipLon": float(src_ecl.lon.deg),
        "SourceEclipLat": float(src_ecl.lat.deg),
        # separations computed in a common frame to avoid cross-frame warnings
        "Moon-Sun_Separation": float(moon_altaz.separation(sun_altaz).deg),
        "Moon-Source_Separation": float(moon_altaz.separation(src_altaz).deg),
    }

# ---------- Config writers (to temp dir) ----------

_INST_BASE = """
# Wavelength grid:
limlam     = 0.36 0.98
# step size [mum]
dlam       = 0.00005
# Line-spread function:
kernrad    = 3
wbox       = 0.8
wgauss     = 0.8
wlorentz   = 0.8
varkern    = 1
kernelfile = output/kernel.dat
""".strip()

_OBS_TMPL = """
sm_h = 2.5
sm_hmin = 2.0
alt      = {alt:.1f}
alpha    = {alpha:.1f}
rho      = {rho:.1f}
altmoon  = {altmoon:.1f}
moondist = {moondist:.2f}
pres     = 744.
ssa      = 0.97
calcds   = N
o3column = 1.
moonscal = 1.0
lon_ecl  = {lon_ecl:.1f}
lat_ecl  = {lat_ecl:.1f}
emis_str = 0.2
temp_str = 290.
msolflux = 101.
season   = 0
time     = 0 
vac_air  = air
pwv      = 3.5
rtcode   = L
resol    = 6e4
filepath = data
incl     = YYYYYYY
""".strip()


@dataclass
class EvalResult:
    """Container for one evaluation result."""
    table: Table
    meta: dict


def _ensure_calcskymodel(eso_sky_model: Optional[str]) -> str:
    """Return path to calcskymodel or raise with clear message."""
    base = eso_sky_model or os.getenv("ESO_SKY_MODEL")
    if not base:
        raise RuntimeError("ESO_SKY_MODEL not set and --eso-sky-model not provided.")
    exe = os.path.join(base, "sm-01_mod2", "bin", "calcskymodel")
    if not os.path.isfile(exe):
        raise RuntimeError(f"calcskymodel not found at {exe}. Ensure ESO_SKY_MODEL is correct.")
    return exe


def _run_one_local(ra: float, dec: float, obstime_iso_ms: str, eso_sky_model: Optional[str]) -> EvalResult:
    """Run calcskymodel in a temp dir, load outputs, build LVM-like table in memory."""
    exe = _ensure_calcskymodel(eso_sky_model)
    base = eso_sky_model or os.getenv("ESO_SKY_MODEL")
    data_dir = os.path.join(base, "sm-01_mod2", "data")
    if not os.path.isdir(data_dir):
        raise RuntimeError(f"Data dir not found: {data_dir}")

    info = _info_la_silla(obstime_iso_ms, ra, dec)
    # sanity
    if info["SourceAlt"] < 0:
        raise RuntimeError(f"Source altitude < 0 deg for RA={ra}, Dec={dec} at {obstime_iso_ms}")

    with tempfile.TemporaryDirectory(prefix="eso_sky_") as td:
        # local working tree
        cfg = os.path.join(td, "config")
        outd = os.path.join(td, "output")
        os.makedirs(cfg, exist_ok=True)
        os.makedirs(outd, exist_ok=True)
        # symlink data
        os.symlink(data_dir, os.path.join(td, "data"))
        # write config files
        with open(os.path.join(cfg, "instrument_etc.par"), "w") as f:
            f.write(_INST_BASE + "\n")
        # longitude difference in [-180, 180]
        lon = info["SourceEclipLon"] - info["SunEclipLon"]
        lon = lon - 360 if lon > 180 else lon
        lon = lon + 360 if lon < -180 else lon
        obs_txt = _OBS_TMPL.format(
            alt=info["SourceAlt"],
            alpha=info["Moon-Sun_Separation"],
            rho=info["Moon-Source_Separation"],
            altmoon=info["MoonAlt"],
            moondist=info["MoonDistanceInMeanUnits"],
            lon_ecl=lon,
            lat_ecl=info["SourceEclipLat"],
        )
        with open(os.path.join(cfg, "skymodel_etc.par"), "w") as f:
            f.write(obs_txt + "\n")

        # run calcskymodel
        env = os.environ.copy()
        cwd = td
        proc = subprocess.run([exe], cwd=cwd, capture_output=True, text=True)
        if proc.returncode != 0 or not os.path.isfile(os.path.join(outd, "radspec.fits")):
            raise RuntimeError(
                "calcskymodel failed. stderr=" + proc.stderr.strip()
            )

        # load outputs (no permanent files)
        rfile = os.path.join(outd, "radspec.fits")
        tfile = os.path.join(outd, "transspec.fits")
        with fits.open(rfile) as _rhdul:
            rtab = Table(_rhdul[1].data)
        with fits.open(tfile) as _thdul:
            ttab = Table(_thdul[1].data)
        ztab = join(rtab, ttab, join_type="left")
        ztab["lam"] *= 1000.0
        # rename + compute LVM-like columns
        ztab.rename_column("lam", "WAVE")
        ztab.rename_column("flux", "FLUX")
        ztab.rename_column("flux_sml", "MOON")
        ztab.rename_column("flux_zl", "ZODI")
        ztab.rename_column("flux_ael", "LINES")
        ztab.rename_column("flux_arc", "DIFFUSE")
        area = np.pi * (37 / 2) ** 2
        q = 1.98644586e-17 * area / ztab["WAVE"]
        for c in ("FLUX", "MOON", "LINES", "ZODI", "DIFFUSE"):
            ztab[c] *= q
        ztab["WAVE"] *= 10.0
        ztab["CONT"] = ztab["MOON"] + ztab["ZODI"] + ztab["DIFFUSE"]
        for c in ("CONT", "MOON", "ZODI", "DIFFUSE"):
            ztab[c] /= ztab["trans"]

    meta = {
        "RA": float(ra),
        "Dec": float(dec),
        "ObsTime": obstime_iso_ms,
        "Observatory": "lasilla",
    }
    return EvalResult(table=ztab, meta=meta)


def evaluate_single(
    ra: float,
    dec: float,
    obstime: Union[str, float, int, "datetime.datetime"],
    *,
    eso_sky_model: Optional[str] = None,
    save_fits: Optional[str] = None,
    return_format: str = "table",
    hdu_extname: Optional[str] = None,
    dtype = np.float32,
) -> Union[Table, Tuple[np.ndarray, dict]]:
    """Evaluate for one (ra, dec, obstime). Optionally save FITS or return arrays.

    return_format: 'table' (Astropy Table) or 'arrays' -> (dict_of_arrays, meta).
    """
    # dtype: accept only np.float32 or np.float64
    if dtype not in (np.float32, np.float64, np.dtype("float32"), np.dtype("float64")):
        raise ValueError("dtype must be np.float32 or np.float64")
    np_dtype = np.float32 if str(np.dtype(dtype)) == "float32" else np.float64
    iso = convert_time(obstime, "iso_ms")
    res = _run_one_local(ra, dec, iso, eso_sky_model)

    if save_fits:
        hdul = fits.HDUList([fits.PrimaryHDU()])
        hdu = fits.BinTableHDU(data=res.table)
        if hdu_extname:
            hdu.name = hdu_extname
        hdul.append(hdu)
        for k, v in res.meta.items():
            hdul[0].header[k[:8]] = v
        hdul.writeto(save_fits, overwrite=True)

    if return_format == "arrays":
        cols = {c: np.asarray(res.table[c], dtype=np_dtype) for c in res.table.colnames}
        return cols, res.meta
    return res.table


def evaluate_batch(
    records: Iterable[Tuple[float, float, Union[str, float, int]]],
    *,
    eso_sky_model: Optional[str] = None,
    ncpu: int = 1,
    save_fits: Optional[str] = None,
    ignore_errors: bool = False,
    show_progress: bool = True,
    columns: Optional[List[str]] = None,
    dtype = np.float32,
) -> List[EvalResult]:
    """Evaluate multiple records; optionally write FITS in per-array layout only."""
    items = list(records)
    if ncpu is None or ncpu < 1:
        ncpu = 1

    # choose columns to write (default minimal set)
    default_cols = ["FLUX", "MOON", "ZODI", "DIFFUSE", "LINES", "CONT"]
    sel_cols = [c.strip().upper() for c in columns] if columns else default_cols
    # dtype handling: accept only np.float32 or np.float64
    if dtype not in (np.float32, np.float64, np.dtype("float32"), np.dtype("float64")):
        raise ValueError("dtype must be np.float32 or np.float64")
    np_dtype = np.float32 if str(np.dtype(dtype)) == "float32" else np.float64

    raw: List[Union[EvalResult, Tuple]]
    if ncpu == 1:
        raw = [None] * len(items)  # type: ignore[assignment]
        bar = tqdm(total=len(items), disable=not show_progress)
        try:
            for i, x in enumerate(items):
                raw[i] = _job_mp(x, eso_sky_model)
                bar.update(1)
        finally:
            bar.close()
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        raw = [None] * len(items)  # type: ignore[assignment]
        with ProcessPoolExecutor(max_workers=ncpu) as ex:
            futures = {ex.submit(_job_mp, items[i], eso_sky_model): i for i in range(len(items))}
            bar = tqdm(total=len(items), disable=not show_progress)
            try:
                for fut in as_completed(futures):
                    i = futures[fut]
                    raw[i] = fut.result()
                    bar.update(1)
            finally:
                bar.close()

    # Collate results and handle errors per-record
    results: List[EvalResult] = []
    errors = 0
    for i, r in enumerate(raw):
        if isinstance(r, tuple) and len(r) >= 2 and r[0] == "error":
            errors += 1
            if not ignore_errors:
                raise RuntimeError(f"Record {i} failed: {r[1]}")
            else:
                # fabricate a NaN table to keep output aligned
                try:
                    # r may be ("error", msg, ra, dec, iso) or just ("error", msg)
                    if len(r) >= 5:
                        _, msg, ra, dec, iso = r[:5]
                    else:
                        # fallback if worker didn't return coords
                        ra = dec = np.nan
                        iso = ""
                        msg = r[1]
                    nan_res = _make_nan_result(ra, dec, iso, error_msg=msg)
                    results.append(nan_res)
                except Exception:
                    pass
                continue
        results.append(r)  # type: ignore[arg-type]

    if save_fits:
        if len(results) == 0:
            raise RuntimeError("No successful records; nothing to write.")
        # Aggregate per-array 2D images
        colnames = results[0].table.colnames
        wave = np.asarray(results[0].table["WAVE"], dtype=float)
        nlam = wave.size
        nrec = len(results)
        # restrict to requested existing numeric columns (exclude WAVE)
        existing_upper = {c.upper(): c for c in colnames if c != "WAVE"}
        numeric_cols = [existing_upper[c] for c in sel_cols if c in existing_upper]
        # Save with numpy shape (nrec, nlam) so fits.info shows Dimensions (nlam, nrec)
        # i.e., (12401, 7947) as requested
        stacks = {c: np.full((nrec, nlam), np.nan, dtype=np_dtype) for c in numeric_cols}
        for i, r in enumerate(results):
            for c in numeric_cols:
                try:
                    vec = np.asarray(r.table[c], dtype=np_dtype)
                    # place vec along axis-1 for record i
                    stacks[c][i, :] = vec
                except Exception:
                    pass
        ra_arr = np.array([r.meta.get("RA", np.nan) for r in results], dtype=float)
        dec_arr = np.array([r.meta.get("Dec", np.nan) for r in results], dtype=float)
        obstime_arr = np.array([str(r.meta.get("ObsTime", "")) for r in results])
        err_arr = np.array([str(r.meta.get("Error", "")) for r in results])

        hdul = fits.HDUList([fits.PrimaryHDU()])
        hdul.append(fits.ImageHDU(data=wave.astype(np_dtype, copy=False), name="WAVE"))
        for c in numeric_cols:
            hdul.append(fits.ImageHDU(data=stacks[c], name=c.upper()))
        meta_tab = Table({"RA": ra_arr, "Dec": dec_arr, "ObsTime": obstime_arr, "Error": err_arr})
        hdul.append(fits.BinTableHDU(data=meta_tab, name="META"))
        hdul[0].header["LAYOUT"] = "per_array"
        hdul.writeto(save_fits, overwrite=True)

    return results


def _job_mp(x: Tuple[float, float, Union[str, float, int]], eso_sky_model: Optional[str]):
    """Top-level worker for multiprocessing (picklable).

    Returns either EvalResult or ("error", message) tuple.
    """
    try:
        ra, dec, t = x
        iso = convert_time(t, "iso_ms")
        return _run_one_local(ra, dec, iso, eso_sky_model)
    except Exception as e:
        # include coords/time for downstream NaN fabrication
        try:
            ra, dec, t = x
            iso = convert_time(t, "iso_ms")
        except Exception:
            ra = dec = np.nan
            iso = ""
        return ("error", str(e), ra, dec, iso)


def _make_nan_result(ra: float, dec: float, obstime_iso_ms: str, error_msg: Optional[str] = None) -> EvalResult:
    """Create a placeholder EvalResult with NaN columns and standard WAVE grid."""
    n = 12401
    wave = np.linspace(3600.0, 9800.0, n)
    cols = {
        "WAVE": wave,
        "FLUX": np.full(n, np.nan),
        "dflux1": np.full(n, np.nan),
        "dflux2": np.full(n, np.nan),
        "MOON": np.full(n, np.nan),
        "flux_ssl": np.full(n, np.nan),
        "ZODI": np.full(n, np.nan),
        "flux_tie": np.full(n, np.nan),
        "flux_tme": np.full(n, np.nan),
        "LINES": np.full(n, np.nan),
        "DIFFUSE": np.full(n, np.nan),
        "trans": np.full(n, np.nan),
        "dtrans1": np.full(n, np.nan),
        "dtrans2": np.full(n, np.nan),
        "trans_ma": np.full(n, np.nan),
        "trans_o3": np.full(n, np.nan),
        "trans_rs": np.full(n, np.nan),
        "trans_ms": np.full(n, np.nan),
        "CONT": np.full(n, np.nan),
    }
    tab = Table(cols)
    meta = {"RA": float(ra), "Dec": float(dec), "ObsTime": obstime_iso_ms, "Observatory": "lasilla"}
    if error_msg:
        meta["Error"] = error_msg
    return EvalResult(table=tab, meta=meta)


# ---------- CLI ----------

@click.group(context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120})
def cli():
    """ESO Sky Model evaluator (local calcskymodel)."""


@cli.command("single", context_settings={"allow_interspersed_args": False})
@click.argument("ra", type=float)
@click.argument("dec", type=float)
@click.argument("obstime", type=str)
@click.option("--eso-sky-model", type=click.Path(exists=True, file_okay=False), default=None, help="Override ESO_SKY_MODEL base path.")
@click.option("--save-fits", type=click.Path(dir_okay=False), default=None, help="Write final FITS (one table).")
@click.option("--extname", type=str, default=None, help="Extension name for saved HDU.")
@click.option("--return-format", type=click.Choice(["table", "arrays"]), default="table")
def cli_single(ra, dec, obstime, eso_sky_model, save_fits, extname, return_format):
    """Run one evaluation.

    Examples:

    1. Return table only

        $ python eso_skymodel_eval.py single 121.75 -29.7 2012-07-17T21:12:14

    2. Save FITS with custom extname

        $ python eso_skymodel_eval.py single --save-fits sky_single.fits --extname SKY 121.75 -29.7 2012-07-17T21:12:14

    3. Override ESO_SKY_MODEL location

        $ python eso_skymodel_eval.py single --eso-sky-model /path/to/ESO_SKY_MODEL 121.75 -29.7 2012-07-17T21:12:14

    4. Return numpy arrays

        $ python eso_skymodel_eval.py single --return-format arrays 121.75 -29.7 2012-07-17T21:12:14
    """
    try:
        tab_or_arrays = evaluate_single(
            ra, dec, obstime, eso_sky_model=eso_sky_model, save_fits=save_fits, return_format=return_format, hdu_extname=extname
        )
        if return_format == "table":
            click.echo(f"Rows: {len(tab_or_arrays)}; Cols: {len(tab_or_arrays.colnames)}")
        else:
            cols, meta = tab_or_arrays
            click.echo(f"Arrays: {list(cols.keys())}; N={len(next(iter(cols.values())))}")
    except Exception as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)


@cli.command("batch", context_settings={"allow_interspersed_args": False})
@click.option("--csv", "csv_path", type=click.Path(exists=True, dir_okay=False), default=None, help="CSV with columns ra,dec,obstime")
@click.option("--ncpu", type=int, default=1, show_default=True)
@click.option("--eso-sky-model", type=click.Path(exists=True, file_okay=False), default=None)
@click.option("--save-fits", type=click.Path(dir_okay=False), default=None, help="Write FITS output (per-array layout only). Columns saved are float32 by default.")
@click.option("--ignore-errors", is_flag=True, default=False, help="Skip records that fail (e.g., below horizon).")
@click.option("--no-progress", is_flag=True, default=False, help="Disable progress bar.")
@click.option("--columns", type=str, default=None, help=f"Comma-separated list of columns to save. Available: {_AVAILABLE_COLS}. Default: FLUX,MOON,ZODI,DIFFUSE,LINES,CONT")
@click.option("--dtype", type=click.Choice(["float32", "float64"]), default="float32", show_default=True, help="Data type for saved arrays.")
def cli_batch(csv_path, ncpu, eso_sky_model, save_fits, ignore_errors, no_progress, columns, dtype):
    """Run batch evaluation from CSV (ra,dec,obstime).

    Examples:

    1. Run on a CSV with 4 processes and save FITS

        $ python eso_skymodel_eval.py batch --csv inputs.csv --ncpu 4 --save-fits sky_batch.fits

    2. Override ESO_SKY_MODEL location

        $ python eso_skymodel_eval.py batch --csv inputs.csv --eso-sky-model /path/to/ESO_SKY_MODEL

    CSV format (header required):
        ra,dec,obstime
        121.75,-29.7,2012-07-17T21:12:14
        200.1,10.5,2019-05-20T08:01:02.123
    """
    try:
        if csv_path is None:
            click.echo("ERROR: --csv is required.", err=True)
            sys.exit(2)
        import csv as _csv

        items: List[Tuple[float, float, Union[str, float, int]]] = []
        with open(csv_path, newline="") as f:
            rdr = _csv.DictReader(f)
            for row in rdr:
                items.append((float(row["ra"]), float(row["dec"]), row["obstime"]))
        # parse columns
        cols = [s.strip() for s in columns.split(",")] if columns else None
        # map CLI dtype string to numpy type
        np_cli_dtype = np.float32 if dtype == "float32" else np.float64
        res = evaluate_batch(
            items,
            eso_sky_model=eso_sky_model,
            ncpu=ncpu,
            save_fits=save_fits,
            ignore_errors=ignore_errors,
            show_progress=not no_progress,
            columns=cols,
            dtype=np_cli_dtype,
        )
        click.echo(f"Completed {len(res)} records.")
        if ignore_errors:
            err_indices = []
            for i, r in enumerate(res):
                err = r.meta.get("Error")
                if err:
                    err_indices.append(i)
            if err_indices:
                click.echo(f"Skipped {len(err_indices)} failed record(s) due to errors:", err=True)
                for i in err_indices:
                    m = res[i].meta
                    click.echo(
                        f"  - idx={i} RA={m.get('RA'):.3f} Dec={m.get('Dec'):.3f} ObsTime={m.get('ObsTime')} | {m.get('Error')}",
                        err=True,
                    )
    except Exception as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
