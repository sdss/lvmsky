"""Stable META row schema consumed by downstream sky-model code."""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


def _float(header: fits.Header, *keys: str) -> float:
    for key in keys:
        try:
            value = float(header.get(key))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            return value
    return np.nan


def _int(header: fits.Header, key: str) -> int:
    try:
        return int(header.get(key))
    except (TypeError, ValueError, OverflowError):
        return -1


def _text(header: fits.Header, *keys: str) -> str:
    for key in keys:
        value: Any = header.get(key)
        if isinstance(value, bytes):
            value = value.decode(errors="replace")
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def median_meta(
    *,
    header: fits.Header,
    path: Path,
    input_index: int,
    expnum: int,
    positions: dict[str, tuple[float, float]],
    separations: dict[str, float],
    near: str,
    far: str,
    counts: dict[str, int],
    sci_percentile: float,
    sky_percentile: float,
    gaia_sigma: int | None,
    gaia_ratio_threshold: float | None,
) -> dict[str, Any]:
    fluxcal = _text(header, "FLUXCAL")
    pwv = _float(header, "PWV_MED")

    def fraction(telescope: str) -> float:
        total = counts[f"{telescope}_good"]
        return counts[f"{telescope}_used"] / total if total else np.nan

    return {
        "path": str(path),
        "input_index": input_index,
        "sci_faint_fiber_percentile": sci_percentile,
        "sky_faint_fiber_percentile": sky_percentile,
        "gaia_sigma": gaia_sigma if gaia_sigma is not None else -1,
        "gaia_ratio_threshold": (
            gaia_ratio_threshold if gaia_ratio_threshold is not None else np.nan
        ),
        "exposure": _int(header, "EXPOSURE"),
        "expnum": expnum,
        "mjd": _int(header, "MJD"),
        "tile_id": _int(header, "TILE_ID"),
        "obstime": _text(header, "OBSTIME", "DATE-OBS"),
        "date_obs": _text(header, "DATE-OBS", "OBSTIME"),
        "fluxcal": fluxcal,
        "pwv_med": pwv,
        "pwv_fallback": fluxcal == "MOD" and (not np.isfinite(pwv) or pwv <= 0),
        "sci_ra": positions["Sci"][0],
        "sci_dec": positions["Sci"][1],
        "skye_ra": positions["SkyE"][0],
        "skye_dec": positions["SkyE"][1],
        "skyw_ra": positions["SkyW"][0],
        "skyw_dec": positions["SkyW"][1],
        "sky_near_label": near,
        "sky_far_label": far,
        "sky_near_ra": positions[near][0],
        "sky_near_dec": positions[near][1],
        "sky_far_ra": positions[far][0],
        "sky_far_dec": positions[far][1],
        "skye_sep_deg": separations["SkyE"],
        "skyw_sep_deg": separations["SkyW"],
        "sky_near_sep_deg": separations[near],
        "sky_far_sep_deg": separations[far],
        "sci_alt": _float(header, "SKY SCI_ALT", "SCIALT", "ALT"),
        "skye_alt": _float(header, "SKY SKYE_ALT", "SKYEALT"),
        "skyw_alt": _float(header, "SKY SKYW_ALT", "SKYWALT"),
        "sci_airmass": _float(header, "SCIAM", "TESCIAM", "AIRMASS"),
        "skye_airmass": _float(header, "SKYEAM", "TESKYEAM"),
        "skyw_airmass": _float(header, "SKYWAM", "TESKYWAM"),
        "sci_moon_sep": _float(header, "SKY SCI_MOON_SEP"),
        "skye_moon_sep": _float(header, "SKY SKYE_MOON_SEP"),
        "skyw_moon_sep": _float(header, "SKY SKYW_MOON_SEP"),
        "moon_alt": _float(header, "SKY MOON_ALT", "MOONALT"),
        "sun_alt": _float(header, "SKY SUN_ALT", "SUNALT"),
        "moon_ra": _float(header, "SKY MOON_RA", "MOONRA"),
        "moon_dec": _float(header, "SKY MOON_DEC", "MOONDEC"),
        "moon_phase": _float(header, "SKY MOON_PHASE", "MOONPHA"),
        "moon_fli": _float(header, "SKY MOON_FLI"),
        "moon_illum": _float(header, "SKY MOON_FLI", "MOONILL"),
        "fibers_sci_good": counts["Sci_good"],
        "fibers_sci_used": counts["Sci_used"],
        "fibers_skye_good": counts["SkyE_good"],
        "fibers_skye_used": counts["SkyE_used"],
        "fibers_skyw_good": counts["SkyW_good"],
        "fibers_skyw_used": counts["SkyW_used"],
        "fibers_sky_near_good": counts[f"{near}_good"],
        "fibers_sky_near_used": counts[f"{near}_used"],
        "fibers_sky_far_good": counts[f"{far}_good"],
        "fibers_sky_far_used": counts[f"{far}_used"],
        "fibers_sci": counts["Sci_used"],
        "fiberfrac_sci": fraction("Sci"),
        "fibers_sky_near": counts[f"{near}_used"],
        "fiberfrac_sky_near": fraction(near),
        "fibers_sky_far": counts[f"{far}_used"],
        "fiberfrac_sky_far": fraction(far),
        "worker_pid": os.getpid(),
    }
