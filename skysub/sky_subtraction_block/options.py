"""Shared defaults and English help; importing this module needs no fitting runtime."""
from __future__ import annotations

import math

DEFAULTS = dict(
    slit_count=4, offset_count=11, wave_count=4,
    offset_roughness=0.0, wave_roughness=0.0, slit_roughness=0.0,
    reference_prior=0.0, centered=False, keep_preparation=False,
)
HELP = {
    "ensemble": (
        "Trained ML model archive (.pt) used for ordinary sky inference, not training. "
        "Required unless reusing same-exposure saved state. Intermediate calculations stay "
        "in memory; no preliminary FITS is written by default. --ensemble is a compatible alias."
    ),
    "gaia_dir": (
        "Direct parent of lvmGAIA-fibers-00011017.fits (eight-digit exposure), "
        "e.g. /cache/gaia/fibers, not /cache/gaia. Requires PRIMARY.EXPNUM and a FIBERS "
        "table whose fiber IDs/coordinates match SLITMAP. Unset/missing/identity-mismatched "
        "cache queries Gaia DR3 at https://gaia.ari.uni-heidelberg.de/tap (3 attempts); "
        "source matches go to OUTPUT_PARENT/.skyml-gaia/. Malformed FITS/schema is an error."
    ),
    "seed": (
        "Advanced: reuse median/decomposition/ML/correction state from an existing compact "
        "SkyML FITS for this exact CFrame. Not a random seed. Input hash, fiber IDs and "
        "wavelengths must match. A snapshot can be made with --keep-preparation. "
        "Unset: compute in memory, or reuse an existing adjacent "
        "lvmSFrame-skyml-preparation-EXPNUM.fits for compatibility. Reused ML predictions "
        "are preserved even if a model file is supplied."
    ),
    "keep_preparation": (
        "Also save newly computed median/ML/correction state beside the final product as "
        "lvmSFrame-skyml-preparation-EXPNUM.fits for later --reuse-state reruns. "
        "Off by default: ordinary runs write only the final SFrame."
    ),
    "jobs": (
        "Worker processes per exposure (--jobs is a compatible alias) on parallel stages; 1 runs in the "
        "exposure process. Each process has one native numerical thread. Serial stages "
        "remain serial. Batch defaults to 1; single/legacy to 8."
    ),
    "slit_count": (
        "Cubic controls across physical fiber position within a spectrograph (integer >=4). "
        "More controls allow finer slit variation and require more Sky-fiber coverage."
    ),
    "offset_count": (
        "Cubic offset M-spline controls for R/Z profile shape (odd integer >=5). "
        "Support remains +/-2.5 Angstrom; this changes flexibility, not support width."
    ),
    "wave_count": (
        "Cubic controls for R/Z LSF variation with wavelength (integer >=4). "
        "More controls allow faster wavelength variation and add computational cost."
    ),
    "offset_roughness": (
        "Dimensionless penalty on squared second differences across offset coefficients. "
        "Larger values favor smoother profile coefficient patterns; 0 disables it. "
        "Finite >=0; scaled by mean native-IVAR data-Hessian diagonal. R/Z only."
    ),
    "wave_roughness": (
        "Dimensionless penalty on squared second differences along wavelength controls. "
        "Larger values favor smoother wavelength variation; 0 disables it. "
        "Finite >=0; scaled by mean data-Hessian diagonal. R/Z only."
    ),
    "slit_roughness": (
        "Dimensionless penalty on squared second differences along slit controls. "
        "Larger values favor smoother fiber-position variation; 0 disables it. "
        "Finite >=0; scaled by mean data-Hessian diagonal. R/Z only."
    ),
    "reference_prior": (
        "Dimensionless squared-distance penalty toward the Science-median LSF projected "
        "onto the requested basis and repeated across slit position. Larger values pull "
        "toward this reference; 0 disables this soft prior. Finite >=0; scaled by mean "
        "data-Hessian diagonal. Separate from automatic reference fallback. Data chi-square "
        "excludes all penalties."
    ),
    "centered": (
        "Constrain fitted R/Z profiles to zero first moment; default allows centroid variation. "
        "This is not a wavelength-shift fit. Unsupported blocks adopting the reference keep "
        "its centroid and record that choice in QA."
    ),
}


def validate_model(options):
    for name, minimum in (("slit_count", 4), ("wave_count", 4), ("offset_count", 5)):
        if options[name] < minimum:
            raise ValueError(f"--{name.replace('_', '-')} must be >= {minimum}")
    if options["offset_count"] % 2 != 1:
        raise ValueError("--offset-count must be odd")
    for name in ("offset_roughness", "wave_roughness", "slit_roughness", "reference_prior"):
        if not math.isfinite(options[name]) or options[name] < 0:
            raise ValueError(f"--{name.replace('_', '-')} must be finite and >= 0")
