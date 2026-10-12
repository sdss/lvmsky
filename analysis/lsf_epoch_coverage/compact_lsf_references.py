"""Separate small usable LSF references from their full calculation archive."""
from pathlib import Path
import os

import numpy as np
from astropy.io import fits

CORE = ("PRIMARY", "REFERENCES", "LSF_COEF", "LSF_KNOTS", "LSF_CONFIG", "WAVE",
        "COEF_SCATTER", "COEF_CONDERR", "SOURCES", "SOURCE_CONFIG")


def write_reference_hdus(hdus, path, destination):
    """Publish the archive and a small standalone consumer file atomically."""
    path, destination = Path(path), Path(destination)
    archive = destination.parent / "lsf_reference_audit" / path.relative_to(destination)
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive_temp = archive.with_suffix(".tmp.fits")
    fits.HDUList(hdus).writeto(archive_temp, overwrite=True, checksum=True)
    compact = fits.HDUList([h.copy() for h in hdus if h.name in CORE])
    compact[0].header["PACKING"] = "compact"
    compact[0].header["AUDIT"] = os.path.relpath(archive, path.parent)
    temporary = path.with_suffix(".tmp.fits")
    compact.writeto(temporary, overwrite=True, checksum=True)
    os.replace(archive_temp, archive)
    os.replace(temporary, path)


def compact_reference(path, destination):
    """Preserve every coefficient exactly; retain full records in audit/."""
    path, destination = Path(path), Path(destination)
    archive = destination.parent / "lsf_reference_audit" / path.relative_to(destination)
    archive.parent.mkdir(parents=True, exist_ok=True)
    with fits.open(path, memmap=True) as source:
        if source[0].header.get("PACKING") == "compact":
            if not archive.exists():
                raise ValueError("Compact reference must retain its calculation archive")
            return archive
        if archive.exists():
            raise FileExistsError(f"Refusing to overwrite an existing archive: {archive}")
        hdus = fits.HDUList([h.copy() for h in source if h.name in CORE])
        hdus[0].header["PACKING"] = "compact"
        hdus[0].header["AUDIT"] = os.path.relpath(archive, path.parent)
        temporary = path.with_suffix(".compact.tmp.fits")
        hdus.writeto(temporary, overwrite=True, checksum=True)
        with fits.open(temporary) as saved:
            for name in ("LSF_COEF", "LSF_KNOTS", "WAVE", "COEF_SCATTER", "COEF_CONDERR"):
                np.testing.assert_array_equal(source[name].data, saved[name].data)
    os.replace(path, archive)
    try:
        os.replace(temporary, path)
    except BaseException:
        os.replace(archive, path)
        raise
    return archive
