from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from typer.main import get_command

from lvm_medians.cli import app
from lvm_medians.gaia import derive_fiber_table
from lvm_medians.stack import (
    build_stack,
    gaia_clean_mask,
    gaia_ratio_clean_mask,
    read_manifest,
)


def make_sframe(path: Path) -> None:
    wave = np.linspace(4000, 8000, 8, dtype=np.float32)
    flux = np.arange(1, 7, dtype=np.float32)[:, None] * np.ones((6, 8), dtype=np.float32)
    mask = np.zeros_like(flux, dtype=np.int32)
    mask[0, 0] = 1
    sky = np.ones_like(flux)
    slitmap = Table(
        {
            "fiberid": np.arange(1, 7),
            "ra": np.linspace(10, 11, 6),
            "dec": np.linspace(-2, -1, 6),
            "fibstatus": np.zeros(6, dtype=np.int16),
            "telescope": ["Sci", "Sci", "SkyE", "SkyE", "SkyW", "SkyW"],
            "targettype": ["science"] * 6,
        }
    )
    header = fits.Header({"EXPNUM": 1, "DRPVER": "test", "BUNIT": "erg / (Angstrom s cm2)"})
    fits.HDUList(
        [
            fits.PrimaryHDU(header=header),
            fits.ImageHDU(wave, name="WAVE"),
            fits.ImageHDU(flux, name="FLUX"),
            fits.ImageHDU(sky, name="SKY"),
            fits.ImageHDU(np.full_like(flux, 4), name="IVAR"),
            fits.ImageHDU(np.full_like(flux, 4), name="SKY_IVAR"),
            fits.ImageHDU(mask, name="MASK"),
            fits.ImageHDU(np.ones_like(flux), name="LSF"),
            fits.BinTableHDU(slitmap, name="SLITMAP"),
        ]
    ).writeto(path)


def test_manifest_masks_and_cli_defaults(tmp_path: Path) -> None:
    sframe = tmp_path / "lvmSFrame-00000001.fits"
    make_sframe(sframe)
    manifest = tmp_path / "sframes.txt"
    manifest.write_text(f"{sframe.name}\n", encoding="utf-8")
    assert read_manifest(manifest) == [(0, sframe)]
    data = Table(
        {
            "n_gaia_sources": [0, 1, 2],
            "gaia_g_n_valid_sources": [0, 1, 1],
            "gaia_g_flux": [0.0, 0.5, 0.5],
            "lvm_flux_g_err": [1.0, 1.0, 1.0],
            "lvm_flux_g_valid_fraction": [1.0, 1.0, 1.0],
        }
    )
    assert gaia_clean_mask(data, 1).tolist() == [True, True, False]
    command = get_command(app)
    assert set(command.commands) == {
        "scan",
        "fetch-gaia",
        "combine-gaia",
        "build-medians",
        "status",
    }
    fetch_params = {
        parameter.name: parameter.default for parameter in command.commands["fetch-gaia"].params
    }
    assert fetch_params["query_workers"] == 5
    build_params = {parameter.name for parameter in command.commands["build-medians"].params}
    assert "gaia_ratio_threshold" in build_params
    ratio_data = Table(
        {
            "n_gaia_sources": [0, 1, 1, 2],
            "gaia_g_n_valid_sources": [0, 1, 1, 1],
            "lvm_flux_plus_sky_g": [10.0] * 4,
            "lvm_flux_plus_sky_g_valid_fraction": [1.0] * 4,
            "ratio_gaia_to_lvm_flux_plus_sky_g": [np.nan, 0.1, 0.11, 0.05],
        }
    )
    assert gaia_ratio_clean_mask(ratio_data, 0.1).tolist() == [True, True, False, False]
    fibers, _ = derive_fiber_table(sframe, Table())
    for name in ("svo_g_photon_coverage", "svo_g_energy_coverage"):
        assert 0 < float(fibers[name][0]) <= 1

    assert fibers["lvm_flux_g_valid_fraction"][0] < fibers["lvm_flux_g_valid_fraction"][1]


def test_build_modes_and_combined_ivar(tmp_path: Path) -> None:
    sframe = tmp_path / "lvmSFrame-00000001.fits"
    make_sframe(sframe)
    manifest = tmp_path / "sframes.txt"
    manifest.write_text(f"{sframe}\n", encoding="utf-8")

    median = tmp_path / "median.fits"
    result = build_stack(manifest, median, workers=1)
    assert result["complete"]
    with fits.open(median, checksum=True) as hdul:
        assert len(hdul["META"].data) == 1
        assert {
            "pwv_med",
            "sci_airmass",
            "skye_airmass",
            "skyw_airmass",
            "sky_near_label",
            "sky_far_label",
            "moon_fli",
            "fibers_sci_used",
        } <= set(hdul["META"].columns.names)
        assert hdul["INPUT_STATUS"].data[0]["status"] == "OK"
        assert all(hdu.verify_checksum() == 1 for hdu in hdul if "CHECKSUM" in hdu.header)

    gaia_dir = tmp_path / "gaia"
    gaia_dir.mkdir()
    Table(
        {
            "expnum": [1] * 6,
            "fiberid": np.arange(1, 7),
            "n_gaia_sources": [1] * 6,
            "gaia_g_n_valid_sources": [1] * 6,
            "lvm_flux_plus_sky_g": [10.0] * 6,
            "lvm_flux_plus_sky_g_valid_fraction": [1.0] * 6,
            "ratio_gaia_to_lvm_flux_plus_sky_g": [0.1] * 6,
        }
    ).write(gaia_dir / "lvmGAIA-fibers-00000001.fits")
    ratio_output = tmp_path / "ratio.fits"
    result = build_stack(
        manifest,
        ratio_output,
        workers=1,
        gaia_dir=gaia_dir,
        gaia_ratio_threshold=0.1,
    )
    assert result["complete"]
    with fits.open(ratio_output) as hdul:
        assert hdul[0].header["GAIATHR"] == 0.1
        assert hdul["META"].data["gaia_ratio_threshold"][0] == 0.1

    faint = tmp_path / "faint.fits"
    result = build_stack(manifest, faint, mode="faint-fibers", workers=1, fibers_per_telescope=1)
    assert result["complete"]
    with fits.open(faint) as hdul:
        assert hdul["FLUX"].data.shape == (3, 8)
        np.testing.assert_allclose(hdul["IVAR"].data, 2.0)
        assert hdul[0].header["NFIBTEL"] == 1


def test_partial_build_records_bad_input(tmp_path: Path) -> None:
    good = tmp_path / "lvmSFrame-00000001.fits"
    bad = tmp_path / "lvmSFrame-00000002.fits"
    make_sframe(good)
    make_sframe(bad)
    with fits.open(bad, mode="update") as hdul:
        hdul["WAVE"].data += 1
    manifest = tmp_path / "sframes.txt"
    manifest.write_text(f"{good}\n{bad}\n", encoding="utf-8")

    output = tmp_path / "partial.fits"
    result = build_stack(manifest, output, workers=1)
    assert result["ok"] == 1 and result["error"] == 1 and not result["complete"]
    with fits.open(output) as hdul:
        assert not hdul[0].header["COMPLETE"]
        assert list(hdul["INPUT_STATUS"].data["status"]) == ["OK", "ERROR"]
