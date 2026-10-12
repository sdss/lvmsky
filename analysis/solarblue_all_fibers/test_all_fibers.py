import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from analysis.solarblue_all_fibers.run_all_fibers import (
    attach_input_context,
    compact_product,
    read_selected_manifest,
    selected_fibers,
    write_immutable,
)


def test_all_science_and_sky_fibers_are_kept_with_nine_references_first():
    rows = []
    for spectrograph in (1, 2, 3):
        for telescope, count in (("Sci", 10), ("SkyE", 2), ("SkyW", 2), ("Spec", 1)):
            for position in range(count):
                rows.append((len(rows), spectrograph, telescope, position, 0))
    table = Table(rows=rows, names=("fiberid", "spectrographid", "telescope", "ypix_b", "fibstatus"))
    indices, reference = selected_fibers(table)
    assert len(reference) == 9
    assert np.array_equal(indices[:9], reference)
    assert len(indices) == 42
    assert set(np.asarray(table[indices]["telescope"])) == {"Sci", "SkyE", "SkyW"}


def test_selection_artifacts_are_immutable(tmp_path):
    path = tmp_path / "selected.csv"
    write_immutable(path, "expnum\n1\n")
    write_immutable(path, "expnum\n1\n")
    with pytest.raises(RuntimeError, match="Immutable selection artifact differs"):
        write_immutable(path, "expnum\n2\n")


def test_frozen_selection_is_validated(tmp_path):
    path = tmp_path / "selected.csv"
    path.write_text(
        "expnum,mjd,moon_alt,moon_fli,source_path,source_available\n"
        "1,60000,61,.96,/tmp/frame.fits,True\n"
    )
    assert read_selected_manifest(path, 60, .95).expnum.tolist() == [1]
    with pytest.raises(ValueError, match="ineligible"):
        read_selected_manifest(path, 62, .95)


def test_fit_product_gets_self_contained_input_context(tmp_path):
    source = tmp_path / "input.fits"
    product = tmp_path / "fit.fits"
    fits.HDUList([
        fits.PrimaryHDU(), fits.ImageHDU(np.arange(3.0), name="WAVE"),
        fits.ImageHDU(np.ones((2, 3)), name="IVAR_SCI"),
        fits.ImageHDU(np.zeros((2, 3), dtype=np.uint8), name="MASK_SCI"),
        fits.BinTableHDU(Table({"expnum": [1, 1], "fiberid": [10, 11]}), name="META"),
    ]).writeto(source)
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.ones((2, 3)), name="BESTFIT")]).writeto(product)
    attach_input_context(source, product, rows=1, include_wave=True)
    with fits.open(product) as hdul:
        assert np.array_equal(hdul["WAVE"].data, np.arange(3.0))
        assert len(hdul["INPUT_META"].data) == 1
        assert hdul["INPUT_META"].data["fiberid"][0] == 10
        assert hdul["INPUT_IVAR"].data.shape == (1, 3)
        assert hdul["INPUT_MASK"].data.dtype == np.uint8


def test_compact_product_keeps_only_requested_float32_arrays(tmp_path):
    product = tmp_path / "fit.fits"
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(np.ones((2, 3), dtype=np.float64), name="BESTFIT_LSF"),
        fits.ImageHDU(np.ones((2, 3), dtype=np.float64), name="BESTFIT"),
        fits.BinTableHDU(Table({"a": np.ones(2, dtype=np.float64)}), name="COEF"),
        fits.BinTableHDU(Table({"a": np.ones(2, dtype=np.float64)}), name="COEF_ERR"),
        fits.ImageHDU(np.ones((2, 2, 2), dtype=np.float64), name="COEF_COV_MOON"),
    ]).writeto(product)
    compact_product(product, reference=True)
    with fits.open(product) as hdul:
        assert hdul[0].header["STORFMT"] == "minimal-f32-v2"
        assert "BESTFIT" not in hdul and "COEF_ERR" not in hdul and "COEF_COV_MOON" not in hdul
        assert hdul["BESTFIT_LSF"].data.dtype.itemsize == 4
        assert hdul["COEF"].data["a"].dtype.itemsize == 4


def test_compact_product_stores_lsf_config_once_and_only_failed_text(tmp_path):
    product = tmp_path / "fit.fits"
    lsf_meta = Table({
        "spectrum_index": [0, 0, 0], "channel": ["B", "R", "Z"],
        "available": [True, True, True], "status": ["Solved", "fallback", "Solved"],
        "fit_status": ["Solved", "failed", "Solved"], "failure_reason": ["", "bad", ""],
        "lower": [3600.0, 5700.0, 7500.0], "degree": [3, 3, 3],
    })
    fits.HDUList([
        fits.PrimaryHDU(), fits.BinTableHDU(lsf_meta, name="LSF_META"),
    ]).writeto(product)
    compact_product(product)
    with fits.open(product) as hdul:
        assert len(hdul["LSF_CONFIG"].data) == 3
        assert len(hdul["LSF_FAILURES"].data) == 1
        assert set(hdul["LSF_META"].columns.names) == {
            "spectrum_index", "channel", "available",
        }
