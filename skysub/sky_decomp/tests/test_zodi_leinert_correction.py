"""The empirical Leinert zodi correction and the decomposition/cache pairing."""

import numpy as np
import pytest
from astropy.io import fits

from skysub.sky_decomp.moon_zodi_model import (
    ZODI_LEINERT_CORRECTIONS,
    zodi_leinert_correction_log10,
)

NAME = "lvm-ecl-2026-09"


def _reference(dlam, beta):
    """Independent evaluation of the stored polynomial."""
    spec = ZODI_LEINERT_CORRECTIONS[NAME]
    u = 2.0 * dlam / 180.0 - 1.0
    v = 2.0 * np.sin(np.deg2rad(abs(beta))) - 1.0
    pu = np.polynomial.legendre.legvander([u], 2)[0]
    pv = np.polynomial.legendre.legvander([v], 2)[0]
    return sum(c * pu[i] * pv[j] for (i, j), c in zip(spec["terms"], spec["coef"]))


def test_none_is_exactly_zero():
    for lon, lat in ((0.0, 0.0), (123.4, -56.7), (-170.0, 89.0)):
        assert zodi_leinert_correction_log10("none", lon, lat) == 0.0


def test_unknown_name_raises_rather_than_defaulting():
    with pytest.raises(KeyError):
        zodi_leinert_correction_log10("no-such-correction", 90.0, 10.0)


def test_longitude_folds_to_absolute_sun_relative_value():
    # The model returns relative longitude wrapped to [-180, 180]; the fit used
    # |lambda - lambda_sun| in [0, 180].
    a = zodi_leinert_correction_log10(NAME, 160.0, 20.0)
    for lon in (-160.0, 200.0, 520.0, -200.0):
        assert zodi_leinert_correction_log10(NAME, lon, 20.0) == pytest.approx(a, abs=1e-12)


def test_symmetric_in_ecliptic_latitude():
    assert (zodi_leinert_correction_log10(NAME, 100.0, -35.0)
            == pytest.approx(zodi_leinert_correction_log10(NAME, 100.0, 35.0), abs=1e-12))


@pytest.mark.parametrize("dlam,beta", [(75.0, 0.0), (120.0, 0.0), (120.0, 90.0), (180.0, 45.0)])
def test_matches_independent_evaluation(dlam, beta):
    assert zodi_leinert_correction_log10(NAME, dlam, beta) == pytest.approx(
        _reference(dlam, beta), abs=1e-12)


def test_measured_shape_pole_brighter_than_ecliptic():
    # The fitted finding: the Leinert profile is too steep in latitude, so the
    # correction raises the poles relative to the ecliptic (+0.30 dex at 120).
    rise = (zodi_leinert_correction_log10(NAME, 120.0, 90.0)
            - zodi_leinert_correction_log10(NAME, 120.0, 0.0))
    assert 0.25 < rise < 0.35


def _write_decomp(path, tag):
    prim = fits.PrimaryHDU()
    if tag is not None:
        prim.header["ZODICORR"] = tag
    fits.HDUList([prim]).writeto(path)


def test_cache_reads_the_decomposition_tag(tmp_path):
    from skysub.mlp_predictor.moon_model_cache import decomposition_zodi_correction
    prefix = str(tmp_path / "stack")
    assert decomposition_zodi_correction(prefix) == "none"          # no products
    _write_decomp(f"{prefix}_decomp_sci_x.fits", None)
    assert decomposition_zodi_correction(prefix) == "none"          # pre-2026-09-24
    _write_decomp(f"{prefix}_decomp_sky1_x.fits", None)
    assert decomposition_zodi_correction(prefix) == "none"


def test_cache_refuses_disagreeing_products(tmp_path):
    from skysub.mlp_predictor.moon_model_cache import decomposition_zodi_correction
    prefix = str(tmp_path / "stack")
    _write_decomp(f"{prefix}_decomp_sci_x.fits", NAME)
    assert decomposition_zodi_correction(prefix) == NAME
    _write_decomp(f"{prefix}_decomp_sky1_x.fits", "none")
    with pytest.raises(RuntimeError, match="disagree"):
        decomposition_zodi_correction(prefix)
