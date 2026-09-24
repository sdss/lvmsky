"""ML context geometry: separations are topocentric, the moon phase follows META MOON_PHASE."""

import numpy as np
import pytest
from astropy.coordinates import AltAz, SkyCoord, get_body, get_sun
from astropy.time import Time
import astropy.units as u

import sys
from pathlib import Path

# mlp_predictor imports `sky_decomp` as a top-level package, the way the
# notebooks run it from skysub/, so import it the same way here.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import data, inference  # noqa: E402


def _one_row(mjd, ra, dec):
    return inference.build_triplet_from_pointings(
        obstime_mjd=mjd, sci_ra=ra, sci_dec=dec,
        sky_near_ra=ra + 3.0, sky_near_dec=dec, sky_far_ra=ra - 8.0, sky_far_dec=dec,
        verbose=False)


@pytest.mark.parametrize("mjd,ra,dec", [(60500.20, 250.0, -30.0), (60650.05, 40.0, -60.0),
                                        (61100.30, 180.0, -10.0)])
def test_sun_and_moon_separations_are_real_angles(mjd, ra, dec):
    tr = _one_row(mjd, ra, dec)
    names = list(tr["ctx_names"])
    ctx = np.asarray(tr["ctx_sci"], float)[0]
    t = Time(mjd, format="mjd", scale="utc")
    p = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    # Independent references: the pointing moved into the Sun's GCRS frame for
    # the elongation; AltAz at LCO for the (parallax-sensitive) moon.
    elong = get_sun(t).separation(p).deg
    lco = data._lco_earth_location()
    aa = AltAz(obstime=t, location=lco)
    moon_sep = get_body("moon", t, location=lco).transform_to(aa).separation(p.transform_to(aa)).deg
    assert ctx[names.index("sun_sep")] == pytest.approx(elong, abs=0.05)
    assert ctx[names.index("moon_sep")] == pytest.approx(moon_sep, abs=0.05)
    # The pre-2026-09-24 bug made moon_sep == 180 - elongation (ICRS
    # separation with the bodies at the barycentre).  Guard it, where the two
    # quantities are far enough apart to tell.
    if abs(moon_sep - (180.0 - elong)) > 5.0:
        assert abs(ctx[names.index("moon_sep")] - (180.0 - elong)) > 1.0


# Known lunar phases (UTC), 2024: first quarter 04-15 19:13, full 04-23 23:49,
# last quarter 05-01 11:27.
@pytest.mark.parametrize("iso,expected", [("2024-04-15T19:13:00", 90.0),
                                          ("2024-04-23T23:49:00", 180.0),
                                          ("2024-05-01T11:27:00", 270.0)])
def test_moon_phase_follows_meta_convention(iso, expected):
    ph = float(inference._compute_moon_phase_deg(np.array([Time(iso).mjd]))[0])
    # Topocentric vs geocentric and elongation vs ecliptic-longitude definitions
    # keep this within a couple of degrees of the nominal event.
    assert ((ph - expected + 180.0) % 360.0) - 180.0 == pytest.approx(0.0, abs=3.0)


def test_moon_fli_consistent_with_phase():
    mjd = np.linspace(60400.0, 60430.0, 13)
    ph = inference._compute_moon_phase_deg(mjd)
    t = Time(mjd, format="mjd", scale="utc")
    elong = get_body("moon", t, location=data._lco_earth_location()).separation(get_sun(t)).deg
    np.testing.assert_allclose((1.0 - np.cos(np.deg2rad(ph))) / 2.0,
                               (1.0 - np.cos(np.deg2rad(elong))) / 2.0, atol=1e-6)
