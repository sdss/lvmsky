"""Round-trip heterogeneous spline histories without spectral image extensions."""

import numpy as np
from astropy.io import fits
from scipy.interpolate import BSpline

from skysub.sky_subtraction.product import history_table


def test_rz_spline_parameters_survive_fits_roundtrip(tmp_path):
    records = []
    expected = []
    for channel, lo, hi in [('R', 5787.0, 7454.0), ('Z', 7454.0, 9800.5)]:
        wave = np.arange(lo, hi, 0.5)
        knots = np.r_[np.repeat(wave[0], 4), np.arange(wave[0] + 20, wave[-1], 20), np.repeat(wave[-1], 4)]
        coef = np.linspace(0.2, 3, len(knots) - 4)
        expected.append((wave, BSpline(knots, coef, 3)(wave)))
        records.append(
            dict(channel=channel, continuum_knots=knots, continuum_coef=coef, variance_continuum_coef=coef * 0.9)
        )
    path = tmp_path / 'history.fits'
    fits.BinTableHDU(history_table(records)).writeto(path)
    with fits.open(path) as h:
        for row, (wave, truth) in zip(h[1].data, expected, strict=True):
            knots = row['continuum_knots'][: row['continuum_knots_n']]
            coef = row['continuum_coef'][: row['continuum_coef_n']]
            np.testing.assert_array_equal(BSpline(knots, coef, 3)(wave), truth)


def test_survey_only_diagnostic_has_empty_fit_without_history():
    from astropy.table import Table

    from skysub.sky_subtraction.reconstruction import SkyMLReconstructor

    r = SkyMLReconstructor.__new__(SkyMLReconstructor)
    r.wave = np.arange(4.0)
    r.h = fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(Table(dict(fiberid=[10])), name='SLITMAP')])
    r._observed_input = (np.ones((1, 4)) * 3, np.ones((1, 4)), np.zeros((1, 4), dtype='u1'))
    r.sky = lambda fiberid: dict(SKY=np.ones(4))
    result = r.fiber(10)
    np.testing.assert_array_equal(result['FLUX'], np.ones(4) * 2)
    assert np.all(np.isnan(result['FIT_MODEL'])) and np.all(result['FIT_MASK'] == 1)


def test_reconstructor_identifies_an_incompatible_product(tmp_path):
    import pytest

    from skysub.sky_subtraction.reconstruction import SkyMLReconstructor

    path = tmp_path / 'ordinary.fits'
    fits.PrimaryHDU().writeto(path)
    with pytest.raises(ValueError, match='coefficient-based'):
        SkyMLReconstructor(path)


def test_survey_reconstruction_uses_recorded_calibration(tmp_path):
    import json

    import pytest

    from skysub.sky_subtraction.io import BUNDLED_CALIBRATION
    from skysub.sky_subtraction.reconstruction import SkyMLReconstructor

    manifest = json.loads((BUNDLED_CALIBRATION / 'manifest.json').read_text())
    r = SkyMLReconstructor.__new__(SkyMLReconstructor)
    r.provenance = dict(calibration=str(BUNDLED_CALIBRATION), calibration_manifest=manifest)
    assert r.calibration == BUNDLED_CALIBRATION
    other = tmp_path / 'lsf'
    other.mkdir()
    (other / 'manifest.json').write_text(json.dumps(dict(manifest, note='different')))
    r = SkyMLReconstructor.__new__(SkyMLReconstructor)
    r.provenance = dict(calibration=str(other), calibration_manifest=manifest)
    with pytest.raises(ValueError, match='differs'):
        _ = r.calibration
