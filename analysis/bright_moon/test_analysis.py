"""Run with the same Python environment as analyze.py."""
import numpy as np
import tempfile
from pathlib import Path
from types import SimpleNamespace
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import analyze
import pandas as pd
from analyze import full_moons, geometry


def test_full_moons():
    # Independent USNO reference, rounded to one minute:
    # https://aa.usno.navy.mil/calculated/moon/phases?date=2024-01-21&format=p&nump=5&submit=Get+Data
    actual = full_moons(Time('2024-01-01').mjd, Time('2024-02-29').mjd)
    expected = Time(['2024-01-25T17:54:00', '2024-02-24T12:30:00']).mjd
    assert len(actual) == 2
    assert np.max(np.abs(actual - expected)) * 86400 < 60
    alt, sun = geometry(Time(['2024-01-25T04:00:00', '2024-01-25T16:00:00']).mjd)
    assert alt[0] > 0 and sun[0] < -18
    assert sun[1] > 0


def test_continuum():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'synthetic.fits'
        flux = np.ones((20, 201), dtype=np.float32) * 1e-15
        flux[0] = 1e-12  # A stellar outlier should not dominate the median.
        mask = np.zeros_like(flux, dtype=np.uint8)
        mask[1] = 1
        fibers = Table({'telescope': ['SkyE'] * 10 + ['SkyW'] * 10,
                        'fibstatus': [0, 0, 1] + [0] * 17})
        hdus = [fits.PrimaryHDU(), fits.ImageHDU(flux, name='FLUX'),
                fits.ImageHDU(np.ones_like(flux) * 2e-15, name='SKY'),
                fits.ImageHDU(np.ones_like(flux), name='IVAR'),
                fits.ImageHDU(mask, name='MASK'),
                fits.ImageHDU(np.linspace(5400, 5500, 201), name='WAVE'),
                fits.BinTableHDU(fibers, name='SLITMAP')]
        hdus[1].header['BUNIT'] = 'erg / (Angstrom s cm2)'
        fits.HDUList(hdus).writeto(path)
        original = analyze.SAS
        try:
            analyze.SAS = Path(directory)
            result = analyze.measure(SimpleNamespace(expnum=1, location=path.name))
        finally:
            analyze.SAS = original
        assert result['measurement_status'] == 'ok'
        assert result['skye_n_fibers'] == 8 and result['skyw_n_fibers'] == 10
        assert np.isclose(result['skye_continuum'] / 1e-15, 3)
        assert np.isclose(result['skyw_continuum'] / 1e-15, 3)


def test_report_tables():
    root = Path(__file__).resolve().parent
    all_rows = pd.read_csv(root / 'all_exposures.csv')
    selected = pd.read_csv(root / 'selected_exposures.csv')
    winners = pd.read_csv(root / 'one_brightest_per_cycle.csv')
    outliers = pd.read_csv(root / 'calibration_outliers.csv')
    assert len(all_rows) == 17260 and all_rows.expnum.is_unique
    assert len(selected) == 1765 and selected.selected.all()
    assert set(outliers.expnum) == {23816}
    assert outliers.sky_continuum_raw.iloc[0] > 0.1 and outliers.sky_continuum.isna().all()
    assert 23816 not in winners.expnum.to_numpy()
    assert len(winners) == 26
    assert np.allclose(winners.sky_continuum,
                       selected.groupby('full_moon_utc').sky_continuum.max().to_numpy())


if __name__ == '__main__':
    test_full_moons()
    test_continuum()
    if (Path(__file__).resolve().parent / 'calibration_outliers.csv').exists():
        test_report_tables()
    print('USNO epochs, LCO night/day, and masked sky-continuum checks passed.')
