"""Self-contained checks of reference weighting, FITS persistence and covariance."""
from dataclasses import asdict
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from analysis.lsf_epoch_coverage.build_lsf_references import dispersion_average, write_product, publish, ROOT
from skysub.sky_decomp.lsf_reference import load_reference_lsf, reference_covariance, validate_coefficients, select_reference_lsf, load_lsf_from_sframe
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig, _wave_fingerprint
from skysub.sky_decomp.lsf_spline2d import LSF_SPLINE2D_REPRESENTATION, evaluate_lsf_density, mspline_bin_integrals
from skysub.sky_decomp.result_io import STATE_CONFIG_COLUMNS
from analysis.lsf_epoch_coverage.compact_lsf_references import compact_reference


def fixtures():
    wave = np.linspace(3600, 9800, 12401).astype(np.float32).astype(float)
    fingerprint = _wave_fingerprint(wave)
    config = asdict(LSFSurfaceIterativeConfig(n_basis=4, n_refinement_cycles=30))
    configs, knots = [], {}
    for channel, lower, upper in (("B", 3600., 5787.), ("R", 5787., 7454.), ("Z", 7454., 9800.)):
        knots[channel] = np.repeat([lower, upper], 4)
        cf = {column: converter(config[key]) for key, column, converter in STATE_CONFIG_COLUMNS}
        cf.update(channel=channel, lower=np.nan if channel == "B" else lower,
                  upper=np.nan if channel == "Z" else upper, n_basis=4, degree=3, n_knots=8,
                  tap_offsets=np.arange(-5, 6), schema_version=1,
                  legacy_kernel_representation=LSF_SPLINE2D_REPRESENTATION,
                  wave_n=fingerprint[0], wave_min=fingerprint[1], wave_max=fingerprint[2], wave_sha256=fingerprint[3],
                  offset_half_width_angstrom=2.5, offset_roughness_fraction=.01)
        configs.append(cf)
    values, rows = [], []
    for ci, channel in enumerate(("B", "R", "Z")):
        for i, width in enumerate((1.2, 1.3, 1.4, 1.3)):
            masses = np.exp(-.5 * (np.arange(-5, 6) / width)**2)
            masses /= masses.sum()
            values.append(np.tile(masses[:, None], (1, 4)))
            rows.append(dict(sample_id=len(rows), expnum=10+i, fiberid=123, channel=channel, telescope="Sci",
                             spectrographid=1, calibration_epoch=60000, full_moon_utc="2023-04-06T00:00:00",
                             phase_bin="waxing_99-100%", mjd=60001 + i//2, observation_mjd_utc=60001+i*.5,
                             source_id=0, fit_row=i, gaia_ratio=.005, gaia_clean_1pct=True, gaia_clean_10pct=True,
                             fwhm_angstrom=1.5 + .01*i, w50_angstrom=1.0 + .01*i,
                             w90_angstrom=2.8 + .01*i, centroid_angstrom=0.0))
    return pd.DataFrame(rows), np.asarray(values), Table(rows=configs), wave, knots


class ReferenceTests(unittest.TestCase):
    def test_constraints_balance_and_small_samples(self):
        _, x, _, _, _ = fixtures()
        same = np.repeat(x[:1], 3, axis=0)
        result = dispersion_average(same, [1, 1, 2])
        np.testing.assert_allclose(result["weights"], [.25, .25, .5])
        np.testing.assert_allclose(result["coefficient"], x[0])
        validate_coefficients(result["coefficient"], tolerance=1e-12)
        single = dispersion_average(x[:1], [1])
        self.assertTrue(np.isnan(single["scatter"]).all())
        self.assertTrue(np.isnan(single["conditional_error"]).all())
        self.assertEqual(single["weights"][0], 1)
        pair = dispersion_average(x[:2], [1, 1])
        np.testing.assert_array_equal(pair["robust"], [1, 1])
        invalid = x[0].copy()
        invalid[0, 0] = 2
        with self.assertRaises(ValueError):
            validate_coefficients(invalid)

    def test_outlier_downweight_and_covariance_constraints(self):
        _, x, _, _, _ = fixtures()
        samples = np.concatenate([x[:4], x[:4], x[:4], np.full((1, 11, 4), 1/11)])
        result = dispersion_average(samples, np.arange(len(samples)))
        self.assertLess(result["weights"][-1], np.median(result["weights"][:-1]))
        np.testing.assert_allclose(result["weights"].sum(), 1)
        covariance = result["factors"].T @ result["factors"] / result["denominator"]
        np.testing.assert_allclose(np.diag(covariance), result["scatter"].ravel()**2)
        direction = np.zeros(44)
        direction[::4] = 1
        np.testing.assert_allclose(covariance @ direction, 0, atol=1e-15)
        self.assertGreaterEqual(np.linalg.eigvalsh(covariance).min(), -1e-14)

    def test_nearest_lunation_and_sframe_provenance(self):
        meta, x, config, wave, knots = fixtures()
        manifest = dict(fingerprint="test", source_fits=[dict(path="source.fits", size=1, mtime_ns=1, sha256="test")])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "publication"
            root.mkdir()
            first = write_product("lunation", (60000, "2023-04-06T00:00:00"), meta, x, config, wave, knots, manifest, root)
            later = meta.copy()
            later["observation_mjd_utc"] += 30
            later["full_moon_utc"] = "2023-05-06T00:00:00"
            second = write_product("lunation", (60000, "2023-05-06T00:00:00"), later, x, config, wave, knots, manifest, root)
            epochs = dict.fromkeys(("B", "R", "Z"), 60000)
            for date, catalog in ((60001., first), (60031., second)):
                state = select_reference_lsf(root, 123, epochs, date)
                for channel in epochs:
                    self.assertEqual(state.metrics[channel]["reference_file"], str(root/catalog.file.iloc[0]))
            with self.assertRaisesRegex(ValueError, "No B LSF"):
                select_reference_lsf(root, 123, dict.fromkeys(epochs, 61000), 60031.)
            frame = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(wave, name="WAVE"),
                                 fits.BinTableHDU(Table(rows=[dict(fiberid=123, spectrographid=1)]), name="SLITMAP")])
            frame[0].header["OBSTIME"] = "2023-03-28T00:00:00"
            for channel in "brz":
                frame[0].header.add_comment(f"/calib/60000/lvm-mlsf-{channel}1.fits")
            path = Path(directory)/"target.fits"
            frame.writeto(path)
            state = load_lsf_from_sframe(path, 123, root)
            self.assertEqual(state.metrics["B"]["reference_file"], str(root/second.file.iloc[0]))
            survey = write_product("survey", (), meta, x, config, wave, knots, manifest, root)
            with fits.open(path, mode="update") as hdus:
                del hdus[0].header["COMMENT"]
                del hdus[0].header["OBSTIME"]
            pooled = load_lsf_from_sframe(path, 123, root, scope="survey")
            expected = load_reference_lsf(root/survey.file.iloc[0], 123, allow_epoch_pooling=True)
            for channel in epochs:
                np.testing.assert_array_equal(pooled.coefficients[channel], expected.coefficients[channel])
                self.assertEqual(pooled.metrics[channel]["reference_file"], str(root/survey.file.iloc[0]))
            with self.assertRaisesRegex(ValueError, "calibration files"):
                load_lsf_from_sframe(path, 123, root)
            with self.assertRaisesRegex(ValueError, "scope"):
                load_lsf_from_sframe(path, 123, root, scope="invalid")
            with self.assertRaisesRegex(ValueError, "physical fiber"):
                load_lsf_from_sframe(path, 999, root, scope="survey")
            with fits.open(path, mode="update") as hdus:
                hdus["WAVE"].data[0] += .01
            with self.assertRaisesRegex(ValueError, "wavelength grid"):
                load_lsf_from_sframe(path, 123, root, scope="survey")

    def test_fits_round_trip_factors_and_guards(self):
        meta, x, config, wave, knots = fixtures()
        manifest = dict(fingerprint="test", source_fits=[dict(path="source.fits", size=1, mtime_ns=1, sha256="test")])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "publication"
            root.mkdir()
            catalog = write_product("epoch", (60000,), meta, x, config, wave, knots, manifest, root)
            path = root / catalog.file.iloc[0]
            state = load_reference_lsf(path, 123)
            for ci, channel in enumerate(("B", "R", "Z")):
                normalized = x[ci*4:ci*4+4].copy()
                normalized /= normalized.sum(axis=1, keepdims=True)
                expected = dispersion_average(normalized, meta.mjd.iloc[ci*4:ci*4+4].to_numpy())
                np.testing.assert_array_equal(state.coefficients[channel], expected["coefficient"])
                covariance = reference_covariance(path, ci)
                np.testing.assert_allclose(np.diag(covariance), expected["scatter"].ravel()**2, rtol=2e-7, atol=1e-15)
                conditional = reference_covariance(path, ci, kind="conditional_mean")
                np.testing.assert_allclose(np.diag(conditional), expected["conditional_error"].ravel()**2, rtol=2e-7, atol=1e-15)
                bins = mspline_bin_integrals(np.array([-2.5]), np.array([2.5]), n_basis=11, half_width=2.5)[0]
                np.testing.assert_allclose(bins @ state.coefficients[channel], 1, atol=1e-12)
            self.assertEqual(state.wave_sha256, _wave_fingerprint(wave)[3])
            self.assertTrue(np.isfinite(evaluate_lsf_density(state, [4700., 6600., 8500.], [-1., 0., 1.])).all())
            with self.assertRaises(ValueError):
                load_reference_lsf(path, 999)
            with fits.open(path, mode="update") as hdus:
                hdus[0].header["APPLYOK"] = False
            with self.assertRaises(ValueError):
                load_reference_lsf(path, 123)
            load_reference_lsf(path, 123, allow_epoch_pooling=True)
            paths = {}
            with fits.open(path) as original:
                for ci, channel in enumerate(("B", "R", "Z")):
                    selected = fits.HDUList([h.copy() for h in original])
                    selected[0].header["APPLYOK"] = True
                    selected["REFERENCES"].data = selected["REFERENCES"].data[[ci]]
                    selected["LSF_COEF"].data = selected["LSF_COEF"].data[[ci]]
                    selected["LSF_KNOTS"].data = selected["LSF_KNOTS"].data[[ci]]
                    paths[channel] = root / f"sparse_{channel}.fits"
                    selected.writeto(paths[channel], checksum=True)
            with self.assertRaises(ValueError):
                load_reference_lsf(paths["B"], 123)
            assembled = load_reference_lsf(paths, 123)
            for channel in ("B", "R", "Z"):
                np.testing.assert_array_equal(assembled.coefficients[channel], state.coefficients[channel])
                self.assertEqual(assembled.metrics[channel]["reference_file"], str(paths[channel]))
            old_covariance = reference_covariance(path, 0)
            compact_reference(path, root)
            compact = load_reference_lsf(path, 123, allow_epoch_pooling=True)
            for channel in ("B", "R", "Z"):
                np.testing.assert_array_equal(compact.coefficients[channel], state.coefficients[channel])
            np.testing.assert_array_equal(reference_covariance(path, 0), old_covariance)
            self.assertLess(path.stat().st_size, (root.parent / "lsf_reference_audit" / path.relative_to(root)).stat().st_size)


if __name__ == "__main__":
    unittest.main()
