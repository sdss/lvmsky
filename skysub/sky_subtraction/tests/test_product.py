"""FITS publication, protected destinations, fiber isolation and correction consistency."""

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
from skysub.sky_subtraction import pipeline
from skysub.sky_subtraction.decomp_adapter import MedianFits
from skysub.sky_subtraction.io import Median
from skysub.sky_subtraction.product import MLCoefficients, ProductInputs, build_product, sha256, write_product
from skysub.sky_subtraction.rendering import render

LSF = Path(__file__).parents[2] / 'sky_decomp/data/lsf'


def test_output_is_reconstructible_and_atomic(tmp_path):
    wave = np.arange(3600.0, 9800.5, 0.5)
    ref = load_survey_ascii_lsf(LSF, 500)
    refs = [ref, deepcopy(ref)]
    n = len(wave)
    flux = np.full((3, n), 1e-14, dtype='f4')
    first = fits.PrimaryHDU()
    first.header['EXPOSURE'] = 1
    image = fits.ImageHDU(flux, name='FLUX')
    image.header['BUNIT'] = 'erg / (Angstrom s cm2)'
    source = tmp_path / 'lvmCFrame-00000001.fits'
    slit = Table(dict(fiberid=[500, 501, 502], telescope=['Sci', 'Sci', 'SkyE']))
    fits.HDUList(
        [
            first,
            image,
            fits.ImageHDU(np.full_like(flux, 1e30), name='IVAR'),
            fits.ImageHDU(np.zeros_like(flux, dtype='u1'), name='MASK'),
            fits.ImageHDU(wave, name='WAVE'),
            fits.ImageHDU(np.ones_like(flux), name='LSF'),
            fits.BinTableHDU(slit, name='SLITMAP'),
        ]
    ).writeto(source)
    raw = np.full_like(flux, 1e-15)
    sky = raw + raw / 2 + raw / 4
    sky[0, wave > 8000] = np.nan
    sky[2] = 3e-15
    sky[2, 5] = 1e-14
    results = {
        k: SimpleNamespace(lsf_state=ref, coef=np.ones(2), fit_status='Solved', reduced_chi2=1.0, t_o2=190.0)
        for k in ('Sci', 'SkyE', 'SkyW')
    }
    median_fits = MedianFits(
        labels=dict(sci='Sci', sky1='SkyE', sky2='SkyW'),
        results=results,
        reliability={},
        science_model=None,
        fit_mask=np.zeros(n, dtype='u1'),
        render_context=dict(settings={}, models=[dict(telescope=k) for k in results]),
        transmission=np.ones((3, n)),
    )
    inputs = ProductInputs(
        sky=sky,
        rows=np.array([0, 1]),
        arm_rows=np.array([2]),
        arm_states=[ref],
        adopted=refs,
        candidates=refs,
        qa=[dict(channel=ch, accepted=True, reason='insitu', fiberid=f) for f in [500, 501] for ch in 'RZ'],
        history=[],
        photometry=slit,
        sources=None,
        medians={k: Median(flux[0], np.ones(n), np.full(n, 2)) for k in results},
        median_fits=median_fits,
        initial_states=[ref, ref, ref],
        stack_meta=Table({'expnum': [1]}),
        ml=MLCoefficients(['a', 'b'], np.ones(2), np.ones(2) * 0.1, np.ones(2)),
        fit_mask=np.zeros(n, dtype='u1'),
    )
    dest = tmp_path / 'lvmSFrame-skyml-00000001.fits'
    before = sha256(source)
    write_product(source, dest, inputs, {})
    assert sha256(source) == before
    with fits.open(dest, checksum=True) as h:
        np.testing.assert_allclose(flux - h['SKY'].data, h['FLUX'].data, rtol=1e-6, equal_nan=True)
        assert all(
            n not in h
            for n in [
                'OBSFLUX',
                'SKY_IVAR',
                'SKY_RAW',
                'SKY_RES',
                'SKY_OH',
                'SKY_EPI',
                'SKY_SURVEY',
                'SKYML_VALID',
                'FIT_MODEL',
                'FIT_SKY',
                'FIT_CONT',
            ]
        )
        assert 'GAIA' in h and 'MED_QA' in h and 'CAND_LSF_COEF' in h and 'INIT_LSF_COEF' in h
        assert h['SKY'].header['BUNIT'] == h['FLUX'].header['BUNIT']
        assert np.all(h['IVAR'].data[0, wave > 8000] == 0)
        assert np.all(h['MASK'].data[0, wave > 8000] != 0)
        np.testing.assert_array_equal(h['SLITMAP'].data['skyml_applied'], [True, True, True])
        assert h['FLUX'].data[2, 5] == 0 and h['SKY'].data[2, 5] == np.float32(1e-14)
    # The in-memory product (used by the block method) equals what is published.
    memory = build_product(source, inputs, {})
    with fits.open(dest) as h:
        assert sorted(x.name for x in memory) == sorted(x.name for x in h)
        for name in ('FLUX', 'SKY', 'MED_COEF'):
            np.testing.assert_array_equal(memory[name].data, h[name].data)
        assert memory['ML_COEF'].data['coef'].dtype.kind == 'f' and memory['ML_COEF'].data['coef'].itemsize == 4
    existing = sha256(dest)
    with pytest.raises(FileExistsError):
        write_product(source, dest, inputs, {})
    assert sha256(dest) == existing
    assert not list(tmp_path.glob('.skyml-*.fits'))


def test_protected_destination_stops_before_input_read(tmp_path):
    output = tmp_path / 'lvmSFrame-skyml-00000001.fits'
    output.touch()
    with pytest.raises(FileExistsError):
        pipeline.run(tmp_path / 'lvmCFrame-00000001.fits', 'unused', output=output)
    with pytest.raises(ValueError, match='protect'):
        pipeline.run(tmp_path / 'lvmCFrame-00000001.fits', 'unused', output=tmp_path / 'lvmSFrame-00000001.fits')


def _fiber_work(reference, wave, survey_only=False):
    n = len(wave)
    return pipeline.FiberWork(
        model=None,
        wave=wave,
        flux=np.ones((1, n)),
        ivar=np.ones((1, n)),
        mask=np.zeros((1, n)),
        fiberids=np.array([500]),
        references={0: reference},
        comparison=reference,
        designs={},
        science_mask=np.zeros(n, bool),
        coef=None,
        oh_coef=None,
        arms=[],
        survey_only=survey_only,
        settings=pipeline.SETTINGS,
    )


def test_failing_fiber_fit_keeps_survey_lsf(monkeypatch):
    reference = load_survey_ascii_lsf(LSF, 500)
    reference.metrics['Z'] = dict(reference.metrics['Z'], quality='warning')
    wave = np.arange(3600.0, 9800.5, 0.5)

    def broken(*args):
        raise RuntimeError('solver crashed')

    monkeypatch.setattr(pipeline, 'fit_fiber', broken)
    monkeypatch.setattr(pipeline, 'sky_parts', lambda *args: SimpleNamespace(total=np.zeros(len(wave))))
    monkeypatch.setattr(pipeline, '_WORK', _fiber_work(reference, wave))
    result = pipeline._fiber_task(0)
    reasons = {r['channel']: r['reason'] for r in result.fit.qa}
    assert reasons == {'R': 'fiber_error:RuntimeError', 'Z': 'unusable_survey_fallback'}
    assert all(r['fiberid'] == 500 and r['input_row'] == 0 for r in result.fit.qa)
    np.testing.assert_array_equal(result.fit.adopted.coefficients['R'], reference.coefficients['R'])
    assert np.isnan(result.sky[wave >= 7454.0]).all() and np.all(result.sky[wave < 7454.0] == 0)


def _small_model():
    from skysub.sky_decomp.lsf_spline2d import SkyDecompLSFSpline2D
    from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig

    data = Path(__file__).parents[2] / 'sky_decomp'
    wave = np.load(data / 'tests/data/lsf_surface_iterative_row837_n5_golden.npz')['wave'].astype(float)
    return SkyDecompLSFSpline2D(
        wave,
        lsf_sigma=0.7,
        base_dir=data / 'data',
        n_spline_knots=9,
        n_zodi_spline_knots=3,
        offset_half_width_angstrom=2.5,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=1, n_basis=4),
    )


def test_full_sky_lsf_tensor_matches_renderer_when_initial_profiles_differ():
    from skysub.sky_subtraction.fiber_lsf import sky_designs

    ref = load_survey_ascii_lsf(LSF, 500)
    model = _small_model()
    coefficient = np.full(len(model.design_names), 0.001)
    _, _, matrices = render(model, coefficient, ref)
    other = deepcopy(ref)
    for ch in 'RZ':
        c = other.coefficients[ch]
        other.coefficients[ch] = 0.8 * c + 0.1 * np.roll(c, 1, axis=0) + 0.1 * np.roll(c, -1, axis=0)
    render(model, coefficient, other)
    designs = sky_designs(model, coefficient, ref, matrices)
    expected, _, _ = render(model, coefficient, ref)
    for ch, (take, design) in designs.items():
        actual = (design['lines'] + design['solar']) @ ref.coefficients[ch].ravel() + design['fixed']
        np.testing.assert_allclose(actual, expected[take], rtol=1e-10, atol=1e-12)


def test_oh_coefficient_correction_equals_template_correction():
    """Rendering corrected coefficients must add exactly the fitted template corrections."""
    ref = load_survey_ascii_lsf(LSF, 500)
    model = _small_model()
    rng = np.random.default_rng(3)
    coefficient = rng.uniform(0.001, 0.01, len(model.design_names))
    _, components, matrices = render(model, coefficient, ref)
    templates = pipeline.line_templates(model, matrices, coefficient, atoms=False, o2=False)
    parameters = {
        name: 1 + rng.uniform(-0.01, 0.01) if name.endswith('_tilt') else 1 + rng.uniform(-0.2, 0.2)
        for name in templates
    }
    corrected = pipeline.corrected_coefficients(model, matrices, coefficient, parameters)
    actual = model._components_from_coef(corrected, matrices)['oh'] - components['oh']
    expected = sum((parameters[name] - 1) * template for name, template in templates.items())
    scale = np.max(np.abs(expected))
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12 * scale)


def test_explicit_overwrite_replaces_only_after_successful_write(tmp_path):
    from skysub.sky_subtraction.product import publish
    output = tmp_path / "result.fits"
    original = b"previous result"
    output.write_bytes(original)

    class NewProduct:
        def writeto(self, path, **kwargs):
            assert output.read_bytes() == original
            fits.HDUList([fits.PrimaryHDU(np.arange(3.0))]).writeto(path, **kwargs)
            assert output.read_bytes() == original

    publish(NewProduct(), output, 0o644, overwrite=True)
    with fits.open(output) as hdus:
        np.testing.assert_array_equal(hdus[0].data, np.arange(3.0))
        assert hdus[0].verify_checksum() == 1
    assert not list(tmp_path.glob(".skyml-*.fits"))


def test_failed_overwrite_retains_previous_output_and_cleans_temporary(tmp_path):
    from skysub.sky_subtraction.product import publish
    output = tmp_path / "result.fits"
    output.write_bytes(b"previous result")

    class BrokenProduct:
        def writeto(self, path, **kwargs):
            Path(path).write_bytes(b"partial replacement")
            raise OSError("Simulated write failure")

    with pytest.raises(OSError, match="Simulated"):
        publish(BrokenProduct(), output, 0o644, overwrite=True)
    assert output.read_bytes() == b"previous result"
    assert not list(tmp_path.glob(".skyml-*.fits"))


def test_default_publication_refuses_a_destination_created_during_write(tmp_path):
    from skysub.sky_subtraction.product import publish
    output = tmp_path / "result.fits"

    class RacingProduct:
        def writeto(self, path, **kwargs):
            Path(path).write_bytes(b"new result")
            output.write_bytes(b"other writer")

    with pytest.raises(FileExistsError):
        publish(RacingProduct(), output, 0o644)
    assert output.read_bytes() == b"other writer"
    assert not list(tmp_path.glob(".skyml-*.fits"))
