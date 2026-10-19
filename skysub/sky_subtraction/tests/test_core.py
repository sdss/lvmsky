"""Focused contracts for selection, masked medians, and continuous LSF recovery."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
from skysub.sky_decomp.lsf_spline2d import _integrated_components, _line_design, evaluate_lsf_diagnostics
from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
from skysub.sky_subtraction import fiber_lsf
from skysub.sky_subtraction.fiber_lsf import fit_channel
from skysub.sky_subtraction.io import median_spectrum, select_fibers

ROOT = Path(__file__).parents[2] / 'sky_decomp/data/lsf'


def lines_only(lines):
    return dict(lines=lines, solar=np.zeros_like(lines), fixed=np.zeros(len(lines)))


def test_selection_join_order_and_science_fallback():
    n = 70
    ids = np.arange(n, 0, -1)
    slit = Table(dict(fiberid=ids, fibstatus=np.zeros(n, int), telescope=['Sci'] * 60 + ['SkyE'] * 5 + ['SkyW'] * 5))
    table = Table(
        dict(
            n_gaia_sources=np.ones(n, int),
            gaia_g_n_valid_sources=np.ones(n, int),
            gaia_g_flux=np.ones(n),
            ratio_gaia_to_lvm_flux_plus_sky_g=np.r_[np.linspace(0.005, 0.1, 60), np.full(10, 0.001)],
            lvm_flux_plus_sky_g=np.ones(n),
            lvm_flux_plus_sky_g_valid_fraction=np.ones(n),
        )
    )
    selection = select_fibers(slit, table)
    selected = selection.indices
    assert len(selected['Sci']) == 50
    assert 'selected' not in table.colnames
    annotated = selection.annotate(table)
    assert annotated['selected'].sum() == 60 and set(annotated['selection_reason']) == {
        'lowest_ratio_50',
        'gaia_clean',
        'not_selected',
    }
    assert len(selected['SkyE']) == len(selected['SkyW']) == 5
    assert np.max(table['ratio_gaia_to_lvm_flux_plus_sky_g'][selected['Sci']]) < 0.085
    table['gaia_g_n_valid_sources'][60:65] = 0
    with pytest.raises(ValueError, match='SkyE'):
        select_fibers(slit, table)


def test_masked_median_and_gaussian_variance():
    rng = np.random.default_rng(44)
    samples = rng.normal(size=(101, 10000))
    flux, ivar, count = median_spectrum(samples, np.ones_like(samples), np.zeros_like(samples), np.arange(101))
    assert np.var(flux) == pytest.approx(np.pi / (2 * 101), rel=0.05)
    assert np.mean(1 / ivar) == pytest.approx(np.pi / (2 * 101))
    mask = np.zeros((3, 3), int)
    mask[:2, 0] = 1
    mask[:, 2] = 1
    f, v, n = median_spectrum(np.array([[100, 1, 1], [200, 2, 2], [3, 3, 3.0]]), np.ones((3, 3)), mask, np.arange(3))
    assert f[0] == 3 and v[0] == 1 and np.isnan(f[2]) and v[2] == 0
    assert list(n) == [1, 3, 0]


@pytest.mark.parametrize('channel,lower,upper', [('R', 5787.0, 7454.0), ('Z', 7454.0, 9800.0)])
@pytest.mark.parametrize('contrast', [0.0, 1.0, 10.0])
def test_known_lsf_under_stellar_structure(channel, lower, upper, contrast):
    rng = np.random.default_rng(12)
    ref = load_survey_ascii_lsf(ROOT, 500)
    truth = deepcopy(ref)
    c = ref.coefficients[channel]
    broadened = 0.5 * c + 0.25 * np.roll(c, 1, axis=0) + 0.25 * np.roll(c, -1, axis=0)
    broadened /= broadened.sum(axis=0)
    truth.coefficients[channel] = broadened
    wave = np.arange(lower, upper, 0.5)
    lines = np.arange(lower + 12, upper - 12, 22.7)
    knots, degree = ref.knot_vectors[channel], ref.degrees[channel]
    b = evaluate_bspline_basis(np.clip(lines, knots[degree], knots[-degree - 1]), knots, degree)
    design = _line_design(
        _integrated_components(wave, lines, np.ones(len(wave), bool), 11, 2.5), np.full(len(lines), 5.0), b, 11
    ).toarray()
    sky = design @ truth.coefficients[channel].ravel()
    # Broad molecular bands plus narrow stellar absorption crossing some sky lines.
    star = contrast * (1.5 + 0.4 * np.sin(wave / 90) - 0.6 * np.exp(-0.5 * ((wave - (lower + upper) / 2) / 35) ** 2))
    star -= contrast * 0.15 * sum(np.exp(-0.5 * ((wave - center) / 0.8) ** 2) for center in lines[::7])
    flux = 1.1 * sky + star + rng.normal(0, 0.01, len(wave))
    fit = fit_channel(
        wave, flux, np.full(len(wave), 1e4), lines_only(design), ref, channel, np.zeros(len(wave), bool), ref
    )
    recovered, info, candidate = fit.adopted, fit.qa, fit.candidate
    if contrast == 0:
        assert info['accepted'], info
    if info['accepted']:
        probes = np.linspace(lower + 1, upper - 1, 9)
        d = evaluate_lsf_diagnostics(candidate, probes)
        t = evaluate_lsf_diagnostics(truth, probes)
        assert np.max(np.abs(d['w50_angstrom'] / t['w50_angstrom'] - 1)) < 0.12
    else:
        np.testing.assert_array_equal(recovered, c)
    assert info['reason'] in ('insitu', 'profile_guard', 'line_information')
    if contrast == 10:
        assert not info['accepted'] and not info['line_information_pass']


def test_insufficient_information_keeps_survey():
    ref = load_survey_ascii_lsf(ROOT, 500)
    wave = np.arange(6000.0, 6400.0, 0.5)
    design = np.zeros((len(wave), ref.coefficients['R'].size))
    fit = fit_channel(
        wave, np.ones(len(wave)), np.ones(len(wave)), lines_only(design), ref, 'R', np.zeros(len(wave), bool), ref
    )
    recovered, info = fit.adopted, fit.qa
    assert not info['accepted'] and info['reason'] == 'insufficient_lines'
    np.testing.assert_array_equal(recovered, ref.coefficients['R'])


def test_cframe_gaia_and_cached_photometry_agree(tmp_path, monkeypatch):
    from astropy.io import fits

    from skysub.medians_computation.src.lvm_medians import gaia
    from skysub.sky_subtraction.io import CFrame, resolve_gaia

    wave = np.arange(3600.0, 9800.0, 2.0)
    flux = np.full((3, len(wave)), 1e-14)
    ivar = np.full_like(flux, 1e30)
    mask = np.zeros_like(flux, int)
    slit = Table(
        dict(
            fiberid=[3, 1, 2],
            ra=[10.0, 11.0, 12.0],
            dec=[-20.0, -21.0, -22.0],
            fibstatus=[0, 0, 0],
            telescope=['Sci', 'SkyE', 'SkyW'],
        )
    )
    path = tmp_path / 'lvmCFrame-00000007.fits'
    h = fits.PrimaryHDU()
    h.header['EXPOSURE'] = 7
    f = fits.ImageHDU(flux, name='FLUX')
    f.header['BUNIT'] = 'erg / (Angstrom s cm2)'
    fits.HDUList(
        [
            h,
            f,
            fits.ImageHDU(ivar, name='IVAR'),
            fits.ImageHDU(mask, name='MASK'),
            fits.ImageHDU(wave, name='WAVE'),
            fits.BinTableHDU(slit, name='SLITMAP'),
        ]
    ).writeto(path)
    sources = Table(dict(fiberid=[3], phot_g_mean_mag=[18.0], phot_g_mean_flux_over_error=[100.0]))
    from_tap, _ = gaia.derive_fiber_table(path, sources, observed_total=True)
    cache = tmp_path / 'cache'
    cache.mkdir()
    prim = fits.PrimaryHDU()
    prim.header['EXPNUM'] = 7
    fits.HDUList([prim, fits.BinTableHDU(from_tap[[2, 0, 1]], name='FIBERS')]).writeto(
        cache / 'lvmGAIA-fibers-00000007.fits'
    )
    frame = CFrame(h.header, slit, wave, flux, ivar, mask, None)
    cached, _, info = resolve_gaia(path, frame, cache, tmp_path / 'tap')
    monkeypatch.setattr(gaia, '_download', lambda *args, **kwargs: sources)
    queried, _, tap = resolve_gaia(path, frame, None, tmp_path / 'tap')
    assert list(cached['fiberid']) == list(slit['fiberid'])
    np.testing.assert_allclose(
        cached['ratio_gaia_to_lvm_flux_plus_sky_g'], queried['ratio_gaia_to_lvm_flux_plus_sky_g']
    )
    np.testing.assert_allclose(queried['lvm_flux_plus_sky_g'], 1e-14)
    assert list(queried['n_gaia_sources']) == [1, 0, 0]
    assert info['mode'] == 'aggregate cache'
    assert tap['mode'] == 'ARI TAP' and tap['cache_miss'] == 'no cache directory'


def test_observed_total_preserves_sframe_semantics(tmp_path):
    from astropy.io import fits

    from skysub.medians_computation.src.lvm_medians import gaia

    wave = np.arange(3600.0, 9800.0, 2.0)
    flux = np.full((1, len(wave)), 1e-14)
    slit = Table(dict(fiberid=[1], ra=[10.0], dec=[20.0]))
    source = Table(dict(fiberid=[], phot_g_mean_mag=[], phot_g_mean_flux_over_error=[]))
    h = fits.PrimaryHDU()
    h.header['EXPOSURE'] = 7
    f = fits.ImageHDU(flux, name='FLUX')
    f.header['BUNIT'] = 'erg / (Angstrom s cm2)'
    path = tmp_path / 'lvmSFrame-00000007.fits'
    fits.HDUList(
        [
            h,
            f,
            fits.ImageHDU(np.full_like(flux, 1e30), name='IVAR'),
            fits.ImageHDU(wave, name='WAVE'),
            fits.ImageHDU(2 * flux, name='SKY'),
            fits.ImageHDU(np.full_like(flux, 1e30), name='SKY_IVAR'),
            fits.BinTableHDU(slit, name='SLITMAP'),
        ]
    ).writeto(path)
    original, _ = gaia.derive_fiber_table(path, source)
    total, _ = gaia.derive_fiber_table(path, source, observed_total=True)
    np.testing.assert_allclose(original['lvm_flux_plus_sky_g'], 3e-14)
    np.testing.assert_allclose(total['lvm_flux_plus_sky_g'], 1e-14)


@pytest.mark.parametrize('channel,lower,upper', [('R', 5787.0, 7454.0), ('Z', 7454.0, 9800.0)])
@pytest.mark.parametrize('contrast', [0.3, 3.0, 10.0])
def test_empirical_m5v_bands(channel, lower, upper, contrast):
    ref = load_survey_ascii_lsf(ROOT, 500)
    wave = np.arange(lower, upper, 0.5)
    lines = np.arange(lower + 12, upper - 12, 22.7)
    knots, degree = ref.knot_vectors[channel], ref.degrees[channel]
    b = evaluate_bspline_basis(np.clip(lines, knots[degree], knots[-degree - 1]), knots, degree)
    design = _line_design(
        _integrated_components(wave, lines, np.ones(len(wave), bool), 11, 2.5), np.full(len(lines), 5.0), b, 11
    ).toarray()
    stellar = np.loadtxt(Path(__file__).parent / 'data/pickles_m5v.csv', delimiter=',')
    star = contrast * np.interp(wave, stellar[:, 0], stellar[:, 1])
    sky = design @ ref.coefficients[channel].ravel()
    flux = sky + star + np.random.default_rng(32).normal(0, 0.01, len(wave))
    fit = fit_channel(
        wave, flux, np.full(len(wave), 1e4), lines_only(design), ref, channel, np.zeros(len(wave), bool), ref
    )
    recovered, info, candidate = fit.adopted, fit.qa, fit.candidate
    if info['accepted']:
        probes = np.linspace(lower + 1, upper - 1, 9)
        d = evaluate_lsf_diagnostics(candidate, probes)
        t = evaluate_lsf_diagnostics(ref, probes)
        assert np.max(np.abs(d['w50_angstrom'] / t['w50_angstrom'] - 1)) < 0.12
    else:
        np.testing.assert_array_equal(recovered, ref.coefficients[channel])


def test_equal_coefficient_distance_has_different_profile_effects():
    ref = load_survey_ascii_lsf(ROOT, 500)
    a, b = deepcopy(ref), deepcopy(ref)
    # Equal movements in coefficient space redistribute very different displacement masses.
    eps = 0.01
    a.coefficients['R'][5] -= eps
    a.coefficients['R'][4] += eps
    b.coefficients['R'][5] -= eps
    b.coefficients['R'][0] += eps
    assert np.linalg.norm(a.coefficients['R'] - ref.coefficients['R']) == pytest.approx(
        np.linalg.norm(b.coefficients['R'] - ref.coefficients['R'])
    )
    probes = np.array([6500.0])
    assert not np.allclose(
        evaluate_lsf_diagnostics(a, probes)['centroid_angstrom'],
        evaluate_lsf_diagnostics(b, probes)['centroid_angstrom'],
    )


def test_narrow_science_emission_is_preserved():
    ref = load_survey_ascii_lsf(ROOT, 500)
    channel = 'R'
    wave = np.arange(5787.0, 7454.0, 0.5)
    lines = np.arange(5800.0, 7440.0, 22.7)
    knots, degree = ref.knot_vectors[channel], ref.degrees[channel]
    basis = evaluate_bspline_basis(np.clip(lines, knots[degree], knots[-degree - 1]), knots, degree)
    design = _line_design(
        _integrated_components(wave, lines, np.ones(len(wave), bool), 11, 2.5), np.full(len(lines), 5.0), basis, 11
    ).toarray()
    sky = design @ ref.coefficients[channel].ravel()
    centre = lines[len(lines) // 2]
    object_line = 30 * np.exp(-0.5 * ((wave - centre) / 1.0) ** 2)
    science_mask = np.abs(wave - centre) < 5
    flux = sky + object_line + 1.0
    fit = fit_channel(wave, flux, np.full(len(wave), 1e4), lines_only(design), ref, channel, science_mask, ref)
    recovered, info = fit.adopted, fit.qa
    assert info['accepted']
    residual = flux - design @ recovered.ravel() - 1.0
    assert np.sum(residual[science_mask]) / np.sum(object_line[science_mask]) == pytest.approx(1.0, abs=0.01)


def test_mean_profile_approximates_median_for_small_fiber_variations():
    from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_density

    ref = load_survey_ascii_lsf(ROOT, 500)
    delta = np.linspace(-2.5, 2.5, 501)
    profiles = []
    coefficients = []
    for fraction in np.linspace(0, 0.12, 21):
        state = deepcopy(ref)
        c = ref.coefficients['B']
        wide = 0.5 * c + 0.25 * np.roll(c, 1, axis=0) + 0.25 * np.roll(c, -1, axis=0)
        state.coefficients['B'] = (1 - fraction) * c + fraction * wide
        profiles.append(evaluate_lsf_density(state, np.array([4500.0]), delta)[0])
        coefficients.append(state.coefficients['B'])
    mean = deepcopy(ref)
    mean.coefficients['B'] = np.mean(coefficients, axis=0)
    profile = evaluate_lsf_density(mean, np.array([4500.0]), delta)[0]
    assert np.max(np.abs(profile - np.median(profiles, axis=0))) / profile.max() < 0.01


@pytest.mark.parametrize('channel,lower,upper', [('R', 5787.0, 7454.0), ('Z', 7454.0, 9800.0)])
def test_full_sky_with_strong_solar_absorption(channel, lower, upper):
    """A bright structured sky continuum must enter the LSF model, not the source spline."""
    ref = load_survey_ascii_lsf(ROOT, 500)
    truth = deepcopy(ref)
    c = ref.coefficients[channel]
    wide = 0.5 * c + 0.25 * np.roll(c, 1, axis=0) + 0.25 * np.roll(c, -1, axis=0)
    truth.coefficients[channel] = 0.7 * c + 0.3 * wide
    wave = np.arange(lower, upper, 0.5)
    knots, degree = ref.knot_vectors[channel], ref.degrees[channel]
    basis = evaluate_bspline_basis(np.clip(wave, knots[degree], knots[-degree - 1]), knots, degree)
    absorption = sum(
        0.4 * np.exp(-0.5 * ((wave - centre) / 0.65) ** 2) for centre in np.arange(lower + 8, upper - 8, 17.3)
    )
    solar_source = 50 * (1 - absorption)
    solar = _line_design(
        _integrated_components(wave, wave, np.ones(len(wave), bool), 11, 2.5), solar_source * 0.5, basis, 11
    ).toarray()
    centres = np.arange(lower + 12, upper - 12, 22.7)
    line_basis = evaluate_bspline_basis(np.clip(centres, knots[degree], knots[-degree - 1]), knots, degree)
    lines = _line_design(
        _integrated_components(wave, centres, np.ones(len(wave), bool), 11, 2.5),
        np.full(len(centres), 5.0),
        line_basis,
        11,
    ).toarray()
    stellar = 3 + np.sin(wave / 90)
    sky = (solar + lines) @ truth.coefficients[channel].ravel()
    flux = sky + stellar + np.random.default_rng(41).normal(0, 0.01, len(wave))
    fit = fit_channel(
        wave,
        flux,
        np.full(len(wave), 1e4),
        dict(lines=lines, solar=solar, fixed=np.zeros(len(wave))),
        ref,
        channel,
        np.zeros(len(wave), bool),
        ref,
    )
    recovered, info, candidate = fit.adopted, fit.qa, fit.candidate
    assert not info['accepted'] and info['reason'] == 'line_information'
    np.testing.assert_array_equal(recovered, ref.coefficients[channel])
    probes = np.linspace(lower + 1, upper - 1, 9)
    expected = evaluate_lsf_diagnostics(truth, probes)
    actual = evaluate_lsf_diagnostics(candidate, probes)
    assert np.max(np.abs(actual['w50_angstrom'] / expected['w50_angstrom'] - 1)) < 0.025
    assert np.all(fit.history[-1]['continuum_coef'] >= 0)


def _synthetic_channel(channel='R', lower=5787.0, upper=7454.0):
    ref = load_survey_ascii_lsf(ROOT, 500)
    wave = np.arange(lower, upper, 0.5)
    lines = np.arange(lower + 12, upper - 12, 22.7)
    knots, degree = ref.knot_vectors[channel], ref.degrees[channel]
    basis = evaluate_bspline_basis(np.clip(lines, knots[degree], knots[-degree - 1]), knots, degree)
    design = _line_design(
        _integrated_components(wave, lines, np.ones(len(wave), bool), 11, 2.5), np.full(len(lines), 5.0), basis, 11
    ).toarray()
    flux = design @ ref.coefficients[channel].ravel() + 1.0
    return ref, wave, flux, lines_only(design)


def test_lsf_solver_failure_keeps_reason_and_survey(monkeypatch):
    ref, wave, flux, design = _synthetic_channel()
    original, calls = fiber_lsf.fit_lsf_surface, []

    def failing_second_call(*args):
        calls.append(1)
        coefficient, metric = original(*args)
        if len(calls) == 2:
            metric.update(status='Failed', reason='solver_failed')
            coefficient = np.full_like(coefficient, np.nan)
        return coefficient, metric

    monkeypatch.setattr(fiber_lsf, 'fit_lsf_surface', failing_second_call)
    fit = fit_channel(wave, flux, np.full(len(wave), 1e4), design, ref, 'R', np.zeros(len(wave), bool), ref)
    assert fit.qa['reason'] == 'lsf_solver_solver_failed' and not fit.qa['accepted']
    assert len(fit.history) == 1 and np.isfinite(fit.candidate.coefficients['R']).all()
    np.testing.assert_array_equal(fit.adopted, ref.coefficients['R'])


def test_invalid_profile_is_a_rejection_not_an_error(monkeypatch):
    ref, wave, flux, design = _synthetic_channel()

    def invalid(_):
        raise ValueError('not unit-sum')

    monkeypatch.setattr(fiber_lsf, 'validate_coefficients', invalid)
    fit = fit_channel(wave, flux, np.full(len(wave), 1e4), design, ref, 'R', np.zeros(len(wave), bool), ref)
    assert fit.qa['reason'] == 'invalid_profile' and not fit.qa['accepted']
    np.testing.assert_array_equal(fit.adopted, ref.coefficients['R'])


def test_channel_slices_follow_open_state_bounds():
    from skysub.sky_subtraction.rendering import channel_slices

    ref = load_survey_ascii_lsf(ROOT, 500)
    wave = np.arange(3600.0, 9900.0, 0.5)
    takes = channel_slices(wave, ref)
    assert sum(t.astype(int) for t in takes.values()).tolist() == [1] * len(wave)
    assert takes['Z'][-1] and takes['B'][0] and not takes['R'][wave >= 7454.0].any()


def test_settings_record_roundtrip_ignores_unknown_keys():
    import json

    from skysub.sky_subtraction.settings import SETTINGS, Settings

    record = json.loads(json.dumps(dict(SETTINGS.asdict(), jobs=4, gaia_threshold=0.01)))
    assert Settings.from_record(record) == SETTINGS
    assert Settings.from_record(dict(variance_floor=0.05)).variance_floor == 0.05
    assert Settings.from_record(None) == SETTINGS


def test_product_layout_follows_lvmsframe(tmp_path):
    from astropy.io import fits

    from skysub.sky_subtraction.product import SHARED_ORDER, lvmsframe_layout

    source = tmp_path / 'lvmCFrame-1.fits'
    cols = fits.ColDefs([fits.Column('x', 'D', unit='mm', array=np.zeros(2))])
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols, name='SLITMAP')]).writeto(source)
    names = ['SKYML_META', 'BLOCK_ONLY', *reversed(SHARED_ORDER)]
    hdus = [fits.PrimaryHDU(), *(fits.ImageHDU(name=n) for n in names if n != 'SLITMAP')]
    hdus.append(fits.BinTableHDU(Table({'x': np.zeros(2)}), name='SLITMAP'))
    product = lvmsframe_layout(fits.HDUList(hdus), tmp_path / 'lvmSFrame-skyml-1.fits', source)
    assert [x.name for x in product[1:]] == [*SHARED_ORDER, 'BLOCK_ONLY', 'SKYML_META']
    assert product[0].header['FILENAME'] == 'lvmSFrame-skyml-1.fits'
    assert product['SLITMAP'].columns['x'].unit == 'mm'


def test_unfitted_qa_has_the_schema_of_a_fitted_channel():
    from skysub.sky_subtraction.fiber_lsf import unfitted_qa

    # Column order of LSF_QA in an in-fiber product (exposure 4108), before fiberid and input_row.
    columns = [
        'channel', 'accepted', 'reason', 'n_line_pixels', 'n_supported_basis', 'iterations', 'chi2_red',
        'width_change', 'centroid_change', 'line_information_pass', 'quality_bin_counts', 'quality_bin_contrast',
        'quality_bin_snr', 'line_amplitude', 'solar_amplitude', 'downweighted', 'converged', 'median_width_change',
    ]
    assert [*unfitted_qa('R', 'sky_block')] == columns
