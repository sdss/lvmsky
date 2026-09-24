"""Contract tests for the fixed-weight A/B scorer.

The scorer is the instrument the decomposition A/B is read off, so the
properties that make its comparisons meaningful are asserted here rather
than assumed: the bands partition the grid, the weighting depends only on
the data (never on the arm being scored), a run scored against itself gives
exactly unity, and a model that is worse everywhere scores worse.
"""

import numpy as np
import pytest
from astropy.io import fits

from skysub.sky_decomp import ab_score
from skysub.sky_decomp import reliability as _rel


N_PIX, N_ROWS = 400, 12
FACTOR = 1e14


def _wave():
    return np.linspace(3700.0, 9700.0, N_PIX)


def _write_pair(tmp_path, *, model_scale=1.0, seed=0):
    """A minimal (base cube, decomp product) pair the scorer can read."""
    rng = np.random.default_rng(seed)
    wave = _wave()
    truth = 1e-14 * (2.0 + np.sin(wave / 300.0))[None, :] * np.linspace(
        0.5, 2.0, N_ROWS)[:, None]
    flux = truth + 1e-16 * rng.normal(size=truth.shape)
    # The model misses the truth by a fixed fraction; model_scale > 1 is a
    # uniformly worse arm.
    model = (truth + model_scale * (truth - flux)) * FACTOR

    base = tmp_path / f"base_{seed}.fits"
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(wave.astype(np.float32), name="WAVE"),
        fits.ImageHDU(flux.astype(np.float32), name="FLUX_SCI"),
    ]).writeto(base, overwrite=True)

    bits = np.zeros(N_ROWS, dtype=np.int64)
    bits[0] = _rel.RELIABILITY_DIFFUSE_COLLAPSED
    bits[1] = _rel.RELIABILITY_ZODI_ANCHOR_PINNED
    meta = fits.BinTableHDU.from_columns([
        fits.Column(name="reliability", format="J", array=bits),
        fits.Column(name="moon_share", format="D",
                    array=np.linspace(0.02, 0.9, N_ROWS)),
        fits.Column(name="zodi_int", format="D",
                    array=np.linspace(1e3, 5e3, N_ROWS)),
    ], name="META")
    dec = tmp_path / f"dec_{seed}_{model_scale}.fits"
    fits.HDUList([
        fits.PrimaryHDU(), meta,
        fits.ImageHDU(model, name="BESTFIT_LSF"),
    ]).writeto(dec, overwrite=True)
    return base, dec


def test_bands_partition_the_grid():
    """Every pixel falls in exactly one band, and the bands are ordered."""
    wave = np.linspace(3600.0, 9799.0, 5000)
    counts = np.zeros(wave.size, dtype=int)
    for _name, lo, hi in ab_score.BANDS:
        counts += ((wave >= lo) & (wave < hi)).astype(int)
    assert np.all(counts == 1)
    edges = [lo for _n, lo, _h in ab_score.BANDS]
    assert edges == sorted(edges)


def test_reference_ivar_ignores_the_model(tmp_path):
    """The weighting is a function of the DATA, so two arms share it exactly.

    This is the property that makes chi2 comparable across arms; if the
    weighting ever picked up anything from the fit, the score would reward
    an arm for its own weighting choice.
    """
    wave = _wave()
    sens = ab_score.absolute_sensitivity(wave)
    flux = 1e-14 * (2.0 + np.sin(wave / 300.0)) * FACTOR
    w1 = ab_score.reference_ivar(flux, wave, sens=sens, flux_scale=FACTOR)
    w2 = ab_score.reference_ivar(flux, wave, sens=sens, flux_scale=FACTOR)
    np.testing.assert_array_equal(w1, w2)
    assert np.isclose(np.mean(w1[w1 > 0]), 1.0, rtol=1e-9)


def test_reference_ivar_is_scale_free_in_factor():
    """FACTOR cancels: the same physical row scores the same at any FACTOR."""
    wave = _wave()
    sens = ab_score.absolute_sensitivity(wave)
    phys = 1e-14 * (2.0 + np.sin(wave / 300.0))
    a = ab_score.reference_ivar(phys * 1e14, wave, sens=sens, flux_scale=1e14)
    b = ab_score.reference_ivar(phys * 1e12, wave, sens=sens, flux_scale=1e12)
    np.testing.assert_allclose(a, b, rtol=1e-10)


def test_self_comparison_is_exactly_unity(tmp_path):
    base, dec = _write_pair(tmp_path)
    s = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR)
    paired = ab_score.compare_paired(s, s)
    for band in ("b", "r", "z", "full"):
        assert paired[f"ratio_{band}_median"] == pytest.approx(1.0, rel=1e-12)
        assert paired[f"n_dropped_{band}"] == 0


def test_a_uniformly_worse_arm_scores_worse(tmp_path):
    base, good = _write_pair(tmp_path, model_scale=1.0)
    _, bad = _write_pair(tmp_path, model_scale=3.0)
    s_good = ab_score.score_decomposition(base, good, kind="sci", factor=FACTOR)
    s_bad = ab_score.score_decomposition(base, bad, kind="sci", factor=FACTOR)
    # resid scales as (1 + model_scale), so chi2 scales as its square.
    paired = ab_score.compare_paired(s_good, s_bad)
    assert paired["ratio_full_median"] == pytest.approx(4.0, rel=1e-3)
    assert paired["frac_improved_full"] == 0.0
    assert s_bad["chi2_full_median"] > s_good["chi2_full_median"]


def test_reliability_rates_are_read_off_the_bits(tmp_path):
    base, dec = _write_pair(tmp_path)
    s = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR)
    assert s["frac_diffuse_collapsed"] == pytest.approx(1.0 / N_ROWS)
    assert s["frac_zodi_anchor_pinned"] == pytest.approx(1.0 / N_ROWS)
    assert s["frac_reversed"] == 0.0
    assert s["frac_any_error"] == pytest.approx(1.0 / N_ROWS)
    assert s["frac_any_warning"] == pytest.approx(1.0 / N_ROWS)


def test_mismatched_row_sets_are_refused(tmp_path):
    base, dec = _write_pair(tmp_path)
    full = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR)
    part = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR,
                                        rows=np.arange(N_ROWS - 2))
    with pytest.raises(ValueError, match="different row sets"):
        ab_score.compare_paired(full, part)


def test_masked_pixels_do_not_enter_the_score(tmp_path):
    """Pixels masked out of the fit carry no information about it."""
    base, dec = _write_pair(tmp_path)
    mask = np.zeros(N_PIX, dtype=bool)
    mask[100:140] = True
    s_all = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR)
    s_msk = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR,
                                         science_line_mask=mask)
    assert s_all["chi2_full_median"] != s_msk["chi2_full_median"]
    # ... and masking a band out entirely leaves that band unscored.
    all_b = np.zeros(N_PIX, dtype=bool)
    wave = _wave()
    all_b[(wave >= ab_score.BANDS[0][1]) & (wave < ab_score.BANDS[0][2])] = True
    s_nob = ab_score.score_decomposition(base, dec, kind="sci", factor=FACTOR,
                                         science_line_mask=all_b)
    assert np.isnan(s_nob["chi2_b_median"])
