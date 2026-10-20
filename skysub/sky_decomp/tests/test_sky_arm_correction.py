"""Sky-arm residual correction: recovers a shared pattern, rejects broken arms."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mlp_predictor.sky_arm_correction import (  # noqa: E402
    SkyArm,
    line_protection_masks,
    protect_science_lines,
    protection_velocities,
    sky_arm_residual_correction,
)

WAVE = np.arange(3600.0, 6000.0, 0.5)
AIRGLOW = 1.0


def _pattern():
    # A solar-line-like fractional residual: narrow dips, zero mean-ish.
    p = np.zeros_like(WAVE)
    for c in (3934.0, 3968.0, 4101.0, 4300.0, 4861.0):
        p -= 0.3 * np.exp(-0.5 * ((WAVE - c) / 1.2) ** 2)
    return p


def _arm(solar, pattern, noise_sigma, rng):
    # The residual is a fraction of the arm's TOTAL model, so the model ratio
    # carries it to the science exactly.
    model = solar + AIRGLOW
    obs = model + pattern * model + rng.normal(0.0, noise_sigma, WAVE.size)
    return SkyArm(obs, model, solar, np.full(WAVE.size, noise_sigma ** 2))


def _sci(solar):
    return solar + AIRGLOW, solar      # (sci_model, sci_solar)


def test_recovers_the_shared_pattern_at_the_science_model_level():
    rng = np.random.default_rng(0)
    p = _pattern()
    s_sci = np.full(WAVE.size, 10.0)
    near = _arm(np.full(WAVE.size, 8.0), p, 1e-3, rng)
    far = _arm(np.full(WAVE.size, 12.0), p, 1e-3, rng)
    corr = sky_arm_residual_correction(WAVE, *_sci(s_sci), [near, far], smoothing=0)
    core = WAVE < 5000                  # inside the requested range
    # The truth the science prediction misses is p * M_sci.
    assert np.allclose(corr[core], (p * (s_sci + AIRGLOW))[core], atol=5e-3)
    assert np.all(corr[WAVE >= 5100] == 0.0)


def test_the_band_edge_is_a_smooth_taper_not_a_step():
    s = np.full(WAVE.size, 10.0)
    arm = SkyArm(np.full(WAVE.size, 12.0), np.full(WAVE.size, 11.0), s, np.full(WAVE.size, 1e-6))
    corr = sky_arm_residual_correction(WAVE, *_sci(s), [arm], smoothing=0, scale=None,
                                       max_wavelength=5000.0, taper_width=100.0)
    # Full correction up to the boundary, a ramp just outside it, zero beyond.
    assert np.allclose(corr[WAVE <= 5000], 1.0)
    assert np.all(corr[WAVE >= 5100] == 0.0)
    ramp = corr[(WAVE >= 5000) & (WAVE <= 5100)]
    assert np.all(np.diff(ramp) <= 0) and ramp[0] > 0.99 and ramp[-1] < 0.01
    assert np.max(np.abs(np.diff(corr))) < 0.02          # no step anywhere
    # The full band leaves it untapered.
    full = sky_arm_residual_correction(WAVE, *_sci(s), [arm], smoothing=0, scale=None,
                                       max_wavelength=np.inf)
    assert np.allclose(full, 1.0)


def test_a_broken_arm_is_gated_out():
    rng = np.random.default_rng(1)
    p = _pattern()
    s = np.full(WAVE.size, 10.0)
    near = _arm(s, p, 1e-3, rng)
    bad = SkyArm(near.observed + 50.0, near.model, s, near.variance)   # residual ~5x continuum
    corr, info = sky_arm_residual_correction(WAVE, *_sci(s), [near, bad], smoothing=0,
                                             return_info=True)
    assert info["arm_used"].tolist() == [True, False]
    corr_near = sky_arm_residual_correction(WAVE, *_sci(s), [near], smoothing=0)
    assert np.allclose(corr, corr_near)


def test_inverse_variance_weighting_favours_the_quiet_arm():
    rng = np.random.default_rng(2)
    s = np.full(WAVE.size, 10.0)
    quiet = _arm(s, _pattern(), 1e-3, rng)
    noisy = _arm(s, _pattern(), 1e-1, rng)
    corr = sky_arm_residual_correction(WAVE, *_sci(s), [quiet, noisy], smoothing=0)
    truth = _pattern() * (s + AIRGLOW)
    core = WAVE < 4900
    err_both = np.std((corr - truth)[core])
    err_noisy = np.std((sky_arm_residual_correction(WAVE, *_sci(s), [noisy], smoothing=0)
                        - truth)[core])
    assert err_both < 0.05 * err_noisy


def test_auto_smoothing_only_touches_faint_rows():
    rng = np.random.default_rng(3)
    p = _pattern()
    bright, faint = np.full(WAVE.size, 100.0), np.full(WAVE.size, 0.01)
    a, b = _arm(bright, p, 1e-2, rng), _arm(faint, p, 1e-2, rng)
    arm = SkyArm(*(np.vstack([getattr(a, f), getattr(b, f)])
                   for f in ("observed", "model", "solar", "variance")))
    sols = np.vstack([bright, faint])
    _, info = sky_arm_residual_correction(WAVE, sols + AIRGLOW, sols, [arm],
                                          auto_snr=100.0, return_info=True)
    assert info["smoothed"].tolist() == [False, True]


def test_single_row_matches_the_batch_row_and_validates_inputs():
    rng = np.random.default_rng(4)
    s = np.full(WAVE.size, 10.0)
    a1, a2 = _arm(s, _pattern(), 1e-2, rng), _arm(s, _pattern(), 1e-2, rng)
    stacked = SkyArm(*(np.vstack([getattr(a1, f), getattr(a2, f)])
                       for f in ("observed", "model", "solar", "variance")))
    m, sol = _sci(s)
    batch = sky_arm_residual_correction(WAVE, np.vstack([m, m]), np.vstack([sol, sol]), [stacked])
    assert np.allclose(batch[1], sky_arm_residual_correction(WAVE, m, sol, [a2]))
    with pytest.raises(ValueError):
        sky_arm_residual_correction(WAVE, m, sol, [a1], smoothing="median")
    with pytest.raises(ValueError):
        sky_arm_residual_correction(WAVE, m, sol, [a1], scale="oh")
    with pytest.raises(ValueError):
        sky_arm_residual_correction(WAVE[:-1], m, sol, [a1])


def test_a_narrow_feature_in_one_arm_is_not_transferred():
    rng = np.random.default_rng(3)
    p = _pattern()
    s = np.full(WAVE.size, 10.0)
    near = _arm(s, p, 1e-3, rng)
    far = _arm(s, p, 1e-3, rng)
    line = 5.0 * np.exp(-0.5 * ((WAVE - 4500.0) / 1.06) ** 2)     # emission in the near field only
    near = SkyArm(near.observed + line, near.model, near.solar, near.variance)
    win = np.abs(WAVE - 4500.0) < 6
    kw = dict(smoothing=0, max_wavelength=np.inf)
    off = sky_arm_residual_correction(WAVE, *_sci(s), [near, far], consistency_clip=None, **kw)
    on, info = sky_arm_residual_correction(WAVE, *_sci(s), [near, far], return_info=True, **kw)
    assert off[win].sum() > 0.4 * line[win].sum()                 # half the line, without the clip
    assert abs(on[win].sum()) < 0.02 * line[win].sum()            # gone with it (the default)
    # The shared pattern is still corrected away from the line.
    away = (np.abs(WAVE - 4500.0) > 10) & (WAVE < 5900)
    assert np.allclose(on[away], (p * (s + AIRGLOW))[away], atol=5e-3)
    assert 0 < info["clipped_fraction"] < 0.05


# --- Science-line protection over the full band ---------------------------

FULL = np.arange(3600.0, 9800.0, 0.5)
C_KM_S = 299792.458


def _gauss(centre, amplitude, sigma=1.0):
    return amplitude * np.exp(-0.5 * ((FULL - centre) / sigma) ** 2)


def test_line_masks_cover_the_lines_at_rest_and_at_the_cloud_velocity():
    windows, cores = line_protection_masks(FULL, 1.5, velocities=(262.0,))
    at = lambda lam: windows[np.argmin(np.abs(FULL - lam))]
    assert at(6562.8) and at(9068.6) and at(4861.3)
    assert at(6562.8 * (1 + 262.0 / C_KM_S))
    assert not at(6450.0)
    assert np.all(windows[cores])                       # every core lies inside its window
    # Interpolated blue windows are narrower (300 km/s) than the red ones (500 km/s).
    rest, _ = line_protection_masks(FULL, 1.5)
    run = lambda lam: rest[np.abs(FULL - lam) < 15].sum() * 0.5
    assert abs(run(4861.33) - 2 * 4861.33 * 300 / C_KM_S) < 1.0
    assert abs(run(6678.15) - 2 * 6678.15 * 500 / C_KM_S) < 1.0
    assert 0.1 < windows.mean() < 0.3
    assert protection_velocities(81.0, -69.7) == (262.0,)
    assert protection_velocities(14.0, -73.0) == (146.0,)
    assert protection_velocities(180.0, 10.0) == ()


def test_emission_in_both_arms_is_not_transferred():
    # Diffuse emission common to both fields passes the consistency clip; the
    # line protection must stop it, also where it sits on an OH line.
    continuum = np.full(FULL.size, 10.0)
    oh = _gauss(9068.6, 30.0, 1.0) + _gauss(6562.0, 5.0, 1.0)
    model = continuum + oh
    lines = sum(_gauss(lam, 10.0, 1.0) for lam in (4861.33, 6562.80, 6583.45, 9068.60))
    arm = lambda: SkyArm(model + lines, model, continuum, np.full(FULL.size, 1e-4))
    kw = dict(smoothing=0, max_wavelength=np.inf, scale=None)
    off = sky_arm_residual_correction(FULL, model, continuum, [arm(), arm()], **kw)
    on = sky_arm_residual_correction(FULL, model, continuum, [arm(), arm()], sci_oh=oh,
                                     line_masks=line_protection_masks(FULL, 1.5), **kw)
    for lam in (4861.33, 6562.80, 6583.45, 9068.60):
        w = np.abs(FULL - lam) < 4
        assert off[w].sum() > 0.9 * lines[w].sum()
        assert abs(on[w].sum()) < 1e-6 * lines[w].sum()


def test_blue_windows_are_interpolated_and_red_windows_keep_the_oh_part():
    windows, cores = line_protection_masks(FULL, 1.5)
    corr = 0.01 * (FULL - 3600.0)                        # a correction linear in wavelength
    oh_frac = np.where(np.abs(FULL - 6575.0) < 1.0, 0.8, 0.0)   # an OH line inside the Halpha window
    out = protect_science_lines(corr, FULL, windows, cores, oh_frac)
    blue = FULL < 5800
    assert np.allclose(out[blue], corr[blue], atol=1e-9)   # a linear correction survives interpolation
    red = (FULL >= 5800) & windows & ~cores
    assert np.allclose(out[red], (corr * oh_frac)[red])
    assert np.all(out[cores & (FULL >= 5800)] == 0.0)
    assert np.array_equal(out[~windows], corr[~windows])
    zeroed = protect_science_lines(corr, FULL, windows, cores, oh_frac, mode="windows")
    assert np.all(zeroed[windows & (FULL >= 5800)] == 0.0)
    with pytest.raises(ValueError):
        sky_arm_residual_correction(FULL, np.ones(FULL.size), np.ones(FULL.size),
                                    [SkyArm(*(np.ones(FULL.size),) * 4)], line_masks=(windows, cores))
