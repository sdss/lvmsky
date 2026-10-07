"""Sky-line scaling: recovers band brightness under a science continuum."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mlp_predictor.sky_line_scaling import sky_line_scaling_correction  # noqa: E402

WAVE = np.arange(6000.0, 9000.0, 0.5)
SIG = 1.1                                # LSF sigma, A


def _comb(centres, amp=10.0):
    t = np.zeros_like(WAVE)
    for c in centres:
        t += amp * np.exp(-0.5 * ((WAVE - c) / SIG) ** 2)
    return t


def _templates():
    rng = np.random.default_rng(0)
    return {"OH_a": _comb(np.sort(rng.uniform(6100, 8900, 60))),
            "OH_b": _comb(np.sort(rng.uniform(6100, 8900, 60)))}


def _science_continuum():
    # Bright, curved: 30x the line peak times a slow wiggle.
    x = (WAVE - 7500.0) / 1500.0
    return 300.0 * (1.0 + 0.4 * x - 0.3 * x ** 2 + 0.05 * np.sin(WAVE / 40.0))


def _observe(T, scales, noise=0.05, extra=None, seed=1):
    rng = np.random.default_rng(seed)
    sky = 5.0 + sum(T.values())                       # predicted sky: lines + flat continuum
    obs = (5.0 + sum(s * T[k] for k, s in scales.items()) + _science_continuum()
           + rng.normal(0.0, noise, WAVE.size))
    if extra is not None:
        obs = obs + extra
    return sky, obs


def _fit(obs, sky, T, **kw):
    return sky_line_scaling_correction(WAVE, obs, sky, T, return_info=True, **kw)


def test_recovers_band_scales_under_a_bright_science_continuum():
    T = _templates()
    truth = {"OH_a": 1.06, "OH_b": 0.95}
    sky, obs = _observe(T, truth)
    corr, info = _fit(obs, sky, T, prior_sigma=1.0)
    for k, s in truth.items():
        assert info["scales"][k] == pytest.approx(s, abs=2e-3)
    # The correction is the line-flux change only; the continuum is untouched.
    expected = sum((truth[k] - 1.0) * T[k] for k in T)
    assert np.max(np.abs(corr - expected)) < 0.05


def test_robust_to_narrow_science_absorption_lines():
    T = _templates()
    truth = {"OH_a": 1.05, "OH_b": 1.0}
    rng = np.random.default_rng(3)
    # Deep stellar-like absorption lines, some landing on sky lines.
    absorb = -_comb(rng.uniform(6100, 8900, 25), amp=60.0)
    sky, obs = _observe(T, truth, extra=absorb)
    _, robust = _fit(obs, sky, T, prior_sigma=1.0)
    _, plain = _fit(obs, sky, T, prior_sigma=1.0, huber_k=1e9)
    err_r = max(abs(robust["scales"][k] - truth[k]) for k in T)
    err_p = max(abs(plain["scales"][k] - truth[k]) for k in T)
    assert err_r < 0.01 and err_r < 0.5 * err_p


def test_a_template_without_signal_stays_at_the_prediction():
    T = _templates()
    T["weak"] = 1e-6 * _comb([7000.0])
    sky, obs = _observe({k: v for k, v in T.items() if k != "weak"}, {"OH_a": 1.0, "OH_b": 1.0}, noise=1.0)
    sky = sky + T["weak"]
    _, info = _fit(obs, sky, T, prior_sigma=0.1)
    assert info["scales"]["weak"] == pytest.approx(1.0, abs=1e-3)


def test_masked_pixels_are_ignored():
    T = _templates()
    spike = np.zeros_like(WAVE)
    sel = np.abs(WAVE - 6563.0) < 4.0
    spike[sel] = 1e4                                 # a science emission line
    sky, obs = _observe(T, {"OH_a": 1.0, "OH_b": 1.0}, extra=spike)
    _, info = _fit(obs, sky, T, mask=sel, prior_sigma=1.0, huber_k=1e9)
    assert all(abs(s - 1.0) < 2e-3 for s in info["scales"].values())


def test_no_templates_means_no_correction():
    sky = np.ones_like(WAVE)
    corr = sky_line_scaling_correction(WAVE, sky, sky, {"none": np.zeros_like(WAVE)})
    assert np.all(corr == 0.0)
    with pytest.raises(ValueError):
        sky_line_scaling_correction(WAVE[:-1], sky, sky, _templates())
