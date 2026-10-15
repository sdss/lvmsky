"""Observed evidence, not model-pixel counts, controls the new LSF gate."""

import numpy as np

from skysub.sky_subtraction.lsf_information import line_information


def spectrum(continuum=1.0, noise=0.05):
    wave = np.arange(6000.0, 6800.0, 0.5)
    lines = sum(5 * np.exp(-0.5 * ((wave - c) / 0.7) ** 2) for c in np.arange(6040.0, 6800.0, 60.0))
    source = np.full_like(wave, continuum)
    return wave, source + lines, np.full_like(wave, noise**-2), lines, source, np.ones(len(wave), bool)


def test_lines_across_channel_with_strong_contrast_and_snr_pass():
    evidence = line_information(*spectrum())
    assert evidence['line_information_pass']
    assert np.all(evidence['quality_bin_counts'] >= 2)


def test_bright_continuum_and_low_snr_each_keep_survey():
    assert not line_information(*spectrum(continuum=20.0))['line_information_pass']
    assert not line_information(*spectrum(noise=10.0))['line_information_pass']


def test_many_lines_in_one_region_cannot_support_whole_channel():
    args = list(spectrum())
    wave = args[0]
    args[-1] = wave < 6400
    info = line_information(*args)
    assert not info['line_information_pass']
    assert np.all(info['quality_bin_counts'][2:] == 0)


def test_stellar_absorption_can_cancel_an_expected_sky_line():
    args = list(spectrum())
    args[1] = args[4].copy()
    info = line_information(*args)
    assert not info['line_information_pass']
    assert np.all(info['quality_bin_counts'] == 0)
