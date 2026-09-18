"""Contract tests for the per-row reliability flags and the reversal retry.

These guard the bits that get PERSISTED in the products (so their meaning must
not drift) and the writer seam that carries them, without needing any data
asset or a real fit.
"""

import numpy as np
import pytest

from skysub.sky_decomp import reliability as rel


def _power_law(wave, slope, amplitude=1.0):
    """Power law normalised to unit total, so amplitude sets the SHARE.

    Without the normalisation a steep law like ``wave ** -3.7`` integrates to
    almost nothing next to a flat one, the moon share collapses, and the row is
    untestable by construction rather than by intent.
    """
    flux = np.asarray(wave, dtype=float) ** float(slope)
    return amplitude * flux / flux.sum()


WAVE = np.linspace(3600.0, 9800.0, 1200)
LOG_WAVE = np.log(WAVE)


def test_loglog_slope_recovers_a_known_power_law():
    for slope in (-3.7, -0.3, 0.0, 1.4):
        got = rel.loglog_slope(_power_law(WAVE, slope), LOG_WAVE)
        assert got == pytest.approx(slope, abs=1e-10)


def test_loglog_slope_ignores_zeros_and_needs_enough_pixels():
    flux = _power_law(WAVE, -2.0)
    switched_off = flux.copy()
    switched_off[: WAVE.size // 2] = 0.0  # the QP leaves exact zeros
    assert rel.loglog_slope(switched_off, LOG_WAVE) == pytest.approx(-2.0, abs=1e-10)

    nearly_all_off = flux.copy()
    nearly_all_off[rel.REVERSAL_MIN_PIXELS - 1:] = 0.0
    assert np.isnan(rel.loglog_slope(nearly_all_off, LOG_WAVE))


def test_reversal_state_flags_swapped_colours_only_when_both_families_carry_flux():
    moon = _power_law(WAVE, -3.7)
    zodi = _power_law(WAVE, -0.3)
    # Physical ordering: zodi is redder, so separation > 0 and nothing is flagged.
    is_reversed, testable, info = rel.reversal_state(
        {"moon": moon, "zodi": zodi}, LOG_WAVE)
    assert testable and not is_reversed
    assert info["separation"] == pytest.approx(3.4, abs=1e-8)

    # Swap the two: the moon is now the redder family -- a reversal.
    is_reversed, testable, info = rel.reversal_state(
        {"moon": zodi, "zodi": moon}, LOG_WAVE)
    assert testable and is_reversed
    assert info["separation"] == pytest.approx(-3.4, abs=1e-8)

    # One family switched off carries no ordering information.
    faint = moon * 1.0e-8
    is_reversed, testable, _info = rel.reversal_state(
        {"moon": faint, "zodi": zodi}, LOG_WAVE)
    assert not testable and not is_reversed


def test_reversal_state_is_indifferent_to_amplitude_only_to_shape():
    moon = _power_law(WAVE, -0.3, amplitude=1.0e3)
    zodi = _power_law(WAVE, -3.7, amplitude=1.0e3)
    a = rel.reversal_state({"moon": moon, "zodi": zodi}, LOG_WAVE)
    b = rel.reversal_state({"moon": moon * 137.0, "zodi": zodi * 137.0}, LOG_WAVE)
    assert a[0] is b[0] is True
    assert a[2]["separation"] == pytest.approx(b[2]["separation"], abs=1e-12)


def test_diffuse_collapsed_only_when_the_block_is_switched_off():
    bestfit = _power_law(WAVE, -1.0, amplitude=1.0e3)
    assert not rel.diffuse_collapsed({"diffuse": bestfit * 0.04}, bestfit)
    assert rel.diffuse_collapsed({"diffuse": bestfit * 1.0e-12}, bestfit)
    assert rel.diffuse_collapsed({"diffuse": np.zeros_like(bestfit)}, bestfit)
    # Falls back to the individual species when no summed plane is present.
    species = {k: bestfit * 1.0e-12 for k in rel.DIFFUSE_COMPONENT_KEYS}
    assert rel.diffuse_collapsed(species, bestfit)
    species["feo"] = bestfit * 0.05
    assert not rel.diffuse_collapsed(species, bestfit)


def test_flag_bits_are_distinct_powers_of_two_and_describe_round_trips():
    values = [value for value, _name in rel.RELIABILITY_BITS]
    assert len(set(values)) == len(values)
    for value in values:
        assert value > 0 and value & (value - 1) == 0
    assert rel.describe(0) == "ok"
    both = rel.RELIABILITY_REVERSAL_RETRIED | rel.RELIABILITY_REVERSAL_RECOVERED
    assert rel.describe(both) == "reversal_retried|reversal_recovered"
    assert "unknown" in rel.describe(1 << 30)


def test_extra_meta_columns_are_validated_and_pivoted():
    from skysub.sky_decomp.result_io import extra_meta_columns

    assert extra_meta_columns(None, 3) == {}
    assert extra_meta_columns([{"a": 1}, {"a": 2}], 2) == {"a": [1, 2]}

    with pytest.raises(ValueError, match="entries for"):
        extra_meta_columns([{"a": 1}], 2)
    with pytest.raises(ValueError, match="disagree on their columns"):
        extra_meta_columns([{"a": 1}, {"b": 2}], 2)
    with pytest.raises(ValueError, match="collides with a result field"):
        extra_meta_columns([{"reduced_chi2": 1.0}], 1,
                           reserved={"reduced_chi2": []})


def test_cached_rows_are_marked_not_evaluated_rather_than_clean():
    from skysub import decompose_parallel

    flags = {"reliability": np.int32(0), "reversal_separation": 1.5,
             "reversal_moon_frac": 0.4, "reversal_retry_bound": float("nan")}
    filled = decompose_parallel._reliability_extra_meta([flags, None], 2)
    assert filled[0] is flags
    # -1, never 0: a row whose flags were not computed must not read as clean.
    assert int(filled[1]["reliability"]) == -1
    assert np.isnan(filled[1]["reversal_separation"])
    assert set(filled[1]) == set(flags)
    assert decompose_parallel._reliability_extra_meta([None, None], 2) is None


def _prior(**overrides):
    prior = {
        "moon_fraction": 0.5, "amp_prior_tol": 3.0, "amp_prior_floor": 0.02,
        "zodi_total": 100.0, "zodi_amp_bound": 2.0,
        "diffuse_oh_amp": 10.0, "diffuse_oh_centre_log10": -0.6489,
        "diffuse_oh_bound_dex": 0.15, "diffuse_oh_relax_dex": 0.0,
        "diffuse_oh_gate_frac": 0.6,
        "diffuse_ratio_nominal": (0.1, 0.6, 0.3),
        "diffuse_ratio_bound_dex": 0.2,
        "moon_ratio_bound": 0.7, "zodi_ratio_bound": 0.7,
    }
    prior.update(overrides)
    return prior


def _flat(total, n=1200):
    return np.full(n, float(total) / n)


def test_zodi_anchor_bit_fires_only_on_the_bracket():
    # v == kappa_z * Z exactly -> on the ceiling.
    bits, info = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(200.0)}, _prior())
    assert bits & rel.RELIABILITY_ZODI_ANCHOR_PINNED
    assert info["zodi_int"] == pytest.approx(200.0)
    # v == Z / kappa_z exactly -> on the floor, also pinned.
    bits, _ = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(50.0)}, _prior())
    assert bits & rel.RELIABILITY_ZODI_ANCHOR_PINNED
    # interior
    bits, _ = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(120.0)}, _prior())
    assert not bits & rel.RELIABILITY_ZODI_ANCHOR_PINNED


def test_moon_share_bit_catches_the_dark_time_floor():
    # f = 0 collapses the bracket to [0, amp_prior_floor]; a dark-time fit sits
    # on 0.02 exactly -- the "ghost" moon block.
    components = {"moon": _flat(2.0), "zodi": _flat(98.0)}
    bits, info = rel.constraint_bits(components, _prior(moon_fraction=0.0))
    assert info["moon_share"] == pytest.approx(0.02)
    assert bits & rel.RELIABILITY_MOON_SHARE_PINNED
    bits, _ = rel.constraint_bits({"moon": _flat(30.0), "zodi": _flat(70.0)},
                                  _prior(moon_fraction=0.0))
    assert not bits & rel.RELIABILITY_MOON_SHARE_PINNED


def test_diffuse_oh_cap_bit_is_gated_on_moon_brightness():
    on_cap = 10.0 ** (-0.6489 + 0.15) * 10.0
    components = {"moon": _flat(100.0), "zodi": _flat(120.0),
                  "diffuse": _flat(on_cap)}
    bits, _ = rel.constraint_bits(components, _prior(moon_fraction=0.9))
    assert bits & rel.RELIABILITY_DIFFUSE_OH_CAP_BINDING
    # Below the gate there is no cap at all, so it cannot bind.
    bits, _ = rel.constraint_bits(components, _prior(moon_fraction=0.3))
    assert not bits & rel.RELIABILITY_DIFFUSE_OH_CAP_BINDING


def test_shape_bound_counts_pairs_and_ignores_switched_off_families():
    beta = 0.7
    on_bound = np.array([1.0, beta, beta ** 2, beta ** 3], dtype=float)
    bits, info = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(120.0)}, _prior(),
        coef_blocks={"moon": on_bound})
    assert info["shape_bound_pairs"] == 3
    assert bits & rel.RELIABILITY_SHAPE_BOUND_ACTIVE
    # An interior shape touches no bound.
    bits, info = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(120.0)}, _prior(),
        coef_blocks={"moon": np.array([1.0, 0.9, 0.85, 0.8])})
    assert info["shape_bound_pairs"] == 0
    assert not bits & rel.RELIABILITY_SHAPE_BOUND_ACTIVE
    # Zeros are a family switched OFF, not a shaped one.
    bits, info = rel.constraint_bits(
        {"moon": _flat(100.0), "zodi": _flat(120.0)}, _prior(),
        coef_blocks={"moon": np.zeros(4)})
    assert info["shape_bound_pairs"] == 0


def test_colour_excess_is_zero_for_identical_shapes_and_scale_free():
    wave = np.linspace(3600.0, 9800.0, 2000)
    sci = _power_law(wave, -1.0, amplitude=1.0e3)
    assert rel.sci_colour_excess(sci, sci, sci, wave) == pytest.approx(0.0, abs=1e-12)
    # A throughput difference cancels; only shape survives.
    assert rel.sci_colour_excess(sci * 7.0, sci, sci, wave) == pytest.approx(
        0.0, abs=1e-12)
    redder = _power_law(wave, -0.5, amplitude=1.0e3)
    assert rel.sci_colour_excess(redder, sci, sci, wave) > 0.0


def test_prediction_reliability_reports_not_evaluated_without_a_reference():
    from skysub.mlp_predictor.inference import (
        prediction_reliability, _upper_bound_by_name)

    names = ["HO2", "FeO", "O2Ac", "Moon_bs00"]
    coef = np.array([[1e-9, 1e-9, 1e-9, 1.0],
                     [0.1, 0.5, 2.0, 1.0]])
    # No reference at all: -1 per row, NEVER 0 -- a test that did not run must
    # not read as a clean row.
    assert np.all(prediction_reliability(coef, names) == -1)

    upper = {"continuum": np.array([0.2, 2.2, 6.0]), "moon": np.array([3.0])}
    groups = {"continuum": [0, 1, 2], "moon": [3]}
    assert _upper_bound_by_name(names, upper, groups)["FeO"] == pytest.approx(2.2)
    bits = prediction_reliability(coef, names, upper, groups)
    assert bits[0] == rel.RELIABILITY_DIFFUSE_COLLAPSED
    assert bits[1] == 0
