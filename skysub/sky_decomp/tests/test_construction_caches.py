"""Contracts behind the fast per-row construction of the telluric models.

`decompose_parallel` builds a fresh decomposer for every fitted row of a
telluric fit model (the DRP transmission and the source airmass are per-row),
where the non-telluric models build one decomposer per worker.  Construction
is therefore on the hot path, and it is made cheap by a windowed line-profile
sum plus caches on the parts that cannot vary from row to row.  These tests
pin the two properties that make that safe: the window is exact to float64,
and nothing row-dependent is ever shared between instances.
"""

from pathlib import Path

import numpy as np
import pytest

from skysub.sky_decomp.fit import (
    GRP2VECTOR_SIGMA_CUTOFF,
    SPLIT_ZODI_CONTINUUM_DEFAULTS,
    grp2vector,
)
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT
from skysub.sky_decomp.telluric_corrected_lines import (
    SkyDecompAdam25kTelluricSplitZodiLSFSpline2D,
)

GOLDEN = (
    Path(__file__).parent / "data/lsf_surface_iterative_row837_n5_golden.npz"
)


class _FakeTelluric:
    """Stand-in for the DRP calculator; the line model never calls it."""

    def match_to_data(self, wave, lsf, pwv, *, airmass, lsf_in_wavelength):
        return np.full_like(wave, 0.85, dtype=float)


def _dense_grp2vector(line_wave, line_amp, wave, lsf_sigma):
    """The unwindowed form this routine replaced, kept as the reference."""
    cent = np.asarray(line_wave, float)
    amp = np.asarray(line_amp, float)
    sig = (
        np.interp(cent, wave, lsf_sigma)
        if np.ndim(lsf_sigma) > 0
        else float(lsf_sigma)
    )
    yy = (wave[:, None] - cent) / sig
    return np.sum(amp[None, :] * np.exp(-0.5 * yy**2), axis=1)


def _wave():
    return np.asarray(np.load(GOLDEN)["wave"], dtype=float)


@pytest.mark.parametrize(
    "n_lines, sigma",
    [(1, 0.5), (7, 0.15), (70, 0.5), (500, 1.2), (3000, None)],
)
def test_windowed_line_sum_matches_the_dense_sum_to_float64(n_lines, sigma):
    wave = _wave()
    rng = np.random.default_rng(20260917)
    cent = rng.uniform(wave[0] - 5.0, wave[-1] + 5.0, n_lines)
    amp = rng.uniform(0.0, 1.0, n_lines)
    # `None` exercises the per-pixel sigma array, which is interpolated per line.
    lsf_sigma = (
        np.linspace(0.4, 1.4, wave.size) if sigma is None else sigma
    )

    dense = _dense_grp2vector(cent, amp, wave, lsf_sigma)
    windowed = grp2vector(cent, amp, wave, lsf_sigma)

    # The two differ only by summation order and by tails the window drops,
    # both of which are below the float64 resolution of the peak.
    assert np.max(np.abs(dense - windowed)) <= 8.0 * np.finfo(float).eps * dense.max()
    assert np.allclose(np.sum(dense), np.sum(windowed), rtol=1.0e-14, atol=0.0)


def test_windowed_line_sum_handles_degenerate_input():
    wave = _wave()
    empty = grp2vector(np.array([]), np.array([]), wave, 0.5)
    assert empty.shape == wave.shape and not np.any(empty)
    # A line well outside the grid reaches no pixel at all.
    assert not np.any(grp2vector(np.array([100.0]), np.array([1.0]), wave, 0.5))
    # A line just off the blue end still reaches the pixels inside the window.
    edge = grp2vector(np.array([wave[0] - 1.0]), np.array([1.0]), wave, 0.5)
    assert edge[0] > 0.0
    assert GRP2VECTOR_SIGMA_CUTOFF >= 8.0


def _telluric_model(pwv_mm, source_airmass, drp_transmission):
    return SkyDecompAdam25kTelluricSplitZodiLSFSpline2D(
        _wave(),
        telluric_calculator=_FakeTelluric(),
        pwv_mm=pwv_mm,
        source_airmass=source_airmass,
        drp_transmission=drp_transmission,
        lsf_sigma=0.5,
        base_dir=DEFAULT_DATA_ROOT,
        moon_smooth_lambda=0.1,
        moon_interline_boost=0.0,
        config=LSFSurfaceIterativeConfig(
            n_refinement_cycles=1, roughness_fraction=1.0e-4
        ),
    )


def test_per_row_reconstruction_shares_only_row_independent_state():
    wave = _wave()
    span = (wave - wave[0]) / np.ptp(wave)
    first_row = (4.0, 1.3, 0.80 + 0.15 * span)
    other_row = (2.0, 1.8, 0.60 + 0.35 * span)

    first = _telluric_model(*first_row)
    # A different row in between: if anything row-dependent were cached, the
    # rebuild below would come back carrying this row's transmission.
    other = _telluric_model(*other_row)
    rebuilt = _telluric_model(*first_row)

    for name in (
        "_line_transmission_values",
        "matrix_oh",
        "matrix_atom",
        "matrix_orc",
        "design_matrix",
        "drp_transmission",
    ):
        np.testing.assert_array_equal(
            getattr(rebuilt, name), getattr(first, name), err_msg=name
        )
        assert not np.allclose(
            getattr(other, name), getattr(first, name)
        ), f"{name} must depend on the row"

    for channel in first._line_components:
        probe = np.ones(first._line_components[channel].shape[1])
        np.testing.assert_array_equal(
            rebuilt._line_components[channel] @ probe,
            first._line_components[channel] @ probe,
        )
        # The continuum M-spline integrals carry no transmission at all, so
        # every row shares the very same matrix.
        shared = other._continuum_components[channel]
        mine = rebuilt._continuum_components[channel]
        assert mine.shape == shared.shape
        np.testing.assert_array_equal(mine.data, shared.data)
        np.testing.assert_array_equal(mine.indices, shared.indices)
        np.testing.assert_array_equal(mine.indptr, shared.indptr)

    # Row-independent families come from the memoised assets unchanged.
    for name in ("vector_moon", "matrix_diffuse", "matrix_moon", "matrix_zodi"):
        np.testing.assert_array_equal(
            getattr(other, name), getattr(first, name), err_msg=name
        )


def test_per_row_reconstruction_gives_a_bit_identical_fit():
    golden = np.load(GOLDEN)
    wave = np.asarray(golden["wave"], dtype=float)
    flux = np.asarray(golden["flux"], dtype=float)
    ivar = np.asarray(golden["ivar"], dtype=float)
    row = (4.0, 1.3, np.full_like(wave, 0.85))

    first = _telluric_model(*row).fit(flux, ivar)
    _telluric_model(2.0, 1.8, np.full_like(wave, 0.7))
    rebuilt = _telluric_model(*row).fit(flux, ivar)

    np.testing.assert_array_equal(rebuilt.coef, first.coef)
    np.testing.assert_array_equal(rebuilt.bestfit, first.bestfit)
    assert rebuilt.r2 == first.r2


def test_cached_line_catalog_requires_a_class_level_amplitude_rule():
    from skysub.sky_decomp.lsf_spline2d import SkyDecompLSFSpline2D

    class _InstanceAmplitude(SkyDecompLSFSpline2D):
        def _oh_amplitude(self, group):  # not a staticmethod
            return np.asarray(group["Aij"] * group["gi"], dtype=float)

    with pytest.raises(TypeError, match="must be a staticmethod"):
        _InstanceAmplitude(
            _wave(),
            **(
                SPLIT_ZODI_CONTINUUM_DEFAULTS
                | {
                    "lsf_sigma": 0.5,
                    "base_dir": DEFAULT_DATA_ROOT,
                    "moon_smooth_lambda": 0.1,
                }
            ),
        )
