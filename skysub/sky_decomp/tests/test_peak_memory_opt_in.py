"""Allocation tracing for ``peak_memory_mb`` is opt-in: it multiplies the fit time."""

import tracemalloc
from pathlib import Path

import numpy as np
import pytest

from skysub.sky_decomp.fit import SkyDecomp, peak_memory_requested
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.moon_zodi_model import DEFAULT_DATA_ROOT
from skysub.sky_decomp.telluric_corrected_lines import SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D

GOLDEN = Path(__file__).parent / "data/lsf_surface_iterative_row837_n5_golden.npz"


class _FakeTelluric:
    """Stand-in for the DRP calculator; the line model never calls it."""

    def match_to_data(self, wave, lsf, pwv, *, airmass, lsf_in_wavelength):
        return np.full_like(wave, 0.85, dtype=float)


def _plain_fit():
    """``SkyDecomp.fit`` on a small synthetic spectrum."""
    wave = np.linspace(4000.0, 9500.0, 400)
    model = SkyDecomp(wave, lsf_sigma=1.5, n_spline_knots=25)
    flux = 1e-16 * (model.matrix_moon.sum(0) + 0.3 * model.matrix_oh.sum(0))
    return model.fit(flux, np.full_like(wave, (0.01 * np.median(flux)) ** -2), verbose=False)


def _iterative_fit():
    """``SkyDecompLSFSurfaceIterative.fit`` on the golden row with one refinement cycle."""
    golden = np.load(GOLDEN)
    wave = np.asarray(golden["wave"], dtype=float)
    model = SkyDecompPalaceCorrAijcVNFSplitZodiLSFSpline2D(
        wave,
        telluric_calculator=_FakeTelluric(),
        pwv_mm=4.0,
        source_airmass=1.3,
        drp_transmission=np.full_like(wave, 0.85),
        lsf_sigma=0.5,
        base_dir=DEFAULT_DATA_ROOT,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=1, roughness_fraction=1.0e-4),
    )
    return model.fit(np.asarray(golden["flux"], dtype=float), np.asarray(golden["ivar"], dtype=float))


@pytest.mark.parametrize("value, expected", [(None, False), ("", False), ("0", False), ("1", True), ("yes", True)])
def test_peak_memory_requested_reads_environment(monkeypatch, value, expected):
    if value is None:
        monkeypatch.delenv("LVMSKY_MEASURE_PEAK_MEMORY", raising=False)
    else:
        monkeypatch.setenv("LVMSKY_MEASURE_PEAK_MEMORY", value)
    assert peak_memory_requested() is expected


@pytest.mark.parametrize("fit", [_plain_fit, _iterative_fit], ids=["plain", "iterative"])
def test_default_fit_leaves_tracemalloc_off(monkeypatch, fit):
    monkeypatch.delenv("LVMSKY_MEASURE_PEAK_MEMORY", raising=False)
    assert not tracemalloc.is_tracing()
    traced = []
    monkeypatch.setattr(tracemalloc, "start", lambda *args: traced.append(args))
    result = fit()
    assert traced == []
    assert not tracemalloc.is_tracing()
    assert np.isnan(result.peak_memory_mb)


@pytest.mark.parametrize("fit", [_plain_fit, _iterative_fit], ids=["plain", "iterative"])
def test_requested_fit_reports_peak_memory(monkeypatch, fit):
    monkeypatch.setenv("LVMSKY_MEASURE_PEAK_MEMORY", "1")
    assert not tracemalloc.is_tracing()
    result = fit()
    # The fit stops the tracing it started.
    assert not tracemalloc.is_tracing()
    assert np.isfinite(result.peak_memory_mb) and result.peak_memory_mb > 0.0
