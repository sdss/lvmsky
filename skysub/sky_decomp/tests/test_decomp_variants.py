"""The ML side reads, and reconstructs with, the one decomposition the pipeline makes.

The solarblue coefficients are nearly identical to the older palace fits, so a
reconstruction on the wrong class or OH file is invisible in coefficient space;
these tests pin the plumbing that prevents it.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mlp_predictor import data  # noqa: E402
from mlp_predictor.config import DataConfig  # noqa: E402
from skysub.sky_decomp.moon_zodi_model import SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX  # noqa: E402
from skysub import decompose_parallel  # noqa: E402


def test_suffix_matches_decompose_parallel():
    assert data.DECOMP_SUFFIX == decompose_parallel.FIT_MODEL_SUFFIX
    assert DataConfig().decomp_suffix == data.DECOMP_SUFFIX
    assert data.DECOMP_SUFFIX.endswith("_solarblue")


def test_oh_file_matches_the_fit():
    assert data.PALACE_OH_SUFFIX == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    assert decompose_parallel.PALACE_OH_SUFFIX == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX


def _bundle():
    return {"telluric_calculator": object(), "pwv_mm": 5.0, "source_airmass": 1.2,
            "drp_transmission": np.ones(3)}


@pytest.fixture
def captured(monkeypatch):
    import sky_decomp.solar_lsf as sl
    calls = []

    class Recorder:
        def __init__(self, wave, **kwargs):
            calls.append(kwargs)

    monkeypatch.setattr(sl, "SkyDecompPalaceCorrSolarBlueLSF", Recorder)
    return calls


def test_reconstruction_uses_the_solarblue_class_and_the_bundle(captured):
    data.make_reconstruction_decomposer(
        np.arange(3.0), n_spline_knots=11, base_dir=".", telluric=_bundle())
    kw = captured[-1]
    assert kw["palace_oh_suffix"] == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    assert kw["pwv_mm"] == 5.0 and kw["source_airmass"] == 1.2


def test_reconstruction_requires_a_telluric_bundle():
    with pytest.raises(ValueError, match="telluric"):
        data.make_reconstruction_decomposer(
            np.arange(3.0), n_spline_knots=11, base_dir=".", telluric=None)
