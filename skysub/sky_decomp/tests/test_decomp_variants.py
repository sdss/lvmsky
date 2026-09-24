"""The ML side picks the OH line file from the decomposition suffix.

The palacecorr coefficients are nearly identical to the palace ones, so a
reconstruction on the wrong OH file is invisible in coefficient space; these
tests pin the plumbing that prevents it.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from mlp_predictor import data  # noqa: E402
from skysub.sky_decomp.moon_zodi_model import SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX  # noqa: E402
from skysub import decompose_parallel  # noqa: E402

PALACE = "_palace_aijc_vnf_split_zodi_lsf_spline2d"
PALACECORR = "_palacecorr_aijc_vnf_split_zodi_lsf_spline2d"


def test_variant_suffixes_match_decompose_parallel():
    fit = decompose_parallel.FIT_MODEL_SUFFIXES
    assert fit[decompose_parallel.PALACE_VNF_SPLIT_ZODI_FIT_MODEL] == PALACE
    assert fit[decompose_parallel.PALACECORR_VNF_SPLIT_ZODI_FIT_MODEL] == PALACECORR
    assert data.DECOMP_VARIANTS["telluric"]["suffix"] == PALACE
    assert data.DECOMP_VARIANTS["telluric-palacecorr"]["suffix"] == PALACECORR


def test_palacecorr_oh_file_matches_the_fit():
    assert data.palace_oh_suffix_for(PALACECORR) == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    assert data.palace_oh_suffix_for(PALACE) is None


def test_legacy_suffix_is_non_telluric_and_unknown_telluric_raises():
    assert data.decomp_variant_spec("_lsf_surface_iterative_split_zodi") is None
    assert data.decomp_variant_spec("") is None
    assert data.decomp_variant_spec(None) is None
    with pytest.raises(ValueError, match="not registered"):
        data.decomp_variant_spec("_palacenew_aijc_vnf_split_zodi_lsf_spline2d")


def test_lookup_requires_a_telluric_suffix():
    with pytest.raises(TypeError):
        data.make_telluric_row_lookup("unused.fits")
    with pytest.raises(ValueError, match="not a telluric"):
        data.make_telluric_row_lookup("unused.fits",
                                      decomp_suffix="_lsf_surface_iterative_split_zodi")


def _bundle(oh):
    return {"telluric_calculator": object(), "pwv_mm": 5.0, "source_airmass": 1.2,
            "drp_transmission": np.ones(3), "palace_oh_suffix": oh}


@pytest.fixture
def captured(monkeypatch):
    import sky_decomp.residual_pca as rp
    calls = []

    class Recorder:
        def __init__(self, wave, **kwargs):
            calls.append(kwargs)

    monkeypatch.setattr(rp, "SkyDecompPalaceAijcVNFSplitZodiLSFSpline2D", Recorder)
    return calls


def test_bundle_oh_file_reaches_the_constructor(captured):
    data.make_reconstruction_decomposer(
        np.arange(3.0), n_spline_knots=11, base_dir=".",
        telluric=_bundle(SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX))
    assert captured[-1]["palace_oh_suffix"] == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    data.make_reconstruction_decomposer(
        np.arange(3.0), n_spline_knots=11, base_dir=".", telluric=_bundle(None))
    assert captured[-1]["palace_oh_suffix"] is None


def test_explicit_oh_file_may_agree_but_not_contradict(captured):
    bundle = _bundle(SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX)
    data.make_reconstruction_decomposer(
        np.arange(3.0), n_spline_knots=11, base_dir=".", telluric=bundle,
        palace_oh_suffix=SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX)
    assert captured[-1]["palace_oh_suffix"] == SKYFAR_LINEAR_RIDGE_PALACE_OH_SUFFIX
    with pytest.raises(ValueError, match="disagrees"):
        data.make_reconstruction_decomposer(
            np.arange(3.0), n_spline_knots=11, base_dir=".", telluric=bundle,
            palace_oh_suffix="_other")
    assert "palace_oh_suffix" in bundle   # the caller's bundle is not mutated
