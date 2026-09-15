from dataclasses import asdict
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits

import skysub.decompose_parallel as parallel
from skysub.sky_decomp.lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    LSFSurfaceState,
)
from skysub.sky_decomp.result_io import load_lsf_surface_state


def _result():
    channels = ("B", "R", "Z")
    state = LSFSurfaceState(
        coefficients={name: np.ones((11, 1)) for name in channels},
        knot_vectors={name: np.array([3600.0, 9800.0]) for name in channels},
        degrees={name: 0 for name in channels},
        channel_bounds={name: (None, None) for name in channels},
        tap_offsets=np.arange(-5, 6),
        config=asdict(LSFSurfaceIterativeConfig()),
        metrics={name: {"status": "ok"} for name in channels},
        requested_cycles=5,
        completed_cycles=5,
        wave_n=12_401,
        wave_min=3600.0,
        wave_max=9800.0,
        wave_sha256="wave",
        fit_status="Solved",
    )
    values = {
        "t_o2": 1.0,
        "t_o2_err": 0.1,
        "o2_prefit_amp": 1.0,
        "reduced_chi2": 1.0,
        "r2": 0.9,
        "rms_resid": 0.2,
        "resid_level": 0.0,
        "fit_status": "Solved",
        "fit_summary": "ok",
        "fit_elapsed_sec": 1.0,
        "peak_memory_mb": 2.0,
        "o2_fit_status": "Solved",
        "o2_fit_summary": "ok",
        "o2_fit_elapsed_sec": 0.1,
        "o2_valid_frac": 1.0,
    }
    return SimpleNamespace(
        **values,
        design_names=["OH_000", "O2_b01"],
        coef=np.array([1.0, 2.0]),
        coef_err=np.array([0.1, 0.2]),
        coef_cov_moon=None,
        coef_cov_zodi=None,
        lsf_state=state,
    )


def test_compact_cache_preserves_rows_and_marks_failures(tmp_path, monkeypatch):
    monkeypatch.setattr(parallel, "_WORKER_COMPACT_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(parallel, "_WORKER_RUN_FINGERPRINT", "run")
    parallel._save_compact_cache("sci", 0, result=_result())
    try:
        raise ValueError("invalid PWV")
    except ValueError as error:
        parallel._save_compact_cache("sci", 1, error=error)

    output = tmp_path / "compact.fits"
    parallel._write_compact_fits(
        tmp_path / "cache", "sci", 2, "run", "test-method", output
    )
    with fits.open(output) as hdul:
        assert [hdu.name for hdu in hdul[:7]] == [
            "PRIMARY",
            "META",
            "COEF",
            "COEF_ERR",
            "LSF_COEF",
            "LSF_KNOTS",
            "LSF_META",
        ]
        np.testing.assert_allclose(hdul["COEF"].data["OH_000"], [1.0, np.nan])
        assert hdul["META"].data["input_valid"].tolist() == [True, False]
        assert "invalid PWV" in hdul["META"].data["error_message"][1]

    assert load_lsf_surface_state(output, 0).completed_cycles == 5
    with pytest.raises(ValueError, match="No fitted LSF state"):
        load_lsf_surface_state(output, 1)
