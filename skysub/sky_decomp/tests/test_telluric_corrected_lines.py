from pathlib import Path

import numpy as np
import pytest

from skysub.sky_decomp.lsf_spline2d import SkyDecompLSFSpline2D
from skysub.sky_decomp.lsf_surface_iterative import LSFSurfaceIterativeConfig
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PALACE_OH_SUFFIX,
)
from skysub.sky_decomp.telluric_corrected_lines import (
    SkyDecompAdam25kNivContinuumLSFSpline2D,
    SkyDecompAdam25kTelluricLSFSpline2D,
    SkyDecompTelluricCorrectedLinesLSFSpline2D,
    SkyDecompTelluricLinesLSFSpline2D,
    calculate_drp_transmission,
    calculate_line_transmission,
    restore_drp_input,
)


GOLDEN = Path(__file__).parent / "data/lsf_surface_iterative_row837_n5_golden.npz"


class _FakeTelluric:
    def __init__(self, wave_air, transmission):
        self.wave_air = np.asarray(wave_air, dtype=float)
        self.transmission = np.asarray(transmission, dtype=float)
        self.last_lsf = None
        self.calc_calls = 0

    def calc_transmission(self, pwv, *, airmass):
        assert pwv > 0.0 and airmass > 0.0
        self.calc_calls += 1
        return self.transmission.copy()

    def match_to_data(self, wave, lsf, pwv, *, airmass, lsf_in_wavelength):
        assert pwv > 0.0 and airmass > 0.0 and lsf_in_wavelength
        self.last_lsf = np.asarray(lsf).copy()
        return np.full_like(wave, 0.8, dtype=float)


def _wave():
    return np.asarray(np.load(GOLDEN)["wave"], dtype=float)


def _model_kwargs():
    return {
        "lsf_sigma": 0.5,
        "base_dir": DEFAULT_DATA_ROOT,
        "n_spline_knots": 25,
        "n_zodi_spline_knots": 3,
        "zodi_smooth_lambda": 0.1,
        "moon_smooth_lambda": 0.1,
        "config": LSFSurfaceIterativeConfig(
            n_refinement_cycles=1,
            roughness_fraction=1.0e-4,
        ),
    }


def test_restore_drp_input_preserves_the_native_objective():
    wave = np.arange(4.0)
    flux = np.array([1.0, 2.0, 3.0, 4.0])
    ivar = np.ones_like(flux)
    lsf = np.array([[1.0, 2.0, np.nan, 4.0], [3.0, 4.0, 6.0, 8.0]])
    telluric = _FakeTelluric(wave, np.ones_like(wave))

    flux_pre, ivar_pre, transmission = restore_drp_input(
        flux, ivar, wave, lsf, 5.0, 1.2, telluric
    )
    np.testing.assert_array_equal(telluric.last_lsf, np.nanmedian(lsf, axis=0))
    model_pre = np.array([0.7, 1.6, 2.5, 3.4])
    np.testing.assert_allclose(
        np.sum((flux_pre - model_pre) ** 2 * ivar_pre),
        np.sum((flux - model_pre / transmission) ** 2 * ivar),
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    np.testing.assert_array_equal(
        calculate_drp_transmission(wave, lsf, 5.0, 1.2, telluric),
        transmission,
    )
    with pytest.raises(ValueError, match="pwv_mm"):
        restore_drp_input(flux, ivar, wave, lsf, 0.0, 1.2, telluric)


def test_palace_line_formula_matches_the_published_equation():
    transmission_ref = np.array([0.98, 0.72, 0.31])
    fraction_h2o = np.array([0.0, 0.4, 1.0])
    tau_ref = -np.log(transmission_ref)
    pwv_mm = 6.0
    airmass = 1.7

    actual = calculate_line_transmission(
        tau_ref * (1.0 - fraction_h2o),
        tau_ref * fraction_h2o,
        pwv_mm,
        airmass,
    )
    expected = transmission_ref ** (
        (1.0 + (pwv_mm / 2.5 - 1.0) * fraction_h2o) * airmass
    )
    np.testing.assert_allclose(actual, expected, rtol=2.0e-15, atol=0.0)


def test_all_line_families_use_r4m_line_coefficients_and_new_oh_default():
    wave = _wave()
    wave_hr = np.linspace(wave[0] - 10.0, wave[-1] + 10.0, 40_000)
    transmission_hr = 0.65 + 0.3 * (wave_hr - wave_hr[0]) / np.ptp(wave_hr)
    telluric = _FakeTelluric(wave_hr, transmission_hr)
    intrinsic = SkyDecompLSFSpline2D(wave, **_model_kwargs())
    candidate = SkyDecompTelluricLinesLSFSpline2D(
        wave,
        **_model_kwargs(),
        telluric_calculator=telluric,
        pwv_mm=5.0,
        line_airmass=1.4,
    )

    assert candidate.palace_oh_suffix == DEFAULT_PALACE_OH_SUFFIX
    assert candidate._pmd_path("pmd_popmodel_OH.dat").name == (
        f"pmd_popmodel_OH{DEFAULT_PALACE_OH_SUFFIX}.dat"
    )
    np.testing.assert_array_equal(candidate._line_wave, intrinsic._line_wave)
    np.testing.assert_allclose(
        candidate._line_weights(),
        intrinsic._line_weights() * candidate._line_transmission_values,
        rtol=1.0e-14,
        atol=0.0,
    )
    assert telluric.calc_calls == 0
    assert np.any(
        candidate._line_transmission_values
        != np.interp(candidate._line_wave, wave_hr, transmission_hr)
    )
    for family in ("oh", "atom", "orc", "o2"):
        family_slice = candidate._group_slices[family]
        line_mask = (candidate._line_group >= family_slice.start) & (
            candidate._line_group < family_slice.stop
        )
        assert np.any(
            candidate._line_weights()[line_mask] != intrinsic._line_weights()[line_mask]
        )


def test_corrected_domain_operator_divides_only_after_lsf():
    wave = _wave()
    wave_hr = np.linspace(wave[0] - 10.0, wave[-1] + 10.0, 40_000)
    telluric = _FakeTelluric(wave_hr, np.linspace(0.7, 0.9, wave_hr.size))
    drp_transmission = np.linspace(0.72, 0.94, wave.size)
    pre = SkyDecompTelluricLinesLSFSpline2D(
        wave,
        **_model_kwargs(),
        telluric_calculator=telluric,
        pwv_mm=5.0,
        line_airmass=1.35,
    )
    corrected = SkyDecompTelluricCorrectedLinesLSFSpline2D(
        wave,
        **_model_kwargs(),
        telluric_calculator=telluric,
        pwv_mm=5.0,
        source_airmass=1.35,
        drp_transmission=drp_transmission,
    )

    for family in ("oh", "atom", "orc"):
        np.testing.assert_allclose(
            getattr(corrected, f"matrix_{family}"),
            getattr(pre, f"matrix_{family}") / drp_transmission,
            rtol=2.0e-15,
            atol=0.0,
        )
    for channel, matrix in pre._line_components.items():
        probe = np.linspace(0.5, 1.5, matrix.shape[1])
        np.testing.assert_allclose(
            corrected._line_components[channel] @ probe,
            (matrix @ probe) / drp_transmission,
            rtol=2.0e-14,
            atol=1.0e-14,
        )


def test_named_adam25k_class_locks_its_oh_asset():
    wave = _wave()
    wave_hr = np.linspace(wave[0] - 10.0, wave[-1] + 10.0, 40_000)
    model = SkyDecompAdam25kTelluricLSFSpline2D(
        wave,
        **_model_kwargs(),
        telluric_calculator=_FakeTelluric(wave_hr, np.ones_like(wave_hr)),
        pwv_mm=2.0,
        source_airmass=1.1,
        drp_transmission=np.ones_like(wave),
    )

    assert model.palace_oh_suffix == DEFAULT_PALACE_OH_SUFFIX
    with pytest.raises(ValueError, match="requires its bundled OH source table"):
        SkyDecompAdam25kTelluricLSFSpline2D(
            wave,
            **(_model_kwargs() | {"palace_oh_suffix": "_other"}),
            telluric_calculator=_FakeTelluric(wave_hr, np.ones_like(wave_hr)),
            pwv_mm=2.0,
            source_airmass=1.1,
            drp_transmission=np.ones_like(wave),
        )


def test_niv_adam_class_locks_the_continuum_contract():
    wave = _wave()
    wave_hr = np.linspace(wave[0] - 10.0, wave[-1] + 10.0, 40_000)
    common = {
        **_model_kwargs(),
        "telluric_calculator": _FakeTelluric(wave_hr, np.ones_like(wave_hr)),
        "pwv_mm": 2.0,
        "source_airmass": 1.1,
        "drp_transmission": np.ones_like(wave),
    }
    common.pop("n_spline_knots", None)
    common.pop("n_zodi_spline_knots", None)
    model = SkyDecompAdam25kNivContinuumLSFSpline2D(wave, **common)

    assert model.split_zodi is True
    assert model.n_spline_knots == 11
    assert model.n_zodi_spline_knots == 1
    assert model.diffuse_oh_scope == "block"
    with pytest.raises(ValueError, match="requires moon_ratio_bound=0.7"):
        SkyDecompAdam25kNivContinuumLSFSpline2D(
            wave,
            **(common | {"moon_ratio_bound": 0.6}),
        )
