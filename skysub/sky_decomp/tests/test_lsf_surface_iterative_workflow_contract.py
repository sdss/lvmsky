from dataclasses import asdict
from types import MethodType

import numpy as np
import pytest

import skysub.sky_decomp.lsf_surface_iterative as iterative_module
from skysub.sky_decomp.lsf_surface_iterative import (
    LSFSurfaceIterativeConfig,
    LSFSurfaceState,
    SkyDecompLSFSurfaceIterative,
)


def _matrix_bundle(scale, size):
    row = np.full((1, size), float(scale))
    return {name: row.copy() for name in ("oh", "moon", "diffuse", "atom", "orc", "o2")}


def _lsf_state(config, scale, *, status="fit", knots=None):
    knots = np.array([0.0, 1.0]) if knots is None else np.asarray(knots, dtype=float)
    return LSFSurfaceState(
        coefficients={"B": np.array([[float(scale)]])},
        knot_vectors={"B": knots.copy()},
        degrees={"B": 0},
        channel_bounds={"B": (None, None)},
        tap_offsets=np.array([0]),
        config=asdict(config),
        metrics={"B": {"status": status}},
        requested_cycles=config.n_refinement_cycles,
        completed_cycles=0,
        wave_n=8,
        wave_min=5000.0,
        wave_max=5007.0,
        wave_sha256="workflow-test",
    )


@pytest.fixture
def workflow_model(monkeypatch):
    """Build a small deterministic substitute for the expensive fit internals."""

    monkeypatch.setattr(
        iterative_module,
        "build_skyline_mask",
        lambda wave, matrices, **kwargs: np.zeros(wave.size, dtype=bool),
    )
    monkeypatch.setattr(
        iterative_module,
        "continuum_fit_weights",
        lambda wave, residual, ivar, skyline_mask, **kwargs: (
            ivar.copy(),
            {"B": 1.0},
        ),
    )

    def build(*, cycles, solver_steps, lsf_states):
        model = object.__new__(SkyDecompLSFSurfaceIterative)
        model.wave = np.arange(5000.0, 5008.0)
        model.config = LSFSurfaceIterativeConfig(n_refinement_cycles=cycles)
        model.lsf_sigma = 0.5
        initial = _matrix_bundle(1.0, model.wave.size)
        for name, matrix in initial.items():
            setattr(model, f"matrix_{name}", matrix)
        for name in ("oh", "atom", "orc", "o2"):
            setattr(model, f"matrix_{name}_stick", initial[name].copy())

        model.lsf_surface_state = None
        model._lsf_surface = None
        model.lsf_kernels = {}
        model.lsf_metrics = {}
        model.design_names = list(initial)
        model.t_o2 = 191.5
        model.t_o2_err = 1.0
        model.o2_fit_status = "Solved"
        model.o2_fit_summary = "stub"
        model.o2_fit_elapsed_sec = 0.0
        model.o2_valid_frac = 1.0
        model.o2_prefit_amp = np.nan
        model.moon_knots_used = np.array([], dtype=float)
        model.moon_boosted_pixels_used = np.array([], dtype=float)

        steps = iter(solver_steps)
        states = iter(lsf_states)

        def fit_design(self, design, flux, ivar, **kwargs):
            status, coefficient = next(steps)
            coefficient = np.asarray(coefficient, dtype=float)
            return {
                "status": status,
                "coef": coefficient,
                "bestfit": np.full_like(flux, 123.0),
                "resid": np.full_like(flux, -5.0),
                "resid_level": -15.0,
            }

        def set_lsf_state(self, state):
            self.lsf_surface_state = state
            if state is None:
                self._lsf_surface = None
                self.lsf_kernels = {}
                self.lsf_metrics = {}
                return
            scale = float(state.coefficients["B"][0, 0])
            self._lsf_surface = np.full((self.wave.size, 1), scale)
            self.lsf_kernels = {"B": np.array([scale])}
            self.lsf_metrics = state.metrics.copy()

        def fit_lsf_channels(self, flux, ivar, source, fixed_background):
            state = next(states)
            self._set_lsf_state(state)
            return state

        def assemble_refined_matrices(self):
            scale = float(self.lsf_surface_state.coefficients["B"][0, 0])
            return _matrix_bundle(scale, self.wave.size)

        def components_from_coef(self, coefficient, matrices):
            slices = self._component_slices(matrices)
            return {name: matrices[name].T @ coefficient[slices[name]] for name in matrices}

        model._prefit_o2 = MethodType(lambda self, flux, ivar: None, model)
        model._fit_design = MethodType(fit_design, model)
        model._set_lsf_state = MethodType(set_lsf_state, model)
        model._fit_lsf_channels = MethodType(fit_lsf_channels, model)
        model._assemble_refined_matrices = MethodType(
            assemble_refined_matrices,
            model,
        )
        model._components_from_coef = MethodType(components_from_coef, model)
        return model

    return build


def test_nominal_seed_failure_skips_refinement_and_final_closure(workflow_model):
    model = workflow_model(
        cycles=1,
        solver_steps=[("PrimalInfeasible", [1, 2, 3, 4, 5, 6])],
        lsf_states=[],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    assert result.fit_status == "failed:nominal_seed:PrimalInfeasible"
    assert result.lsf_state.completed_cycles == 0
    assert result.lsf_state.failure_reason == "nominal_seed:PrimalInfeasible"
    assert result.lsf_state.final_continuum_status == "not_run"
    assert result.lsf_state.final_line_status == "not_run"
    assert [metric["stage"] for metric in model.stage_metrics] == ["00_nominal_joint_seed"]


def test_continuum_failure_keeps_last_cycle_and_runs_final_closure(workflow_model):
    config = LSFSurfaceIterativeConfig(n_refinement_cycles=2)
    accepted = _lsf_state(config, 2.0)
    model = workflow_model(
        cycles=2,
        solver_steps=[
            ("Solved", [1, 1, 1, 1, 1, 1]),
            ("Solved", [2, 3]),
            ("Solved", [4, 5, 6, 7]),
            ("PrimalInfeasible", [0, 0]),
            ("Solved", [8, 9]),
            ("Solved", [10, 11, 12, 13]),
        ],
        lsf_states=[accepted],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    assert result.lsf_state is accepted
    assert result.lsf_state.completed_cycles == 1
    assert result.lsf_state.failure_reason == "cycle_2_continuum:PrimalInfeasible"
    assert [metric["stage"] for metric in model.stage_metrics][-2:] == [
        "final_continuum",
        "final_lines",
    ]
    np.testing.assert_array_equal(result.bestfit_lsf, np.full(8, 126.0))


def test_line_failure_rolls_back_cycle_state_and_refined_matrices(workflow_model):
    config = LSFSurfaceIterativeConfig(n_refinement_cycles=2)
    accepted = _lsf_state(config, 2.0)
    rejected = _lsf_state(config, 3.0)
    model = workflow_model(
        cycles=2,
        solver_steps=[
            ("Solved", [1, 1, 1, 1, 1, 1]),
            ("Solved", [2, 3]),
            ("Solved", [4, 5, 6, 7]),
            ("Solved", [20, 30]),
            ("PrimalInfeasible", [0, 0, 0, 0]),
            ("Solved", [8, 9]),
            ("Solved", [10, 11, 12, 13]),
        ],
        lsf_states=[accepted, rejected],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    assert result.lsf_state is accepted
    assert result.lsf_state.completed_cycles == 1
    assert result.lsf_state.failure_reason == "cycle_2_lines:PrimalInfeasible"
    assert [metric["stage"] for metric in model.stage_metrics][-1] == "final_lines"
    np.testing.assert_array_equal(result.bestfit_lsf, np.full(8, 126.0))


def test_later_channel_fallback_is_kept_after_a_successful_cycle(workflow_model):
    config = LSFSurfaceIterativeConfig(n_refinement_cycles=2)
    fixed_knots = np.array([5000.0, 5003.5, 5007.0])
    first = _lsf_state(config, 2.0, knots=fixed_knots)
    fallback = _lsf_state(
        config,
        2.0,
        status="fallback",
        knots=fixed_knots,
    )
    fallback.metrics["B"]["reason"] = "previous_surface:PrimalInfeasible"
    fallback.metrics["B"]["knots_fixed"] = True
    model = workflow_model(
        cycles=2,
        solver_steps=[
            ("Solved", [1, 1, 1, 1, 1, 1]),
            ("Solved", [2, 3]),
            ("Solved", [4, 5, 6, 7]),
            ("Solved", [8, 9]),
            ("Solved", [10, 11, 12, 13]),
            ("Solved", [14, 15]),
            ("Solved", [16, 17, 18, 19]),
        ],
        lsf_states=[first, fallback],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    assert result.fit_status == "Solved"
    assert result.lsf_state is fallback
    assert result.lsf_state.completed_cycles == 2
    assert result.lsf_state.metrics["B"]["status"] == "fallback"
    assert result.lsf_state.metrics["B"]["knots_fixed"] is True
    np.testing.assert_array_equal(
        result.lsf_state.knot_vectors["B"],
        first.knot_vectors["B"],
    )


@pytest.mark.parametrize(
    ("final_steps", "expected_continuum", "expected_line", "failure"),
    [
        (
            [("PrimalInfeasible", [0, 0])],
            "PrimalInfeasible",
            "not_run",
            "final_continuum:PrimalInfeasible",
        ),
        (
            [("Solved", [20, 30]), ("PrimalInfeasible", [0, 0, 0, 0])],
            "Solved",
            "PrimalInfeasible",
            "final_lines:PrimalInfeasible",
        ),
    ],
)
def test_final_stage_failures_record_status_and_keep_last_complete_model(
    workflow_model,
    final_steps,
    expected_continuum,
    expected_line,
    failure,
):
    config = LSFSurfaceIterativeConfig(n_refinement_cycles=1)
    model = workflow_model(
        cycles=1,
        solver_steps=[
            ("Solved", [1, 1, 1, 1, 1, 1]),
            ("Solved", [2, 3]),
            ("Solved", [4, 5, 6, 7]),
            *final_steps,
        ],
        lsf_states=[_lsf_state(config, 2.0)],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    assert result.fit_status == f"failed:{failure}"
    assert result.lsf_state.final_continuum_status == expected_continuum
    assert result.lsf_state.final_line_status == expected_line
    assert result.lsf_state.failure_reason == failure
    np.testing.assert_array_equal(result.bestfit_lsf, np.full(8, 54.0))


def test_result_keeps_seed_bestfit_and_residual_but_exposes_final_lsf_model(
    workflow_model,
):
    config = LSFSurfaceIterativeConfig(n_refinement_cycles=1)
    model = workflow_model(
        cycles=1,
        solver_steps=[
            ("Solved", [1, 1, 1, 1, 1, 1]),
            ("Solved", [2, 3]),
            ("Solved", [4, 5, 6, 7]),
            ("Solved", [8, 9]),
            ("Solved", [10, 11, 12, 13]),
        ],
        lsf_states=[_lsf_state(config, 2.0)],
    )

    result = model.fit(np.zeros(8), np.ones(8))

    np.testing.assert_array_equal(result.bestfit, np.full(8, 123.0))
    np.testing.assert_array_equal(result.resid, np.full(8, -5.0))
    np.testing.assert_array_equal(result.bestfit_lsf, np.full(8, 126.0))
    assert not np.array_equal(result.bestfit, result.bestfit_lsf)
