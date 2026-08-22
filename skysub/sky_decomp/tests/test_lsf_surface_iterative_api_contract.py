"""Stable API contract for the iterative LSF implementation facade."""

from dataclasses import MISSING, fields
import inspect
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from skysub import decompose_parallel
import skysub.sky_decomp.lsf_surface_iterative as iterative


def _signature_without_annotations(callable_object):
    signature = inspect.signature(callable_object)
    parameters = [
        parameter.replace(annotation=inspect.Parameter.empty)
        for parameter in signature.parameters.values()
    ]
    return str(
        signature.replace(
            parameters=parameters,
            return_annotation=inspect.Signature.empty,
        )
    )


def test_public_exports_and_pickle_qualified_names():
    assert iterative.__all__ == [
        "LSFSurfaceIterativeConfig",
        "LSFSurfaceIterativeResult",
        "LSFSurfaceState",
        "SkyDecompLSFSurfaceIterative",
        "apply_lsf_channelwise",
        "apply_lsf_surface",
        "build_bspline_basis",
        "build_lsf_operator",
        "build_skyline_mask",
        "continuum_fit_weights",
        "evaluate_bspline_basis",
        "evaluate_lsf_kernel",
        "evaluate_lsf_surface",
        "fit_bspline_channel",
        "fit_lsf_surface",
        "kernel_moments",
        "load_lsf_surface_state",
    ]

    for class_name in (
        "LSFSurfaceIterativeConfig",
        "LSFSurfaceState",
        "LSFSurfaceIterativeResult",
        "SkyDecompLSFSurfaceIterative",
    ):
        class_object = getattr(iterative, class_name)
        assert class_object.__module__ == "skysub.sky_decomp.lsf_surface_iterative"
        assert class_object.__qualname__ == class_name


def test_dataclass_field_order_and_defaults():
    config_fields = fields(iterative.LSFSurfaceIterativeConfig)
    assert tuple(field.name for field in config_fields) == (
        "line_weight",
        "skyline_cumulative_fraction",
        "skyline_half_width_angstrom",
        "huber_transition_sigma",
        "n_refinement_cycles",
        "n_basis",
        "degree",
        "roughness_fraction",
        "fallback_prior_fraction",
        "information_prior_max_boost",
        "background_degree",
        "blue_fit_lower",
    )
    assert tuple(field.default for field in config_fields) == (
        5.0e-4,
        0.99,
        2.0,
        3.0,
        5,
        6,
        3,
        0.01,
        1.0e-4,
        50.0,
        3,
        5500.0,
    )

    state_fields = fields(iterative.LSFSurfaceState)
    assert tuple(field.name for field in state_fields) == (
        "coefficients",
        "knot_vectors",
        "degrees",
        "channel_bounds",
        "tap_offsets",
        "config",
        "metrics",
        "requested_cycles",
        "completed_cycles",
        "wave_n",
        "wave_min",
        "wave_max",
        "wave_sha256",
        "fit_status",
        "failure_reason",
        "final_continuum_status",
        "final_line_status",
        "knot_strategy",
        "legacy_kernel_representation",
        "schema_version",
    )
    assert all(field.default is MISSING for field in state_fields[:13])
    assert tuple(field.default for field in state_fields[13:]) == (
        "",
        "",
        "not_run",
        "not_run",
        "first_cycle_adaptive_then_fixed",
        "native_grid_channel_median",
        1,
    )

    assert tuple(field.name for field in fields(iterative.LSFSurfaceIterativeResult)) == (
        "coef",
        "coef_err",
        "bestfit",
        "resid",
        "resid_level",
        "fit_status",
        "fit_summary",
        "reduced_chi2",
        "fit_elapsed_sec",
        "components",
        "design_names",
        "t_o2",
        "t_o2_err",
        "r2",
        "rms_resid",
        "peak_memory_mb",
        "o2_fit_status",
        "o2_fit_summary",
        "o2_fit_elapsed_sec",
        "o2_valid_frac",
        "lsf_kernels",
        "lsf_metrics",
        "bestfit_lsf",
        "moon_knots",
        "moon_boosted_pixels",
        "vector_o2",
        "o2_prefit_amp",
        "bestfit_lsf_sigma",
        "zodi_names",
        "zodi_knots",
        "lsf_state",
    )


def test_public_callable_signatures():
    expected = {
        "apply_lsf_channelwise": ("(state, wave, values, *, allow_different_grid=False)"),
        "apply_lsf_surface": (
            "(wave, values, kernel_surface, *, channel_bounds=None, tap_offsets=None)"
        ),
        "build_bspline_basis": "(wave, n_basis, degree, information)",
        "build_skyline_mask": (
            "(wave, line_stick_matrices, *, cumulative_fraction=0.99, half_width_angstrom=2.0)"
        ),
        "continuum_fit_weights": (
            "(wave, residual, ivar, skyline_mask, *, line_weight=0.0005, "
            "huber_transition_sigma=3.0)"
        ),
        "evaluate_bspline_basis": "(wave, knots, degree)",
        "evaluate_lsf_kernel": "(state, wavelength)",
        "evaluate_lsf_surface": "(state, wave)",
        "fit_bspline_channel": (
            "(wave, source, target, ivar, fallback_kernel, *, n_basis=6, "
            "degree=3, roughness_fraction=0.01, fallback_prior_fraction=0.0001, "
            "information_prior_max_boost=50.0, background_degree=3, "
            "free_amplitude=False, knot_vector=None)"
        ),
        "fit_lsf_surface": (
            "(wave, flux, ivar, source, fixed_background, fallback_kernels, "
            "config, *, previous_state=None)"
        ),
        "kernel_moments": "(kernel_surface)",
        "load_lsf_surface_state": "(filename, spectrum_index=0)",
    }
    for name, signature in expected.items():
        assert _signature_without_annotations(getattr(iterative, name)) == signature

    assert _signature_without_annotations(iterative.SkyDecompLSFSurfaceIterative.__init__) == (
        "(self, *args, config=None, **kwargs)"
    )
    assert _signature_without_annotations(iterative.SkyDecompLSFSurfaceIterative.fit) == (
        "(self, flux, ivar, *, verbose=False, moon_amp_prior=None, "
        "zodi_amp_prior=None, moon_amp_prior_lambda=0.0, "
        "zodi_amp_prior_lambda=0.0)"
    )


def test_batch_and_cli_defaults_resolve_to_five_cycles(monkeypatch):
    assert (
        inspect.signature(decompose_parallel.init_worker).parameters["n_refinement_cycles"].default
        == 5
    )
    assert inspect.signature(decompose_parallel.run).parameters["n_refinement_cycles"].default == 5

    captured = {}
    monkeypatch.setattr(
        decompose_parallel,
        "run",
        lambda **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        decompose_parallel,
        "extract_meta_and_coef_products",
        lambda **kwargs: None,
    )
    monkeypatch.setattr(
        decompose_parallel,
        "thin_fits_every_n",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(sys, "argv", ["decompose_parallel.py", "input.fits", "palace"])

    decompose_parallel.main()

    assert captured["n_refinement_cycles"] == 5


def test_moon_zodi_cli_does_not_require_legacy_palace_path(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        decompose_parallel,
        "run",
        lambda **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        decompose_parallel,
        "extract_meta_and_coef_products",
        lambda **kwargs: None,
    )
    monkeypatch.setattr(
        decompose_parallel,
        "thin_fits_every_n",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "decompose_parallel.py",
            "input.fits",
            "--fit-model",
            decompose_parallel.MOON_ZODI_FIT_MODEL,
            "--moon-zodi-data-root",
            "/external/data",
        ],
    )

    decompose_parallel.main()

    assert captured["palace_dir"] is None
    assert str(captured["moon_zodi_data_root"]) == "/external/data"


def test_legacy_cli_still_requires_palace_path(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        ["decompose_parallel.py", "input.fits", "--fit-model", "baseline"],
    )

    with pytest.raises(SystemExit) as error:
        decompose_parallel.main()

    assert error.value.code == 2


def test_chunk_worker_preserves_indices_and_reports_each_row(monkeypatch):
    progress = []

    class Decomposer:
        def fit(self, flux, ivar, **kwargs):
            return float(flux[0])

    monkeypatch.setattr(decompose_parallel, "_WORKER_DECOMPOSER", Decomposer())
    monkeypatch.setattr(decompose_parallel, "_WORKER_FACTOR", 2.0)
    monkeypatch.setattr(decompose_parallel, "_WORKER_FIT_MODEL", "baseline")
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_FLUX",
        {"sci": np.arange(12, dtype=float).reshape(4, 3)},
    )
    monkeypatch.setattr(
        decompose_parallel,
        "_WORKER_PROGRESS_QUEUE",
        SimpleNamespace(put=progress.append),
    )

    kind, results = decompose_parallel.fit_chunk_worker(("sci", 1, 3))

    assert kind == "sci"
    assert results == [(1, 6.0), (2, 12.0)]
    assert progress == [1, 1]


def test_batch_run_restores_row_order_and_writes_three_outputs(monkeypatch, tmp_path):
    class Future:
        def __init__(self, task):
            self.task = task

        def result(self):
            kind, first, last = self.task
            return kind, [(index, f"{kind}-{index}") for index in range(first, last)]

    class Executor:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, worker, task):
            return Future(task)

    class ProgressQueue:
        def get_nowait(self):
            raise decompose_parallel.queue_mod.Empty

        def close(self):
            pass

        def join_thread(self):
            pass

    class Counter:
        def __init__(self, *_args):
            self.value = 0

    class ProgressBar:
        def set_postfix(self, **kwargs):
            pass

        def refresh(self):
            pass

        def update(self, value):
            pass

        def close(self):
            pass

    class FitsFile:
        def __enter__(self):
            return {"FLUX_SCI": SimpleNamespace(data=np.empty((3, 4)))}

        def __exit__(self, *args):
            pass

    written = []
    monkeypatch.setattr(decompose_parallel, "resolve_base_dir", lambda path: tmp_path)
    monkeypatch.setattr(decompose_parallel.fits, "getdata", lambda *args: np.arange(4.0))
    monkeypatch.setattr(decompose_parallel.fits, "open", lambda *args, **kwargs: FitsFile())
    monkeypatch.setattr(decompose_parallel.mp, "get_context", lambda *args: SimpleNamespace(
        Queue=ProgressQueue,
        Value=Counter,
    ))
    monkeypatch.setattr(decompose_parallel, "ProcessPoolExecutor", Executor)
    monkeypatch.setattr(
        decompose_parallel,
        "wait",
        lambda pending, **kwargs: (set(pending), set()),
    )
    monkeypatch.setattr(decompose_parallel, "tqdm", lambda **kwargs: ProgressBar())
    monkeypatch.setattr(
        decompose_parallel,
        "results_to_fits",
        lambda results, path: written.append((results, path.name)),
    )

    decompose_parallel.run(
        "input.fits",
        tmp_path,
        n_workers=2,
        lsf_sigma=0.5,
        factor=1.0,
        output_dir=tmp_path,
        chunk_size=2,
        max_in_flight=2,
    )

    assert written == [
        (["sci-0", "sci-1", "sci-2"], "input_decomp_sci.fits"),
        (["sky1-0", "sky1-1", "sky1-2"], "input_decomp_sky1.fits"),
        (["sky2-0", "sky2-1", "sky2-2"], "input_decomp_sky2.fits"),
    ]
