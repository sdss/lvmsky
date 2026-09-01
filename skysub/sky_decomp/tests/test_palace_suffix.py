import importlib.util
import inspect
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from skysub import decompose_parallel
from skysub.sky_decomp.fit import SkyDecomp
from skysub.sky_decomp.lsf_surface_iterative import SkyDecompLSFSurfaceIterative
from skysub.sky_decomp.moon_zodi_model import (
    DEFAULT_PALACE_DIFFUSE_SUFFIX,
    DEFAULT_PALACE_OH_SUFFIX,
)


def test_decompose_parallel_direct_script_imports_package_from_any_cwd(tmp_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    completed = subprocess.run(
        [sys.executable, str(Path(decompose_parallel.__file__).resolve()), "--help"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert completed.returncode == 0, completed.stderr
    assert "LVM sky spectral decomposition" in completed.stdout


def _path_only_model(tmp_path, palace_suffix):
    model = object.__new__(SkyDecomp)
    model.pmd_dir = tmp_path / "palace" / "PMD"
    model.pmd_dir.mkdir(parents=True)
    model.palace_suffix = SkyDecomp._validate_palace_suffix(palace_suffix)
    model.palace_oh_suffix = None
    model.palace_diffuse_suffix = None
    return model


def test_default_palace_paths_are_unchanged(tmp_path):
    model = _path_only_model(tmp_path, None)
    expected = model.pmd_dir / "pmd_refcont.dat"
    expected.touch()

    assert model._pmd_path("pmd_refcont.dat") == expected


def test_suffix_selects_only_versioned_oh_and_diffuse_tables(tmp_path):
    model = _path_only_model(tmp_path, "_adam_100_v1")
    expected = {
        "pmd_refcont.dat": model.pmd_dir / "pmd_refcont_adam_100_v1.dat",
        "pmd_popmodel_OH.dat": model.pmd_dir / "pmd_popmodel_OH_adam_100_v1.dat",
        "pmd_popmodel_O2.dat": model.pmd_dir / "pmd_popmodel_O2.dat",
        "pmd_intdata_atom.dat": model.pmd_dir / "pmd_intdata_atom.dat",
        "pmd_intmodel_Orc.dat": model.pmd_dir / "pmd_intmodel_Orc.dat",
    }
    for path in expected.values():
        path.touch()

    assert model.palace_suffix == "_adam_100_v1"
    for canonical_name, selected_path in expected.items():
        assert model._pmd_path(canonical_name) == selected_path


def test_suffix_is_appended_exactly_without_adding_a_separator(tmp_path):
    model = _path_only_model(tmp_path, "revision1")
    expected = model.pmd_dir / "pmd_refcontrevision1.dat"
    expected.touch()

    assert model.palace_suffix == "revision1"
    assert model._pmd_path("pmd_refcont.dat") == expected


def test_table_specific_suffixes_override_the_common_suffix(tmp_path):
    model = _path_only_model(tmp_path, "_common")
    model.palace_oh_suffix = "_oh_v2"
    model.palace_diffuse_suffix = "_diffuse_adam"
    oh_path = model.pmd_dir / "pmd_popmodel_OH_oh_v2.dat"
    diffuse_path = model.pmd_dir / "pmd_refcont_diffuse_adam.dat"
    oh_path.touch()
    diffuse_path.touch()

    assert model._pmd_path("pmd_popmodel_OH.dat") == oh_path
    assert model._pmd_path("pmd_refcont.dat") == diffuse_path


def test_native_diffuse_table_is_loaded_row_for_row_without_interpolation(
    tmp_path, monkeypatch
):
    wave = np.array([3600.0, 3600.5, 3601.0])
    diffuse = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    )
    table = {
        "lam": wave / 1.0e4,
        "trans": np.ones(3),
        "fH2O": np.ones(3),
        "fcHO2": diffuse[0],
        "fcFeO": diffuse[1],
        "fcO2Ac": diffuse[2],
        "isX": np.zeros(3, dtype=int),
    }
    path = tmp_path / "pmd_refcont_native.dat"
    from astropy.table import Table

    Table(table).write(path, format="ascii.basic")
    model = _path_only_model(tmp_path, None)
    model.wave = wave
    monkeypatch.setattr(model, "_pmd_path", lambda name: path)
    monkeypatch.setattr(
        np,
        "interp",
        lambda *args, **kwargs: pytest.fail("np.interp must not be called"),
    )

    matrix, names = model._build_diffuse()

    assert names == ["HO2", "FeO", "O2Ac"]
    assert np.array_equal(matrix, diffuse)


@pytest.mark.parametrize("suffix", ["", "_", "../v1", "folder/v1", "_v1.dat"])
def test_invalid_palace_suffix_is_rejected(suffix):
    with pytest.raises(ValueError, match="palace_suffix"):
        SkyDecomp._validate_palace_suffix(suffix)


def test_lsf_surface_constructor_propagates_suffix_to_skydecomp(monkeypatch):
    captured = {}

    def fake_init(self, *args, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(SkyDecomp, "__init__", fake_init)

    SkyDecompLSFSurfaceIterative(
        np.arange(4.0),
        palace_suffix="_common",
        palace_oh_suffix="_oh_v2",
        palace_diffuse_suffix="_diffuse_adam",
    )

    assert captured["palace_suffix"] == "_common"
    assert captured["palace_oh_suffix"] == "_oh_v2"
    assert captured["palace_diffuse_suffix"] == "_diffuse_adam"


def test_parallel_worker_propagates_suffix_to_lsf_surface_model(monkeypatch):
    captured = {}

    class FakeDecomposer:
        def __init__(self, *args, **kwargs):
            captured.update(kwargs)

    fake_module = SimpleNamespace(
        LSFSurfaceIterativeConfig=lambda **kwargs: kwargs,
        SkyDecompLSFSurfaceIterative=FakeDecomposer,
    )
    fake_hdul = {
        "FLUX_SCI": SimpleNamespace(data=np.empty((1, 2))),
        "FLUX_SKY_NEAR": SimpleNamespace(data=np.empty((1, 2))),
        "FLUX_SKY_FAR": SimpleNamespace(data=np.empty((1, 2))),
    }
    monkeypatch.setitem(
        sys.modules,
        "skysub.sky_decomp.lsf_surface_iterative",
        fake_module,
    )
    monkeypatch.setattr(decompose_parallel, "_clamp_native_threads", lambda *args: None)
    monkeypatch.setattr(decompose_parallel.fits, "open", lambda *args, **kwargs: fake_hdul)

    decompose_parallel.init_worker(
        np.arange(2.0),
        0.5,
        "/custom/base",
        1.0,
        "input.fits",
        fit_model="lsf-surface-iterative",
        palace_oh_suffix="_custom_oh",
        palace_diffuse_suffix="_custom_diffuse",
    )

    assert captured["base_dir"] == "/custom/base"
    assert captured["palace_oh_suffix"] == "_custom_oh"
    assert captured["palace_diffuse_suffix"] == "_custom_diffuse"


def test_bundled_skydecomp_owns_default_oh_and_diffuse_tables(monkeypatch, tmp_path):
    (tmp_path / "bundle_manifest.json").touch()
    pmd_dir = tmp_path / "palace" / "PMD"
    pmd_dir.mkdir(parents=True)
    oh_path = pmd_dir / f"pmd_popmodel_OH{DEFAULT_PALACE_OH_SUFFIX}.dat"
    diffuse_path = pmd_dir / f"pmd_refcont{DEFAULT_PALACE_DIFFUSE_SUFFIX}.dat"
    oh_path.touch()
    diffuse_path.touch()
    monkeypatch.setattr(SkyDecomp, "_build_static_basis", lambda self: None)

    model = SkyDecomp(np.arange(2.0), base_dir=tmp_path)

    assert model.palace_suffix is None
    assert model.palace_oh_suffix == DEFAULT_PALACE_OH_SUFFIX
    assert model.palace_diffuse_suffix == DEFAULT_PALACE_DIFFUSE_SUFFIX
    assert model._pmd_path("pmd_popmodel_OH.dat") == oh_path
    assert model._pmd_path("pmd_refcont.dat") == diffuse_path


def test_explicit_palace_suffixes_override_bundled_defaults(monkeypatch, tmp_path):
    (tmp_path / "bundle_manifest.json").touch()
    monkeypatch.setattr(SkyDecomp, "_build_static_basis", lambda self: None)

    table_specific = SkyDecomp(
        np.arange(2.0),
        base_dir=tmp_path,
        palace_oh_suffix="_custom_oh",
        palace_diffuse_suffix="_custom_diffuse",
    )
    common = SkyDecomp(
        np.arange(2.0),
        base_dir=tmp_path,
        palace_suffix="_custom_common",
    )

    assert table_specific.palace_oh_suffix == "_custom_oh"
    assert table_specific.palace_diffuse_suffix == "_custom_diffuse"
    assert common.palace_suffix == "_custom_common"
    assert common.palace_oh_suffix is None
    assert common.palace_diffuse_suffix is None


def test_nonbundled_skydecomp_keeps_legacy_unsuffixed_default(monkeypatch, tmp_path):
    monkeypatch.setattr(SkyDecomp, "_build_static_basis", lambda self: None)

    model = SkyDecomp(np.arange(2.0), base_dir=tmp_path)

    assert model.palace_suffix is None
    assert model.palace_oh_suffix is None
    assert model.palace_diffuse_suffix is None


def test_split_zodi_without_legacy_path_uses_packaged_bundle(monkeypatch, tmp_path):
    bundle_root = tmp_path / "bundle"
    validated = []
    monkeypatch.setattr(
        decompose_parallel,
        "DEFAULT_MOON_ZODI_DATA_ROOT",
        bundle_root,
    )
    monkeypatch.setattr(
        decompose_parallel,
        "validate_decomposition_data_root",
        validated.append,
    )
    monkeypatch.setattr(
        decompose_parallel,
        "resolve_base_dir",
        lambda *_args: pytest.fail("legacy PALACE resolution must not run"),
    )

    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        decompose_parallel.SPLIT_ZODI_FIT_MODEL,
    )

    assert base_dir == bundle_root.resolve()
    assert data_root == bundle_root.resolve()
    assert validated == [str(bundle_root.resolve())]


def test_split_zodi_positional_bundle_uses_bundle_contract(monkeypatch, tmp_path):
    bundle_root = tmp_path / "bundle"
    bundle_root.mkdir()
    (bundle_root / "bundle_manifest.json").touch()
    validated = []
    monkeypatch.setattr(
        decompose_parallel,
        "validate_decomposition_data_root",
        validated.append,
    )
    monkeypatch.setattr(
        decompose_parallel,
        "resolve_base_dir",
        lambda *_args: pytest.fail("legacy PALACE resolution must not run"),
    )

    base_dir, data_root = decompose_parallel.resolve_runtime_data_roots(
        decompose_parallel.SPLIT_ZODI_FIT_MODEL,
        palace_dir=bundle_root,
    )
    assert base_dir == bundle_root.resolve()
    assert data_root == bundle_root.resolve()
    assert validated == [str(bundle_root.resolve())]


def test_cli_forwards_optional_palace_suffix(monkeypatch):
    captured = {}
    monkeypatch.setattr(decompose_parallel, "run", lambda **kwargs: captured.update(kwargs))
    monkeypatch.setattr(decompose_parallel, "extract_meta_and_coef_products", lambda **kwargs: None)
    monkeypatch.setattr(decompose_parallel, "thin_fits_every_n", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "decompose_parallel.py",
            "input.fits",
            "/custom/base",
            "--fit-model",
            "lsf-surface-iterative",
            "--palace-suffix",
            "_adam_100_v1",
            "--palace-oh-suffix",
            "_custom_oh",
            "--palace-diffuse-suffix",
            "_custom_diffuse",
        ],
    )

    decompose_parallel.main()

    assert captured["palace_dir"] == "/custom/base"
    assert captured["palace_suffix"] == "_adam_100_v1"
    assert captured["palace_oh_suffix"] == "_custom_oh"
    assert captured["palace_diffuse_suffix"] == "_custom_diffuse"


def test_split_zodi_cli_does_not_require_legacy_palace_path(monkeypatch):
    captured = {}
    monkeypatch.setattr(decompose_parallel, "run", lambda **kwargs: captured.update(kwargs))
    monkeypatch.setattr(decompose_parallel, "extract_meta_and_coef_products", lambda **kwargs: None)
    monkeypatch.setattr(decompose_parallel, "thin_fits_every_n", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "decompose_parallel.py",
            "input.fits",
            "--fit-model",
            decompose_parallel.SPLIT_ZODI_FIT_MODEL,
        ],
    )

    decompose_parallel.main()

    assert captured["palace_dir"] is None
    assert captured["moon_zodi_data_root"] is None
    assert captured["palace_suffix"] is None
    assert captured["palace_oh_suffix"] is None
    assert captured["palace_diffuse_suffix"] is None


def test_mlp_gaussian_reconstruction_forwards_skydecomp_selection(
    monkeypatch, tmp_path
):
    skysub_root = Path(__file__).resolve().parents[2]
    monkeypatch.syspath_prepend(str(skysub_root))
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "matplotlib"))
    module_path = skysub_root / "mlp_predictor" / "data.py"
    spec = importlib.util.spec_from_file_location("mlp_data_for_test", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    signature = inspect.signature(module.reconstruct_with_lsf)
    assert signature.parameters["palace_oh_suffix"].default is None
    assert signature.parameters["palace_diffuse_suffix"].default is None

    captured = {}
    expected = {"total": np.arange(2.0)}
    monkeypatch.setattr(
        module,
        "reconstruct_component_spectra",
        lambda **kwargs: captured.update(kwargs) or expected,
    )

    result = module.reconstruct_with_lsf(
        np.arange(2.0),
        np.ones(3),
        np.ones(2),
        base_dir=tmp_path,
        split_zodi=True,
        n_zodi_spline_knots=7,
        palace_oh_suffix="_custom_oh",
        palace_diffuse_suffix="_custom_diffuse",
    )

    assert result is expected
    assert captured["split_zodi"] is True
    assert captured["n_zodi_spline_knots"] == 7
    assert captured["palace_oh_suffix"] == "_custom_oh"
    assert captured["palace_diffuse_suffix"] == "_custom_diffuse"
