from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from skysub import decompose_parallel
from skysub.sky_decomp.fit import SkyDecomp
from skysub.sky_decomp.lsf_surface_iterative import SkyDecompLSFSurfaceIterative


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
    monkeypatch.setitem(sys.modules, "sky_decomp.lsf_surface_iterative", fake_module)
    monkeypatch.setattr(decompose_parallel, "_clamp_native_threads", lambda *args: None)
    monkeypatch.setattr(decompose_parallel.fits, "open", lambda *args, **kwargs: fake_hdul)

    decompose_parallel.init_worker(
        np.arange(2.0),
        0.5,
        "/custom/base",
        1.0,
        "input.fits",
        fit_model="lsf-surface-iterative",
        palace_oh_suffix="_joint_v2_updated",
        palace_diffuse_suffix="_joint_native_adam_invsky_p2_10000iter_lvm",
    )

    assert captured["base_dir"] == "/custom/base"
    assert captured["palace_oh_suffix"] == "_joint_v2_updated"
    assert captured["palace_diffuse_suffix"] == (
        "_joint_native_adam_invsky_p2_10000iter_lvm"
    )


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
            "_joint_v2_updated",
            "--palace-diffuse-suffix",
            "_joint_native_adam_invsky_p2_10000iter_lvm",
        ],
    )

    decompose_parallel.main()

    assert captured["palace_dir"] == "/custom/base"
    assert captured["palace_suffix"] == "_adam_100_v1"
    assert captured["palace_oh_suffix"] == "_joint_v2_updated"
    assert captured["palace_diffuse_suffix"] == (
        "_joint_native_adam_invsky_p2_10000iter_lvm"
    )
