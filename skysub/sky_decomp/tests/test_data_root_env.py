"""LVMSKY_DATA_ROOT redirects the data root used by sky_decomp and mlp_predictor.

The variable is read when moon_zodi_model is imported, so every check runs in a
fresh interpreter with the environment it needs.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

SKYSUB = Path(__file__).resolve().parents[2]
REPO = SKYSUB.parent
BUNDLED = (SKYSUB / "sky_decomp" / "data").resolve()

PROBE_SKY_DECOMP = """
import json
from sky_decomp import moon_zodi_model, pixel_weights
print(json.dumps({
    "root": str(moon_zodi_model.DEFAULT_DATA_ROOT),
    "dir": str(moon_zodi_model.DEFAULT_DATA_DIR),
    "sens": str(pixel_weights.sensitivity_dir()),
    "sens_arg": str(pixel_weights.sensitivity_dir("/explicit")),
}))
"""

PROBE_MLP = """
import json
from mlp_predictor import data
print(json.dumps({"recon": str(data._infer_base_dir_for_reconstruction())}))
"""


def _probe(code, root=None):
    env = {k: v for k, v in os.environ.items() if k != "LVMSKY_DATA_ROOT"}
    if root is not None:
        env["LVMSKY_DATA_ROOT"] = str(root)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO), str(SKYSUB), env.get("PYTHONPATH", "")])
    out = subprocess.run([sys.executable, "-c", code], env=env, cwd=str(REPO),
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])


@pytest.fixture
def relocated_root(tmp_path):
    """A complete data root outside the code tree: a real directory whose
    entries link to the bundled data (so it resolves to itself, not the bundle)."""
    root = tmp_path / "lvmsky_data"
    root.mkdir()
    for entry in BUNDLED.iterdir():
        (root / entry.name).symlink_to(entry, target_is_directory=entry.is_dir())
    return root.resolve()


def test_unset_uses_bundled_data():
    paths = _probe(PROBE_SKY_DECOMP)
    assert Path(paths["root"]) == BUNDLED
    assert Path(paths["dir"]) == BUNDLED / "moon_zodi"
    assert Path(paths["sens"]) == BUNDLED / "sensitivity"


def test_env_redirects_sky_decomp(relocated_root):
    paths = _probe(PROBE_SKY_DECOMP, relocated_root)
    assert Path(paths["root"]) == relocated_root
    assert Path(paths["dir"]) == relocated_root / "moon_zodi"
    assert Path(paths["sens"]) == relocated_root / "sensitivity"


def test_explicit_argument_still_wins(relocated_root):
    paths = _probe(PROBE_SKY_DECOMP, relocated_root)
    assert paths["sens_arg"] == "/explicit"


def test_env_redirects_mlp_predictor_reconstruction(relocated_root):
    pytest.importorskip("torch")  # mlp_predictor imports torch
    paths = _probe(PROBE_MLP, relocated_root)
    assert Path(paths["recon"]) == relocated_root
