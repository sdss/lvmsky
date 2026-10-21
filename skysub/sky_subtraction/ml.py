"""The single import point of the ML predictor package.

The trained archive pickles objects of the top-level ``mlp_predictor`` package, and its
modules import ``sky_decomp`` as a top-level package as well; both need ``skysub/`` on
sys.path. ``mlp_predictor/__init__`` loads the whole package, torch included, so
reconstruction loads it too although it never runs inference.
"""

import sys
from pathlib import Path

_SKYSUB = str(Path(__file__).resolve().parents[1])
if _SKYSUB not in sys.path:
    sys.path.insert(0, _SKYSUB)

from mlp_predictor import data, inference, serialization  # noqa: E402
from mlp_predictor.sky_arm_correction import (  # noqa: E402
    SkyArm,
    line_protection_masks,
    protection_velocities,
    sky_arm_residual_correction,
)
from mlp_predictor.sky_line_scaling import (  # noqa: E402
    _E_SCALE,
    NEBULAR_LINES_RED,
    UNMEASURED_VELOCITY_KM_S,
    line_templates,
    oh_group_labels,
    sky_line_scaling_correction,
)

__all__ = [
    'data',
    'inference',
    'serialization',
    'SkyArm',
    'line_protection_masks',
    'protection_velocities',
    'sky_arm_residual_correction',
    'E_SCALE',
    'NEBULAR_LINES_RED',
    'UNMEASURED_VELOCITY_KM_S',
    'line_templates',
    'oh_group_labels',
    'sky_line_scaling_correction',
]
E_SCALE = _E_SCALE
