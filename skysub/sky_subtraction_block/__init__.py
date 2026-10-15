"""Sky-only block LSF reconstruction, independent of Science in-situ fitting."""
from pathlib import Path
import sys

# The trained archive and numerical models retain top-level sky_decomp/mlp_predictor imports.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
