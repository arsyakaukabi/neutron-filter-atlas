"""Make `src/` importable when a pipeline is run as `python pipelines/NN_name.py` from the repository root."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
