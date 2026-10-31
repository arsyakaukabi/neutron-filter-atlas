"""Shared fixtures. Tests that need data skip with a clear message when data/ is empty (fresh clone)."""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from filteratlas.common.config import cfg  # noqa: E402
from filteratlas.common.paths import DATA  # noqa: E402

FILES = cfg("paths")["files"]


def _have(*names):
    return all((DATA / FILES[n]).exists() for n in names)


def need(*names):
    """Skip unless the registered data files exist (data/README.md explains where they come from)."""
    missing = [n for n in names if not (DATA / FILES[n]).exists()]
    if missing:
        pytest.skip(f"needs data files {missing} (place the distributed data under data/; see data/README.md)")


def need_dir(key):
    p = DATA / cfg("paths")["dirs"][key]
    if not p.exists() or not any(p.iterdir()):
        pytest.skip(f"needs data directory {p} (see data/README.md)")


@pytest.fixture
def derived_json():
    def load(name):
        need(name)
        return json.loads((DATA / FILES[name]).read_text())
    return load


@pytest.fixture
def atlas():
    need("atlas.parquet")
    import pandas as pd
    return pd.read_parquet(DATA / FILES["atlas.parquet"])
