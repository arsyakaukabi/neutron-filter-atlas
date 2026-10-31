"""Every data path used by the code, resolved from `configs/paths.yaml`.

- `D(key)`  -> a data directory (e.g. D("mp") = data/upstream/mp).
- `F(name)` -> a data file by its registered name (e.g. F("atlas.parquet") = data/derived/atlas.parquet).

File names are kept from earlier versions of the code; only the directory layout is defined here.
The data root is `data/` in the repository; set FA_DATA to point at another data tree. Parent directories are
created on demand.
"""
import os
from pathlib import Path

from .config import REPO, cfg

DATA = Path(os.environ.get("FA_DATA", REPO / cfg("paths")["data_root"]))


def D(key):
    p = DATA / cfg("paths")["dirs"][key]
    p.mkdir(parents=True, exist_ok=True)
    return p


def F(name):
    files = cfg("paths")["files"]
    if name not in files:
        raise KeyError(f"data file '{name}' is not registered in configs/paths.yaml")
    p = DATA / files[name]
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def env_file():
    """Path of the .env file with API keys (never inside the repository's tracked files)."""
    return Path(os.environ.get("FA_ENV", REPO / ".env"))
