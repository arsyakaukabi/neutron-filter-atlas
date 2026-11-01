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
    """Path of a data directory registered in `configs/paths.yaml` (created if absent).

    Parameters:
        key: directory key under `dirs` in `configs/paths.yaml`, e.g. "mp" or "nucdata".
    Returns:
        pathlib.Path under the data root (`data/`, or `FA_DATA` if set).
    Example:
        >>> from filteratlas.common.paths import D
        >>> D("mp")
    """
    p = DATA / cfg("paths")["dirs"][key]
    p.mkdir(parents=True, exist_ok=True)
    return p


def F(name):
    """Path of a data file registered in `configs/paths.yaml` (its parent directory is created if absent).

    Parameters:
        name: file name under `files` in `configs/paths.yaml`, e.g. "atlas.parquet".
    Returns:
        pathlib.Path under the data root (`data/`, or `FA_DATA` if set).
    Raises:
        KeyError if the name is not registered.
    Example:
        >>> from filteratlas.common.paths import F
        >>> F("atlas.parquet")
    """
    files = cfg("paths")["files"]
    if name not in files:
        raise KeyError(f"data file '{name}' is not registered in configs/paths.yaml")
    p = DATA / files[name]
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def env_file():
    """Path of the `.env` file with API keys (`FA_ENV` overrides; the file is never tracked in git).

    Returns:
        pathlib.Path, by default `.env` in the repository root.
    Example:
        >>> from filteratlas.common.paths import env_file
        >>> env_file().name
    """
    return Path(os.environ.get("FA_ENV", REPO / ".env"))
