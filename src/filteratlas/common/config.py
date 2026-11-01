"""Load the YAML configuration in `configs/`.

Every tunable number of the pipeline lives in one of these files, each with a comment on where the value comes
from. Code reads them through `cfg(name)`.

The repository root is located from this file; set FA_REPO to use a different checkout's configs.
"""
import os
import re
from functools import lru_cache
from pathlib import Path

import yaml

REPO = Path(os.environ.get("FA_REPO", Path(__file__).resolve().parents[3]))
CONFIGS = REPO / "configs"


class _Loader(yaml.SafeLoader):
    """SafeLoader that also reads 1e5, 0.988e6 and 2.249e-6 as floats (YAML 1.1 would read the first two as strings)."""


_Loader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(r"^[-+]?(?:[0-9][0-9_]*\.[0-9_]*(?:[eE][-+]?[0-9]+)?|[0-9][0-9_]*[eE][-+]?[0-9]+|\.[0-9_]+(?:[eE][-+]?[0-9]+)?"
               r"|[-+]?\.(?:inf|Inf|INF)|\.(?:nan|NaN|NAN))$"),
    list("-+0123456789."))


@lru_cache(maxsize=None)
def cfg(name):
    """Parsed `configs/<name>.yaml`, cached.

    Values are returned exactly as written in the file (lists stay lists); numbers such as 1e5 are read as floats.

    Parameters:
        name: file name without extension, e.g. "figures_of_merit".
    Returns:
        dict with the content of the file.
    Example:
        >>> from filteratlas.common.config import cfg
        >>> cfg("figures_of_merit")["band_A"]
    """
    return yaml.load((CONFIGS / f"{name}.yaml").read_text(encoding="utf-8"), Loader=_Loader)


def tup(x):
    """Convert YAML lists to tuples, recursively.

    Parameters:
        x: any value read from a configuration file.
    Returns:
        the same value with every list replaced by a tuple.
    Example:
        >>> from filteratlas.common.config import tup
        >>> tup([[0.5, 3.0], 1])
    """
    if isinstance(x, list):
        return tuple(tup(v) for v in x)
    return x
