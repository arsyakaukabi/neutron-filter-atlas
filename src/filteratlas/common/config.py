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
    """Parsed `configs/<name>.yaml` (cached). Values are returned exactly as written; lists stay lists."""
    return yaml.load((CONFIGS / f"{name}.yaml").read_text(encoding="utf-8"), Loader=_Loader)


def tup(x):
    """YAML lists -> tuples, recursively (the original code used tuples; equality and hashing are unchanged)."""
    if isinstance(x, list):
        return tuple(tup(v) for v in x)
    return x
