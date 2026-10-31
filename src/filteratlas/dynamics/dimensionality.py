"""Structural dimensionality (Cheon et al. 2017, pymatgen) of every compound evaluated OK with MP phonons.

Used for a sensitivity analysis added after inspection of the results (0-D = molecular / noble-gas crystals).
Output: stage1_dimensionality.parquet
"""
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from filteratlas.common.paths import D, F
from filteratlas.common.constants import N_WORKERS  # noqa: E402


def _dim_one(args):
    mid, path = args
    import warnings
    from pymatgen.analysis.dimensionality import get_dimensionality_cheon
    from pymatgen.core import Structure
    from filteratlas.ingest import mp as mq
    warnings.filterwarnings("ignore")
    try:
        s = Structure.from_dict(mq.load_dos(path)["structure"])
        return mid, str(get_dimensionality_cheon(s))
    except Exception as e:  # noqa: BLE001
        return mid, f"error:{type(e).__name__}"


def main():
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    t = t[t.status == "ok"]
    jobs = [(r.material_id, D("mp") / "phonon_dos" / r.phonon_method / f"{r.phonon_id}.json.gz")
            for r in t.itertuples()]
    with Pool(N_WORKERS) as p:
        rows = p.map(_dim_one, jobs, chunksize=8)
    d = pd.DataFrame(rows, columns=["material_id", "dim"])
    d.to_parquet(F("stage1_dimensionality.parquet"))
    print(d.dim.value_counts().to_dict())


if __name__ == "__main__":
    main()
