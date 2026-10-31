"""Atlas data fetch (network only, cached under D("mp"); reruns are offline):
1. MP summary structures for every proxy-bound survivor without MP phonons (Debye score, supercell size, SevenNet);
2. MP force constants for the MP-phonon survivors that pass the DOS rule and were not screened in step 02 (needed
   for the dynamical-stability screen on force constants).
"""
import sys
from pathlib import Path

import pandas as pd

from filteratlas.ingest import mp as mq  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402
from filteratlas.dynamics.stability_screen import fetch_fc  # noqa: E402


def d7a_todo():
    p = pd.read_parquet(F("stage3_proxy.parquet"))
    s = p[~p.excluded_by_proxy_bound & (p.phonon_id != "")]
    t = pd.read_parquet(F("stage1_triviality.parquet"), columns=["material_id", "status", "stable"])
    done = set(pd.read_parquet(F("stage1r_stability_full.parquet")).material_id)
    m = s.merge(t, on="material_id", how="left")
    return m[(m.status == "ok") & m.stable.fillna(False) & ~m.material_id.isin(done)]


def main(fc_shard=None):
    """fc_shard="i/n" (or --fc-shard i/n) fetches only the i-th of n interleaved parts of the force-constant list
    (network-bound; parallel shards; files already on disk are skipped)."""
    if fc_shard is None and "--fc-shard" in sys.argv:
        fc_shard = sys.argv[sys.argv.index("--fc-shard") + 1]
    if fc_shard is None:
        p = pd.read_parquet(F("stage3_proxy.parquet"))
        ids = sorted(p[~p.excluded_by_proxy_bound & (p.phonon_id == "")].material_id)
        got = mq.structures(ids)
        print(f"structures: {len(got)} / {len(ids)}", flush=True)
    si, sn = map(int, fc_shard.split("/")) if fc_shard else (0, 1)
    todo = d7a_todo()
    rows = todo[todo.phonon_method == "pheasy"].sort_values("material_id").to_dict("records")
    fetch_fc([r for k, r in enumerate(rows) if k % sn == si])
    print("force constants fetched", flush=True)


if __name__ == "__main__":
    main()
