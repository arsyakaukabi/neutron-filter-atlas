"""Fix the MLIP benchmark set before any MLIP run (called by step 03, MLIP dynamics).

150 usable compounds with pheasy MP phonons that are dynamically clean (commensurate-q minimum >= -0.5 THz and not
flagged anywhere on the 40 A mesh, so the DFT reference itself is clean), stratified by cell size
(60 x <=10 sites, 60 x 11-20, 30 x 21-40), seed 20261024; plus the EXFOR benchmark crystals.
-> the MLIP benchmark-set JSON in data/derived
"""
import json
import sys
from pathlib import Path

import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402

from filteratlas.common.config import cfg  # noqa: E402

_MS = cfg("mlip")["benchmark_set"]
SEED = _MS["seed"]
STRATA = [((a, b), n) for a, b, n in _MS["strata_nsites_n"]]
# EXFOR benchmark crystals (phases fixed before the computation): MP id with phonons where available, else the summary entry
EXFOR_CRYSTALS = {f: v[0] for f, v in cfg("benchmark")["crystals"].items() if f != "ZrO2"}   # same phases as the benchmark


def main():
    u = pd.read_parquet(F("stage1_usable.parquet"))
    sf = pd.read_parquet(F("stage1r_stability_full.parquet"))
    m = u.merge(sf[["material_id", "status", "unstable_05", "unstable_comm_05"]].rename(columns={"status": "fc_status"}),
                on="material_id")
    pool = m[(m.phonon_method == "pheasy") & (m.fc_status == "ok") & ~m.unstable_05.astype(bool)
             & ~m.unstable_comm_05.astype(bool)]
    sel = []
    for i, ((lo, hi), n) in enumerate(STRATA):
        g = pool[(pool.nsites >= lo) & (pool.nsites <= hi)]
        take = g.sort_values("material_id").sample(min(n, len(g)), random_state=SEED + i)
        sel += [{"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol, "nsites": int(r.nsites),
                 "phonon_id": r.phonon_id, "phonon_method": r.phonon_method, "arm": f"random_{lo}-{hi}",
                 "fc_npz": f"cache/mp/force_constants/{r.phonon_id}.npz"} for r in take.itertuples()]
    s = pd.read_parquet(D("mp") / "summary_exp_nearhull.parquet")
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    for f, sg in EXFOR_CRYSTALS.items():
        e = s[(s.formula == f) & (s.sg_symbol == sg)].sort_values("e_hull")
        ph = t[t.material_id.isin(e.material_id) & (t.status == "ok") & (t.phonon_method == "pheasy")]
        if len(ph):
            r = ph.merge(e[["material_id", "e_hull"]], on="material_id").sort_values("e_hull").iloc[0]
            npz = D("mp") / "force_constants" / f"{r.phonon_id}.npz"
            sel.append({"material_id": r.material_id, "formula": f, "sg_symbol": sg, "nsites": int(r.nsites),
                        "phonon_id": r.phonon_id, "phonon_method": "pheasy", "arm": "exfor_crystal",
                        "fc_npz": f"cache/mp/force_constants/{npz.name}" if npz.exists() else None})   # older-layout label, resolved by file name
        else:   # no MP phonons: MLIP-only (structure from the summary entry, supercell chosen by the MLIP script)
            r = e.iloc[0]
            sel.append({"material_id": r.material_id, "formula": f, "sg_symbol": sg, "nsites": int(r.nsites),
                        "phonon_id": None, "phonon_method": None, "arm": "exfor_crystal_mlip_only", "fc_npz": None})
    out = {"seed": SEED, "strata": [[list(k), n] for k, n in STRATA], "pool_size": len(pool), "n": len(sel),
           "compounds": sel}
    (F("stage2_mlip_set.json")).write_text(json.dumps(out, indent=1))
    print(len(pool), "eligible;", len(sel), "selected;",
          pd.Series([c["arm"] for c in sel]).value_counts().to_dict())


if __name__ == "__main__":
    main()
