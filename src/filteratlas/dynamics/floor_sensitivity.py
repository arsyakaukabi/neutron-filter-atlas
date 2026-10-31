"""Triviality-test statistic at VDOS floors of 1 and 3 meV (the adopted value is 2 meV).

Recomputes full-model L0.8 for every usable compound (current usable set) at floors 1 and 3 meV; the proxy does
not depend on the floor. Also lists any compound whose conversion fails at either floor.
-> stage1r_floor_sensitivity.json (+ .parquet with per-compound values)
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from filteratlas.common.paths import D, F
from filteratlas.common.config import cfg, tup
from filteratlas.common.constants import N_WORKERS, TEMPS  # noqa: E402

FLOORS_SENS = tup(cfg("vdos")["floor_sensitivity_meV"])   # triviality test at 1 and 3 meV floors (adopted MP DOS floor: 2 meV)


def job(args):
    mid, meth, pid, floor = args
    from filteratlas.ingest import mp as mq
    from filteratlas.xs import ncmat as nb
    from filteratlas.xs import evaluate as xe
    try:
        dos = mq.load_dos(D("mp") / "phonon_dos" / meth / f"{pid}.json.gz")
        ev = xe.evaluate(nb.build_vdos_material(dos, f"fs{floor}_{mid}", floor_meV=floor), TEMPS)
        return {"material_id": mid, "floor": floor, "L08_293K": ev["L08_293K"], "L08_77K": ev["L08_77K"], "ok": True}
    except Exception as e:  # noqa: BLE001
        return {"material_id": mid, "floor": floor, "ok": False, "error": repr(e)[:150]}


def main():
    u = pd.read_parquet(F("stage1_usable.parquet"))
    jobs = [(r.material_id, r.phonon_method, r.phonon_id, fl) for fl in FLOORS_SENS for r in u.itertuples()]
    with Pool(N_WORKERS) as p:
        rows = pd.DataFrame(p.map(job, jobs, chunksize=8))
    rows = rows.reindex(columns=["material_id", "floor", "L08_293K", "L08_77K", "ok", "error"])
    rows.to_parquet(F("stage1r_floor_sensitivity.parquet"))
    res = {"n_usable": len(u), "n_failures": int((~rows.ok).sum()),
           "failures": rows[~rows.ok][["material_id", "floor", "error"]].to_dict("records")}
    for fl in FLOORS_SENS:
        r = rows[(rows.floor == fl) & rows.ok].merge(u[["material_id", "L08_proxy_293K", "L08_proxy_77K"]],
                                                     on="material_id")
        res[f"floor_{fl:g}meV"] = {t: {"rho": round(float(spearmanr(r[f"L08_{t}"], r[f"L08_proxy_{t}"]).statistic), 4),
                                       "n": len(r)} for t in ("293K", "77K")}
    res["floor_2meV_adopted"] = {t: {"rho": round(float(spearmanr(u[f"L08_{t}"], u[f"L08_proxy_{t}"]).statistic), 4),
                                     "n": len(u)} for t in ("293K", "77K")}
    (F("stage1r_floor_sensitivity.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
