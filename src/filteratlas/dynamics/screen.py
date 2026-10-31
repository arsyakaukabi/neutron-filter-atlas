"""Atlas funnel: dynamical-stability screen and Debye scores of the survivors of the proxy bound.

Two entry points:
- the dynamical-stability screen on force constants (40 A Gamma-centred mesh + supercell-commensurate q,
  symmetrised MP force constants) of the MP-phonon survivors that pass the DOS rule and were not screened
  before; merged with the earlier screen for every MP-phonon survivor -> stage3_stability.parquet (step 02);
- debye: for every survivor without MP phonons: SevenNet supercell size (diagonal, >= 12 A per axis, as
  mlip_run.start_cell) and, with K_VRH > 0 and G_VRH > 0, the Debye-from-elastic L0.8 at 293 / 77 K
  (MP summary structure) -> stage3_debye.parquet (step 03).
"""
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.constants import N_WORKERS, TEMPS  # noqa: E402
from filteratlas.common.paths import F

SC_MIN_A = cfg("stability")["supercell_min_length_A"]   # minimum supercell length per axis
UNSTABLE = cfg("stability")["fc_unstable_THz"]          # dynamical-stability threshold on force-constant frequencies


def d7a_one(r):
    from filteratlas.dynamics.stability_screen import min_freq
    rec = {"material_id": r["material_id"]}
    try:
        fmin, fgam, _, fcomm = min_freq(r["phonon_id"])
        rec.update({"fc_min_freq_THz": fmin, "fc_commensurate_min_THz": fcomm, "fc_status": "ok"})
    except FileNotFoundError:
        rec["fc_status"] = "no_force_constants"
    except Exception as e:  # noqa: BLE001
        rec.update({"fc_status": "error", "fc_error": repr(e)[:200]})
    return rec


def d7a():
    from filteratlas.ingest.structures_fc import d7a_todo
    todo = d7a_todo()
    rows = todo[todo.phonon_method == "pheasy"].to_dict("records")
    t0 = time.time()
    with Pool(N_WORKERS) as pool:
        new = pd.DataFrame(pool.map(d7a_one, rows, chunksize=4))
    print(f"screened {len(new)} in {time.time() - t0:.0f}s", flush=True)
    p = pd.read_parquet(F("stage3_proxy.parquet"))
    s = p[~p.excluded_by_proxy_bound & (p.phonon_id != "")][["material_id", "formula", "phonon_method", "phonon_id"]]
    t = pd.read_parquet(F("stage1_triviality.parquet"), columns=["material_id", "status", "stable", "unstable_weight"])
    old = pd.read_parquet(F("stage1r_stability_full.parquet"))
    old = old.rename(columns={"status": "fc_status"})[["material_id", "fc_status", "fc_min_freq_THz", "fc_commensurate_min_THz"]]
    old["d7a_source"] = "stage1r_stability_full"
    new["d7a_source"] = "stage3_screen"
    scr = pd.concat([old, new], ignore_index=True).drop_duplicates("material_id", keep="last")
    m = s.merge(t, on="material_id", how="left").merge(scr, on="material_id", how="left")
    m["dos_convert_ok"] = m.status == "ok"
    m["dos_rule_pass"] = m.stable.fillna(False).astype(bool)
    # no force constants available (DFPT-only entries, pheasy entries MP returns none for, and the 8 step-02
    # "FileNotFoundError" entries, which are the same case) -> unscreenable
    m["unscreenable"] = m.dos_convert_ok & m.dos_rule_pass & (m.fc_status != "ok")
    ok = m.fc_status == "ok"
    m["unstable_comm"] = ok & (m.fc_commensurate_min_THz < UNSTABLE)
    m["interp_only"] = ok & (m.fc_min_freq_THz < UNSTABLE) & ~m.unstable_comm
    m["unstable_dft"] = m.dos_convert_ok & (~m.dos_rule_pass | m.unstable_comm)
    m["d7a_pass"] = m.dos_convert_ok & m.dos_rule_pass & ~m.unstable_comm
    m = m.sort_values("material_id").reset_index(drop=True)
    m.to_parquet(F("stage3_stability.parquet"))
    print(m[["dos_convert_ok", "dos_rule_pass", "unscreenable", "unstable_comm", "interp_only", "unstable_dft",
             "d7a_pass"]].sum().to_dict(), m.fc_status.value_counts(dropna=False).to_dict())


def debye_one(r):
    from filteratlas.ingest import mp as mq
    from filteratlas.xs import ncmat as nb
    from filteratlas.xs import evaluate as xe
    rec = {"material_id": r["material_id"]}
    s = mq.structures([r["material_id"]]).get(r["material_id"])
    if s is None:
        rec["structure_status"] = "missing"
        return rec
    n = [max(1, int(np.ceil(SC_MIN_A / L))) for L in s.lattice.abc]
    rec.update({"structure_status": "ok", "nsites_struct": len(s), "sc_diag": n, "supercell_atoms": len(s) * int(np.prod(n))})
    if r["K_vrh"] and r["G_vrh"] and r["K_vrh"] > 0 and r["G_vrh"] > 0:
        try:
            th = nb.debye_temperature(r["K_vrh"], r["G_vrh"], r["density"], r["nsites"], r["volume"])
            rec["theta_D_elastic"] = th
            cfg = nb.build_debye_material(s, th, f"s3deb_{r['material_id']}")
            for T in TEMPS:
                mu, info = xe.mu_offbragg(cfg, T, xe.LAM)
                tag = "293K" if T > 100 else "77K"
                rec[f"debye_L08_{tag}"] = xe.l_target(mu)
            rec["debye_status"] = "ok"
        except Exception as e:  # noqa: BLE001
            rec.update({"debye_status": "error", "debye_error": repr(e)[:200]})
    else:
        rec["debye_status"] = "no_elastic"
    return rec


def _init():
    # structures are read once per worker from the cache (mq.structures re-reads the file on every call otherwise)
    import gzip
    import json
    from filteratlas.ingest import mp as mq
    from pymatgen.core import Structure
    have = {}
    with gzip.open(mq.MPC / "structures.jsonl.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            have[d["material_id"]] = d["structure"]
    mq.structures = lambda ids, refresh=False: {i: Structure.from_dict(have[i]) for i in ids if i in have}


def debye():
    p = pd.read_parquet(F("stage3_proxy.parquet"))
    c = pd.read_parquet(F("stage0_candidates.parquet"))[["material_id", "density", "volume"]]
    s = p[~p.excluded_by_proxy_bound & (p.phonon_id == "")].merge(c, on="material_id")
    t0 = time.time()
    with Pool(N_WORKERS, initializer=_init) as pool:
        d = pd.DataFrame(pool.map(debye_one, s.to_dict("records"), chunksize=8))
    print(f"debye/supercell for {len(d)} in {time.time() - t0:.0f}s", flush=True)
    d = d.sort_values("material_id").reset_index(drop=True)
    d.to_parquet(F("stage3_debye.parquet"))
    print(d.structure_status.value_counts().to_dict(), d.debye_status.value_counts().to_dict(),
          "supercell>480:", int((d.supercell_atoms > 480).sum()))


if __name__ == "__main__":
    {"d7a": d7a, "debye": debye}[sys.argv[1]]()   # the pipeline steps call the two functions directly
