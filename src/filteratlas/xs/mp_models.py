"""MP-phonon models for the triviality test: full model vs elemental proxy vs Debye-from-elastic, per compound.

For every candidate with an MP phonon DOS: stability check, MP-VDOS NCrystal model, mu at 1.0/1.8 A,
band-weighted L0.8 (thermal and hot band) at 293.15 K and 77 K, proxy analogues, and the Debye-from-elastic
ablation where VRH moduli exist. Positive and negative controls are forced through the same code.
Output: stage1_triviality.parquet, stage1_controls.parquet
"""
import json
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.config import cfg
from filteratlas.common.constants import N_WORKERS, TEMPS  # noqa: E402
from filteratlas.common.paths import D, F

DOSDIR = D("mp") / "phonon_dos"
UNIV_THETA = cfg("figures_of_merit")["universal_debye_K"]
CTRL = cfg("controls")["stage1"]


def one(row):
    rec = {k: row[k] for k in ("material_id", "formula", "sg_symbol", "spacegroup", "nelements", "nsites",
                               "flags", "phonon_method", "phonon_id", "K_vrh", "G_vrh")}
    t0 = time.time()
    try:
        dos = mq.load_dos(DOSDIR / row["phonon_method"] / f"{row['phonon_id']}.json.gz")
        w, ok = nb.stability(dos)
        rec.update({"unstable_weight": w, "stable": ok})
        cfg = nb.build_vdos_material(dos, f"vdos_{row['material_id']}")
        rec.update(xe.evaluate(cfg, TEMPS))
        rec["mean_phonon_meV"] = mean_phonon(dos)
        for T in TEMPS:   # stored curves for band / moderator sensitivity
            tag = "293K" if T > 100 else "77K"
            rec[f"curve_full_{tag}"], rec[f"curve_proxy_{tag}"] = xe.curves(cfg, T)
        for th in UNIV_THETA:   # universal-theta_D nulls: composition + density + one fixed stiffness
            ucfg = nb.build_debye_material(dos["structure"], th, f"univ{th:.0f}_{row['material_id']}")
            rec.update({f"univ{th:.0f}_{k}": v for k, v in xe.evaluate(ucfg, TEMPS).items() if k.startswith("L08")})
        if row["K_vrh"] and row["G_vrh"] and row["K_vrh"] > 0 and row["G_vrh"] > 0:
            th = nb.debye_temperature(row["K_vrh"], row["G_vrh"], row["density"], row["nsites"], row["volume"])
            rec["theta_D_elastic"] = th
            dcfg = nb.build_debye_material(dos["structure"], th, f"deb_{row['material_id']}")
            rec.update({f"debye_{k}": v for k, v in xe.evaluate(dcfg, TEMPS).items() if k.startswith(("L08", "mu"))})
        rec["status"] = "ok"
    except Exception as e:  # noqa: BLE001
        rec["status"] = "error"
        rec["error"] = (repr(e) + " | " + traceback.format_exc().splitlines()[-2])[:300]
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def mean_phonon(dos):
    """Mean phonon energy of a DOS (first moment over positive frequencies).

    Parameters:
        dos: phonon DOS dict with "frequencies" (THz) and "densities".
    Returns:
        mean energy in meV, or nan if the DOS has no positive-frequency weight.
    Example:
        >>> from filteratlas.xs.mp_models import mean_phonon
        >>> mean_phonon(dos)
    """
    import numpy as np
    f = np.asarray(dos["frequencies"], float)
    d = np.clip(np.asarray(dos["densities"], float), 0, None) * (f > 0)
    return float(np.trapezoid(d * f, f) / np.trapezoid(d, f) * 4.135667696) if d.sum() > 0 else float("nan")


def controls():
    """Positive (sapphire, Si, alpha-quartz, Bi) and negative (CdS, LiF, Gd2O3) controls."""
    summ = mq.summary_table()
    rows = []
    for name, ncm in CTRL["positive_curated"]:
        rows.append({"control": name, "kind": "positive", "source": "ncrystal_stdlib", **xe.evaluate(ncm, TEMPS)})
    # Positive controls through the MP pipeline under test
    cand = pd.read_parquet(F("stage0_candidates.parquet"))
    for name, f, sg in CTRL["positive_mp"]:
        r = cand[(cand.formula == f) & (cand.sg_symbol == sg) & (cand.phonon_id != "")].sort_values("e_hull").iloc[0]
        dos = mq.load_dos(DOSDIR / r.phonon_method / f"{r.phonon_id}.json.gz")
        cfg = nb.build_vdos_material(dos, f"posmp_{r.material_id}")
        rows.append({"control": name, "kind": "positive", "source": f"mp_{r.phonon_method}",
                     "material_id": r.material_id, **xe.evaluate(cfg, TEMPS)})
    neg = summ[summ.formula.isin(CTRL["negative"]) & ((summ.phonon_pheasy != "") | (summ.phonon_dfpt != ""))]
    neg = neg.sort_values("e_hull").groupby("formula").head(1)
    pairs = [((r.phonon_pheasy or r.phonon_dfpt).split(",")[0], "pheasy" if r.phonon_pheasy else "dfpt")
             for r in neg.itertuples()]
    paths = mq.phonon_dos_bulk(pairs)
    # Negative controls with no MP phonons (e.g. Gd2O3): forced through with a Debye model (absorption dominates,
    # so the dynamics model is irrelevant to where they rank); labelled as such.
    for f in sorted(set(CTRL["negative"]) - set(neg.formula)):
        r = summ[summ.formula == f].sort_values("e_hull").iloc[0]
        th = nb.debye_temperature(r.K_vrh, r.G_vrh, r.density, r.nsites, r.volume) if r.G_vrh else float("nan")
        th = th if th == th else CTRL["negative_debye_fallback_K"]
        s = mq.structures([r.material_id])[r.material_id]
        cfg = nb.build_debye_material(s, th, f"neg_{r.material_id}")
        rows.append({"control": f, "kind": "negative", "source": f"debye_{th:.0f}K (no MP phonons)",
                     "material_id": r.material_id, **xe.evaluate(cfg, TEMPS)})
    for r, (pid, meth) in zip(neg.itertuples(), pairs):
        if (pid, meth) not in paths:
            rows.append({"control": r.formula, "kind": "negative", "source": "missing_dos"})
            continue
        dos = mq.load_dos(paths[(pid, meth)])
        cfg = nb.build_vdos_material(dos, f"neg_{r.material_id}")
        rows.append({"control": r.formula, "kind": "negative", "source": f"mp_{meth}", "material_id": r.material_id,
                     **xe.evaluate(cfg, TEMPS)})
    return pd.DataFrame(rows)


def main():
    cand = pd.read_parquet(F("stage0_candidates.parquet"))
    cand = cand[cand.phonon_id != ""]
    have = {(p.parent.name, p.name.replace(".json.gz", "")) for p in DOSDIR.glob("*/*.json.gz")}
    cand = cand[[(m, i) in have for m, i in zip(cand.phonon_method, cand.phonon_id)]]
    rows = cand.to_dict("records")
    t0 = time.time()
    with Pool(N_WORKERS) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(one, rows, chunksize=4)):
            out.append(r)
            if i % 250 == 0:
                print(f"{i}/{len(rows)} {time.time() - t0:.0f}s", flush=True)
    # imap_unordered returns rows in completion order; sort so every downstream sample/bootstrap is reproducible
    df = pd.DataFrame(out).sort_values("material_id").reset_index(drop=True)
    df.to_parquet(F("stage1_triviality.parquet"))
    ctrl = controls()
    ctrl.to_parquet(F("stage1_controls.parquet"))
    print(json.dumps({"n": len(df), "ok": int((df.status == "ok").sum()), "stable": int(df.stable.fillna(False).sum()),
                      "errors": df[df.status == "error"].error.str[:80].value_counts().head(8).to_dict(),
                      "wall_s": round(time.time() - t0)}, indent=1))


if __name__ == "__main__":
    main()
