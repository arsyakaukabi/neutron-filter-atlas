"""Funnel selection, fixed before the computation: the SevenNet list (budget 1,500, priorities a-d), the
100-compound drift set, and the funnel counts.

Inputs (data/derived): the proxy bound, the force-constant stability screen, the Debye scores, the tier table,
the triviality-test table and the MLIP benchmark set.
-> the SevenNet set JSON, the funnel table (one row per candidate: funnel label, priority),
   data/derived/atlas_counts.json (funnel part; the atlas evaluation adds the rest)
Usage: python pipelines/03_dynamics_mlip.py --only funnel             (selection)
       python -m filteratlas.atlas.funnel --not-computed [N]   (best not_computed_budget compounds by proxy and Debye)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F  # noqa: E402

_FU = cfg("funnel")
BUDGET, N_DRIFT, SEED, MAX_SC = _FU["sevennet_budget"], _FU["drift_set_size"], _FU["seed"], _FU["max_supercell_atoms"]
FCDIR = _FU["fc_dir_label"]   # directory label stored in the SevenNet set JSON (older directory layout)
P80 = _FU["debye_priority_percentile"]
SC_MIN_A = cfg("stability")["supercell_min_length_A"]


def load():
    p = pd.read_parquet(F("stage3_proxy.parquet"))
    st = pd.read_parquet(F("stage3_stability.parquet"))
    db = pd.read_parquet(F("stage3_debye.parquet"))
    tiers = pd.read_parquet(F("stage3_tiers.parquet"))
    t1 = pd.read_parquet(F("stage1_triviality.parquet"), columns=["material_id", "L08_293K", "L08_77K"])
    return p, st, db, tiers, t1


def main():
    p, st, db, tiers, t1 = load()
    tierA = set(tiers.material_id)
    surv = p[~p.excluded_by_proxy_bound]
    # step 3: MP DFT survivors after the force-constant screen; P80 of their full-model L0.8 per temperature
    dft = st.merge(t1, on="material_id", how="left")
    dpass = dft[dft.d7a_pass]
    p80 = {T: float(np.percentile(dpass[f"L08_{T}"], P80)) for T in ("293K", "77K")}
    # step 4-5: survivors without MP phonons
    nph = surv[surv.phonon_id == ""].merge(db, on="material_id", how="left")
    nph["tierA"] = nph.material_id.isin(tierA)
    nph["debye_ok"] = nph.debye_status == "ok"
    nph["debye_score"] = np.where(nph.debye_ok, np.maximum(nph.debye_L08_293K / p80["293K"],
                                                           nph.debye_L08_77K / p80["77K"]), np.nan)
    nph["too_large"] = nph.supercell_atoms > MAX_SC
    pr = pd.Series("d", index=nph.index)
    pr[~nph.debye_ok] = "c"
    pr[nph.debye_ok & (nph.debye_score >= 1.0)] = "b"
    pr[nph.tierA] = "a"
    nph["priority"] = pr
    nph["order_key"] = np.where(nph.priority == "c", -nph.proxy_L08,
                                np.where(nph.priority == "a", 0.0, -nph.debye_score))
    nph = nph.sort_values(["priority", "order_key", "material_id"]).reset_index(drop=True)
    eligible = nph[~nph.too_large]
    chosen = eligible.head(BUDGET)
    nph["funnel_label"] = np.where(nph.too_large, "supercell_too_large",
                                   np.where(nph.material_id.isin(chosen.material_id), "sevennet_selected",
                                            "not_computed_budget"))
    # drift set: MP-DFT survivors passing the force-constant screen, not in the MLIP benchmark set, MP supercell <= 480 atoms, FC cached
    s2 = {c["material_id"] for c in json.loads((F("stage2_mlip_set.json")).read_text())["compounds"]}
    pool = dpass[(dpass.phonon_method == "pheasy") & ~dpass.material_id.isin(s2) & (dpass.fc_status == "ok")]
    pool = pool.sort_values("material_id")
    sizes = []
    for r in pool.itertuples():
        z = np.load(D("mp") / "force_constants" / f"{r.phonon_id}.npz")
        sizes.append(int(round(len(z["symbols"]) * abs(np.linalg.det(z["sc"])))))
    pool = pool.assign(supercell_atoms=sizes)
    pool = pool[pool.supercell_atoms <= MAX_SC]
    drift = pool.sample(N_DRIFT, random_state=SEED).sort_values("material_id")
    comp = []
    for r in chosen.itertuples():
        comp.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol,
                     "nsites": int(r.nsites_struct), "supercell_atoms": int(r.supercell_atoms),
                     "arm": f"atlas_{r.priority}", "fc_npz": None, "tierA": bool(r.tierA)})
    # Tier-A compounds above the 480-atom cap (the four garnets) run last, alone,
    # with the same protocol and supercell rule; they stay flagged supercell_too_large if they cannot run
    for r in nph[nph.tierA & nph.too_large].itertuples():
        comp.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol,
                     "nsites": int(r.nsites_struct), "supercell_atoms": int(r.supercell_atoms),
                     "arm": "atlas_a_large_D15", "fc_npz": None, "tierA": True})
    # Tier-A compounds whose MP DOS cannot be converted get SevenNet too (outside the budget)
    conv = dft[(~dft.dos_convert_ok) & dft.material_id.isin(tierA)].merge(
        pd.read_parquet(F("stage3_debye.parquet"))[["material_id"]], on="material_id", how="left")
    from filteratlas.ingest import mp as mq
    st16 = mq.structures(sorted(conv.material_id))
    for r in conv.itertuples():
        s_ = st16[r.material_id]
        n16 = [max(1, int(np.ceil(SC_MIN_A / L))) for L in s_.lattice.abc]
        comp.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": None, "nsites": len(s_),
                     "supercell_atoms": int(len(s_) * np.prod(n16)), "arm": "atlas_a_D16", "fc_npz": None, "tierA": True})
    for r in drift.itertuples():
        comp.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": None, "nsites": None,
                     "supercell_atoms": int(r.supercell_atoms), "arm": "drift", "phonon_id": r.phonon_id,
                     "fc_npz": f"{FCDIR}/{r.phonon_id}.npz", "tierA": r.material_id in tierA})
    meta = {"budget": BUDGET, "n_drift": N_DRIFT, "seed": SEED, "max_supercell_atoms": MAX_SC,
            "p80_dft_L08_cm": p80, "n": len(comp),
            "by_arm": pd.Series([c["arm"] for c in comp]).value_counts().to_dict(),
            "supercell_atoms_selected": pd.Series([c["supercell_atoms"] for c in comp]).describe().round(1).to_dict()}
    (F("stage3_mlip_set.json")).write_text(json.dumps({**meta, "compounds": comp}, indent=1))
    # funnel table and counts
    f = p[["material_id", "formula", "sg_symbol", "excluded_by_proxy_bound", "proxy_L08", "proxy_T_1cm"]].copy()
    f = f.merge(dft[["material_id", "dos_convert_ok", "d7a_pass", "unstable_dft", "interp_only", "unscreenable"]],
                on="material_id", how="left")
    f = f.merge(nph[["material_id", "priority", "debye_score", "debye_L08_293K", "debye_L08_77K", "supercell_atoms",
                     "funnel_label"]], on="material_id", how="left")
    lab = np.where(f.excluded_by_proxy_bound, "excluded_by_proxy_bound",
                   np.where(f.funnel_label.notna(), f.funnel_label,
                            np.where(f.d7a_pass == True, "mp_dft",  # noqa: E712
                                     np.where(f.unstable_dft == True, "unstable_dft", "mp_dos_conversion_failed"))))  # noqa: E712
    f["funnel_step"] = lab
    f["tierA"] = f.material_id.isin(tierA)
    f = f.sort_values("material_id").reset_index(drop=True)
    f.to_parquet(F("stage3_funnel.parquet"))
    counts = {"step1_practical_set": len(p),
              "step2_excluded_by_proxy_bound": int(p.excluded_by_proxy_bound.sum()),
              "step2_survivors": len(surv),
              "step3_with_mp_phonons": int((surv.phonon_id != "").sum()),
              "step3_dos_conversion_failed": int((~dft.dos_convert_ok).sum()),
              "step3_unstable_dft_excluded": int(dft.unstable_dft.sum()),
              "step3_unstable_dft_by_dos_rule": int((dft.dos_convert_ok & ~dft.dos_rule_pass).sum()),
              "step3_unstable_dft_by_commensurate_q": int((dft.dos_convert_ok & dft.dos_rule_pass & dft.unstable_comm).sum()),
              "step3_d7a_pass": int(dft.d7a_pass.sum()), "step3_interp_only": int((dft.d7a_pass & dft.interp_only).sum()),
              "step3_unscreenable": int((dft.d7a_pass & dft.unscreenable).sum()),
              "step3_p80_L08_cm": p80,
              "step4_without_mp_phonons": len(nph), "step4_with_debye_score": int(nph.debye_ok.sum()),
              "step4_debye_errors": int((nph.debye_status == "error").sum()),
              "step5_supercell_too_large": int(nph.too_large.sum()),
              "step5_priority_available": nph[~nph.too_large].priority.value_counts().sort_index().to_dict(),
              "step5_selected": chosen.priority.value_counts().sort_index().to_dict(),
              "step5_selected_total": len(chosen), "step5_not_computed_budget": int((nph.funnel_label == "not_computed_budget").sum()),
              "step5_tierA_without_mp_phonons": int(nph.tierA.sum()),
              "step5_tierA_too_large": int((nph.tierA & nph.too_large).sum()),
              "drift_set": len(drift), "drift_pool": len(pool)}
    old = json.loads((F("atlas_counts.json")).read_text()) if (F("atlas_counts.json")).exists() else {}
    old["funnel"] = counts
    (F("atlas_counts.json")).write_text(json.dumps(old, indent=1))
    print(json.dumps({**counts, "set": meta}, indent=1))


def not_computed(n=20):
    f = pd.read_parquet(F("stage3_funnel.parquet"))
    nc = f[f.funnel_step == "not_computed_budget"]
    cols = ["material_id", "formula", "sg_symbol", "proxy_L08", "debye_L08_293K", "debye_L08_77K", "debye_score", "priority"]
    out = {"by_proxy_L08": nc.nlargest(n, "proxy_L08")[cols].round(3).to_dict("records"),
           "by_debye_score": nc.dropna(subset=["debye_score"]).nlargest(n, "debye_score")[cols].round(3).to_dict("records")}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    if "--not-computed" in sys.argv:
        not_computed(int(sys.argv[-1]) if sys.argv[-1].isdigit() else 20)
    else:
        main()
