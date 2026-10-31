"""Atlas summary numbers that are derived (lists, counts, splits) rather than read directly from one output file.
Everything is computed here from data/derived and the run logs; nothing is typed by hand.
-> the summary-numbers JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402

TAGS = {"293K": 293.15, "77K": 77.0}


def names(df):
    return [f"{r.formula} {r.sg_symbol}" for r in df.itertuples()]


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    one = at[at["T"] > 100]
    fm = pd.read_parquet(F("stage3_fronts.parquet"))
    fj = json.loads((F("stage3_fronts.json")).read_text())
    sset = json.loads((F("stage3_mlip_set.json")).read_text())["compounds"]
    out = {}
    # --- compute record (SevenNet run log) ------------------------------------------------------------------------
    log = [json.loads(x) for x in open(D("logs") / "mlip_atlas.jsonl", encoding="utf-8") if x.strip()]
    last = {}
    for r in log:
        last[r["material_id"]] = r
    ok = [r for r in last.values() if r.get("status") == "ok"]
    outs = {p.name.replace(".json.gz", "") for p in (D("mlip_atlas") / "sevennet_mf_ompa").glob("*.json.gz")}
    drift = {c["material_id"] for c in sset if c["arm"] == "drift"}
    out["compute"] = {
        "n_outputs": len(outs), "n_outputs_atlas": len(outs - drift), "n_outputs_drift": len(outs & drift),
        "n_set": len(sset), "n_ok_gpu": sum(r["device"] == "cuda" for r in ok), "n_ok_cpu": sum(r["device"] == "cpu" for r in ok),
        "median_wall_s_gpu": round(float(np.median([r["wall_s"] for r in ok if r["device"] == "cuda"])), 1),
        "median_wall_s_cpu": round(float(np.median([r["wall_s"] for r in ok if r["device"] == "cpu"])), 1),
        "max_supercell_atoms_ok_gpu": max(r["supercell_atoms"] for r in ok if r["device"] == "cuda"),
        "n_gpu_out_of_memory_last_try": sum(r.get("status") == "error" and "CUDA" in r.get("error", "") for r in last.values()
                                            if r["material_id"] not in outs),
        "n_logged_tries": len(log)}
    # host label: the "machine" field where a record carries one, otherwise a label from the device
    mach = lambda r: r.get("machine", {"cuda": "laptop GPU", "cpu": "laptop CPU"}[r["device"]])  # noqa: E731
    out["compute"]["ok_by_machine"] = pd.Series([mach(r) for r in ok]).value_counts().to_dict()
    out["compute"]["ok_gpu_raised_cap"] = sum("gpu_alloc" in r and r["device"] == "cuda" for r in ok)
    out["compute"]["median_wall_s_by_machine"] = {m: round(float(np.median([r["wall_s"] for r in ok if mach(r) == m])), 1)
                                                  for m in out["compute"]["ok_by_machine"]}
    for k, m in (("laptop_gpu", "laptop GPU"), ("laptop_cpu", "laptop CPU"), ("vm_gpu", "VM T4 GPU, llm-06"),
                 ("vm_cpu", "VM CPU, llm-06")):   # plain keys for the report templates
        out["compute"][f"ok_{k}"] = out["compute"]["ok_by_machine"].get(m, 0)
        out["compute"][f"median_wall_s_{k}"] = out["compute"]["median_wall_s_by_machine"].get(m)
    out["compute"]["max_supercell_atoms_ok_by_machine"] = {m: max(r["supercell_atoms"] for r in ok if mach(r) == m)
                                                           for m in out["compute"]["ok_by_machine"]}
    # comparison of the same compound's relaxation across hosts (every try logs its relaxation, also before a phonon
    # out-of-memory)
    tries = {}
    for r in log:
        if r.get("volume_ratio") is not None and r.get("relax_converged") is not None:
            tries.setdefault(r["material_id"], []).append(r)
    pairs, flips = {}, []
    dv = []
    for m, rs in tries.items():
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                a_, b_ = mach(rs[i]), mach(rs[j])
                if a_ == b_:
                    continue
                key = " vs ".join(sorted((a_, b_)))
                pairs[key] = pairs.get(key, 0) + 1
                if rs[i]["relax_converged"] and rs[j]["relax_converged"]:
                    dv.append(abs(np.log(rs[i]["volume_ratio"] / rs[j]["volume_ratio"])))
                elif rs[i]["relax_converged"] != rs[j]["relax_converged"]:
                    flips.append({"compound": f"{rs[i]['formula']} ({m})", "supercell_atoms": rs[i].get("supercell_atoms"),
                                  "converged_on": [mach(x) for x in (rs[i], rs[j]) if x["relax_converged"]],
                                  "not_converged_on": [mach(x) for x in (rs[i], rs[j]) if not x["relax_converged"]],
                                  "steps": [x.get("relax_steps") for x in (rs[i], rs[j])]})
    out["compute"]["cross_machine_relaxation"] = {
        "pairs_by_machines": pairs, "n_pairs": sum(pairs.values()), "n_both_converged": len(dv),
        "median_abs_dlnV": float(f"{np.median(dv):.2g}") if dv else None, "max_abs_dlnV": float(f"{max(dv):.2g}") if dv else None,
        "n_convergence_disagree": len(flips), "convergence_flips": flips}
    # --- funnel outcome ---------------------------------------------------------------------------------------------
    nrm = one[one.funnel_step == "not_run_memory"]
    out["not_run_memory"] = {"n": len(nrm), "tierA": names(nrm[nrm.tier == "A"].sort_values("formula")),
                             "n_tierA": int((nrm.tier == "A").sum()), "n_tierB": int((nrm.tier == "B").sum()),
                             "drift_not_run": len(drift - outs)}
    a = one[one.tier == "A"]
    out["tierA"] = {"n": len(a), "excluded_by_proxy_bound": names(a[a.funnel_step == "excluded_by_proxy_bound"].sort_values("formula")),
                    "unstable_dft": names(a[a.funnel_step == "unstable_dft"].sort_values("formula")),
                    "full_model": int(a.vdos_source.ne("none").sum()), "sevennet": int((a.vdos_source == "sevennet").sum()),
                    "mp_dft": int((a.vdos_source == "mp_dft").sum()),
                    "outside_domain": names(a[a.outside_domain].sort_values("formula")),
                    "untested_heavy_main": int((a.main_front_eligible & a.untested_heavy).sum()),
                    "phase_mismatch_full_model": names(a[a.vdos_source.ne("none") & a.tierA_phase_mismatch].sort_values("formula"))}
    # --- fronts: members and within-1-sigma sets with flags ----------------------------------------------------------
    fr = {}
    for key, v in fj["fronts"].items():
        T = TAGS[key.split("|")[0]]
        w = fm[(fm.key == key) & fm.within_1sigma][["material_id"]].merge(at[at["T"] == T], on="material_id")
        mem = [x["material_id"] for x in v["front"]]
        fr[key] = {"front": [f"{x['formula']} {x['sg_symbol']}" for x in v["front"]],
                   "front_sevennet": [f"{x['formula']} {x['sg_symbol']}" for x in v["front"] if x["vdos_source"] == "sevennet"],
                   "front_flags": {f"{r.formula} {r.sg_symbol}": [c for c in ("interp_only", "tierA_phase_mismatch",
                                   "untested_heavy", "hygroscopic", "flag_Be", "flag_Pb", "flag_Bi", "mlip_unstable") if bool(getattr(r, c))]
                                   for r in at[(at["T"] == T) & at.material_id.isin(mem)].itertuples()},
                   "within_1sigma_not_front": names(w[~w.material_id.isin(mem)].sort_values("A_fast_1E", ascending=False)),
                   "n_within_1sigma": len(w)}
    out["fronts"] = {k.replace("|", "_"): v for k, v in fr.items()}   # "|" is a format separator of the report templates
    out["front_sizes"] = {k.replace("|", "_"): {"n_eligible": v["n_eligible"], "n_front": v["n_front"],
                                                "n_within_1sigma": v["n_within_1sigma"]} for k, v in fj["fronts"].items()}
    out["ladder"] = {k.replace("|", "_"): v for k, v in fj["ladder"].items()}
    out["sensitivity"] = fj["sensitivity"]
    out["kill3"] = fj["kill3_preview"]
    out["who_moves"] = fj["who_moves"]
    out["who_moves_text"] = {k: [f"{r['formula']} {r['sg_symbol']}" for r in fj["who_moves"][k][:4]]
                             for k in ("up_proxy_to_full_293K", "down_proxy_to_full_293K", "up_293_to_77", "down_293_to_77")}
    # --- best compounds without a full model --------------------------------------------------------
    fun = pd.read_parquet(F("stage3_funnel.parquet"))
    def top(step, col, n=5):
        d = fun[fun.funnel_step == step] if step != "not_run_memory" else fun[fun.material_id.isin(nrm.material_id)]
        d = d.dropna(subset=[col]).nlargest(n, col)
        return [{"compound": f"{r.formula} {r.sg_symbol}", "proxy_L08": round(float(r.proxy_L08), 2),
                 "debye_L08_293K": None if r.debye_L08_293K != r.debye_L08_293K else round(float(r.debye_L08_293K), 2),
                 "debye_L08_77K": None if r.debye_L08_77K != r.debye_L08_77K else round(float(r.debye_L08_77K), 2)}
                for r in d.itertuples()]
    out["referee_query"] = {"not_computed_budget_by_proxy": top("not_computed_budget", "proxy_L08"),
                            "not_computed_budget_by_debye": top("not_computed_budget", "debye_L08_293K"),
                            "not_run_memory_by_proxy": top("not_run_memory", "proxy_L08"),
                            "not_run_memory_by_debye": top("not_run_memory", "debye_L08_293K")}
    full = one[one.main_front_eligible]
    for t, T in TAGS.items():
        ft = at[(at["T"] == T) & at.main_front_eligible]
        out["referee_query"][f"max_full_L08_main_{t}"] = round(float(ft.L08.max()), 2)
        failed = one[one.funnel_step == "mlip_failed"].material_id
        nc = fun[fun.funnel_step.isin(["not_computed_budget"]) | fun.material_id.isin(nrm.material_id)
                 | fun.material_id.isin(failed)]
        col = f"debye_L08_{t}"
        out["referee_query"][f"n_without_full_model_debye_above_max_full_{t}"] = int((nc[col] > ft.L08.max()).sum())
        out["referee_query"][f"n_without_full_model_debye_above_P80_full_{t}"] = int(
            (nc[col] > np.percentile(ft.L08, 80)).sum())
    out["referee_query"]["n_full_model_main"] = len(full)
    # --- coherent check, Bragg, OpenMC summaries -------------------------------------------------------------------------
    coh = json.loads((F("stage3_coherent.json")).read_text())
    out["coherent"] = {x: {t: {"dlnL08": round(float(coh["crystals"][x][t]["dlnL08"]), 5),
                               "band_mean_rel_change_mu": round(float(coh["crystals"][x][t]["band_mean_rel_change_mu"]), 5)}
                           for t in TAGS} for x in coh["crystals"]}
    out["coherent"]["max_abs_dlnL08"] = round(max(abs(v[t]["dlnL08"]) for v in out["coherent"].values() for t in TAGS), 4)
    om = json.loads((F("stage3_openmc_check.json")).read_text())
    rd = [abs(c[s]["rel_diff_uncollided"]) for c in om["compounds"] for s in ("1E", "watt")]
    bu = [c[s]["buildup_ratio_T_all_over_T_uncollided"] for c in om["compounds"] for s in ("1E", "watt")]
    out["openmc"] = {"compounds": [c["name"] for c in om["compounds"]], "max_abs_rel_diff_percent": round(100 * max(rd), 2),
                     "buildup_min": round(min(bu), 1), "buildup_max": round(max(bu), 1),
                     "all_pass": all(c[s]["pass"] for c in om["compounds"] for s in ("1E", "watt"))}
    br = json.loads((F("stage3_bragg.json")).read_text())
    lc = [r["lambda_c_A"] for r in br["rows"]]
    out["bragg"] = {"n": br["n"], "lambda_c_min": round(min(lc), 2), "lambda_c_max": round(max(lc), 2),
                    **{k: br[k] for k in ("n_empty_band_lambda_c_ge_10A", "rho_lambda_c_vs_T_pass", "rho_lambda_c_vs_cold_fom", "rho_flux_fraction_vs_cold_fom", "top5_cold_fom")}}
    ctl = json.loads((F("stage3_controls.json")).read_text())
    out["controls"] = {"bi_known_failing": ctl["positive"]["Bi_known_failing"], "all_pass": ctl["all_pass"]}
    from filteratlas.atlas.assemble import SIG_DFT, SIG_MLIP_ADD, SIG_SOFT
    from filteratlas.atlas.fronts import SIG_D
    out["sigma"] = {"dft": SIG_DFT, "mlip_293K": round(float(np.hypot(SIG_DFT, SIG_MLIP_ADD["293K"])), 3),
                    "mlip_77K": round(float(np.hypot(SIG_DFT, SIG_MLIP_ADD["77K"])), 3), "soft_floor": SIG_SOFT,
                    "construction_D": SIG_D, "openmc_pass_percent": round(100 * SIG_DFT, 1)}
    out["nucdata_max_abs_rel_diff_percent"] = round(100 * json.loads((F("stage3_nucdata_check.json")).read_text())
                                                    ["max_abs_rel_diff"], 4)
    st1 = json.loads((F("stage1_stats.json")).read_text())["gate"]
    out["stage1_rho"] = {t: round(st1[t]["primary_L08_vs_proxy"]["rho"], 2) for t in TAGS}
    # --- what the funnel could have missed, and further checks ------------------------------------------------------------------------------------
    mi = json.loads((F("stage3_missed.json")).read_text())
    tbl = pd.read_csv(F("stage3_not_computed_table.csv"))
    out["missed"] = {t: {"debye_front_contenders": mi[t]["debye_front_contenders"],
                         "n_not_excluded_by_proxy_bound_B": mi[t]["n_not_excluded_by_proxy_bound"]["B"],
                         "tierA_not_excluded_by_proxy_bound": names(tbl[(tbl.tier == "A") & tbl[f"not_excluded_by_proxy_bound_{t}"]])}
                     for t in TAGS}
    out["missed"]["n_without_full_model"] = mi["n_without_full_model"]
    out["missed_tol"] = mi["proxy_bound_tolerance"]
    out["missed"]["n_without_debye"] = mi["n_without_debye"]
    out["missed"]["n_with_debye"] = mi["n_with_debye"]
    out["missed"]["density"] = {k: {"not_run_memory": v.get("not_run_memory", {}).get("median_n_A3"),
                                    "sevennet": v.get("sevennet", {}).get("median_n_A3")}
                                for k, v in mi["density_by_supercell_bin"].items()}
    # OpenMC margins SiC vs sapphire
    oc = {c["name"]: c for c in om["compounds"]}
    sic = [k for k in oc if "SiC" in k][0]
    out["openmc"]["margin"] = {s: {"uncollided_SiC": round(oc[sic][s]["A_uncollided"], 3),
                                   "uncollided_sapphire": round(oc["sapphire"][s]["A_uncollided"], 3),
                                   "all_SiC": round(oc[sic][s]["A_all"], 3), "all_sapphire": round(oc["sapphire"][s]["A_all"], 3),
                                   "all_lead_percent": round(100 * (oc[sic][s]["A_all"] / oc["sapphire"][s]["A_all"] - 1), 1),
                                   "uncollided_lead_percent": round(100 * (oc[sic][s]["A_uncollided"] / oc["sapphire"][s]["A_uncollided"] - 1), 1)}
                               for s in ("1E", "watt")}
    # Tier-A ranks of the known filters on both axes, 293 K main set
    d = at[(at["T"] == 293.15) & (at.tier == "A") & at.main_front_eligible].copy()
    d["rank_fast"], d["rank_gamma"] = d.A_fast_1E.rank(ascending=False), d.A_gamma_1MeV.rank(ascending=False)
    out["ranks_293K_tierA"] = {"n": len(d), **{f"{r.formula}_{r.sg_symbol}": [int(r.rank_fast), int(r.rank_gamma)]
                                               for r in d[d.formula.isin(["Al2O3", "Si", "MgO", "SiO2", "SiC", "BeAl2O4",
                                                                          "Y3Al5O12", "Ge3(BiO3)4"])].itertuples()}}
    # the compounds whose log records carry a "machine" field: where they land
    vm_ok = {r["material_id"] for r in ok if "machine" in r}
    b = at[(at["T"] == 293.15) & (at.tier == "B") & at.main_front_eligible & ~at.unstable_dft].copy()
    b["rank_fast"], b["rank_gamma"] = b.A_fast_1E.rank(ascending=False), b.A_gamma_1MeV.rank(ascending=False)
    v1 = one[one.material_id.isin(vm_ok)]
    out["d19"] = {
        "n_ok": len(vm_ok), "n_tierA": int((v1.tier == "A").sum()), "n_tierB": int((v1.tier == "B").sum()),
        "by_funnel_step": v1.funnel_step.value_counts().to_dict(),
        "mlip_unstable_tierA": names(v1[(v1.tier == "A") & v1.mlip_unstable].sort_values("formula")),
        "ranks_293K_tierB": {"n": len(b), **{f"{r.formula}_{r.sg_symbol.replace('/', '_')}": [int(r.rank_fast), int(r.rank_gamma)]
                                             for r in b[b.material_id.isin(vm_ok) & b.formula.isin(["ZrAl3", "BaRuO3"])].itertuples()}},
        "on_front": sorted({f"{x['formula']} {x['sg_symbol']} ({k.replace('|', ' ')})" for k, v in fj["fronts"].items()
                            for x in v["front"] if x["material_id"] in vm_ok}),
        "within_1sigma": sorted({f"{n_} ({k.replace('|', ' ')})" for k in fj["fronts"]
                                 for n_ in names(fm[(fm.key == k) & fm.within_1sigma & fm.material_id.isin(vm_ok)]
                                                 [["material_id", "T"]].merge(at, on=["material_id", "T"]))})}
    out["d19"]["within_1sigma_compounds"] = sorted({x.split(" (")[0] for x in out["d19"]["within_1sigma"]})
    out["d19"]["on_front_compounds"] = sorted({x.split(" (")[0] for x in out["d19"]["on_front"]})
    # per compound, the front keys whose within-1sigma set it joined; and the Tier-A ranks of all these Tier-A
    # compounds in the 293 K main set
    w1, ids = {}, {}
    for k in fj["fronts"]:
        sel = fm[(fm.key == k) & fm.within_1sigma & fm.material_id.isin(vm_ok)][["material_id", "T"]].merge(at, on=["material_id", "T"])
        for r in sel.itertuples():
            c = f"{r.formula} {r.sg_symbol}"
            w1.setdefault(c, set()).add(k.replace("|", " "))
            ids.setdefault(c, set()).add(r.material_id)
    out["d19"]["within_1sigma_by_compound"] = {c: sorted(v) for c, v in sorted(w1.items())}
    out["d19"]["within_1sigma_text"] = [
        f"{c}{f' ({len(ids[c])} entries)' if len(ids[c]) > 1 else ''}: {len(v)} sets ({', '.join(sorted(v))})"
        for c, v in sorted(w1.items())]
    dA = d[d.material_id.isin(vm_ok)].sort_values("rank_fast")
    out["d19"]["tierA_ranks_293K"] = {f"{r.formula} {r.sg_symbol}": [int(r.rank_fast), int(r.rank_gamma)]
                                      for r in dA.itertuples()}
    out["d19"]["tierA_best_A_fast_293K"] = f"{dA.iloc[0].formula} {dA.iloc[0].sg_symbol}" if len(dA) else None
    out["d19"]["n_tierA_main_293K"] = len(dA)
    # composition of the Tier-B within-1sigma set
    w = out["fronts"]["293K_B_main_main"]["within_1sigma_not_front"]
    out["tierB_within_293K_composition"] = {"n": len(w), "SiC": sum(x.startswith("SiC ") for x in w),
                                            "Ce": sum("Ce" in x.split()[0] for x in w),
                                            "Zr_no_Ce": sum("Zr" in x.split()[0] and "Ce" not in x.split()[0] for x in w)}
    # flags of individually discussed compounds
    flagged = {}
    for nm in ("Ge(Bi3O5)4", "BeAl2O4", "Si(Bi3O5)4", "C", "ZrO2", "ZrBi", "Zr2PbC", "CeS"):
        r = at[(at.formula == nm) & (at["T"] > 100) & at.vdos_source.ne("none")]
        if len(r):
            r = r.iloc[0]
            flagged[nm] = [c for c in ("untested_heavy", "flag_Be", "flag_Bi", "flag_Pb", "interp_only", "tierA_phase_mismatch",
                                       "not_off_bragg_single_crystal", "unscreenable") if bool(r[c])]
    out["flags_named"] = flagged
    # SiC polytypes and the unscreenable outlier
    sc = at[(at.formula == "SiC") & (at["T"] > 100) & at.vdos_source.ne("none")]
    out["sic_polytypes_L08_293K"] = {"n": len(sc), "min": round(float(sc.L08.min()), 2), "max": round(float(sc.L08.max()), 2),
                                     "screened_min": round(float(sc[~sc.unscreenable].L08.min()), 2),
                                     "unscreenable": [f"{r.sg_symbol} {r.L08:.2f} cm" for r in sc[sc.unscreenable].itertuples()],
                                     "dln_sigma": round(float(abs(np.log(sc[sc.unscreenable].L08.min() / sc[~sc.unscreenable].L08.min()))
                                                              / out["sigma"]["dft"]), 1)}
    out["bi_proxy_L08"] = round(float(at[(at.formula == "Bi") & (at.tier == "A")].proxy_L08.iloc[0]), 0)
    out["kill2_n_evaluated"] = len(json.loads((F("stage2_exfor.json")).read_text())["kill2"]["evaluated_crystals"])
    g = json.loads((F("stage3_laptop_gpu_check.json")).read_text())
    out["gpu_check"] = {"n": len(g["rows"]), "atoms_min": min(r["supercell_atoms"] for r in g["rows"]),
                        "atoms_max": max(r["supercell_atoms"] for r in g["rows"])}
    out["coherent"]["max_abs_dsigma_b_by_qcut"] = {
        q: round(max(abs(a - b) for x in coh["crystals"].values() for t in TAGS
                     for a, b in zip(x[t][q]["coh_exact_b"], x[t][q]["coh_incapprox_b"])), 4) for q in ("8.0", "11.0", "14.0")}
    out["referee_query"]["not_run_memory_by_debye"] = list({r["compound"]: r for r in
                                                           out["referee_query"]["not_run_memory_by_debye"]}.values())
    out["tierB_unscreenable"] = fj["tierB_exclude_unscreenable"]
    out["ladder_extra"] = {t: {"proxy_vs_debye_rho": fj["ladder"][f"{t}|all"]["proxy_vs_debye"]["rho"],
                               "proxy_vs_full_on_debye_subset_rho": fj["ladder"][f"{t}|all"]["proxy_vs_full_on_debye_subset"]["rho"],
                               "kendall_debye_vs_full_top10": fj["ladder"][f"{t}|all"]["kendall_debye_vs_full_top10_full"],
                               "debye_level_front": fj["ladder"][f"{t}|all"]["debye_level_front"],
                               "full_front_on_debye_subset": fj["ladder"][f"{t}|all"]["full_front_on_debye_subset"]} for t in TAGS}
    # budget outcome per priority arm and the extra arms (large-cell garnets, Tier-A compounds without a usable MP DOS)
    arm = pd.DataFrame(sset)[["material_id", "arm"]].merge(one[["material_id", "funnel_step"]], on="material_id")
    bud = arm[arm.arm.isin(["atlas_a", "atlas_b", "atlas_c"])]
    out["budget"] = {"n_slots": len(bud), "by_outcome": {k: int((bud.funnel_step == k).sum()) for k in
                                                         sorted({"sevennet", "mlip_failed", "not_run_memory"} | set(bud.funnel_step))},
                     "extras_D15": int((arm.arm == "atlas_a_large_D15").sum()), "extras_D16": int((arm.arm == "atlas_a_D16").sum())}
    cpu = [r for r in log if r.get("device") == "cpu"]
    out["compute"]["max_atoms_cpu_job"] = max(r["supercell_atoms"] for r in cpu)
    full1 = one[one.vdos_source.ne("none")]
    out["domain_overlap"] = {"soft": int(full1.soft_lattice.sum()), "heavy_soft_subset_of_soft": int(full1.heavy_soft.sum()),
                             "molecular_not_soft": int((full1.molecular & ~full1.soft_lattice).sum()),
                             "outside_domain_total": int(full1.outside_domain.sum())}
    (F("stage3_text_numbers.json")).write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
