"""Derived numbers of the validation and uncertainty analysis (step 06), so that no quoted number is hand-typed.

Every value is computed here from the JSON and parquet files written by steps 05 (validation) and 06 (uncertainty).
-> the validation text-numbers JSON in data/derived
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402

THZ_TO_MEV = 4.135667696


def main():
    ex = json.loads((F("stage2_exfor.json")).read_text())
    res = pd.DataFrame(json.loads((F("stage2_residuals.json")).read_text())["datasets"])
    bn = pd.read_parquet(F("stage2_exfor_binned.parquet"))
    mlip = json.loads((F("stage2_mlip.json")).read_text())["models"]
    unc = json.loads((F("stage2_uncertainty_model.json")).read_text())
    orb = json.loads((F("stage2_orb_diagnostic.json")).read_text())
    aud = {(a["formula"], a["subentry"]): a for a in ex["audit"]}
    r = lambda f, s, m: float(res[(res.formula == f) & (res.subentry == s) & (res.model == m)].r_d.iloc[0])  # noqa: E731
    out = {}
    out["si_harvey_13733004_r_d"] = {m: round(r("Si", 13733004, m), 3)
                                     for m in ("mp_dft", "mlip_sevennet_mf_ompa", "curated")}
    lo = bn[(bn.formula == "Al2O3") & bn.subentry.isin([13761007, 14302002]) & (bn.model == "mp_dft") & (bn.E_c < 0.01)]
    out["sapphire_harvey_mildner_mean_ln_ratio_below_10meV"] = round(float(lo.lr.mean()), 3)
    pg = res[(res.formula == "C") & (res.model == "mp_dft")]
    out["pyrolytic_graphite_mean_abs_bias_factor"] = round(float(np.exp(pg["mean"].abs().mean())), 2)
    out["pyrolytic_graphite_free_atom_ratio"] = [aud[("C", s)]["free_atom_ratio"] for s in (13761002, 13761003)]
    out["sapphire_spread_percent"] = round(100 * (np.exp(ex["crystals"]["Al2O3"]["S_c"]) - 1))
    out["si_spread_percent"] = round(100 * (np.exp(ex["crystals"]["Si"]["S_c"]) - 1))
    out["cantargi_brugger_max_r_d_mp_dft"] = round(max(r("Al2O3", 31743002, "mp_dft"),
                                                       *(r("Si", s, "mp_dft") for s in (10943002, 10943003, 10943004))), 3)
    out["mace_softening_percent"] = round(100 * (1 - mlip["mace_mp0_medium"]["softening_median"]))
    out["sevennet_softening_percent"] = round(100 * (1 - mlip["sevennet_mf_ompa"]["softening_median"]))
    out["debye300_sapphire_factor"] = round(float(np.exp(ex["crystals"]["Al2O3"]["per_model"]["debye300"]["R_c"])), 1)
    out["model_vs_data_band_rms_DFT"] = unc["band_DFT"]["rms"]
    out["model_vs_data_over_mlip_vs_dft"] = round(unc["band_DFT"]["rms"] / unc["mlip_vs_dft_dlnL08"]["293K"]["p68_abs_dlog"], 1)
    # ORB diagnostic: compounds whose failure depends on the displacement amplitude (base fails, disp 0.03 A differs
    # by more than a factor 2 in the commensurate minimum or passes)
    runs = pd.DataFrame(orb["runs"])
    piv = runs.pivot_table(index="formula", columns="variant", values="min_freq_comm_THz")
    dep = piv[(piv["disp03"] > -0.5) | (piv["disp03"] / piv["base"] < 0.5)]
    out["orb_diag_n_compounds"] = int(len(piv))
    out["orb_diag_displacement_sensitive"] = int(len(dep))
    out["orb_diag_displacement_sensitive_list"] = sorted(dep.index)
    out["orb_diag_rescued_by_disp03"] = sorted(piv[piv["disp03"] > -0.5].index)
    out["zro2_duplicate_pair"] = ex["crystals"]["ZrO2"]["spread_pairs"]
    out["min_pair_spread_excluding_ZrO2"] = min(p["rms"] for f, c in ex["crystals"].items() if f != "ZrO2"
                                                for p in c["spread_pairs"])
    out["sevennet_stratum_max_over_min"] = round(max(v for v in mlip["sevennet_mf_ompa"]["p68_by_stratum_293K"].values() if v)
                                                 / min(v for v in mlip["sevennet_mf_ompa"]["p68_by_stratum_293K"].values() if v), 2)
    # robustness variants and the Pb/Bi quasi-harmonic run
    rob = json.loads((F("stage2_d11_robustness.json")).read_text())
    V = rob["variants"]
    short = {"as_reported": "D11 as reported (q20, all-form S_c)",
             "sc_floor": "D11, single-crystal S_c (q20)",
             "revised_q10": "D11 + reclassification, single-crystal S_c, q10 (C kept failing)",
             "revised": "D11 + reclassification, single-crystal S_c, q20 (C kept failing)",
             "revised_q50": "D11 + reclassification, single-crystal S_c, q50 (C kept failing)",
             "no_floor": "D11 + reclassification, flat 0.15 (no floor), q20",
             "no_floor_no_reclass": "D11, flat 0.15 (no floor), q20, no reclassification",
             "without_13761": "D11 + reclassification, without Harvey & Smith entry 13761 (C kept failing)",
             "T_minus20": "D11 + reclassification, assumed temperatures -20 K",
             "T_plus20": "D11 + reclassification, assumed temperatures +20 K"}
    out["d11_variants"] = {k: {"n": V[v]["n_failing"], "of": V[v]["n_evaluable"], "failing": V[v]["failing"],
                               "fires": V[v]["fires"]} for k, v in short.items()}
    lr = V["leave-one-dataset-out (reclassified, single-crystal S_c, q20): range"]
    out["d11_loo_range"] = [lr["min"], lr["max"]]
    out["d11_single_dataset_crystals"] = rob["single_dataset_crystals"]
    out["d11_entry_13761_only_crystals"] = rob["entry_13761_only_crystals"]
    out["reclassified"] = [x[0] for x in rob["reclassified"]]
    rd = pd.DataFrame(rob["datasets"])
    out["sapphire_signed_bias_q20"] = [round(x, 3) for x in rd[rd.formula == "Al2O3"].sort_values("subentry").bias_atlas_q20]
    out["sc_temp_assumed"] = [int(rd.temp_assumed.sum()), len(rd)]
    sp = pd.DataFrame(rob["polycrystal_cutoff_split"])
    sp = sp[sp.model == "mp_dft"]
    for f in ("BeO", "C", "Cu"):
        g = sp[sp.formula == f]
        out[f"cutoff_split_{f}"] = {"below_mean_ln": [float(g.mean_ln_below.min()), float(g.mean_ln_below.max())],
                                    "above_r_d": [float(g.r_d_above.min()), float(g.r_d_above.max())]}
    g = sp[sp.formula == "Ge"]
    out["cutoff_split_Ge_above_r_d"] = [float(g.r_d_above.min()), float(g.r_d_above.max())]
    qha = json.loads((F("stage2_pbbi_qha.json")).read_text())
    out["pb_static_volume_A3_per_atom"] = round(qha["Pb"]["V0_static_A3"], 2)
    out["pb_static_fcc_a_A"] = round((4 * qha["Pb"]["V0_static_A3"]) ** (1 / 3), 3)   # 1-atom primitive cell
    out["bi_qha_min_mesh_freq_THz"] = [round(min(qha["Bi"]["min_freq_mesh_THz"]), 3), round(max(qha["Bi"]["min_freq_mesh_THz"]), 3)]
    pj = json.loads((F("stage2_pbbi_diag.json")).read_text())["crystals"]
    out["pb_scale_closing_range"] = [min(pj["Pb"]["scale_closes_gap"]), max(pj["Pb"]["scale_closes_gap"])]
    out["pg_poly_model_r_d_curated"] = [round(float(x), 2) for x in rd[rd.formula == "C"].r_poly_curated.agg(["min", "max"])]
    out["pb_softening_percent_range"] = [round(100 * (1 - max(pj["Pb"]["scale_closes_gap"]))),
                                         round(100 * (1 - min(pj["Pb"]["scale_closes_gap"])))]
    spb = sp[sp.formula.isin(["BeO", "C", "Cu"])]
    out["cutoff_split_n_bins_below"] = [int(spb.n_below.min()), int(spb.n_below.max())]
    out["n_single_dataset_crystals"] = len(rob["single_dataset_crystals"])
    # per-crystal table after reclassification, sapphire envelope dependence, cut-off split, MLIP sensitivity, GPU run
    recl = {(int(x[0]), float(x[1])) for x in rob["reclassified"]}
    keep = rd[[(int(a_), float(b_)) not in recl for a_, b_ in zip(rd.subentry, rd["T"])]]
    fl = rob["spread_floors"]["q20"]["S_c_single_crystal_without_reclassified"]
    tab = {}
    for f in sorted(keep.formula.unique()):
        thr = max(0.15, fl[f]) if fl.get(f) is not None else 0.15
        R = float(np.sqrt(np.mean(keep[keep.formula == f].r_atlas_q20 ** 2)))
        tab[f] = {"threshold": round(thr, 3), "R_c_atlas": round(R, 3), "fails": "yes" if R > thr else "no",
                  "n_datasets": int((keep.formula == f).sum())}
    out["revised_table"] = tab
    bm = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    x = bm[(bm.subentry == 13761005) & (bm.model == "mlip_sevennet_mf_ompa") & (bm.envelope_q == 20)]
    out["mgo_005_sevennet_r_d"] = round(float(x.r_d.iloc[0]), 3)
    clean = [("Ge", None), ("Cu", None), ("SiO2", None), ("Zn", None), ("Al", None)]
    vals = [tab[f]["R_c_atlas"] for f, _ in clean]
    vals += [round(float(rd[(rd.subentry == s_)].r_atlas_q20.iloc[0]), 3) for s_ in (10943002, 10943003, 10943004, 31743002)]
    out["clean_single_crystal_r_range"] = [min(vals), max(vals)]
    sa = rd[rd.formula == "Al2O3"].sort_values("subentry")
    out["sapphire_bias_by_envelope"] = {q: [round(v, 3) for v in sa[f"bias_atlas_{q}"]] for q in ("q10", "q20", "q50")}
    out["sapphire_subentries"] = [int(v) for v in sa.subentry]
    spa = pd.DataFrame(rob["polycrystal_cutoff_split"])
    out["slope_below_BeO_mp_dft"] = float(spa[(spa.formula == "BeO") & (spa.model == "mp_dft")].slope_excess_vs_ln_lambda_below.iloc[0])
    cs = spa[spa.formula == "C"].groupby("model").mean_ln_below.agg(["min", "max"]).round(2)
    out["graphite_below_mean_ln_by_model"] = {m: [float(r["min"]), float(r["max"])] for m, r in cs.iterrows()}
    unc = json.loads((F("stage2_uncertainty_model.json")).read_text())
    out["unc_band_dft_crystals"] = sorted(unc["band_per_crystal_DFT"])
    sens = json.loads((F("stage2_mlip_sensitivity.json")).read_text())
    so = sens["orb_disp030"]
    out["orb_disp030_fail_percent"] = round(100 * so["random_failed"] / so["random_n"], 1)
    out["orb_disp030_still_excluded"] = "yes" if so["random_failed"] / so["random_n"] > 0.10 else "no"
    out["orb_ref_fail_percent"] = round(100 * so["random_failed_ref"] / so["random_n"], 1)
    out["gpu_speedup_per_compound"] = round(1 / sens["sevennet_cuda"]["wall_time"]["median_ratio_new_over_ref"], 1)
    log = [json.loads(x) for x in open((D("logs") / "mlip_sevennet_mf_ompa_cuda.jsonl"), encoding="utf-8") if x.strip()]
    last, ever_err = {}, set()
    for x in log:
        last[x["material_id"]] = x
        if x.get("status") == "error":
            ever_err.add(x["material_id"])
    out["gpu_oom_unresolved"] = sorted(last[m]["formula"] for m in ever_err if last[m].get("status") != "ok")
    out["gpu_oom_rescued_by_single_process"] = sorted(last[m]["formula"] for m in ever_err if last[m].get("status") == "ok")
    # memory held by our own process at the failure, from the PyTorch message of the final (single-process) attempts
    held = [float(re.search(r"this process has ([\d.]+) GiB", last[m]["error"]).group(1))
            for m in ever_err if last[m].get("status") != "ok"]
    out["gpu_oom_unresolved_own_process_GiB"] = [min(held), max(held)]
    out["gpu_total_GiB"] = float(re.search(r"total capacity of ([\d.]+) GiB", next(last[m]["error"] for m in sorted(ever_err) if last[m].get("status") != "ok")).group(1))
    (F("stage2_text_numbers.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
