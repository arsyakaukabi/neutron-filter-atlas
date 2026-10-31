"""Summary numbers of the shortlist, with template-friendly keys.

Reads the shortlist derived files (plus the atlas, the candidate table and the uncertainty model). Compound keys are slugs: formula and space group with every character outside
[A-Za-z0-9] replaced by "_" (e.g. "SiC P6_3mc" -> "SiC_P6_3mc"); lists are rendered with |join.
-> derived/stage4_text_numbers.json
"""
import json
import re

import numpy as np
import pandas as pd

from filteratlas.common.paths import F

KEY = {"sic6h": "SiC P6_3mc [alhn]", "sapphire": "Al2O3 R-3c [abrz]", "chrysoberyl": "BeAl2O4 Pnma [aeon]",
       "zro2": "ZrO2 P2_1/c [aefy]", "bi12geo20": "Ge(Bi3O5)4 I23 [bioe]", "bgo": "Ge3(BiO3)4 I-43d [biwe]",
       "yag": "Y3Al5O12 Ia-3d [aeni]", "mgo": "MgO Fm-3m [abwr]", "si": "Si Fd-3m [aaft]", "quartz": "SiO2 P3_221 [akgo]",
       "bi": "Bi R-3m [bigm]"}


def slug(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name.split(" [")[0]).strip("_")


def short(name):
    return name.split(" [")[0]


def fmt_p(q, p=0.90):
    """Two decimals, three within 0.01 of the statement threshold."""
    return f"{q:.3f}" if min(abs(q - p), abs(q - (1 - p))) < 0.01 else f"{q:.2f}"


def fmt_bq(x):
    return "< 1e-6" if x is not None and x < 1e-6 else (f"{x:.3g}" if x is not None else None)


def main():
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    # SiC entries are named with their polytype or site count, never by space group alone
    disp_ = {f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]":
             f"{r.formula} {r.sg_symbol}" + (f" ({r.polytype})" if r.polytype else "") for r in sl.itertuples()}

    def disp(name):
        return disp_.get(name, short(name))
    rk = json.loads(F("stage4_rank_probabilities.json").read_text())
    tr = json.loads(F("stage4_transport.json").read_text())
    dp = json.loads(F("stage4_dips.json").read_text())
    ts = json.loads(F("stage4_temperature.json").read_text())
    ds = json.loads(F("stage4_dossier.json").read_text())
    ct = json.loads(F("stage4_controls.json").read_text())
    m = rk["main"]
    out = {"shortlist": {"n": len(sl), "n_tierA": int((sl.tier == "A").sum()), "n_tierB": int((sl.tier == "B").sum()),
                         "n_SiC": int((sl.formula == "SiC").sum()),
                         "tierA": [f"{r.formula} {r.sg_symbol}" for r in sl[sl.tier == "A"].itertuples()],
                         "tierB_front": [f"{r.formula} {r.sg_symbol}" for r in sl[(sl.tier == "B") & sl.reasons.str.contains("front")].itertuples()]}}
    out["mc"] = {"n_draws": rk["settings"]["n_draws"], "n_compounds": rk["n_compounds"]}
    k3 = {}
    for sc in ("main", "independent_T", "common", "sigma_D"):
        k = rk[sc]["kill3"]
        k3[sc] = {"P_K_293K": k["P_K"]["293K"], "P_K_77K": k["P_K"]["77K"], "P_R": k["P_R"], "tau_median": k["tau_median"],
                  "tau_lo": k["tau_p2.5_p97.5"][0], "tau_hi": k["tau_p2.5_p97.5"][1], "fires": k["fires"]}
    out["kill3"] = k3
    pf = {}
    for t in ("293K", "77K"):
        pf[t] = [f"{disp(n)} {p:.2f}" for n, p in m[t]["p_front"]["main_A"].items() if p >= 0.05]
        pf[f"{t}_B"] = [f"{disp(n)} {p:.2f}" for n, p in list(m[t]["p_front"]["main_B"].items())[:8]]
    out["p_front"] = pf
    beats = {}
    for k, n in KEY.items():
        if k == "sapphire":
            continue
        for t in ("293K", "77K"):
            for ax, a in (("fast", "A_fast_1E"), ("gamma", "A_gamma_1MeV")):
                beats[f"{k}_{ax}_{t}"] = m[t]["p_beats_sapphire"][a].get(n)
    out["p_beats_sapphire"] = beats
    st = m["statements"]
    main_axes = [s for s in st if s["axis"] in ("A_fast_1E", "A_gamma_1MeV")]
    out["statements"] = {"n_all_axes": len(st), "n_main_axes": len(main_axes),
                         "beats_sapphire_main_axes": [f"{disp(s['compound'])} ({s['axis'].replace('A_', '')}, {s['T']}, P = {fmt_p(s['P'])})"
                                                      for s in main_axes if s["statement"] == "beats sapphire"],
                         "tierA_nonknown_beats": [f"{short(s['compound'])} ({s['axis']}, {s['T']})" for s in main_axes
                                                  if s["statement"] == "beats sapphire" and s["compound"] in
                                                  {f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]" for r in sl.itertuples()
                                                   if r.tier == "A" and not r.known_filter}]}
    out["sensitivity_changes"] = {sc: {"lost": len(v["lost"]), "gained": len(v["gained"])}
                                  for sc, v in rk["sensitivity_changes"].items()}
    fr = tr["framing"]
    out["framing"] = {"framing": fr["framing"], "framing_ratio_model": tr["sensitivity_ratio_model"]["framing"]["framing"],
                      "n_qualifying_before_transport": len(fr["qualifying_before_transport"]),
                      "transport_induced": [f"{disp(s['compound'])} ({s['axis']}, {s['T']}, P = {s['P']:.2f}; flags: {s['flags']})"
                                            for s in fr["transport_induced_not_counted"]],
                      "transport_induced_ratio_model": [f"{disp(s['compound'])} ({s['axis']}, {s['T']}, P = {s['P']:.2f}; flags: {s['flags']})"
                                                        for s in tr["sensitivity_ratio_model"]["framing"]["transport_induced_not_counted"]],
                      "candidates": [f"{disp(c)} (P on front {v['293K']:.2f} at 293 K, {v['77K']:.2f} at 77 K)"
                                     for c, v in fr["candidates"].items()]}
    # the framing under the common-component sensitivity run (fixed before the computation)
    cc = tr["sensitivity_common_component"]
    fmt_s = lambda s: f"{disp(s['compound'])} ({s['axis']}, {s['T']}, P = {fmt_p(s['P'])})"  # noqa: E731
    out["framing"]["common"] = {
        "qualifying_before_transport": [fmt_s(s) for s in cc["shift"]["framing"]["qualifying_before_transport"]],
        "framing_shift": cc["shift"]["framing"]["framing"], "framing_ratio": cc["ratio"]["framing"]["framing"],
        "discovery_shift": [fmt_s(s) for s in cc["shift"]["framing"]["discovery_statements"]],
        "discovery_ratio": [fmt_s(s) for s in cc["ratio"]["framing"]["discovery_statements"]],
        "chrysoberyl_gamma_293K_common": rk["common"]["293K"]["p_beats_sapphire"]["A_gamma_1MeV"].get(KEY["chrysoberyl"])}
    # untransported statements on these compounds, with their status after transport
    ptr = tr["p_beats_sapphire_transported"]
    p_ = rk["settings"]["p_statement"]
    ba = []
    for s in st:
        if s["axis"] not in ("A_fast_1E", "A_gamma_1MeV") or short(s["compound"]) not in (
                "Y3Al5O12 Ia-3d", "Ge3(BiO3)4 I-43d", "Be2C Fm-3m", "BeO P6_3mc"):
            continue
        q = ptr[s["T"]][s["axis"]].get(s["compound"])
        after = ("beats sapphire" if q > p_ else "sapphire beats it" if q < 1 - p_ else "not distinguishable")
        ba.append(f"{disp(s['compound'])}, {s['axis']}, {s['T']}: {s['statement']} before transport (P = {fmt_p(s['P'])}); "
                  f"after transport {after} (P = {fmt_p(q)})")
    out["statements"]["before_after_transport"] = ba
    # sapphire conditioned on its measured residual (scenario added after inspection of the results)
    so = rk.get("sapphire_offset")
    if so:
        tso = tr.get("sensitivity_sapphire_offset", {})
        out["sapphire_offset"] = {
            "offset": round(so["offset_lnL08"], 3), "P_K_293K": so["kill3"]["P_K"]["293K"], "P_K_77K": so["kill3"]["P_K"]["77K"],
            "fires": so["kill3"]["fires"],
            **{f"{k}_{ax}_{t}": so[t]["p_beats_sapphire"][a].get(KEY[k]) for k in ("sic6h", "chrysoberyl")
               for t in ("293K", "77K") for ax, a in (("fast", "A_fast_1E"), ("gamma", "A_gamma_1MeV"))},
            "p_front_293K": [f"{disp(n)} {p:.2f}" for n, p in so["293K"]["p_front"]["main_A"].items() if p >= 0.05],
            "p_front_77K": [f"{disp(n)} {p:.2f}" for n, p in so["77K"]["p_front"]["main_A"].items() if p >= 0.05],
            "framing_after_transport": tso.get("framing", {}).get("framing"),
            "n_qualifying_before_transport": len(tso.get("framing", {}).get("qualifying_before_transport", []))}
    unc = json.loads(F("stage2_uncertainty_model.json").read_text())["band_DFT"]
    m_, s_ = unc["common_mode_mean"], unc["compound_specific_sd"]
    out["common_component_from_validation"] = {"common_mode_mean": m_, "compound_specific_sd": s_,
                                               "c_hat": round(m_ ** 2 / (m_ ** 2 + s_ ** 2), 2)}
    mx = tr["max_abs_check_rel_diff"]
    out["transport"] = {"n_runs": tr["n_runs"], "check_photon_1MeV_pct": round(100 * mx["g1MeV"], 2),
                        "check_photon_6MeV_pct": round(100 * mx["g6MeV"], 2), "check_thermal_pct": round(100 * mx["thermal"], 2),
                        "check_1E_pct": round(100 * mx["1E"], 2), "check_watt_pct": round(100 * mx["watt"], 2),
                        "n_failed_preregistered": tr["check_preregistered_setting"]["n_failed"],
                        "check_1E_preregistered_pct": round(100 * tr["check_preregistered_setting"]["max_abs_check_rel_diff"]["1E"], 2),
                        "n_failed_coh_elas": sum(v for k, v in tr["checks_failed_by_cause"].items() if k.startswith("coherent elastic")),
                        "n_failed_ncrystal_model": sum(v for k, v in tr["checks_failed_by_cause"].items() if k.startswith("NCrystal's model")),
                        "n_failed_be_n2n": sum(v for k, v in tr["checks_failed_by_cause"].items() if "(n,2n)" in k),
                        "n_unexplained": len(tr["checks_unexplained"]),
                        "n_runs_1E_all": tr["n_runs_1E_all_components"],
                        "free_atom_max_pct": round(100 * max(abs(d["rel_diff"]) for d in tr["free_atom_diagnostic"] if d["case"].endswith("_1E")), 2),
                        "thermalised_1E_max_pct": round(100 * tr["thermalised_fraction_range"]["1E"][1], 1),
                        "thermalised_watt_max_pct": round(100 * tr["thermalised_fraction_range"]["watt"][1], 1),
                        "front_A_293K_ratio_model": [disp(x) for x in tr["shortlist_fronts"]["293K_A"]["transported_ratio_model"]],
                        "front_A_77K_ratio_model": [disp(x) for x in tr["shortlist_fronts"]["77K_A"]["transported_ratio_model"]],
                        "correction_fast_min": round(tr["correction_range"]["A_fast_1E"][0], 2),
                        "correction_fast_max": round(tr["correction_range"]["A_fast_1E"][1], 2),
                        "correction_gamma_min": round(tr["correction_range"]["A_gamma_1MeV"][0], 2),
                        "correction_gamma_max": round(tr["correction_range"]["A_gamma_1MeV"][1], 2),
                        "front_A_293K_transported": [disp(x) for x in tr["shortlist_fronts"]["293K_A"]["transported"]],
                        "front_A_77K_transported": [disp(x) for x in tr["shortlist_fronts"]["77K_A"]["transported"]],
                        "front_A_293K_uncollided": [disp(x) for x in tr["shortlist_fronts"]["293K_A"]["uncollided"]],
                        "front_A_77K_uncollided": [disp(x) for x in tr["shortlist_fronts"]["77K_A"]["uncollided"]],
                        "statements_lost": len(tr["statements_changed_by_transport"]["lost"]),
                        "statements_gained": len(tr["statements_changed_by_transport"]["gained"])}
    for k in ("sapphire", "sic6h", "chrysoberyl", "bi12geo20"):
        for t in ("293K", "77K"):
            e = tr["at_L08"][KEY[k]][t]
            out["transport"][f"{k}_{t}"] = {"fast_unc": e["A_fast_1E"]["uncollided"], "fast_tr": e["A_fast_1E"]["transport"],
                                            "gamma_unc": e["A_gamma_1MeV"]["uncollided"], "gamma_tr": e["A_gamma_1MeV"]["transport"],
                                            "T_th_all": round(e["T_th_all_thermal"], 3)}
    dd = {}
    for n, r in dp["compounds"].items():
        e = r["mosaic"]["0.3deg"]
        dd[slug(n)] = {"best_0.3": round(e["blocked_fraction"]["best"], 2), "median_0.3": round(e["blocked_fraction"]["median"], 2),
                       "best_0.1": round(r["mosaic"]["0.1deg"]["blocked_fraction"]["best"], 2),
                       "best_1.0": round(r["mosaic"]["1.0deg"]["blocked_fraction"]["best"], 2),
                       "best_low_index_0.3": round(e["low_index_directions"]["best_blocked_fraction"], 2),
                       "best_low_index_dir": e["low_index_directions"]["best"],
                       "median_low_index_0.3": round(e["low_index_directions"]["median_blocked_fraction"], 2),
                       "best_all_0.3": round(e["best_blocked_fraction_grid_and_low_index"], 2),
                       "best_low_index_0.1": round(r["mosaic"]["0.1deg"]["low_index_directions"]["best_blocked_fraction"], 2)}
    sap = dp["compounds"][KEY["sapphire"]]["mosaic"]["0.3deg"]["literature_directions"]["beam || c [001]"]
    si = dp["compounds"][KEY["si"]]["mosaic"]["0.3deg"]["literature_directions"]
    out["dips"] = {"compounds": dd, "sapphire_c_blocked": round(sap["blocked_fraction"], 2),
                   "si_100_blocked": round(si["<100>"]["blocked_fraction"], 2), "si_111_blocked": round(si["<111>"]["blocked_fraction"], 2),
                   "si_110_blocked": round(si["<110>"]["blocked_fraction"], 2),
                   "b5_matched_0002": round(sap["mildner_b5_check"]["tol_0.0002"]["fraction_of_model_dips_matched"], 2),
                   "b5_chance_0002": round(sap["mildner_b5_check"]["tol_0.0002"]["chance_level"], 3),
                   "b5_n_dips": sap["mildner_b5_check"]["n_model_dips"], "n_directions": dp["settings"]["n_directions"],
                   "b5_n_unmatched_0002": round(sap["mildner_b5_check"]["n_model_dips"]
                                                * (1 - sap["mildner_b5_check"]["tol_0.0002"]["fraction_of_model_dips_matched"])),
                   "n_low_index": dp["compounds"][KEY["sapphire"]]["mosaic"]["0.3deg"]["low_index_directions"]["n"]}
    # sapphire tilt scan about c (Born et al. / Nieman et al. as quoted by Mildner et al. 1993)
    ts_ = dp["compounds"][KEY["sapphire"]]["mosaic"]["0.3deg"]["tilt_scan_about_c"]
    out["dips"]["tilt"] = {d: {f"{r['tilt_deg']:g}": [round(r["blocked_fraction_min"], 2), round(r["median"], 2), round(r["max"], 2)]
                               for r in rows} for d, rows in ts_.items()}
    # is the compound ranking (grid-best) the same at the three mosaics?
    order = {m_: [n for n, _ in sorted(dd.items(), key=lambda kv: kv[1][f"best_{m_}"])] for m_ in ("0.1", "0.3", "1.0")}
    out["dips"]["ranking_same_all_mosaics"] = order["0.1"] == order["0.3"] == order["1.0"]
    names_ = {slug(k): short(k) for k in dp["compounds"]}

    def swaps(a, b):
        """Adjacent pairs that trade places between two rankings, each pair once."""
        return [f"{names_[a[i]]} and {names_[a[i + 1]]}" for i in range(len(a) - 1)
                if a[i] != b[i] and a[i] == b[i + 1] and a[i + 1] == b[i]]
    out["dips"]["ranking_swaps_vs_0.3"] = {m_: swaps(order["0.3"], order[m_]) for m_ in ("0.1", "1.0")}
    tt = {}
    for n, r in ts["compounds"].items():
        s = r["series"]
        tt[slug(n)] = {"L08_293": round(s[0]["L08"], 2), "L08_77": round(s[-1]["L08"], 2), "gain": round(r["L08_gain_293_to_77"], 2),
                       "best_interval": f"{r['largest_gain_per_K']['from_K']:.0f}-{r['largest_gain_per_K']['to_K']:.0f} K",
                       "best_relative_interval": f"{r['largest_relative_gain_per_K']['from_K']:.0f}-{r['largest_relative_gain_per_K']['to_K']:.0f} K",
                       "best_margin_pct": round(100 * (lambda g: g[0] / g[1] - 1)(sorted((x["dL08_per_K"] for x in r["intervals"]), reverse=True)), 1)}
    out["temperature"] = tt
    dz = {}
    for n, r in ds["compounds"].items():
        a = r.get("activation", {})
        dz[slug(n)] = {"melting_point": r.get("melting_point", {}).get("value", "not found"),
                       "largest": r.get("largest_stated_crystal", {}).get("value", "not found"),
                       "growth": r.get("growth_method", {}).get("value", "not found"),
                       "water": r.get("water", {}).get("value", "not found"),
                       "price_per_cm3": (r.get("price", {}) or {}).get("price_per_cm3") or "not found",
                       "toxicity": "; ".join(r["toxicity"]),
                       "activity_end_Bq": a.get("0d_after_end", {}).get("total_Bq"),
                       "activity_1d_Bq": a.get("1d_after_end", {}).get("total_Bq"),
                       "activity_7d_Bq": a.get("7d_after_end", {}).get("total_Bq"),
                       "activity_1d_text": fmt_bq(a.get("1d_after_end", {}).get("total_Bq")),
                       "activity_7d_text": fmt_bq(a.get("7d_after_end", {}).get("total_Bq")),
                       "dominant_7d_shares": {x["nuclide"]: x["share"] for x in a.get("7d_after_end", {}).get("dominant", [])},
                       "dominant_7d": ", ".join(x["nuclide"] for x in a.get("7d_after_end", {}).get("dominant", [])) or "none",
                       "measurement_search": r.get("measurement_search", {}).get("value", "")}
    out["dossier"] = dz
    out["controls"] = {"sapphire_pct_293": ct["positive"]["sapphire"]["293K"]["L08_percentile_tierA"],
                       "sapphire_pct_77": ct["positive"]["sapphire"]["77K"]["L08_percentile_tierA"],
                       "si_pct_293": ct["positive"]["Si"]["293K"]["L08_percentile_tierA"],
                       "si_pct_77": ct["positive"]["Si"]["77K"]["L08_percentile_tierA"],
                       "quartz_pct_293": ct["positive"]["quartz"]["293K"]["L08_percentile_tierA"],
                       "quartz_pct_77": ct["positive"]["quartz"]["77K"]["L08_percentile_tierA"],
                       "top_quartile_pass": ct["pass"]["sapphire_Si_top_quartile_both_T"],
                       "quartz_pass": ct["pass"]["quartz_intermediate_both_T"],
                       "bi_overstatement": ct["Bi_known_failing"]["L08_overstatement_percent_stage2"],
                       "n_tierA": ct["positive"]["sapphire"]["293K"]["n_tierA"],
                       "nieman": ct["nieman_ratio"], "nieman_pass": ct["pass"]["sapphire_largest_nieman_ratio_both_T"]}
    out["anchors"] = ct["anchors"]
    # SiC polytypes in the atlas
    at = pd.read_parquet(F("atlas.parquet"))
    cand = pd.read_parquet(F("stage0_candidates.parquet"))
    sic = at[(at.formula == "SiC") & (at["T"] > 100) & at.vdos_source.ne("none")].merge(sl[["material_id", "polytype"]], on="material_id")
    is15 = sic.polytype == "15R"
    out["sic"] = {"has_4H_entry": bool(((cand.formula == "SiC") & (cand.sg_symbol == "P6_3mc") & (cand.nsites == 8)).any()),
                  "n_full_model": len(sic), "n_excl_15R": int((~is15).sum()),
                  "L08_min_excl_15R": round(float(sic[~is15].L08.min()), 2),
                  "L08_max_excl_15R": round(float(sic[~is15].L08.max()), 2),
                  "lnL08_range_excl_15R": round(float(np.log(sic[~is15].L08.max() / sic[~is15].L08.min())), 3),
                  "L08_15R": round(float(sic[is15].L08.iloc[0]), 2)}
    out["mc"]["n_tierA_eligible"] = rk["main"]["kill3"]["n_tierA_both_T"]
    F("stage4_text_numbers.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("kill3", "framing")}, indent=1, default=str))


if __name__ == "__main__":
    main()
