"""Slab-transport results of the shortlist and the framing, from the OpenMC runs.

Per case: the uncollided result against the analytic one (check: within 2 %), the transported rejection
A_transport = -ln(transmitted fraction above 0.625 eV for neutron sources; above 10 keV for photons), the in-scatter
correction A_transport - A_uncollided, the thermalised fraction from fast sources, and the thermal-source transmission
(uncollided vs analytic T_th; all transmitted thermal neutrons).
Effect: fronts within the shortlist with transported values in place of uncollided ones; the pairwise P(X beats
sapphire) recomputed with each compound's in-scatter correction at nominal L0.8 as a fixed shift in every draw
(same draws as the main Monte Carlo); the framing rule (fixed before the computation) applied to the survivors.
The 1/E source is taken from the reruns with all NCrystal components (coherent elastic kept; raw_1E_all/), because
"Bragg off" applied up to 5 eV removed the coherent-elastic part of the free-atom cross section; the runs with the
original setting are kept and reported. Causes of check failures are verified case by case, not assigned by rule. The
statements and the framing are also evaluated for the common-component scenario (fixed before the computation), and
the fronts for the ratio model (added after inspection of the results).
-> derived/stage4_transport.json
"""
import json

import numpy as np
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F
from filteratlas.shortlist import rank_mc as mc

C = cfg("shortlist")
KEYS = {"1E": "A_fast_1E", "watt": "A_fast_watt", "g1MeV": "A_gamma_1MeV", "g6MeV": "A_gamma_6MeV"}


def load_cases(sub="raw"):
    rows = []
    for f in sorted((D("stage4_transport") / sub).glob("*.json")):
        r = json.loads(f.read_text())
        ff = r["fom_fractions"]
        rec = {k: r[k] for k in ("material_id", "formula", "sg_symbol", "tag", "thickness", "L_cm", "source")}
        rec["rel_sd_target_met"] = r["rel_sd_target_met"]
        rec["particles"] = r["raw"]["particles"]
        if r["source"] == "thermal":
            rec["T_th_uncollided"], rec["T_th_analytic"] = ff["uncollided"][0], r["analytic"]["T_th"]
            rec["T_th_all_thermal"] = float(r["raw"]["all"][0])
            rec["check_rel_diff"] = rec["T_th_uncollided"] / rec["T_th_analytic"] - 1
        else:
            k = KEYS[r["source"]]
            a_unc, a_tr = -np.log(ff["uncollided"][0]), -np.log(ff["transported"][0])
            rec.update({"axis": k, "A_analytic": r["analytic"][k], "A_uncollided": a_unc, "A_transport": a_tr,
                        "correction": a_tr - a_unc, "check_rel_diff": a_unc / r["analytic"][k] - 1})
            if not r["source"].startswith("g"):
                rec["thermalised_fraction"] = float(r["raw"]["all"][0])
        rec["check_pass"] = abs(rec["check_rel_diff"]) < C["transport"]["check_tolerance"]
        rows.append(rec)
    return pd.DataFrame(rows)


CASE = ["material_id", "tag", "thickness", "source"]


def main():
    pre = load_cases()                       # the original setting (Bragg off below 5 eV for every source)
    alt = load_cases("raw_1E_all")           # 1/E source with all NCrystal components
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    at = mc.load()
    if len(alt) != (pre.source == "1E").sum():
        raise RuntimeError(f"raw_1E_all has {len(alt)} runs, expected {(pre.source == '1E').sum()}")
    df = pd.concat([pre[pre.source != "1E"], alt], ignore_index=True)
    tol = C["transport"]["check_tolerance"]
    out = {"n_runs": len(df), "n_runs_1E_all_components": len(alt),
           "check_preregistered_setting": {
               "n_failed": int((~pre.check_pass).sum()), "n_runs": len(pre),
               "failed": pre[~pre.check_pass][CASE + ["check_rel_diff"]].to_dict("records"),
               "max_abs_check_rel_diff": {s: float(pre[pre.source == s].check_rel_diff.abs().max()) for s in sorted(pre.source.unique())}},
           "all_checks_pass": bool(df.check_pass.all()),
           "checks_failed": df[~df.check_pass][CASE + ["check_rel_diff"]].to_dict("records"),
           "rel_sd_target_met_all": bool(df.rel_sd_target_met.all()),
           "max_abs_check_rel_diff": {s: float(df[df.source == s].check_rel_diff.abs().max()) for s in sorted(df.source.unique())}}
    # causes verified case by case. A 1/E failure in the original setting is explained if the same case passes
    # with all NCrystal components (coherent elastic restored); a Watt failure on a Be compound is attributed to (n,2n)
    # secondaries (which also start with zero collisions) only if the free-atom rerun shows the same difference
    diag_f = D("stage4_transport") / "diagnostic_free_atom.json"
    diag = json.loads(diag_f.read_text()) if diag_f.exists() else []
    out["free_atom_diagnostic"] = diag
    dg = {d["case"]: d for d in diag}
    alt_pass = {tuple(r[k] for k in CASE): (r["check_pass"], r["check_rel_diff"]) for r in alt.to_dict("records")}
    fail = pre[~pre.check_pass]
    def why(r):
        if r.source == "1E":
            ok, rd = alt_pass[(r.material_id, r.tag, r.thickness, "1E")]
            if ok:
                return "coherent elastic removed up to 5 eV; passes with all NCrystal components"
            d = dg.get(f"{r.material_id}_{r.tag}_{r.thickness}_1E")
            if d and abs(d["rel_diff"]) < tol:
                return ("NCrystal's model below 5 eV differs from the free-atom data even with all components "
                        "(free-atom rerun passes)")
            return "unexplained"
        if r.source == "watt" and "Be" in r.formula:
            d = dg.get(f"{r.material_id}_{r.tag}_{r.thickness}_watt")
            if d and abs(d["rel_diff"] - r.check_rel_diff) < tol / 2:
                return "consistent with (n,2n) secondaries counted as zero-collision (free-atom rerun shows the same deviation)"
        return "unexplained"
    cause = fail.apply(why, axis=1) if len(fail) else pd.Series(dtype=str)
    out["checks_failed_by_cause"] = fail.assign(cause=cause).groupby("cause").size().to_dict()
    out["checks_unexplained"] = fail[cause == "unexplained"][CASE + ["check_rel_diff"]].to_dict("records")
    fm = df[~df.check_pass]   # what still fails in the setting used for the results
    out["checks_failed_by_cause_used_setting"] = (fm.assign(cause=fm.apply(why, axis=1)).groupby("cause").size().to_dict()
                                                  if len(fm) else {})
    out["check_1E_all_components_max_abs_rel_diff"] = float(alt.check_rel_diff.abs().max())
    out["check_explained_case_by_case"] = not out["checks_unexplained"]
    # per-compound table at L0.8
    L = df[df.thickness == "L08"]
    lab = {r.material_id: f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]" for r in sl.itertuples()}
    table = {}
    for (mid, tag), g in L.groupby(["material_id", "tag"]):
        e = {}
        for r in g.itertuples():
            if r.source == "thermal":
                e["T_th_uncollided"], e["T_th_all_thermal"] = r.T_th_uncollided, r.T_th_all_thermal
            else:
                e[r.axis] = {"uncollided": round(r.A_uncollided, 4), "transport": round(r.A_transport, 4),
                             "correction": round(r.correction, 4), "analytic": round(r.A_analytic, 4)}
                if r.source in ("1E", "watt"):
                    e[r.axis]["thermalised_fraction"] = round(r.thermalised_fraction, 5)
        table.setdefault(lab[mid], {})[tag] = e
    out["at_L08"] = table
    five = df[df.thickness == "5cm"]
    out["at_5cm"] = {f"{lab[m]} {t}": {r.axis if r.source != "thermal" else "thermal":
                                       (round(r.correction, 4) if r.source != "thermal" else round(r.T_th_uncollided, 4))
                                       for r in g.itertuples()} for (m, t), g in five.groupby(["material_id", "tag"])}
    out["correction_range"] = {k: [float(L[L.axis == k].correction.min()), float(L[L.axis == k].correction.max())]
                               for k in KEYS.values()}
    # fronts within the shortlist (main axes), uncollided vs transported, per tier and T
    corr = {(r.material_id, r.tag, r.axis): r.correction for r in L[L.source != "thermal"].itertuples()}
    ratio = {(r.material_id, r.tag, r.axis): r.A_transport / r.A_uncollided for r in L[L.source != "thermal"].itertuples()}
    fronts = {}
    for tag, T in mc.TAGS.items():
        for tier in ("A", "B"):
            ids = [m for m in sl[(sl.tier == tier) & sl.main_front_eligible].material_id]
            d = at[(at["T"] == T) & at.material_id.isin(ids)].set_index("material_id").loc[ids]
            x0, y0 = d.A_fast_1E.values, d.A_gamma_1MeV.values
            x1 = x0 + np.array([corr[(m, tag, "A_fast_1E")] for m in ids])
            y1 = y0 + np.array([corr[(m, tag, "A_gamma_1MeV")] for m in ids])
            x2 = x0 * np.array([ratio[(m, tag, "A_fast_1E")] for m in ids])
            y2 = y0 * np.array([ratio[(m, tag, "A_gamma_1MeV")] for m in ids])
            f0 = {lab[m] for m, k in zip(ids, mc.front_mask(x0[None], y0[None])[0]) if k}
            f1 = {lab[m] for m, k in zip(ids, mc.front_mask(x1[None], y1[None])[0]) if k}
            f2 = {lab[m] for m, k in zip(ids, mc.front_mask(x2[None], y2[None])[0]) if k}
            fronts[f"{tag}_{tier}"] = {"uncollided": sorted(f0), "transported": sorted(f1), "changed": f0 != f1,
                                       "transported_ratio_model": sorted(f2)}
    out["shortlist_fronts"] = fronts
    # from the original setting (Bragg off below 5 eV, the atlas convention); the all-components 1/E
    # reruns add powder Bragg scattering in the thermal band, which changes this fraction by -6 to +11 %
    tf = pre[(pre.thickness == "L08") & pre.source.isin(["1E", "watt"])]
    out["thermalised_fraction_range"] = {s: [float(tf[tf.source == s].thermalised_fraction.min()),
                                             float(tf[tf.source == s].thermalised_fraction.max())] for s in ("1E", "watt")}
    # pairwise statements with the corrections as fixed shifts (main scenario draws)
    ids = mc.members(at, sl)
    u, gr = mc.grids(at, ids)
    pos = {m: i for i, m in enumerate(ids)}
    sapph = at[(at.formula == "Al2O3") & (at.sg_symbol == "R-3c") & (at.tier == "A")].material_id.iloc[0]
    N = C["rank_mc"]["n_draws"]
    sl_idx = np.array([pos[m] for m in sl.material_id])

    def beats_for(scenario):
        """P(X beats sapphire) with the corrections as fixed shifts and (ratio model) as ratios, for the draws of the
        given rank_mc scenario (the same random sequence as rank_mc.run, so the untransported probabilities match)."""
        rng = np.random.default_rng(C["rank_mc"]["seed"])
        beats = {t: {k: np.zeros(len(sl)) for k in KEYS.values()} for t in mc.TAGS}
        beats_r = {t: {k: np.zeros(len(sl)) for k in KEYS.values()} for t in mc.TAGS}
        done = 0
        while done < N:
            nd = min(mc.CHUNK, N - done)
            z = rng.standard_normal((nd, len(ids)))
            if scenario == "common":
                c = C["rank_mc"]["common_component"]
                z = np.sqrt(c) * rng.standard_normal((nd, 1)) + np.sqrt(1 - c) * z
            if scenario == "sapphire_offset":   # as in rank_mc
                z[:, pos[sapph]] += mc.sapphire_offset() / gr["293K"]["sB"][pos[sapph]]
            for t in mc.TAGS:
                A = mc.values(gr[t], u, z, gr[t]["sB"][None, :], gr[t]["sB"][None, :])
                for k in KEYS.values():
                    shift = np.array([corr.get((m, t, k), 0.0) for m in sl.material_id])
                    cs = corr.get((sapph, t, k), 0.0)
                    beats[t][k] += ((A[k][:, sl_idx] + shift[None, :]) > (A[k][:, [pos[sapph]]] + cs)).sum(0)
                    rat = np.array([ratio.get((m, t, k), 1.0) for m in sl.material_id])
                    rs = ratio.get((sapph, t, k), 1.0)
                    beats_r[t][k] += ((A[k][:, sl_idx] * rat[None, :]) > (A[k][:, [pos[sapph]]] * rs)).sum(0)
            done += nd
        return beats, beats_r

    p = C["rank_mc"]["p_statement"]
    rk = json.loads(F("stage4_rank_probabilities.json").read_text())
    key4 = lambda s_: (s_["T"], s_["axis"], s_["compound"], s_["statement"])  # noqa: E731
    base = {key4(s_) for s_ in rk["main"]["statements"]}
    tierA_nonknown = {lab[r.material_id] for r in sl.itertuples() if r.tier == "A" and not r.known_filter}
    flags = {lab[r.material_id]: r.flags for r in sl.itertuples()}
    k3 = rk["main"]["kill3"]

    def statements(bt):
        stm = []
        for t in mc.TAGS:
            for k in KEYS.values():
                for m, q in zip(sl.material_id, bt[t][k] / N):
                    if m != sapph and (q > p or q < 1 - p):
                        stm.append({"T": t, "axis": k, "compound": lab[m], "P": round(float(q), 4),
                                    "statement": "beats sapphire" if q > p else "sapphire beats it"})
        return stm

    def framing(stm, bt, scenario="main"):
        # discovery needs a non-known Tier-A compound with P(beats sapphire) > 0.90 on A_fast,1/E or
        # A_gamma,1MeV at either T, not beaten by sapphire (P >= 0.10) on the other axis at that T, after transport
        # "the statement must survive transport": it has to qualify before transport (main Monte Carlo) and still
        # qualify after it; statements that only appear after transport are reported, not counted
        def qualifying(st):
            pt = {(s_["T"], s_["axis"], s_["compound"]): s_ for s_ in st}
            q = {}
            for s_ in st:
                if s_["compound"] in tierA_nonknown and s_["statement"] == "beats sapphire" and s_["axis"] in ("A_fast_1E", "A_gamma_1MeV"):
                    other = "A_gamma_1MeV" if s_["axis"] == "A_fast_1E" else "A_fast_1E"
                    o = pt.get((s_["T"], other, s_["compound"]))
                    if not (o and o["statement"] == "sapphire beats it"):
                        q[(s_["T"], s_["axis"], s_["compound"])] = s_
            return q
        before, after = qualifying(rk[scenario]["statements"]), qualifying(stm)
        def rec(key, s_):
            i = next(j for j, m in enumerate(sl.material_id) if lab[m] == s_["compound"])
            other = "A_gamma_1MeV" if s_["axis"] == "A_fast_1E" else "A_fast_1E"
            return {**s_, "flags": flags[s_["compound"]], "P_other_axis_after_transport": round(float(bt[s_["T"]][other][i] / N), 4)}
        disc = [rec(k, after[k]) for k in after if k in before]
        return {"kill3_fires": k3["fires"], "K_holds": k3["K_holds"], "R_holds": k3["R_holds"],
                "qualifying_before_transport": [before[k] for k in before],
                "discovery_statements": disc,
                "transport_induced_not_counted": [rec(k, after[k]) for k in after if k not in before],
                "lost_by_transport": [before[k] for k in before if k not in after],
                "framing": "atlas + method" if k3["fires"] else ("discovery" if disc else "atlas + method with candidates")}

    beats, beats_r = beats_for("main")
    # the full P(X beats sapphire) after transport, so every untransported statement gets its status
    out["p_beats_sapphire_transported"] = {t: {k: {lab[m]: round(float(q), 4) for m, q in zip(sl.material_id, beats[t][k] / N)
                                                   if m != sapph} for k in KEYS.values()} for t in mc.TAGS}
    stm = statements(beats)
    new = {(s_["T"], s_["axis"], s_["compound"], s_["statement"]) for s_ in stm}
    out["statements_transported"] = stm
    out["statements_changed_by_transport"] = {"lost": sorted(map(list, base - new)), "gained": sorted(map(list, new - base))}
    out["framing"] = framing(stm, beats)
    out["framing"]["candidates"] = {c: {t: rk["main"][t]["p_front"]["main_A"].get(c, 0.0) for t in mc.TAGS}
                                    for c in sorted(tierA_nonknown)
                                    if max(rk["main"][t]["p_front"]["main_A"].get(c, 0.0) for t in mc.TAGS) >= 0.05}
    stm_r = statements(beats_r)
    new_r = {(s_["T"], s_["axis"], s_["compound"], s_["statement"]) for s_ in stm_r}
    out["sensitivity_ratio_model"] = {   # added after inspection of the results: A_transport / A_uncollided at nominal L0.8 as a fixed ratio
        "statements": stm_r, "changed_vs_untransported": {"lost": sorted(map(list, base - new_r)),
                                                          "gained": sorted(map(list, new_r - base))},
        "framing": framing(stm_r, beats_r)}
    # the common-component scenario (fixed before the computation), before and after transport (shift and ratio models)
    cb, cbr = beats_for("common")
    base_c = {key4(s_) for s_ in rk["common"]["statements"]}
    sc_ = {}
    for nm, bt in (("shift", cb), ("ratio", cbr)):
        st_ = statements(bt)
        new_c = {key4(s_) for s_ in st_}
        sc_[nm] = {"statements": st_, "changed_vs_untransported": {"lost": sorted(map(list, base_c - new_c)),
                                                                   "gained": sorted(map(list, new_c - base_c))},
                   "framing": framing(st_, bt, "common")}
    out["sensitivity_common_component"] = sc_
    # sapphire conditioned on its measured residual (added after inspection of the results), before and after transport
    if "sapphire_offset" in rk:
        ob, _ = beats_for("sapphire_offset")
        st_ = statements(ob)
        out["sensitivity_sapphire_offset"] = {"offset_lnL08": mc.sapphire_offset(), "statements": st_,
                                              "framing": framing(st_, ob, "sapphire_offset")}
    F("stage4_transport.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_runs", "all_checks_pass", "max_abs_check_rel_diff", "correction_range",
                                          "framing")}, indent=1, default=str))
    print("ratio model framing:", json.dumps(out["sensitivity_ratio_model"]["framing"], default=str))


if __name__ == "__main__":
    main()
