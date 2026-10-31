"""Consequences of the full dense-mesh stability screen.

Reads stage1r_stability_full.parquet (from stability_screen.run(all_usable=True)) and reports:
  (a) prevalence of dense-mesh instability (< -0.5 THz) by cell size, with Wilson 95 % CIs, and how many of the
      unstable entries the MP-DOS weight rule (fixed before the computation) missed;
  (b) stability of the top-30 lists (full and proxy; one polymorph per formula, non-molecular), anchors, controls;
  (c) instability rate among the elastic-Debye outliers (outside x1.5) vs the rest;
  (d) pheasy-vs-DFPT reproducibility restricted to pairs whose pheasy entry is dense-mesh clean, and a
      SIMEX that draws noise from the empirical pair-difference distribution (heavy tails kept) instead of a
      Gaussian, on the full usable set and on the dynamically clean set;
  (e) large-cell floor validation: FC-tetrahedron VDOS vs MP DOS (2 meV floor) for 20 random dense-mesh-clean
      usable compounds with 21-60 sites (fc_validation covers <= 20 sites only).
-> stage1r_stability_summary.json, stage1r_fc_validation_largecell.parquet
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.config import cfg
from filteratlas.common.constants import TEMPS, UNSTABLE_WEIGHT_MAX  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.dynamics.fc_validation import fc_vdos  # noqa: E402
from filteratlas.dynamics.stability_screen import wilson  # noqa: E402

_S = cfg("stability")
STRATA = [tuple(x) for x in _S["summary_strata_nsites"]]
RNG = np.random.default_rng(_S["summary_seed"])
LARGE_N, LARGE_SEED = _S["large_cell_check"]["n"], _S["large_cell_check"]["seed"]
FC_UNSTABLE = _S["fc_unstable_THz"]
NEAR_ZERO = _S["fc_near_zero_THz"]          # "clean" in the floor checks: no frequency below -0.05 THz
DEBYE_OUTLIER = cfg("analysis")["stage1"]["debye_outlier_factor"]   # elastic-Debye outlier criterion (x1.5)


def simex_empirical(x, y, resid, lambdas=(0.5, 1.0, 1.5, 2.0), B=100):
    """SIMEX with noise drawn (with replacement) from an empirical per-entry error sample `resid` (log units),
    scaled by sqrt(lambda); quadratic extrapolation to lambda = -1."""
    lx = np.log(np.asarray(x, float))
    r0 = spearmanr(lx, y).statistic
    pts = [(0.0, r0)]
    for lam in lambdas:
        pts.append((lam, float(np.mean([spearmanr(lx + np.sqrt(lam) * RNG.choice(resid, len(lx)), y).statistic
                                        for _ in range(B)]))))
    L, R = np.array(pts).T
    return {"rho_obs": round(float(r0), 4), "rho_simex": round(float(np.polyval(np.polyfit(L, R, 2), -1.0)), 4)}


def main():
    u = pd.read_parquet(F("stage1_usable.parquet"))
    sf = pd.read_parquet(F("stage1r_stability_full.parquet"))
    m = u.merge(sf[["material_id", "status", "fc_min_freq_THz", "unstable_05", "unstable_comm_05"]]
                .rename(columns={"status": "fc_status"}),
                on="material_id")
    ok = m[m.fc_status == "ok"].copy()
    ok["unstable_05"] = ok.unstable_05.astype(bool)
    ok["unstable_comm_05"] = ok.unstable_comm_05.astype(bool)
    tag = lambda r: f"{r.formula} {r.sg_symbol}" + ("" if r.unstable_comm_05 == True else " (interpolation-only)")  # noqa: E712,E731
    res = {}

    # (a) prevalence by cell size
    prev = {}
    for lo, hi in STRATA:
        g = ok[(ok.nsites >= lo) & (ok.nsites <= hi)]
        k = int(g.unstable_05.sum())
        prev[f"{lo}-{hi if hi < 10_000 else 'max'}"] = {"n": len(g), "unstable": k,
                                                         "unstable_commensurate": int(g.unstable_comm_05.sum()),
                                                         "frac": round(k / len(g), 4) if len(g) else None,
                                                         "wilson95": wilson(k, len(g))}
    k = int(ok.unstable_05.sum())
    prev["all"] = {"n": len(ok), "unstable": k, "frac": round(k / len(ok), 4), "wilson95": wilson(k, len(ok))}
    prev["missed_by_mp_dos_rule"] = int((ok.unstable_05 & (ok.unstable_weight <= UNSTABLE_WEIGHT_MAX)).sum())
    prev["unstable_with_zero_mp_dos_weight"] = int((ok.unstable_05 & (ok.unstable_weight == 0)).sum())
    res["prevalence_by_nsites"] = prev
    # stratified-sample screen vs full screen: are the per-compound minima identical (same method, same inputs)?
    st = pd.read_parquet(F("stage1r_stability_screen.parquet"))
    st = st[st.status == "ok"].merge(sf[sf.status == "ok"][["material_id", "fc_min_freq_THz"]], on="material_id",
                                     suffixes=("", "_full"))
    strat = json.loads((F("stage1r_stability_screen.json")).read_text())["prevalence"]["all_random_unweighted"]
    full = prev["all"]["frac"]
    res["full_vs_stratified_sample"] = {
        "full_frac": full, "stratified_frac": strat["frac"], "stratified_wilson95": strat["wilson95"],
        "full_within_stratified_ci": bool(strat["wilson95"][0] <= full <= strat["wilson95"][1]),
        "relation": "within" if strat["wilson95"][0] <= full <= strat["wilson95"][1] else "outside"}
    res["stratified_vs_full"] = {"n_shared": len(st), "n_identical_1e-6": int((np.abs(st.fc_min_freq_THz
                                                                                  - st.fc_min_freq_THz_full) < 1e-6).sum())}

    # (b) top lists, anchors, controls
    dd = m[~m.molecular].sort_values("e_hull").groupby("formula").head(1)
    tops = {}
    for t in ("293K", "77K"):
        for lab, col in (("full", f"L08_{t}"), ("proxy", f"L08_proxy_{t}")):
            top = dd.nlargest(30, col)
            tops[f"top30_{lab}_{t}"] = {"n_screened": int((top.fc_status == "ok").sum()),
                                        "n_unstable": int((top.unstable_05 == True).sum()),  # noqa: E712
                                        "n_unscreenable": int((top.fc_status != "ok").sum()),
                                        "n_unstable_commensurate": int((top.unstable_comm_05 == True).sum()),  # noqa: E712
                                        "unstable": [tag(r) for r in top[top.unstable_05 == True].itertuples()]}  # noqa: E712
    res["top_lists"] = tops
    anchors = [("Al2O3", "R-3c"), ("Si", "Fd-3m"), ("MgO", "Fm-3m"), ("MgF2", "P4_2/mnm"), ("BeF2", "P3_121"),
               ("BeO", "P6_3mc"), ("PbF2", "Fm-3m"), ("SiO2", "P3_121"), ("SiO2", "P3_221"), ("AlPO4", "P3_121"),
               ("AlPO4", "P3_221"), ("TeO2", "P4_12_12"), ("TeO2", "P4_32_12")]
    an = []
    for f, sg in anchors:
        g = m[(m.formula == f) & (m.sg_symbol == sg)]
        for r in g.itertuples():
            an.append({"compound": f"{f} {sg}", "status": r.fc_status,
                       "fc_min_freq_THz": None if pd.isna(r.fc_min_freq_THz) else round(float(r.fc_min_freq_THz), 2),
                       "unstable": bool(r.unstable_05 == True),  # noqa: E712
                       "commensurate_unstable": bool(r.unstable_comm_05 == True)})  # noqa: E712
    res["anchors_controls"] = an
    res["anchors_controls_unstable"] = [f"{a['compound']} ({a['fc_min_freq_THz']} THz" +
                                        ("" if a["commensurate_unstable"] else ", interpolation-only") + ")"
                                        for a in an if a["unstable"]]
    res["anchors_controls_unscreenable"] = [a["compound"] for a in an if a["status"] != "ok"]

    # (c) elastic-Debye outliers
    d = ok.dropna(subset=["debye_L08_293K"])
    out = {}
    for t in ("293K", "77K"):
        r = np.abs(np.log(d[f"debye_L08_{t}"] / d[f"L08_{t}"])) > np.log(DEBYE_OUTLIER)
        out[t] = {"n_outliers": int(r.sum()), "unstable_among_outliers": int((r & d.unstable_05).sum()),
                  "frac_unstable_outliers": round(float(d[r].unstable_05.mean()), 3) if r.any() else None,
                  "frac_unstable_rest": round(float(d[~r].unstable_05.mean()), 3),
                  "frac_outside_1.5_clean_only": round(float(r[~d.unstable_05].mean()), 4),
                  "n_clean_with_elastic": int((~d.unstable_05).sum())}
    res["elastic_debye_outliers_vs_stability"] = out

    # (d) reproducibility on the dense-clean population
    pr = pd.read_parquet(F("stage1r_pheasy_vs_dfpt.parquet"))
    pr = pr[(pr.status == "ok") & pr.stable_pheasy.fillna(False) & pr.stable_dfpt.fillna(False)
            & pr.material_id.isin(u.material_id)].merge(m[["material_id", "nsites", "fc_status", "unstable_05"]],
                                                           on="material_id")
    rep = {"n_pairs_usable": len(pr), "nsites_max": int(pr.nsites.max()), "nsites_median": float(pr.nsites.median()),
           "usable_frac_nsites_gt20": round(float((u.nsites > 20).mean()), 3)}
    for lab, g in (("all_pairs", pr), ("pheasy_dense_clean", pr[(pr.fc_status == "ok") & (pr.unstable_05 == False)])):  # noqa: E712
        rep[lab] = {"n": len(g)}
        for t in ("293K", "77K"):
            dl = np.log(g[f"pheasy_L08_{t}"] / g[f"dfpt_L08_{t}"])
            rep[lab][t] = {"median_abs_dlog": round(float(dl.abs().median()), 4),
                           "p90_abs_dlog": round(float(dl.abs().quantile(0.9)), 4),
                           "p99_abs_dlog": round(float(dl.abs().quantile(0.99)), 4),
                           "max_abs_dlog": round(float(dl.abs().max()), 4),
                           "spearman": round(float(spearmanr(g[f"pheasy_L08_{t}"], g[f"dfpt_L08_{t}"]).statistic), 4)}
    fv = pd.read_parquet(F("stage1r_fc_validation.parquet"))
    fv = fv[(fv.status == "ok") & (fv.fc_min_freq_THz < NEAR_ZERO)]
    for t in ("293K", "77K"):
        a = np.abs(np.log(fv[f"floor2_L08_{t}"] / fv[f"fc_L08_{t}"]))
        rep[f"fc_unclean_abs_dlog_range_{t}"] = [round(float(a.min()), 2), round(float(a.max()), 2)]
    rep["fc_unclean_n"] = len(fv)
    # Does a dense-mesh instability corrupt the MP-DOS L0.8 itself? Two direct tests:
    # (i) pheasy-vs-DFPT pairs whose pheasy entry is unstable (same compound, independent workflow);
    # (ii) enantiomorph pairs with one unstable member (physically identical crystals).
    uns = pr[(pr.fc_status == "ok") & (pr.unstable_05 == True)]  # noqa: E712
    rep["pheasy_dense_unstable"] = {"n": len(uns)}
    for t in ("293K", "77K"):
        dl = np.log(uns[f"pheasy_L08_{t}"] / uns[f"dfpt_L08_{t}"]).abs()
        rep["pheasy_dense_unstable"][t] = {"median_abs_dlog": round(float(dl.median()), 4),
                                           "p90_abs_dlog": round(float(dl.quantile(0.9)), 4),
                                           "max_abs_dlog": round(float(dl.max()), 4)}
    # how many DFPT partners of these pairs share the instability (any MP-DOS weight below
    # -0.5 THz)? Where both workflows are unstable, their agreement shows reproducibility, not correctness.
    cand = pd.read_parquet(F("stage0_candidates.parquet")).set_index("material_id")
    shared = []
    for mid in uns.material_id:
        dos = mq.load_dos(D("mp") / "phonon_dos" / "dfpt" / f"{cand.loc[mid, 'phonon_dfpt'].split(',')[0]}.json.gz")
        shared.append(nb.stability(dos)[0] > 0)
    rep["pheasy_dense_unstable"]["n_dfpt_partner_also_negative_weight"] = int(np.sum(shared))
    # ... and how strongly their L0.8 depends on the near-zero regularisation (floor 1 vs 3 meV)
    fsr = pd.read_parquet(F("stage1r_floor_sensitivity.parquet"))
    fsr = fsr[fsr.ok.astype(bool)].pivot(index="material_id", columns="floor", values="L08_293K")
    dfl = np.log(fsr[3.0] / fsr[1.0]).abs().rename("dfl").reset_index().merge(ok[["material_id", "unstable_05"]],
                                                                              on="material_id")
    rep["floor_1_vs_3meV_abs_dlog_293K"] = {
        lab_: {"n": len(g), "median": round(float(g.dfl.median()), 4), "frac_gt_0.05": round(float((g.dfl > 0.05).mean()), 3),
               "max": round(float(g.dfl.max()), 3)}
        for lab_, g in (("unstable", dfl[dfl.unstable_05]), ("clean", dfl[~dfl.unstable_05]))}
    en = json.loads((F("stage1_enantiomorphs.json")).read_text())
    fmin = sf.set_index("material_id").fc_min_freq_THz
    ep = []
    for p in en["pairs"]:
        fa, fb = fmin.get(p["id_a"]), fmin.get(p["id_b"])
        if fa is None or fb is None or pd.isna(fa) or pd.isna(fb) or ((fa < FC_UNSTABLE) == (fb < FC_UNSTABLE)):
            continue
        ep.append({"pair": f"{p['formula']} {p['sg_a']}/{p['sg_b']}",
                   "unstable_member": p["sg_a"] if fa < FC_UNSTABLE else p["sg_b"],
                   "min_freq_THz": round(float(min(fa, fb)), 2),
                   "abs_dlog_293K": round(float(p["abs_log_ratio_293K"]), 3),
                   "abs_dlog_77K": round(float(p["abs_log_ratio_77K"]), 3)})
    rep["enantiomorph_one_unstable"] = ep
    rep["enantiomorph_one_unstable_text"] = [f"{e['pair']}: {e['abs_dlog_293K']}" for e in ep]
    # FC reference for unstable compounds: soft/imaginary modes put spurious weight near zero frequency (no floor
    # is applied to the FC DOS), so the FC L0.8 is itself unreliable there; report the MP/FC ratio direction.
    rep["fc_unclean_frac_mp_above_fc"] = round(float((fv["floor2_L08_293K"] > fv["fc_L08_293K"]).mean()), 3)
    m60 = F("stage1r_stability_full_mesh60.json")
    if m60.exists():
        j = json.loads(m60.read_text())
        res["mesh60_sensitivity"] = {"unstable_05": j["unstable_05"],
                                     "screened_clean_only_0.5": j["gate_variants"]["screened_clean_only_0.5"],
                                     "top10_overlap_clean_dedup_293K": j["top10_overlap_clean_dedup_293K"],
                                     "top10_overlap_clean_dedup_77K": j["top10_overlap_clean_dedup_77K"]}
    clean_set = ok[~ok.unstable_05]
    sim = {}
    for t in ("293K", "77K"):
        dl = np.log(pr[f"pheasy_L08_{t}"] / pr[f"dfpt_L08_{t}"]).values / np.sqrt(2)   # per-entry error sample
        dl = dl - np.median(dl)
        sim[t] = {"primary_all_usable": simex_empirical(m[f"L08_{t}"], m[f"L08_proxy_{t}"], dl),
                  "dense_clean_screened": simex_empirical(clean_set[f"L08_{t}"], clean_set[f"L08_proxy_{t}"], dl)}
    rep["simex_empirical_heavy_tail"] = sim
    res["reproducibility_restricted"] = rep

    # (e) large-cell floor validation
    pool = ok[(~ok.unstable_05) & (ok.nsites > 20) & (ok.fc_min_freq_THz >= NEAR_ZERO) & (ok.phonon_method == "pheasy")]
    sel = pool.sort_values("material_id").sample(min(LARGE_N, len(pool)), random_state=LARGE_SEED)
    rows = []
    for row in sel.itertuples():
        rec = {"material_id": row.material_id, "formula": row.formula, "nsites": row.nsites}
        try:
            dos = mq.load_dos(D("mp") / "phonon_dos" / row.phonon_method / f"{row.phonon_id}.json.gz")
            vd, fmin = fc_vdos(None, row)
            ev_fc = xe.evaluate(nb.build_vdos_material(dos, f"lc_fc_{row.material_id}", vdos=vd), TEMPS)
            ev_2 = xe.evaluate(nb.build_vdos_material(dos, f"lc_f2_{row.material_id}", floor_meV=2.0), TEMPS)
            rec.update({"fc_min_freq_THz": fmin, "status": "ok",
                        **{f"fc_L08_{t}": ev_fc[f"L08_{t}"] for t in ("293K", "77K")},
                        **{f"floor2_L08_{t}": ev_2[f"L08_{t}"] for t in ("293K", "77K")}})
        except Exception as e:  # noqa: BLE001
            rec.update({"status": "error", "error": repr(e)[:200]})
        print(rec["formula"], rec["status"], flush=True)
        rows.append(rec)
    lc = pd.DataFrame(rows)
    lc.to_parquet(F("stage1r_fc_validation_largecell.parquet"))
    lok = lc[lc.status == "ok"]
    val = {"n_selected": len(lc), "n_ok": len(lok), "nsites_range": [int(lok.nsites.min()), int(lok.nsites.max())]}
    for t in ("293K", "77K"):
        dl = np.log(lok[f"floor2_L08_{t}"] / lok[f"fc_L08_{t}"])
        val[t] = {"rms_dlog": round(float(np.sqrt((dl ** 2).mean())), 4), "median_dlog": round(float(dl.median()), 4),
                  "max_abs_dlog": round(float(dl.abs().max()), 4)}
    res["largecell_floor_validation"] = val
    (F("stage1r_stability_summary.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
