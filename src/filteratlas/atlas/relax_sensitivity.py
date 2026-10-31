"""Relaxation-limit sensitivity, added after inspection of the results: SevenNet relaxations that hit the
500-step FIRE limit, rerun from the same starting cell with 2,000 steps (intermediate/mlip_relax_sensitivity/).
The main atlas keeps the 500-step rule fixed before the computation; no code of steps 04-08 reads these outputs.
Reported:
1. how many reruns converge within 2,000 steps (the rest stay mlip_failed here too);
2. for those, the difference ln L0.8(2,000 steps) - ln L0.8(500 steps), both phonon sets through the identical code path;
3. the fronts with them added as full-model compounds (assemble.build with the 2,000-step records, the same rows,
   flags and sigma as any SevenNet compound): which compounds enter a Tier-A or Tier-B front or its 1-sigma band;
4. Kendall tau of the top of the ranking against the main atlas: the 20 highest A_fast,1/E of the augmented main-panel
   set, ranked in the augmented and in the main atlas, where a compound absent from the main atlas takes the rank after
   the last (tau-b handles the ties);
5. the temperature-reordering preview (Tier A, 293 vs 77 K; the preview function of fronts) with and without the reruns.
-> the relaxation-sensitivity JSON in data/derived
"""
import json
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from filteratlas.atlas import assemble
from filteratlas.atlas.assemble import FLAG_COLS, MLIP_DIR, eval_mlip
from filteratlas.atlas.fronts import TAGS, fronts, kill3_preview
from filteratlas.common.config import cfg
from filteratlas.common.constants import N_WORKERS
from filteratlas.common.paths import D, F

RUN_DIR = D("mlip_relax_sensitivity") / cfg("mlip")["atlas_model"]
TOP = 20


def last_records(name):
    out = {}
    for line in open(D("logs") / name, encoding="utf-8"):
        if line.strip():
            r = json.loads(line)
            out[r["material_id"]] = r
    return out


def name(r):
    return f"{r.formula} {r.sg_symbol}"


def stats(x):
    x = np.abs(np.asarray(x, float))
    if not len(x):
        return {"n": 0}
    return {"n": int(len(x)), "median_abs": round(float(np.median(x)), 4), "p68_abs": round(float(np.percentile(x, 68)), 4),
            "p95_abs": round(float(np.percentile(x, 95)), 4), "max_abs": round(float(x.max()), 4)}


def compare(at0, at1, added):
    """Front and 1-sigma membership changes, per front key, and the top-of-ranking Kendall tau."""
    s0, m0 = fronts(at0)
    s1, m1 = fronts(at1)
    nm = at1.drop_duplicates("material_id").set_index("material_id")
    lab = lambda ids: sorted(f"{nm.loc[i, 'formula']} {nm.loc[i, 'sg_symbol']}" for i in ids)  # noqa: E731
    out = {}
    for key in sorted(s1):
        f0 = {x["material_id"] for x in s0.get(key, {"front": []})["front"]}
        f1 = {x["material_id"] for x in s1[key]["front"]}
        w0 = set(m0[(m0.key == key) & m0.within_1sigma].material_id)
        w1 = set(m1[(m1.key == key) & m1.within_1sigma].material_id)
        out[key.replace("|", "_")] = {
            "front_main": lab(f0), "front_with_d20": lab(f1), "front_changed": f0 != f1,
            "d20_compounds_on_front": lab(f1 & added), "left_front": lab(f0 - f1),
            "d20_compounds_within_1sigma": lab((w1 - f1) & added), "left_within_1sigma": lab(w0 - w1),
            "n_eligible_main": s0.get(key, {}).get("n_eligible", 0), "n_eligible_with_d20": s1[key]["n_eligible"]}
    tau = {}
    for tag, T in TAGS.items():
        for tier in ("A", "B"):
            sel = lambda a: a[(a["T"] == T) & (a.tier == tier) & a.main_front_eligible & ~a.unstable_dft]  # noqa: E731
            d0, d1 = sel(at0), sel(at1)
            top = d1.nlargest(TOP, "A_fast_1E")
            r1 = d1.A_fast_1E.rank(ascending=False, method="min")
            r0 = d0.set_index("material_id").A_fast_1E.rank(ascending=False, method="min")
            a = r1.loc[top.index].values
            b = np.array([r0.get(m, len(d0) + 1) for m in top.material_id])
            t = kendalltau(a, b).statistic
            tau[f"{tag}_{tier}"] = {"kendall_tau_b": round(float(t), 4), "n_top": len(top),
                                    "d20_compounds_in_top": lab(set(top.material_id) & added)}
    return out, tau


def main():
    at0 = pd.read_parquet(F("atlas.parquet"))
    for c in FLAG_COLS:
        at0[c] = at0[c].fillna(False).astype(bool)
    one0 = at0[at0["T"] > 100]
    failed = set(one0[one0.funnel_step == "mlip_failed"].material_id)
    main_log, sens_log = last_records("mlip_atlas.jsonl"), last_records("mlip_relax_sensitivity.jsonl")
    vm = {json.loads(x)["material_id"] for x in open(D("logs") / "mlip_atlas.jsonl", encoding="utf-8")
          if x.strip() and "machine" in json.loads(x)}
    files = sorted(str(p) for p in RUN_DIR.glob("*.json.gz"))
    with Pool(N_WORKERS) as pool:
        ev = pd.DataFrame(pool.map(partial(eval_mlip, cache_dir=D("relax_sensitivity_eval")), files, chunksize=4))
    conv = ev[(ev.status == "ok") & (ev.relax_converged == True) & ev.material_id.isin(failed)].copy()  # noqa: E712
    # 500-step phonons of the same compounds (the atlas cache; their outputs exist, only the relaxation flag fails them)
    p500 = {m: eval_mlip(str(MLIP_DIR / f"{m}.json.gz")) for m in conv.material_id if (MLIP_DIR / f"{m}.json.gz").exists()}
    rows = []
    for r in conv.to_dict("records"):
        p = p500.get(r["material_id"])
        if p and p.get("status") == "ok":
            rows.append({"material_id": r["material_id"], **{f"dlnL08_{t}": float(np.log(r[f"L08_{t}"] / p[f"L08_{t}"]))
                                                            for t in ("293K", "77K")}})
    dl = pd.DataFrame(rows)
    tried = failed & set(sens_log)
    status = {m: ("converged" if m in set(conv.material_id) else
                  "not_converged" if sens_log[m].get("status") == "ok" else "error") for m in tried}
    origin = lambda m: "D19" if m in vm else "Stage 3"  # noqa: E731
    at1, _ = assemble.build(ml_override=conv)
    added = set(conv.material_id)
    fr, tau = compare(at0, at1, added)
    for a in (at0, at1):
        for c in FLAG_COLS:
            a[c] = a[c].fillna(False).astype(bool)
    k3 = {"main": kill3_preview(at0), "with_d20": kill3_preview(at1)}
    one1 = at1[at1["T"] > 100]
    names = one1.set_index("material_id")
    steps = [sens_log[m]["relax_steps"] for m in added]
    out = {
        "note": "D20 post-hoc sensitivity; the main atlas keeps the 500-step rule (PREREG_STAGE3 §1.3)",
        "max_steps": cfg("mlip")["relax_sensitivity_max_steps"],
        "n_mlip_failed_main": len(failed), "n_rerun": len(tried), "n_not_rerun": len(failed - tried),
        "not_rerun": sorted(f"{names.loc[m, 'formula']} {names.loc[m, 'sg_symbol']}" for m in failed - tried),
        "by_outcome": pd.Series(status).value_counts().to_dict(),
        "by_outcome_and_origin": pd.Series({m: f"{origin(m)}|{s}" for m, s in status.items()}).value_counts().to_dict(),
        "relax_steps_converged": {"median": float(np.median(steps)) if steps else None,
                                  "max": int(max(steps)) if steps else None},
        "errors": {m: sens_log[m].get("error", "")[:120] for m, s in status.items() if s == "error"},
        # n_above_sigma_sevennet_total: |shift| above the whole SevenNet sigma, hypot(sigma_DFT, sigma_MLIP,add);
        # n_above_sigma_mlip_add: above the MLIP term alone (the SevenNet-DFT difference)
        "dlnL08_2000_minus_500": {t: {**(stats(dl[f"dlnL08_{t}"]) if len(dl) else {"n": 0}),
                                      "median_signed": round(float(dl[f"dlnL08_{t}"].median()), 4) if len(dl) else None,
                                      "sigma_sevennet_total": round(float(np.hypot(assemble.SIG_DFT, assemble.SIG_MLIP_ADD[t])), 4),
                                      "sigma_mlip_add": assemble.SIG_MLIP_ADD[t],
                                      "n_above_sigma_sevennet_total": int((dl[f"dlnL08_{t}"].abs() > np.hypot(
                                          assemble.SIG_DFT, assemble.SIG_MLIP_ADD[t])).sum()) if len(dl) else 0,
                                      "n_above_sigma_mlip_add": int((dl[f"dlnL08_{t}"].abs() > assemble.SIG_MLIP_ADD[t]).sum())
                                      if len(dl) else 0}
                                  for t in ("293K", "77K")},
        "largest_shifts_293K": [{"compound": f"{names.loc[r.material_id, 'formula']} {names.loc[r.material_id, 'sg_symbol']}",
                                 "dlnL08_293K": round(r.dlnL08_293K, 4)}
                                for r in dl.reindex(dl.dlnL08_293K.abs().sort_values(ascending=False).index).head(5).itertuples()]
        if len(dl) else [],
        "added_main_front_eligible": {tier: int(one1[one1.material_id.isin(added) & (one1.tier == tier)].main_front_eligible.sum())
                                      for tier in ("A", "B")},
        "added_tierA": sorted(f"{names.loc[m, 'formula']} {names.loc[m, 'sg_symbol']}" for m in added
                              if names.loc[m, "tier"] == "A"),
        "funnel_with_d20": one1.funnel_step.value_counts().to_dict(),
        "fronts": fr, "top_ranking_tau": tau, "top_n": TOP,
        "kill3_preview": {k: {x: v[x] for x in ("n_tierA_both_T", "kendall_tau_A_fast_1E", "tau_ci95",
                                                "n_pareto_status_changes", "changed")} for k, v in k3.items()},
    }
    # what the converged and the non-converged reruns are
    lab = lambda m: f"{names.loc[m, 'formula']} {names.loc[m, 'sg_symbol']}"  # noqa: E731
    vol = pd.Series({m: sens_log[m]["volume_ratio"] for m in added})
    elig = set(one1[one1.material_id.isin(added) & one1.main_front_eligible].material_id)
    dm = dl[dl.material_id.isin(elig)]
    softmol = {m for m in added if bool(names.loc[m, "soft_lattice"]) or bool(names.loc[m, "molecular"])}
    out["converged_reruns"] = {
        "n": len(added), "n_main_front_eligible": len(elig), "n_soft_or_molecular": len(softmol),
        "volume_ratio": {"median": round(float(vol.median()), 4), "n_above_1.1": int((vol > 1.1).sum()),
                         "n_above_1.2": int((vol > 1.2).sum()), "max": round(float(vol.max()), 3),
                         "max_compound": lab(vol.idxmax())},
        "largest_shifts_293K_soft_or_molecular": [x["compound"] for x in out["largest_shifts_293K"]
                                                  if any(lab(m) == x["compound"] for m in softmol)],
        "dlnL08_main_front_eligible": {t: {**(stats(dm[f"dlnL08_{t}"]) if len(dm) else {"n": 0}),
                                           "n_above_sigma_mlip_add": int((dm[f"dlnL08_{t}"].abs() >
                                                                          assemble.SIG_MLIP_ADD[t]).sum())}
                                       for t in ("293K", "77K")},
        "n_in_any_top": len({x for v in tau.values() for x in v["d20_compounds_in_top"]})}
    out["per_compound"] = [{"compound": lab(r.material_id), "dlnL08_293K": round(float(r.dlnL08_293K), 5),
                            "dlnL08_77K": round(float(r.dlnL08_77K), 5), "volume_ratio": round(float(vol[r.material_id]), 4),
                            "main_front_eligible": r.material_id in elig, "soft_or_molecular": r.material_id in softmol}
                           for r in dl.itertuples()]   # per-compound shifts for plotting
    out["compute"] = {"n_reruns_logged": len(sens_log),   # step 08 must not read the rerun log, so this count lives here
                      "n_raised_cap": sum("gpu_alloc" in r for r in sens_log.values())}
    nc = sorted(m for m, s in status.items() if s == "not_converged")
    tab = pd.read_csv(F("stage3_not_computed_table.csv")).set_index("material_id")
    out["not_converged_reruns"] = {
        "n": len(nc), "by_tier": pd.Series([names.loc[m, "tier"] for m in nc]).value_counts().to_dict() if nc else {},
        "relax_fmax_median": round(float(np.median([sens_log[m]["relax_fmax"] for m in nc])), 4) if nc else None,
        "n_volume_ratio_above_1.1": int(sum(sens_log[m]["volume_ratio"] > 1.1 for m in nc)),
        "n_not_excluded_by_proxy_bound": {t: int(tab.reindex(nc)[f"not_excluded_by_proxy_bound_{t}"].fillna(False)
                                                 .astype(bool).sum()) for t in ("293K", "77K")},
        "compounds": [lab(m) for m in nc]}
    out["mlip_failed_not_excluded_covered_by_d20"] = {
        t: int(tab[(tab.funnel_step == "mlip_failed") & tab[f"not_excluded_by_proxy_bound_{t}"].fillna(False)
                   .astype(bool)].index.isin(list(added)).sum()) for t in ("293K", "77K")}
    out["summary"] = {   # roll-up for the reports
        "any_front_changed": any(v["front_changed"] for v in fr.values()),
        "entered_a_front": sorted({x for v in fr.values() for x in v["d20_compounds_on_front"]}),
        "entered_a_within_1sigma_set": sorted({x for v in fr.values() for x in v["d20_compounds_within_1sigma"]}),
        "left_a_within_1sigma_set": sorted({x for v in fr.values() for x in v["left_within_1sigma"]}),
        "min_top_tau": min(v["kendall_tau_b"] for v in tau.values()),
        "kill3_tau_main": k3["main"]["kendall_tau_A_fast_1E"], "kill3_tau_with_d20": k3["with_d20"]["kendall_tau_A_fast_1E"],
        "kill3_changes_main": k3["main"]["n_pareto_status_changes"],
        "kill3_changes_with_d20": k3["with_d20"]["n_pareto_status_changes"],
        "n_converged": int(out["by_outcome"].get("converged", 0)),
        "n_not_converged": int(out["by_outcome"].get("not_converged", 0)),
        "n_error": int(out["by_outcome"].get("error", 0))}
    F("stage3_relax_sensitivity.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_rerun", "by_outcome", "dlnL08_2000_minus_500", "top_ranking_tau")}, indent=1))


if __name__ == "__main__":
    main()
