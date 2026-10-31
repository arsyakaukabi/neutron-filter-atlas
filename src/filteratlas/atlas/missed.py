"""What could the funnel have missed? Added after inspection of the results.

1. Front-axis test at the Debye level: every survivor without a full model (not_computed_budget, not_run_memory,
   mlip_failed, supercell_too_large) that has a Debye L0.8 gets A_fast,1/E and A_gamma,1MeV at its Debye L0.8, with the
   MP summary-structure density. A compound is a "Debye-level front contender" if this point is not dominated by the
   main front of its tier at that temperature (the main fronts are full-model).
2. Upper-bound test: mu_proxy <= mu_th up to the tolerance measured in step 02 (L0.8,full / L0.8,proxy <= 1.0167 at
   293 K and 1.0195 at 77 K, triviality-test statistics -> proxy_lower_bound), so the proxy L0.8 times that factor bounds
   L0.8 from above and the rejection there bounds the front position from above. Compounds not dominated at that bound
   cannot be ruled out.
3. Density bias of not_run_memory: median atom number density by supercell-size bin, per SevenNet outcome.
4. Supplementary table of every compound without a full model.
-> the missed-compounds JSON and the not-computed table (CSV) in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.fom import rejection as rj  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402
from filteratlas.atlas.fronts import pareto  # noqa: E402

TAGS = {"293K": 293.15, "77K": 77.0}
NOFULL = ("not_computed_budget", "not_run_memory", "mlip_failed", "supercell_too_large")


def rejection(formula, n, L):
    if not (L == L) or L <= 0:
        return np.nan, np.nan
    return float(rj.a_fast(rj.sigma_tot(formula), n, L, "1E")[0]), rj.mu_gamma(formula, n)[0] * L


def dominated(x, y, fx, fy):
    return bool(np.any((fx >= x) & (fy >= y) & ((fx > x) | (fy > y))))


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    fun = pd.read_parquet(F("stage3_funnel.parquet"))
    prox = pd.read_parquet(F("stage3_proxy.parquet"))[["material_id", "n_summary_A3"]]
    one = at[at["T"] > 100][["material_id", "formula", "sg_symbol", "tier", "funnel_step", "proxy_L08",
                             "flag_Be", "flag_Pb", "flag_Bi", "flag_Tl", "flag_As"]]
    tol = {t: json.loads((F("stage1_stats.json")).read_text())["proxy_lower_bound"][f"L08_{t}"]["max_ratio_full_over_proxy"]
           for t in TAGS}
    nf = one[one.funnel_step.isin(NOFULL)].merge(fun[["material_id", "debye_L08_293K", "debye_L08_77K"]],
                                                  on="material_id", how="left").merge(prox, on="material_id", how="left")
    out = {"proxy_bound_tolerance": tol, "n_without_full_model": len(nf), "by_label": nf.funnel_step.value_counts().to_dict(),
           "n_with_debye": int(nf.debye_L08_293K.notna().sum()),
           "n_without_debye": int(nf.debye_L08_293K.isna().sum())}
    rows = []
    for tag, T in TAGS.items():
        fr = {}
        for tier in ("A", "B"):
            d = at[(at["T"] == T) & (at.tier == tier) & at.main_front_eligible & ~at.unstable_dft]
            m = pareto(d.A_fast_1E, d.A_gamma_1MeV)
            fr[tier] = (d.A_fast_1E.values[m], d.A_gamma_1MeV.values[m])
        a_d, g_d, a_p, g_p, cont, bound = [], [], [], [], [], []
        for r in nf.itertuples():
            fx, fy = fr[r.tier]
            a1, g1 = rejection(r.formula, r.n_summary_A3, getattr(r, f"debye_L08_{tag}"))
            p1, q1 = rejection(r.formula, r.n_summary_A3, r.proxy_L08 * tol[tag])
            a_d.append(a1); g_d.append(g1); a_p.append(p1); g_p.append(q1)
            cont.append(a1 == a1 and not dominated(a1, g1, fx, fy))
            bound.append(p1 == p1 and not dominated(p1, q1, fx, fy))
        nf[f"A_fast_1E_debye_{tag}"], nf[f"A_gamma_1MeV_debye_{tag}"] = a_d, g_d
        nf[f"A_fast_1E_proxybound_{tag}"], nf[f"A_gamma_1MeV_proxybound_{tag}"] = a_p, g_p
        nf[f"debye_front_contender_{tag}"], nf[f"not_excluded_by_proxy_bound_{tag}"] = cont, bound
        c = nf[nf[f"debye_front_contender_{tag}"]]
        b = nf[nf[f"not_excluded_by_proxy_bound_{tag}"]]
        out[tag] = {"debye_front_contenders": [f"{r.formula} {r.sg_symbol} (Tier {r.tier}, {r.funnel_step})"
                                               for r in c.itertuples()],
                    "n_debye_front_contenders": {"A": int((c.tier == "A").sum()), "B": int((c.tier == "B").sum())},
                    "n_not_excluded_by_proxy_bound": {"A": int((b.tier == "A").sum()), "B": int((b.tier == "B").sum())},
                    "not_excluded_by_proxy_bound_by_label": b.funnel_step.value_counts().to_dict()}
    # density bias: SevenNet set only, drift arm excluded
    sset = pd.DataFrame(json.loads((F("stage3_mlip_set.json")).read_text())["compounds"])
    s = sset[sset.arm != "drift"][["material_id", "supercell_atoms"]].merge(one[["material_id", "funnel_step"]], on="material_id") \
        .merge(prox, on="material_id")
    dens = {}
    for lo, hi in ((0, 100), (100, 200), (200, 250), (250, 300), (300, 330), (330, 400), (400, 480), (480, 700)):   # bins above 330 atoms resolve the largest supercells
        g = s[(s.supercell_atoms > lo) & (s.supercell_atoms <= hi)]
        dens[f"{lo}-{hi}"] = {k: {"median_n_A3": round(float(v.n_summary_A3.median()), 4), "n": len(v)}
                              for k, v in g.groupby("funnel_step")}
    out["density_by_supercell_bin"] = dens
    out["density_all"] = {k: round(float(v.n_summary_A3.median()), 4) for k, v in s.groupby("funnel_step")}
    (F("stage3_missed.json")).write_text(json.dumps(out, indent=1, default=str))
    cols = ["material_id", "formula", "sg_symbol", "tier", "funnel_step", "n_summary_A3", "proxy_L08",
            "debye_L08_293K", "debye_L08_77K"] + [c for c in nf.columns if c.startswith(("A_fast", "A_gamma", "debye_front",
                                                                                           "not_excluded"))]
    nf[cols].sort_values(["tier", "funnel_step", "proxy_L08"], ascending=[True, True, False]).round(4) \
        .to_csv(F("stage3_not_computed_table.csv"), index=False)
    print(json.dumps({k: v for k, v in out.items() if k not in ("density_by_supercell_bin",)}, indent=1, default=str)[:3500])
    print(json.dumps(out["density_by_supercell_bin"], indent=0)[:1500])


if __name__ == "__main__":
    main()
