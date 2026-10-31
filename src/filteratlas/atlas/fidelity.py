"""Statistics of the triviality test (does a composition-only proxy rank the crystals as well as the full model).

Fixed before the computation: gate statistic (Spearman rho of L0.8 full vs proxy, usable = full mu_th(1.8 A, 293 K)
< 0.3 /cm), secondary statistics, usable-cut / E_hull sensitivity, Debye-from-elastic ablation, lower-bound check.
Added after the first comparison: selection-symmetric variants, band-edge and moderator temperature variants from
stored curves, universal-theta_D nulls and Debye top-k overlap, formula-clustered bootstrap, partial correlations for
'who moves', molecular definition excluding metal-containing compounds, magnitude of lower-bound violations.
Output: the statistics JSON, the usable-set table and the controls JSON in data/derived.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.core import Composition, Element
from scipy.stats import kendalltau, rankdata, spearmanr

from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import KILL_RHO, N_BOOT, USABLE_MU_CUT, USABLE_SENSITIVITY  # noqa: E402

from filteratlas.common.config import cfg

_S1 = cfg("analysis")["stage1"]
RNG = np.random.default_rng(_S1["seed_stats"])
DEBYE_OUTLIER = _S1["debye_outlier_factor"]
TAGS = ("293K", "77K")
NOBLE = {"He", "Ne", "Ar", "Kr", "Xe"}


def rho_ci(x, y, n=N_BOOT, clusters=None):
    """Spearman rho with a percentile bootstrap CI (sampling-only). clusters: resample whole groups."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    r = spearmanr(x, y).statistic
    if clusters is None:
        idx = RNG.integers(0, len(x), size=(n, len(x)))
    else:
        g = pd.Series(np.arange(len(x))).groupby(np.asarray(clusters)).apply(list).tolist()
        idx = [np.concatenate([g[j] for j in RNG.integers(0, len(g), len(g))]) for _ in range(n)]
    boots = np.array([spearmanr(x[i], y[i]).statistic for i in idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return {"rho": round(float(r), 4), "ci95": [round(float(lo), 4), round(float(hi), 4)], "n": int(len(x))}


def overlap(a, b, k, ascending=False):
    ta = set(a.sort_values(ascending=ascending).index[:k])
    tb = set(b.sort_values(ascending=ascending).index[:k])
    return len(ta & tb)


def gate_block(u):
    out = {}
    for t in TAGS:
        out[t] = {
            "primary_L08_vs_proxy": rho_ci(u[f"L08_{t}"], u[f"L08_proxy_{t}"]),
            "primary_L08_vs_proxy_cluster_bootstrap": rho_ci(u[f"L08_{t}"], u[f"L08_proxy_{t}"], clusters=u.formula),
            "kendall_L08": round(float(kendalltau(u[f"L08_{t}"], u[f"L08_proxy_{t}"]).statistic), 4),
            "top10_overlap_L08": overlap(u[f"L08_{t}"], u[f"L08_proxy_{t}"], 10),
            "top30_overlap_L08": overlap(u[f"L08_{t}"], u[f"L08_proxy_{t}"], 30),
            "mu1.0_vs_proxy": rho_ci(u[f"mu_1.0_{t}"], u[f"proxy_1.0_{t}"]),
            "mu1.8_vs_proxy": rho_ci(u[f"mu_1.8_{t}"], u[f"proxy_1.8_{t}"]),
            "kendall_mu1.0": round(float(kendalltau(u[f"mu_1.0_{t}"], u[f"proxy_1.0_{t}"]).statistic), 4),
            "top10_overlap_mu1.0": overlap(u[f"mu_1.0_{t}"], u[f"proxy_1.0_{t}"], 10, True),
            "hot_band_L08_vs_proxy": rho_ci(u[f"L08_hot_{t}"], u[f"L08_hot_proxy_{t}"]),
            "median_frac_mu1.8_not_in_proxy": round(float(np.median(1 - u[f"proxy_1.8_{t}"] / u[f"mu_1.8_{t}"])), 4),
            "median_frac_mu1.0_not_in_proxy": round(float(np.median(1 - u[f"proxy_1.0_{t}"] / u[f"mu_1.0_{t}"])), 4),
        }
    out["kill_primary"] = bool(all(out[t]["primary_L08_vs_proxy"]["rho"] > KILL_RHO for t in TAGS))
    out["secondary_above_kill_both_T"] = [k for k in ("mu1.0_vs_proxy", "mu1.8_vs_proxy", "hot_band_L08_vs_proxy")
                                          if all(out[t][k]["rho"] > KILL_RHO for t in TAGS)]
    out["secondary_above_kill_any_T"] = sorted({f"{k}@{t}" for t in TAGS for k in
                                                ("mu1.0_vs_proxy", "mu1.8_vs_proxy", "hot_band_L08_vs_proxy")
                                                if out[t][k]["rho"] > KILL_RHO})
    return out


def mean_mass(formula):
    c = Composition(formula)
    return float(sum(el.atomic_mass * n for el, n in c.items()) / c.num_atoms)


def sigma_over_mass(formula):
    """Atom-averaged bound scattering cross section per unit mass [b/u] (the proxy's blind spot)."""
    import NCrystal as NC
    c = Composition(formula)
    return float(sum(NC.atomDB(el.symbol).scatteringXS() * n for el, n in c.items())
                 / sum(el.atomic_mass * n for el, n in c.items()))


def partial_spearman(y, x, controls):
    """Spearman partial correlation of y and x given controls (rank-transform, OLS residuals, Pearson)."""
    R = lambda v: rankdata(v)  # noqa: E731
    Z = np.column_stack([np.ones(len(y))] + [R(c) for c in controls])
    ry = R(y) - Z @ np.linalg.lstsq(Z, R(y), rcond=None)[0]
    rx = R(x) - Z @ np.linalg.lstsq(Z, R(x), rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


def is_molecular(formula, dim):
    """0-D / intercalated-molecule AND no metal atom, or noble-gas-only."""
    els = {e.symbol for e in Composition(formula).elements}
    if els <= NOBLE:
        return True
    metal = any(Element(e).is_metal for e in els)
    return dim in ("0D", "intercalated molecule") and not metal


def main(tag=""):
    df = pd.read_parquet(F("stage1_triviality.parquet"))
    cand = pd.read_parquet(F("stage0_candidates.parquet"))[["material_id", "e_hull", "density",
                                                                     "crystal_system"]]
    df = df.merge(cand, on="material_id", how="left")
    counts = {"computed": len(df), "status_ok": int((df.status == "ok").sum()),
              "unstable": int(((df.status == "ok") & ~df.stable.fillna(False)).sum())}
    base = df[(df.status == "ok") & df.stable.fillna(False)].copy().set_index("material_id")
    counts["sample"] = len(base)
    u = base[base["mu_1.8_293K"] < USABLE_MU_CUT].copy()
    counts["usable"] = len(u)
    res = {"counts": counts, "gate": gate_block(u)}

    # --- sensitivity fixed before the computation --------------------------------------------------------------------------
    sens = {}
    for cut in USABLE_SENSITIVITY:
        uu = base[base["mu_1.8_293K"] < cut]
        sens[f"usable_cut_{cut}"] = {t: {**rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500),
                                         "top10_overlap": overlap(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 10)}
                                     for t in TAGS}
    for eh in (0.0, 0.025):
        uu = u[u.e_hull <= eh + 1e-9]
        sens[f"e_hull_le_{eh}"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}
    uu = u[u["flags"].fillna("") == ""]
    sens["no_flagged_elements"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}
    uu = u[u.phonon_method == "pheasy"]
    sens["pheasy_only"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}

    # --- selection-symmetric variants ---------------------------------------------------------------------------
    sel = {}
    uu = base[base["proxy_1.8_293K"] < USABLE_MU_CUT]
    sel["select_on_proxy_mu1.8_293K"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}
    uu = base[(base["mu_1.8_293K"] < USABLE_MU_CUT) | (base["proxy_1.8_293K"] < USABLE_MU_CUT)]
    sel["select_union_full_or_proxy"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}
    uu = base[base["mu_1.8_77K"] < USABLE_MU_CUT]
    sel["select_77K_on_own_mu1.8"] = {"77K": rho_ci(uu["L08_77K"], uu["L08_proxy_77K"], 500)}
    for cut in (_S1["filter_grade_mu_cut_per_cm"],):
        uu = base[base[f"proxy_1.8_293K"] < cut]
        sel[f"select_on_proxy_cut_{cut}"] = {t: rho_ci(uu[f"L08_{t}"], uu[f"L08_proxy_{t}"], 500) for t in TAGS}
    res["selection_symmetric"] = sel

    # --- molecular (metal-containing compounds excluded from the definition) / dedup ------------------------------
    dim = pd.read_parquet(F("stage1_dimensionality.parquet")).set_index("material_id")
    u = u.join(dim, how="left")
    u["molecular"] = [is_molecular(f, d) for f, d in zip(u.formula, u.dim.fillna(""))]
    nonmol = u[~u.molecular]
    sens["exclude_molecular"] = {t: rho_ci(nonmol[f"L08_{t}"], nonmol[f"L08_proxy_{t}"], 500) for t in TAGS}
    dedup = u.sort_values("e_hull").groupby("formula").head(1)
    sens["one_polymorph_per_formula"] = {t: {**rho_ci(dedup[f"L08_{t}"], dedup[f"L08_proxy_{t}"], 500),
                                             "top10_overlap": overlap(dedup[f"L08_{t}"], dedup[f"L08_proxy_{t}"], 10)}
                                         for t in TAGS}
    both = dedup[~dedup.molecular]
    sens["exclude_molecular_and_dedup"] = {t: {**rho_ci(both[f"L08_{t}"], both[f"L08_proxy_{t}"], 500),
                                               "top10_overlap": overlap(both[f"L08_{t}"], both[f"L08_proxy_{t}"], 10)}
                                           for t in TAGS}
    res["n_molecular_in_usable"] = int(u.molecular.sum())
    res["examples_molecular"] = sorted(u[u.molecular].formula.unique().tolist())
    res["dimensionality_failures_in_usable"] = int(u.dim.fillna("error").str.startswith("error").sum())
    res["sensitivity"] = sens

    # --- band edges and moderator temperature from stored curves -----------------------------------------------
    bands = {k: ((v[0], v[1]), v[2]) for k, v in _S1["band_sensitivity"].items()}   # label: (band, moderator T)
    bres = {}
    for name, (band, Tm) in bands.items():
        lt = xe.lambda_T(Tm)
        bres[name] = {}
        for t in TAGS:
            lf = np.array([xe.l_target_band(np.asarray(c), band, lt) for c in u[f"curve_full_{t}"]])
            lp = np.array([xe.l_target_band(np.asarray(c), band, lt) for c in u[f"curve_proxy_{t}"]])
            bres[name][t] = {**rho_ci(lf, lp, 300), "top10_overlap": len(set(np.argsort(-lf)[:10]) &
                                                                          set(np.argsort(-lp)[:10]))}
    res["band_moderator_sensitivity"] = bres

    # --- whole-sample view ---------------------------------------------------------------------------------------
    res["all_stable_compounds"] = {t: rho_ci(base[f"L08_{t}"], base[f"L08_proxy_{t}"], 500) for t in TAGS}

    # --- temperature ---------------------------------------------------------------------------------------------
    res["temperature"] = {"L08_293_vs_77": rho_ci(u["L08_293K"], u["L08_77K"]),
                          "top10_overlap_293_vs_77": overlap(u["L08_293K"], u["L08_77K"], 10),
                          "median_gain_L08_77_over_293": round(float(np.median(u["L08_77K"] / u["L08_293K"])), 3)}

    # --- lower bound (with magnitudes) ----------------------------------------------------------------------------
    lb = {}
    for t in TAGS:
        for w in ("1.0", "1.8"):
            r = base[f"proxy_{w}_{t}"] / base[f"mu_{w}_{t}"]
            lb[f"{w}_{t}"] = {"frac": round(float((r > 1).mean()), 4), "max_ratio": round(float(r.max()), 4),
                              "worst": base.loc[r.idxmax(), "formula"]}
        r = base[f"L08_{t}"] / base[f"L08_proxy_{t}"]
        lb[f"L08_{t}"] = {"frac": round(float((r > 1).mean()), 4), "max_ratio_full_over_proxy": round(float(r.max()), 4),
                          "worst": base.loc[r.idxmax(), "formula"]}
    res["proxy_lower_bound"] = lb

    # --- nulls: universal theta_D and elastic Debye --------------------------------------------------------------
    nulls = {}
    for th in (300, 600):
        nulls[f"universal_theta_{th}K"] = {t: {**rho_ci(u[f"L08_{t}"], u[f"univ{th}_L08_{t}"], 500),
                                              "top10_overlap_with_full": overlap(u[f"L08_{t}"], u[f"univ{th}_L08_{t}"], 10),
                                              "top30_overlap_with_full": overlap(u[f"L08_{t}"], u[f"univ{th}_L08_{t}"], 30)}
                                          for t in TAGS}
    d = u.dropna(subset=["debye_L08_293K"])
    deb = {"n_usable_with_elastic": len(d)}
    for t in TAGS:
        ratio = d[f"debye_L08_{t}"] / d[f"L08_{t}"]
        deb[t] = {**rho_ci(d[f"L08_{t}"], d[f"debye_L08_{t}"], 500),
                  "frac_outside_factor_1.5": round(float(((ratio > DEBYE_OUTLIER) | (ratio < 1 / DEBYE_OUTLIER)).mean()), 4),
                  "median_ratio_debye_over_full": round(float(ratio.median()), 3),
                  "top10_overlap_with_full": overlap(d[f"L08_{t}"], d[f"debye_L08_{t}"], 10),
                  "top30_overlap_with_full": overlap(d[f"L08_{t}"], d[f"debye_L08_{t}"], 30),
                  "proxy_top10_overlap_with_full_same_subset": overlap(d[f"L08_{t}"], d[f"L08_proxy_{t}"], 10),
                  "proxy_rho_same_subset": rho_ci(d[f"L08_{t}"], d[f"L08_proxy_{t}"], 300)["rho"]}
    pre = _S1["debye_prescreen"]
    deb["usable_as_exclusion_prescreen"] = bool(all(deb[t]["rho"] >= pre["min_rho"]
                                                    and deb[t]["frac_outside_factor_1.5"] < pre["max_frac_outside"]
                                                    for t in TAGS))
    nulls["debye_from_elastic"] = deb
    res["nulls"] = nulls

    # --- who moves (partial correlations) -------------------------------------------------------------------------
    u["mean_mass"] = [mean_mass(f) for f in u.formula]
    u["sigma_over_mass"] = [sigma_over_mass(f) for f in u.formula]
    moves = {}
    for t in TAGS:
        u[f"rank_full_{t}"] = u[f"L08_{t}"].rank(ascending=False)
        u[f"rank_proxy_{t}"] = u[f"L08_proxy_{t}"].rank(ascending=False)
        u[f"shift_{t}"] = u[f"rank_full_{t}"] - u[f"rank_proxy_{t}"]  # >0: full ranks it worse than proxy
        moves[t] = {f"raw_{f}": round(float(spearmanr(u[f"shift_{t}"], u[f]).statistic), 3)
                    for f in ("mean_phonon_meV", "mean_mass", "sigma_over_mass")}
        moves[t]["partial_mean_phonon_given_proxyrank_sigmaovermass"] = round(partial_spearman(
            u[f"shift_{t}"], u["mean_phonon_meV"], [u[f"rank_proxy_{t}"], u["sigma_over_mass"]]), 3)
        moves[t]["partial_sigmaovermass_given_proxyrank"] = round(partial_spearman(
            u[f"shift_{t}"], u["sigma_over_mass"], [u[f"rank_proxy_{t}"]]), 3)
    res["who_moves"] = moves
    cols = ["formula", "sg_symbol", "flags", "L08_293K", "L08_proxy_293K", "rank_full_293K", "rank_proxy_293K",
            "shift_293K", "mean_phonon_meV", "sigma_over_mass"]
    res["top_movers_down"] = u.sort_values("shift_293K", ascending=False).head(10)[cols].round(3).reset_index().to_dict(
        "records")
    for t in TAGS:
        res[f"top10_full_dedup_nonmolecular_{t}"] = both.sort_values(f"L08_{t}", ascending=False).head(10)[
            ["formula", "sg_symbol", f"L08_{t}"]].round(2).reset_index().to_dict("records")
        res[f"top10_proxy_dedup_nonmolecular_{t}"] = both.sort_values(f"L08_proxy_{t}", ascending=False).head(10)[
            ["formula", "sg_symbol", f"L08_proxy_{t}"]].round(2).reset_index().to_dict("records")

    keep = [c for c in u.columns if not c.startswith("curve_")]
    u[keep].reset_index().to_parquet(F("stage1_usable.parquet"))
    (F("stage1_stats.json")).write_text(json.dumps(res, indent=1, default=str))
    ctrl = pd.read_parquet(F("stage1_controls.parquet"))
    (F("stage1_controls.json")).write_text(json.dumps(
        {r["control"].replace(" ", "_"): {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}
         for r in ctrl.to_dict("records")}, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ("counts", "gate")}, indent=1))


if __name__ == "__main__":
    main()
