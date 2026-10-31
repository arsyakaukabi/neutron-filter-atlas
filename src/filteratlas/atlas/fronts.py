"""Analysis of the atlas table (data/derived/atlas.parquet):
Pareto fronts (main and secondary, per temperature and tier, plus the outside-domain panel), the "within 1 sigma"
sets, the Tier-A fidelity ladder, who moves, the sensitivity runs, and a preview of the test of whether the growable
front contains anything beyond the known filters and whether temperature reorders it (decided in step 10).
-> the fronts JSON and the front-membership table (per compound x T) in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.core import Composition
from scipy.stats import kendalltau, spearmanr

from filteratlas.fom import rejection as rj
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import N_BOOT  # noqa: E402

TAGS = {"293K": 293.15, "77K": 77.0}
MAIN = ("A_fast_1E", "A_gamma_1MeV")
SECONDARY = ("A_fast_watt", "A_gamma_6MeV")
from filteratlas.common.config import cfg  # noqa: E402
from filteratlas.atlas.assemble import SIG_MLIP_ADD, SIG_SOFT  # noqa: E402

_U = cfg("uncertainty")
SEED = cfg("funnel")["seed"]   # a fresh generator per statistic, so intervals do not depend on call order
KNOWN = set(_U["known_filters"])
SIG_D = _U["sigma_construction_D"]
F_BLOCK = set(_U["f_block_4f"])   # partly filled 4f shells
TAU_REORDER = _U["kill3_tau_threshold"]


def pareto(x, y):
    """Boolean mask of non-dominated points (maximise both)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    order = np.lexsort((-y, -x))
    keep = np.zeros(len(x), bool)
    best = -np.inf
    for i in order:
        if y[i] > best:
            keep[i] = True
            best = y[i]
    # exact ties in x with lower y are dominated; of two identical points only the first is kept
    return keep


def within(xhi, yhi, fx, fy):
    """+1 sigma point not dominated by any front member's central point."""
    out = np.ones(len(xhi), bool)
    for a, b in zip(fx, fy):
        out &= ~((a >= xhi) & (b >= yhi) & ((a > xhi) | (b > yhi)))
    return out


def fronts(at):
    rows, summ = [], {}
    for tag, T in TAGS.items():
        d = at[(at["T"] == T) & at.vdos_source.ne("none") & ~at.unstable_dft]
        f4 = d.formula.map(lambda x: bool(set(Composition(x).as_dict()) & F_BLOCK))
        vc = d.main_front_eligible & ~d.untested_heavy & ~f4      # no Z >= 80, no partly filled 4f shell
        for tier in ("A", "B"):
            for panel, sub in (("main", d[(d.tier == tier) & d.main_front_eligible]),
                               ("outside_domain", d[(d.tier == tier) & d.outside_domain]),
                               ("no_heavy_no_4f", d[(d.tier == tier) & vc])):
                for kind, (ax, ay) in (("main", MAIN), ("secondary", SECONDARY)):
                    if sub.empty:
                        continue
                    m = pareto(sub[ax], sub[ay])
                    w = within(sub[f"{ax}_hi"], sub[f"{ay}_hi"], sub[ax][m], sub[ay][m])
                    key = f"{tag}|{tier}|{panel}|{kind}"
                    f = sub[m].sort_values(ax, ascending=False)
                    summ[key] = {"n_eligible": len(sub), "n_front": int(m.sum()), "n_within_1sigma": int(w.sum()),
                                 "front": f[["material_id", "formula", "sg_symbol", ax, ay, "L08", "vdos_source"]]
                                 .round(4).to_dict("records")}
                    rows.append(pd.DataFrame({"material_id": sub.material_id, "T": T, "key": key, "front": m,
                                              "within_1sigma": w}))
    return summ, pd.concat(rows, ignore_index=True)


def rho_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    r = spearmanr(x, y).statistic
    rng = np.random.default_rng(SEED)
    bs = [spearmanr(x[i], y[i]).statistic for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    return {"rho": round(float(r), 4), "ci95": [round(float(np.nanpercentile(bs, 2.5)), 4),
                                                 round(float(np.nanpercentile(bs, 97.5)), 4)], "n": len(x)}


def overlap(d, a, b, k):
    return len(set(d.nlargest(k, a).material_id) & set(d.nlargest(k, b).material_id))


def ladder(at):
    out = {}
    for tag, T in TAGS.items():
        d = at[(at["T"] == T) & (at.tier == "A") & at.vdos_source.ne("none") & ~at.unstable_dft].copy()
        for sub_name, sub in (("all", d), ("domain", d[d.main_front_eligible])):
            res = {"n": len(sub)}
            res["proxy_vs_full"] = rho_ci(sub.proxy_L08, sub.L08)
            res["proxy_vs_full_top10_overlap"] = overlap(sub, "proxy_L08", "L08", 10)
            res["proxy_vs_full_top30_overlap"] = overlap(sub, "proxy_L08", "L08", 30)
            e = sub.dropna(subset=["debye_L08"])
            res["n_with_debye"] = len(e)
            res["n_with_debye_sevennet"] = int((e.vdos_source == "sevennet").sum())   # elastic Debye value from the Debye table
            if len(e) >= 10:
                res["debye_vs_full"] = rho_ci(e.debye_L08, e.L08)
                res["proxy_vs_debye"] = rho_ci(e.proxy_L08, e.debye_L08)
                res["debye_vs_full_top10_overlap"] = overlap(e, "debye_L08", "L08", 10)
                res["debye_vs_full_top30_overlap"] = overlap(e, "debye_L08", "L08", 30)
                res["proxy_vs_full_top10_overlap_on_debye_subset"] = overlap(e, "proxy_L08", "L08", 10)
            if len(e) >= 10:   # like-for-like proxy on the Debye subset; order within the full top 10
                res["proxy_vs_full_on_debye_subset"] = rho_ci(e.proxy_L08, e.L08)
                top = e.nlargest(10, "L08")
                res["kendall_debye_vs_full_top10_full"] = round(float(kendalltau(top.debye_L08, top.L08).statistic), 3)
                e_a1 = [float(rj.a_fast(rj.sigma_tot(f), n, L, "1E")[0]) for f, n, L in zip(e.formula, e.numberdensity_A3, e.debye_L08)]
                e_g1 = [rj.mu_gamma(f, n)[0] * L for f, n, L in zip(e.formula, e.numberdensity_A3, e.debye_L08)]
                m_e = e.main_front_eligible.values
                res["debye_level_front"] = sorted(f"{r.formula} {r.sg_symbol}" for r in
                                                  e[m_e][pareto(np.array(e_a1)[m_e], np.array(e_g1)[m_e])].itertuples())
                res["full_front_on_debye_subset"] = sorted(f"{r.formula} {r.sg_symbol}" for r in
                                                           e[m_e][pareto(e.A_fast_1E[m_e], e.A_gamma_1MeV[m_e])].itertuples())
            # secondary: A_fast,1/E rankings, proxy A at the proxy L0.8 (summary density)
            sub = sub.assign(A_fast_1E_proxy=[float(rj.a_fast(rj.sigma_tot(f), n, L, "1E")[0]) for f, n, L in
                                              zip(sub.formula, sub.n_summary_A3, sub.proxy_L08)])
            res["A_fast_1E_proxy_vs_full"] = rho_ci(sub.A_fast_1E_proxy, sub.A_fast_1E)
            res["A_fast_1E_top10_overlap"] = overlap(sub, "A_fast_1E_proxy", "A_fast_1E", 10)
            res["top10_full"] = sub.nlargest(10, "L08")[["formula", "sg_symbol", "L08", "vdos_source"]].round(3).to_dict("records")
            res["top10_proxy"] = sub.nlargest(10, "proxy_L08")[["formula", "sg_symbol", "proxy_L08"]].round(3).to_dict("records")
            out[f"{tag}|{sub_name}"] = res
    return out


def who_moves(at):
    a = at[(at.tier == "A") & at.vdos_source.ne("none") & ~at.unstable_dft]
    w = a.pivot_table(index=["material_id", "formula", "sg_symbol"], columns="T", values=["L08", "proxy_L08"]).dropna()
    w.columns = [f"{v}_{'293K' if t > 100 else '77K'}" for v, t in w.columns]
    w = w.reset_index()
    for t in TAGS:
        w[f"rank_full_{t}"] = w[f"L08_{t}"].rank(ascending=False)
        w[f"rank_proxy_{t}"] = w[f"proxy_L08_{t}"].rank(ascending=False)
        w[f"move_proxy_to_full_{t}"] = w[f"rank_proxy_{t}"] - w[f"rank_full_{t}"]
    w["move_293_to_77"] = w["rank_full_293K"] - w["rank_full_77K"]
    w["L08_ratio_77_over_293"] = w["L08_77K"] / w["L08_293K"]
    cols = ["formula", "sg_symbol", "rank_proxy_293K", "rank_full_293K", "rank_full_77K", "L08_ratio_77_over_293"]
    return {"n": len(w),
            "up_proxy_to_full_293K": w.nlargest(8, "move_proxy_to_full_293K")[cols + ["move_proxy_to_full_293K"]].round(2).to_dict("records"),
            "down_proxy_to_full_293K": w.nsmallest(8, "move_proxy_to_full_293K")[cols + ["move_proxy_to_full_293K"]].round(2).to_dict("records"),
            "up_293_to_77": w.nlargest(8, "move_293_to_77")[cols + ["move_293_to_77"]].round(2).to_dict("records"),
            "down_293_to_77": w.nsmallest(8, "move_293_to_77")[cols + ["move_293_to_77"]].round(2).to_dict("records"),
            "median_abs_move_proxy_to_full": {t: float(w[f"move_proxy_to_full_{t}"].abs().median()) for t in TAGS},
            "median_abs_move_293_to_77": float(w.move_293_to_77.abs().median())}


def tierA_main(at, T):
    return at[(at["T"] == T) & (at.tier == "A") & at.main_front_eligible & ~at.unstable_dft]


def recompute(d, L):
    """A_fast,1E and A_gamma,1MeV at new thicknesses L (same density and composition)."""
    a1 = [float(rj.a_fast(rj.sigma_tot(f), n, x, "1E")[0]) for f, n, x in zip(d.formula, d.numberdensity_A3, L)]
    g1 = [rj.mu_gamma(f, n)[0] * x for f, n, x in zip(d.formula, d.numberdensity_A3, L)]
    return np.array(a1), np.array(g1)


def front_ids(d, x, y):
    return set(d.material_id[pareto(x, y)])


def sensitivity(at):
    cur = pd.read_parquet(F("stage3_curves.parquet")).set_index("material_id")
    out = {}
    for tag, T in TAGS.items():
        d = tierA_main(at, T)
        base = front_ids(d, d.A_fast_1E, d.A_gamma_1MeV)
        res = {"base_front": sorted(d[d.material_id.isin(base)].formula + " " + d[d.material_id.isin(base)].sg_symbol)}

        def cmp(name, ids):
            res[name] = {"front": sorted(at[at.material_id.isin(ids)].drop_duplicates("material_id").apply(
                lambda r: f"{r.formula} {r.sg_symbol}", axis=1)), "added": len(ids - base), "removed": len(base - ids),
                "changed": ids != base}
        for name, band, lt in (("band_0.4-3.5A", (0.4, 3.5), xe.LAMBDA_T), ("cold_40K_0.5-4.5A", (0.5, 4.5), xe.lambda_T(40.0))):
            dd = d[d.material_id.isin(cur.index)]
            L = np.array([xe.l_target_band(np.asarray(cur.loc[m, f"curve_full_{tag}"]), band, lt) for m in dd.material_id])
            a1, g1 = recompute(dd, L)
            cmp(name, front_ids(dd, a1, g1))
        for name, mask in (("exclude_interp_only", ~d.interp_only), ("exclude_untested_heavy", ~d.untested_heavy.fillna(False)),
                           ("exclude_mlip_unstable", ~d.mlip_unstable.fillna(False)),
                           ("exclude_unscreenable", ~d.unscreenable.fillna(False))):   # added after inspection of the results
            dd = d[mask]
            cmp(name, front_ids(dd, dd.A_fast_1E, dd.A_gamma_1MeV))
        # sigma construction D: the "within 1 sigma" set changes (point fronts do not)
        m = pareto(d.A_fast_1E, d.A_gamma_1MeV)
        wB = set(d.material_id[within(d.A_fast_1E_hi, d.A_gamma_1MeV_hi, d.A_fast_1E[m], d.A_gamma_1MeV[m])])
        sD = np.where(d.vdos_source == "mp_dft", SIG_D, np.hypot(SIG_D, SIG_MLIP_ADD[tag]))
        sD = np.maximum(sD, np.where(d.soft_lattice.fillna(False), SIG_SOFT, 0.0))
        Lhi = d.L08 * np.exp(sD)
        a1, g1 = recompute(d, Lhi)
        wD = set(d.material_id[within(a1, g1, d.A_fast_1E[m], d.A_gamma_1MeV[m])])
        res["sigma_D_within_1sigma"] = {"n_B": len(wB), "n_D": len(wD), "added": sorted(wD - wB), "removed": sorted(wB - wD)}
        out[tag] = res
    return out


def kill3_preview(at):
    a, b = tierA_main(at, 293.15), tierA_main(at, 77.0)
    m = a.merge(b, on="material_id", suffixes=("_293", "_77"))
    tau = kendalltau(m.A_fast_1E_293, m.A_fast_1E_77).statistic
    bs, rng = [], np.random.default_rng(SEED)
    for _ in range(N_BOOT):
        i = rng.integers(0, len(m), len(m))
        bs.append(kendalltau(m.A_fast_1E_293.values[i], m.A_fast_1E_77.values[i]).statistic)
    fa = front_ids(a, a.A_fast_1E, a.A_gamma_1MeV)
    fb = front_ids(b, b.A_fast_1E, b.A_gamma_1MeV)
    changed = sorted((fa ^ fb) & set(m.material_id))
    names = at.drop_duplicates("material_id").set_index("material_id")
    front_formulas = {t: sorted({names.loc[i, "formula"] for i in f}) for t, f in (("293K", fa), ("77K", fb))}
    only_known = {t: all(f in KNOWN for f in v) for t, v in front_formulas.items()}
    return {"n_tierA_both_T": len(m), "kendall_tau_A_fast_1E": round(float(tau), 4),
            "tau_ci95": [round(float(np.nanpercentile(bs, 2.5)), 4), round(float(np.nanpercentile(bs, 97.5)), 4)],
            "n_pareto_status_changes": len(changed),
            "changed": [f"{names.loc[i, 'formula']} {names.loc[i, 'sg_symbol']}" for i in changed],
            "front_formulas": front_formulas, "front_only_known_filters": only_known,
            "temperature_reorders_front": bool(len(changed) >= 1 or tau < TAU_REORDER),
            "rule": "reorders if >= 1 Pareto-status change or tau < 0.9 (PREREG_STAGE3 §6); decided in Stage 4"}


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    from filteratlas.atlas.assemble import FLAG_COLS
    for c in FLAG_COLS:
        at[c] = at[c].fillna(False).astype(bool)
    prox = pd.read_parquet(F("stage3_proxy.parquet"))[["material_id", "n_summary_A3"]]
    at = at.merge(prox, on="material_id", how="left")
    summ, mem = fronts(at)
    mem.to_parquet(F("stage3_fronts.parquet"))
    out = {"fronts": summ, "ladder": ladder(at), "who_moves": who_moves(at), "sensitivity": sensitivity(at),
           "kill3_preview": kill3_preview(at)}
    tb = {}   # Tier-B main front with `unscreenable` entries removed
    for tag, T in TAGS.items():
        d = at[(at["T"] == T) & (at.tier == "B") & at.main_front_eligible & ~at.unstable_dft]
        base, dd = front_ids(d, d.A_fast_1E, d.A_gamma_1MeV), d[~d.unscreenable]
        ex = front_ids(dd, dd.A_fast_1E, dd.A_gamma_1MeV)
        tb[tag] = {"changed": ex != base, "added": len(ex - base), "removed": len(base - ex),
                   "n_unscreenable_eligible": int(d.unscreenable.sum())}
    out["tierB_exclude_unscreenable"] = tb
    (F("stage3_fronts.json")).write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({"fronts": {k: (v["n_eligible"], v["n_front"], [f"{x['formula']} {x['sg_symbol']}" for x in v["front"]][:12])
                                 for k, v in summ.items()}, "kill3": out["kill3_preview"]}, indent=1, default=str))


if __name__ == "__main__":
    main()
