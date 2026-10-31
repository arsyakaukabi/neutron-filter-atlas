"""Phonon-input noise model from a large reproducibility set (rather than a handful of pairs) + SIMEX.

(a) pheasy (PBEsol, 2026) vs dfpt (PBE, ~2018) phonons for the SAME Materials Project compound (621 candidates
    have both). Their L0.8 spread includes functional and method differences, so it is a conservative estimate
    of the per-compound phonon-input uncertainty.
(b) enantiomorphic pairs and same-formula+space-group duplicates under the current VDOS conversion.
(c) SIMEX (Cook & Stefanski 1994): add noise of variance lambda*s^2 (lambda = 0.5..2) to log L0.8_full,
    extrapolate the mean Spearman rho to lambda = -1. Applied to the usable set of the triviality test and the
    filter-grade subset.
Output: stage1r_reproducibility.json, stage1r_pheasy_vs_dfpt.parquet
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from filteratlas.common.paths import D, F
from filteratlas.common.config import cfg
from filteratlas.common.constants import N_WORKERS, TEMPS, USABLE_MU_CUT  # noqa: E402

DOSDIR = D("mp") / "phonon_dos"
RNG = np.random.default_rng(cfg("analysis")["stage1"]["seed_reproducibility"])
FILTER_GRADE = cfg("analysis")["stage1"]["filter_grade_mu_cut_per_cm"]


def pair_job(args):
    mid, pid_p, pid_d = args
    from filteratlas.ingest import mp as mq
    from filteratlas.xs import ncmat as nb
    from filteratlas.xs import evaluate as xe
    out = {"material_id": mid}
    try:
        for meth, pid in (("pheasy", pid_p), ("dfpt", pid_d)):
            dos = mq.load_dos(DOSDIR / meth / f"{pid}.json.gz")
            w, ok = nb.stability(dos)
            out[f"stable_{meth}"] = ok
            ev = xe.evaluate(nb.build_vdos_material(dos, f"rep_{meth}_{mid}"), TEMPS)
            out.update({f"{meth}_{k}": v for k, v in ev.items() if k.startswith("L08")})
        out["status"] = "ok"
    except Exception as e:  # noqa: BLE001
        out["status"] = "error"
        out["error"] = repr(e)[:200]
    return out


def simex(x, y, s, lambdas=(0.5, 1.0, 1.5, 2.0), B=100):
    """SIMEX-corrected Spearman rho for noise sd s on log(x)."""
    lx = np.log(np.asarray(x, float))
    r0 = spearmanr(lx, y).statistic
    pts = [(0.0, r0)]
    for lam in lambdas:
        pts.append((lam, float(np.mean([spearmanr(lx + RNG.normal(0, np.sqrt(lam) * s, len(lx)), y).statistic
                                        for _ in range(B)]))))
    L, R = np.array(pts).T
    c = np.polyfit(L, R, 2)
    return {"rho_obs": round(float(r0), 4), "rho_simex": round(float(np.polyval(c, -1.0)), 4),
            "curve": [[round(a, 2), round(b, 4)] for a, b in pts]}


def main():
    cand = pd.read_parquet(F("stage0_candidates.parquet"))
    both = cand[(cand.phonon_pheasy != "") & (cand.phonon_dfpt != "")]
    jobs = [(r.material_id, r.phonon_pheasy.split(",")[0], r.phonon_dfpt.split(",")[0]) for r in both.itertuples()]
    with Pool(N_WORKERS) as p:
        rows = p.map(pair_job, jobs, chunksize=4)
    pdf = pd.DataFrame(rows).merge(both[["material_id", "formula", "sg_symbol"]], on="material_id")
    pdf.to_parquet(F("stage1r_pheasy_vs_dfpt.parquet"))
    ok = pdf[(pdf.status == "ok") & pdf.stable_pheasy.fillna(False) & pdf.stable_dfpt.fillna(False)]
    st = pd.read_parquet(F("stage1_triviality.parquet")).set_index("material_id")
    usable_ids = set(st.index[(st.status == "ok") & st.stable.fillna(False) & (st["mu_1.8_293K"] < USABLE_MU_CUT)])
    res = {"pheasy_vs_dfpt": {"n_pairs_computed": int((pdf.status == "ok").sum()), "n_both_stable": len(ok),
                              "n_both_stable_and_usable": int(ok.material_id.isin(usable_ids).sum())}}
    for t in ("293K", "77K"):
        okk = ok[ok.material_id.isin(usable_ids)]
        d = np.log(okk[f"pheasy_L08_{t}"] / okk[f"dfpt_L08_{t}"])
        res["pheasy_vs_dfpt"][t] = {"median_dlog": round(float(d.median()), 4),
                                    "median_abs_dlog": round(float(d.abs().median()), 4),
                                    "p68_abs_dlog": round(float(d.abs().quantile(0.68)), 4),
                                    "p90_abs_dlog": round(float(d.abs().quantile(0.90)), 4),
                                    "sd_dlog": round(float(d.std()), 4),
                                    "per_entry_sigma_est": round(float(d.std() / np.sqrt(2)), 4),
                                    "spearman_pheasy_vs_dfpt": round(float(spearmanr(okk[f"pheasy_L08_{t}"],
                                                                                     okk[f"dfpt_L08_{t}"]).statistic), 4)}

    # (b) enantiomorphs and duplicates under the current conversion
    stt = st[(st.status == "ok") & st.stable.fillna(False)]
    en = json.loads((F("stage1_enantiomorphs.json")).read_text())
    pairs = [(p["id_a"], p["id_b"]) for p in en.get("pairs", [])]
    res["enantiomorph_pairs"] = {}
    for t in ("293K", "77K"):
        d = [abs(np.log(stt.loc[a, f"L08_{t}"] / stt.loc[b, f"L08_{t}"])) for a, b in pairs if a in stt.index and b in stt.index]
        res["enantiomorph_pairs"][t] = {"n": len(d), "median_abs_dlog": round(float(np.median(d)), 4),
                                        "max_abs_dlog": round(float(np.max(d)), 4)}
    g = stt.reset_index().groupby(["formula", "sg_symbol"])
    spreads = {t: [float(np.log(x[f"L08_{t}"].max() / x[f"L08_{t}"].min())) for _, x in g if len(x) > 1]
               for t in ("293K", "77K")}
    res["same_formula_sg_duplicates"] = {t: {"n_groups": len(v), "median_log_spread": round(float(np.median(v)), 4),
                                             "max_log_spread": round(float(np.max(v)), 4)} for t, v in spreads.items()}

    # (c) SIMEX with the conservative per-entry sigma from (a)
    u = stt[stt["mu_1.8_293K"] < USABLE_MU_CUT]
    fg = stt[stt["mu_1.8_293K"] < FILTER_GRADE]
    res["simex"] = {}
    for t in ("293K", "77K"):
        s = res["pheasy_vs_dfpt"][t]["per_entry_sigma_est"]
        res["simex"][t] = {"sigma": s, "primary_usable_0.3": simex(u[f"L08_{t}"], u[f"L08_proxy_{t}"], s),
                           "filter_grade_0.1": simex(fg[f"L08_{t}"], fg[f"L08_proxy_{t}"], s)}
    (F("stage1r_reproducibility.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
