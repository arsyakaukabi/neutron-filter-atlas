"""Score the MLIPs against the DFT reference (MLIP benchmark, step 05, validation).

Every MLIP DOS (<mlip dir>/<model>/<material_id>.json.gz, MP DOS format, written by filteratlas.dynamics.mlip_run) is
converted with the identical code path as the MP DFT DOS (2 meV floor, NCrystal off-Bragg baseline) and compared with
the MP DFT values of the same compound. "dft_fc" (MP force constants through the same phonopy DOS settings) is the
DOS-method control: its error is the floor that no MLIP can beat.
Failures: run error, relaxation not converged, or a commensurate-q frequency < -0.5 THz.
-> the MLIP benchmark parquet table and summary JSON in data/derived
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from filteratlas.common.paths import D, F
from filteratlas.common.constants import N_WORKERS, TEMPS  # noqa: E402

THZ_TO_MEV = 4.135667696
from filteratlas.common.config import cfg  # noqa: E402

FAIL_FREQ = cfg("stability")["fc_unstable_THz"]   # failure threshold, fixed before the computation
RNG = np.random.default_rng(cfg("mlip")["benchmark_eval_seed"])


def mean_energy_meV(dos):
    f, g = np.asarray(dos["frequencies"], float), np.clip(np.asarray(dos["densities"], float), 0, None)
    m = f > 0
    return float(np.trapezoid(f[m] * g[m], f[m]) / np.trapezoid(g[m], f[m]) * THZ_TO_MEV)


def job(path):
    from filteratlas.ingest import mp as mq
    from filteratlas.xs import ncmat as nb
    from filteratlas.xs import evaluate as xe
    p = Path(path)
    dos = mq.load_dos(p)
    meta = dos.get("meta", {})
    rec = {"model": p.parent.name, "material_id": p.name.replace(".json.gz", ""), **{k: meta.get(k) for k in (
        "formula", "arm", "relax_converged", "relax_steps", "volume_ratio", "min_freq_mesh_THz", "min_freq_comm_THz",
        "n_disp", "wall_s")}}
    try:
        ev = xe.evaluate(nb.build_vdos_material(dos, f"ev_{rec['model']}_{rec['material_id']}"), TEMPS)
        rec.update({f"mlip_{k}": v for k, v in ev.items() if k.startswith(("L08_2", "L08_7", "mu_1"))})
        rec["mlip_mean_phonon_meV"] = mean_energy_meV(dos)
        rec["eval_status"] = "ok"
    except Exception as e:  # noqa: BLE001
        rec.update({"eval_status": "error", "eval_error": repr(e)[:200]})
    return rec


def p68_ci(x, B=2000):
    x = np.asarray(x, float)
    bs = [np.percentile(RNG.choice(x, len(x)), 68) for _ in range(B)]
    return round(float(np.percentile(x, 68)), 4), [round(float(np.percentile(bs, 2.5)), 4),
                                                   round(float(np.percentile(bs, 97.5)), 4)]


def main():
    files = sorted(str(p) for p in (D("mlip")).glob("*/*.json.gz"))
    with Pool(N_WORKERS) as pool:
        df = pd.DataFrame(pool.map(job, files, chunksize=4))
    # a compound in two arms (Al2O3 mp-aaaaabrz: random_1-10 and exfor_crystal) has one output file whose meta carries
    # only the arm it was first run under; take the arms from the set itself, one row per (compound, arm)
    arms = pd.DataFrame(json.loads((F("stage2_mlip_set.json")).read_text())["compounds"])[["material_id", "arm"]]
    df = df.drop(columns="arm").merge(arms.drop_duplicates(), on="material_id", how="left")
    # wall time is appended to the run log after the DOS file is written, so it is not in the file's meta
    logs = [json.loads(line) for m in df.model.unique() if (D("logs") / f"mlip_{m}.jsonl").exists()
            for line in open((D("logs") / f"mlip_{m}.jsonl"), encoding="utf-8")]
    wall = pd.DataFrame([r for r in logs if r.get("status") == "ok"])[["model", "material_id", "wall_s"]]
    df = df.drop(columns="wall_s").merge(wall.drop_duplicates(["model", "material_id"], keep="last"),
                                         on=["model", "material_id"], how="left")
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    ref = t[["material_id", "L08_293K", "L08_77K", "mu_1.0_293K", "mu_1.8_293K", "nsites", "phonon_method",
             "phonon_id"]].rename(columns=lambda c: c if c in ("material_id", "nsites", "phonon_method", "phonon_id")
                                  else f"dft_{c}")
    df = df.merge(ref, on="material_id", how="left")
    # DFT mean phonon energy from the MP DOS itself (same definition as for the MLIP DOS)
    from filteratlas.ingest import mp as mq
    df["dft_mean_phonon_meV"] = [mean_energy_meV(mq.load_dos(D("mp") / "phonon_dos" / m / f"{i}.json.gz"))
                                 if isinstance(i, str) else np.nan for m, i in zip(df.phonon_method, df.phonon_id)]
    df["failed"] = ((df.eval_status != "ok") | (df.relax_converged == False)  # noqa: E712
                    | (df.min_freq_comm_THz.fillna(0) < FAIL_FREQ))
    for tg in ("293K", "77K"):
        df[f"dlog_L08_{tg}"] = np.log(df[f"mlip_L08_{tg}"] / df[f"dft_L08_{tg}"])
    df["softening"] = df.mlip_mean_phonon_meV / df.dft_mean_phonon_meV
    df.to_parquet(F("stage2_mlip_benchmark.parquet"))

    res = {"fail_rule": "eval error | relax not converged | commensurate-q min freq < -0.5 THz", "models": {}}
    rnd = df[df.arm.fillna("").str.startswith("random")]
    for mname, g in rnd.groupby("model"):
        ok = g[~g.failed & g.dft_L08_293K.notna()]
        r = {"n_attempted": len(g), "n_failed": int(g.failed.sum()), "failure_rate": round(float(g.failed.mean()), 4),
             "n_scored": len(ok)}
        if len(ok) >= 10:
            for tg in ("293K", "77K"):
                a = ok[f"dlog_L08_{tg}"].abs()
                p68, ci = p68_ci(a)
                r[tg] = {"p68_abs_dlog": p68, "p68_ci95": ci, "median_abs_dlog": round(float(a.median()), 4),
                         "p90_abs_dlog": round(float(a.quantile(0.9)), 4),
                         "median_dlog": round(float(ok[f"dlog_L08_{tg}"].median()), 4),
                         "spearman_vs_dft": round(float(spearmanr(ok[f"mlip_L08_{tg}"], ok[f"dft_L08_{tg}"]).statistic), 4),
                         "top10_overlap": len(set(ok.nlargest(10, f"mlip_L08_{tg}").material_id)
                                              & set(ok.nlargest(10, f"dft_L08_{tg}").material_id))}
            r["softening_median"] = round(float(ok.softening.median()), 4)
            r["spearman_softening_vs_dlog293"] = round(float(spearmanr(ok.softening, ok.dlog_L08_293K).statistic), 4)
            r["wall_s_median"] = round(float(g.wall_s.median()), 1)
            r["wall_s_total_h"] = round(float(g.wall_s.sum() / 3600), 2)
            # stratification: by cell size and by presence of a first-row element (Li-F)
            from pymatgen.core import Composition
            first_row = ok.formula.map(lambda f: any(e.row == 2 for e in Composition(f).elements))   # Li..F
            strata = {"nsites<=10": ok.nsites <= 10, "nsites11-20": (ok.nsites > 10) & (ok.nsites <= 20),
                      "nsites21-40": ok.nsites > 20, "first_row": first_row, "no_first_row": ~first_row}
            r["p68_by_stratum_293K"] = {k: (round(float(ok[v].dlog_L08_293K.abs().quantile(0.68)), 4) if v.sum() >= 10
                                            else None) for k, v in strata.items()}
        res["models"][mname] = r
    # selection rule fixed before the computation (a recommendation; the atlas model is set in the config)
    cands = {m: r for m, r in res["models"].items() if m != "dft_fc" and "293K" in r and r["failure_rate"] <= 0.10}
    order = sorted(cands, key=lambda m: cands[m]["293K"]["p68_abs_dlog"])
    rec = order[0] if order else None
    if len(order) >= 2:
        a, b = cands[order[0]]["293K"], cands[order[1]]["293K"]
        if a["p68_ci95"][1] >= b["p68_ci95"][0]:   # CIs overlap -> Spearman, then wall time
            sa, sb = a["spearman_vs_dft"], b["spearman_vs_dft"]
            if abs(sa - sb) > 0.01:
                rec = order[0] if sa > sb else order[1]
            else:
                rec = min(order[:2], key=lambda m: cands[m]["wall_s_median"])
            res["decision_note"] = "best two p68 CIs overlap; tie broken by Spearman then wall time (B5)"
    res["excluded_failure_rate_gt_10pct"] = [m for m, r in res["models"].items()
                                             if m != "dft_fc" and r.get("failure_rate", 1) > 0.10]
    res["recommended_by_rule"] = rec
    exf = df[df.arm.fillna("").str.startswith("exfor")]
    res["exfor_crystals"] = exf[["model", "formula", "failed", "volume_ratio", "min_freq_comm_THz", "mlip_L08_293K",
                                 "dft_L08_293K", "dlog_L08_293K", "wall_s"]].round(4).to_dict("records")
    (F("stage2_mlip.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "exfor_crystals"}, indent=1, default=str))


if __name__ == "__main__":
    main()
