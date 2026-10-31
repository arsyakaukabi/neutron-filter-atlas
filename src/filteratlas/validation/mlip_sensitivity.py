"""MLIP sensitivity runs (the mlip_sensitivity data directory), added after inspection of the MLIP benchmark and
scored with the identical code path as the benchmark (filteratlas.validation.mlip_eval).

1. ORB v3 with 0.03 A finite displacements (the benchmark fixes 0.01 A): failure rate (same failure rule) and
   p68 |dln L0.8| vs DFT on the random 150, next to the benchmark ORB run. Does not change the MLIP choice.
2. SevenNet on GPU (CUDA; preparation of the atlas run): per-compound reproducibility vs the CPU run (|dln L0.8| GPU
   vs CPU, agreement of the failure flag) and the wall-time ratio from the run logs.
-> the MLIP sensitivity JSON in data/derived
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F
from filteratlas.common.constants import N_WORKERS  # noqa: E402
from filteratlas.validation.mlip_eval import FAIL_FREQ, job, p68_ci  # noqa: E402

LOGS = D("logs")
RUNS = {"orb_disp030": ("orb_v3_cons_inf_omat_disp030", "orb_v3_cons_inf_omat"),
        "sevennet_cuda": ("sevennet_mf_ompa_cuda", "sevennet_mf_ompa")}


def wall(tag):
    f = LOGS / f"mlip_{tag}.jsonl"
    if not f.exists():
        return pd.DataFrame(columns=["material_id", "wall_s"])
    rows = [json.loads(x) for x in open(f, encoding="utf-8") if x.strip()]
    return pd.DataFrame([r for r in rows if r.get("status") == "ok"])[["material_id", "wall_s"]].drop_duplicates(
        "material_id", keep="last")


def failed(d):
    return (d.eval_status != "ok") | (d.relax_converged == False) | (d.min_freq_comm_THz.fillna(0) < FAIL_FREQ)  # noqa: E712


def main():
    base = pd.read_parquet(F("stage2_mlip_benchmark.parquet"))
    arms = pd.DataFrame(json.loads((F("stage2_mlip_set.json")).read_text())["compounds"])[["material_id", "arm"]]
    out = {"note": "post hoc / Stage 3 preparation; the preregistered Stage 2 result is results/stage2_mlip.json"}
    for key, (tag, ref) in RUNS.items():
        files = sorted(str(p) for p in (D("mlip_sensitivity") / tag).glob("*.json.gz"))
        if not files:
            out[key] = "no outputs"
            continue
        with Pool(N_WORKERS) as pool:
            d = pd.DataFrame(pool.map(job, files, chunksize=4))
        # wall time is not in the DOS meta (it is logged after the file is written); take it from the run log below
        d = d.drop(columns=["arm", "wall_s"]).merge(arms.drop_duplicates(), on="material_id", how="left")
        d["failed"] = failed(d)
        r = base[base.model == ref].drop_duplicates(["material_id", "arm"])
        m = d.merge(r[["material_id", "arm", "failed", "mlip_L08_293K", "mlip_L08_77K", "dft_L08_293K", "wall_s"]],
                    on=["material_id", "arm"], suffixes=("", "_ref"))
        m = m.merge(wall(tag), on="material_id", how="left", suffixes=("", "_new"))
        rnd = m[m.arm.str.startswith("random")]
        rec = {"n_outputs": len(d), "n_matched": len(m), "random_n": len(rnd),
               "random_failed": int(rnd.failed.sum()), "random_failed_ref": int(rnd.failed_ref.sum()),
               "failure_flag_agreement": round(float((m.failed == m.failed_ref).mean()), 4)}
        ok = rnd[~rnd.failed & rnd.dft_L08_293K.notna()]
        if len(ok) >= 10:
            p, ci = p68_ci(np.abs(np.log(ok.mlip_L08_293K / ok.dft_L08_293K)))
            rec["random_p68_vs_dft_293K"] = {"p68": p, "ci95": ci, "n": len(ok)}
        both = m[~m.failed & ~m.failed_ref & m.mlip_L08_293K_ref.notna()]
        if len(both):
            dl = np.abs(np.log(both.mlip_L08_293K / both.mlip_L08_293K_ref))
            rec["vs_reference_run_abs_dlnL08_293K"] = {"median": round(float(dl.median()), 5),
                                                       "p68": round(float(dl.quantile(0.68)), 5),
                                                       "max": round(float(dl.max()), 5), "n": len(both)}
        if "wall_s_new" in m and m.wall_s_new.notna().any():
            w = m[m.wall_s_new.notna() & m.wall_s.notna()]
            rec["wall_time"] = {"total_h_new": round(float(w.wall_s_new.sum() / 3600), 2),
                                "total_h_ref": round(float(w.wall_s.sum() / 3600), 2),
                                "median_ratio_new_over_ref": round(float((w.wall_s_new / w.wall_s).median()), 3), "n": len(w)}
        out[key] = rec
    (F("stage2_mlip_sensitivity.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
