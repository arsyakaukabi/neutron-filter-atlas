"""Check added after the first comparison: reproducibility of MP DFT phonon DOS, from enantiomorphic pairs.

Enantiomorphs (e.g. P3_121 / P3_221 alpha-quartz) are physically identical crystals for neutron scattering, so
any difference in our figures of merit between the two MP entries measures the noise of the MP phonon DOS +
our pipeline. Pairs: same reduced formula, enantiomorphic space-group pair, lattice parameters within 1 %.
Also records the DOS weight below 5 meV, the suspected driver. -> the enantiomorph JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import THZ_TO_EV  # noqa: E402

PAIRS = [("P3_121", "P3_221"), ("P3_112", "P3_212"), ("P6_122", "P6_522"), ("P6_222", "P6_422"),
         ("P6_1", "P6_5"), ("P6_2", "P6_4"), ("P4_122", "P4_322"), ("P4_12_12", "P4_32_12"), ("P4_1", "P4_3"),
         ("P4_132", "P4_332"), ("P3_1", "P3_2")]


def low_weight(row, emax_meV=5.0):
    d = mq.load_dos(D("mp") / "phonon_dos" / row.phonon_method / f"{row.phonon_id}.json.gz")
    f, g = np.asarray(d["frequencies"]), np.clip(np.asarray(d["densities"]), 0, None)
    e = f * THZ_TO_EV * 1e3
    return float(np.trapezoid(np.where((e > 0) & (e < emax_meV), g, 0), f) / np.trapezoid(g, f))


def abc(row):
    from pymatgen.core import Structure
    d = mq.load_dos(D("mp") / "phonon_dos" / row.phonon_method / f"{row.phonon_id}.json.gz")
    s = Structure.from_dict(d["structure"])
    return np.array(sorted(s.lattice.abc)), s.volume / len(s)


def main():
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    t = t[(t.status == "ok") & t.stable.fillna(False)]
    rows = []
    for a, b in PAIRS:
        if not b:
            continue
        A, B = t[t.sg_symbol == a], t[t.sg_symbol == b]
        for f in sorted(set(A.formula) & set(B.formula)):
            for ra in A[A.formula == f].itertuples():
                la, va = abc(ra)
                for rb in B[B.formula == f].itertuples():
                    lb, vb = abc(rb)
                    if np.max(np.abs(la / lb - 1)) < 0.01 and abs(va / vb - 1) < 0.01:
                        rec = {"formula": f, "sg_a": a, "sg_b": b, "id_a": ra.material_id, "id_b": rb.material_id,
                               "lowE_weight_a": low_weight(ra), "lowE_weight_b": low_weight(rb)}
                        for tag in ("293K", "77K"):
                            rec[f"L08_a_{tag}"] = getattr(ra, f"L08_{tag}")
                            rec[f"L08_b_{tag}"] = getattr(rb, f"L08_{tag}")
                            rec[f"abs_log_ratio_{tag}"] = abs(np.log(rec[f"L08_a_{tag}"] / rec[f"L08_b_{tag}"]))
                        rows.append(rec)
    df = pd.DataFrame(rows)
    res = {"n_pairs": len(df)}
    if len(df):
        for tag in ("293K", "77K"):
            x = df[f"abs_log_ratio_{tag}"]
            res[tag] = {"median_abs_log_ratio": round(float(x.median()), 4), "p68": round(float(x.quantile(0.68)), 4),
                        "p90": round(float(x.quantile(0.90)), 4), "max": round(float(x.max()), 4),
                        "frac_gt_0.2": round(float((x > 0.2).mean()), 3)}
        dw = np.log(df.lowE_weight_a.clip(1e-5) / df.lowE_weight_b.clip(1e-5))
        dl = np.log(df.L08_a_293K / df.L08_b_293K)
        from scipy.stats import spearmanr
        res["spearman_dlog_lowE_weight_vs_dlog_L08_293K"] = round(float(spearmanr(dw, dl).statistic), 3)
        res["pairs"] = df.round(4).to_dict("records")
    (F("stage1_enantiomorphs.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "pairs"}, indent=1))
    print(df[["formula", "sg_a", "L08_a_293K", "L08_b_293K", "lowE_weight_a", "lowE_weight_b"]].round(4).to_string()
          if len(df) else "no pairs")


if __name__ == "__main__":
    main()
