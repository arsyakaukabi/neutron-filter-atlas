"""Pb / Bi diagnostics of the EXFOR benchmark (DIAGNOSTIC ONLY; the atlas keeps the static model fixed before the
computation).

Against the qualifying single-crystal datasets (q20 envelope, 1-100 meV, as in the benchmark):
  - SevenNet static (atlas), curated NCrystal, elastic Debye (from the per dataset x model benchmark table);
  - SevenNet VDOS with all frequencies scaled by s = 0.80 ... 1.20 (step 0.02): which s minimises R_c, and does any
    s in the 10-20 % range bring R_c below the threshold;
  - SevenNet quasi-harmonic VDOS at the dataset temperature (filteratlas.validation.pbbi_qha; used if the files exist).
-> the Pb/Bi diagnostic JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.validation import exfor_benchmark as S  # noqa: E402
from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

CRYST = ("Pb", "Bi")
SCALES = np.round(np.arange(0.80, 1.2001, 0.02), 2)
THZ_TO_MEV = 4.135667696


def r_d(cfg, d):
    g = d["data"].sort_values("E_eV")
    E, D = g.E_eV.values, g.sigma_per_atom_b.values
    b = S.binned(E, D, S.sigma(cfg, d["T"], E, True), True, S.ENV_Q)
    k = b[(b.E_c >= S.E_LO) & (b.E_c <= S.E_KILL)]
    lr = np.log(k.model / k.data)
    return float(np.sqrt(np.mean(lr ** 2))), float(lr.mean())


def mean_meV(dos):
    f, g = np.asarray(dos["frequencies"], float), np.clip(np.asarray(dos["densities"], float), 0, None)
    m = f > 0
    return float(np.trapezoid(f[m] * g[m], f[m]) / np.trapezoid(g[m], f[m]) * THZ_TO_MEV)


def main():
    full = json.loads((F("stage2_exfor.json")).read_text())
    bench = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    ds = [d for d in S.datasets() if d["formula"] in CRYST and d["qualifies"] and d["form"] in S.SC_FORMS]
    out = {"note": "diagnostic only", "crystals": {}}
    for f in CRYST:
        mid = full["mp_material_ids"].get(f) or S.mp_entry(f, S.CRYSTALS[f][0])[1].material_id
        dos = mq.load_dos(D("mlip") / S.CHOSEN_MLIP / f"{mid}.json.gz")
        dd = [d for d in ds if d["formula"] == f]
        thr = full["crystals"][f]["threshold"]
        rec = {"material_id": mid, "threshold": thr, "datasets": [[d["subentry"], d["T"]] for d in dd],
               "mean_phonon_meV_static": round(mean_meV(dos), 3), "reference_models": {}, "scale_scan": []}
        sub = bench[bench.formula.eq(f) & bench.subentry.isin([d["subentry"] for d in dd]) & (bench.envelope_q == 20)]
        for m in ("curated", "debye_elastic", f"mlip_{S.CHOSEN_MLIP}"):
            x = sub[sub.model == m]
            if len(x):
                rec["reference_models"][m] = {"R_c": round(float(np.sqrt(np.mean(x.r_d ** 2))), 4),
                                              "mean_bias": round(float(x.bias.mean()), 4)}
        for s in SCALES:
            d2 = {**dos, "frequencies": (np.asarray(dos["frequencies"], float) * s).tolist()}
            cfg = nb.build_vdos_material(d2, f"pbbi_{f}_s{int(round(s * 100))}")
            rs = [r_d(cfg, d) for d in dd]
            rec["scale_scan"].append({"s": float(s), "R_c": round(float(np.sqrt(np.mean([r[0] ** 2 for r in rs]))), 4),
                                      "mean_bias": round(float(np.mean([r[1] for r in rs])), 4)})
        best = min(rec["scale_scan"], key=lambda r: r["R_c"])
        rec["best_scale"] = best
        rec["scale_closes_gap"] = [r["s"] for r in rec["scale_scan"] if r["R_c"] <= thr]
        qrs, qmeta = [], {}
        for d in dd:
            p = D("mlip_sensitivity") / f"{S.CHOSEN_MLIP}_qha{int(d['T'])}K" / f"{mid}.json.gz"
            if p.exists():
                qd = mq.load_dos(p)
                qrs.append(r_d(nb.build_vdos_material(qd, f"pbbi_{f}_qha{int(d['T'])}"), d))
                qmeta[str(int(d["T"]))] = {"V_over_V0": round(qd["meta"]["V_over_V0"], 5),
                                           "mean_phonon_meV": round(mean_meV(qd), 3),
                                           "mean_phonon_ratio_to_static": round(mean_meV(qd) / mean_meV(dos), 4)}
        if qrs and len(qrs) == len(dd):
            rec["qha"] = {"R_c": round(float(np.sqrt(np.mean([r[0] ** 2 for r in qrs]))), 4),
                          "mean_bias": round(float(np.mean([r[1] for r in qrs])), 4),
                          "closes_gap": bool(np.sqrt(np.mean([r[0] ** 2 for r in qrs])) <= thr), "per_T": qmeta}
        else:
            rec["qha"] = "not available (QHA DOS files missing)"
        out["crystals"][f] = rec
    (F("stage2_pbbi_diag.json")).write_text(json.dumps(out, indent=1))
    for f, r in out["crystals"].items():
        print(f, "thr", r["threshold"], r["reference_models"], "best", r["best_scale"], "closes", r["scale_closes_gap"],
              "qha", r["qha"])


if __name__ == "__main__":
    main()
