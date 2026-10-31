"""Sensitivity of the acceptance criterion against measured total cross-sections, NARROW variant (pyrolytic graphite
excluded); chosen after the first comparison.

The single-crystal restriction actually used includes pyrolytic graphite and is computed in its own step of step 05
(validation); this narrower variant is kept for comparison.

Recomputes the per-crystal RMS R_c using only single-crystal datasets (form single_crystal / single_crystal_assumed;
pyrolytic graphite, polycrystals and powders dropped), with the same envelope (q = 20), the same qualifying rule and the
same threshold max(0.15, S_c), fixed before the computation and taken from the benchmark summary JSON. Done for the
atlas model under each candidate MLIP (MP DFT where available) and, for reference, the curated NCrystal library files.
Motivation (stated before computing): the failing crystals' residuals are dominated by polycrystal / powder samples
whose small-angle and texture effects the off-Bragg single-crystal model does not describe.
-> the single-crystal-only sensitivity JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402

from filteratlas.common.config import cfg  # noqa: E402
from filteratlas.validation.exfor_benchmark import ENV_Q, KILL_N  # noqa: E402

SC_ONLY = set(cfg("benchmark")["single_crystal_forms_narrow"])   # single-crystal forms without pyrolytic graphite


def main():
    full = json.loads((F("stage2_exfor.json")).read_text())
    b = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    b = b[b.qualifies & (b.envelope_q == ENV_Q) & b.form.isin(SC_ONLY)]
    cands = sorted(m for m in b.model.unique() if m.startswith("mlip_"))
    out = {"note": "POST-HOC sensitivity, NARROW variant without pyrolytic graphite; the D11 amended count is results/stage2_kill2_D11.json; preregistered result is results/stage2_exfor.json",
           "thresholds_from": "stage2_exfor.json (max(0.15, S_c), all qualifying datasets)", "crystals": {}, "scenarios": {}}
    for f, c in full["crystals"].items():
        if not c["counts_for_kill"]:
            continue
        g = b[b.formula == f]
        out["crystals"][f] = {"threshold": c["threshold"], "n_sc_datasets": int(g.subentry.nunique()),
                              "R_c": {m: round(float(np.sqrt(np.mean(gm.r_d ** 2))), 4) for m, gm in g.groupby("model")}}
    for m in cands + ["curated"]:
        fails, missing = [], []
        for f, c in out["crystals"].items():
            atlas = "mp_dft" if "mp_dft" in full["crystals"][f]["per_model"] else m
            if m == "curated":
                atlas = "curated"
            R = c["R_c"].get(atlas)
            if R is None:
                missing.append(f)
            elif R > c["threshold"]:
                fails.append(f)
        out["scenarios"][m] = {"failing": fails, "n_failing": len(fails), "no_single_crystal_data": missing,
                               "n_evaluated": len(out["crystals"]) - len(missing), "fires": len(fails) >= KILL_N}
    (F("stage2_exfor_sc_only.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
