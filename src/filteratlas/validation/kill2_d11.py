"""Acceptance criterion against measured total cross-sections, restricted to single-crystal data (step 05, validation).

The restriction to single-crystal data, the only sample form the model describes, was decided after the first
comparison. The test uses single-crystal datasets only (single_crystal, single_crystal_assumed, pyrolytic_graphite);
thresholds unchanged (R_c > max(0.15, S_c), S_c and thresholds taken from the benchmark summary JSON); the criterion
is not met at >= 6 failing single-crystal-evaluable crystals. Polycrystal / powder datasets are reported separately as
a validation of NCrystal's Bragg model, not as part of the criterion.
Reported side by side: the count over all sample forms, as fixed before the computation, and the single-crystal count,
for the atlas model (MP DFT where available, else the chosen MLIP) and, as scenarios, every candidate MLIP and the
curated NCrystal files; the narrower variant without pyrolytic graphite (as in filteratlas.validation.exfor_sc_only);
and the envelope sensitivity (10th and 50th percentile envelopes, thresholds held fixed).
-> the single-crystal acceptance JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402
from filteratlas.validation.exfor_benchmark import CHOSEN_MLIP, KILL_N, PC_FORMS, SC_FORMS  # noqa: E402

from filteratlas.validation.exfor_benchmark import ENV_Q  # noqa: E402

NARROW = SC_FORMS - {"pyrolytic_graphite"}


def rc_table(b, forms, q):
    g = b[b.qualifies & b.form.isin(forms) & (b.envelope_q == q)]
    return {(f, m): (round(float(np.sqrt(np.mean(x.r_d ** 2))), 4), int(x.subentry.nunique()))
            for (f, m), x in g.groupby(["formula", "model"])}


def count(full, rc, model):
    """Failing / evaluable crystals when `model` is the model for the crystals without MP DFT (as in the atlas), or
    for every crystal when model == 'curated'."""
    fails, evaluable, missing = [], [], []
    for f, c in full["crystals"].items():
        if not c["counts_for_kill"]:
            continue
        m = "curated" if model == "curated" else ("mp_dft" if "mp_dft" in c["per_model"] else model)
        if (f, m) not in rc:
            missing.append(f)
            continue
        evaluable.append(f)
        if rc[(f, m)][0] > c["threshold"]:
            fails.append(f)
    return {"failing": fails, "n_failing": len(fails), "evaluable": evaluable, "n_evaluable": len(evaluable),
            "no_data": missing, "fires": len(fails) >= KILL_N}


def main():
    full = json.loads((F("stage2_exfor.json")).read_text())
    b = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    b = b[b.get("error").isna()] if "error" in b else b
    atlas = f"mlip_{CHOSEN_MLIP}"
    cands = sorted(m for m in b.model.unique() if m.startswith("mlip_"))
    rc = rc_table(b, SC_FORMS, ENV_Q)
    out = {"note": "POST HOC (D11). Preregistered result: kill2 in results/stage2_exfor.json.",
           "atlas_mlip": atlas, "kill_n": KILL_N, "thresholds": {f: c["threshold"] for f, c in full["crystals"].items()},
           "preregistered": {k: full["kill2"][k] for k in ("failing", "n_failing", "evaluated_crystals", "fires_already")},
           "amended_D11": count(full, rc, atlas),
           "amended_D11_scenarios": {m: count(full, rc, m) for m in cands + ["curated"]},
           "narrow_no_pyrolytic_graphite": count(full, rc_table(b, NARROW, ENV_Q), atlas),
           "envelope_sensitivity": {}, "single_crystal_R_c": {}, "polycrystal_validation": {}}
    for q in (10, 20, 50):
        r = rc_table(b, SC_FORMS, q)
        out["envelope_sensitivity"][f"q{q}"] = {
            "atlas": count(full, r, atlas), "curated": count(full, r, "curated"),
            "R_c_atlas": {f: r.get((f, "mp_dft" if "mp_dft" in c["per_model"] else atlas), (None,))[0]
                          for f, c in full["crystals"].items() if c["counts_for_kill"]}}
    for (f, m), (R, n) in sorted(rc.items()):
        out["single_crystal_R_c"].setdefault(f, {"threshold": full["crystals"][f]["threshold"], "R_c": {}, "n_datasets": n})
        out["single_crystal_R_c"][f]["R_c"][m] = R
    for (f, m), (R, n) in sorted(rc_table(b, PC_FORMS, 50).items()):
        thr = full["crystals"][f]["threshold"]
        out["polycrystal_validation"].setdefault(f, {"threshold": thr, "n_datasets": n, "R_c": {}, "exceeds": []})
        out["polycrystal_validation"][f]["R_c"][m] = R
        if R > thr:
            out["polycrystal_validation"][f]["exceeds"].append(m)
    (F("stage2_kill2_D11.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("preregistered", "amended_D11", "amended_D11_scenarios",
                                          "narrow_no_pyrolytic_graphite")}, indent=1))
    for q, v in out["envelope_sensitivity"].items():
        print(q, v["atlas"]["failing"], v["R_c_atlas"])


if __name__ == "__main__":
    main()
