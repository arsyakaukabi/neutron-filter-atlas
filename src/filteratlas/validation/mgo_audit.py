"""MgO dataset conflict in the EXFOR benchmark: Harvey & Smith EXFOR 13761004 vs 13761005, both 'single crystal'.

Model-independent facts: sample text, stated target thickness (atoms/b) and the length it implies at the MP number
density, the free-atom ratio over 1-20 eV from the data audit, and the ratio of the two datasets per energy bin
(lower envelope).
Model tests: each dataset against the single-crystal off-Bragg model (lower envelope, as in the benchmark) and against
NCrystal's full polycrystal model incl. Bragg scattering (binned median), for MP DFT and curated. A dataset that
agrees with the polycrystal model and not with the off-Bragg one behaves as a polycrystal / multi-grain sample.
-> the MgO audit JSON in data/derived
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.validation import exfor_benchmark as S  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

SUBS = tuple(S.cfg("benchmark")["mgo_audit_subentries"])   # MgO 13761004 vs 13761005


def x4_sample(sub):
    txt = (D("exfor") / "13761.x4").read_text(errors="replace").splitlines()
    i0 = next(i for i, l in enumerate(txt) if l.startswith("SUBENT") and str(sub) in l)
    blk = []
    for l in txt[i0 + 1:]:
        if l.startswith(("ENDBIB", "ENDSUBENT")):
            break
        blk.append(l.rstrip())
    sample = [l[11:].strip() for l in blk if l.startswith("SAMPLE") or (blk and l.startswith(" " * 11) and "thickness" in l)]
    run = next((l.split("Run ID:")[1].strip() for l in blk if "Run ID:" in l), None)
    thick = next((float(m.group(1)) for l in blk for m in [re.search(r"([\d.]+) atoms/b", l)] if m), None)
    return {"sample_text": sample, "run_id": run, "thickness_atoms_per_b": thick}


def main():
    ds = {d["subentry"]: d for d in S.datasets() if d["subentry"] in SUBS}
    models, _ = S.build_models("MgO")
    _, summ = S.mp_entry("MgO", "Fm-3m")
    n_atoms_A3 = float(summ.nsites / summ.volume)            # atoms per A^3 from the MP structure
    out = {"note": "diagnostic; model-independent facts first", "mp_number_density_atoms_per_A3": round(n_atoms_A3, 5),
           "datasets": {}}
    env = {}
    for sub, d in ds.items():
        g = d["data"].sort_values("E_eV")
        E, D = g.E_eV.values, g.sigma_per_atom_b.values
        info = x4_sample(sub)
        t = info["thickness_atoms_per_b"]
        rec = {**info, "audit": d["audit"], "renorm_applied": d["renorm"], "free_atom_ratio": d["free_atom_ratio"],
               "models": {}}
        # t atoms/b = t * 1e24 atoms/cm^2; n = n_atoms_A3 * 1e24 atoms/cm^3  ->  L = t / n_atoms_A3  [cm]
        rec["implied_length_cm_if_atoms"] = round(t / n_atoms_A3, 3) if t else None
        rec["implied_length_cm_if_formula_units"] = round(2 * t / n_atoms_A3, 3) if t else None
        for mname in ("mp_dft", "curated"):
            for label, single in (("off_bragg_single_crystal", True), ("full_polycrystal_incl_bragg", False)):
                M = S.sigma(models[mname], d["T"], E, single)
                b = S.binned(E, D, M, single, 20 if single else 50)
                k = b[(b.E_c >= S.E_LO) & (b.E_c <= S.E_KILL)]
                lr = np.log(k.model / k.data)
                rec["models"][f"{mname}:{label}"] = {"r_d": round(float(np.sqrt(np.mean(lr ** 2))), 4),
                                                     "bias": round(float(lr.mean()), 4), "n_bins": len(k)}
        b = S.binned(E, D, D, True, S.ENV_Q)
        env[sub] = b.set_index("bin")
        out["datasets"][str(sub)] = rec
    a, c = env[SUBS[0]], env[SUBS[1]]
    com = a.index.intersection(c.index)
    out["ratio_004_over_005_lower_envelope"] = [{"E_meV": round(1e3 * a.loc[i, "E_c"], 3),
                                                 "ratio": round(float(a.loc[i, "data"] / c.loc[i, "data"]), 3)}
                                                for i in com]
    lo = [r["ratio"] for r in out["ratio_004_over_005_lower_envelope"] if r["E_meV"] <= 100]
    out["median_ratio_1_100meV"] = round(float(np.median(lo)), 3) if lo else None
    d4, d5 = (out["datasets"][str(x)] for x in SUBS)
    # epithermal (1-20 eV) ratio of the two datasets after the audit renormalisation, from their free-atom ratios
    out["epithermal_ratio_004_over_005_after_D9"] = round(d4["free_atom_ratio"] * d4["renorm_applied"]
                                                          / (d5["free_atom_ratio"] * d5["renorm_applied"]), 3)
    # sensitivity: 13761004 removed as a polycrystal-behaving sample -> MgO has no spread pair, threshold 0.15
    full = json.loads((F("stage2_exfor.json")).read_text())
    r5 = d5["models"]["mp_dft:off_bragg_single_crystal"]["r_d"]
    out["sensitivity_without_13761004"] = {"R_c_mp_dft": r5, "threshold": S.KILL_RMS, "fails": r5 > S.KILL_RMS,
                                           "threshold_with_004": full["crystals"]["MgO"]["threshold"]}
    best = {s: min(v["models"], key=lambda m: v["models"][m]["r_d"]) for s, v in out["datasets"].items()}
    out["best_model_per_dataset"] = best
    (F("stage2_mgo_audit.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
