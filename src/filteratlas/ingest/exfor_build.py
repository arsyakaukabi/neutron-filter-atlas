"""Fetch + parse EXFOR filter-crystal total cross sections into a tidy table (step 00).

Output: exfor_benchmark.csv (one row per data point) and stage0_exfor_summary.json
(one record per subentry with material, sample form, temperature, energy range, n points).
The sample form is classified from the SAMPLE/TITLE text by the rules in FORM_RULES (checked by hand);
'room temperature' without a number is recorded as 293 K with temp_assumed=True.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

from filteratlas.ingest import exfor as exfor  # noqa: E402
from filteratlas.common.config import cfg
from filteratlas.common.paths import F

_B = cfg("benchmark")


def entries():
    """Entry list fixed before the computation plus the verified extra entries (curated/exfor_candidates_extra.csv)."""
    ents = list(_B["exfor_entries"])
    extra = F("exfor_candidates_extra.csv")   # entries found in a literature search (verified ones are used)
    if extra.exists():
        ex = pd.read_csv(extra)
        ents += [int(e) for e in ex.loc[ex["verified"].astype(str).str.lower().str.startswith("yes"), "entry"]]
    return sorted(set(ents))


TARGET = {k: tuple(v) for k, v in _B["target_codes"].items()}   # EXFOR target code -> (formula, atoms per formula unit)
# Per-subentry overrides where the EXFOR target code would mislabel the material (checked against SAMPLE text):
OVERRIDE = {k: tuple(v) for k, v in _B["subentry_overrides"].items()}
FORM_RULES = [tuple(r) for r in _B["sample_form_rules"]]


def main():
    rows, summ = [], []
    for e in entries():
        try:
            df = exfor.parse(e)
        except Exception as ex_:  # noqa: BLE001
            summ.append({"entry": e, "error": repr(ex_)[:200]})
            continue
        if df.empty:
            summ.append({"entry": e, "error": "no parseable (N,TOT) DATA block (energy/sigma units not recognised)"})
            continue
        for (sub, rx), g in df.groupby(["subentry", "reaction"]):
            code = re.match(r"\(([^()]+)\(N,TOT\)", rx)
            code = code.group(1) if code else ""
            formula, apf = TARGET.get(code, (None, None))
            note = ""
            if sub in OVERRIDE:
                formula, apf, note = OVERRIDE[sub]
            text = (g["sample"].iloc[0] + " " + g["title"].iloc[0]).lower()
            form = next((f for k, f in FORM_RULES if k in text), "unknown")
            temp = g.temp_K.dropna()
            t_val, assumed = (float(temp.iloc[0]), False) if len(temp) else (293.0, True)
            if "room temp" not in text and not len(temp):
                assumed = True
            derived = "DERIV" in rx.upper() or ",AV)" in rx.upper()
            rec = {"entry": e, "subentry": sub, "target_code": code, "formula": formula, "atoms_per_formula": apf,
                   "sample_form": form, "temp_K": t_val, "temp_assumed": assumed, "n": len(g),
                   "E_min_eV": float(g.E_eV.min()), "E_max_eV": float(g.E_eV.max()),
                   "derived_or_average": derived, "reference": g.reference.iloc[0], "author": g.author.iloc[0][:80],
                   "sample_text": g["sample"].iloc[0][:200], "status": g.status.iloc[0][:80],
                   "reaction": rx[:120], "override_note": note}
            summ.append(rec)
            if formula and not derived:
                gg = g.copy()
                gg["sigma_per_atom_b"] = gg.sigma_b / apf
                gg["err_per_atom_b"] = gg.err_b / apf
                for k in ("target_code", "formula", "atoms_per_formula", "sample_form", "temp_assumed"):
                    gg[k] = rec[k]
                gg["temp_K"] = gg.temp_K.fillna(t_val)
                rows.append(gg)
    tidy = pd.concat(rows, ignore_index=True)
    cols = ["entry", "subentry", "formula", "target_code", "sample_form", "temp_K", "temp_assumed", "E_eV", "sigma_b",
            "err_b", "atoms_per_formula", "sigma_per_atom_b", "err_per_atom_b", "reference", "author", "status"]
    tidy[cols].to_csv(F("exfor_benchmark.csv"), index=False)
    (F("stage0_exfor_summary.json")).write_text(json.dumps(summ, indent=1, default=str))
    s = pd.DataFrame([r for r in summ if "error" not in r])
    print(s[s.formula.notna()].groupby(["formula", "sample_form"]).agg(subentries=("subentry", "count"),
                                                                        points=("n", "sum")).to_string())
    print("errors:", [r for r in summ if "error" in r])


if __name__ == "__main__":
    main()
