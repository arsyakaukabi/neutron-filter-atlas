"""Tiers: map data/curated/growable_crystals.csv onto candidate material_ids (Tier A).

Rule (fixed before any ranking): reduced-formula match; keep candidates whose MP crystal system
matches the row's crystal_system_or_phase text; a stated Hermann-Mauguin symbol that matches MP sg_symbol wins
(tierA_match = sg); else a structure-type keyword (hcp, rutile, CsCl type, ...) mapped to its space group(s)
(structure_type); else the lowest-E_hull crystal-system match (crystal_system_lowest_ehull); else the lowest-E_hull
entry of the formula, flagged tierA_phase_mismatch. Solid solutions are listed as unmapped; growable formulas absent
from the candidate set are listed with the reason.
-> the tier table (one row per Tier-A material_id) and its JSON summary in data/derived
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

from filteratlas.common.paths import D, F
from filteratlas.common.constants import E_HULL_MAX, RADIOACTIVE, STRONG_ABSORBERS  # noqa: E402

from filteratlas.common.config import cfg

_T = cfg("funnel")["tier_a_matching"]
SYSTEMS = tuple(_T["crystal_systems"])
# structure-type keywords -> MP space-group symbols (standard crystallography). Order matters: longer phrases first.
STRUCTURE_TYPES = [(kw, set(sgs)) for kw, sgs in _T["structure_types"]]
SOLID_SOLUTION = re.compile(_T["solid_solution_regex"])


def systems_in(text):
    t = text.lower()
    s = {w for w in SYSTEMS if w in t}
    if "rhombohedral" in s:
        s = (s - {"rhombohedral"}) | {"trigonal"}
    if "trigonal/hexagonal" in t:
        s |= {"trigonal", "hexagonal"}
    return {x.capitalize() for x in s}


def type_symbols(text):
    t = text.lower()
    out = set()
    for kw, sgs in STRUCTURE_TYPES:
        if re.search(r"(?<![a-z0-9])" + re.escape(kw) + r"(?![a-z0-9])", t):
            out |= sgs
            t = t.replace(kw, " ")
    return out


def hm_symbols(text):
    """Hermann-Mauguin symbols written in the text, normalised to MP style (P3121 -> P3_121, I41/a -> I4_1/a)."""
    out = set()
    for tok in re.findall(r"\b([PIFCRA][-0-9a-z_]{1,10}(?:/[a-z0-9]+)?)", text):   # "P3121/P3221" -> two symbols
        if not re.search(r"\d", tok):
            continue
        t = re.sub(r"^([PIFCRA])(-?)([2346])([1-5])(?=[a-z0-9/]|$)", r"\1\2\3_\4", tok)
        out |= {tok, t}
    return out


def main():
    g = pd.read_csv(F("growable_crystals.csv"))
    c = pd.read_parquet(F("stage0_candidates.parquet"))
    c["reduced"] = [Composition(f).reduced_formula for f in c.formula]
    rows, unmapped, notcand = [], [], []
    for f, grp in g.groupby("formula"):
        text = " | ".join(sorted(set(grp.crystal_system_or_phase.astype(str))))
        notes = " | ".join(sorted(set(grp.notes.dropna().astype(str))))[:300]
        if SOLID_SOLUTION.search(f):
            unmapped.append({"formula": f, "phase_text": text, "reason": "unmapped_solid_solution"})
            continue
        try:
            comp = Composition(f)
        except Exception:  # noqa: BLE001
            unmapped.append({"formula": f, "phase_text": text, "reason": "unparseable_formula"})
            continue
        red = comp.reduced_formula
        cand = c[c.reduced == red]
        if cand.empty:
            els = {str(e) for e in comp.elements}
            why = ("H-bearing" if "H" in els else "excluded absorber " + ",".join(sorted(els & STRONG_ABSORBERS))
                   if els & STRONG_ABSORBERS else "radioactive" if els & RADIOACTIVE
                   else f"no observed MP entry with E_hull <= {E_HULL_MAX} eV/atom")
            notcand.append({"formula": f, "reduced": red, "reason": why})
            continue
        sysm = systems_in(text)
        hm = hm_symbols(text)
        bys = cand[cand.crystal_system.isin(sysm)] if sysm else cand.iloc[0:0]
        sg = bys[bys.sg_symbol.isin(hm)] if hm else bys.iloc[0:0]
        st = bys[bys.sg_symbol.isin(type_symbols(text))]
        graphite = red == "C" and "graphite" in text.lower()
        diamond = red == "C" and "diamond" in text.lower()
        picks = []
        if red == "C":   # two distinct Tier-A materials share the formula: handle each phase text separately
            if diamond:
                d = cand[cand.sg_symbol == "Fd-3m"].sort_values("e_hull").head(1)
                picks += [(r, "structure_type", "diamond") for r in d.itertuples()]
                if d.empty:
                    notcand.append({"formula": "C (diamond)", "reduced": red,
                                    "reason": f"phase not in candidate set (no Fd-3m C with E_hull <= {E_HULL_MAX} eV/atom)"})
            if graphite:
                h = cand[cand.crystal_system == "Hexagonal"].sort_values("e_hull").head(1)
                picks += [(r, "crystal_system_lowest_ehull", "graphite") for r in h.itertuples()]
        elif len(sg):   # one entry per matched space group (lowest E_hull); enantiomorph pairs keep both
            picks = [(r, "sg", "") for r in sg.sort_values("e_hull").groupby("sg_symbol").head(1).itertuples()]
        elif len(st):
            picks = [(r, "structure_type", "") for r in st.sort_values("e_hull").groupby("sg_symbol").head(1).itertuples()]
        elif len(bys):
            picks = [(r, "crystal_system_lowest_ehull", "") for r in bys.sort_values("e_hull").head(1).itertuples()]
        else:
            picks = [(r, "lowest_ehull_phase_mismatch", "") for r in cand.sort_values("e_hull").head(1).itertuples()]
        for r, how, sub in picks:
            rows.append({"material_id": r.material_id, "growable_formula": f, "formula": r.formula,
                         "sg_symbol": r.sg_symbol, "crystal_system": r.crystal_system, "e_hull": r.e_hull,
                         "phase_text": text, "tierA_match": how, "n_vendor_rows": len(grp),
                         "n_polymorph_candidates": len(cand),
                         "tierA_phase_mismatch": how == "lowest_ehull_phase_mismatch",
                         "not_off_bragg_single_crystal": sub == "graphite", "notes": notes})
    t = pd.DataFrame(rows).sort_values("material_id").reset_index(drop=True)
    t.to_parquet(F("stage3_tiers.parquet"))
    out = {"n_growable_rows": len(g), "n_growable_formulas": int(g.formula.nunique()),
           "n_tierA_material_ids": len(t), "n_tierA_formulas_mapped": int(t.growable_formula.nunique()),
           "match_counts": t.tierA_match.value_counts().to_dict(),
           "phase_mismatch": t[t.tierA_phase_mismatch][["growable_formula", "sg_symbol", "phase_text"]].to_dict("records"),
           "unmapped_solid_solution": unmapped, "growable_not_candidate": notcand}
    (F("stage3_tiers.json")).write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: v for k, v in out.items() if k != "growable_not_candidate"}, indent=1, default=str))
    print("not candidate:", len(notcand), pd.DataFrame(notcand).reason.str[:20].value_counts().to_dict())


if __name__ == "__main__":
    main()
