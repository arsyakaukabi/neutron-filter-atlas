"""Input of the OpenMC check of the atlas: sapphire, Si and the top Tier-A compound of the 293 K main front by
A_fast,1/E (or the next one if that is sapphire or Si), with the atlas density, L0.8(293 K) and analytic A_fast.
-> the OpenMC-check input JSON in data/derived (read by filteratlas.transport.openmc_check in the OpenMC environment)
"""
import json

import pandas as pd

from filteratlas.common.paths import F
from filteratlas.fom import rejection as rj

LANDMARKS = (("sapphire", "Al2O3", "R-3c"), ("Si", "Si", "Fd-3m"))


def row_of(at, mid):
    return at[(at.material_id == mid) & (at["T"] > 100)].iloc[0]


def record(name, r):
    return {"name": name, "formula": r.formula, "material_id": r.material_id, "fractions": rj.fractions(r.formula),
            "n_A3": float(r.numberdensity_A3), "L_cm": float(r.L08), "A_fast_1E": float(r.A_fast_1E),
            "A_fast_watt": float(r.A_fast_watt)}


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    one = at[(at["T"] > 100) & (at.tier == "A")]
    out, used = [], set()
    for name, f, sg in LANDMARKS:
        r = one[(one.formula == f) & (one.sg_symbol == sg)].iloc[0]
        out.append(record(name, r))
        used.add(r.material_id)
    fronts = json.loads(F("stage3_fronts.json").read_text())["fronts"]["293K|A|main|main"]["front"]
    top = [m for m in sorted(fronts, key=lambda x: -x["A_fast_1E"]) if m["material_id"] not in used][0]
    r = row_of(at, top["material_id"])
    out.append(record(f"top_tierA_{r.formula}", r))
    F("stage3_openmc_input.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
