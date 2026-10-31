"""Inputs of the OpenMC slab-transport runs of the shortlist, written in the main environment (OpenMC runs in its own).

For every shortlisted compound: its NCMAT (atlas structure and per-element VDOS: the MP DFT DOS or the SevenNet
output, as in the atlas), and per temperature the slab thicknesses (L0.8, and 5 cm where L0.8 < 5 cm) with the
analytic uncollided A_fast (1/E, Watt), A_gamma (1, 6 MeV) and T_th at each thickness, from the atlas code.
-> intermediate/stage4_transport/ncmat/<material_id>.ncmat, intermediate/stage4_transport/cases.json
"""
import json

import numpy as np
import pandas as pd

from filteratlas.atlas.assemble import MLIP_DIR
from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F
from filteratlas.fom import rejection as rj
from filteratlas.ingest import mp as mq
from filteratlas.xs import evaluate as xe
from filteratlas.xs import ncmat as nb

C = cfg("shortlist")["transport"]
TAGS = {"293K": 293.15, "77K": 77.0}


def dos_for(mid, vdos_source, cand):
    if vdos_source == "sevennet":
        return mq.load_dos(MLIP_DIR / f"{mid}.json.gz")
    r = cand.loc[mid]
    return mq.load_dos(D("mp") / "phonon_dos" / r.phonon_method / f"{r.phonon_id}.json.gz")


def ncmat_text(dos):
    c, labels = nb.composer(dos["structure"])
    vd = nb.element_vdos(dos)
    for lab, el in labels.items():
        eg, dd = vd[el]
        c.set_dyninfo_vdos(lab, eg, dd)
    text = c.create_ncmat()
    if "vdosdebye" in text or "debye_temp" in text:
        raise RuntimeError("fallback Debye dynamics in the generated NCMAT")
    return text


def main():
    out_dir = D("stage4_transport")
    (out_dir / "ncmat").mkdir(parents=True, exist_ok=True)
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    at = pd.read_parquet(F("atlas.parquet")).set_index(["material_id", "T"])
    cand = pd.read_parquet(F("stage0_candidates.parquet")).set_index("material_id")
    cases = []
    for r in sl.itertuples():
        text = ncmat_text(dos_for(r.material_id, r.vdos_source, cand))
        path = out_dir / "ncmat" / f"{r.material_id}.ncmat"
        path.write_text(text, encoding="utf-8")
        name = nb.build_vdos_material(dos_for(r.material_id, r.vdos_source, cand), f"s4_{r.material_id}")
        sig = rj.sigma_tot(r.formula)
        for tag, T in TAGS.items():
            a = at.loc[(r.material_id, T)]
            # the NCMAT reproduces the atlas: T_th at the atlas L0.8 must be 0.80
            mu, _ = xe.mu_offbragg(name, T, xe.LAM)
            tth = float(xe.t_band(mu, a.L08))
            if abs(tth - 0.8) > 1e-3:
                raise RuntimeError(f"{r.material_id} {tag}: T_th(L0.8) = {tth}")
            g1, g6 = a.A_gamma_1MeV / a.L08, a.A_gamma_6MeV / a.L08
            for kind, L in (("L08", float(a.L08)), ("5cm", C["min_thickness_cm"])):
                if kind == "5cm" and a.L08 >= C["min_thickness_cm"]:
                    continue
                cases.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol,
                              "T": T, "tag": tag, "thickness": kind, "L_cm": L, "ncmat": path.name,
                              "analytic": {"A_fast_1E": float(rj.a_fast(sig, a.numberdensity_A3, L, "1E")[0]),
                                           "A_fast_watt": float(rj.a_fast(sig, a.numberdensity_A3, L, "watt")[0]),
                                           "A_gamma_1MeV": float(g1 * L), "A_gamma_6MeV": float(g6 * L),
                                           "T_th": float(xe.t_band(mu, L))}})
    (out_dir / "cases.json").write_text(json.dumps({"lambda_T_A": xe.LAMBDA_T, "band_A": list(map(float, xe.BAND)),
                                                    "cases": cases}, indent=1), encoding="utf-8")
    print(len(sl), "compounds,", len(cases), "slab cases")


if __name__ == "__main__":
    main()
