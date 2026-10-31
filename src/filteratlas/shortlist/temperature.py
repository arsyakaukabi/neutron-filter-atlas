"""Temperature behaviour of the Tier-A shortlist.

T_th(1 cm), L0.8, A_fast (1/E, Watt) and A_gamma (1, 6 MeV) at 293, 250, 200, 150, 100 and 77 K, with the atlas code and
the compound's atlas VDOS (NCrystal temperature only; the free-atom fast and photon data are temperature-independent,
as in the atlas). Cooling gain per interval: dL0.8/dT (fixed before the computation) and d ln L0.8/dT (secondary reading).
-> derived/stage4_temperature.json
"""
import json

import numpy as np
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import F
from filteratlas.fom import rejection as rj
from filteratlas.shortlist.transport_inputs import dos_for
from filteratlas.xs import evaluate as xe
from filteratlas.xs import ncmat as nb

TS = cfg("shortlist")["temperature_series_K"]


def main():
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    sl = sl[sl.tier == "A"]
    at = pd.read_parquet(F("atlas.parquet")).set_index(["material_id", "T"])
    cand = pd.read_parquet(F("stage0_candidates.parquet")).set_index("material_id")
    out = {"temperatures_K": TS, "compounds": {}}
    for r in sl.itertuples():
        name = nb.build_vdos_material(dos_for(r.material_id, r.vdos_source, cand), f"tser_{r.material_id}")
        sig = rj.sigma_tot(r.formula)
        rows = []
        for T in TS:
            mu, info = xe.mu_offbragg(name, T, xe.LAM)
            n = info.numberdensity
            L = float(xe.l_target(mu))
            g1, g6 = rj.mu_gamma(r.formula, n)
            rows.append({"T": T, "T_th_1cm": float(xe.t_band(mu, 1.0)), "L08": L,
                         "A_fast_1E": float(rj.a_fast(sig, n, L, "1E")[0]), "A_fast_watt": float(rj.a_fast(sig, n, L, "watt")[0]),
                         "A_gamma_1MeV": g1 * L, "A_gamma_6MeV": g6 * L})
        for T, tag in ((293.15, "293K"), (77.0, "77K")):   # the series reproduces the atlas at the atlas temperatures
            a = at.loc[(r.material_id, T)]
            row = next(x for x in rows if x["T"] == T)
            if abs(row["L08"] / a.L08 - 1) > 1e-6:
                raise RuntimeError(f"{r.material_id} {tag}: L0.8 {row['L08']} vs atlas {a.L08}")
        gains = []
        for a_, b_ in zip(rows[:-1], rows[1:]):
            gains.append({"from_K": a_["T"], "to_K": b_["T"], "dL08_cm": b_["L08"] - a_["L08"],
                          "dL08_per_K": (b_["L08"] - a_["L08"]) / (a_["T"] - b_["T"]),
                          "dlnL08_per_K": (np.log(b_["L08"]) - np.log(a_["L08"])) / (a_["T"] - b_["T"])})
        # the cooling gain is dL0.8/dT; the relative gain d ln L0.8/dT is a secondary reading
        best = max(gains, key=lambda g: g["dL08_per_K"])
        best_ln = max(gains, key=lambda g: g["dlnL08_per_K"])
        out["compounds"][f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]"] = {
            "series": rows, "intervals": gains,
            "largest_gain_per_K": {"from_K": best["from_K"], "to_K": best["to_K"], "dL08_per_K": best["dL08_per_K"]},
            "largest_relative_gain_per_K": {"from_K": best_ln["from_K"], "to_K": best_ln["to_K"],
                                            "dlnL08_per_K": best_ln["dlnL08_per_K"]},
            "L08_gain_293_to_77": rows[-1]["L08"] / rows[0]["L08"]}
        print(r.formula, [round(x["L08"], 3) for x in rows], flush=True)
    F("stage4_temperature.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
