"""Atlas funnel, first steps (step 01): candidate-set checks and the elemental proxy bound.

For every candidate: mu_proxy(lambda) = n * sum_i x_i (sigma_abs,i(lambda) + sigma_inc,i) with NCrystal 4.4.6 atom data,
n = nsites / volume of the MP summary structure (uniform for all compounds), T_th,proxy(L = 1 cm) and L0.8,proxy over
the thermal band. Compounds with T_th,proxy(1 cm) < 0.5 are labelled excluded_by_proxy_bound.
-> stage3_proxy.parquet, stage3_proxy.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import NCrystal as NC
import pandas as pd
from pymatgen.core import Composition

from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import RADIOACTIVE, STRONG_ABSORBERS  # noqa: E402

from filteratlas.common.config import cfg

T_BOUND, L_BOUND = cfg("funnel")["proxy_bound"]["T_th_min"], cfg("funnel")["proxy_bound"]["L_cm"]     # proxy-bound rule (configs/funnel.yaml)


def atom_xs(sym, cache={}):
    """(sigma_abs at 1.798 A, sigma_inc) [barn] for the natural element."""
    if sym not in cache:
        a = NC.atomDB(sym)
        cache[sym] = (a.captureXS(), a.incoherentXS())
    return cache[sym]


def main():
    c = pd.read_parquet(F("stage0_candidates.parquet"))
    els = c.elements.str.split(",")
    checks = {"n_candidates": len(c),
              "n_with_H": int(els.apply(lambda e: "H" in e).sum()),
              "n_with_radioactive": int(els.apply(lambda e: bool(RADIOACTIVE & set(e))).sum()),
              "n_with_excluded_absorber": int(els.apply(lambda e: bool(STRONG_ABSORBERS & set(e))).sum()),
              "n_theoretical": int(c.theoretical.sum()), "max_e_hull": float(c.e_hull.max())}
    rows = []
    for r in c.itertuples():
        comp = Composition(r.formula).fractional_composition
        n = r.nsites / r.volume                                     # atoms / A^3 -> mu in cm^-1 with barn
        sabs = sum(float(x) * atom_xs(str(e))[0] for e, x in comp.items())
        sinc = sum(float(x) * atom_xs(str(e))[1] for e, x in comp.items())
        mu = n * (sabs * xe.LAM / 1.798 + sinc)
        rows.append({"material_id": r.material_id, "n_summary_A3": n, "sigma_abs_1798_b": sabs, "sigma_inc_b": sinc,
                     "proxy_mu_1.8": float(n * (sabs + sinc)), "proxy_T_1cm": xe.t_band(mu, L_BOUND),
                     "proxy_L08": xe.l_target(mu)})
    p = c[["material_id", "formula", "sg_symbol", "nsites", "e_hull", "flags", "phonon_method", "phonon_id",
           "K_vrh", "G_vrh"]].merge(pd.DataFrame(rows), on="material_id")
    p["excluded_by_proxy_bound"] = p.proxy_T_1cm < T_BOUND
    p = p.sort_values("material_id").reset_index(drop=True)
    p.to_parquet(F("stage3_proxy.parquet"))
    out = {"checks": checks, "rule": f"T_th,proxy(L={L_BOUND} cm) < {T_BOUND}",
           "n_excluded_by_proxy_bound": int(p.excluded_by_proxy_bound.sum()),
           "n_survivors": int((~p.excluded_by_proxy_bound).sum()),
           "survivors_with_mp_phonons": int((~p.excluded_by_proxy_bound & (p.phonon_id != "")).sum()),
           "survivors_without_mp_phonons": int((~p.excluded_by_proxy_bound & (p.phonon_id == "")).sum())}
    (F("stage3_proxy.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
