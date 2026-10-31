"""Model control (fixed before the computation): the Debye-300 K fallback on sapphire must reproduce the reference
~2.85x overestimate (configs/controls.yaml) vs EXFOR 31743 (median model/data, 1 meV-1 eV, off-Bragg baseline, 293 K).
Also reports the MP-pheasy VDOS and NCrystal curated models for context. -> stage0_model_control.json
"""
import json
import sys
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb  # noqa: E402
from filteratlas.common.config import cfg as _cfg
from filteratlas.common.paths import D, F  # noqa: E402

MC = _cfg("controls")["model_control"]


def main():
    ex = pd.read_csv(F("exfor_benchmark.csv"))
    (e_lo, e_hi) = MC["energy_window_eV"]
    d = ex[(ex.entry == MC["exfor_entry"]) & (ex.E_eV >= e_lo) & (ex.E_eV <= e_hi)]
    dos = mq.load_dos(D("mp") / "phonon_dos" / MC["phonon_method"] / f"{MC['phonon_id']}.json.gz")   # mp-aaaaabrz, R-3c sapphire
    models = {
        "debye300_fallback": nb.build_debye_material(dos["structure"], MC["debye_K"], "ctrl_sapphire_d300"),
        "mp_pheasy_vdos": nb.build_vdos_material(dos, "ctrl_sapphire_mpvdos"),
        "ncrystal_curated": "stdlib::Al2O3_sg167_Corundum.ncmat",
    }
    res = {"exfor": "31743 (Cantargi et al. 2015, 293 K, single crystal)", "n_points": int(len(d))}
    for k, cfg in models.items():
        c = f"{cfg};temp=293.15K"
        s = (NC.createScatter(f"{c};comp=inelas,incoh_elas").xsect(ekin=d.E_eV.values)
             + NC.createAbsorption(c).xsect(ekin=d.E_eV.values))
        r = s / d.sigma_per_atom_b.values
        res[k] = {"median_ratio": round(float(np.median(r)), 3),
                  "rms_log_ratio": round(float(np.sqrt(np.mean(np.log(r) ** 2))), 3)}
    ref, tol = MC["phase0_reference_ratio"], MC["tolerance"]
    res["phase0_reference_median_ratio_debye300"] = ref
    res["control_pass"] = bool(abs(res["debye300_fallback"]["median_ratio"] - ref) / ref < tol)
    (F("stage0_model_control.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
