"""MLIP benchmark diagnostic, added after inspection of the results: are ORB v3's commensurate imaginary modes a
property of the model or of our calculation settings?

ORB v3 conservative-inf-omat failed 41 / 149 random compounds (commensurate-q min freq < -0.5 THz, all relaxations
converged; MACE-MP-0 failed 7 / 149 with the same pipeline). For 7 small failed compounds this re-runs the
filteratlas.dynamics.mlip_run pipeline with one setting changed at a time:
  base    - identical settings (reproducibility check)
  disp03  - finite displacement 0.03 A instead of 0.01 A (force-noise / anharmonicity check)
  nn_inf  - max_num_neighbors 10000 instead of the adapter's 120 (kNN truncation of the 6 A graph)
  tight   - relaxation fmax 0.0005 eV/A, <= 2000 steps (residual-force check)
  norelax - MP DFT structure, no relaxation (geometry vs force-constant check)
  mpa     - orb_v3_conservative_inf_mpa weights (MPtrj + Alexandria training set instead of OMat24)
The benchmark ORB result, fixed before the computation, is NOT replaced by anything here; this only decides how the
failure is described.
-> the ORB diagnostic JSON in data/derived
Usage (MLIP environment): python pipelines/05_validation.py --with mlip
  or: PYTHONPATH=src python -m filteratlas.validation.orb_diagnostic [--threads N]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from filteratlas.common.paths import D, F

from filteratlas.dynamics import mlip_run as R

COMPOUNDS = ("P2Ru", "As2Ru", "Ca2Sb", "ZnPdF6", "CaSi2", "GeTe", "Cs")
VARIANTS = ("base", "disp03", "nn_inf", "tight", "norelax", "mpa")


def orb_calc(weights="omat", max_nn=None):
    from orb_models.forcefield import pretrained
    from orb_models.forcefield.inference.calculator import ORBCalculator
    load = pretrained.orb_v3_conservative_inf_omat if weights == "omat" else pretrained.orb_v3_conservative_inf_mpa
    orbff, adapter = load(device="cpu", precision="float64", compile=False)
    return ORBCalculator(orbff, atoms_adapter=adapter, device="cpu", max_num_neighbors=max_nn)


def run_one(c, structs, calc, variant):
    from ase import Atoms
    R.DISP, R.FMAX, R.MAXSTEP = 0.01, 0.005, 500
    if variant == "disp03":
        R.DISP = 0.03
    if variant == "tight":
        R.FMAX, R.MAXSTEP = 0.0005, 2000
    syms, lat, frac, sc, prim, _ = R.start_cell(c, structs)
    rec = {"formula": c["formula"], "material_id": c["material_id"], "variant": variant, "disp_A": R.DISP}
    t0 = time.time()
    if variant != "norelax":
        at, conv, nstep, fmax = R.relax(Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True), calc)
        rec.update({"relax_converged": conv, "relax_steps": nstep, "relax_fmax": fmax,
                    "volume_ratio": float(at.get_volume() / abs(np.linalg.det(lat)))})
        syms, lat, frac = list(at.get_chemical_symbols()), np.asarray(at.cell), at.get_scaled_positions()
    else:
        a = Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True)
        a.calc = calc
        rec.update({"volume_ratio": 1.0, "residual_fmax_at_mp_structure": float(np.abs(a.get_forces()).max())})
    _, fmesh, fcomm, nd = R.phonons(syms, lat, frac, sc, prim, calc=calc)
    rec.update({"min_freq_mesh_THz": fmesh, "min_freq_comm_THz": fcomm, "n_disp": nd,
                "failed": bool(fcomm < -0.5), "wall_s": round(time.time() - t0, 1)})
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    cset = {c["formula"]: c for c in json.loads(F("stage2_mlip_set.json").read_text())["compounds"]
            if c["arm"].startswith("random")}
    structs = R.load_structures()
    calcs = {"omat": orb_calc("omat"), "omat_nn": orb_calc("omat", 10000), "mpa": orb_calc("mpa")}
    out = []
    for f in COMPOUNDS:
        for v in VARIANTS:
            calc = calcs["omat_nn" if v == "nn_inf" else "mpa" if v == "mpa" else "omat"]
            try:
                rec = run_one(cset[f], structs, calc, v)
            except Exception as e:  # noqa: BLE001
                rec = {"formula": f, "variant": v, "error": repr(e)[:200]}
            out.append(rec)
            print(json.dumps({k: (round(x, 3) if isinstance(x, float) else x) for k, x in rec.items()}), flush=True)
    summary = {v: {"n_failed": sum(r.get("failed", False) for r in out if r["variant"] == v),
                   "n": sum(1 for r in out if r["variant"] == v and "failed" in r)} for v in VARIANTS}
    res = {"compounds": COMPOUNDS, "variants": VARIANTS, "fail_rule": "commensurate-q min freq < -0.5 THz",
           "summary": summary, "runs": out}
    F("stage2_orb_diagnostic.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
