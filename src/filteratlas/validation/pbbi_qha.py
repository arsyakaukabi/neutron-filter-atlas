"""Pb / Bi diagnostic of the EXFOR benchmark (DIAGNOSTIC ONLY; the atlas keeps the static model fixed before the
computation).

Quasi-harmonic test: does the thermal expansion at the measurement temperature close the Pb / Bi gap to the
single-crystal data? With SevenNet-MF-ompa (the atlas MLIP, MLIP-benchmark settings otherwise):
  1. relax the MP structure (benchmark settings), V0 = static volume;
  2. 7 isotropic strains V/V0 = 0.97 ... 1.06; at each, relax internal positions only (FixSymmetry, cell fixed;
     for Bi the static c/a is kept: an approximation, Bi's thermal expansion is anisotropic), electronic energy and
     phonopy (0.01 A, diagonal supercell >= 12 A as in filteratlas.dynamics.mlip_run), thermal properties 0-400 K;
  3. phonopy QHA (Vinet EOS) -> V(T) at the dataset temperatures;
  4. phonons at V(T) -> site-projected DOS in the MP format, <mlip_sensitivity dir>/sevennet_mf_ompa_qha<T>K/<id>.json.gz.
-> the Pb/Bi QHA JSON in data/derived (V(T)/V0, energies, free energies), the DOS files above.
Usage (MLIP environment with SevenNet): python pipelines/05_validation.py --with mlip
  or: PYTHONPATH=src python -m filteratlas.validation.pbbi_qha [--threads N]
"""
import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np

from filteratlas.common.paths import D, F

from filteratlas.dynamics import mlip_run as R  # noqa: E402

MODEL = "sevennet_mf_ompa"
STRAINS = (0.97, 0.985, 1.0, 1.015, 1.03, 1.045, 1.06)
TEMPS = (77.0, 293.0, 297.0)


def relax_positions(atoms, calc):
    from ase.constraints import FixSymmetry
    from ase.optimize import FIRE
    atoms = atoms.copy()
    atoms.calc = calc
    atoms.set_constraint(FixSymmetry(atoms))
    FIRE(atoms, logfile=None).run(fmax=R.FMAX, steps=R.MAXSTEP)
    atoms.set_constraint()
    return atoms, float(atoms.get_potential_energy())


def thermal(atoms, sc, calc):
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    ph = Phonopy(PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.asarray(atoms.cell),
                              scaled_positions=atoms.get_scaled_positions()), supercell_matrix=sc)
    ph.generate_displacements(distance=R.DISP)
    f = []
    for s in ph.supercells_with_displacements:
        a = Atoms(symbols=s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        a.calc = calc
        f.append(a.get_forces())
    ph.forces = np.array(f)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_mesh(R.MESH, is_gamma_center=True)
    fmin = float(np.min(ph.mesh.frequencies))
    ph.run_thermal_properties(t_min=0, t_max=400, t_step=1)
    tp = ph.get_thermal_properties_dict()
    return tp, fmin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.threads)
    from ase import Atoms
    from phonopy import PhonopyQHA
    calc, _ = R.calculator(MODEL)
    structs = R.load_structures()
    cset = {c["formula"]: c for c in json.loads(F("stage2_mlip_set.json").read_text())["compounds"]
            if c["formula"] in ("Pb", "Bi")}
    out = {"model": MODEL, "strains_V_over_V0": STRAINS, "note": "diagnostic only; isotropic strain (Bi c/a fixed)"}
    for f, c in cset.items():
        syms, lat, frac, sc, prim, _ = R.start_cell(c, structs)
        at0, conv, nstep, fmax = R.relax(Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True), calc)
        V0 = float(at0.get_volume())
        vols, ens, fe, cv, en, fmins = [], [], [], [], [], []
        for s in STRAINS:
            at = at0.copy()
            at.set_cell(np.asarray(at0.cell) * s ** (1 / 3), scale_atoms=True)
            at, E = relax_positions(at, calc)
            tp, fmin = thermal(at, sc, calc)
            vols.append(float(at.get_volume())); ens.append(E); fmins.append(fmin)
            fe.append(tp["free_energy"]); cv.append(tp["heat_capacity"]); en.append(tp["entropy"])
            print(f, s, round(E, 5), "min freq", round(fmin, 3), flush=True)
        T = np.asarray(tp["temperatures"])
        q = PhonopyQHA(volumes=vols, electronic_energies=ens, temperatures=T, free_energy=np.array(fe).T,
                       cv=np.array(cv).T, entropy=np.array(en).T, eos="vinet", t_max=400)
        VT = np.asarray(q.volume_temperature)
        rec = {"material_id": c["material_id"], "relax_converged": conv, "V0_static_A3": V0, "volumes_A3": vols,
               "energies_eV": ens, "min_freq_mesh_THz": fmins,
               "V_over_V0_at_T": {}, "dos_files": {}}
        for t in TEMPS:
            vt = float(np.interp(t, T[:len(VT)], VT))
            rec["V_over_V0_at_T"][str(int(t))] = round(vt / V0, 5)
            at = at0.copy()
            at.set_cell(np.asarray(at0.cell) * (vt / V0) ** (1 / 3), scale_atoms=True)
            at, _ = relax_positions(at, calc)
            dos, fm, fcm, nd = R.phonons(list(at.get_chemical_symbols()), np.asarray(at.cell), at.get_scaled_positions(),
                                         sc, prim, calc=calc)
            od = D("mlip_sensitivity") / f"{MODEL}_qha{int(t)}K"
            od.mkdir(parents=True, exist_ok=True)
            p = od / f"{c['material_id']}.json.gz"
            with gzip.open(p, "wt", encoding="utf-8") as fh:
                json.dump({**dos, "meta": {"material_id": c["material_id"], "formula": f, "model": MODEL, "T_K": t,
                                           "V_over_V0": vt / V0, "min_freq_mesh_THz": fm,
                                           "min_freq_comm_THz": fcm}}, fh)
            rec["dos_files"][str(int(t))] = f"cache/mlip_sensitivity/{od.name}/{p.name}"   # older-layout label, not a live path
            print(f, "T", t, "V/V0", round(vt / V0, 4), flush=True)
        out[f] = rec
    F("stage2_pbbi_qha.json").write_text(json.dumps(out, indent=1))
    print("QHA_DONE", flush=True)


if __name__ == "__main__":
    main()
