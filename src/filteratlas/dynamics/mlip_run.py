"""MLIP phonons for the benchmark set (protocol fixed before the computation). No NCrystal needed.
Environments: MACE-MP-0 and ORB v3 in requirements-mlip-benchmark.txt (e3nn 0.4.4, required by the pickled MACE-MP-0
model); SevenNet in requirements-mlip.txt (sevenn pulls e3nn 0.6.0).

For each compound in stage2_mlip_set.json and one model:
  1. start from the MP phonon-calculation cell (force-constant npz) or, if none, the MP summary structure;
  2. relax positions + cell with symmetry fixed (FixSymmetry, FrechetCellFilter, FIRE, fmax 0.005 eV/A, <= 500 steps);
  3. phonopy finite displacements (0.01 A), same supercell matrix as MP pheasy (else diagonal, >= 12 A per axis),
     symmetrised force constants;
  4. minimum frequency on a 40 A Gamma-centred mesh and at supercell-commensurate q (dynamical-stability rule);
  5. site-projected DOS (40 A mesh, tetrahedron) written in the SAME format as the cached MP DOS
     (frequencies THz, densities, projected_densities per site, structure) so that the main pipeline converts it
     with the identical code path (2 meV floor, NCrystal).
Model "dft_fc" skips 1-3 and uses the MP force constants (same phonopy DOS settings) as a DOS-method control.
Resumable: existing outputs are skipped. -> D("mlip")/<model>/<material_id>.json.gz, D("logs")/mlip_<model>.jsonl
Usage (MLIP environment): python -m filteratlas.dynamics.mlip_run <model> [--threads N] [--only arm-prefix]
    [--shard i/n] [--ids id1,id2]
--shard i/n runs only compounds with index % n == i, so n processes can share one model run.
--ids runs only the listed material_ids (one process per compound).
--disp D (sensitivity only; the protocol value is 0.01 A) writes to D("mlip_sensitivity")/<model>_disp<1000*D>/
and D("logs")/mlip_<model>_disp<1000*D>.jsonl, so it never enters the benchmark evaluation (D("mlip")).
--device cuda (GPU test) likewise writes to D("mlip_sensitivity")/<model>_cuda/.
"""
import argparse
import gzip
import json
import time
import traceback
from pathlib import Path

import numpy as np

from filteratlas.common.config import cfg, tup
from filteratlas.common.paths import D, F

_M = cfg("mlip")
MODELS = tup(_M["models"])
DISP, FMAX, MAXSTEP, MESH, SC_MIN_A = (_M["displacement_A"], _M["relax_fmax_eV_per_A"], _M["relax_max_steps"],
                                       _M["mesh_length_A"], cfg("stability")["supercell_min_length_A"])
PROTOCOL_DISP = DISP   # the protocol displacement; --disp overrides DISP for the sensitivity runs


def fc_path(rel):
    """stage2_mlip_set.json stores relative paths of an older layout; resolve the force-constant npz by file name."""
    return D("mp") / "force_constants" / Path(rel).name


def calculator(name, device="cpu"):
    if name == "mace_mp0_medium":
        from mace.calculators import mace_mp
        return mace_mp(model="medium", default_dtype="float64", device=device), "float64"
    if name == "orb_v3_cons_inf_omat":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator
        orbff, adapter = pretrained.orb_v3_conservative_inf_omat(device=device, precision="float64", compile=False)
        return ORBCalculator(orbff, atoms_adapter=adapter, device=device), "float64"
    if name == "sevennet_mf_ompa":
        from sevenn.calculator import SevenNetCalculator
        return SevenNetCalculator("7net-mf-ompa", modal="mpa", device=device), "float32 (package default)"
    raise ValueError(name)


def load_structures():
    out = {}
    with gzip.open(D("mp") / "structures.jsonl.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            out[r["material_id"]] = r["structure"]
    return out


def start_cell(c, structs):
    """(symbols, lattice, frac, supercell_matrix, primitive_matrix, fc or None)."""
    if c["fc_npz"]:
        z = np.load(fc_path(c["fc_npz"]), allow_pickle=False)
        return list(z["symbols"]), z["lattice"], z["frac"], z["sc"], z["prim"], z["fc"]
    from pymatgen.core import Structure
    s = Structure.from_dict(structs[c["material_id"]])
    n = [max(1, int(np.ceil(SC_MIN_A / L))) for L in s.lattice.abc]
    return [str(x.specie) for x in s], s.lattice.matrix, s.frac_coords, np.diag(n), np.eye(3), None


def relax(atoms, calc, steps=None):
    """steps: None = the protocol limit (relax_max_steps); larger only for the relaxation-sensitivity rerun of failed
    relaxations (kept separate from the atlas)."""
    from ase.constraints import FixSymmetry
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    atoms = atoms.copy()
    atoms.calc = calc
    atoms.set_constraint(FixSymmetry(atoms))
    opt = FIRE(FrechetCellFilter(atoms), logfile=None)
    conv = bool(opt.run(fmax=FMAX, steps=MAXSTEP if steps is None else steps))
    fmax = float(np.abs(atoms.get_forces()).max())
    atoms.set_constraint()
    return atoms, conv, int(opt.get_number_of_steps()), fmax


def phonons(syms, lat, frac, sc, prim, calc=None, fc=None):
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.harmonic.dynmat_to_fc import get_commensurate_points
    from phonopy.structure.atoms import PhonopyAtoms
    ph = Phonopy(PhonopyAtoms(symbols=syms, cell=lat, scaled_positions=frac), supercell_matrix=sc, primitive_matrix=prim)
    ndisp = 0
    if fc is None:
        ph.generate_displacements(distance=DISP)
        forces = []
        for s in ph.supercells_with_displacements:
            a = Atoms(symbols=s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
            a.calc = calc
            forces.append(a.get_forces())
        ndisp = len(forces)
        ph.forces = np.array(forces)
        ph.produce_force_constants()
    else:
        ph.force_constants = fc
    ph.symmetrize_force_constants()
    ph.run_mesh(MESH, is_gamma_center=True)
    fmesh = float(np.min(ph.mesh.frequencies))
    smat = np.rint(np.linalg.inv(ph.primitive_matrix) @ ph.supercell_matrix).astype(int)
    qc = get_commensurate_points(smat)
    ph.run_qpoints(qc)
    fq = np.asarray(ph.get_qpoints_dict()["frequencies"])
    gam = np.all(np.abs(qc) < 1e-8, axis=1)
    fcomm = float(min((fq[i][np.argsort(np.abs(fq[i]))[3:]] if gam[i] else fq[i]).min() for i in range(len(qc))
                      if len(fq[i]) > 3 or not gam[i]))
    ph.run_mesh(MESH, with_eigenvectors=True, is_mesh_symmetry=False, is_gamma_center=True)
    ph.run_projected_dos(freq_pitch=0.05)
    pdos = ph.projected_dos
    p = ph.primitive
    struct = {"@module": "pymatgen.core.structure", "@class": "Structure", "charge": 0,
              "lattice": {"matrix": np.asarray(p.cell).tolist()},
              "sites": [{"species": [{"element": e, "occu": 1}], "abc": list(map(float, x)), "label": e,
                         "properties": {}} for e, x in zip(p.symbols, p.scaled_positions)]}
    proj = np.asarray(pdos.projected_dos)
    return {"frequencies": np.asarray(pdos.frequency_points).tolist(), "densities": proj.sum(axis=0).tolist(),
            "projected_densities": proj.tolist(), "structure": struct}, fmesh, fcomm, ndisp


def main():
    global DISP
    ap = argparse.ArgumentParser()
    ap.add_argument("model", choices=MODELS)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--only", default="")
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--ids", default="")
    ap.add_argument("--disp", type=float, default=DISP)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    sens = abs(a.disp - DISP) > 1e-12 or a.device != "cpu"
    DISP = a.disp
    tag = (a.model + (f"_disp{round(DISP * 1000):03d}" if abs(DISP - PROTOCOL_DISP) > 1e-12 else "")
           + ("" if a.device == "cpu" else f"_{a.device}"))
    import torch
    torch.set_num_threads(a.threads)
    cset = [c for c in json.loads(F("stage2_mlip_set.json").read_text())["compounds"]
            if c["arm"].startswith(a.only)]
    si, sn = map(int, a.shard.split("/"))
    cset = [c for k, c in enumerate(cset) if k % sn == si]
    if a.ids:
        cset = [c for c in cset if c["material_id"] in a.ids.split(",")]
    outdir = D("mlip_sensitivity" if sens else "mlip") / tag
    outdir.mkdir(parents=True, exist_ok=True)
    logf = D("logs") / f"mlip_{tag}.jsonl"
    calc, prec = (None, "n/a") if a.model == "dft_fc" else calculator(a.model, a.device)
    structs = load_structures()
    for i, c in enumerate(cset):
        out = outdir / f"{c['material_id']}.json.gz"
        if out.exists():
            continue
        rec = {"material_id": c["material_id"], "formula": c["formula"], "arm": c["arm"], "model": a.model,
               "precision": prec, "disp_A": DISP, "device": a.device}
        t0 = time.time()
        try:
            syms, lat, frac, sc, prim, fc = start_cell(c, structs)
            if a.model == "dft_fc":
                if fc is None:
                    raise ValueError("no MP force constants")
                rec.update({"relax_converged": None})
            else:
                from ase import Atoms
                at, conv, nstep, fmax = relax(Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True), calc)
                rec.update({"relax_converged": conv, "relax_steps": nstep, "relax_fmax": fmax,
                            "volume_ratio": float(at.get_volume() / abs(np.linalg.det(lat)))})
                syms, lat, frac, fc = list(at.get_chemical_symbols()), np.asarray(at.cell), at.get_scaled_positions(), None
            dos, fmesh, fcomm, ndisp = phonons(syms, lat, frac, sc, prim, calc=calc, fc=fc)
            rec.update({"status": "ok", "min_freq_mesh_THz": fmesh, "min_freq_comm_THz": fcomm, "n_disp": ndisp})
            with gzip.open(out, "wt", encoding="utf-8") as fh:
                json.dump({**dos, "meta": rec}, fh)
        except Exception as e:  # noqa: BLE001
            rec.update({"status": "error", "error": (repr(e) + " | " + traceback.format_exc().splitlines()[-2])[:300]})
        rec["wall_s"] = round(time.time() - t0, 1)
        with open(logf, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"[{i + 1}/{len(cset)}] {c['formula']} {rec['status']} {rec['wall_s']}s", flush=True)


if __name__ == "__main__":
    main()
