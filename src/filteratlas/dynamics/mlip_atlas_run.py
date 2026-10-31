"""SevenNet phonons for the atlas. The protocol is that of the MLIP benchmark, unchanged
(filteratlas.dynamics.mlip_run: relax with FixSymmetry/FrechetCellFilter/FIRE fmax 0.005 <= 500 steps; phonopy 0.01 A
displacements; diagonal supercell >= 12 A per axis, or the MP supercell matrix when MP force constants exist (drift
set); 40 A mesh; commensurate-q minimum; site-projected DOS in the MP DOS format).

Compounds come from stage3_mlip_set.json (step 03 funnel). CPU and GPU outputs go to the same folder, because the GPU
run reproduces the CPU one (stage2_mlip_sensitivity.json); the device is recorded per compound.
-> D("mlip_atlas")/sevennet_mf_ompa/<material_id>.json.gz, D("logs")/mlip_atlas.jsonl
Usage (MLIP environment; normally launched by filteratlas.dynamics.scheduler):
    python -m filteratlas.dynamics.mlip_atlas_run --ids id1,id2 [--threads N] [--device cpu|cuda]
"""
import argparse
import gzip
import json
import time
import traceback
from pathlib import Path

import numpy as np

from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F
from filteratlas.dynamics import mlip_run as s2

MODEL = cfg("mlip")["atlas_model"]   # selected by the MLIP-benchmark rule fixed before the computation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--relax-sensitivity", action="store_true",
                    help="D20 post-hoc: relax with relax_sensitivity_max_steps; separate folder and log, never read by the atlas")
    a = ap.parse_args()
    steps = cfg("mlip")["relax_sensitivity_max_steps"] if a.relax_sensitivity else None
    import os
    import torch
    torch.set_num_threads(a.threads)
    if a.device == "cuda":   # under the Windows GPU driver a job above GPU memory spills into shared system memory
        # and crawls; capping PyTorch's allocator turns that into an immediate out-of-memory error instead
        torch.cuda.set_per_process_memory_fraction(float(os.environ.get("FA_GPU_MEM_FRACTION",
                                                                        cfg("mlip")["gpu_memory_fraction"])))
    sel = {c["material_id"]: c for c in json.loads(F("stage3_mlip_set.json").read_text())["compounds"]}
    ids = [i for i in a.ids.split(",") if i in sel]
    outdir = D("mlip_relax_sensitivity" if a.relax_sensitivity else "mlip_atlas") / MODEL
    outdir.mkdir(parents=True, exist_ok=True)
    logf = D("logs") / ("mlip_relax_sensitivity.jsonl" if a.relax_sensitivity else "mlip_atlas.jsonl")
    calc, prec = s2.calculator(MODEL, a.device)
    structs = s2.load_structures()
    for k, mid in enumerate(ids):
        c = sel[mid]
        out = outdir / f"{mid}.json.gz"
        if out.exists():
            continue
        rec = {"material_id": mid, "formula": c["formula"], "arm": c["arm"], "model": MODEL, "precision": prec,
               "disp_A": s2.DISP, "device": a.device, "threads": a.threads}
        if steps is not None:
            rec["relax_max_steps"] = steps
        t0 = time.time()
        try:
            from ase import Atoms
            syms, lat, frac, sc, prim, _ = s2.start_cell(c, structs)
            rec["supercell_atoms"] = int(len(syms) * round(abs(np.linalg.det(sc))))
            at, conv, nstep, fmax = s2.relax(Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True), calc, steps)
            rec.update({"relax_converged": conv, "relax_steps": nstep, "relax_fmax": fmax,
                        "volume_ratio": float(at.get_volume() / abs(np.linalg.det(lat)))})
            dos, fmesh, fcomm, ndisp = s2.phonons(list(at.get_chemical_symbols()), np.asarray(at.cell),
                                                  at.get_scaled_positions(), sc, prim, calc=calc)
            rec.update({"status": "ok", "min_freq_mesh_THz": fmesh, "min_freq_comm_THz": fcomm, "n_disp": ndisp})
            with gzip.open(out, "wt", encoding="utf-8") as fh:
                json.dump({**dos, "meta": rec}, fh)
        except Exception as e:  # noqa: BLE001
            rec.update({"status": "error", "error": (repr(e) + " | " + traceback.format_exc().splitlines()[-2])[:300]})
        rec["wall_s"] = round(time.time() - t0, 1)
        with open(logf, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"[{k + 1}/{len(ids)}] {c['formula']} {rec['status']} {rec['wall_s']}s", flush=True)


if __name__ == "__main__":
    main()
