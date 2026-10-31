"""GPU reproducibility sample: 20 random compounds of the MLIP-benchmark random set (seed in configs/mlip.yaml), whose
CPU SevenNet outputs are cached (D("mlip")/sevennet_mf_ompa/), recomputed on the GPU with the identical protocol
(filteratlas.dynamics.mlip_run; 0.9 allocator cap). One compound at a time; resumable.
-> gpu20/<material_id>.json.gz in the GPU-check data directory (compared by filteratlas.dynamics.gpu_check)
Usage (MLIP environment with a CUDA build of torch): python -m filteratlas.dynamics.gpu_check_run
"""
import gzip
import json
import time
from pathlib import Path

import numpy as np
import torch

from filteratlas.dynamics import mlip_run as s2

from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F

OUT = D("laptop_gpu_check") / "gpu20"
_G = cfg("mlip")["gpu_check"]


def pick():
    S = json.loads(F("stage2_mlip_set.json").read_text())["compounds"]
    c = [x for x in S if x["arm"].startswith("random")
         and (D("mlip") / "sevennet_mf_ompa" / f"{x['material_id']}.json.gz").exists()]
    rng = np.random.default_rng(_G["seed"])
    return [c[i] for i in sorted(rng.choice(len(c), _G["n"], replace=False))]


def main():
    torch.cuda.set_per_process_memory_fraction(float(cfg("mlip")["gpu_memory_fraction"]))
    OUT.mkdir(parents=True, exist_ok=True)
    calc, _ = s2.calculator("sevennet_mf_ompa", "cuda")
    structs = s2.load_structures()
    from ase import Atoms
    for c in pick():
        f = OUT / f"{c['material_id']}.json.gz"
        if f.exists():
            continue
        t0 = time.time()
        try:
            syms, lat, frac, sc, prim, _ = s2.start_cell(c, structs)
            at, conv, n, fm = s2.relax(Atoms(symbols=syms, cell=lat, scaled_positions=frac, pbc=True), calc)
            dos, fmesh, fcomm, nd = s2.phonons(list(at.get_chemical_symbols()), np.asarray(at.cell),
                                               at.get_scaled_positions(), sc, prim, calc=calc)
            meta = {"relax_converged": conv, "wall_s": round(time.time() - t0, 1),
                    "supercell_atoms": int(len(syms) * round(abs(np.linalg.det(sc))))}
            with gzip.open(f, "wt", encoding="utf-8") as fh:
                json.dump({**dos, "meta": meta}, fh)
            print(c["formula"], meta, flush=True)
        except Exception as e:  # noqa: BLE001
            print(c["formula"], "error", repr(e)[:200], flush=True)
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
