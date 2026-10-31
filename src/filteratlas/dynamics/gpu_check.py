"""GPU vs CPU reproducibility of SevenNet-MF-ompa phonons (2 compounds of a first GPU test plus a random 20).
The first two GPU DOS (gpu_check_cuda.json.gz) were computed with the MLIP-benchmark protocol for two compounds of the
benchmark random set; the CPU reference is the cached benchmark CPU output of the same compounds
(D("mlip")/sevennet_mf_ompa/). Both go through the atlas code path (NCrystal, 2 meV floor).
-> stage3_laptop_gpu_check.json
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import TEMPS  # noqa: E402


def l08(dos, tag):
    cfg = nb.build_vdos_material(dos, f"gpuchk_{tag}")
    return xe.evaluate(cfg, TEMPS)


def main():
    gpu = json.load(gzip.open(D("laptop_gpu_check") / "gpu_check_cuda.json.gz", "rt"))
    rows = []
    for mid, g in gpu.items():
        ref = mq.load_dos(D("mlip") / "sevennet_mf_ompa" / f"{mid}.json.gz")
        a, b = l08(g["dos"], mid + "_gpu"), l08(ref, mid + "_cpu")
        rows.append({"material_id": mid, "supercell_atoms": g["natoms_sc"], "gpu_wall_s": round(g["wall"], 1),
                     **{f"dlnL08_{t}": round(float(np.log(a[f"L08_{t}"] / b[f"L08_{t}"])), 5) for t in ("293K", "77K")}})
    # 20 random compounds of the MLIP-benchmark set, recomputed on the GPU by filteratlas.dynamics.gpu_check_run
    for f in sorted((D("laptop_gpu_check") / "gpu20").glob("*.json.gz")):
        mid = f.name.replace(".json.gz", "")
        g = mq.load_dos(f)
        ref = mq.load_dos(D("mlip") / "sevennet_mf_ompa" / f"{mid}.json.gz")
        a, b = l08(g, mid + "_g20"), l08(ref, mid + "_c20")
        rows.append({"material_id": mid, "supercell_atoms": g["meta"]["supercell_atoms"], "gpu_wall_s": g["meta"]["wall_s"],
                     "sample": "random20", **{f"dlnL08_{t}": round(float(np.log(a[f"L08_{t}"] / b[f"L08_{t}"])), 5)
                                              for t in ("293K", "77K")}})
    d = np.array([abs(r[k]) for r in rows for k in ("dlnL08_293K",)])
    out = {"rows": rows, "n": len(rows), "median_abs_dlnL08_293K": round(float(np.median(d)), 5),
           "max_abs_dlnL08": max(abs(r[k]) for r in rows for k in ("dlnL08_293K", "dlnL08_77K")),
           "stage2_gpu_vs_cpu_max_abs_dlnL08": json.loads((F("stage2_mlip_sensitivity.json")).read_text())
           ["sevennet_cuda"]["vs_reference_run_abs_dlnL08_293K"]["max"]}
    (F("stage3_laptop_gpu_check.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
