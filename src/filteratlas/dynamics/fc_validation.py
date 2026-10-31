"""Principled low-energy VDOS treatment, validated against force constants.

For a validation set (all enantiomorph entries + classic filters + 30 random usable compounds with <= 20 sites,
fixed seed), rebuild the phonon DOS from the MP force constants with phonopy on a dense q-mesh and the
tetrahedron method. That DOS has the correct g ~ E^2 acoustic onset, so it is the reference. Compare L0.8
(293 K, 77 K) from: (i) the FC-tetrahedron VDOS, (ii) the MP DOS conversion with floors 0/1/2/3 meV.
Output: stage1r_fc_validation.parquet + .json
"""
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.config import cfg, tup
from filteratlas.common.constants import THZ_TO_EV, TEMPS, VDOS_NPTS  # noqa: E402
from filteratlas.common.paths import D, F

FLOORS = tup(cfg("stability")["fc_validation_floors_meV"])
FCDIR = D("mp") / "force_constants"
FCDIR.mkdir(parents=True, exist_ok=True)
MESH_LENGTH = cfg("stability")["fc_validation_mesh_length_A"]   # Angstrom; phonopy q-mesh density


def fc_vdos(rester, row):
    """Per-element VDOS from MP force constants via phonopy (dense mesh, tetrahedron). Cached."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    out = FCDIR / f"{row.phonon_id}.npz"
    if not out.exists():
        legacy = str(rester.materials.summary.search(material_ids=[row.material_id], fields=["material_id"])[0]
                     .material_id)
        docs = rester.materials.phonon.search(material_ids=[legacy], fields=[
            "supercell_matrix", "primitive_matrix", "structure", "phonon_method", "identifier"])
        doc = next(d for d in docs if str(d.phonon_method).endswith(row.phonon_method)
                   and str(d.identifier) == row.phonon_id)
        fc = np.array(rester.materials.phonon.get_forceconstants_from_phonon_id(row.phonon_id, row.phonon_method))
        s = doc.structure
        np.savez_compressed(out, fc=fc, sc=np.array(doc.supercell_matrix), prim=np.eye(3) if doc.primitive_matrix is None else np.array(doc.primitive_matrix),
                            lattice=s.lattice.matrix, frac=s.frac_coords, symbols=np.array([str(x.specie) for x in s]))
    z = np.load(out, allow_pickle=False)
    unit = PhonopyAtoms(symbols=list(z["symbols"]), cell=z["lattice"], scaled_positions=z["frac"])
    ph = Phonopy(unit, supercell_matrix=z["sc"], primitive_matrix=z["prim"])
    ph.force_constants = z["fc"]
    ph.run_mesh(MESH_LENGTH, with_eigenvectors=True, is_mesh_symmetry=False)
    ph.run_projected_dos(freq_pitch=0.02)
    pd_ = ph.projected_dos
    f, proj = np.asarray(pd_.frequency_points), np.asarray(pd_.projected_dos)
    syms = ph.primitive.symbols
    vd = {}
    for el in sorted(set(syms)):
        dens = np.clip(proj[[i for i, x in enumerate(syms) if x == el]].sum(axis=0), 0, None)
        m = f * THZ_TO_EV >= 1e-5
        grid = np.linspace(f[m][0], f[m][-1], VDOS_NPTS)
        vd[el] = (grid * THZ_TO_EV, np.interp(grid, f[m], dens[m]))
    fmin = float(np.min(ph.mesh.frequencies))
    return vd, fmin


def main():
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    t = t[(t.status == "ok") & t.stable.fillna(False) & (t.phonon_method == "pheasy")]
    en = json.loads((F("stage1_enantiomorphs.json")).read_text())
    ids = {p["id_a"] for p in en["pairs"]} | {p["id_b"] for p in en["pairs"]}
    classic = {("Al2O3", "R-3c"), ("Si", "Fd-3m"), ("MgO", "Fm-3m"), ("MgF2", "P4_2/mnm"), ("BeF2", "P3_121"),
               ("BeO", "P6_3mc"), ("PbF2", "Fm-3m"), ("AlP", "F-43m"), ("AlPO4", "P6_422")}
    ids |= set(t[[(f, s) in classic for f, s in zip(t.formula, t.sg_symbol)]].material_id)
    pool = t[(t["mu_1.8_293K"] < 0.3) & (t.nsites <= 20) & ~t.material_id.isin(ids)]
    ids |= set(pool.sort_values("material_id").sample(30, random_state=20261021).material_id)
    sel = t[t.material_id.isin(ids)]
    rows = []
    with mq.rester() as r:
        for row in sel.itertuples():
            rec = {"material_id": row.material_id, "formula": row.formula, "sg": row.sg_symbol, "nsites": row.nsites}
            t0 = time.time()
            try:
                dos = mq.load_dos(D("mp") / "phonon_dos" / row.phonon_method / f"{row.phonon_id}.json.gz")
                vd, fmin = fc_vdos(r, row)
                rec["fc_min_freq_THz"] = fmin
                ev = xe.evaluate(nb.build_vdos_material(dos, f"fc_{row.material_id}", vdos=vd), TEMPS)
                rec.update({f"fc_{k}": v for k, v in ev.items() if k.startswith("L08")})
                for fl in FLOORS:
                    ev = xe.evaluate(nb.build_vdos_material(dos, f"fl{fl}_{row.material_id}", floor_meV=fl), TEMPS)
                    rec.update({f"floor{fl:g}_{k}": v for k, v in ev.items() if k.startswith("L08")})
                rec["status"] = "ok"
            except Exception as e:  # noqa: BLE001
                rec["status"] = "error"
                rec["error"] = (repr(e) + " | " + traceback.format_exc().splitlines()[-2])[:300]
            rec["seconds"] = round(time.time() - t0, 1)
            print(rec["formula"], rec["sg"], rec["status"], rec["seconds"], flush=True)
            rows.append(rec)
    pd.DataFrame(rows).to_parquet(F("stage1r_fc_validation.parquet"))
    summarise()


def summarise():
    """FC-tetrahedron DOS is a valid reference only where the FC reconstruction is dynamically clean
    (min frequency on the dense mesh >= -0.05 THz); entries with imaginary branches are reported separately."""
    df = pd.read_parquet(F("stage1r_fc_validation.parquet"))
    ok = df[df.status == "ok"]
    clean = ok[ok.fc_min_freq_THz >= -0.05]
    res = {"n_selected": len(df), "n_ok": len(ok), "n_fc_clean": len(clean), "mesh_length_A": MESH_LENGTH,
           "n_errors": int((df.status != "ok").sum())}
    for subset, d0 in (("fc_clean", clean), ("all_ok", ok)):
        res[subset] = {}
        for fl in FLOORS:
            for tag in ("293K", "77K"):
                d = np.log(d0[f"floor{fl:g}_L08_{tag}"] / d0[f"fc_L08_{tag}"])
                res[subset][f"floor{fl:g}_{tag}"] = {"median_dlog": round(float(d.median()), 4),
                                                     "median_abs_dlog": round(float(d.abs().median()), 4),
                                                     "p90_abs_dlog": round(float(d.abs().quantile(0.9)), 4),
                                                     "max_abs_dlog": round(float(d.abs().max()), 4),
                                                     "rms_dlog": round(float(np.sqrt((d ** 2).mean())), 4)}
    (F("stage1r_fc_validation.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    summarise() if "--summary" in sys.argv else main()
