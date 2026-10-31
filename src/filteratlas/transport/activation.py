"""Activation of the Tier-A shortlist: OpenMC 0.16.0 depletion with IndependentOperator.

1 cm3 of crystal at its atlas number density, irradiated 7 days in a constant thermal flux of 1e8 n cm^-2 s^-1 (a 293 K
Maxwellian, one-group cross sections from MicroXS.from_multigroup_flux), then cooled 1 and 7 days. Chain: the
ENDF/B-VIII.0 thermal-spectrum chain distributed by OpenMC (openmc.org/data; SHA-256 in the output).
Two steps, both run by step 10 (shortlist):
  prepare (main environment):   python pipelines/10_shortlist.py
  run (OpenMC environment, Linux; OPENMC_CROSS_SECTIONS=data/upstream/endfb80_hdf5/cross_sections.xml):
      python pipelines/10_shortlist.py --with openmc
      or: PYTHONPATH=src python -m filteratlas.transport.activation run <activation dir> <work dir>
-> <activation dir>/inputs.json, results.json (read into the practicality dossier)
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from filteratlas.common.config import cfg

C = cfg("shortlist")["dossier"]
K_B_EV = 8.617333262e-5


def prepare():
    import pandas as pd
    from pymatgen.core import Composition

    from filteratlas.common.paths import D, F
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    sl = sl[sl.tier == "A"]
    at = pd.read_parquet(F("atlas.parquet")).set_index(["material_id", "T"])
    comps = []
    for r in sl.itertuples():
        a = at.loc[(r.material_id, 293.15)]
        comps.append({"material_id": r.material_id, "name": f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]",
                      "fractions": {str(e): float(x) for e, x in Composition(r.formula).fractional_composition.items()},
                      "n_atom_per_b_cm": float(a.numberdensity_A3)})
    out = D("stage4_activation")
    (out / "inputs.json").write_text(json.dumps({"compounds": comps}, indent=1), encoding="utf-8")
    print(len(comps), "compounds")


def run(adir, wdir):
    import openmc
    import openmc.deplete as od
    adir, wdir = Path(adir).resolve(), Path(wdir).resolve()
    chain_file = adir / "chain_endfb80_thermal.xml"
    chain = od.Chain.from_xml(chain_file)
    half = {n.name: n.half_life for n in chain.nuclides}
    # 293 K Maxwellian flux per unit energy, phi(E) ~ E exp(-E / kT), on 400 log groups from 1e-5 eV to 10 eV
    edges = np.logspace(-5, 1, 401)
    e_mid = np.sqrt(edges[1:] * edges[:-1])
    kt = K_B_EV * 293.15
    phi = e_mid * np.exp(-e_mid / kt) * np.diff(edges)
    phi = phi / phi.sum()
    comps = json.loads((adir / "inputs.json").read_text())["compounds"]
    days = C["irradiation_days"]
    cool = C["cooling_days"]                     # [0, 1, 7]: end of irradiation, +1 d, +7 d
    steps = [days] + [b - a for a, b in zip(cool[:-1], cool[1:])]
    rates = [1.0] + [0.0] * (len(steps) - 1)
    results = []
    for c in comps:
        wd = wdir / c["material_id"]
        wd.mkdir(parents=True, exist_ok=True)
        mat = openmc.Material(name=c["name"])
        for el, x in c["fractions"].items():
            mat.add_element(el, x, "ao")
        mat.set_density("atom/b-cm", c["n_atom_per_b_cm"])
        mat.volume = C["volume_cm3"]
        mat.depletable = True
        nucs = [n for n, _ in mat.get_nuclide_atom_densities().items()]
        micro = od.MicroXS.from_multigroup_flux(energies=edges, multigroup_flux=phi, chain_file=chain_file,
                                                nuclides=nucs, temperature=294.0)
        op = od.IndependentOperator(openmc.Materials([mat]), fluxes=[C["flux_n_cm2_s"] * C["volume_cm3"]],
                                    micros=[micro], chain_file=chain_file, normalization_mode="source-rate")
        integ = od.PredictorIntegrator(op, steps, source_rates=rates, timestep_units="d")
        integ.integrate(path=str(wd / "depletion_results.h5"))
        res = od.Results(wd / "depletion_results.h5")
        times, act = res.get_activity(mat, units="Bq", by_nuclide=True)
        rec = {"material_id": c["material_id"], "name": c["name"], "times_d": [float(t) / 86400 for t in times], "at": {}}
        for i, t in enumerate(times):
            if i == 0:
                continue
            a = act[i]
            tot = float(sum(a.values()))
            top = sorted(((k, float(v)) for k, v in a.items() if v > 0), key=lambda kv: -kv[1])[:3]
            label = f"{cool[i - 1]}d_after_end" if i >= 1 else "start"
            rec["at"][label] = {"total_Bq": tot, "top3": [{"nuclide": k, "Bq": v, "half_life_s": half.get(k)} for k, v in top]}
        results.append(rec)
        print(c["name"], {k: round(v["total_Bq"], 3) for k, v in rec["at"].items()}, flush=True)
    sha = hashlib.sha256(chain_file.read_bytes()).hexdigest()
    (adir / "results.json").write_text(json.dumps({"openmc": openmc.__version__, "chain_sha256": sha,
                                                   "flux_n_cm2_s": C["flux_n_cm2_s"], "irradiation_days": days,
                                                   "cooling_days": cool, "volume_cm3": C["volume_cm3"],
                                                   "compounds": results}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare()
    else:
        run(sys.argv[2], sys.argv[3])
