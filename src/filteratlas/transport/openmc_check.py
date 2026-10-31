"""OpenMC check of the atlas: uncollided vs transported A_fast on three slabs.
Runs in the OpenMC environment (Linux; OpenMC 0.16.0, same ENDF/B-VIII.0 NNDC data, free atom), with
OPENMC_CROSS_SECTIONS=data/upstream/endfb80_hdf5/cross_sections.xml:
  python pipelines/08_atlas.py --with openmc
  or: PYTHONPATH=src python -m filteratlas.transport.openmc_check <in.json> <out.json> <workdir>
The input is written by filteratlas.transport.openmc_input.

Slab of thickness L0.8(293 K), laterally 1e4 cm (vacuum), pencil beam along +x; sources 1/E (1 eV - 100 keV) and
Watt (a = 0.988 MeV, b = 2.249 /MeV, 100 keV - 10 MeV, tabulated). Tally: neutron current through the back face,
with CollisionFilter([0]) (uncollided) and without (all transmitted neutrons). URR probability tables are off so
the uncollided result uses the same smooth cross sections as the analytic A_fast.
in.json: [{"name", "fractions": {element: atom fraction}, "n_A3", "L_cm", "A_fast_1E", "A_fast_watt"}]
"""
import json
import sys
from pathlib import Path

import numpy as np
import openmc

from filteratlas.common.config import cfg

_OC = cfg("figures_of_merit")["openmc_check"]
_FOM = cfg("figures_of_merit")
N_PART, N_BATCH = _OC["particles_per_batch"], _OC["batches"]


def spectrum(kind):
    if kind == "1E":   # tabulated 1/E (PowerLaw with exponent -1 is the singular case of power-law sampling)
        e = np.logspace(0, 5, _OC["source_points_1E"])
        return openmc.stats.Tabular(e, 1.0 / e, "linear-linear")
    e = np.logspace(5, 7, _OC["source_points_watt"])
    p = np.exp(-e / _FOM["watt_a_eV"]) * np.sinh(np.sqrt(_FOM["watt_b_per_eV"] * e))
    return openmc.stats.Tabular(e, p, "linear-linear")


def run(c, kind, wd):
    wd.mkdir(parents=True, exist_ok=True)
    mat = openmc.Material(name=c["name"])
    for el, x in c["fractions"].items():
        mat.add_element(el, float(x), "ao")
    mat.set_density("atom/b-cm", float(c["n_A3"]))   # 1 atom/A^3 = 1 atom/(b cm)
    L = float(c["L_cm"])
    x0, x1 = openmc.XPlane(0.0), openmc.XPlane(L)
    xl, xr = openmc.XPlane(-1.0, boundary_type="vacuum"), openmc.XPlane(L + 1.0, boundary_type="vacuum")
    ylo, yhi = openmc.YPlane(-1e4, boundary_type="vacuum"), openmc.YPlane(1e4, boundary_type="vacuum")
    zlo, zhi = openmc.ZPlane(-1e4, boundary_type="vacuum"), openmc.ZPlane(1e4, boundary_type="vacuum")
    lat = +ylo & -yhi & +zlo & -zhi
    cells = [openmc.Cell(region=+xl & -x0 & lat), openmc.Cell(fill=mat, region=+x0 & -x1 & lat),
             openmc.Cell(region=+x1 & -xr & lat)]
    geom = openmc.Geometry(openmc.Universe(cells=cells))
    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.particles, s.batches = N_PART, N_BATCH
    s.source = openmc.IndependentSource(space=openmc.stats.Point((-0.5, 0, 0)),
                                        angle=openmc.stats.Monodirectional((1, 0, 0)), energy=spectrum(kind))
    s.temperature = {"default": 294.0}
    s.ptables = False
    s.output = {"tallies": False, "summary": False}
    sf = openmc.SurfaceFilter([x1])
    t0 = openmc.Tally(name="uncollided")
    t0.filters = [sf, openmc.CollisionFilter([0])]
    t0.scores = ["current"]
    t1 = openmc.Tally(name="all")
    t1.filters = [sf]
    t1.scores = ["current"]
    model = openmc.Model(geometry=geom, materials=openmc.Materials([mat]), settings=s, tallies=openmc.Tallies([t0, t1]))
    sp = model.run(cwd=str(wd), output=False)
    with openmc.StatePoint(sp) as st:
        u = st.get_tally(name="uncollided")
        a = st.get_tally(name="all")
        tu, su = float(u.mean.ravel()[0]), float(u.std_dev.ravel()[0])
        ta, sa = float(a.mean.ravel()[0]), float(a.std_dev.ravel()[0])
    return {"T_uncollided": tu, "T_uncollided_sd": su, "A_uncollided": -np.log(tu), "A_uncollided_sd": su / tu,
            "T_all": ta, "T_all_sd": sa, "A_all": -np.log(ta)}


def main(in_json=None, out_json=None, workdir=None):
    in_json, out_json, workdir = in_json or sys.argv[1], out_json or sys.argv[2], workdir or sys.argv[3]
    comps = json.loads(Path(in_json).read_text())
    wd = Path(workdir)
    out = []
    for c in comps:
        rec = dict(c)
        for kind in ("1E", "watt"):
            r = run(c, kind, wd / f"{c['name']}_{kind}")
            ana = c[f"A_fast_{kind}"]
            rec[kind] = {**r, "A_analytic": ana, "rel_diff_uncollided": (r["A_uncollided"] - ana) / ana,
                         "pass": bool(abs(r["A_uncollided"] - ana) / ana < cfg("uncertainty")["sigma_model_dft"]),
                         "buildup_ratio_T_all_over_T_uncollided": r["T_all"] / r["T_uncollided"]}
            print(c["name"], kind, rec[kind], flush=True)
        out.append(rec)
    Path(out_json).write_text(json.dumps({"openmc": openmc.__version__, "particles": N_PART * N_BATCH,
                                             "compounds": out}, indent=1, default=float))


if __name__ == "__main__":
    main()
