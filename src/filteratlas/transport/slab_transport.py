"""Energy-resolved slab transport of the shortlist. Runs in the OpenMC environment (Linux; OpenMC 0.16.0,
NCrystal 4.4.6, ENDF/B-VIII.0 NNDC HDF5), needs no pymatgen; OPENMC_CROSS_SECTIONS=data/upstream/endfb80_hdf5/cross_sections.xml:
  python pipelines/10_shortlist.py --with openmc
  or: PYTHONPATH=src python -m filteratlas.transport.slab_transport <transport dir> <work dir>

Per case (compound, temperature, thickness): a slab in vacuum, laterally 1e4 cm, pencil beam at normal incidence.
Material: the compound's atlas NCMAT through NCrystal with coherent elastic (Bragg) scattering off, as in the atlas
mu_th; above the NCrystal range ENDF/B-VIII.0 free atom at 294 K (nearest temperature, as in the atlas).
Sources: thermal (300 K Maxwellian per unit wavelength over the atlas band), 1/E (1 eV - 100 keV), Watt
(100 keV - 10 MeV), photons at 1 and 6 MeV.
Tallies: current through the back face, uncollided (CollisionFilter([0])) and all, in energy groups
[0, 0.625 eV, 100 keV, 20 MeV] for neutrons and [10 keV, 20 MeV] for photons. Histories are raised by a factor 4 from
1e6 (up to 6.4e7) until every group that enters a figure of merit has a relative SD below 2 %.
-> <transport dir>/raw/<case>.json (resumable), read by filteratlas.shortlist.transport_results
"python -m filteratlas.transport.slab_transport all_components_1E <transport dir> <work dir>" reruns the 1/E
source with all NCrystal components (coherent elastic kept) -> <transport dir>/raw_1E_all/.
"""
import json
import sys
from pathlib import Path

import numpy as np
import openmc

from filteratlas.common.config import cfg

C = cfg("shortlist")["transport"]
FOM = cfg("figures_of_merit")
H, M_N, EV = 6.62607015e-34, 1.67492749804e-27, 1.602176634e-19


def source_energy(kind, lam_t, band):
    if kind == "thermal":   # phi(lambda) ~ lambda^-5 exp(-(lambda_T/lambda)^2) over the band, as p(E) on an energy grid
        lam = np.linspace(band[0], band[1], 4000)
        e = (H / (lam * 1e-10)) ** 2 / (2 * M_N) / EV
        p = lam ** -5 * np.exp(-(lam_t / lam) ** 2) * lam ** 3   # |d lambda / dE| ~ lambda^3
        o = np.argsort(e)
        return openmc.stats.Tabular(e[o], p[o], "linear-linear")
    if kind == "1E":
        e = np.logspace(0, 5, 2000)
        return openmc.stats.Tabular(e, 1.0 / e, "linear-linear")
    if kind == "watt":
        e = np.logspace(5, 7, 2000)
        return openmc.stats.Tabular(e, np.exp(-e / FOM["watt_a_eV"]) * np.sinh(np.sqrt(FOM["watt_b_per_eV"] * e)),
                                    "linear-linear")
    return openmc.stats.Discrete([{"g1MeV": 1e6, "g6MeV": 6e6}[kind]], [1.0])


def run_once(case, kind, ncmat_dir, wd, n_particles, lam_t, band, free_atom=False, all_components=False):
    wd.mkdir(parents=True, exist_ok=True)
    # all_components: coherent elastic kept, for the 1/E source; with "Bragg off" extended to the whole NCrystal
    # range (up to 5 eV) the material misses the coherent-elastic part of the free-atom cross section above 1 eV
    comp = "" if all_components else ";comp=inelas,incoh_elas"
    mat = openmc.Material.from_ncrystal(f"{ncmat_dir / case['ncmat']};temp={case['T']}K{comp}")
    if free_atom:   # diagnostic: the same nuclides and density, free atom everywhere, as the analytic A_fast
        free = openmc.Material(name="free")
        dens = mat.get_nuclide_atom_densities()
        for nuc, d in dens.items():
            free.add_nuclide(nuc, d)
        free.set_density("atom/b-cm", sum(dens.values()))
        mat = free
    L = float(case["L_cm"])
    lat_ = C["lateral_cm"]
    x0, x1 = openmc.XPlane(0.0), openmc.XPlane(L)
    xl, xr = openmc.XPlane(-1.0, boundary_type="vacuum"), openmc.XPlane(L + 1.0, boundary_type="vacuum")
    lat = (+openmc.YPlane(-lat_, boundary_type="vacuum") & -openmc.YPlane(lat_, boundary_type="vacuum")
           & +openmc.ZPlane(-lat_, boundary_type="vacuum") & -openmc.ZPlane(lat_, boundary_type="vacuum"))
    cells = [openmc.Cell(region=+xl & -x0 & lat), openmc.Cell(fill=mat, region=+x0 & -x1 & lat),
             openmc.Cell(region=+x1 & -xr & lat)]
    photon = kind.startswith("g")
    s = openmc.Settings()
    s.run_mode = "fixed source"
    s.batches = C["batches"]
    s.particles = int(n_particles // C["batches"])
    s.source = openmc.IndependentSource(space=openmc.stats.Point((-0.5, 0, 0)), angle=openmc.stats.Monodirectional((1, 0, 0)),
                                        energy=source_energy(kind, lam_t, band), particle="photon" if photon else "neutron")
    s.photon_transport = photon
    s.temperature = {"default": 294.0, "method": "nearest", "tolerance": 300.0}
    s.ptables = False
    s.output = {"tallies": False, "summary": False}
    sf = openmc.SurfaceFilter([x1])
    groups = [C["photon_cut_eV"], 2e7] if photon else [0.0, C["e_thermal_cut_eV"], C["e_fast_cut_eV"], 2e7]
    pf = openmc.ParticleFilter(["photon" if photon else "neutron"])
    t0 = openmc.Tally(name="uncollided")
    # photons: secondary photons (annihilation, bremsstrahlung, fluorescence) also start with zero collisions;
    # the uncollided primaries are those with zero collisions at the source energy (+-0.01 %)
    e0 = {"g1MeV": 1e6, "g6MeV": 6e6}.get(kind)
    t0.filters = [sf, pf, openmc.CollisionFilter([0]),
                  openmc.EnergyFilter([e0 * 0.9999, e0 * 1.0001] if photon else groups)]
    t0.scores = ["current"]
    t1 = openmc.Tally(name="all")
    t1.filters = [sf, pf, openmc.EnergyFilter(groups)]
    t1.scores = ["current"]
    model = openmc.Model(geometry=openmc.Geometry(openmc.Universe(cells=cells)), materials=openmc.Materials([mat]),
                         settings=s, tallies=openmc.Tallies([t0, t1]))
    sp = model.run(cwd=str(wd), output=False)
    with openmc.StatePoint(sp) as st:
        u, a = st.get_tally(name="uncollided"), st.get_tally(name="all")
        um, us = u.mean.ravel(), u.std_dev.ravel()
        am, as_ = a.mean.ravel(), a.std_dev.ravel()
    return {"groups_eV": groups, "uncollided": um.tolist(), "uncollided_sd": us.tolist(), "all": am.tolist(),
            "all_sd": as_.tolist(), "particles": int(s.particles * s.batches)}


def fom_groups(kind, r):
    """(value, sd) of the transmitted fractions that enter a figure of merit."""
    u, us, a, as_ = map(np.array, (r["uncollided"], r["uncollided_sd"], r["all"], r["all_sd"]))
    if kind.startswith("g"):
        return {"uncollided": (u.sum(), np.sqrt((us ** 2).sum())), "transported": (a.sum(), np.sqrt((as_ ** 2).sum()))}
    if kind == "thermal":
        return {"uncollided": (u.sum(), np.sqrt((us ** 2).sum()))}
    return {"uncollided": (u.sum(), np.sqrt((us ** 2).sum())),
            "transported": (a[1:].sum(), np.sqrt((as_[1:] ** 2).sum()))}   # all transmitted above 0.625 eV


def main(tdir=None, workdir=None, all_components_1E=False):
    """all_components_1E: rerun only the 1/E source with all NCrystal components -> <transport dir>/raw_1E_all/"""
    tdir, workdir = Path(tdir or sys.argv[1]).resolve(), Path(workdir or sys.argv[2]).resolve()   # OpenMC runs in its own cwd
    meta = json.loads((tdir / "cases.json").read_text())
    raw = tdir / ("raw_1E_all" if all_components_1E else "raw")
    raw.mkdir(exist_ok=True)
    for case in meta["cases"]:
        for kind in (("1E",) if all_components_1E else ("thermal", "1E", "watt", "g1MeV", "g6MeV")):
            key = f"{case['material_id']}_{case['tag']}_{case['thickness']}_{kind}"
            f = raw / f"{key}.json"
            if f.exists():
                continue
            n = C["particles_start"]
            while True:
                r = run_once(case, kind, tdir / "ncmat", workdir / (key + ("_all" if all_components_1E else "")), n,
                             meta["lambda_T_A"], meta["band_A"], all_components=all_components_1E)
                g = fom_groups(kind, r)
                rel = {k: (sd / v if v > 0 else np.inf) for k, (v, sd) in g.items()}
                if max(rel.values()) < C["rel_sd_target"] or n * C["particles_factor"] > C["particles_cap"]:
                    break
                n *= C["particles_factor"]
            out = {**{k: case[k] for k in ("material_id", "formula", "sg_symbol", "T", "tag", "thickness", "L_cm")},
                   "source": kind, "analytic": case["analytic"], "raw": r,
                   "fom_fractions": {k: [float(v), float(sd)] for k, (v, sd) in g.items()},
                   "rel_sd": {k: float(x) for k, x in rel.items()},
                   "rel_sd_target_met": bool(max(rel.values()) < C["rel_sd_target"]), "openmc": openmc.__version__,
                   "ncrystal_components": "all" if all_components_1E else "inelas,incoh_elas"}
            f.write_text(json.dumps(out, indent=1))
            print(key, {k: round(v[0], 5) for k, v in out["fom_fractions"].items()}, out["rel_sd"], flush=True)




def diagnose(tdir, workdir, keys):
    """Free-atom diagnostic: rerun the given cases (material_id_tag_thickness_source) with a free-atom material, uncollided
    only, and compare with the analytic value. -> <transport dir>/diagnostic_free_atom.json"""
    tdir, workdir = Path(tdir).resolve(), Path(workdir).resolve()
    meta = json.loads((tdir / "cases.json").read_text())
    by = {f"{c['material_id']}_{c['tag']}_{c['thickness']}": c for c in meta["cases"]}
    f = tdir / "diagnostic_free_atom.json"
    out = [d for d in (json.loads(f.read_text()) if f.exists() else []) if d["case"] not in keys]   # merge with earlier runs
    for k in keys:
        base, kind = k.rsplit("_", 1)
        c = by[base]
        r = run_once(c, kind, tdir / "ncmat", workdir / f"diag_{k}", 4 * C["particles_start"], meta["lambda_T_A"],
                     meta["band_A"], free_atom=True)
        u = float(sum(r["uncollided"]))
        ana = c["analytic"]["A_fast_1E" if kind == "1E" else "A_fast_watt"]
        out.append({"case": k, "A_uncollided_free_atom": -float(np.log(u)), "A_analytic": ana,
                    "rel_diff": -float(np.log(u)) / ana - 1})
        print(out[-1], flush=True)
    (tdir / "diagnostic_free_atom.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[1] == "diagnose":
        diagnose(sys.argv[2], sys.argv[3], sys.argv[4].split(","))
    elif len(sys.argv) > 3 and sys.argv[1] == "all_components_1E":
        main(sys.argv[2], sys.argv[3], all_components_1E=True)
    else:
        main()
