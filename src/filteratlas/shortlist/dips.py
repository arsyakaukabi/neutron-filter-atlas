"""Bragg-dip density and best orientation for the Tier-A shortlist.

NCrystal single-crystal model (Gaussian mosaicity) built from the atlas structure and VDOS at 293 K, coherent elastic
only. The crystal is fixed in the lab frame: NCrystal's Cartesian crystal frame (a along x, b in the xy plane, z along
c*) is mapped onto the lab frame (dir1=@crys:0,0,1@lab:0,0,1, dir2=@crys:1,0,0@lab:1,0,0); orientations
are sampled as beam directions on a Fibonacci grid over the upper hemisphere, plus the literature directions (sapphire
beam || c; Si along <111>, <100>, <110>).
Per direction: T_B(lambda) = exp(-n sigma_Bragg(lambda) L0.8(293 K)) on 6,000 wavelengths over 0.5-3.0 A; a dip is a
contiguous interval with T_B < 0.9; metrics: dips per A, flux-weighted blocked fraction (300 K Maxwellian), longest
dip-free interval. Best orientation: smallest blocked fraction, ties broken by the longest dip-free interval.
Sapphire check: dip centres for beam || c against Mildner et al. (1993) eq. (B5).
The same metrics along every low-index lattice direction [uvw] (|u|,|v|,|w| <= 2) for every
compound, so the zone axes are sampled for all crystals and not only for sapphire and Si; the best orientation is
reported over the grid alone and over grid + low-index directions. Sapphire tilt scan about c (0-3 deg, 12 azimuths)
to test the orientation statements of Nieman et al. and Born et al. as quoted by Mildner et al. (1993), for a pencil
beam and averaged over a uniform cone of beam divergence (0.25 and 0.5 deg half-angle).
-> derived/stage4_dips.json
"""
import json
import time

import numpy as np
import NCrystal as NC
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import F
from filteratlas.shortlist.transport_inputs import dos_for
from filteratlas.xs import evaluate as xe
from filteratlas.xs import ncmat as nb

C = cfg("shortlist")["dips"]
LIT_DIRS = {"Al2O3": {"beam || c [001]": "c"}, "Si": {"<111>": (1, 1, 1), "<100>": (1, 0, 0), "<110>": (1, 1, 0)}}


def lattice_matrix(si):
    """Rows: the direct lattice vectors a, b, c in the lab frame of scatter(), which is NCrystal's Cartesian crystal
    frame (a along x, b in the xy plane, z along c*; this differs from c along z only for cells with alpha or
    beta != 90 deg, here monoclinic ZrO2)."""
    a, b, c = si["a"], si["b"], si["c"]
    al, be, ga = np.radians([si["alpha"], si["beta"], si["gamma"]])
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    return np.array([va, vb, [cx, cy, np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))]])


def low_index_directions(si, nmax=2):
    """Unit lab vectors of the lattice directions [uvw], |u|,|v|,|w| <= nmax, gcd 1, one of each +- pair."""
    from math import gcd
    M = lattice_matrix(si)
    out = {}
    r = range(-nmax, nmax + 1)
    for u in r:
        for v in r:
            for w in r:
                if (u, v, w) == (0, 0, 0) or gcd(gcd(abs(u), abs(v)), abs(w)) != 1:
                    continue
                if next(x for x in (u, v, w) if x) < 0:   # one of each +- pair: first non-zero index positive
                    continue
                vec = u * M[0] + v * M[1] + w * M[2]
                out[(u, v, w)] = vec / np.linalg.norm(vec)
    return {f"[{u} {v} {w}]": d for (u, v, w), d in out.items()}


def cone(axis, half_deg, n_ring=3, n_az=12):
    """Directions filling a uniform cone (by solid angle) about a unit axis; the axis itself first."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    e1 = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    out = [axis]
    for k in range(1, n_ring + 1):
        th = np.radians(half_deg) * np.sqrt(k / n_ring)
        for j in range(n_az):
            ph = 2 * np.pi * (j + 0.5 * (k % 2)) / n_az
            out.append(np.cos(th) * axis + np.sin(th) * (np.cos(ph) * e1 + np.sin(ph) * e2))
    return np.array(out)


def tilt_scan(sc, ek, n, L, lam, phi, tilts=(0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0), n_az=12, divergences=(0.0, 0.25, 0.5)):
    """Blocked fraction for beam directions tilted by theta from c (azimuths around c), per beam divergence."""
    out = {}
    for div in divergences:
        rows = []
        for th in tilts:
            bfs = []
            for j in range(1 if th == 0 else n_az):
                ph = 2 * np.pi * j / n_az
                d = np.array([np.sin(np.radians(th)) * np.cos(ph), np.sin(np.radians(th)) * np.sin(ph), np.cos(np.radians(th))])
                dirs = [d] if div == 0 else cone(d, div)
                T = np.mean([np.exp(-n * np.array(sc.crossSection(ek, tuple(x))) * L) for x in dirs], axis=0)
                bfs.append(metrics(T, lam, phi)["blocked_fraction"])
            rows.append({"tilt_deg": th, "blocked_fraction_min": float(min(bfs)), "median": float(np.median(bfs)),
                         "max": float(max(bfs))})
        out[f"divergence_{div}deg"] = rows
    return out


def fibonacci_hemisphere(n):
    i = np.arange(n) + 0.5
    z = 1 - i / n                            # z in (0, 1]: upper hemisphere
    phi = np.pi * (1 + 5 ** 0.5) * i
    r = np.sqrt(1 - z ** 2)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def scatter(name, mosaic):
    cfgstr = (f"{name};temp={C['temperature_K']}K;mos={mosaic}deg;dir1=@crys:0,0,1@lab:0,0,1;"
              f"dir2=@crys:1,0,0@lab:1,0,0;dirtol=180deg;comp=coh_elas")
    return NC.createScatter(cfgstr)


def metrics(Tb, lam, phi):
    blocked = Tb < C["dip_transmission"]
    edges = np.diff(blocked.astype(int))
    n_dips = int((edges == 1).sum() + int(blocked[0]))
    free, best, run = ~blocked, 0.0, 0
    dl = lam[1] - lam[0]
    for f in free:
        run = run + 1 if f else 0
        best = max(best, run * dl)
    return {"n_dips": n_dips, "dips_per_A": n_dips / (lam[-1] - lam[0]),
            "blocked_fraction": float((phi * blocked).sum() / phi.sum()), "dip_free_A": float(best)}


def dip_centres(Tb, lam):
    blocked = Tb < C["dip_transmission"]
    out, start = [], None
    for i, b in enumerate(blocked):
        if b and start is None:
            start = i
        if not b and start is not None:
            out.append(float(lam[start:i][np.argmin(Tb[start:i])]))
            start = None
    return out


def mildner_b5(info, lam_band):
    """Mildner et al. 1993 eq. (B5), beam || c: lambda = 2 c l / [l^2 + (4 c^2 / 3 a^2)(h^2 + hk + k^2)]."""
    a, c = info.structure_info["a"], info.structure_info["c"]
    out = []
    for h in range(-8, 9):
        for k in range(-8, 9):
            for l in range(1, 40):
                lam = 2 * c * l / (l ** 2 + (4 * c ** 2 / (3 * a ** 2)) * (h * h + h * k + k * k))
                if lam_band[0] <= lam <= lam_band[1]:
                    out.append(lam)
    return np.array(sorted(set(np.round(out, 4))))


def main():
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    sl = sl[sl.tier == "A"]
    at = pd.read_parquet(F("atlas.parquet")).set_index(["material_id", "T"])
    cand = pd.read_parquet(F("stage0_candidates.parquet")).set_index("material_id")
    lam = np.linspace(xe.BAND[0], xe.BAND[1], C["n_wavelengths"])
    phi = xe.maxwellian_flux(lam)
    ek = NC.wl2ekin(lam)
    dirs = fibonacci_hemisphere(C["n_directions"])
    out = {"settings": {k: C[k] for k in C}, "compounds": {}}
    for r in sl.itertuples():
        t0 = time.time()
        name = nb.build_vdos_material(dos_for(r.material_id, r.vdos_source, cand), f"dip_{r.material_id}")
        L = float(at.loc[(r.material_id, 293.15)].L08)
        info = NC.createInfo(f"{name};temp={C['temperature_K']}K")
        n = info.numberdensity
        rec = {"formula": r.formula, "sg_symbol": r.sg_symbol, "L08_293K_cm": L, "mosaic": {}}
        for mos in [C["mosaic_deg"]] + list(C["mosaic_sensitivity_deg"]):
            sc = scatter(name, mos)
            res = []
            for d in dirs:
                xs = np.array(sc.crossSection(ek, tuple(d)))
                res.append(metrics(np.exp(-n * xs * L), lam, phi))
            bf = np.array([x["blocked_fraction"] for x in res])
            dfree = np.array([x["dip_free_A"] for x in res])
            order = np.lexsort((-dfree, bf))
            b = order[0]
            entry = {"n_directions": len(dirs),
                     "blocked_fraction": {"best": float(bf[b]), "median": float(np.median(bf)), "worst": float(bf.max())},
                     "dips_per_A": {"best_orientation": res[b]["dips_per_A"],
                                    "median": float(np.median([x["dips_per_A"] for x in res]))},
                     "dip_free_A": {"best_orientation": float(dfree[b]), "median": float(np.median(dfree))},
                     "best_direction_lab": [round(float(v), 4) for v in dirs[b]]}
            lit = {}
            for lab_name, v in LIT_DIRS.get(r.formula, {}).items():
                vec = (0.0, 0.0, 1.0) if v == "c" else tuple(np.array(v, float) / np.linalg.norm(v))
                xs = np.array(sc.crossSection(ek, vec))
                Tb = np.exp(-n * xs * L)
                lit[lab_name] = {**metrics(Tb, lam, phi), "rank_of_blocked_fraction_among_grid":
                                 int((bf < metrics(Tb, lam, phi)["blocked_fraction"]).sum()) + 1}
                if r.formula == "Al2O3" and v == "c" and mos == C["mosaic_deg"]:
                    cen = np.array(dip_centres(Tb, lam))
                    b5 = mildner_b5(info, (lam[0], lam[-1]))
                    near = np.array([float(np.min(np.abs(b5 - x)) / x) for x in cen]) if len(b5) else np.array([])
                    chk = {"n_model_dips": len(cen), "n_b5_lines_in_band": len(b5),
                           "model_dip_centres_A": [round(x, 4) for x in cen[:40]]}
                    # B5 gives a dense set of lines, so a match means something only against the chance level: the
                    # fraction of the band that lies within the same tolerance of some B5 line
                    for tol in (0.005, 0.001, 0.0002):
                        chance = float(np.mean([np.min(np.abs(b5 - x)) / x < tol for x in lam]))
                        chk[f"tol_{tol}"] = {"fraction_of_model_dips_matched": float(np.mean(near < tol)) if len(near) else None,
                                             "chance_level": chance}
                    lit[lab_name]["mildner_b5_check"] = chk
            entry["literature_directions"] = lit
            # low-index lattice directions for every compound
            low = {}
            for nm, vec in low_index_directions(info.structure_info).items():
                low[nm] = metrics(np.exp(-n * np.array(sc.crossSection(ek, tuple(vec))) * L), lam, phi)
            lb = min(low, key=lambda k: (low[k]["blocked_fraction"], -low[k]["dip_free_A"]))
            entry["low_index_directions"] = {"n": len(low), "best": lb, "best_blocked_fraction": low[lb]["blocked_fraction"],
                                             "median_blocked_fraction": float(np.median([x["blocked_fraction"] for x in low.values()])),
                                             "per_direction": {k: round(v["blocked_fraction"], 4) for k, v in low.items()}}
            entry["best_blocked_fraction_grid_and_low_index"] = float(min(bf[b], low[lb]["blocked_fraction"]))
            if r.formula == "Al2O3" and mos == C["mosaic_deg"]:
                entry["tilt_scan_about_c"] = tilt_scan(sc, ek, n, L, lam, phi)
            rec["mosaic"][f"{mos}deg"] = entry
        rec["seconds"] = round(time.time() - t0, 1)
        out["compounds"][f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]"] = rec
        print(r.formula, rec["seconds"], "s; best blocked", rec["mosaic"][f"{C['mosaic_deg']}deg"]["blocked_fraction"], flush=True)
    F("stage4_dips.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
