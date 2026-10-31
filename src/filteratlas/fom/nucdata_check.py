"""Convergence check of the A_fast energy grid: the 5,000-points-per-decade grid vs the exact
union grid of the isotopes, for Al2O3 (sapphire), Bi and CaWO4, at L = 5 cm and 20 cm, both spectra. Pass: < 0.1 %.
-> the nuclear-data check JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np

from filteratlas.fom import rejection as rj  # noqa: E402
from filteratlas.common.config import cfg
from filteratlas.common.paths import D, F  # noqa: E402

_NC = cfg("figures_of_merit")["nucdata_check"]
CASES = _NC["cases_atoms_per_A3"]   # rough atom densities (atoms/A^3); only for the check


def exact(formula, n, L, spectrum):
    z = np.load(D("nucdata") / "endfb80_elements.npz")
    fr = rj.fractions(formula)
    iso = {}
    for k in z.files:
        if k.startswith("check::") and k.endswith("|E"):
            el, nuc = k.split("::")[1].split("|")[:2]
            iso.setdefault(el, []).append((z[f"check::{el}|{nuc}|E"], z[f"check::{el}|{nuc}|xs"], float(z[f"check::{el}|{nuc}|a"])))
    lo, hi = cfg("figures_of_merit")["spectrum_1E_eV" if spectrum == "1E" else "spectrum_watt_eV"]
    E = np.unique(np.concatenate([e for v in iso.values() for e, _, _ in v] + [np.array([lo, hi])]))
    E = E[(E >= lo) & (E <= hi)]
    sig = sum(x * sum(a * np.interp(E, e, s) for e, s, a in iso[el]) for el, x in fr.items())
    S = 1 / E if spectrum == "1E" else np.exp(-E / rj.WATT_A) * np.sinh(np.sqrt(rj.WATT_B * E))
    return -np.log(np.trapezoid(S * np.exp(-n * sig * L), E) / np.trapezoid(S, E)), len(E)


def main():
    out, ok = [], True
    for f, n in CASES.items():
        for L in _NC["L_cm"]:
            for sp in ("1E", "watt"):
                a_grid = float(rj.a_fast(rj.sigma_tot(f), n, L, sp)[0])
                a_ex, npts = exact(f, n, L, sp)
                rel = (a_grid - a_ex) / a_ex
                ok &= abs(rel) < _NC["tolerance"]
                out.append({"formula": f, "L_cm": L, "spectrum": sp, "A_grid": a_grid, "A_exact": float(a_ex),
                            "rel_diff": float(rel), "n_union_points": npts})
    res = {"cases": out, "max_abs_rel_diff": max(abs(r["rel_diff"]) for r in out), "pass": bool(ok)}
    (F("stage3_nucdata_check.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
