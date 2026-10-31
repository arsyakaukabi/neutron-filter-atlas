"""Nuclear and photon data for the fast-neutron and photon rejection figures of merit. Needs openmc
(environment.yml); constants from configs/figures_of_merit.yaml. Run by step 07:
  python pipelines/07_figures_of_merit.py --with openmc
or directly: python -m filteratlas.ingest.nucdata data/upstream/endfb80_hdf5 data/intermediate/nucdata

From the ENDF/B-VIII.0 (NNDC build, openmc_data_downloader) HDF5 files in <xs_dir>:
- elemental free-atom total cross section sigma_tot(E) at 294 K on a log grid of 5,000 points per decade over
  1 eV - 10 MeV: natural-abundance-weighted sum over isotopes (openmc.data.NATURAL_ABUNDANCE), lin-lin
  interpolation of each isotope's pointwise MT=1; missing isotopes -> abundances renormalised, element flagged
  nucdata_incomplete if > 1 % of the abundance is missing;
- the exact (union-grid) isotope data for the convergence-check elements;
- total photon cross section (coherent + incoherent + photoelectric + pair, nuclear and electron fields) at 1 and
  6 MeV from the photoatomic data.
-> <out_dir>/endfb80_elements.npz, <out_dir>/endfb80_meta.json
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import openmc.data as od

from filteratlas.common.config import cfg, tup

_F = cfg("figures_of_merit")
EMIN, EMAX, PER_DECADE = _F["nucdata_grid"]["emin_eV"], _F["nucdata_grid"]["emax_eV"], _F["nucdata_grid"]["per_decade"]
TEMP = _F["nucdata_temperature"]
CHECK_ELEMENTS = tup(_F["nucdata_check_elements"])
PHOTON_E = tup(_F["photon_energies_eV"])
PHOTON_MT = tup(_F["photon_mt"])   # coherent, incoherent, pair (nuclear), pair (electron), photoelectric


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def isotope_total(path):
    nuc = od.IncidentNeutron.from_hdf5(str(path))
    e = np.asarray(nuc.energy[TEMP])
    xs = nuc[1].xs[TEMP]
    y = np.asarray(xs(e) if callable(xs) else xs.y)
    return e, y


def main(xs_dir=None, out_dir=None):
    xsdir = Path(xs_dir or sys.argv[1]).expanduser()
    out = Path(out_dir or sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    grid = np.logspace(np.log10(EMIN), np.log10(EMAX), int(PER_DECADE * np.log10(EMAX / EMIN)) + 1)
    files = {p.stem.replace("ENDFB-8.0-NNDC_", ""): p for p in xsdir.glob("ENDFB-8.0-NNDC_*.h5")}
    elements = sorted(k for k in files if k.isalpha())   # photon files are named by element symbol (one per element)
    sig, meta, check = {}, {"library": "ENDFB-8.0-NNDC (openmc_data_downloader)", "temperature": TEMP,
                            "grid": {"emin_eV": EMIN, "emax_eV": EMAX, "per_decade": PER_DECADE, "n": len(grid)},
                            "elements": {}, "files_sha256": {}}, {}
    photon = {}
    for el in elements:
        iso = {n: a for n, a in od.NATURAL_ABUNDANCE.items() if "".join(c for c in n if c.isalpha()) == el}
        have = {n: a for n, a in iso.items() if n in files}
        missing = sum(a for n, a in iso.items() if n not in files)
        rec = {"isotopes": sorted(have), "missing_isotopes": sorted(set(iso) - set(have)),
               "missing_abundance": round(missing, 6), "nucdata_incomplete": missing > 0.01}
        if have:
            tot = sum(have.values())
            s = np.zeros_like(grid)
            for n, a in have.items():
                e, y = isotope_total(files[n])
                s += a / tot * np.interp(grid, e, y)
                meta["files_sha256"][n] = sha(files[n])
                if el in CHECK_ELEMENTS:
                    m = (e >= EMIN) & (e <= EMAX)
                    check[f"{el}|{n}|E"], check[f"{el}|{n}|xs"], check[f"{el}|{n}|a"] = e[m], y[m], np.array(a / tot)
            sig[el] = s
            rec["sigma_tot_1eV_b"] = float(s[0])
        if el in files:   # photoatomic file for the element
            ph = od.IncidentPhoton.from_hdf5(str(files[el]))
            vals = []
            for E in PHOTON_E:
                v = 0.0
                for mt in PHOTON_MT:
                    if mt in ph.reactions:
                        v += float(ph[mt].xs(E))
                vals.append(v)
            photon[el] = vals
            rec["photon_sigma_b_1MeV_6MeV"] = vals
            meta["files_sha256"][el] = sha(files[el])
        meta["elements"][el] = rec
    els = sorted(sig)
    np.savez_compressed(out / "endfb80_elements.npz", grid=grid, elements=np.array(els),
                        sigma_tot=np.array([sig[e] for e in els]),
                        photon_elements=np.array(sorted(photon)), photon_sigma=np.array([photon[e] for e in sorted(photon)]),
                        photon_E=np.array(PHOTON_E), **{f"check::{k}": v for k, v in check.items()})
    (out / "endfb80_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"neutron elements {len(els)}, photon elements {len(photon)}, incomplete "
          f"{[e for e, r in meta['elements'].items() if r['nucdata_incomplete']]}")


if __name__ == "__main__":
    main()
