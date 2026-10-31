"""Uncollided fast-neutron and gamma rejection, from the elemental ENDF/B-VIII.0 tables
written by filteratlas.ingest.nucdata (data/intermediate/nucdata/endfb80_elements.npz).

A_fast(L) = -ln[ int S(E) exp(-Sigma_tot(E) L) dE / int S(E) dE ], Sigma_tot = n sum_i x_i sigma_tot,i(E)
  spectra: 1/E over 1 eV - 100 keV; Watt (a = 0.988 MeV, b = 2.249 /MeV) over 100 keV - 10 MeV.
A_gamma = mu_gamma(E) L at 1 and 6 MeV. mu_epi = n sum_i x_i sigma_tot,i(1 eV).
Units: n in atoms/A^3 and sigma in barn give cm^-1.
"""
from functools import lru_cache

import numpy as np
from pymatgen.core import Composition

from filteratlas.common.config import cfg
from filteratlas.common.paths import D

_FOM = cfg("figures_of_merit")
NUC = D("nucdata") / "endfb80_elements.npz"
WATT_A, WATT_B = _FOM["watt_a_eV"], _FOM["watt_b_per_eV"]     # eV, 1/eV


@lru_cache(maxsize=1)
def tables():
    z = np.load(NUC)
    E = z["grid"]
    el = {str(e): i for i, e in enumerate(z["elements"])}
    ph = {str(e): i for i, e in enumerate(z["photon_elements"])}
    spectra = {}
    (lo1, hi1), (low, hiw) = _FOM["spectrum_1E_eV"], _FOM["spectrum_watt_eV"]
    for name, lo, hi, s in (("1E", lo1, hi1, 1.0 / E), ("watt", low, hiw,
                                                        np.exp(-E / WATT_A) * np.sinh(np.sqrt(WATT_B * E)))):
        m = (E >= lo) & (E <= hi)
        e, w = E[m], s[m]
        # trapezoid weights so that sum(wt * f) = int S f dE on the sub-grid
        dw = np.zeros_like(e)
        de = np.diff(e)
        dw[:-1] += 0.5 * de * w[:-1]
        dw[1:] += 0.5 * de * w[1:]
        spectra[name] = (m, dw / dw.sum())
    return {"E": E, "sig": z["sigma_tot"], "el": el, "phsig": z["photon_sigma"], "ph": ph,
            "phE": z["photon_E"], "spectra": spectra}


def fractions(formula):
    return {str(e): float(x) for e, x in Composition(formula).fractional_composition.items()}


def sigma_tot(formula):
    """Mean free-atom total cross section per atom [b] on the energy grid."""
    t = tables()
    return sum(x * t["sig"][t["el"][e]] for e, x in fractions(formula).items())


def a_fast(sig, n, L, spectrum):
    """-ln uncollided transmission for spectrum '1E' or 'watt' through L cm (n atoms/A^3, sig barn on grid)."""
    m, w = tables()["spectra"][spectrum]
    L = np.atleast_1d(np.asarray(L, float))
    T = np.array([float(np.dot(w, np.exp(-n * sig[m] * x))) for x in L])
    return -np.log(T)


def mu_gamma(formula, n):
    """(mu at 1 MeV, mu at 6 MeV) [cm^-1], total photon cross section incl. coherent (narrow-beam)."""
    t = tables()
    s = sum(x * t["phsig"][t["ph"][e]] for e, x in fractions(formula).items())
    return float(n * s[0]), float(n * s[1])


def mu_epi(formula, n):
    """Free-atom total attenuation at 1 eV [cm^-1] (Nieman-type ratio numerator)."""
    return float(n * sigma_tot(formula)[0])
