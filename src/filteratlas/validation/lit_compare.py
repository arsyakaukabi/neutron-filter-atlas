"""Literature cross-check of the benchmark, sapphire (step 05, validation).

Literature inputs (transcribed with page references):
  Mildner, Arif, Stone & Crawford, J. Appl. Cryst. 26 (1993) 438, hemex grade, ambient temperature,
    Eq. (3) Table 2:  Sigma = A lam + C (1 - exp[-(B/lam^2 + D/lam^4)])
    Eq. (4) Table 3:  Sigma = A lam + C [1 - (lam^2/2B)(1 - exp[-(2B/lam^2 + D/lam^4)])]
    (Sigma in m^-1, lam in nm), N = 2.35e19 Al2O3 units per mm^3 (p. 441);
  Tennant (1988) as quoted by Mildner et al. p. 438: transmission at 0.12 nm through 101.6 mm "increases by only 18%"
    on cooling from 300 to 77 K;
  Born et al. (1987) as quoted p. 439: a 90 mm filter "has a transmission of about 0.8" for 0.12-0.24 nm.
Model: MP DFT VDOS (atlas model), off-Bragg baseline, per atom -> per Al2O3 unit (x5); n from the MP structure.
-> the sapphire literature-comparison JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np

from filteratlas.validation import exfor_benchmark as S  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

HEMEX = {"eq3": {"A": 5.96, "C": 25.2, "B": 1.17e-3, "D": 3.42e-6},      # Table 2, p. 442
         "eq4": {"A": 6.03, "C": 28.4, "B2": 5.78e-3, "D": 13.14e-6}}    # Table 3, p. 443
N_MILDNER_PER_M3 = 2.35e19 * 1e9                                          # p. 441, Al2O3 units
LAMS = (0.10, 0.12, 0.1798, 0.20, 0.30, 0.40, 0.60)                       # nm
E_meV = lambda lam: 0.0818042 / lam ** 2 * 1e3 / 100  # noqa: E731   E[meV] = 81.8042 / lam[A]^2


def sigma_mildner(lam, p, form):
    if form == "eq3":
        Sig = p["A"] * lam + p["C"] * (1 - np.exp(-(p["B"] / lam ** 2 + p["D"] / lam ** 4)))
    else:
        Sig = p["A"] * lam + p["C"] * (1 - lam ** 2 / p["B2"] * (1 - np.exp(-(p["B2"] / lam ** 2 + p["D"] / lam ** 4))))
    return Sig / N_MILDNER_PER_M3 * 1e28        # barn per Al2O3 unit


def main():
    models, _ = S.build_models("Al2O3")
    cfg = models["mp_dft"]
    _, summ = S.mp_entry("Al2O3", "R-3c")
    n_units_cm3 = float(summ.nsites / summ.volume) / 5 * 1e24       # Al2O3 units per cm^3 from the MP cell
    E = np.array([81.8042 / (10 * l) ** 2 for l in LAMS]) * 1e-3     # eV
    out = {"note": "literature inputs transcribed in notes/stage2_anchors.md section 4", "per_lambda": []}
    s293 = S.sigma(cfg, 293.0, E, True) * 5
    for l, e, sm in zip(LAMS, E, s293):
        out["per_lambda"].append({"lambda_nm": l, "E_meV": round(e * 1e3, 2), "model_b_per_unit": round(float(sm), 3),
                                  "mildner_eq3_b": round(float(sigma_mildner(l, HEMEX["eq3"], "eq3")), 3),
                                  "mildner_eq4_b": round(float(sigma_mildner(l, HEMEX["eq4"], "eq4")), 3),
                                  "ratio_model_eq4": round(float(sm / sigma_mildner(l, HEMEX["eq4"], "eq4")), 3)})
    win = [r for r in out["per_lambda"] if 0.12 <= r["lambda_nm"] <= 0.40]
    out["model_over_mildner_eq4_0.12_0.40nm"] = [min(r["ratio_model_eq4"] for r in win), max(r["ratio_model_eq4"] for r in win)]
    out["model_vs_mildner_eq4_rms_ln_0.10_0.60nm"] = round(float(np.sqrt(np.mean(
        [np.log(r["ratio_model_eq4"]) ** 2 for r in out["per_lambda"]]))), 3)
    # Tennant: T(0.12 nm, 101.6 mm) ratio 77 K / 300 K
    e12 = np.array([81.8042 / 1.2 ** 2 * 1e-3])
    T = {t: float(np.exp(-n_units_cm3 * S.sigma(cfg, t, e12, True)[0] * 5 * 1e-24 * 10.16)) for t in (77.0, 300.0)}
    out["tennant_0.12nm_101.6mm"] = {"T_300K": round(T[300.0], 4), "T_77K": round(T[77.0], 4),
                                     "model_increase_percent": round(100 * (T[77.0] / T[300.0] - 1), 1),
                                     "literature_increase_percent": 18}
    # Born: 90 mm filter, 0.12-0.24 nm, 293 K
    eb = np.array([81.8042 / (10 * l) ** 2 for l in (0.12, 0.18, 0.24)]) * 1e-3
    tb = np.exp(-n_units_cm3 * S.sigma(cfg, 293.0, eb, True) * 5 * 1e-24 * 9.0)
    out["born_90mm_0.12_0.24nm"] = {"model_T": [round(float(x), 3) for x in tb], "literature_T_about": 0.8}
    # transcribed literature constants quoted in the text (Mildner et al. 1993, pp. 442-443)
    out["mildner_fitted_A_b_per_nm"] = [2.537, 2.565]
    out["mildner_absorption_b_per_nm"] = 2.575
    out["model_number_density_units_per_cm3"] = n_units_cm3
    out["mildner_number_density_units_per_cm3"] = 2.35e19 * 1e3
    (F("stage2_lit_compare.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
