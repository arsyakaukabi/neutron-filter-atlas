"""Sapphire literature anchors of the atlas (step 08, atlas), computed with the atlas models. The parameters below are
transcribed from the cited sources.

1. Sapphire quality factor (Mildner & Lamaze 1998, J. Appl. Cryst. 31, 835, p. 839): sigma at the attenuation minimum
   divided by the free-atom (high-energy) cross section. Measured 0.076(1) (hemex, ambient). Nieman et al. 1980 as
   quoted there: 0.06 at room temperature, 0.05 at 80 K. Model: min over 0.3-10 A of the atlas off-Bragg mu(lambda, T)
   divided by n * sum x_i sigma_free,i (NCrystal free scattering + absorption at 1 eV, negligible).
2. Sapphire attenuation fits Sigma(lambda) = A lambda + C [1 - exp(-B/lambda^2 - D/lambda^4)] (lambda in nm, Sigma in 1/m):
   Mildner & Lamaze Table 1, mean unirradiated hemex (ambient): A 8.08, C 21.0, B 1.64e-3, D 12.3e-6;
   Table 3, Tennant (1988) data: 300 K A 8.88, C 20.6, B 0.57e-3, D 19.4e-6; 80 K A 6.24, C 19.2, B -0.15e-3,
   D 15.8e-6. Model/fit ratio of mu at 0.12, 0.2, 0.3, 0.4 nm, and T(101.6 mm) at 0.12 nm both temperatures.
3. Tennant cooling gain (Mildner 1993 p. 438): transmission of 0.12 nm neutrons through 101.6 mm sapphire increases by
   18 % from 300 to 77 K.
4. Born et al. 1987 (Mildner 1993 p. 439): 90 mm sapphire, T about 0.8 for 0.12-0.24 nm.
-> the atlas anchors JSON in data/derived
"""
import json
import sys
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd
from pymatgen.core import Composition

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D as _D, F  # noqa: E402   # D is also a literature fit parameter below

FITS = {"mildner_lamaze_hemex_ambient": (293.15, 8.08, 21.0, 1.64e-3, 12.3e-6),
        "tennant_300K": (293.15, 8.88, 20.6, 0.57e-3, 19.4e-6),
        "tennant_80K": (77.0, 6.24, 19.2, -0.15e-3, 15.8e-6)}


def sigma_fit(lam_nm, A, C, B, D):
    return A * lam_nm + C * (1 - np.exp(-B / lam_nm ** 2 - D / lam_nm ** 4))   # 1/m


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    row = at[(at.formula == "Al2O3") & (at.sg_symbol == "R-3c") & (at.tier == "A")].iloc[0]
    t1 = pd.read_parquet(F("stage1_triviality.parquet")).set_index("material_id").loc[row.material_id]
    cfg = nb.build_vdos_material(mq.load_dos(_D("mp") / "phonon_dos" / t1.phonon_method / f"{t1.phonon_id}.json.gz"),
                                 "anc_sapphire")
    n = float(row.numberdensity_A3)
    comp = Composition("Al2O3").fractional_composition
    s_free = sum(float(x) * NC.atomDB(str(e)).freeScatteringXS() for e, x in comp.items())
    out = {"material_id": row.material_id, "n_A3": n, "sigma_free_b_per_atom": s_free}
    lam = np.linspace(0.3, 10.0, 971)
    qf = {}
    for T in (293.15, 77.0):
        mu, _ = xe.mu_offbragg(cfg, T, lam)
        i = int(np.argmin(mu))
        qf[f"{T:g}K"] = {"quality_factor": float(mu[i] / (n * s_free)), "lambda_min_A": float(lam[i])}
    out["quality_factor"] = {"model": qf, "mildner_lamaze_1998_measured_ambient": 0.076,
                             "nieman_1980_room_temperature": 0.06, "nieman_1980_80K": 0.05}
    fits = {}
    for name, (T, A, C, B, D) in FITS.items():
        lams_nm = np.array([0.12, 0.2, 0.3, 0.4])
        mu_model, _ = xe.mu_offbragg(cfg, T, lams_nm * 10)                    # cm^-1
        mu_fit = sigma_fit(lams_nm, A, C, B, D) / 100.0                        # 1/m -> 1/cm
        fits[name] = {"T_model_K": T, "lambda_nm": lams_nm.tolist(), "mu_fit_cm": mu_fit.round(5).tolist(),
                      "mu_model_cm": mu_model.round(5).tolist(), "model_over_fit": (mu_model / mu_fit).round(3).tolist(),
                      "T_101.6mm_0.12nm_fit": float(np.exp(-mu_fit[0] * 10.16)),
                      "T_101.6mm_0.12nm_model": float(np.exp(-mu_model[0] * 10.16))}
    out["attenuation_fits"] = fits
    g_fit = fits["tennant_80K"]["T_101.6mm_0.12nm_fit"] / fits["tennant_300K"]["T_101.6mm_0.12nm_fit"] - 1
    g_mod = fits["tennant_80K"]["T_101.6mm_0.12nm_model"] / fits["tennant_300K"]["T_101.6mm_0.12nm_model"] - 1
    out["tennant_cooling_gain_0.12nm_101.6mm"] = {"literature_percent": 18, "tennant_fits_percent": round(100 * g_fit, 1),
                                                  "model_percent": round(100 * g_mod, 1)}
    mu, _ = xe.mu_offbragg(cfg, 293.15, np.linspace(1.2, 2.4, 25))
    out["born_90mm_0.12_0.24nm"] = {"literature_T_about": 0.8, "model_T_range": [float(np.exp(-mu.max() * 9.0)),
                                                                                  float(np.exp(-mu.min() * 9.0))]}
    (F("stage3_anchors.json")).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
