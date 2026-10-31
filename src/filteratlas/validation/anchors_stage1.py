"""Literature cross-checks of the MP-DFT model (step 05, validation). The parameters below are transcribed from the
cited sources.

1. Freund-type model (as written by Barker, Mildner, Rodriguez & Thiyagarajan 2008, J. Appl. Cryst. 41:1003,
   Eqs. 1-5, parameters from their Table 2) vs our full model, sapphire and Si, 293 K and 77 K.
   Self-check: B0 and B_T for MgF2 must reproduce Barker's printed 2.29e-3 nm^2 and 2.62/0.24e-3 nm^2.
2. Measured-fit curves (Barker 2008 Table 3: Si 77 K [Freund et al. 1985], Al2O3 300 K [Mildner & Lamaze 1998],
   MgO 77 K [Thiyagarajan et al. 1998], MgF2 77/293 K [Barker 2008]) vs our full model, 1-3 A.
3. Holmryd & Connor 1969 abstract: 15 cm MgO at 83 K transmits >70 % at 2 A; SiO2 nearly as good;
   PbF2 better than Bi at 1-3.5 A; Bi ~50 % at >4 A (15 cm).
4. Al-Qasir & Qteish 2017 abstract: BeF2 better filter than MgF2 and MgO; BeO transmits low-energy neutrons better.
5. Nieman et al. 1980 abstract: sapphire more effective (mu_epi/mu_th) than quartz or Si at 300 K and 80 K.
Output: the anchors JSON in data/derived.
"""
import json
import sys
from math import factorial
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import bernoulli

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D as _D, F  # noqa: E402   # D is also a literature fit parameter below

H, KB, MN, MP_ = 6.62607015e-34, 1.380649e-23, 1.67492749804e-27, 1.67262192e-27
BN = bernoulli(22)


def R_series(x):
    if x > 6:
        return 3.3 * x ** -3.5
    return sum(BN[n] * x ** (n - 1) / (factorial(n) * (n + 2.5)) for n in range(23))


def freund_sigma(lam_nm, T, N, M, A_abs, theta, sigma_fa, sigma_b):
    """Barker 2008 Eqs (1)-(5); returns macroscopic Sigma in m^-1."""
    lam_m = lam_nm * 1e-9
    s_abs = A_abs * lam_nm
    s_sph = (3 * N * sigma_b / M) * np.sqrt(2 * MN * KB * theta / H ** 2) * lam_m * R_series(theta / T)
    B0 = 3 * H ** 2 / (2 * KB * theta * M * MP_) * 1e18                      # nm^2
    integ = quad(lambda z: z / np.expm1(z), 0, theta / T)[0]
    BT = 4 * B0 * T ** 2 / theta ** 2 * integ
    s2 = 0.5 - 1 / (3 * M)
    s_mph = sigma_fa * (1 - np.exp(-(B0 + BT) * s2 / lam_nm ** 2))
    return s_abs + s_sph + s_mph, B0, BT


def fit_curve(lam_nm, A, sfa, B, D):
    return A * lam_nm + sfa * (1 - np.exp(-(B / lam_nm ** 2 + D / lam_nm ** 4)))


def sigma_bound_per_fu(formula):
    from pymatgen.core import Composition
    return sum(n * NC.atomDB(el.symbol).scatteringXS() for el, n in Composition(formula).items())


def model_mu(cfg, T, lam_A):
    return xe.mu_offbragg(cfg, T, lam_A)[0] * 100.0   # cm^-1 -> m^-1


def main():
    res = {}
    # --- self-check of the transcription: Barker's MgF2 B0 and B_T ---------------------------------------------
    _, b0, bt293 = freund_sigma(0.2, 293.0, 3.071, 20.405, 1.40, 610.0, 32.8, 11.75)
    _, _, bt77 = freund_sigma(0.2, 77.0, 3.071, 20.405, 1.40, 610.0, 32.8, 11.75)
    res["selfcheck_MgF2_B0_BT_nm2"] = {"B0": round(b0, 5), "BT_293": round(bt293, 5), "BT_77": round(bt77, 5),
                                       "barker_printed": {"B0": 2.29e-3, "BT_293": 2.62e-3, "BT_77": 0.24e-3}}

    # --- our models (MP DFT VDOS, lowest-e_hull polymorph with the classical structure) -----------------------
    s1 = pd.read_parquet(F("stage1_triviality.parquet"))
    PICK = {"Al2O3": "R-3c", "Si": "Fd-3m", "MgO": "Fm-3m", "MgF2": "P4_2/mnm", "SiO2": None, "PbF2": None,
            "BeF2": None, "BeO": None}
    cfgs, picked = {}, {}
    cand = pd.read_parquet(F("stage0_candidates.parquet"))
    for f, sg in PICK.items():
        rows = s1[(s1.formula == f) & (s1.status == "ok") & s1.stable.fillna(False)]
        if sg:
            rows = rows[rows.sg_symbol == sg]
        if f == "SiO2":
            # alpha-quartz: the P3_221 entry is dynamically clean on a dense mesh; P3_121 (mp-aaaaakjg) is not
            # (force-constant validation of the MP phonons). P3_121 is reported for information only.
            rows = rows[rows.sg_symbol == "P3_221"]
        rows = rows.merge(cand[["material_id", "e_hull"]], on="material_id").sort_values("e_hull")
        if rows.empty:
            picked[f] = None
            continue
        r = rows.iloc[0]
        dos = mq.load_dos(_D("mp") / "phonon_dos" / r.phonon_method / f"{r.phonon_id}.json.gz")
        cfgs[f] = nb.build_vdos_material(dos, f"anchor_{f}")
        picked[f] = {"material_id": r.material_id, "sg": r.sg_symbol, "e_hull": float(r.e_hull)}
    cfgs["Bi"] = "stdlib::Bi_sg166.ncmat"
    picked["Bi"] = {"source": "NCrystal stdlib (no MP phonons for Bi)"}
    res["models_used"] = picked

    # --- 1. Freund (Barker form, Table 2 parameters) vs full model --------------------------------------------
    lam_nm = np.linspace(0.05, 0.30, 26)       # 0.5-3 A thermal band
    T2 = {"Al2O3": (2.345, 17.351, 6.05, 1032.0, 32.99), "Si": (4.9939, 28.086, 4.76, 543.0, 10.10)}
    fr = {}
    for f, (N, M, A, th, sfa) in T2.items():
        sb = sigma_bound_per_fu(f)
        for T in (293.15, 77.0):
            fre = np.array([freund_sigma(l, T, N, M, A, th, sfa, sb)[0] for l in lam_nm])
            ours = model_mu(cfgs[f], T, lam_nm * 10)
            lr = np.log(ours / fre)
            fr[f"{f}_{int(T)}K"] = {"median_ratio_ours_over_freund": round(float(np.exp(np.median(lr))), 3),
                                    "rms_log_ratio": round(float(np.sqrt(np.mean(lr ** 2))), 3),
                                    "ratio_at_1A": round(float(ours[np.argmin(abs(lam_nm - 0.1))] /
                                                               fre[np.argmin(abs(lam_nm - 0.1))]), 3),
                                    "ratio_at_2.5A": round(float(ours[np.argmin(abs(lam_nm - 0.25))] /
                                                                 fre[np.argmin(abs(lam_nm - 0.25))]), 3)}
    res["freund_barker_form"] = fr

    # --- 2. measured-fit curves (Barker Table 3) vs full model, 1-3 A ------------------------------------------
    lam2 = np.linspace(0.10, 0.30, 21)
    T3 = [("Si", 77.0, 4.83, 10.10, 2.11e-3, 0.0, "Freund et al. 1985"),
          ("Al2O3", 293.15, 8.11, 21.2, 1.61e-3, 1.29e-5, "Mildner & Lamaze 1998 (300 K)"),
          ("MgO", 77.0, 2.62, 38.4, 1.3e-3, 0.0, "Thiyagarajan et al. 1998"),
          ("MgF2", 77.0, 1.49, 32.8, 0.35e-3, 2.85e-5, "Barker et al. 2008"),
          ("MgF2", 293.15, 6.32, 32.8, 3.68e-3, 0.72e-5, "Barker et al. 2008")]
    mf = {}
    for f, T, A, sfa, B, D, src in T3:
        meas = fit_curve(lam2, A, sfa, B, D)
        ours = model_mu(cfgs[f], T, lam2 * 10)
        lr = np.log(ours / meas)
        mf[f"{f}_{int(T)}K"] = {"source": src, "median_ratio_ours_over_measfit": round(float(np.exp(np.median(lr))), 3),
                                "rms_log_ratio": round(float(np.sqrt(np.mean(lr ** 2))), 3)}
    res["measured_fit_curves_barker_table3"] = mf

    # --- 3. Holmryd & Connor 1969 (abstract) -------------------------------------------------------------------
    L = 15.0
    hc = {}
    # informational: the dynamically unstable alpha-quartz P3_121 entry (not used for verdicts)
    _r121 = s1[(s1.formula == "SiO2") & (s1.sg_symbol == "P3_121") & (s1.status == "ok")]
    if len(_r121):
        _d = mq.load_dos(_D("mp") / "phonon_dos" / _r121.iloc[0].phonon_method / f"{_r121.iloc[0].phonon_id}.json.gz")
        cfgs["SiO2_P3_121_unstable"] = nb.build_vdos_material(_d, "anchor_SiO2_P3_121")
    for f in ("MgO", "SiO2", "SiO2_P3_121_unstable", "PbF2", "Bi"):
        if f not in cfgs:
            hc[f] = "not in sample"
            continue
        mu2 = xe.mu_offbragg(cfgs[f], 83.0, [2.0])[0][0]
        band = xe.mu_offbragg(cfgs[f], 83.0, np.linspace(1.0, 3.5, 26))[0]
        long = xe.mu_offbragg(cfgs[f], 83.0, np.linspace(4.0, 6.0, 11))[0]
        hc[f] = {"T_2A_83K_15cm": round(float(np.exp(-mu2 * L)), 3),
                 "mean_T_1-3.5A_83K_15cm": round(float(np.mean(np.exp(-band * L))), 3),
                 "mean_T_4-6A_83K_15cm": round(float(np.mean(np.exp(-long * L))), 3)}
    hc["literature"] = "MgO >70 % at 2 A (15 cm, 83 K); SiO2 nearly as transparent; PbF2 > any Bi at 1-3.5 A; Bi ~50 % at >4 A"
    res["holmryd_connor_1969"] = hc

    # --- 4. Al-Qasir & Qteish 2017 (abstract) ------------------------------------------------------------------
    aq = {}
    for f in ("BeF2", "BeO", "MgF2", "MgO"):
        if f not in cfgs:
            aq[f] = "not in sample"
            continue
        r = xe.evaluate(cfgs[f], (293.15, 77.0))
        lowE = {f"mu_4A_{int(T)}K": round(float(xe.mu_offbragg(cfgs[f], T, [4.0])[0][0]), 4) for T in (293.15, 77.0)}
        aq[f] = {"L08_293K": round(r["L08_293K"], 2), "L08_77K": round(r["L08_77K"], 2), **lowE}
    aq["literature"] = "BeF2 better filter than MgF2 and MgO; BeO transmits low-energy neutrons better than BeF2"
    res["al_qasir_2017"] = aq

    # --- 5. Nieman et al. 1980 (abstract): mu_epi/mu_th, sapphire best ------------------------------------------
    nm = {}
    for f in ("Al2O3", "SiO2", "Si"):
        info = NC.createInfo(cfgs[f])
        sfree = sum(x * ad.scatteringXS() * (ad.averageMassAMU() / (ad.averageMassAMU() + 1.00866)) ** 2
                    for x, ad in info.composition)
        mu_epi = info.numberdensity * sfree           # free-atom limit, cm^-1
        for T in (293.15, 80.0):
            mu_th = xe.mu_offbragg(cfgs[f], T, [1.0])[0][0]
            nm[f"{f}_{int(T)}K"] = {"mu_epi_over_mu_th_1A": round(float(mu_epi / mu_th), 2)}
    nm["literature"] = "sapphire significantly more effective than quartz or Si (300 K and 80 K)"
    res["nieman_1980"] = nm

    # --- verdict table. The criteria were fixed before the model with the 2 meV floor on the MP vibrational DOS was
    # run. Caveat: values for two anchors (SiO2 T = 0.743 for the P3_121 entry; MgF2 0.901 / 0.997) were already
    # known when the criteria were written, so for those two the criterion was not blind; the SiO2 threshold is
    # therefore labelled as chosen after the first comparison. Criteria:
    #  Freund / measured fits: RMS log ratio <= 0.15 -> pass, 0.15-0.30 -> marginal, > 0.30 -> fail.
    #  H&C MgO: T(2 A) >= 0.70; H&C SiO2 'nearly as transparent': T_SiO2/T_MgO >= 0.85; H&C PbF2 > Bi (1-3.5 A);
    #  H&C Bi > 4 A '~50 %': 0.35 <= T <= 0.65; Nieman: sapphire ratio highest at both T;
    #  Nieman via Habib 2007 (secondary, attribution uncertain): quartz vs Si ratios within 20 % at 80 K;
    #  Al-Qasir: each stated ordering checked separately in each regime.
    def grade(rms):
        return "pass" if rms <= 0.15 else ("marginal" if rms <= 0.30 else "fail")


    V = []
    for k, v in fr.items():
        V.append({"anchor": f"Freund (Barker 2008 form) {k}", "value": v["median_ratio_ours_over_freund"],
                  "metric": "median ours/Freund", "rms_log": v["rms_log_ratio"], "verdict": grade(v["rms_log_ratio"])})
    for k, v in mf.items():
        V.append({"anchor": f"Measured fit {k} ({v['source']})", "value": v["median_ratio_ours_over_measfit"],
                  "metric": "median ours/fit", "rms_log": v["rms_log_ratio"], "verdict": grade(v["rms_log_ratio"])})
    if isinstance(hc.get("MgO"), dict):
        V.append({"anchor": "Holmryd&Connor MgO T(2A,83K,15cm) > 0.70", "value": hc["MgO"]["T_2A_83K_15cm"],
                  "verdict": "pass" if hc["MgO"]["T_2A_83K_15cm"] >= 0.70 else "fail"})
        if isinstance(hc.get("SiO2"), dict):
            r = hc["SiO2"]["T_2A_83K_15cm"] / hc["MgO"]["T_2A_83K_15cm"]
            V.append({"anchor": "Holmryd&Connor SiO2 'nearly as transparent' as MgO (ratio >= 0.85; post-hoc threshold)",
                      "value": round(r, 3), "verdict": "pass" if r >= 0.85 else "fail",
                      "note": "threshold not blind: value disclosed in referee round 1; dynamically clean P3_221 entry"})
    if isinstance(hc.get("PbF2"), dict) and isinstance(hc.get("Bi"), dict):
        V.append({"anchor": "Holmryd&Connor PbF2 > Bi at 1-3.5 A", "value": [hc["PbF2"]["mean_T_1-3.5A_83K_15cm"],
                                                                              hc["Bi"]["mean_T_1-3.5A_83K_15cm"]],
                  "verdict": "pass" if hc["PbF2"]["mean_T_1-3.5A_83K_15cm"] > hc["Bi"]["mean_T_1-3.5A_83K_15cm"] else "fail"})
        b = hc["Bi"]["mean_T_4-6A_83K_15cm"]
        V.append({"anchor": "Holmryd&Connor Bi > 4 A ~0.5 (0.35-0.65)", "value": b,
                  "verdict": "pass" if 0.35 <= b <= 0.65 else "fail"})
    for T in ("293K", "80K"):
        s, q, si = (nm[f"{f}_{T}"]["mu_epi_over_mu_th_1A"] for f in ("Al2O3", "SiO2", "Si"))
        V.append({"anchor": f"Nieman: sapphire most effective at {T}", "value": [s, q, si],
                  "verdict": "pass" if s > q and s > si else "fail"})
    q, si = nm["SiO2_80K"]["mu_epi_over_mu_th_1A"], nm["Si_80K"]["mu_epi_over_mu_th_1A"]
    V.append({"anchor": "Nieman via Habib 2007 (secondary): quartz ~ Si within 20 % at 80 K", "value": round(q / si, 3),
              "verdict": "pass" if abs(np.log(q / si)) <= np.log(1.2) else "fail", "note": "attribution uncertain"})
    if all(isinstance(aq.get(f), dict) for f in ("BeF2", "MgF2", "MgO", "BeO")):
        for key, better in (("L08_293K", "high"), ("L08_77K", "high"), ("mu_4A_293K", "low"), ("mu_4A_77K", "low")):
            bf = aq["BeF2"][key]
            others = [aq["MgF2"][key], aq["MgO"][key]]
            ok = all(bf > o for o in others) if better == "high" else all(bf < o for o in others)
            V.append({"anchor": f"Al-Qasir: BeF2 better than MgF2 and MgO ({key})", "value": [bf] + others,
                      "verdict": "pass" if ok else "fail"})
        for key in ("mu_4A_293K", "mu_4A_77K"):
            V.append({"anchor": f"Al-Qasir: BeO transmits low-energy neutrons better than BeF2 ({key})",
                      "value": [aq["BeO"][key], aq["BeF2"][key]],
                      "verdict": "pass" if aq["BeO"][key] < aq["BeF2"][key] else "fail"})
    res["verdicts"] = V
    res["verdict_counts"] = {k: sum(v["verdict"] == k for v in V) for k in ("pass", "marginal", "fail")}
    res["mgo_note"] = ("MgO 77 K vs Thiyagarajan et al. 1998 fit: unexplained. Barker et al.'s impurity remark concerns "
                       "only the long-wavelength linear term A, which changes the fit by at most ~17 % over 1-3 A "
                       "(referee estimate); the residual sits in the multiphonon/Debye-Waller regime.")
    (F("stage1_anchors.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
