"""Residual structure of the EXFOR benchmark and the model-uncertainty table used by the atlas (step 06, uncertainty).

1. Residual structure. For every qualifying single-crystal dataset (single-crystal forms incl. pyrolytic graphite;
   ZrO2 excluded, phase unknown) and the models mp_dft / SevenNet / curated: binned ln(sigma_model/sigma_data), q20
   envelope, 1-100 meV. Classes fixed before computing: agrees | offset | curved (dip / hump) | rising | falling |
   irregular. Envelope dependence: r_d(q10) - r_d(q50).
2. Uncertainty model. 68th percentile of |ln(sigma_model/sigma_data)| over 1-100 meV bins, single-crystal datasets,
   each crystal weighted equally (bins weighted 1/n_bins of their crystal), bootstrap over crystals (seed fixed):
   sigma_model(DFT) over the crystals with MP DFT; sigma_model(SevenNet) over all evaluable crystals (MP DFT where it
   exists, otherwise SevenNet = the atlas rule) and over the SevenNet-only crystals; heavy_soft_metal (Pb, Bi).
   Variant without MgO 13761004 (polycrystal-behaving sample, see filteratlas.validation.mgo_audit). Plus the
   MLIP-vs-DFT interval from the MLIP benchmark. Since L0.8 ~ 1/sigma in the band, ln-sigma errors carry over
   to ln L0.8 one to one.
-> the residuals JSON and the uncertainty-model JSON in data/derived
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.paths import D, F  # noqa: E402
from filteratlas.common.constants import BAND  # noqa: E402
from filteratlas.xs.evaluate import maxwellian_flux  # noqa: E402
from filteratlas.validation.exfor_benchmark import CHOSEN_MLIP, E_KILL, E_LO, KILL_RMS, SC_FORMS  # noqa: E402

ATLAS = f"mlip_{CHOSEN_MLIP}"
MODELS = ("mp_dft", ATLAS, "curated")
from filteratlas.common.config import cfg  # noqa: E402

HEAVY = tuple(cfg("benchmark")["heavy_soft_metals"])
RNG = np.random.default_rng(cfg("benchmark")["residuals_seed"])


def classify(E, lr):
    x = np.log10(E / 1e-2)
    r = float(np.sqrt(np.mean(lr ** 2)))
    mean, sd = float(lr.mean()), float(lr.std())
    p1 = np.polyfit(x, lr, 1)
    res1 = float(np.sqrt(np.mean((lr - np.polyval(p1, x)) ** 2)))
    p2 = np.polyfit(x, lr, 2) if len(x) >= 4 else np.array([0.0, *p1])
    res2 = float(np.sqrt(np.mean((lr - np.polyval(p2, x)) ** 2)))
    if r <= KILL_RMS:
        cls = "agrees"
    elif abs(mean) >= 0.10 and sd <= 0.10:
        cls = "offset"
    elif abs(p2[0]) >= 0.15 and res1 > 0 and res2 <= 0.7 * res1:
        cls = "curved_dip" if p2[0] > 0 else "curved_hump"
    elif p1[0] >= 0.10:
        cls = "rising"
    elif p1[0] <= -0.10:
        cls = "falling"
    else:
        cls = "irregular"
    return {"class": cls, "r_d": round(r, 4), "mean": round(mean, 4), "sd": round(sd, 4),
            "slope_per_decade": round(float(p1[0]), 4), "curvature": round(float(p2[0]), 4), "n_bins": len(x)}


def p68_weighted(a, w):
    o = np.argsort(a)
    c = np.cumsum(w[o]) / w.sum()
    return float(a[o][np.searchsorted(c, 0.68)])


def pooled(bins, crystals, B=2000):
    g = bins[bins.formula.isin(crystals)]
    if g.empty:
        return None
    a = np.abs(g.lr.values)
    w = 1.0 / g.groupby("formula").lr.transform("size").values
    est = p68_weighted(a, w)
    cs = sorted(g.formula.unique())
    per = {c: np.abs(g[g.formula == c].lr.values) for c in cs}
    bs = []
    for _ in range(B):
        pick = RNG.choice(cs, len(cs))   # resample crystals; each drawn crystal keeps weight 1 over its bins
        bs.append(p68_weighted(np.concatenate([per[c] for c in pick]),
                               np.concatenate([np.full(len(per[c]), 1.0 / len(per[c])) for c in pick])))
    return {"p68_abs_ln": round(est, 4), "ci95": [round(float(np.percentile(bs, 2.5)), 4),
                                                   round(float(np.percentile(bs, 97.5)), 4)],
            "crystals": cs, "n_bins": len(g)}


def main():
    full = json.loads((F("stage2_exfor.json")).read_text())
    qual = {(a["subentry"], a["T"]) for a in full["audit"] if a["qualifies"]}
    bn = pd.read_parquet(F("stage2_exfor_binned.parquet"))
    bn = bn[np.array([(s, t) in qual for s, t in zip(bn.subentry, bn["T"])]) & bn.form.isin(SC_FORMS)
            & (bn.formula != "ZrO2") & (bn.E_c >= E_LO) & (bn.E_c <= E_KILL)]
    bench = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    bench = bench[bench.qualifies & bench.form.isin(SC_FORMS)]
    out = {"note": "classes fixed in NOTEBOOK 2026-10-22 before computing", "datasets": [], "class_counts": {}}
    for (f, sub, T, m), g in bn[bn.model.isin(MODELS)].groupby(["formula", "subentry", "T", "model"]):
        g = g.sort_values("E_c")
        rec = {"formula": f, "subentry": int(sub), "T": float(T), "model": m, **classify(g.E_c.values, g.lr.values)}
        e = bench[(bench.subentry == sub) & (bench["T"] == T) & (bench.model == m)].set_index("envelope_q").r_d
        rec["envelope_dependence_q10_minus_q50"] = round(float(e.get(10) - e.get(50)), 4) if 10 in e and 50 in e else None
        out["datasets"].append(rec)
    d = pd.DataFrame(out["datasets"])
    out["class_counts"] = {m: d[d.model == m]["class"].value_counts().to_dict() for m in MODELS}
    (F("stage2_residuals.json")).write_text(json.dumps(out, indent=1))

    dft_cr = sorted(f for f, c in full["crystals"].items() if "mp_dft" in c["per_model"] and f != "ZrO2")
    atlas_bins = pd.concat([bn[(bn.model == "mp_dft") & bn.formula.isin(dft_cr)],
                            bn[(bn.model == ATLAS) & ~bn.formula.isin(dft_cr)]])
    no_004 = lambda x: x[x.subentry != 13761004]  # noqa: E731
    mlip = json.loads((F("stage2_mlip.json")).read_text())["models"][CHOSEN_MLIP]
    round1 = {"note": "SUPERSEDED (referee round 1, M4): per-bin p68 over 1-100 meV, incl. out-of-scope data",
              "sigma_model_DFT": pooled(bn[bn.model == "mp_dft"], dft_cr),
              "sigma_model_DFT_without_MgO_13761004": pooled(no_004(bn[bn.model == "mp_dft"]), dft_cr),
              "sigma_model_SevenNet_all_crystals": pooled(bn[bn.model == ATLAS], sorted(bn[bn.model == ATLAS].formula.unique())),
              "sigma_model_heavy_soft_metal": pooled(bn[bn.model == ATLAS], list(HEAVY)),
              "sigma_model_curated": pooled(bn[bn.model == "curated"], sorted(bn.formula.unique()))}

    # --- revised model (replaces the per-bin p68 above) ----------------------------------------------------------
    # scope: single-crystal datasets that are off-Bragg measurements: no pyrolytic graphite, no dataset
    # reclassified as polycrystal-behaving by the general rule (robustness analysis of step 05), no ZrO2
    rob = json.loads((F("stage2_d11_robustness.json")).read_text())
    recl = {(int(a), float(b)) for a, b in rob["reclassified"]}
    bf = pd.read_parquet(F("stage2_exfor_binned.parquet"))
    bf = bf[np.array([(s_, t_) in qual and (int(s_), float(t_)) not in recl for s_, t_ in zip(bf.subentry, bf["T"])])
            & bf.form.isin(SC_FORMS - {"pyrolytic_graphite"}) & (bf.formula != "ZrO2")]
    lam = np.sqrt(81.8042 / (bf.E_c * 1e3))                                   # Angstrom
    bf = bf.assign(lam=lam, w=maxwellian_flux(lam) * lam * np.log(10) / 20 / 2)  # flux x d(lambda) of a log-E bin
    band = bf[(bf.lam >= BAND[0]) & (bf.lam <= BAND[1])]
    per_ds = (band.groupby(["formula", "subentry", "T", "model"])
              .apply(lambda g: pd.Series({"band_mean_ln": float(np.sum(g.w * g.lr) / np.sum(g.w)), "n_bins": len(g)}),
                     include_groups=False).reset_index())
    per_cr = per_ds.groupby(["formula", "model"]).band_mean_ln.mean().unstack()

    def summary(vals):
        v = np.asarray([x for x in vals if np.isfinite(x)], float)
        return {"n_crystals": len(v), "rms": round(float(np.sqrt(np.mean(v ** 2))), 4),
                "common_mode_mean": round(float(v.mean()), 4), "compound_specific_sd": round(float(v.std(ddof=1)), 4),
                "p68_abs": round(float(np.percentile(np.abs(v), 68)), 4)}

    dft_in = [f for f in dft_cr if f in per_cr.index and np.isfinite(per_cr.loc[f, "mp_dft"])]
    sc_rc = pd.DataFrame(rob["datasets"])
    sc_rc = sc_rc[np.array([(int(a), float(b)) not in recl for a, b in zip(sc_rc.subentry, sc_rc["T"])]) & (sc_rc.form != "pyrolytic_graphite")]
    rc_in = {f: float(np.sqrt(np.mean(g.r_atlas_q20 ** 2))) for f, g in sc_rc.groupby("formula") if f in dft_in}
    unc = {"note": ("Revised after referee round 1 (M4). Scope: off-Bragg single-crystal datasets (no pyrolytic graphite, "
                    "no reclassified polycrystal-behaving sets, no ZrO2). Primary quantity: per crystal, the mean over its "
                    "datasets of the flux-weighted (300 K Maxwellian, 0.5-3.0 A) mean ln(sigma_model/sigma_data); "
                    "ln L0.8 errors equal minus these to first order. PREREGISTRATION s7 'RMS log-ratio' reported as "
                    "the primary statistic; p68 alongside."),
           "preregistered_rms_all_forms_DFT": round(float(np.sqrt(np.mean(
               [full["crystals"][f]["per_model"]["mp_dft"]["R_c"] ** 2 for f in dft_cr]))), 4),
           "preregistered_rms_in_scope_DFT_1_100meV": round(float(np.sqrt(np.mean([v ** 2 for v in rc_in.values()]))), 4),
           "in_scope_R_c_DFT_1_100meV": {f: round(v, 4) for f, v in rc_in.items()},
           "band_per_crystal_DFT": {f: round(float(per_cr.loc[f, "mp_dft"]), 4) for f in dft_in},
           "band_DFT": summary([per_cr.loc[f, "mp_dft"] for f in dft_in]),
           "band_per_crystal_SevenNet": {f: round(float(per_cr.loc[f, ATLAS]), 4) for f in per_cr.index
                                         if ATLAS in per_cr and np.isfinite(per_cr.loc[f, ATLAS])},
           "band_SevenNet_all": summary(per_cr[ATLAS].values) if ATLAS in per_cr else None,
           "band_curated_all": summary(per_cr["curated"].values),
           "heavy_soft_metal_band_values_no_ci": {f: round(float(per_cr.loc[f, ATLAS]), 4) for f in HEAVY
                                                  if f in per_cr.index},
           "heavy_soft_metal_abs_max": round(float(max(abs(per_cr.loc[f, ATLAS]) for f in HEAVY if f in per_cr.index)), 4),
           "heavy_soft_metal_L08_overstatement_percent": [round(100 * (np.exp(-per_cr.loc[f, ATLAS]) - 1), 1)
                                                          for f in HEAVY if f in per_cr.index],
           "mlip_vs_dft_dlnL08": {t: {"p68_abs_dlog": mlip[t]["p68_abs_dlog"], "ci95": mlip[t]["p68_ci95"],
                                      "median_signed_dlog": mlip[t]["median_dlog"]} for t in ("293K", "77K")},
           "per_dataset_band": per_ds.round(4).to_dict("records"),
           "round1_superseded": round1}
    # the per-crystal-mean construction is one of several defensible ones; report all, each with a bootstrap CI over
    # crystals (seed fixed), so that the choice for the atlas can be made by principle, not by size
    def boot(groups, stat, B=4000):
        keys = sorted(groups)
        est = stat([groups[k] for k in keys])
        bs = [stat([groups[k] for k in RNG.choice(keys, len(keys))]) for _ in range(B)]
        return {"value": round(float(est), 4), "ci95": [round(float(np.percentile(bs, 2.5)), 4),
                                                         round(float(np.percentile(bs, 97.5)), 4)], "n_crystals": len(keys)}
    dsd = per_ds[(per_ds.model == "mp_dft") & per_ds.formula.isin(dft_in)]
    bands = band[(band.model == "mp_dft") & band.formula.isin(dft_in)]
    g_ds = {f: g.band_mean_ln.values for f, g in dsd.groupby("formula")}
    g_bin = {f: (g.lr.values, g.w.values) for f, g in bands.groupby("formula")}
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))  # noqa: E731
    unc["band_DFT_constructions"] = {
        "A_per_crystal_mean_then_rms": {**boot(g_ds, lambda L: rms([np.mean(x) for x in L])),
                                        "definition": "mean of dataset band means per crystal; RMS over crystals"},
        "B_per_crystal_rms_then_rms": {**boot(g_ds, lambda L: rms([rms(x) for x in L])),
                                       "definition": "RMS of dataset band means per crystal (as R_c is built); RMS over crystals"},
        "C_per_dataset": {**boot(g_ds, lambda L: rms(np.concatenate(L))),
                          "definition": "RMS over all in-scope datasets of their band means (datasets weighted equally)"},
        "D_within_band_bins": {**boot(g_bin, lambda L: float(np.sqrt(np.sum(np.concatenate([x[1] * x[0] ** 2 for x in L]))
                                                                     / np.sum(np.concatenate([x[1] for x in L]))))),
                               "definition": "flux-weighted RMS of per-bin ln ratios within the band, pooled"}}
    unc["band_DFT_constructions"]["note"] = (
        "A is smallest partly because opposite-signed datasets of one crystal cancel (e.g. Si). The Stage 3 choice is to "
        "be fixed in PREREG_STAGE3 by principle (what the atlas value represents), not by size.")
    sd = unc["band_DFT"]["rms"]
    unc["proposed_for_prereg_stage3"] = {
        "DFT_sigma_ln": sd,
        "MLIP_sigma_ln_293K": round(float(np.hypot(sd, mlip["293K"]["p68_abs_dlog"])), 4),
        "MLIP_bias_dlnL08_293K": mlip["293K"]["median_dlog"],
        "common_mode_fraction_of_variance": round(float(unc["band_DFT"]["common_mode_mean"] ** 2 / sd ** 2), 3),
        "rule": "DFT: sigma = RMS of in-scope band-weighted crystal means (PREREGISTRATION s7 statistic). MLIP: that "
                "(+) the MLIP-vs-DFT p68 (s7), with the MLIP median signed bias reported, not subtracted. "
                "For the Stage 4 Monte Carlo split sigma^2 into a common-mode part (mean^2) shared by all compounds and "
                "a compound-specific remainder. heavy_soft_metal: definition and width to be fixed in PREREG_STAGE3."}
    (F("stage2_uncertainty_model.json")).write_text(json.dumps(unc, indent=1))
    print(json.dumps(out["class_counts"], indent=1))
    print(json.dumps({k: v for k, v in unc.items() if k not in ("note", "per_dataset_band", "round1_superseded")}, indent=1))


if __name__ == "__main__":
    main()
