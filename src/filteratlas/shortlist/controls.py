"""Controls of the shortlist after transport, and the literature anchors.

Positive controls: sapphire and Si must stay good thermal filters (L0.8 in the top quartile of the Tier-A full-model
set at both temperatures), quartz intermediate (rank percentile 25-75), with their transport corrections and the
thermal-source transmission; Bi stays the known failing control (L0.8 overstated in the EXFOR benchmark). The Nieman
ratio of sapphire, Si and quartz tests the literature's effectiveness order. Anchors are tabulated from the anchor
file of the atlas step (literature values as quoted, with source and page).
-> derived/stage4_controls.json
"""
import json

import pandas as pd

from filteratlas.common.paths import F

REF = {"sapphire": ("Al2O3", "R-3c"), "Si": ("Si", "Fd-3m"), "quartz": ("SiO2", "P3_221"), "MgO": ("MgO", "Fm-3m")}


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    tr = json.loads(F("stage4_transport.json").read_text())
    anc = json.loads(F("stage3_anchors.json").read_text())
    unc = json.loads(F("stage2_uncertainty_model.json").read_text())
    out = {"positive": {}, "pass": {}}
    for tag, T in (("293K", 293.15), ("77K", 77.0)):
        a = at[(at["T"] == T) & (at.tier == "A") & at.vdos_source.ne("none")].copy()
        a["pct"] = a.L08.rank(pct=True) * 100           # 100 = most transparent
        for k, (f, sg) in REF.items():
            r = a[(a.formula == f) & (a.sg_symbol == sg)].iloc[0]
            lab = next((n for n in tr["at_L08"] if n.startswith(f"{f} {sg} ")), None)
            e = tr["at_L08"].get(lab, {}).get(tag, {})
            out["positive"].setdefault(k, {})[tag] = {
                "L08_cm": round(float(r.L08), 3), "L08_percentile_tierA": round(float(r.pct), 1), "n_tierA": len(a),
                "T_th_uncollided_openmc": e.get("T_th_uncollided"), "T_th_all_thermal_openmc": e.get("T_th_all_thermal"),
                "correction_A_fast_1E": e.get("A_fast_1E", {}).get("correction"),
                "correction_A_gamma_1MeV": e.get("A_gamma_1MeV", {}).get("correction")}
    pos = out["positive"]
    out["pass"]["sapphire_Si_top_quartile_both_T"] = all(pos[k][t]["L08_percentile_tierA"] >= 75
                                                         for k in ("sapphire", "Si") for t in ("293K", "77K"))
    out["pass"]["quartz_intermediate_both_T"] = all(25 <= pos["quartz"][t]["L08_percentile_tierA"] <= 75
                                                    for t in ("293K", "77K"))
    # the literature's own effectiveness metric is the Nieman ratio mu_epi(1 eV) / mu_th(1.8 A);
    # "sapphire significantly more effective than either quartz or silicon" holds if sapphire's ratio is the largest
    nie = {}
    for tag, T in (("293K", 293.15), ("77K", 77.0)):
        a = at[(at["T"] == T) & (at.tier == "A")]
        nie[tag] = {k: round(float(a[(a.formula == f) & (a.sg_symbol == sg)].nieman_ratio.iloc[0]), 2)
                    for k, (f, sg) in REF.items() if k in ("sapphire", "Si", "quartz")}
    out["nieman_ratio"] = nie
    out["pass"]["sapphire_largest_nieman_ratio_both_T"] = all(v["sapphire"] > max(v["Si"], v["quartz"]) for v in nie.values())
    bi = next((n for n in tr["at_L08"] if n.startswith("Bi R-3m ")), None)
    out["Bi_known_failing"] = {"L08_overstatement_percent_stage2": unc["heavy_soft_metal_L08_overstatement_percent"][1],
                               "transport": tr["at_L08"].get(bi)}
    q, af = anc["quality_factor"], anc["attenuation_fits"]
    out["anchors"] = [
        {"anchor": "Nieman et al. 1980 quality factor, sapphire, room temperature (as quoted by Mildner & Lamaze 1998 p. 839)",
         "literature": q["nieman_1980_room_temperature"], "model": round(q["model"]["293.15K"]["quality_factor"], 4),
         "other_measurement": {"Mildner & Lamaze 1998 hemex, ambient": q["mildner_lamaze_1998_measured_ambient"]}},
        {"anchor": "Nieman et al. 1980 quality factor, sapphire, 80 K (as quoted by Mildner & Lamaze 1998 p. 839)",
         "literature": q["nieman_1980_80K"], "model": round(q["model"]["77K"]["quality_factor"], 4)},
        *[{"anchor": f"attenuation vs wavelength, {k} (Mildner et al. 1993 / Mildner & Lamaze 1998 fits)",
           "lambda_nm": v["lambda_nm"], "literature_mu_cm": v["mu_fit_cm"], "model_mu_cm": v["mu_model_cm"],
           "model_over_fit": v["model_over_fit"], "T_101.6mm_0.12nm": {"fit": round(v["T_101.6mm_0.12nm_fit"], 4),
                                                                      "model": round(v["T_101.6mm_0.12nm_model"], 4)}}
          for k, v in af.items()],
        {"anchor": "Tennant cooling gain at 0.12 nm through 101.6 mm (percent)",
         "literature": anc["tennant_cooling_gain_0.12nm_101.6mm"]["literature_percent"],
         "tennant_fits": anc["tennant_cooling_gain_0.12nm_101.6mm"]["tennant_fits_percent"],
         "model": anc["tennant_cooling_gain_0.12nm_101.6mm"]["model_percent"]},
        {"anchor": "Born et al. 1987 90 mm sapphire filter transmission, 0.12-0.24 nm (as quoted by Mildner et al. 1993)",
         "literature": anc["born_90mm_0.12_0.24nm"]["literature_T_about"],
         "model_range": [round(x, 3) for x in anc["born_90mm_0.12_0.24nm"]["model_T_range"]]},
    ]
    out["all_pass"] = all(out["pass"].values())
    F("stage4_controls.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps({"pass": out["pass"], "positive": out["positive"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
