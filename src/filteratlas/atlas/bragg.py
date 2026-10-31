"""Bragg cut-off section for Tier A:
lambda_c = NCrystal braggthreshold of the polycrystal built from the atlas VDOS structure, and the cold pass-band
transmission of a 5 cm polycrystal at 77 K for a 40 K Maxwellian over lambda_c < lambda <= 10 A, where Bragg
scattering vanishes so the off-Bragg mu applies. The band is empty when lambda_c >= 10 A (T_pass NaN, cold FoM 0).
Cold figure of merit: flux fraction beyond lambda_c x pass-band transmission.

-> the Bragg-section JSON in data/derived
"""
import json
import sys
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb
from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

from filteratlas.common.config import cfg  # noqa: E402

_BR = cfg("figures_of_merit")["bragg_section"]
L_CM, T_K, T_MOD, LAM_MAX = _BR["thickness_cm"], _BR["temperature_K"], _BR["moderator_T_K"], _BR["lambda_max_A"]


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    t1 = pd.read_parquet(F("stage1_triviality.parquet"), columns=["material_id", "phonon_method", "phonon_id"])
    a = at[(at.tier == "A") & at.vdos_source.ne("none") & (at["T"] == T_K)].merge(t1, on="material_id", how="left")
    rows = []
    for r in a.itertuples():
        try:
            path = (D("mp") / "phonon_dos" / r.phonon_method / f"{r.phonon_id}.json.gz" if r.vdos_source == "mp_dft"
                    else D("mlip_atlas") / "sevennet_mf_ompa" / f"{r.material_id}.json.gz")
            cfg = nb.build_vdos_material(mq.load_dos(path), f"brg_{r.material_id}")
            info = NC.createInfo(f"{cfg};temp={T_K}K")
            lc = float(info.braggthreshold)
            if lc < LAM_MAX:
                lam = np.linspace(lc * 1.0001, LAM_MAX, 400)
                mu, _ = xe.mu_offbragg(cfg, T_K, lam)
                w = xe.maxwellian_flux(lam, xe.lambda_T(T_MOD))
                tp = float(np.trapezoid(w * np.exp(-mu * L_CM), lam) / np.trapezoid(w, lam))
                mu_lc = float(mu[0])
            else:   # the band lambda_c < lambda <= 10 A is empty
                tp, mu_lc = np.nan, np.nan
            lw = np.linspace(0.5, LAM_MAX, 800)
            ww = xe.maxwellian_flux(lw, xe.lambda_T(T_MOD))
            frac = float(np.trapezoid(ww[lw > lc], lw[lw > lc]) / np.trapezoid(ww, lw))
            rows.append({"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol,
                         "lambda_c_A": lc, "T_pass_5cm_77K": tp, "flux_fraction_beyond_lambda_c": frac,
                         "cold_fom": frac * tp if lc < LAM_MAX else 0.0,   # transmitted fraction of the cold flux
                         "mu_at_lambda_c_77K": mu_lc, "vdos_source": r.vdos_source})
        except Exception as e:  # noqa: BLE001
            rows.append({"material_id": r.material_id, "formula": r.formula, "error": repr(e)[:200]})
    d = pd.DataFrame(rows)
    ok = d.dropna(subset=["lambda_c_A"]) if "lambda_c_A" in d else d
    from scipy.stats import spearmanr
    band = ok[ok.lambda_c_A < LAM_MAX]
    res = {"definition": __doc__.split("\n\n")[0], "n": len(d), "n_ok": len(ok),
           "n_empty_band_lambda_c_ge_10A": int((ok.lambda_c_A >= LAM_MAX).sum()),
           "rho_lambda_c_vs_T_pass": round(float(spearmanr(band.lambda_c_A, band.T_pass_5cm_77K).statistic), 3),
           "rho_lambda_c_vs_cold_fom": round(float(spearmanr(ok.lambda_c_A, ok.cold_fom).statistic), 3),
           "rho_flux_fraction_vs_cold_fom": round(float(spearmanr(ok.flux_fraction_beyond_lambda_c, ok.cold_fom).statistic), 3),
           "top5_cold_fom": [f"{r.formula} {r.sg_symbol}" for r in ok.nlargest(5, "cold_fom").itertuples()],
           "rows": d.round(4).to_dict("records")}
    (F("stage3_bragg.json")).write_text(json.dumps(res, indent=1, default=str))
    print({k: v for k, v in res.items() if k != "rows"})


if __name__ == "__main__":
    main()
