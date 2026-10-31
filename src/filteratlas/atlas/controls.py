"""Atlas controls, fixed before the computation.

positive: sapphire, Si, MgO, alpha-quartz -- the atlas entries' MP DOS through the atlas code path
          (nb.build_vdos_material) into the benchmark code of step 05 (exfor_benchmark.sigma / binned, q20, audited
          data); per-dataset r_d and the single-crystal R_c (single-crystal data only, reclassification rule applied)
          must reproduce the benchmark to 3 decimals.
          Bi: the known failing control (atlas row and the 15 % overstatement from the uncertainty model of step 06).
negative: CdS, LiF (natural Li), Gd2O3 forced through (control table of step 02, same code) must lie below the
          1st percentile of the atlas L0.8 at both temperatures; no H-bearing compound in the atlas.
model:    Debye-300 K fallback on sapphire vs EXFOR 31743, recomputed now, reproduces the model control of step 02 to
          3 significant digits.
drift (i): the atlas SevenNet code path (assemble.eval_mlip) on the cached benchmark random-set SevenNet outputs
          reproduces the benchmark p68 |dln L0.8| at 293 K (0.0462).
drift (ii): the 100-compound drift set (data/intermediate/mlip_atlas/, arm "drift") vs MP DFT: p68 within the
          benchmark CI.
-> the controls JSON in data/derived
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb  # noqa: E402
from filteratlas.common.paths import D, F
from filteratlas.common.constants import N_WORKERS  # noqa: E402

from filteratlas.common.config import cfg  # noqa: E402

_C3 = cfg("controls")["stage3"]
POS = dict(_C3["positive"])   # formula -> space group


def positive(at):
    from filteratlas.validation import exfor_benchmark as s2
    bench = pd.read_parquet(F("stage2_exfor_benchmark.parquet"))
    rob = json.loads((F("stage2_d11_robustness.json")).read_text())
    recl = {int(x[0]) for x in rob["reclassified"]}
    tn = json.loads((F("stage2_text_numbers.json")).read_text())["revised_table"]
    ds = [d for d in s2.datasets() if d["formula"] in POS and d["qualifies"] and d["form"] in s2.SC_FORMS]
    t1 = pd.read_parquet(F("stage1_triviality.parquet")).set_index("material_id")
    out = {}
    for f, sg in POS.items():
        row = at[(at.formula == f) & (at.sg_symbol == sg) & (at.tier == "A")].drop_duplicates("material_id")
        mid = row.material_id.iloc[0]
        ph, _ = s2.mp_entry(f, sg)
        r = t1.loc[mid]
        cfg = nb.build_vdos_material(mq.load_dos(D("mp") / "phonon_dos" / r.phonon_method / f"{r.phonon_id}.json.gz"),
                                     f"ctl3_{mid}")
        per = []
        for d in [d for d in ds if d["formula"] == f]:
            g = d["data"].sort_values("E_eV")
            b = s2.binned(g.E_eV.values, g.sigma_per_atom_b.values, s2.sigma(cfg, d["T"], g.E_eV.values, True), True, 20)
            b["lr"] = np.log(b.model / b.data)
            k = b[(b.E_c >= s2.E_LO) & (b.E_c <= s2.E_KILL)]
            r_atlas = float(np.sqrt(np.mean(k.lr ** 2)))
            ref = bench[(bench.subentry == d["subentry"]) & (bench["T"] == d["T"]) & (bench.model == "mp_dft")
                        & (bench.envelope_q == 20)].r_d
            per.append({"subentry": d["subentry"], "T": d["T"], "r_d_atlas": round(r_atlas, 4),
                        "r_d_stage2": round(float(ref.iloc[0]), 4), "reclassified_D12": d["subentry"] in recl,
                        "match_3dp": round(r_atlas, 3) == round(float(ref.iloc[0]), 3)})
        kept = [p["r_d_atlas"] for p in per if not p["reclassified_D12"]]
        Rc = float(np.sqrt(np.mean(np.square(kept))))
        out[f] = {"atlas_material_id": mid, "stage2_material_id": ph.material_id, "same_entry": mid == ph.material_id,
                  "datasets": per, "R_c_sc_atlas": round(Rc, 3), "R_c_sc_stage2": tn[f]["R_c_atlas"],
                  "R_c_match_3dp": round(Rc, 3) == round(tn[f]["R_c_atlas"], 3),
                  "atlas_L08_293K": float(row.L08.iloc[0]) if "L08" in row else None}
    out["pass"] = all(v["R_c_match_3dp"] and v["same_entry"] and all(p["match_3dp"] for p in v["datasets"])
                      for v in out.values() if isinstance(v, dict))
    unc = json.loads((F("stage2_uncertainty_model.json")).read_text())
    bi = at[(at.formula == "Bi") & (at.tier == "A")]
    out["Bi_known_failing"] = {"atlas_rows": bi[["T", "L08", "vdos_source", "soft_lattice", "heavy_soft", "sigma_model"]]
                               .to_dict("records"),
                               "L08_overstatement_percent_stage2": unc["heavy_soft_metal_L08_overstatement_percent"]}
    return out


def negative(at):
    c = pd.read_parquet(F("stage1_controls.parquet"))
    c = c[c.kind == "negative"]
    full = at[at.vdos_source.ne("none")]
    res = {}
    for tag, T in (("293K", 293.15), ("77K", 77.0)):
        p1 = float(np.percentile(full[full["T"] == T].L08, _C3["negative_below_percentile"]))
        res[tag] = {"atlas_L08_p1_cm": p1,
                    "controls": {r.control: {"L08": float(getattr(r, f"L08_{tag}")), "below_p1": bool(getattr(r, f"L08_{tag}") < p1),
                                             "source": r.source} for r in c.itertuples()}}
    from pymatgen.core import Composition
    nH = int(at.formula.map(lambda x: "H" in {e.symbol for e in Composition(x).elements}).sum())
    res["n_H_bearing_in_atlas"] = nH
    res["pass"] = nH == 0 and all(v["below_p1"] for t in ("293K", "77K") for v in res[t]["controls"].values())
    return res


def model_control():
    ref = json.loads((F("stage0_model_control.json")).read_text())
    ex = pd.read_csv(F("exfor_benchmark.csv"))
    mc = cfg("controls")["model_control"]
    d = ex[(ex.entry == mc["exfor_entry"]) & (ex.E_eV >= mc["energy_window_eV"][0]) & (ex.E_eV <= mc["energy_window_eV"][1])]
    dos = mq.load_dos(D("mp") / "phonon_dos" / mc["phonon_method"] / f"{mc['phonon_id']}.json.gz")
    ncm = nb.build_debye_material(dos["structure"], mc["debye_K"], "ctl3_sapphire_d300")
    c = f"{ncm};temp=293.15K"
    s = (NC.createScatter(f"{c};comp=inelas,incoh_elas").xsect(ekin=d.E_eV.values)
         + NC.createAbsorption(c).xsect(ekin=d.E_eV.values))
    med = float(np.median(s / d.sigma_per_atom_b.values))
    return {"median_ratio_now": round(med, 3), "stage0": ref["debye300_fallback"]["median_ratio"],
            "pass": f"{med:.3g}" == f"{ref['debye300_fallback']['median_ratio']:.3g}"}


def _p68(x):
    return float(np.percentile(np.abs(x), 68))


def drift():
    from filteratlas.atlas.assemble import eval_mlip
    from filteratlas.validation.mlip_eval import p68_ci
    t1 = pd.read_parquet(F("stage1_triviality.parquet"), columns=["material_id", "L08_293K", "L08_77K"])
    bench = pd.read_parquet(F("stage2_mlip_benchmark.parquet"))
    s2 = bench[(bench.model == "sevennet_mf_ompa") & bench.arm.str.startswith("random") & ~bench.failed]
    files = [str(D("mlip") / "sevennet_mf_ompa" / f"{m}.json.gz") for m in s2.material_id]
    with Pool(N_WORKERS) as pool:
        ev = pd.DataFrame(pool.map(eval_mlip, files, chunksize=4))
    m = ev[ev.status == "ok"].merge(t1, on="material_id", suffixes=("_mlip", "_dft"))
    p, ci = p68_ci(np.abs(np.log(m.L08_293K_mlip / m.L08_293K_dft)))
    ref = json.loads((F("stage2_mlip.json")).read_text())["models"]["sevennet_mf_ompa"]["293K"]["p68_abs_dlog"]
    out = {"i_code_path": {"n": len(m), "p68_now": p, "ci95": ci, "stage2": ref, "pass": round(p, 4) == round(ref, 4)}}
    sset = json.loads((F("stage3_mlip_set.json")).read_text())["compounds"]
    dids = [c["material_id"] for c in sset if c["arm"] == "drift"]
    dfiles = [str(D("mlip_atlas") / "sevennet_mf_ompa" / f"{i}.json.gz") for i in dids
              if (D("mlip_atlas") / "sevennet_mf_ompa" / f"{i}.json.gz").exists()]
    if dfiles:
        with Pool(N_WORKERS) as pool:
            dv = pd.DataFrame(pool.map(eval_mlip, dfiles, chunksize=4))
        dv = dv[dv.status == "ok"]
        fail = (~dv.stable.astype(bool)) | (dv.min_freq_comm_THz < -0.5) | (dv.relax_converged == False)  # noqa: E712
        ok = dv[~fail].merge(t1, on="material_id", suffixes=("_mlip", "_dft"))
        res = {"n_set": len(dids), "n_outputs": len(dv), "n_failed_B3": int(fail.sum()), "n_scored": len(ok)}
        s2u = json.loads((F("stage2_uncertainty_model.json")).read_text())["mlip_vs_dft_dlnL08"]
        lo, hi = s2u["293K"]["ci95"]                       # MLIP-benchmark interval (drift control (ii))
        p, ci = p68_ci(np.abs(np.log(ok.L08_293K_mlip / ok.L08_293K_dft)))
        res["293K"] = {"p68": p, "ci95": ci, "stage2_ci95": [lo, hi], "pass": lo <= p <= hi}
        p77, ci77 = p68_ci(np.abs(np.log(ok.L08_77K_mlip / ok.L08_77K_dft)))
        res["77K"] = {"p68": p77, "ci95": ci77, "stage2_p68": s2u["77K"]["p68_abs_dlog"], "stage2_ci95": s2u["77K"]["ci95"]}
        res["median_signed_dlog_293K"] = round(float(np.median(np.log(ok.L08_293K_mlip / ok.L08_293K_dft))), 4)
        out["ii_population"] = res
    return out


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    res = {"positive": positive(at), "negative": negative(at), "model_control": model_control(), "drift": drift()}
    res["all_pass"] = bool(res["positive"]["pass"] and res["negative"]["pass"] and res["model_control"]["pass"]
                           and res["drift"]["i_code_path"]["pass"]
                           and res["drift"].get("ii_population", {}).get("293K", {}).get("pass", False))
    (F("stage3_controls.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: (v.get("pass") if isinstance(v, dict) else v) for k, v in res.items()}, indent=1, default=str))


if __name__ == "__main__":
    main()
