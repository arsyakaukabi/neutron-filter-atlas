"""Robustness of the EXFOR benchmark restricted to single-crystal data (step 05, validation). Everything here was
added after inspection of the first comparison.

1. Reclassification rule (one model-light rule applied to EVERY qualifying single-crystal dataset): a dataset
   labelled single crystal is "polycrystal-behaving" if the off-Bragg baseline fails it (r_off > 0.15, lower envelope)
   and NCrystal's full polycrystal model fits it at least twice as well (r_poly <= r_off / 2, bin median). Uses the
   curated NCrystal files (independent of the atlas VDOS); MP DFT reported alongside.
2. Single-crystal-only dataset-spread floors S_c^SC at the same envelope percentile (q10/q20/q50), with and without the
   reclassified datasets.
3. Robustness of the count of crystals failing the acceptance criterion against measured total cross-sections
   (single-crystal data only): variants listed in VARIANTS; pyrolytic graphite always stays counted as failing
   (removing it would be a second exclusion chosen after the first comparison, in the model's favour).
4. Assumed temperatures: datasets with temp_assumed re-evaluated at T +/- 20 K with the atlas model.
5. Leave-one-dataset-out and without the Harvey & Smith entry 13761.
6. Polycrystal residuals split at the Bragg cut-off E_c = 81.8042 / lambda_c^2 meV (lambda_c = NCrystal braggthreshold;
   bins straddling E_c are dropped): mean ln(model/data) and slope vs ln(lambda) below E_c; r_d above E_c up to 100 meV.
7. Near-duplicate datasets: same formula and T, pairwise RMS of binned data < 0.01.
-> the robustness JSON in data/derived
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.validation import exfor_benchmark as S  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

ATLAS = f"mlip_{S.CHOSEN_MLIP}"
RECLASS_FACTOR = S.cfg("benchmark")["d12_reclassification"]["factor"]      # reclassification rule
DUP_TOL = S.cfg("benchmark")["d12_reclassification"]["duplicate_tolerance"]   # near-duplicate tolerance
DT = S.cfg("benchmark")["temperature_shift_K"]                                 # assumed-temperature shift: +-20 K


def r_and_bias(cfg, d, single, q, T=None):
    g = d["data"].sort_values("E_eV")
    E, D = g.E_eV.values, g.sigma_per_atom_b.values
    b = S.binned(E, D, S.sigma(cfg, T or d["T"], E, single), single, q)
    k = b[(b.E_c >= S.E_LO) & (b.E_c <= S.E_KILL)]
    lr = np.log(k.model / k.data)
    return float(np.sqrt(np.mean(lr ** 2))), float(lr.mean())


def envelopes(ds, q):
    return {(d["subentry"], d["T"]): S.binned(d["data"].E_eV.values, d["data"].sigma_per_atom_b.values,
                                               d["data"].sigma_per_atom_b.values, True, q).set_index("bin")
            for d in ds}


def spread_sc(ds, env, formula, exclude=()):
    q = [d for d in ds if d["formula"] == formula and (d["subentry"], d["T"]) not in exclude]
    rms = []
    for a, b in combinations(q, 2):
        if abs(a["T"] - b["T"]) > 10:
            continue
        x, y = env[(a["subentry"], a["T"])], env[(b["subentry"], b["T"])]
        com = [c for c in x.index.intersection(y.index) if S.E_LO <= x.loc[c, "E_c"] <= S.E_KILL]
        if len(com) >= 5:
            rms.append(float(np.sqrt(np.mean(np.log(x.loc[com, "data"].values / y.loc[com, "data"].values) ** 2))))
    return float(np.median(rms)) if rms else None


def count(rc, thr, crystals, force_fail=("C",)):
    """rc: {crystal: R_c or None}; thr: {crystal: threshold}. Crystals with no data are 'no_data' unless forced."""
    fails, ev = [], []
    for f in crystals:
        if rc.get(f) is None:
            if f in force_fail:
                fails.append(f)
                ev.append(f)
            continue
        ev.append(f)
        if rc[f] > thr[f]:
            fails.append(f)
    return {"failing": fails, "n_failing": len(fails), "n_evaluable": len(ev), "fires": len(fails) >= S.KILL_N}


def main():
    full = json.loads((F("stage2_exfor.json")).read_text())
    counted = [f for f, c in full["crystals"].items() if c["counts_for_kill"]]
    ds_all = [d for d in S.datasets() if d["qualifies"]]
    sc = [d for d in ds_all if d["form"] in S.SC_FORMS and d["formula"] != "ZrO2"]
    pc = [d for d in ds_all if d["form"] in S.PC_FORMS]
    models = {f: S.build_models(f)[0] for f in counted}
    atlas_of = {f: ("mp_dft" if "mp_dft" in models[f] else ATLAS) for f in counted}
    out = {"note": "POST HOC, referee round 1 (review/stage2_referee.md)", "reclass_rule":
           f"r_off > {S.KILL_RMS} and r_poly <= r_off/{RECLASS_FACTOR:g} (curated NCrystal)", "datasets": []}

    # per single-crystal dataset: off-Bragg r at q10/20/50 (atlas, curated) and polycrystal r (curated, mp_dft)
    rows = []
    for d in sc:
        f = d["formula"]
        rec = {"formula": f, "subentry": d["subentry"], "T": d["T"], "form": d["form"], "temp_assumed": d["temp_assumed"],
               "author": d["author"], "entry": d["entry"]}
        for q in tuple(sorted(S.ENV_Q_ALL)):
            rec[f"r_atlas_q{q}"], rec[f"bias_atlas_q{q}"] = r_and_bias(models[f][atlas_of[f]], d, True, q)
        rec["r_off_curated"], _ = r_and_bias(models[f]["curated"], d, True, S.ENV_Q)
        rec["r_poly_curated"], rec["bias_poly_curated"] = r_and_bias(models[f]["curated"], d, False, 50)
        if "mp_dft" in models[f]:
            rec["r_poly_mp_dft"], _ = r_and_bias(models[f]["mp_dft"], d, False, 50)
        rec["polycrystal_behaving"] = bool(rec["r_off_curated"] > S.KILL_RMS
                                           and rec["r_poly_curated"] <= rec["r_off_curated"] / RECLASS_FACTOR)
        if d["temp_assumed"]:
            for dt in (-DT, DT):
                rec[f"r_atlas_q20_T{dt:+d}"], _ = r_and_bias(models[f][atlas_of[f]], d, True, S.ENV_Q, T=d["T"] + dt)
        rows.append(rec)
        print(f, d["subentry"], round(rec["r_atlas_q20"], 3), "off", round(rec["r_off_curated"], 3),
              "poly", round(rec["r_poly_curated"], 3), rec["polycrystal_behaving"], flush=True)
    t = pd.DataFrame(rows)
    out["datasets"] = rows
    recl = {(r.subentry, r.T) for r in t.itertuples() if r.polycrystal_behaving}
    out["reclassified"] = sorted([[int(s), float(T)] for s, T in recl])

    def rc_of(col, excl=(), drop=None):
        g = t[[(s, T) not in excl for s, T in zip(t.subentry, t["T"])]]
        if drop is not None:
            g = g[~drop(g)]
        return {f: (float(np.sqrt(np.mean(g[g.formula == f][col] ** 2))) if (g.formula == f).any() else None)
                for f in counted}

    thr_pre = {f: full["crystals"][f]["threshold"] for f in counted}
    thr_sc, thr_sc_rc, smap = {}, {}, {}
    for q in tuple(sorted(S.ENV_Q_ALL)):
        env = envelopes(sc, q)
        s1 = {f: spread_sc(sc, env, f) for f in counted}
        s2 = {f: spread_sc(sc, env, f, exclude=recl) for f in counted}
        thr_sc[q] = {f: max(S.KILL_RMS, s1[f]) if s1[f] is not None else S.KILL_RMS for f in counted}
        thr_sc_rc[q] = {f: max(S.KILL_RMS, s2[f]) if s2[f] is not None else S.KILL_RMS for f in counted}
        smap[f"q{q}"] = {"S_c_single_crystal": s1, "S_c_single_crystal_without_reclassified": s2}
    out["spread_floors"] = smap
    flat = {f: S.KILL_RMS for f in counted}
    V = {}
    V["D11 as reported (q20, all-form S_c)"] = count(rc_of("r_atlas_q20"), thr_pre, counted, force_fail=())
    V["D11, single-crystal S_c (q20)"] = count(rc_of("r_atlas_q20"), thr_sc[S.ENV_Q], counted, force_fail=())
    for q in tuple(sorted(S.ENV_Q_ALL)):
        V[f"D11 + reclassification, single-crystal S_c, q{q} (C kept failing)"] = count(
            rc_of(f"r_atlas_q{q}", excl=recl), thr_sc_rc[q], counted)
    V["D11 + reclassification, flat 0.15 (no floor), q20"] = count(rc_of("r_atlas_q20", excl=recl), flat, counted)
    V["D11, flat 0.15 (no floor), q20, no reclassification"] = count(rc_of("r_atlas_q20"), flat, counted, force_fail=())
    # floors recomputed without the removed data
    env20 = envelopes(sc, S.ENV_Q)
    ex13761 = recl | {(d["subentry"], d["T"]) for d in sc if d["entry"] == 13761}
    thr_no13761 = {f: (max(S.KILL_RMS, v) if (v := spread_sc(sc, env20, f, exclude=ex13761)) is not None else S.KILL_RMS)
                   for f in counted}
    V["D11 + reclassification, without Harvey & Smith entry 13761 (C kept failing)"] = count(
        rc_of("r_atlas_q20", excl=recl, drop=lambda g: g.entry == 13761), thr_no13761, counted)
    for dt in (-DT, DT):
        col = f"r_atlas_q20_T{dt:+d}"
        t[col + "_eff"] = t[col].fillna(t["r_atlas_q20"]) if col in t else t["r_atlas_q20"]
        V[f"D11 + reclassification, assumed temperatures {dt:+d} K"] = count(rc_of(col + "_eff", excl=recl),
                                                                             thr_sc_rc[20], counted)
    loo = []
    keep = t[[(s, T) not in recl for s, T in zip(t.subentry, t["T"])]]
    for r in keep.itertuples():
        g = keep[~((keep.subentry == r.subentry) & (keep["T"] == r.T))]
        rc = {f: (float(np.sqrt(np.mean(g[g.formula == f].r_atlas_q20 ** 2))) if (g.formula == f).any() else None)
              for f in counted}
        exl = recl | {(r.subentry, r.T)}
        thr_l = {f: (max(S.KILL_RMS, v) if (v := spread_sc(sc, env20, f, exclude=exl)) is not None else S.KILL_RMS)
                 for f in counted}
        c = count(rc, thr_l, counted)
        loo.append({"left_out": [int(r.subentry), float(r.T), r.formula], **c})
    V["leave-one-dataset-out (reclassified, single-crystal S_c, q20): range"] = {
        "min": min(x["n_failing"] for x in loo), "max": max(x["n_failing"] for x in loo),
        "cases_changing_count": [x for x in loo if x["n_failing"] != V[
            "D11 + reclassification, single-crystal S_c, q20 (C kept failing)"]["n_failing"]]}
    out["variants"] = V
    out["single_dataset_crystals"] = sorted(f for f in counted if (keep.formula == f).sum() == 1)
    out["entry_13761_only_crystals"] = sorted(f for f in counted if (keep.formula == f).any()
                                              and set(keep[keep.formula == f].entry) == {13761})

    # 6. polycrystal datasets split at the Bragg cut-off
    bn = pd.read_parquet(F("stage2_exfor_binned.parquet"))
    split = []
    for d in pc:
        f = d["formula"]
        if f not in models or "curated" not in models[f]:
            continue
        lam_c = float(NC.createInfo(models[f]["curated"]).braggthreshold)
        Ec = 81.8042 / lam_c ** 2 * 1e-3                          # eV
        for m in ("mp_dft", "curated", ATLAS):
            g = bn[(bn.subentry == d["subentry"]) & (bn["T"] == d["T"]) & (bn.model == m)].sort_values("E_c")
            if g.empty:
                continue
            lo_edge, hi_edge = S.EDGES[g["bin"].values], S.EDGES[g["bin"].values + 1]
            lo = g[(hi_edge <= Ec) & (g.E_c >= S.E_LO)]           # bins entirely below the cut-off
            hi = g[(lo_edge >= Ec) & (g.E_c <= S.E_KILL)]         # bins entirely above; straddling bins dropped
            rec = {"formula": f, "subentry": d["subentry"], "T": d["T"], "form": d["form"], "model": m,
                   "E_cutoff_meV": round(Ec * 1e3, 3), "n_below": len(lo), "n_above": len(hi),
                   "mean_ln_below": round(float(lo.lr.mean()), 3) if len(lo) else None,
                   "r_d_above": round(float(np.sqrt(np.mean(hi.lr ** 2))), 3) if len(hi) else None}
            if len(lo) >= 3:   # slope of ln(sigma_data/sigma_model) vs ln(lambda): > 0 means excess rising with lambda
                lam = 0.286 / np.sqrt(lo.E_c.values)
                rec["slope_excess_vs_ln_lambda_below"] = round(float(np.polyfit(np.log(lam), -lo.lr.values, 1)[0]), 3)
            split.append(rec)
    out["polycrystal_cutoff_split"] = split

    # 7. near-duplicates among all qualifying datasets
    env = envelopes(ds_all, 50)
    dups = []
    for a, b in combinations(ds_all, 2):
        if a["formula"] != b["formula"] or abs(a["T"] - b["T"]) > 0.5:
            continue
        x, y = env[(a["subentry"], a["T"])], env[(b["subentry"], b["T"])]
        com = x.index.intersection(y.index)
        if len(com) >= 5:
            r = float(np.sqrt(np.mean(np.log(x.loc[com, "data"].values / y.loc[com, "data"].values) ** 2)))
            if r < DUP_TOL:
                dups.append({"formula": a["formula"], "a": a["subentry"], "b": b["subentry"], "T": a["T"], "rms": round(r, 5)})
    out["near_duplicates"] = dups
    (F("stage2_d11_robustness.json")).write_text(json.dumps(out, indent=1, default=str))
    for k, v in V.items():
        print(k, "->", v if "range" in k else (v["n_failing"], v["failing"], v["fires"]))
    print("reclassified", out["reclassified"], "single-dataset", out["single_dataset_crystals"],
          "13761-only", out["entry_13761_only_crystals"], "dups", dups)


if __name__ == "__main__":
    main()
