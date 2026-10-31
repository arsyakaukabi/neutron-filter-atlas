"""Model vs EXFOR benchmark and the acceptance criterion against measured total cross-sections (step 05, validation).

For every EXFOR dataset (subentry x temperature) and every available model of its crystal:
  curated NCrystal stdlib | MP DFT VDOS (2 meV floor, the atlas model) | Debye theta_D from the MP elastic tensor |
  Debye-300 K fallback | MLIP VDOS per candidate (<mlip dir>/<model>/<material_id>.json.gz, from
  filteratlas.dynamics.mlip_run)
compute sigma_tot per atom at the dataset temperature, with the physics of the sample form (single crystal:
off-Bragg baseline vs lower data envelope; polycrystal/powder: full NCrystal incl. Bragg vs binned median),
on a 20-bin/decade log-E grid. Metrics: r_d (RMS ln ratio over 1-100 meV), bias, failing ranges; per crystal
R_c, dataset-spread floor S_c, fail flag; the acceptance criterion. Also the wiring check (curated sapphire vs the
three EXFOR sets of NCrystal's datacmp plot) as a figure input.
-> per dataset x model and binned parquet tables in data/derived, and the benchmark summary JSON
"""
import json
import sys
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.ingest import mp as mq
from filteratlas.xs import ncmat as nb  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

from filteratlas.common.config import cfg  # noqa: E402

_BM = cfg("benchmark")
E_LO, E_HI, E_KILL = _BM["energy_window_eV"][0], _BM["energy_window_eV"][1], _BM["kill_energy_split_eV"]
EDGES = 10 ** np.arange(np.log10(E_LO), np.log10(E_HI) + 1e-9, 1 / _BM["bins_per_decade"])
KILL_RMS, KILL_N = _BM["kill_rms"], _BM["kill_n_failing"]
SC_FORMS = set(_BM["single_crystal_forms"])
PC_FORMS = set(_BM["polycrystal_forms"])

# per crystal, fixed before the computation: phase; curated file; whether the crystal counts for the acceptance criterion
CRYSTALS = {k: tuple(v) for k, v in _BM["crystals"].items()}
SENSITIVITY_PHASES = dict(_BM["sensitivity_phases"])
ENV_Q = _BM["envelope_percentile"]                                 # single-crystal lower envelope (q20), fixed before the computation
ENV_Q_ALL = (ENV_Q, *_BM["envelope_percentile_sensitivity"])       # loop order (20, 10, 50)
CHOSEN_MLIP = cfg("mlip")["atlas_model"]   # the MLIP chosen by the MLIP benchmark for the atlas


AUDIT = "--no-audit" not in sys.argv
SUFFIX = "" if AUDIT else "_preaudit"
FREE_WIN, FREE_TOL = tuple(_BM["free_atom_window_eV"]), _BM["free_atom_tolerance"]


def run(audit=True):
    """Entry point (replaces --no-audit): audit=False writes the superseded pre-audit variant (*_preaudit)."""
    global AUDIT, SUFFIX
    AUDIT, SUFFIX = audit, ("" if audit else "_preaudit")
    return main()


def free_atom_sigma(formula, E):
    """Model-independent epithermal limit: sum_i x_i [sigma_bound,i (A_i/(A_i+1))^2 + sigma_abs,i(E)] per atom."""
    from pymatgen.core import Composition
    comp = Composition(formula).fractional_composition
    tot = np.zeros_like(E, dtype=float)
    for el, x in comp.items():
        a = NC.atomDB(str(el))
        A = a.averageMassAMU() / 1.00866491595
        tot += x * ((a.coherentXS() + a.incoherentXS()) * (A / (A + 1)) ** 2 + a.captureXS() * np.sqrt(0.0253 / E))
    return tot


def audit(ex_all):
    """Data audit (model-independent): (1) superseded subentries (EXFOR STATUS SPSDD), (2) exact duplicates,
    (3) normalisation vs the free-atom limit over 1-20 eV (factor 1, n or 1/n atoms per formula; else excluded)."""
    rows = {}
    seen = {}
    for (sub, T), g in ex_all.groupby(["subentry", "temp_K"]):
        key = (int(sub), float(T))
        if "SPSDD" in str(g.status.iloc[0]):
            rows[key] = {"audit": "superseded", "renorm": 1.0}
            continue
        sig = (g.formula.iloc[0], tuple(np.round(g.E_eV.values, 9)), tuple(np.round(g.sigma_b.values, 6)))
        if sig in seen:
            rows[key] = {"audit": f"duplicate_of_{seen[sig][0]}", "renorm": 1.0}
            continue
        seen[sig] = key
        w = g[(g.E_eV >= FREE_WIN[0]) & (g.E_eV <= FREE_WIN[1])]
        if len(w) < 3:
            rows[key] = {"audit": "normalisation_unverifiable", "renorm": 1.0}
            continue
        r = float(np.median(w.sigma_per_atom_b.values / free_atom_sigma(g.formula.iloc[0], w.E_eV.values)))
        n = float(g.atoms_per_formula.iloc[0])
        fits = [f for f in sorted({1.0, n, 1.0 / n}) if abs(r * f - 1) <= FREE_TOL]
        rows[key] = ({"audit": "ok" if fits[0] == 1.0 else f"renormalised_x{fits[0]:g}", "renorm": fits[0],
                      "free_atom_ratio": round(r, 3)} if fits else
                     {"audit": "normalisation_failed", "renorm": 1.0, "free_atom_ratio": round(r, 3)})
    return rows


def datasets():
    ex_all = pd.read_csv(F("exfor_benchmark.csv"))
    ex_all = ex_all[ex_all.formula.isin(CRYSTALS) & (ex_all.sigma_per_atom_b > 0)]
    aud = audit(ex_all) if AUDIT else {}
    ex = ex_all[(ex_all.E_eV >= E_LO) & (ex_all.E_eV <= E_HI)]
    out = []
    for (sub, T), g in ex.groupby(["subentry", "temp_K"]):
        a = aud.get((int(sub), float(T)), {"audit": "not_audited", "renorm": 1.0})
        g = g.assign(sigma_per_atom_b=g.sigma_per_atom_b * a["renorm"])
        n_kill = int(((g.E_eV >= E_LO) & (g.E_eV <= E_KILL)).sum())
        form = g.sample_form.iloc[0]
        qual = ((n_kill >= 10) and (15 <= T <= 600) and (form in SC_FORMS | PC_FORMS)
                and a["audit"] not in ("superseded", "normalisation_failed") and not a["audit"].startswith("duplicate"))
        out.append({"subentry": int(sub), "entry": int(g.entry.iloc[0]), "formula": g.formula.iloc[0], "T": float(T),
                    "form": form, "temp_assumed": bool(g.temp_assumed.iloc[0]), "author": str(g.author.iloc[0]),
                    "n_points": len(g), "n_kill_range": n_kill, "qualifies": bool(qual), "audit": a["audit"],
                    "renorm": a["renorm"], "free_atom_ratio": a.get("free_atom_ratio"),
                    "counts_for_kill": bool(qual and CRYSTALS[g.formula.iloc[0]][2]), "data": g})
    return out


def mp_entry(formula, sg):
    """(phonon row from the MP phonon table of step 02 or None, summary row) for the lowest-E_hull entry of formula/sg."""
    s = pd.read_parquet(D("mp") / "summary_exp_nearhull.parquet")
    s = s[(s.formula == formula) & (s.sg_symbol == sg)].sort_values("e_hull")
    if s.empty:
        return None, None
    t = pd.read_parquet(F("stage1_triviality.parquet"))
    t = t[t.material_id.isin(s.material_id) & (t.status == "ok")]
    ph = t.merge(s[["material_id", "e_hull"]], on="material_id").sort_values("e_hull").head(1)
    return (ph.iloc[0] if len(ph) else None), s.iloc[0]


def build_models(formula):
    sg, curated, _ = CRYSTALS[formula]
    models = {"curated": f"stdlib::{curated}"}
    ph, summ = mp_entry(formula, sg)
    struct = None
    if ph is not None:
        dos = mq.load_dos(D("mp") / "phonon_dos" / ph.phonon_method / f"{ph.phonon_id}.json.gz")
        models["mp_dft"] = nb.build_vdos_material(dos, f"s2_{formula}_mpdft")
        struct = dos["structure"]
    for alt_name, alt_sg in ([("mp_dft_alt_" + SENSITIVITY_PHASES[formula], SENSITIVITY_PHASES[formula])]
                             if formula in SENSITIVITY_PHASES else []):
        aph, _ = mp_entry(formula, alt_sg)
        if aph is not None:
            dos = mq.load_dos(D("mp") / "phonon_dos" / aph.phonon_method / f"{aph.phonon_id}.json.gz")
            models[alt_name] = nb.build_vdos_material(dos, f"s2_{formula}_alt")
    if summ is not None:
        if struct is None:
            struct = mq.structures([summ.material_id])[summ.material_id]
        models["debye300"] = nb.build_debye_material(struct, 300.0, f"s2_{formula}_d300")
        if pd.notna(summ.K_vrh) and pd.notna(summ.G_vrh) and summ.K_vrh > 0 and summ.G_vrh > 0:
            th = nb.debye_temperature(summ.K_vrh, summ.G_vrh, summ.density, summ.nsites, summ.volume)
            models["debye_elastic"] = nb.build_debye_material(struct, th, f"s2_{formula}_delastic")
    # MLIP phonon DOS written by filteratlas.dynamics.mlip_run in the MP DOS format (relaxed structure + site-projected DOS),
    # converted by the identical code path as the DFT DOS (2 meV floor); keyed by the MP material id of the phase
    mid = (ph.material_id if ph is not None else (summ.material_id if summ is not None else None))
    for f in sorted((D("mlip")).glob(f"*/{mid}.json.gz")) if mid else []:
        if f.parent.name == "dft_fc":
            continue
        models[f"mlip_{f.parent.name}"] = nb.build_vdos_material(mq.load_dos(f), f"s2_{formula}_{f.parent.name}")
    return models, (None if ph is None else ph.material_id)


def sigma(cfg, T, E, single):
    c = f"{cfg};temp={T}K"
    sc = NC.createScatter(f"{c};comp=inelas,incoh_elas" if single else c)
    return sc.xsect(ekin=E) + NC.createAbsorption(c).xsect(ekin=E)


def binned(E, data, model, single, q=ENV_Q):
    idx = np.digitize(E, EDGES) - 1
    rows = []
    for b in np.unique(idx):
        m = idx == b
        if b < 0 or b >= len(EDGES) - 1:
            continue
        d = data[m]
        dv = (np.percentile(d, q) if m.sum() >= 3 else d.min()) if single else np.median(d)
        rows.append({"bin": int(b), "E_c": float(np.sqrt(EDGES[b] * EDGES[b + 1])), "n": int(m.sum()),
                     "data": float(dv), "model": float(np.median(model[m]))})
    return pd.DataFrame(rows)


def main():
    ds = datasets()
    models_by = {}
    for f in CRYSTALS:
        try:
            models_by[f] = build_models(f)
        except Exception as e:  # noqa: BLE001
            print("model build failed", f, repr(e)[:200], flush=True)
            models_by[f] = ({}, None)
    rows, bins = [], []
    for d in ds:
        g = d["data"].sort_values("E_eV")
        E, S = g.E_eV.values, g.sigma_per_atom_b.values
        single = d["form"] in SC_FORMS
        for mname, cfg in models_by[d["formula"]][0].items():
            try:
                M = sigma(cfg, d["T"], E, single)
            except Exception as e:  # noqa: BLE001
                rows.append({**{k: v for k, v in d.items() if k != "data"}, "model": mname, "error": repr(e)[:200]})
                continue
            for q in (ENV_Q_ALL if single else (50,)):
                b = binned(E, S, M, single, q)
                b["lr"] = np.log(b.model / b.data)
                k = b[(b.E_c >= E_LO) & (b.E_c <= E_KILL)]
                bad = k[np.abs(k.lr) > KILL_RMS]
                rec = {**{kk: v for kk, v in d.items() if kk != "data"}, "model": mname, "envelope_q": q,
                       "n_bins_kill": len(k), "r_d": float(np.sqrt(np.mean(k.lr ** 2))) if len(k) else np.nan,
                       "bias": float(k.lr.mean()) if len(k) else np.nan,
                       "r_d_1meV_1eV": float(np.sqrt(np.mean(b.lr ** 2))),
                       "fail_E_ranges_meV": [[round(1e3 * x, 2) for x in (bad.E_c.min(), bad.E_c.max())]] if len(bad) else []}
                rows.append(rec)
                if q in (ENV_Q, 50) and (q == ENV_Q or not single):
                    bins.append(b.rename(columns={"model": "sigma_model", "data": "sigma_data"})
                                .assign(subentry=d["subentry"], T=d["T"], formula=d["formula"], model=mname,
                                        form=d["form"]))
        print(d["formula"], d["subentry"], d["T"], d["form"], "done", flush=True)
    df = pd.DataFrame(rows)
    df["fail_E_ranges_meV"] = df["fail_E_ranges_meV"].astype(str)
    df.to_parquet(F(f"stage2_exfor_benchmark{SUFFIX}.parquet"))
    pd.concat(bins, ignore_index=True).to_parquet(F(f"stage2_exfor_binned{SUFFIX}.parquet"))
    summarise(ds, df, pd.concat(bins, ignore_index=True), models_by)


def spread(ds, binsdf, formula):
    """S_c: median over same-form pairs (|dT| <= 10 K) of the RMS ln ratio of their binned data over common bins."""
    q = [d for d in ds if d["formula"] == formula and d["qualifies"]]
    any_model = binsdf[binsdf.formula == formula].model.iloc[0] if (binsdf.formula == formula).any() else None
    if any_model is None:
        return None, []
    bd = binsdf[(binsdf.formula == formula) & (binsdf.model == any_model)]
    pairs = []
    for i in range(len(q)):
        for j in range(i + 1, len(q)):
            a, b = q[i], q[j]
            same_form = (a["form"] in SC_FORMS) == (b["form"] in SC_FORMS)
            if not same_form or abs(a["T"] - b["T"]) > 10:
                continue
            x = bd[(bd.subentry == a["subentry"]) & (bd["T"] == a["T"])].set_index("bin")
            y = bd[(bd.subentry == b["subentry"]) & (bd["T"] == b["T"])].set_index("bin")
            com = x.index.intersection(y.index)
            com = [c for c in com if E_LO <= x.loc[c, "E_c"] <= E_KILL]
            if len(com) >= 5:
                r = np.log(x.loc[com, "sigma_data"].values / y.loc[com, "sigma_data"].values)
                pairs.append({"a": a["subentry"], "b": b["subentry"], "T": a["T"], "n_bins": len(com),
                              "rms": round(float(np.sqrt(np.mean(r ** 2))), 4)})
    return (float(np.median([p["rms"] for p in pairs])) if pairs else None), pairs


def summarise(ds, df, binsdf, models_by):
    ok = df[df.get("error").isna() if "error" in df else slice(None)]
    prim = ok[(ok.envelope_q == ENV_Q) | ~ok.form.isin(SC_FORMS)]
    prim = prim[prim.qualifies]
    res = {"grid_bins_per_decade": 20, "kill_threshold_rms_ln": KILL_RMS, "kill_n_crystals": KILL_N,
           "datasets": {"n_total": len(ds), "n_qualifying": sum(d["qualifies"] for d in ds),
                        "n_counting": sum(d["counts_for_kill"] for d in ds)},
           "models_available": {f: sorted(m[0]) for f, m in models_by.items()},
           "mp_material_ids": {f: m[1] for f, m in models_by.items()}, "crystals": {}}
    for f in CRYSTALS:
        S, pairs = spread(ds, binsdf, f)
        thr = max(KILL_RMS, S) if S is not None else KILL_RMS
        g = prim[prim.formula == f]
        per_model = {}
        for mname, gm in g.groupby("model"):
            R = float(np.sqrt(np.mean(gm.r_d ** 2)))
            per_model[mname] = {"R_c": round(R, 4), "n_datasets": len(gm), "fails": bool(R > thr),
                                "mean_bias": round(float(gm.bias.mean()), 4)}
        atlas = ("mp_dft" if "mp_dft" in per_model else
                 f"mlip_{CHOSEN_MLIP}" if CHOSEN_MLIP and f"mlip_{CHOSEN_MLIP}" in per_model else None)
        res["crystals"][f] = {"phase": CRYSTALS[f][0], "counts_for_kill": CRYSTALS[f][2], "S_c": None if S is None else round(S, 4),
                              "spread_pairs": pairs, "threshold": round(thr, 4), "atlas_model": atlas or "pending_mlip",
                              "per_model": per_model}
    fails = [f for f, c in res["crystals"].items() if c["counts_for_kill"] and c["atlas_model"] != "pending_mlip"
             and c["per_model"].get(c["atlas_model"], {}).get("fails")]
    evaluated = [f for f, c in res["crystals"].items() if c["counts_for_kill"] and c["atlas_model"] != "pending_mlip"
                 and c["atlas_model"] in c["per_model"]]
    pending = [f for f, c in res["crystals"].items() if c["counts_for_kill"] and c["atlas_model"] == "pending_mlip"]
    res["kill2"] = {"evaluated_crystals": evaluated, "failing": fails, "n_failing": len(fails),
                    "pending_mlip": pending, "fires_already": len(fails) >= KILL_N,
                    "can_still_fire": len(fails) + len(pending) >= KILL_N}
    # scenario per candidate MLIP (the crystals without MP DFT phonons take that MLIP; DFT crystals are unchanged)
    cands = sorted({m for c in res["crystals"].values() for m in c["per_model"] if m.startswith("mlip_")})
    no_dft = [f for f, c in res["crystals"].items() if c["counts_for_kill"] and "mp_dft" not in c["per_model"]]
    dft_fails = [f for f in fails if f not in no_dft]
    res["kill2_if_mlip"] = {}
    for m in cands:
        mf = dft_fails + [f for f in no_dft if res["crystals"][f]["per_model"].get(m, {}).get("fails")]
        ne = [f for f in no_dft if m not in res["crystals"][f]["per_model"]]
        res["kill2_if_mlip"][m] = {"failing": mf, "n_failing": len(mf), "not_evaluable": ne,
                                   "fires": len(mf) >= KILL_N,
                                   "R_c_no_dft": {f: res["crystals"][f]["per_model"].get(m, {}).get("R_c") for f in no_dft}}
    res["audit_applied"] = AUDIT
    lab = [("renormalised" if d["audit"].startswith("renormalised") else "duplicate" if d["audit"].startswith("duplicate")
            else d["audit"]) for d in ds]
    res["audit_summary"] = {k: lab.count(k) for k in sorted(set(lab))}
    res["audit"] = [{k: d[k] for k in ("formula", "subentry", "T", "form", "audit", "renorm", "free_atom_ratio", "qualifies")}
                    for d in ds]
    (F(f"stage2_exfor{SUFFIX}.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({"kill2": res["kill2"], "kill2_if_mlip": res["kill2_if_mlip"]}, indent=1))


if __name__ == "__main__":
    main()
