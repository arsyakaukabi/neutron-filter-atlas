"""Dense-mesh dynamical-stability screen from MP force constants.

Sample (fixed seed): usable compounds stratified by nsites (<=10, 11-20, 21-40, >40; 50 random each), plus the
top-30 of the full-model and proxy rankings (one polymorph per formula, non-molecular) at both temperatures,
plus all anchor / control compounds. For each: phonopy from MP force constants (symmetrized), frequencies on a
40 A q-mesh plus Gamma; unstable if min frequency < -0.5 THz (the frequency threshold fixed before the computation;
-0.1 THz reported as well). Outputs: prevalence per stratum with Wilson 95 % CI, the triviality-test rho on the random
sample with and without dense-mesh-unstable compounds, and stability of the top lists.
-> stage1r_stability_screen.parquet + .json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from filteratlas.ingest import mp as mq  # noqa: E402
from filteratlas.common.config import cfg, tup
from filteratlas.common.constants import UNSTABLE_WEIGHT_MAX
from filteratlas.common.paths import D, F  # noqa: E402

_S = cfg("stability")
FCDIR = D("mp") / "force_constants"
STRATA = [tuple(x) for x in _S["screen_strata_nsites"]]
N_PER = _S["screen_n_per_stratum"]
DEFAULT_MESH = _S["fc_mesh_length_A"]
MESH = DEFAULT_MESH
SUFFIX = ""   # mesh-density sensitivity runs get their own files (set by run(mesh=...))
SEED = _S["screen_seed"]
FC_UNSTABLE = _S["fc_unstable_THz"]
FC_SECONDARY = _S["fc_secondary_THz"]   # reported alongside (stricter)
FCDIR.mkdir(parents=True, exist_ok=True)


def run(all_usable=False, mesh=None):
    """Entry point (replaces the --all / --mesh command-line flags): the stratified screen, or all usable compounds;
    mesh = q-mesh length in Angstrom (default 40, the dynamical-stability rule; 60 is the mesh-density sensitivity
    run)."""
    global MESH, SUFFIX
    MESH = float(mesh) if mesh is not None else DEFAULT_MESH
    SUFFIX = "" if MESH == DEFAULT_MESH else f"_mesh{MESH:g}"
    return main_all() if all_usable else main()


def wilson(k, n, z=1.96):
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(float(c - h), 3), round(float(c + h), 3)]


def fetch_fc(rows, batch=10):
    """Bulk-fetch force constants + supercell/primitive matrices + structure; cache as npz per phonon id."""
    from emmet.core.mpid import AlphaID
    need = [r for r in rows if not (FCDIR / f"{r['phonon_id']}.npz").exists()]
    if not need:
        return
    from filteratlas.common.constants import MP_KEY
    if not MP_KEY:   # offline rerun: ids MP returned no force constants for stay absent -> recorded as no_force_constants
        print(f"fetch_fc: {len(need)} force-constant files absent and no MP_API_KEY; not fetched", flush=True)
        return
    with mq.rester() as m:
        ph = m.materials.phonon
        lbl, _ = ph._get_delta_table("materialsproject-parsed", "phonon/force-constants/", label="ph_force_constants")
        # supercell/primitive matrices and the calculation structure from the phonon docs (legacy ids)
        docs = {}
        for i in range(0, len(need), 200):
            chunk = need[i:i + 200]
            leg = {str(d.model_dump()["material_id"]): str(d.material_id) for d in
                   m.materials.summary.search(material_ids=[r["material_id"] for r in chunk], fields=["material_id"])}
            for d in ph.search(material_ids=list(leg.values()), fields=["supercell_matrix", "primitive_matrix",
                                                                         "structure", "phonon_method", "identifier"]):
                docs[(str(d.identifier), str(getattr(d.phonon_method, "value", d.phonon_method)))] = d
        for meth in sorted({r["phonon_method"] for r in need}):
            ids = [r["phonon_id"] for r in need if r["phonon_method"] == meth]
            for i in range(0, len(ids), batch):
                chunk = ids[i:i + batch]
                key = {str(AlphaID(p, padlen=8)): p for p in chunk}
                q = (f"SELECT identifier, force_constants FROM {lbl} WHERE identifier IN "
                     f"({','.join(repr(k) for k in key)}) AND phonon_method='{meth}'")
                for attempt in range(4):
                    try:
                        rows_ = ph._query_delta_single(q, label=lbl).to_pylist(maps_as_pydicts="strict")
                        break
                    except Exception:  # noqa: BLE001
                        time.sleep(10 * (attempt + 1))
                else:
                    continue
                for r in rows_:
                    pid = key.get(str(r["identifier"]))
                    d = docs.get((pid, meth))
                    if pid is None or d is None or r.get("force_constants") is None:
                        continue
                    s = d.structure
                    np.savez_compressed(FCDIR / f"{pid}.npz", fc=np.array(r["force_constants"], dtype=float),
                                        sc=np.array(d.supercell_matrix), prim=np.eye(3) if d.primitive_matrix is None else np.array(d.primitive_matrix),
                                        lattice=s.lattice.matrix, frac=s.frac_coords,
                                        symbols=np.array([str(x.specie) for x in s], dtype="U3"))
                print(f"[fc] {meth} {i + len(chunk)}/{len(ids)}", flush=True)


def min_freq(pid):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    z = np.load(FCDIR / f"{pid}.npz", allow_pickle=False)
    unit = PhonopyAtoms(symbols=list(z["symbols"]), cell=z["lattice"], scaled_positions=z["frac"])
    ph = Phonopy(unit, supercell_matrix=z["sc"], primitive_matrix=z["prim"])
    ph.force_constants = z["fc"]
    ph.symmetrize_force_constants()
    ph.run_mesh(MESH, is_gamma_center=True)
    fmesh = float(np.min(ph.mesh.frequencies))
    # Gamma: drop the three acoustic modes = the three of smallest |nu| (an imaginary optical mode sorts first,
    # so slicing [3:] of the ascending list would discard it)
    ph.run_qpoints([[0, 0, 0]])
    fg = np.asarray(ph.get_qpoints_dict()["frequencies"][0])
    fgam = float(np.min(fg[np.argsort(np.abs(fg))[3:]])) if len(fg) > 3 else 0.0
    # supercell-commensurate q-points: frequencies there follow directly from the DFT forces (no Fourier
    # interpolation); an instability seen only between them is likely an interpolation artefact
    from phonopy.harmonic.dynmat_to_fc import get_commensurate_points
    smat = np.rint(np.linalg.inv(ph.primitive_matrix) @ ph.supercell_matrix).astype(int)
    qc = get_commensurate_points(smat)
    ph.run_qpoints(qc)
    fc_ = np.asarray(ph.get_qpoints_dict()["frequencies"])
    gam = np.all(np.abs(qc) < 1e-8, axis=1)
    vals = [fc_[i][np.argsort(np.abs(fc_[i]))[3:]] if gam[i] else fc_[i] for i in range(len(qc))]
    fcomm = float(min(v.min() for v in vals if len(v)))
    return min(fmesh, fgam), fgam, len(ph.primitive), fcomm


def main():
    u = pd.read_parquet(F("stage1_usable.parquet"))
    rng = np.random.default_rng(SEED)
    sel = []
    for lo, hi in STRATA:
        pool = u[(u.nsites >= lo) & (u.nsites <= hi)]
        take = pool.sort_values("material_id").sample(min(N_PER, len(pool)), random_state=int(rng.integers(1e9)))
        sel += [{**r, "arm": f"random_{lo}-{hi if hi < 10_000 else 'max'}"} for r in take.to_dict("records")]
    st = json.loads((F("stage1_stats.json")).read_text())
    tops = set()
    for t in ("293K", "77K"):
        both = u[~u.molecular].sort_values("e_hull").groupby("formula").head(1)
        tops |= set(both.nlargest(30, f"L08_{t}").material_id) | set(both.nlargest(30, f"L08_proxy_{t}").material_id)
    anchors = {"Al2O3": "R-3c", "Si": "Fd-3m", "MgO": "Fm-3m", "MgF2": "P4_2/mnm", "BeF2": "P3_121", "BeO": "P6_3mc",
               "PbF2": "Fm-3m", "SiO2": None}
    for f, sg in anchors.items():
        m_ = (u.formula == f) & ((u.sg_symbol == sg) if sg else u.sg_symbol.isin(["P3_121", "P3_221"]))
        tops |= set(u[m_].material_id)
    have = {r["material_id"] for r in sel}
    sel += [{**r, "arm": "top_or_anchor"} for r in u[u.material_id.isin(tops - have)].to_dict("records")]
    for r in sel:
        r["in_top_or_anchor"] = r["material_id"] in tops
    fetch_fc(sel)
    rows = []
    for r in sel:
        rec = {k: r[k] for k in ("material_id", "formula", "sg_symbol", "nsites", "arm", "in_top_or_anchor",
                                 "unstable_weight", "L08_293K", "L08_77K", "L08_proxy_293K", "L08_proxy_77K")}
        try:
            fmin, fgam, nprim, fcomm = min_freq(r["phonon_id"])
            rec.update({"fc_min_freq_THz": fmin, "fc_gamma_min_optical_THz": fgam, "nprim": nprim,
                        "fc_commensurate_min_THz": fcomm, "status": "ok"})
        except Exception as e:  # noqa: BLE001
            rec.update({"status": "error", "error": repr(e)[:200]})
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["unstable_05"] = df.fc_min_freq_THz < FC_UNSTABLE
    df["unstable_01"] = df.fc_min_freq_THz < FC_SECONDARY
    df.to_parquet(F("stage1r_stability_screen.parquet"))
    ok = df[df.status == "ok"]
    rnd = ok[ok.arm.str.startswith("random")]
    res = {"mesh_length_A": MESH, "n_selected": len(df), "n_ok": len(ok), "n_errors": int((df.status != "ok").sum()),
           "errors": df[df.status != "ok"].get("error", pd.Series(dtype=str)).str[:80].value_counts().to_dict()}
    prev = {}
    for arm, g in rnd.groupby("arm"):
        k, n = int(g.unstable_05.sum()), len(g)
        prev[arm] = {"n": n, "unstable_lt_-0.5THz": k, "frac": round(k / n, 3), "wilson95": wilson(k, n),
                     "unstable_lt_-0.1THz": int(g.unstable_01.sum())}
    k, n = int(rnd.unstable_05.sum()), len(rnd)
    # weights so the stratified sample estimates the prevalence in the whole usable set
    w = {arm: ((u.nsites >= int(arm.split("_")[1].split("-")[0])) &
               (u.nsites <= (10_000 if arm.endswith("max") else int(arm.split("-")[1])))).sum() / len(u)
         for arm in prev}
    prev["all_random_unweighted"] = {"n": n, "unstable": k, "frac": round(k / n, 3), "wilson95": wilson(k, n)}
    prev["usable_set_weighted_estimate"] = round(float(sum(w[a] * prev[a]["frac"] for a in w)), 3)
    prev["mp_dos_rule_missed"] = int((rnd.unstable_05 & (rnd.unstable_weight <= UNSTABLE_WEIGHT_MAX)).sum())
    res["prevalence"] = prev
    gate = {}
    for t in ("293K", "77K"):
        allr = spearmanr(rnd[f"L08_{t}"], rnd[f"L08_proxy_{t}"]).statistic
        cl = rnd[~rnd.unstable_05]
        clr = spearmanr(cl[f"L08_{t}"], cl[f"L08_proxy_{t}"]).statistic
        cl1 = rnd[~rnd.unstable_01]
        gate[t] = {"rho_random_sample_all": round(float(allr), 4), "n_all": len(rnd),
                   "rho_random_sample_stable_0.5": round(float(clr), 4), "n_stable_0.5": len(cl),
                   "rho_random_sample_stable_0.1": round(float(spearmanr(cl1[f"L08_{t}"], cl1[f"L08_proxy_{t}"]).statistic), 4),
                   "n_stable_0.1": len(cl1),
                   "delta_rho_excluding_unstable": round(float(clr - allr), 4)}
        # bootstrap the paired difference
        rngb = np.random.default_rng(5)
        diffs = []
        for _ in range(1000):
            b = rnd.iloc[rngb.integers(0, len(rnd), len(rnd))]
            bc = b[~b.unstable_05]
            diffs.append(spearmanr(bc[f"L08_{t}"], bc[f"L08_proxy_{t}"]).statistic
                         - spearmanr(b[f"L08_{t}"], b[f"L08_proxy_{t}"]).statistic)
        gate[t]["delta_ci95"] = [round(float(np.percentile(diffs, 2.5)), 4), round(float(np.percentile(diffs, 97.5)), 4)]
    res["gate_on_random_sample"] = gate
    topok = ok[ok.in_top_or_anchor]
    res["top_and_anchor_unstable"] = topok[topok.unstable_05][["formula", "sg_symbol", "fc_min_freq_THz"]].round(2).to_dict("records")
    res["n_top_and_anchor_screened"] = len(topok)
    (F("stage1r_stability_screen.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


def main_all():
    """Screen ALL usable compounds with force constants (pheasy); DFPT-only entries have no FC in the MP Delta
    table and are flagged 'unscreenable'. Gives the triviality test on the dynamically clean set directly (no
    extrapolation). -> stage1r_stability_full.parquet + .json (with mesh L: suffix _meshL)"""
    u = pd.read_parquet(F("stage1_usable.parquet"))
    rows = u.to_dict("records")
    fetch_fc([r for r in rows if r["phonon_method"] == "pheasy"])
    out = []
    for r in rows:
        rec = {"material_id": r["material_id"], "formula": r["formula"], "sg_symbol": r["sg_symbol"],
               "phonon_method": r["phonon_method"]}
        if r["phonon_method"] != "pheasy":
            rec["status"] = "unscreenable_dfpt_no_fc"
        else:
            try:
                fmin, fgam, nprim, fcomm = min_freq(r["phonon_id"])
                rec.update({"fc_min_freq_THz": fmin, "fc_gamma_min_optical_THz": fgam,
                            "fc_commensurate_min_THz": fcomm, "status": "ok"})
            except Exception as e:  # noqa: BLE001
                rec.update({"status": "error", "error": repr(e)[:200]})
        out.append(rec)
    df = pd.DataFrame(out)
    df["unstable_05"] = df.fc_min_freq_THz < FC_UNSTABLE
    df["unstable_01"] = df.fc_min_freq_THz < FC_SECONDARY
    df["unstable_comm_05"] = df.fc_commensurate_min_THz < FC_UNSTABLE      # evidenced at DFT-supercell q-points
    df["interp_only_05"] = df.unstable_05 & ~df.unstable_comm_05      # imaginary only between them
    df.to_parquet(F(f"stage1r_stability_full{SUFFIX}.parquet"))
    m = u.merge(df[["material_id", "status", "fc_min_freq_THz", "unstable_05", "unstable_01", "unstable_comm_05"]]
                .rename(columns={"status": "fc_status"}), on="material_id")   # u has its own 'status' column
    ok = m[m.fc_status == "ok"]
    res = {"mesh_length_A": MESH, "n_usable": len(u), "n_screened_ok": len(ok),
           "n_unscreenable_dfpt": int((m.fc_status == "unscreenable_dfpt_no_fc").sum()),
           "n_error": int((m.fc_status == "error").sum()),
           "errors": (df[df.status == "error"].error.str[:80].value_counts().to_dict()
                      if "error" in df else {}),
           "unstable_05": {"k": int(ok.unstable_05.sum()), "frac": round(float(ok.unstable_05.mean()), 4),
                           "wilson95": wilson(int(ok.unstable_05.sum()), len(ok))},
           "unstable_01": {"k": int(ok.unstable_01.sum()), "frac": round(float(ok.unstable_01.mean()), 4)},
           "unstable_commensurate_05": {"k": int(ok.unstable_comm_05.sum()), "frac": round(float(ok.unstable_comm_05.mean()), 4),
                                        "wilson95": wilson(int(ok.unstable_comm_05.sum()), len(ok))},
           "interpolation_only_05": {"k": int((ok.unstable_05 & ~ok.unstable_comm_05).sum()),
                                     "frac_of_flagged": round(float((ok.unstable_05 & ~ok.unstable_comm_05).sum()
                                                                    / max(ok.unstable_05.sum(), 1)), 4)}}
    variants = {
        "all_usable_(preregistered)": m,
        "exclude_dense_unstable_0.5_keep_unscreened": m[~(m.unstable_05 == True)],  # noqa: E712
        "exclude_dense_unstable_0.1_keep_unscreened": m[~(m.unstable_01 == True)],  # noqa: E712
        "exclude_commensurate_unstable_0.5_keep_unscreened": m[~(m.unstable_comm_05 == True)],  # noqa: E712
        "screened_exclude_commensurate_unstable_0.5": ok[~ok.unstable_comm_05],
        "screened_clean_only_0.5": ok[~ok.unstable_05],
        "screened_clean_only_0.1": ok[~ok.unstable_01],
    }
    g = {}
    for name, d in variants.items():
        g[name] = {"n": len(d)}
        for t in ("293K", "77K"):
            g[name][t] = round(float(spearmanr(d[f"L08_{t}"], d[f"L08_proxy_{t}"]).statistic), 4)
    # clustered bootstrap CI for the main corrected variant
    d = variants["screened_clean_only_0.5"]
    rngb = np.random.default_rng(11)
    groups = d.groupby("formula").indices
    keys = list(groups)
    for t in ("293K", "77K"):
        bs = []
        x, y = d[f"L08_{t}"].values, d[f"L08_proxy_{t}"].values
        for _ in range(1000):
            idx = np.concatenate([groups[keys[j]] for j in rngb.integers(0, len(keys), len(keys))])
            bs.append(spearmanr(x[idx], y[idx]).statistic)
        g["screened_clean_only_0.5"][f"ci95_cluster_{t}"] = [round(float(np.percentile(bs, 2.5)), 4),
                                                            round(float(np.percentile(bs, 97.5)), 4)]
    res["gate_variants"] = g
    clean = d[~d.molecular].sort_values("e_hull").groupby("formula").head(1)
    for t in ("293K", "77K"):
        res[f"top10_overlap_clean_dedup_{t}"] = len(set(clean.nlargest(10, f"L08_{t}").material_id)
                                                    & set(clean.nlargest(10, f"L08_proxy_{t}").material_id))
        res[f"top10_full_clean_{t}"] = clean.nlargest(10, f"L08_{t}")[["formula", "sg_symbol"]].to_dict("records")
    (F(f"stage1r_stability_full{SUFFIX}.json")).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    run("--all" in sys.argv, float(sys.argv[sys.argv.index("--mesh") + 1]) if "--mesh" in sys.argv else None)
