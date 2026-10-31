"""Atlas table: one row per candidate x temperature, all figures of merit, proxy,
VDOS source, sigma_model, tier, flags, funnel_step.

Full-model values:
- mp_dft: the triviality-test values of step 02 (MP DOS, 2 meV floor), the compound passing the force-constant screen;
- sevennet: SevenNet DOS (data/intermediate/mlip_atlas/sevennet_mf_ompa/) through the identical conversion and NCrystal
  code path (evaluated here; workers cache per-compound results in data/intermediate/atlas_eval/ so an interrupted
  run loses little).
Rejection: filteratlas.fom.rejection at L0.8 and at L0.8*exp(+-sigma). sigma_model from configs/uncertainty.yaml.
-> data/derived/atlas.parquet, the curve table (LAM_WIDE curves of every full-model Tier-A compound),
   data/derived/atlas_counts.json ("atlas" part)
"""
import gzip
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.core import Composition, Element

from filteratlas.fom import rejection as rj  # noqa: E402
from filteratlas.common.config import cfg
from filteratlas.common.constants import N_WORKERS, TEMPS  # noqa: E402
from filteratlas.common.paths import D, F

TAGS = ("293K", "77K")
# sigma values: configs/uncertainty.yaml (construction B; MLIP-vs-DFT p68; heavy-soft band maximum)
_U = cfg("uncertainty")
SIG_DFT = _U["sigma_model_dft"]
SIG_MLIP_ADD = dict(_U["sigma_mlip_add"])
SIG_SOFT = _U["sigma_soft_floor"]
SOFT_MEV, Z_HEAVY = _U["soft_lattice_mean_phonon_meV"], _U["heavy_Z"]
MLIP_DIR = D("mlip_atlas") / cfg("mlip")["atlas_model"]
EVAL_DIR = D("atlas_eval")
NOBLE = set(_U["noble_gases"])
_H = pd.read_csv(F("hygroscopic_families.csv"))   # data/curated
HYGRO_CATIONS = set(_H[_H.kind == "cation"].value)
HYGRO_FORMULAS = set(_H[_H.kind == "formula"].value)
HYGRO_ANIONS = set(_H[_H.kind == "anion"].value)
HYGRO_HALIDES = set(_H[_H.kind == "halide"].value)
FC_UNSTABLE = cfg("stability")["fc_unstable_THz"]
NUC_INCOMPLETE = {e: v["nucdata_incomplete"] for e, v in
                  json.loads((D("nucdata") / "endfb80_meta.json").read_text())["elements"].items()}
FLAG_COLS = ("flag_Be", "flag_Tl", "flag_As", "flag_Pb", "flag_Bi", "hygroscopic", "melting_point_unknown", "interp_only",
             "unscreenable", "unstable_dft", "not_computed_budget", "supercell_too_large", "mlip_failed",
             "tierA_phase_mismatch", "not_off_bragg_single_crystal", "soft_lattice", "heavy_soft", "untested_heavy",
             "molecular", "mlip_unstable", "not_run_memory", "nucdata_incomplete")


def hygroscopic(formula):
    """Rule-based family list (approximate)."""
    c = Composition(formula)
    if c.reduced_formula in HYGRO_FORMULAS:
        return True
    els = {e.symbol for e in c.elements}
    an = {e for e in els if e in HYGRO_ANIONS}
    return bool(an) and an <= HYGRO_HALIDES and bool(els & HYGRO_CATIONS)


def is_molecular(formula, dim):
    """Same definition as the triviality-test statistics of step 02, reused unchanged."""
    els = {e.symbol for e in Composition(formula).elements}
    if els <= NOBLE:
        return True
    metal = any(Element(e).is_metal for e in els)
    return dim in ("0D", "intercalated molecule") and not metal


def eval_mlip(path, cache_dir=None):
    """SevenNet DOS -> the same record as for MP DOS (stability, L0.8, hot band, mu, curves, mean phonon, dimensionality).
    cache_dir: per-compound cache (default intermediate/atlas_eval; the relaxation-sensitivity run uses its own)."""
    from filteratlas.ingest import mp as mq
    from filteratlas.xs import ncmat as nb
    from filteratlas.xs import evaluate as xe
    from filteratlas.xs.mp_models import mean_phonon
    path = Path(path)
    cache_dir = Path(cache_dir) if cache_dir else EVAL_DIR
    cache = cache_dir / path.name.replace(".json.gz", ".json")
    if cache.exists() and cache.stat().st_mtime > path.stat().st_mtime:
        return json.loads(cache.read_text())
    mid = path.name.replace(".json.gz", "")
    rec = {"material_id": mid}
    try:
        dos = mq.load_dos(path)
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            meta = json.load(fh)["meta"]
        rec.update({k: meta.get(k) for k in ("relax_converged", "min_freq_comm_THz", "min_freq_mesh_THz", "volume_ratio",
                                              "device", "wall_s")})
        w, ok = nb.stability(dos)
        rec.update({"unstable_weight": w, "stable": ok})
        cfg = nb.build_vdos_material(dos, f"atl_{mid}")
        rec.update(xe.evaluate(cfg, TEMPS))
        rec["mean_phonon_meV"] = mean_phonon(dos)
        for T in TEMPS:
            tag = "293K" if T > 100 else "77K"
            rec[f"curve_full_{tag}"], rec[f"curve_proxy_{tag}"] = xe.curves(cfg, T)
        from pymatgen.analysis.dimensionality import get_dimensionality_cheon
        from pymatgen.core import Structure
        try:   # pymatgen's bonding analysis can recurse too deep on large cells; that must not drop the compound
            rec["dim"] = str(get_dimensionality_cheon(Structure.from_dict(dos["structure"])))
        except RecursionError as e:
            rec["dim"] = f"error:{type(e).__name__}"
        rec["status"] = "ok"
    except Exception as e:  # noqa: BLE001
        rec.update({"status": "error", "error": repr(e)[:300]})
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(rec))
    return rec


def mlip_log_status():
    """Last logged status per compound (relaxation / phonon failures leave no DOS file)."""
    f = (D("logs") / "mlip_atlas.jsonl")
    if not f.exists():
        return {}
    out = {}
    for line in open(f, encoding="utf-8"):
        if line.strip():
            r = json.loads(line)
            out[r["material_id"]] = r
    return out


def not_run_memory_ids():
    """Compounds the scheduler logged as not_run_memory (fit neither the GPU memory nor the RAM)."""
    f = (D("logs") / "laptop_scheduler_stage3.jsonl")
    if not f.exists():
        return set()
    return {json.loads(x)["material_id"] for x in open(f, encoding="utf-8")
            if x.strip() and json.loads(x).get("status") == "not_run_memory"}


def main():
    at, full = build()
    write(at, full)


def build(ml_override=None):
    """Atlas rows and full-model records. ml_override (relaxation sensitivity only): SevenNet evaluation records that
    replace the atlas records of the same compounds before the relaxation-convergence rule is applied."""
    fun = pd.read_parquet(F("stage3_funnel.parquet"))
    prox = pd.read_parquet(F("stage3_proxy.parquet"))
    cand = pd.read_parquet(F("stage0_candidates.parquet"))[["material_id", "elements", "crystal_system", "e_hull"]]
    tiers = pd.read_parquet(F("stage3_tiers.parquet"))
    sset = json.loads((F("stage3_mlip_set.json")).read_text())
    drift_ids = {c["material_id"] for c in sset["compounds"] if c["arm"] == "drift"}
    # --- full-model values -------------------------------------------------------------------------------------
    t1 = pd.read_parquet(F("stage1_triviality.parquet"))
    dims = pd.read_parquet(F("stage1_dimensionality.parquet"))
    dft = t1[t1.material_id.isin(fun[fun.funnel_step == "mp_dft"].material_id)].merge(dims, on="material_id", how="left")
    dft["vdos_source"] = "mp_dft"
    files = sorted(str(p) for p in MLIP_DIR.glob("*.json.gz") if p.name.replace(".json.gz", "") not in drift_ids)
    t0 = time.time()
    with Pool(N_WORKERS) as pool:
        ml = pd.DataFrame(pool.map(eval_mlip, files, chunksize=4)) if files else pd.DataFrame(columns=["material_id", "status"])
    print(f"evaluated {len(ml)} SevenNet DOS in {time.time() - t0:.0f}s", flush=True)
    if ml_override is not None and len(ml_override):
        ml = pd.concat([ml[~ml.material_id.isin(ml_override.material_id)], ml_override], ignore_index=True)
    # relaxation not converged -> mlip_failed, no full-model values
    ml_ok = ml[(ml.status == "ok") & (ml.relax_converged != False)].copy()   # noqa: E712
    ml_ok["vdos_source"] = "sevennet"
    full = pd.concat([dft, ml_ok], ignore_index=True)
    # --- one row per candidate, base columns ---------------------------------------------------------------------
    base = fun.merge(prox[["material_id", "n_summary_A3", "proxy_mu_1.8"]], on="material_id").merge(cand, on="material_id")
    base = base.merge(tiers[["material_id", "tierA_match", "tierA_phase_mismatch", "not_off_bragg_single_crystal"]],
                      on="material_id", how="left")
    base["tier"] = np.where(base.tierA, "A", "B")
    # the elastic-tensor Debye temperature of compounds without MP phonons (Debye table of step 03), so that SevenNet
    # and no-model rows carry it too; MP-DFT rows keep the value from step 02
    deb = pd.read_parquet(F("stage3_debye.parquet"))[["material_id", "theta_D_elastic"]]
    base = base.merge(deb.rename(columns={"theta_D_elastic": "theta_D_debye_table"}), on="material_id", how="left")
    base = base.rename(columns={"proxy_mu_1.8": "proxy_mu_18"})   # itertuples cannot address a dotted name
    logst = mlip_log_status()
    extra = {c["material_id"] for c in sset["compounds"] if c["arm"] in ("atlas_a_large_D15", "atlas_a_D16")}
    have = {p.name.replace(".json.gz", "") for p in MLIP_DIR.glob("*.json.gz")}
    # the added large-cell and Tier-A arms count as SevenNet compounds once they have run or been tried
    sel = (fun.funnel_step == "sevennet_selected") | (fun.material_id.isin(extra) & (fun.material_id.isin(have)
                                                                                     | fun.material_id.isin(logst)))
    ok_ids = set(ml_ok.material_id)
    base.loc[sel & base.material_id.isin(ok_ids), "funnel_step"] = "sevennet"
    # out of GPU memory (or never tried because the supercell exceeds the GPU memory cap and the
    # RAM) -> not_run_memory; any other logged failure (relaxation, phonons, DOS conversion) -> mlip_failed
    oom = base.material_id.map(lambda m: m in logst and logst[m].get("status") == "error" and any(
        k in logst[m].get("error", "") for k in ("CUDA", "OutOfMemory", "cuda")))
    nrm = base.material_id.isin(not_run_memory_ids()) & ~base.material_id.isin(have)   # a later try may have run
    oom = oom & ~base.material_id.isin(have)
    pend = sel & ~base.material_id.isin(ok_ids)
    base.loc[pend & (oom | nrm), "funnel_step"] = "not_run_memory"
    failed = pend & ~(oom | nrm) & base.material_id.map(lambda m: m in logst or m in set(ml.material_id))
    base.loc[failed, "funnel_step"] = "mlip_failed"
    base.loc[pend & ~(oom | nrm) & ~failed, "funnel_step"] = "sevennet_pending"
    rows = []
    fullx = full.set_index("material_id")
    for r in base.itertuples():
        has = r.material_id in fullx.index and r.funnel_step in ("mp_dft", "sevennet")
        f = fullx.loc[r.material_id] if has else None
        els = set(r.elements.split(","))
        zmax = max(Element(e).Z for e in els)
        common = {"material_id": r.material_id, "formula": r.formula, "sg_symbol": r.sg_symbol,
                  "crystal_system": r.crystal_system, "e_hull": r.e_hull, "tier": r.tier,
                  "tierA_match": r.tierA_match, "funnel_step": r.funnel_step,
                  "vdos_source": f["vdos_source"] if has else "none",
                  "proxy_mu_1.8": r.proxy_mu_18,
                  "proxy_T_1cm": r.proxy_T_1cm, "proxy_L08": r.proxy_L08, "priority": r.priority,
                  "flag_Be": "Be" in els, "flag_Tl": "Tl" in els, "flag_As": "As" in els, "flag_Pb": "Pb" in els,
                  "flag_Bi": "Bi" in els, "hygroscopic": hygroscopic(r.formula), "melting_point_unknown": True,
                  "interp_only": bool(r.interp_only) if r.interp_only == r.interp_only and r.interp_only is not None else False,
                  "unscreenable": bool(r.unscreenable) if r.unscreenable == r.unscreenable and r.unscreenable is not None else False,
                  "unstable_dft": r.funnel_step == "unstable_dft",
                  "not_computed_budget": r.funnel_step == "not_computed_budget",
                  "supercell_too_large": r.funnel_step == "supercell_too_large",
                  "mlip_failed": r.funnel_step == "mlip_failed",
                  "not_run_memory": r.funnel_step == "not_run_memory",
                  "tierA_phase_mismatch": bool(r.tierA_phase_mismatch) if r.tierA_phase_mismatch == r.tierA_phase_mismatch else False,
                  "not_off_bragg_single_crystal": bool(r.not_off_bragg_single_crystal) if r.not_off_bragg_single_crystal == r.not_off_bragg_single_crystal else False,
                  "zmax": zmax, "nucdata_incomplete": any(NUC_INCOMPLETE.get(e, False) for e in els)}
        if has:
            mean_ph = float(f["mean_phonon_meV"])
            soft = mean_ph < SOFT_MEV
            mlip_unst = (f["vdos_source"] == "sevennet") and (not bool(f["stable"]) or
                                                               (f["min_freq_comm_THz"] is not None and f["min_freq_comm_THz"] < FC_UNSTABLE))
            ml_interp = (f["vdos_source"] == "sevennet") and not mlip_unst and f["min_freq_mesh_THz"] is not None \
                and f["min_freq_mesh_THz"] < FC_UNSTABLE
            common.update({"mean_phonon_meV": mean_ph, "soft_lattice": soft, "heavy_soft": soft and zmax >= Z_HEAVY,
                           "untested_heavy": (not soft) and zmax >= Z_HEAVY,
                           "molecular": is_molecular(r.formula, f["dim"] if isinstance(f["dim"], str) else ""),
                           "mlip_unstable": bool(mlip_unst), "numberdensity_A3": float(f["numberdensity_A3"]),
                           "theta_D_elastic": f.get("theta_D_elastic") if f["vdos_source"] == "mp_dft"
                           else r.theta_D_debye_table})
            if ml_interp:
                common["interp_only"] = True
        sig = rj.sigma_tot(r.formula)
        for tag in TAGS:
            row = dict(common, T=293.15 if tag == "293K" else 77.0)
            if has:
                n = float(f["numberdensity_A3"])
                L = float(f[f"L08_{tag}"])
                s = SIG_DFT if f["vdos_source"] == "mp_dft" else float(np.hypot(SIG_DFT, SIG_MLIP_ADD[tag]))
                if common["soft_lattice"]:
                    s = max(s, SIG_SOFT)
                Ls = np.array([L, L * np.exp(-s), L * np.exp(s)])
                a1 = rj.a_fast(sig, n, Ls, "1E")
                aw = rj.a_fast(sig, n, Ls, "watt")
                g1, g6 = rj.mu_gamma(r.formula, n)
                row.update({"L08": L, "L08_hot": float(f[f"L08_hot_{tag}"]), "sigma_model": s,
                            "mu_1.0": float(f[f"mu_1.0_{tag}"]), "mu_1.8": float(f[f"mu_1.8_{tag}"]),
                            "L08_proxy_vdos_n": float(f[f"L08_proxy_{tag}"]),
                            "A_fast_1E": a1[0], "A_fast_1E_lo": a1[1], "A_fast_1E_hi": a1[2],
                            "A_fast_watt": aw[0], "A_fast_watt_lo": aw[1], "A_fast_watt_hi": aw[2],
                            "A_gamma_1MeV": g1 * L, "A_gamma_1MeV_lo": g1 * Ls[1], "A_gamma_1MeV_hi": g1 * Ls[2],
                            "A_gamma_6MeV": g6 * L, "A_gamma_6MeV_lo": g6 * Ls[1], "A_gamma_6MeV_hi": g6 * Ls[2],
                            "mu_epi_1eV": rj.mu_epi(r.formula, n),
                            "nieman_ratio": rj.mu_epi(r.formula, n) / float(f[f"mu_1.8_{tag}"])})
                dl = f.get(f"debye_L08_{tag}")
                if dl is None or dl != dl:   # SevenNet records carry no Debye value; take the Debye table's
                    dl = getattr(r, f"debye_L08_{tag}")
                row["debye_L08"] = float(dl) if dl is not None and dl == dl else np.nan
            else:
                row["debye_L08"] = getattr(r, f"debye_L08_{tag}")
                row["theta_D_elastic"] = r.theta_D_debye_table
            rows.append(row)
    at = pd.DataFrame(rows)
    for c in FLAG_COLS:                       # flags absent for compounds without a full model -> False
        at[c] = at[c].fillna(False).astype(bool)
    at["outside_domain"] = at.soft_lattice | at.heavy_soft | at.molecular
    at["main_front_eligible"] = at.vdos_source.ne("none") & ~at.outside_domain
    at = at.sort_values(["material_id", "T"], ascending=[True, False]).reset_index(drop=True)
    return at, full


def write(at, full):
    at.to_parquet(F("atlas.parquet"))
    # LAM_WIDE curves for the Tier-A full-model compounds (band / moderator sensitivity runs)
    tA = set(at[(at.tier == "A") & at.vdos_source.ne("none")].material_id)
    cur = full[full.material_id.isin(tA)][["material_id", "vdos_source", "curve_full_293K", "curve_full_77K",
                                           "curve_proxy_293K", "curve_proxy_77K"]]
    cur.to_parquet(F("stage3_curves.parquet"))
    one = at[at["T"] > 100]
    counts = {"rows": len(at), "compounds": int(one.material_id.nunique()),
              "by_funnel_step": one.funnel_step.value_counts().to_dict(),
              "by_vdos_source": one.vdos_source.value_counts().to_dict(),
              "tierA_by_vdos_source": one[one.tier == "A"].vdos_source.value_counts().to_dict(),
              "tierA_by_funnel_step": one[one.tier == "A"].funnel_step.value_counts().to_dict(),
              "full_model_flags": {k: int(one[one.vdos_source.ne("none")][k].fillna(False).sum()) for k in
                                   ("soft_lattice", "heavy_soft", "untested_heavy", "molecular", "mlip_unstable",
                                    "interp_only", "unscreenable", "hygroscopic")},
              "main_front_eligible": {"A": int(one[(one.tier == "A")].main_front_eligible.sum()),
                                      "B": int(one[(one.tier == "B")].main_front_eligible.sum())},
              "with_H_in_atlas": int(one.formula.map(lambda x: "H" in {e.symbol for e in Composition(x).elements}).sum())}
    old = json.loads((F("atlas_counts.json")).read_text())
    old["atlas"] = counts
    (F("atlas_counts.json")).write_text(json.dumps(old, indent=1, default=str))
    print(json.dumps(counts, indent=1, default=str))


if __name__ == "__main__":
    main()
