"""Rank probabilities of the shortlist, and the test of whether the growable front contains anything beyond the known
filters (K) and whether temperature reorders it (R).

Monte Carlo over sigma_model: in every draw each compound's ln L0.8 moves by sigma_i z_i (z_i ~ N(0,1), independent
between compounds, the same z_i at 293 and 77 K). A_gamma moves as exp(sigma z); A_fast is interpolated from a grid
of thicknesses L0.8 exp(sigma u), u in [-5, 5], computed with the atlas rejection code. Fronts are recomputed in every
draw on the full main-front-eligible set of the tier.
Outputs per temperature and tier: P(on the main front), P(best on each axis); for every shortlisted compound
P(X beats sapphire) on each axis; the probabilities P_K(T) and P_R of that test; sensitivity runs (independent
z at the two temperatures, a common component between compounds, sigma construction D). Added after inspection of the
results: "sapphire_offset" shifts sapphire's ln L0.8 by its measured band residual (stage2_uncertainty_model.json ->
band_per_crystal_DFT/Al2O3; the model cross section lies above the sapphire data, so the model L0.8 is too short) in
every draw, the rest as in the main run. Run it alone with "python -m filteratlas.shortlist.rank_mc sapphire_offset";
it is merged into the existing output.
-> derived/stage4_rank_probabilities.json
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from filteratlas.atlas.assemble import FLAG_COLS, SIG_MLIP_ADD, SIG_SOFT
from filteratlas.common.config import cfg
from filteratlas.common.paths import F
from filteratlas.fom import rejection as rj

C = cfg("shortlist")["rank_mc"]
U = cfg("uncertainty")
KNOWN = set(U["known_filters"])
SIG_D = U["sigma_construction_D"]
TAGS = {"293K": 293.15, "77K": 77.0}
AXES = {"main": ("A_fast_1E", "A_gamma_1MeV"), "secondary": ("A_fast_watt", "A_gamma_6MeV")}
CHUNK = 1000


def label(r):
    return f"{r.formula} {r.sg_symbol} [{r.material_id[-4:]}]"


def load():
    at = pd.read_parquet(F("atlas.parquet"))
    for c in FLAG_COLS:
        at[c] = at[c].fillna(False).astype(bool)
    return at


def members(at, sl):
    """Compounds in the draws: the main-front-eligible sets of both tiers (at either temperature) and the shortlist."""
    elig = at[at.main_front_eligible & ~at.unstable_dft & at.vdos_source.ne("none")]
    ids = sorted(set(elig.material_id) | set(sl.material_id))
    return ids


def grids(at, ids):
    """Per compound and temperature: L0.8, sigma (B and D), and A values on the u grid [n_ids, n_grid]."""
    u = np.linspace(-C["grid_sigmas"], C["grid_sigmas"], C["n_grid"])
    out = {}
    for tag, T in TAGS.items():
        d = at[at["T"] == T].set_index("material_id").loc[ids]
        n = len(ids)
        g = {k: np.full((n, len(u)), np.nan) for k in ("A_fast_1E", "A_fast_watt")}
        sB, sD = d.sigma_model.values.astype(float), np.zeros(n)
        g1 = (d.A_gamma_1MeV / d.L08).values
        g6 = (d.A_gamma_6MeV / d.L08).values
        for i, (mid, r) in enumerate(d.iterrows()):
            s = sD[i] = max(SIG_D if r.vdos_source == "mp_dft" else float(np.hypot(SIG_D, SIG_MLIP_ADD[tag])),
                            SIG_SOFT if r.soft_lattice else 0.0)
            sig = rj.sigma_tot(r.formula)
            L = r.L08 * np.exp(sB[i] * u)
            for spec, key in (("1E", "A_fast_1E"), ("watt", "A_fast_watt")):
                g[key][i] = rj.a_fast(sig, r.numberdensity_A3, L, spec)
            mid_i = len(u) // 2
            if abs(g["A_fast_1E"][i, mid_i] - r.A_fast_1E) > 1e-6 * max(1.0, r.A_fast_1E):
                raise RuntimeError(f"{mid}: grid centre {g['A_fast_1E'][i, mid_i]} != atlas {r.A_fast_1E}")
        out[tag] = {"L08": d.L08.values, "sB": sB, "sD": sD, "g1": g1, "g6": g6, **g,
                    "tier": d.tier.values, "elig": (d.main_front_eligible & ~d.unstable_dft).values,
                    "known": d.formula.isin(KNOWN).values,
                    "mismatch": d.tierA_phase_mismatch.fillna(False).astype(bool).values}
    return u, out


def values(gr, u, z, sigma, sB):
    """A values for draws z [n_draw, n_ids] with the given sigma; A_fast interpolated on the sigma_B grid."""
    s = z * sigma / sB                      # position on the u grid (the grid was built with sigma_B)
    pos = (np.clip(s, u[0], u[-1]) - u[0]) / (u[1] - u[0])
    i0 = np.clip(np.floor(pos).astype(int), 0, len(u) - 2)
    w = pos - i0
    out = {}
    for k in ("A_fast_1E", "A_fast_watt"):
        G = gr[k]
        rows = np.arange(G.shape[0])[None, :]
        out[k] = G[rows, i0] * (1 - w) + G[rows, i0 + 1] * w
    Lz = gr["L08"][None, :] * np.exp(z * sigma)
    out["A_gamma_1MeV"] = gr["g1"][None, :] * Lz
    out["A_gamma_6MeV"] = gr["g6"][None, :] * Lz
    return out


def front_mask(X, Y):
    """Per draw (rows) non-dominated columns, maximising both."""
    o = np.argsort(-X, axis=1, kind="stable")
    ys = np.take_along_axis(Y, o, axis=1)
    prev = np.maximum.accumulate(np.concatenate([np.full((ys.shape[0], 1), -np.inf), ys[:, :-1]], axis=1), axis=1)
    keep_sorted = ys > prev
    m = np.zeros_like(keep_sorted)
    np.put_along_axis(m, o, keep_sorted, axis=1)
    return m


def sapphire_offset():
    """Measured band residual ln(sigma_model / sigma_data) of sapphire (in-scope single-crystal data, DFT VDOS)."""
    return float(json.loads(F("stage2_uncertainty_model.json").read_text())["band_per_crystal_DFT"]["Al2O3"])


def run(at, ids, u, gr, sl, scenario):
    """One Monte Carlo scenario; returns probabilities and the K / R tallies."""
    rng = np.random.default_rng(C["seed"])
    nid = len(ids)
    pos = {m: i for i, m in enumerate(ids)}
    sapph = pos[at[(at.formula == "Al2O3") & (at.sg_symbol == "R-3c") & (at.tier == "A")].material_id.iloc[0]]
    sl_idx = np.array([pos[m] for m in sl.material_id])
    tiers = {t: {tier: np.where(gr[t]["elig"] & (gr[t]["tier"] == tier))[0] for tier in ("A", "B")} for t in TAGS}
    if scenario == "no_phase_mismatch":   # Tier-A entries whose database phase is not the grown one left out
        tiers = {t: {"A": v["A"][~gr[t]["mismatch"][v["A"]]], "B": v["B"]} for t, v in tiers.items()}
    both_A = np.intersect1d(tiers["293K"]["A"], tiers["77K"]["A"])
    acc = {t: {"front": {ax: {tier: np.zeros(nid) for tier in ("A", "B")} for ax in AXES},
               "best": {k: {tier: np.zeros(nid) for tier in ("A", "B")} for k in sum(AXES.values(), ())},
               "beats": {k: np.zeros(len(sl_idx)) for k in sum(AXES.values(), ())},
               "K": 0} for t in TAGS}
    taus, reorder = [], 0
    n_done = 0
    while n_done < C["n_draws"]:
        nd = min(CHUNK, C["n_draws"] - n_done)
        z0 = rng.standard_normal((nd, nid))
        if scenario == "common":
            c = C["common_component"]
            z0 = np.sqrt(c) * rng.standard_normal((nd, 1)) + np.sqrt(1 - c) * z0
        if scenario == "sapphire_offset":   # ln L0.8 of sapphire + its measured residual, in units of its sigma
            z0 = z0.copy()
            z0[:, sapph] += sapphire_offset() / gr["293K"]["sB"][sapph]
        z = {"293K": z0, "77K": rng.standard_normal((nd, nid)) if scenario == "independent_T" else z0}
        A_by_T, frontA = {}, {}
        for t in TAGS:
            sig = gr[t]["sD"] if scenario == "sigma_D" else gr[t]["sB"]
            A = values(gr[t], u, z[t], sig[None, :], gr[t]["sB"][None, :])
            A_by_T[t] = A
            for ax, (kx, ky) in AXES.items():
                for tier, idx in tiers[t].items():
                    fm = front_mask(A[kx][:, idx], A[ky][:, idx])
                    acc[t]["front"][ax][tier][idx] += fm.sum(0)
                    if ax == "main" and tier == "A":
                        frontA[t] = fm
                        acc[t]["K"] += int(np.all(~fm | gr[t]["known"][idx][None, :], axis=1).sum())
            for k in sum(AXES.values(), ()):
                for tier, idx in tiers[t].items():
                    b = idx[np.argmax(A[k][:, idx], axis=1)]
                    np.add.at(acc[t]["best"][k][tier], b, 1)
                acc[t]["beats"][k] += (A[k][:, sl_idx] > A[k][:, [sapph]]).sum(0)
        # R: tau 293 vs 77 K on A_fast,1/E over the Tier-A set eligible at both T, and the count of
        # Pareto-status changes
        ia = {t: {m: j for j, m in enumerate(tiers[t]["A"])} for t in TAGS}
        c293 = np.array([ia["293K"][m] for m in both_A])
        c77 = np.array([ia["77K"][m] for m in both_A])
        for r in range(nd):
            tau = kendalltau(A_by_T["293K"]["A_fast_1E"][r, both_A], A_by_T["77K"]["A_fast_1E"][r, both_A]).statistic
            ch = int((frontA["293K"][r, c293] != frontA["77K"][r, c77]).sum())
            taus.append(tau)
            reorder += int(ch >= 1 or tau < U["kill3_tau_threshold"])
        n_done += nd
    N = C["n_draws"]
    names = {m: label(r) for m, r in at[at["T"] > 100].set_index("material_id", drop=False).loc[ids].iterrows()}
    res = {"n_draws": N}
    for t in TAGS:
        res[t] = {"p_front": {}, "p_best": {}, "p_beats_sapphire": {}, "P_K": round(acc[t]["K"] / N, 4)}
        for ax in AXES:
            for tier, idx in tiers[t].items():
                p = acc[t]["front"][ax][tier][idx] / N
                res[t]["p_front"][f"{ax}_{tier}"] = {names[ids[i]]: round(float(q), 4)
                                                     for q, i in sorted(zip(p, idx), reverse=True) if q >= 0.001}
        for k in sum(AXES.values(), ()):
            for tier, idx in tiers[t].items():
                p = acc[t]["best"][k][tier][idx] / N
                res[t]["p_best"][f"{k}_{tier}"] = {names[ids[i]]: round(float(q), 4)
                                                   for q, i in sorted(zip(p, idx), reverse=True) if q >= 0.001}
            res[t]["p_beats_sapphire"][k] = {names[ids[i]]: round(float(q / N), 4)
                                             for q, i in zip(acc[t]["beats"][k], sl_idx) if i != sapph}
    taus = np.array(taus)
    res["kill3"] = {"P_K": {t: res[t]["P_K"] for t in TAGS}, "P_R": round(reorder / N, 4),
                    "tau_median": round(float(np.median(taus)), 4),
                    "tau_p2.5_p97.5": [round(float(np.percentile(taus, 2.5)), 4), round(float(np.percentile(taus, 97.5)), 4)],
                    "n_tierA_both_T": len(both_A)}
    k3 = res["kill3"]
    k3["K_holds"] = all(v >= C["kill3_p"] for v in k3["P_K"].values())
    k3["R_holds"] = k3["P_R"] >= C["kill3_p"]
    k3["fires"] = k3["K_holds"] and not k3["R_holds"]
    return res


def statements(res, sl_names, p):
    """'X beats sapphire' (P > p) and 'sapphire beats X' (P < 1 - p), per temperature and axis."""
    out = []
    for t in TAGS:
        for k, d in res[t]["p_beats_sapphire"].items():
            for name, q in d.items():
                if q > p:
                    out.append({"T": t, "axis": k, "compound": name, "P": q, "statement": "beats sapphire"})
                elif q < 1 - p:
                    out.append({"T": t, "axis": k, "compound": name, "P": q, "statement": "sapphire beats it"})
    return out


def main(only=None):
    at = load()
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    ids = members(at, sl)
    u, gr = grids(at, ids)
    out = {"settings": {k: C[k] for k in ("n_draws", "seed", "n_grid", "grid_sigmas", "p_statement",
                                          "common_component", "kill3_p")},
           "n_compounds": len(ids)}
    p = C["p_statement"]
    sl_names = [label(r) for r in sl.merge(at[at["T"] > 100][["material_id"]], on="material_id").itertuples()]
    if only:   # merge one scenario into the existing output (the others are deterministic and unchanged)
        out = json.loads(F("stage4_rank_probabilities.json").read_text())
    for sc in (only or ("main", "independent_T", "common", "sigma_D", "sapphire_offset")):
        r = run(at, ids, u, gr, sl, sc)
        r["statements"] = statements(r, sl_names, p)
        if sc == "sapphire_offset":
            r["offset_lnL08"] = sapphire_offset()
        out[sc] = r
        print(sc, "kill3", r["kill3"], "statements", len(r["statements"]), flush=True)
    base = {(s["T"], s["axis"], s["compound"], s["statement"]) for s in out["main"]["statements"]}
    out["sensitivity_changes"] = {}
    for sc in [s for s in ("independent_T", "common", "sigma_D", "sapphire_offset", "no_phase_mismatch") if s in out]:
        other = {(s["T"], s["axis"], s["compound"], s["statement"]) for s in out[sc]["statements"]}
        out["sensitivity_changes"][sc] = {"lost": sorted(map(list, base - other)), "gained": sorted(map(list, other - base))}
    F("stage4_rank_probabilities.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] or None)
