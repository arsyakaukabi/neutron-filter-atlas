"""Size and sign of the kernel systematic from NCrystal's incoherent approximation for coherent inelastic
scattering, for sapphire and Si only (MP force constants, the atlas entries).

Quantity: the one-phonon coherent inelastic total cross section, averaged over crystal orientations, computed exactly
(mode sum with interference) and in the incoherent approximation, with identical constants, q-mesh, Debye-Waller
factors and kinematics. The difference Delta sigma(E, T) = sigma_coh,1ph(exact) - sigma_coh,1ph(incoherent approx.)
is the correction to the atlas sigma_tot (the incoherent part sigma_inc and multiphonon terms are unchanged; coherent
multiphonon interference is neglected, stated as a limitation).

Formulas (per primitive cell; phonopy eigenvectors, type-II phase convention u_ld ~ e_d exp[iq.(R_l + r_d)]):
  orientation average of a single-crystal mode term a delta^3(Q - Q_m) delta(eps - eps_m) gives a / (2 k^2 Q_m)
  if the scattering is kinematically allowed (|k - k'| <= Q_m <= k + k').
  coherent:   a = (2pi)^3/v0 * (1/N_q) * C/eps * |sum_d b_d/sqrt(M_d) e^{-W_d} e^{i tau.r_d} (Q.e_d)|^2 * (n + 1/2 +- 1/2)
  incoherent approximation (no momentum conservation, angular average of (Q.e_d)^2 = Q^2 |e_d|^2 / 3):
              sigma = (2pi/k^2) (1/N_q) sum_modes C/eps (n + 1/2 +- 1/2) int_{Q-}^{Q+} Q sum_d b_d^2/M_d e^{-2W_d} Q^2 |e_d|^2/3 dQ
  C = hbar^2/(2 m_n) / (1 amu / m_n) [meV A^2], b in sqrt(barn), M in amu, Q in 1/A, eps in meV -> barn per cell.
  2W_d = Q^2 (hbar^2/(2 M_d)) <coth(eps/2kT)/eps>_{g_d} (isotropic, partial DOS from the same mesh).
Only Q <= Q_C is summed for both; beyond it the comb average equals the incoherent approximation (checked by varying
Q_C). The correction is propagated to L0.8 on the atlas mu(lambda) and to the EXFOR benchmark residuals.
-> stage3_coherent.json
"""
import json
import sys
import time
from pathlib import Path

import NCrystal as NC
import numpy as np
import pandas as pd

from filteratlas.xs import evaluate as xe  # noqa: E402
from filteratlas.common.config import cfg, tup
from filteratlas.common.paths import D, F  # noqa: E402

HB2_2MN = 2.072124      # hbar^2/(2 m_n) in meV A^2
MN_AMU = 1.00866491595
KB = 0.08617333         # meV/K
THZ_MEV = 4.135667696
_CC = cfg("figures_of_merit")["coherent_check"]
CRYSTALS = _CC["crystals_phonon_id"]
MESH = _CC["mesh"]
TEMPS = tup(cfg("figures_of_merit")["temperatures_K"])
LAM_GRID = np.linspace(_CC["lambda_A"]["min"], _CC["lambda_A"]["max"], _CC["lambda_A"]["n"])


def phonons(pid, mesh):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    z = np.load(D("mp") / "force_constants" / f"{pid}.npz", allow_pickle=False)
    unit = PhonopyAtoms(symbols=list(z["symbols"]), cell=z["lattice"], scaled_positions=z["frac"])
    ph = Phonopy(unit, supercell_matrix=z["sc"], primitive_matrix=z["prim"])
    ph.force_constants = z["fc"]
    ph.symmetrize_force_constants()
    ph.run_mesh([mesh] * 3, is_gamma_center=True, with_eigenvectors=True, is_mesh_symmetry=False)
    d = ph.get_mesh_dict()
    p = ph.primitive
    return {"q_frac": d["qpoints"], "freq_meV": d["frequencies"] * THZ_MEV, "eig": d["eigenvectors"],
            "lat": np.asarray(p.cell), "frac": np.asarray(p.scaled_positions), "symbols": list(p.symbols),
            "masses": np.asarray(p.masses)}


def setup(P, T):
    lat = P["lat"]
    rec = 2 * np.pi * np.linalg.inv(lat).T            # rows: reciprocal vectors (1/A)
    v0 = abs(np.linalg.det(lat))
    nat = len(P["symbols"])
    b = np.array([NC.atomDB(s).coherentScatLen() for s in P["symbols"]])     # sqrt(barn)
    M = P["masses"]
    eps = P["freq_meV"]                                  # (Nq, 3 nat)
    ok = eps > 0.5                                       # drop the acoustic Gamma modes / numerical zeros
    e = P["eig"].reshape(len(eps), nat, 3, 3 * nat)      # [q, atom, xyz, mode]
    w2 = np.abs(e) ** 2                                  # |e_d,alpha|^2
    ed2 = w2.sum(axis=2)                                 # [q, atom, mode] |e_d|^2
    beta = 1.0 / (KB * T)
    # partial-DOS-based isotropic mean square displacement factor per atom: 2W = Q^2 * u2[d]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        coth = np.where(ok, 1.0 / np.tanh(0.5 * beta * eps), 0.0)
        inv = np.where(ok, 1.0 / eps, 0.0)
    # isotropic <(Q^.u_d)^2> = (hbar^2/2M_d) (1/N_q) sum_{q,m} (|e_d|^2 / 3) coth(eps/2kT) / eps ; 2W_d = Q^2 * u2[d]
    u2 = np.array([HB2_2MN * MN_AMU / M[d] * np.sum(ed2[:, d, :] / 3.0 * coth * inv) / len(eps) for d in range(nat)])
    nB = np.where(ok, 1.0 / np.expm1(beta * np.where(ok, eps, 1.0)), 0.0)
    return {"rec": rec, "v0": v0, "nat": nat, "b": b, "M": M, "eps": eps, "ok": ok, "e": e, "ed2": ed2, "u2": u2,
            "nB": nB, "pos": P["frac"] @ P["lat"], "q": P["q_frac"] @ rec, "Nq": len(eps)}


def kin(k, eps, sign):
    """k' for emission (sign=+1, neutron loses eps) or absorption (-1); nan where not allowed."""
    kp2 = k ** 2 - sign * eps / HB2_2MN
    return np.sqrt(np.where(kp2 > 0, kp2, np.nan))


def sigma_coh_all(S, lams, qcs):
    """Orientation-averaged exact one-phonon coherent cross section [barn per cell], {qc: array over lams}.
    Mode terms (q, tau, mode) with 0 < |Q| <= qc, weight w = (2pi)^3/v0/N_q * C/eps * |F|^2, accumulated per q-point
    so memory stays small; each term contributes w * occ / (2 k^2 Q) where kinematically allowed."""
    qmax = max(qcs)
    rn = np.linalg.norm(S["rec"], axis=1)
    nmax = int(np.ceil((qmax + rn.max()) / rn.min())) + 1
    rng = np.arange(-nmax, nmax + 1)
    taus = np.array(np.meshgrid(rng, rng, rng, indexing="ij")).reshape(3, -1).T @ S["rec"]
    taus = taus[np.linalg.norm(taus, axis=1) <= qmax + np.linalg.norm(S["q"], axis=1).max()]
    phase = np.exp(1j * taus @ S["pos"].T)               # [tau, atom]
    bm = S["b"] / np.sqrt(S["M"])
    ks = 2 * np.pi / np.asarray(lams)
    acc = {qc: np.zeros(len(ks)) for qc in qcs}
    for iq in range(S["Nq"]):
        Q = S["q"][iq][None, :] + taus
        Qn = np.linalg.norm(Q, axis=1)
        sel = (Qn <= qmax) & (Qn > 1e-9)
        if not sel.any():
            continue
        Q, Qn, ph = Q[sel], Qn[sel], phase[sel]
        dw = np.exp(-0.5 * Qn[:, None] ** 2 * S["u2"][None, :])          # e^{-W_d}
        Qe = np.einsum("tx,dxm->tdm", Q, S["e"][iq])
        F2 = np.abs(np.einsum("td,tdm->tm", bm[None, :] * dw * ph, Qe)) ** 2
        mk = S["ok"][iq]
        eps = S["eps"][iq][mk]
        w = (2 * np.pi) ** 3 / S["v0"] / S["Nq"] * HB2_2MN * MN_AMU / eps[None, :] * F2[:, mk]   # [tau, mode]
        nb = S["nB"][iq][mk]
        for j, k in enumerate(ks):
            for sign, occ in ((1, nb + 1.0), (-1, nb)):
                kp = kin(k, eps, sign)                       # [mode]
                good = ~np.isnan(kp)
                if not good.any():
                    continue
                kpg = kp[good]
                allowed = (Qn[:, None] >= np.abs(k - kpg)[None, :]) & (Qn[:, None] <= (k + kpg)[None, :])
                term = np.where(allowed, w[:, good] * occ[good][None, :] / (2 * k ** 2 * Qn[:, None]), 0.0)
                for qc in qcs:
                    acc[qc][j] += float(term[Qn <= qc].sum())
    return acc


def sigma_inc(S, lams, qc):
    """Incoherent approximation of the same one-phonon coherent term [barn per cell], Q <= qc, closed-form Q integral:
    int Q^3 e^{-a Q^2} dQ = -(1 + a Q^2) e^{-a Q^2} / (2 a^2)."""
    ok = S["ok"]
    eps = np.where(ok, S["eps"], 1.0)
    out = []
    for lam in lams:
        k = 2 * np.pi / lam
        tot = 0.0
        for sign, occ in ((1, S["nB"] + 1.0), (-1, S["nB"])):
            kp = kin(k, eps, sign)
            good = ok & ~np.isnan(kp)
            lo = np.where(good, np.abs(k - np.nan_to_num(kp)), 0.0)
            hi = np.where(good, np.minimum(k + np.nan_to_num(kp), qc), 0.0)
            good &= hi > lo
            integ = np.zeros_like(eps)
            for d in range(S["nat"]):
                a = S["u2"][d]
                F = lambda x: -(1 + a * x ** 2) * np.exp(-a * x ** 2) / (2 * a ** 2)  # noqa: E731
                integ += S["b"][d] ** 2 / S["M"][d] * S["ed2"][:, d, :] * (F(hi) - F(lo)) / 3.0
            tot += float(np.sum(np.where(good, 2 * np.pi / k ** 2 / S["Nq"] * HB2_2MN * MN_AMU / eps * occ * integ, 0.0)))
        out.append(tot)
    return np.array(out)


def main(qcs=None):
    qcs = qcs or [float(x) for x in sys.argv[1:]] or list(_CC["q_cut_per_A"])
    at = pd.read_parquet(F("atlas.parquet"))
    cur = pd.read_parquet(F("stage3_curves.parquet")).set_index("material_id")
    out = {"method": __doc__.split("\n\n")[1], "q_cut_1_per_A": qcs, "lambda_A": LAM_GRID.tolist(), "crystals": {}}
    for f, pid in CRYSTALS.items():
        t0 = time.time()
        P = phonons(pid, MESH[f])
        res = {"phonon_id": pid, "mesh": MESH[f], "n_atoms_cell": len(P["symbols"])}
        row = at[(at.formula == f) & (at.tier == "A")].drop_duplicates("material_id").iloc[0]
        n = float(row.numberdensity_A3)
        for T in TEMPS:
            S = setup(P, T)
            tag = "293K" if T > 100 else "77K"
            per = {}
            cohs = sigma_coh_all(S, LAM_GRID, qcs)
            for qc in qcs:
                c = cohs[qc] / S["nat"]
                i = sigma_inc(S, LAM_GRID, qc) / S["nat"]
                per[str(qc)] = {"coh_exact_b": c.tolist(), "coh_incapprox_b": i.tolist(), "delta_b": (c - i).tolist()}
            d = np.array(per[str(qcs[-1])]["delta_b"])
            mu = np.interp(xe.LAM, xe.LAM_WIDE, np.asarray(cur.loc[row.material_id, f"curve_full_{tag}"]))
            mu_c = mu + n * np.interp(xe.LAM, LAM_GRID, d)
            L0, L1 = xe.l_target(mu), xe.l_target(mu_c)
            per["dlnL08"] = float(np.log(L1 / L0))
            per["L08_atlas"], per["L08_corrected"] = L0, L1
            per["band_mean_rel_change_mu"] = float(np.trapezoid(xe.maxwellian_flux(xe.LAM) * (mu_c / mu - 1), xe.LAM)
                                                   / np.trapezoid(xe.maxwellian_flux(xe.LAM), xe.LAM))
            res[tag] = per
        res["wall_s"] = round(time.time() - t0, 1)
        out["crystals"][f] = res
        print(f, {t: (res[t]["dlnL08"], res[t]["band_mean_rel_change_mu"]) for t in ("293K", "77K")}, res["wall_s"], flush=True)
    (F("stage3_coherent.json")).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
