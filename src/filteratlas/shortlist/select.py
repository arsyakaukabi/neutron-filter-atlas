"""The shortlist: who is shortlisted, and why.

Tier A: main-front members and compounds within 1 sigma of the Tier-A main front at 293 K or 77 K; the reference
filters; YAG and BGO; every SiC entry with a full model. Tier B: the union of the Tier-B main-front members at both
temperatures, without heavy_soft, molecular and mlip_unstable compounds (untested_heavy kept, flagged).
-> derived/stage4_shortlist.parquet (one row per compound, with its reasons)
"""
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.paths import F

C = cfg("shortlist")["shortlist"]
FLAGS = ["flag_Be", "flag_Tl", "flag_As", "flag_Pb", "flag_Bi", "untested_heavy", "interp_only", "unscreenable",
         "tierA_phase_mismatch", "mlip_unstable", "hygroscopic"]


def polytype(r):
    """SiC polytype label from space group and sites per cell (2H/4H/6H share P6_3mc)."""
    if r.sg_symbol == "P6_3mc":
        return {4: "2H", 8: "4H", 12: "6H"}.get(r.nsites, f"hexagonal, {r.nsites} sites")
    if r.sg_symbol == "F-43m":
        return "3C"
    if r.sg_symbol == "R3m" and r.nsites == 10:
        return "15R"
    return f"{r.sg_symbol}, {r.nsites} sites"


def pick(one, formula, sg):
    """The atlas entry of a reference compound: the Tier-A one, else the lowest-E_hull one with a full model."""
    d = one[(one.formula == formula) & (one.sg_symbol == sg)]
    a = d[d.tier == "A"]
    d = a if len(a) else d[d.vdos_source != "none"]
    return d.sort_values("e_hull").material_id.iloc[0] if len(d) else None


def main():
    at = pd.read_parquet(F("atlas.parquet"))
    one = at[at["T"] > 100].copy()
    nsites = pd.read_parquet(F("stage0_candidates.parquet"))[["material_id", "nsites"]]
    one = one.merge(nsites, on="material_id", how="left").set_index("material_id", drop=False)
    fm = pd.read_parquet(F("stage3_fronts.parquet"))
    why = {}

    def add(m, reason):
        why.setdefault(m, []).append(reason)

    for T in ("293K", "77K"):
        k = fm[fm.key == f"{T}|A|main|main"]
        for m in k[k.front].material_id:
            add(m, f"Tier-A front {T}")
        for m in k[k.within_1sigma & ~k.front].material_id:
            add(m, f"within 1 sigma of the Tier-A front {T}")
    for formula, sg in C["reference_filters"]:
        m = C["quartz_material_id"] if formula == "SiO2" else pick(one, formula, sg)
        add(m, "reference filter")
    for formula, sg in C["hidden_by_memory_gap"]:
        add(pick(one, formula, sg), "hidden by the v0.3.0 memory gap")
    for m in one[(one.formula == "SiC") & (one.vdos_source != "none")].material_id:
        add(m, "SiC polytype")
    excl = C["tierB_exclude_flags"]
    tierB_dropped = {}
    for T in ("293K", "77K"):
        k = fm[(fm.key == f"{T}|B|main|main") & fm.front]
        for m in k.material_id:
            bad = [f for f in excl if bool(one.loc[m, f])]
            if bad:
                tierB_dropped[m] = bad
            else:
                add(m, f"Tier-B front {T}")
    rows = []
    for m, reasons in why.items():
        r = one.loc[m]
        rows.append({"material_id": m, "formula": r.formula, "sg_symbol": r.sg_symbol, "tier": r.tier,
                     "vdos_source": r.vdos_source, "reasons": "; ".join(dict.fromkeys(reasons)),
                     "polytype": polytype(r) if r.formula == "SiC" else "",
                     "flags": ", ".join(f for f in FLAGS if bool(r[f])),
                     "outside_domain": bool(r.outside_domain), "main_front_eligible": bool(r.main_front_eligible),
                     "known_filter": r.formula in set(cfg("uncertainty")["known_filters"])})
    sl = pd.DataFrame(rows).sort_values(["tier", "formula", "sg_symbol"]).reset_index(drop=True)
    sl.to_parquet(F("stage4_shortlist.parquet"), index=False)
    print(sl[["formula", "sg_symbol", "tier", "vdos_source", "polytype", "reasons", "flags"]].to_string())
    print("Tier-B front members dropped:", {f"{one.loc[m, 'formula']} {one.loc[m, 'sg_symbol']}": v
                                           for m, v in tierB_dropped.items()})


if __name__ == "__main__":
    main()
