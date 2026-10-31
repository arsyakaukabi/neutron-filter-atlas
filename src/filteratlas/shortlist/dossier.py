"""Practicality dossier for the Tier-A shortlist.

Merges the sourced facts of data/curated/stage4_dossier.csv (melting point, growth method, largest stated crystal,
water/hygroscopicity, price, measurement search), the toxicity flags of the atlas (Be, Tl, As, Pb, Bi), and the
activation run (intermediate/stage4_activation/results.json; nuclides below 1 % of the total activity are dropped as
numerical noise). Price per cm3 is computed from each catalogue price and its stated product size and given as a range over the items;
these are thin substrates, so the figure is not the price of a filter-sized block.
-> derived/stage4_dossier.json
"""
import json

import pandas as pd

from filteratlas.common.paths import D, F

TOX = {"flag_Be": "Be: toxic dust when machined", "flag_Tl": "Tl: toxic", "flag_As": "As: toxic",
       "flag_Pb": "Pb: toxic", "flag_Bi": "Bi: activates to Po-210 (alpha, 138 d) under neutrons"}


def volume_cm3(size):
    if not size:
        return None
    if size.startswith("disc:"):
        d, t = (float(x) for x in size[5:].split("x"))
        return 3.141592653589793 * (d / 20) ** 2 * (t / 10)
    a, b, c = (float(x) for x in size.split("x"))
    return a * b * c / 1000.0


def main():
    sl = pd.read_parquet(F("stage4_shortlist.parquet"))
    sl = sl[sl.tier == "A"]
    at = pd.read_parquet(F("atlas.parquet"))
    one = at[at["T"] > 100].set_index("material_id")
    cur = pd.read_csv(F("stage4_dossier.csv"), dtype=str).fillna("")
    act = json.loads((D("stage4_activation") / "results.json").read_text())
    act_by = {c["material_id"]: c for c in act["compounds"]}
    out = {"activation_settings": {k: act[k] for k in ("openmc", "chain_sha256", "flux_n_cm2_s", "irradiation_days",
                                                        "cooling_days", "volume_cm3")}, "compounds": {}}
    for r in sl.itertuples():
        key = f"{r.formula} {r.sg_symbol}"
        rows = cur[cur.compound == key]
        rec = {"material_id": r.material_id}
        prices = []
        for q in rows.itertuples():
            e = {"value": q.value, "source": q.source_url, "access_date": q.access_date}
            if q.note:
                e["note"] = q.note
            if q.quantity == "price":
                v = volume_cm3(q.size_or_unit)
                try:
                    amount, cur_ = q.value.split()[0], q.value.split()[1]
                    e["per_cm3"] = float(amount) / v if v else None
                    e["currency"] = cur_
                except (ValueError, IndexError):
                    pass
                prices.append(e)
                continue
            rec[q.quantity] = e
        if prices:   # every catalogue item, and the range per cm3 (all thin substrates)
            per = [(x["per_cm3"], x["currency"]) for x in prices if x.get("per_cm3")]
            cur_ = {c for _, c in per}
            rng = None
            if per and len(cur_) == 1:
                lo, hi = min(p_ for p_, _ in per), max(p_ for p_, _ in per)
                c_ = cur_.pop()
                rng = (f"{lo:.0f} {c_}/cm3 (thin substrate)" if len(per) == 1 or round(lo) == round(hi)
                       else f"{lo:.0f}-{hi:.0f} {c_}/cm3 ({len(per)} thin-substrate items)")
            rec["price"] = {"value": prices[0]["value"] if not per else "; ".join(x["value"] for x in prices),
                            "items": prices, "price_per_cm3": rng}
        rec["toxicity"] = [TOX[f] for f in TOX if bool(one.loc[r.material_id, f])] or ["none flagged"]
        a = act_by.get(r.material_id)
        if a:
            rec["activation"] = {}
            for t, v in a["at"].items():
                tot = v["total_Bq"]
                rec["activation"][t] = {"total_Bq": float(f"{tot:.3g}"),
                                        "dominant": [{"nuclide": x["nuclide"], "share": round(x["Bq"] / tot, 3),
                                                      "half_life_h": round(x["half_life_s"] / 3600, 3) if x["half_life_s"] else None}
                                                     for x in v["top3"] if tot > 0 and x["Bq"] / tot >= 0.01]}
        out["compounds"][f"{key} [{r.material_id[-4:]}]"] = rec
    F("stage4_dossier.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(len(out["compounds"]), "compounds")


if __name__ == "__main__":
    main()
