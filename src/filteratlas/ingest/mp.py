"""Materials Project access: practical candidate set, structures, phonon DOS (bulk, via the public Delta table).

All MP responses are cached under D("mp") (configs/paths.yaml) so reruns are offline.
"""
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from filteratlas.common.config import cfg
from filteratlas.common.constants import E_HULL_MAX, FLAG_ELEMENTS, MP_KEY, RADIOACTIVE, STRONG_ABSORBERS
from filteratlas.common.paths import D

MPC = D("mp")
SUMMARY_FIELDS = cfg("candidates")["mp_summary_fields"]


def rester():
    from mp_api.client import MPRester
    return MPRester(MP_KEY, mute_progress_bars=True)


def summary_table(refresh=False):
    """Experimental, near-hull summary rows (all elements; filtering happens in candidate_set)."""
    out = MPC / "summary_exp_nearhull.parquet"
    if out.exists() and not refresh:
        return pd.read_parquet(out)
    rows = []
    with rester() as m:
        meta = {"db_version": m.db_version, "fetched": time.strftime("%Y-%m-%d %H:%M")}
        docs = m.materials.summary.search(energy_above_hull=(0, E_HULL_MAX), theoretical=False,
                                          fields=SUMMARY_FIELDS)
        for d in docs:
            dd = d if isinstance(d, dict) else d.model_dump()
            sym = dd.get("symmetry") or {}
            k, g = dd.get("bulk_modulus") or {}, dd.get("shear_modulus") or {}
            ph = dd.get("phonon_IDs") or {}
            rows.append({
                "material_id": str(dd["material_id"]), "formula": dd["formula_pretty"],
                "elements": ",".join(str(e) for e in dd["elements"]), "nelements": dd["nelements"],
                "nsites": dd["nsites"], "density": dd["density"], "volume": dd["volume"],
                "e_hull": dd["energy_above_hull"], "is_stable": dd["is_stable"],
                "theoretical": dd["theoretical"],
                "spacegroup": sym.get("number"), "sg_symbol": sym.get("symbol"),
                "crystal_system": str(sym.get("crystal_system")),
                "K_vrh": k.get("vrh") if isinstance(k, dict) else None,
                "G_vrh": g.get("vrh") if isinstance(g, dict) else None,
                "phonon_pheasy": ",".join(str(x) for x in (ph.get("pheasy") or [])) if isinstance(ph, dict) else "",
                "phonon_dfpt": ",".join(str(x) for x in (ph.get("dfpt") or [])) if isinstance(ph, dict) else "",
            })
    df = pd.DataFrame(rows)
    df.to_parquet(out)
    (MPC / "summary_meta.json").write_text(json.dumps({**meta, "n": len(df)}, indent=1))
    return df


def candidate_set(df=None):
    """Apply the element filters fixed before the computation. Returns (kept, exclusion_counts)."""
    df = summary_table() if df is None else df
    els = df.elements.str.split(",").apply(set)
    reasons = pd.Series("", index=df.index)
    reasons[els.apply(lambda s: "H" in s)] = "H"
    reasons[(reasons == "") & els.apply(lambda s: bool(s & RADIOACTIVE))] = "radioactive"
    reasons[(reasons == "") & els.apply(lambda s: bool(s & STRONG_ABSORBERS))] = "strong_absorber"
    kept = df[reasons == ""].copy()
    kept["flags"] = els[reasons == ""].apply(lambda s: ",".join(sorted(s & FLAG_ELEMENTS)))
    kept["phonon_method"] = np.where(kept.phonon_pheasy != "", "pheasy",
                                     np.where(kept.phonon_dfpt != "", "dfpt", ""))
    kept["phonon_id"] = np.where(kept.phonon_pheasy != "", kept.phonon_pheasy.str.split(",").str[0],
                                 np.where(kept.phonon_dfpt != "", kept.phonon_dfpt.str.split(",").str[0], ""))
    counts = {"summary_exp_nearhull": len(df), **reasons[reasons != ""].value_counts().to_dict(),
              "kept": len(kept), "kept_with_phonon": int((kept.phonon_id != "").sum()),
              "kept_with_elastic": int(kept.G_vrh.notna().sum())}
    return kept, counts


def structures(material_ids, refresh=False):
    """Structures (pymatgen dicts) for material ids, cached in D("mp")/structures.jsonl.gz."""
    from pymatgen.core import Structure
    path = MPC / "structures.jsonl.gz"
    have = {}
    if path.exists() and not refresh:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                have[r["material_id"]] = r["structure"]
    need = [i for i in material_ids if i not in have]
    if need:
        with rester() as m, gzip.open(path, "at", encoding="utf-8") as fh:
            for i in range(0, len(need), 500):
                for d in m.materials.summary.search(material_ids=need[i:i + 500],
                                                    fields=["material_id", "structure"]):
                    # key by the model_dump() form (AlphaID, as in summary_table); str(d.material_id) can
                    # give the legacy numeric form for the same entry
                    mid = str(d.model_dump()["material_id"])
                    s = d.structure.as_dict()
                    have[mid] = s
                    fh.write(json.dumps({"material_id": mid, "structure": s}) + "\n")
    return {i: Structure.from_dict(have[i]) for i in material_ids if i in have}


def phonon_dos_bulk(id_method_pairs, batch=200, refresh=False):
    """Fetch phonon DOS rows from the public MP Delta table for (phonon_id, method) pairs.

    Cached per identifier in D("mp")/phonon_dos/<method>/<identifier>.json.gz with keys
    frequencies (THz), densities, projected_densities (sites x freq), structure (dict).
    Returns dict {(phonon_id, method): path}; failures are logged to D("mp")/phonon_dos/failures.json.
    """
    from emmet.core.mpid import AlphaID
    root = MPC / "phonon_dos"
    out, todo = {}, {}
    for pid, meth in id_method_pairs:
        p = root / meth / f"{pid}.json.gz"
        if p.exists() and not refresh:
            out[(pid, meth)] = p
        else:
            todo.setdefault(meth, []).append(pid)
    failures = {}
    if todo:
        with rester() as m:
            ph = m.materials.phonon
            lbl, _ = ph._get_delta_table("materialsproject-parsed", "phonon/electronic-structure/dos/",
                                         label="ph_dos")
            for meth, pids in todo.items():
                (root / meth).mkdir(parents=True, exist_ok=True)
                for i in range(0, len(pids), batch):
                    chunk = pids[i:i + batch]
                    key = {str(AlphaID(p.split("-")[-1], padlen=8)): p for p in chunk}
                    ids = ",".join(f"'{k}'" for k in key)
                    q = f"SELECT * FROM {lbl} WHERE identifier IN ({ids}) AND phonon_method='{meth}'"
                    for attempt in range(4):
                        try:
                            rows = ph._query_delta_single(q, label=lbl).to_pylist(maps_as_pydicts="strict")
                            break
                        except Exception as e:  # noqa: BLE001
                            err = repr(e)[:300]
                            time.sleep(10 * (attempt + 1))
                    else:
                        for p in chunk:
                            failures[p] = f"query failed: {err}"
                        continue
                    got = set()
                    for r in rows:
                        pid = key.get(str(r.get("identifier")))
                        dos = r.get("dos")
                        if pid is None or not dos:
                            continue
                        got.add(pid)
                        rec = {"frequencies": dos.get("frequencies"), "densities": dos.get("densities"),
                               "projected_densities": dos.get("projected_densities"),
                               "structure": dos.get("structure")}
                        p = root / meth / f"{pid}.json.gz"
                        with gzip.open(p, "wt", encoding="utf-8") as fh:
                            json.dump(rec, fh)
                        out[(pid, meth)] = p
                    for p in set(chunk) - got:
                        failures[p] = "no row returned"
                    print(f"[phonon_dos] {meth} {i + len(chunk)}/{len(pids)} got {len(got)}", flush=True)
    fpath = root / "failures.json"
    old = json.loads(fpath.read_text()) if fpath.exists() else {}
    old.update(failures)
    fpath.parent.mkdir(parents=True, exist_ok=True)
    fpath.write_text(json.dumps(old, indent=1))
    return out


def load_dos(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)
