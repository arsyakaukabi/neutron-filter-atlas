"""Practical candidate set + bulk MP phonon DOS for every candidate that has one (steps 00 and 01).

Writes stage0_candidates.parquet and stage0_mp_meta.json.
"""
import json
import sys
import time
from pathlib import Path

from filteratlas.ingest import mp as mq  # noqa: E402
from filteratlas.common.paths import D, F  # noqa: E402

def main():
    t0 = time.time()
    df = mq.summary_table()
    kept, counts = mq.candidate_set(df)
    kept.to_parquet(F("stage0_candidates.parquet"))
    pairs = [(r.phonon_id, r.phonon_method) for r in kept.itertuples() if r.phonon_id]
    got = mq.phonon_dos_bulk(pairs)
    fails = json.loads((D("mp") / "phonon_dos" / "failures.json").read_text())
    meta = {**json.loads((D("mp") / "summary_meta.json").read_text()), "filter_counts": counts,
            "phonon_requested": len(pairs), "phonon_fetched": len(got),
            "phonon_fetch_failures": len([p for p, _ in pairs if p in fails]),
            "wall_s": round(time.time() - t0)}
    (F("stage0_mp_meta.json")).write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
