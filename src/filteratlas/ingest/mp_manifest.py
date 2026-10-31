"""SHA-256 manifest of the frozen Materials Project inputs (summary table, phonon DOS, force constants) and the
EXFOR entries, so a rerun can prove it used the same data. Output: cache_manifest.json
"""
import hashlib
import json
import sys
from pathlib import Path

from filteratlas.common.paths import D, F  # noqa: E402


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    files = sorted([*(D("mp")).glob("*.parquet"), *(D("mp")).glob("*.json"), *(D("mp")).glob("*.jsonl.gz"),
                    *(D("mp") / "phonon_dos").glob("*/*.json.gz"), *(D("mp") / "force_constants").glob("*.npz"),
                    *(D("exfor")).glob("*.x4")])
    up = D("mp").parent   # data/upstream (manifest keys "mp/...", "exfor/...")
    manifest = {str(p.relative_to(up)).replace("\\", "/"): sha(p) for p in files}
    combined = hashlib.sha256("".join(f"{k}:{v}\n" for k, v in manifest.items()).encode()).hexdigest()
    meta = json.loads((D("mp") / "summary_meta.json").read_text())
    (F("cache_manifest.json")).write_text(json.dumps({"mp_db_version": meta.get("db_version"), "n_files": len(manifest),
                                                             "combined_sha256": combined, "files": manifest}, indent=1))
    print(len(manifest), "files; combined", combined[:16])


if __name__ == "__main__":
    main()
