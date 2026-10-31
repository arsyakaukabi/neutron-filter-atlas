"""Verify data/ against data/MANIFEST.csv: missing, extra and mismatched files (size and SHA-256).

Exit code 0 when every manifest file is present with the right checksum (extra files are reported, not failed);
1 otherwise. --fast compares sizes only. --class upstream,derived restricts the check (e.g. a collaborator who only
reads the atlas checks --class derived).
Usage: python scripts/check_data.py [--fast] [--class derived]
"""
import argparse
import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from filteratlas.common.hashing import sha256  # noqa: E402


def check(data_root=REPO / "data", classes=None, fast=False, manifest=None):
    manifest = manifest or data_root / "MANIFEST.csv"
    rows = list(csv.DictReader(open(manifest, encoding="utf-8")))
    if classes:
        rows = [r for r in rows if r["class"] in classes]
    want = {r["path"]: r for r in rows}
    missing, bad = [], []
    for path, r in want.items():
        p = data_root.parent / path
        if not p.exists():
            missing.append(path)
        elif p.stat().st_size != int(r["bytes"]) or (not fast and sha256(p) != r["sha256"]):
            bad.append(path)
    have = set()
    for cls in (classes or ("upstream", "intermediate", "derived", "curated")):
        have |= {f"data/{p.relative_to(data_root).as_posix()}" for p in (data_root / cls).rglob("*") if p.is_file()}
    extra = sorted(have - set(want) - {"data/README.md", "data/MANIFEST.csv"})
    return missing, bad, extra, len(want)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--class", dest="classes", default="")
    a = ap.parse_args()
    classes = [c for c in a.classes.split(",") if c] or None
    missing, bad, extra, n = check(classes=classes, fast=a.fast)
    print(f"manifest files: {n}; missing: {len(missing)}; mismatched: {len(bad)}; extra (not in manifest): {len(extra)}")
    for lab, lst in (("MISSING", missing), ("MISMATCH", bad), ("EXTRA", extra)):
        for x in lst[:50]:
            print(f"  {lab} {x}")
        if len(lst) > 50:
            print(f"  ... {len(lst) - 50} more")
    if bad:
        print("Mismatched upstream files usually mean a re-fetch from the current Materials Project release (MP serves "
              "only its current release); derived results built from them can differ from the published atlas.")
    sys.exit(0 if not missing and not bad else 1)


if __name__ == "__main__":
    main()
