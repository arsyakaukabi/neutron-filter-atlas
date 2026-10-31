"""scripts/check_data.py: passes on a complete data tree; reports missing files on an empty one."""
import csv
import importlib.util
import shutil

import pytest

from conftest import DATA, REPO

spec = importlib.util.spec_from_file_location("check_data", REPO / "scripts" / "check_data.py")
check_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_data)


def test_check_data_complete_tree():
    if not (DATA / "MANIFEST.csv").exists() or not (DATA / "derived").exists():
        pytest.skip("needs the distributed data under data/ (see data/README.md)")
    missing, bad, extra, n = check_data.check(DATA, fast=True)
    assert n > 0 and not missing and not bad, (missing[:5], bad[:5])


def test_check_data_reports_missing(tmp_path):
    """On a tree with only the curated files, every manifest file outside curated/ is reported missing."""
    manifest = DATA / "MANIFEST.csv"
    if not manifest.exists():
        pytest.skip("data/MANIFEST.csv not present")
    root = tmp_path / "data"
    shutil.copytree(DATA / "curated", root / "curated")
    missing, bad, extra, n = check_data.check(root, fast=True, manifest=manifest)
    rows = list(csv.DictReader(open(manifest, encoding="utf-8")))
    assert len(missing) == sum(r["class"] != "curated" for r in rows)
    assert not bad
