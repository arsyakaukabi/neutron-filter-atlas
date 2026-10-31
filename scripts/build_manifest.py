"""Regenerate data/MANIFEST.csv from the data tree (run before distributing the data).

One row per file under data/upstream, data/intermediate, data/derived and data/curated:
    path, bytes, sha256, class, produced_by, source, notes
- produced_by: the pipeline step that writes the file, or "upstream" / "curated".
- source: for upstream files, the database and version or retrieval date; for computed files, the git commit given
  with --commit, or "computed by this pipeline". A file whose SHA-256 equals its row in the previous manifest
  (--previous) keeps that row's source.
- notes: free text (empty unless set by hand).
make check-data reads only path, bytes, sha256 and class.
Usage: python scripts/build_manifest.py [--commit <sha>] [--previous <old MANIFEST.csv>]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from filteratlas.common.hashing import sha256  # noqa: E402

DATA = REPO / "data"
PATHS = yaml.safe_load((REPO / "configs" / "paths.yaml").read_text(encoding="utf-8"))

# registered file -> pipeline step that writes it
PRODUCER = {
    "stage0_candidates.parquet": "01_candidates/candidates", "stage0_mp_meta.json": "01_candidates/candidates",
    "stage0_exfor_summary.json": "00_fetch_upstream/exfor", "exfor_benchmark.csv": "00_fetch_upstream/exfor",
    "cache_manifest.json": "00_fetch_upstream/mp_manifest", "stage0_model_control.json": "02_dynamics_mp/model_control",
    "stage1_triviality.parquet": "02_dynamics_mp/mp_models", "stage1_controls.parquet": "02_dynamics_mp/mp_models",
    "stage1_dimensionality.parquet": "02_dynamics_mp/dimensionality",
    "stage1_usable.parquet": "02_dynamics_mp/fidelity", "stage1_stats.json": "02_dynamics_mp/fidelity",
    "stage1_controls.json": "02_dynamics_mp/fidelity", "stage1_enantiomorphs.json": "02_dynamics_mp/enantiomorphs",
    "stage1r_fc_validation.parquet": "02_dynamics_mp/fc_validation", "stage1r_fc_validation.json": "02_dynamics_mp/fc_validation",
    "stage1r_floor_sensitivity.parquet": "02_dynamics_mp/floor_sensitivity",
    "stage1r_floor_sensitivity.json": "02_dynamics_mp/floor_sensitivity",
    "stage1r_pheasy_vs_dfpt.parquet": "02_dynamics_mp/reproducibility", "stage1r_reproducibility.json": "02_dynamics_mp/reproducibility",
    "stage1r_stability_screen.parquet": "02_dynamics_mp/stability_screen_sample",
    "stage1r_stability_screen.json": "02_dynamics_mp/stability_screen_sample",
    "stage1r_stability_full.parquet": "02_dynamics_mp/stability_screen_all",
    "stage1r_stability_full.json": "02_dynamics_mp/stability_screen_all",
    "stage1r_stability_full_mesh60.parquet": "02_dynamics_mp/stability_screen_mesh60",
    "stage1r_stability_full_mesh60.json": "02_dynamics_mp/stability_screen_mesh60",
    "stage1r_stability_summary.json": "02_dynamics_mp/stability_summary",
    "stage1r_fc_validation_largecell.parquet": "02_dynamics_mp/stability_summary",
    "stage3_stability.parquet": "02_dynamics_mp/d7a",
    "stage3_proxy.parquet": "01_candidates/proxy_bound", "stage3_proxy.json": "01_candidates/proxy_bound",
    "stage3_tiers.parquet": "01_candidates/tiers", "stage3_tiers.json": "01_candidates/tiers",
    "stage2_mlip_set.json": "03_dynamics_mlip/mlip_benchmark_set", "stage3_debye.parquet": "03_dynamics_mlip/debye_scores",
    "stage3_funnel.parquet": "03_dynamics_mlip/funnel", "stage3_mlip_set.json": "03_dynamics_mlip/funnel",
    "stage1_anchors.json": "05_validation/anchors_stage1", "stage2_mlip.json": "05_validation/mlip_eval",
    "stage2_mlip_benchmark.parquet": "05_validation/mlip_eval", "stage2_mlip_sensitivity.json": "05_validation/mlip_sensitivity",
    "stage2_orb_diagnostic.json": "05_validation/orb_diagnostic", "stage2_exfor.json": "05_validation/exfor_benchmark",
    "stage2_exfor_benchmark.parquet": "05_validation/exfor_benchmark",
    "stage2_exfor_binned.parquet": "05_validation/exfor_benchmark", "stage2_kill2_D11.json": "05_validation/kill2_d11",
    "stage2_exfor_sc_only.json": "05_validation/exfor_sc_only", "stage2_d11_robustness.json": "05_validation/d11_robustness",
    "stage2_mgo_audit.json": "05_validation/mgo_audit", "stage2_pbbi_qha.json": "05_validation/pbbi_qha",
    "stage2_pbbi_diag.json": "05_validation/pbbi_diag", "stage2_lit_compare.json": "05_validation/lit_compare",
    "stage2_residuals.json": "06_uncertainty/residuals", "stage2_uncertainty_model.json": "06_uncertainty/residuals",
    "stage2_text_numbers.json": "06_uncertainty/stage2_summary", "stage3_nucdata_check.json": "07_figures_of_merit/nucdata_check",
    "atlas.parquet": "08_atlas/assemble", "stage3_curves.parquet": "08_atlas/assemble", "atlas_counts.json": "03_dynamics_mlip/funnel + 08_atlas/assemble",
    "stage3_fronts.json": "08_atlas/fronts", "stage3_fronts.parquet": "08_atlas/fronts",
    "stage3_controls.json": "08_atlas/controls", "stage3_bragg.json": "08_atlas/bragg",
    "stage3_anchors.json": "08_atlas/anchors_stage3", "stage3_coherent.json": "08_atlas/coherent",
    "stage3_openmc_input.json": "08_atlas/openmc_input", "stage3_openmc_check.json": "08_atlas/openmc_check",
    "stage3_missed.json": "08_atlas/missed", "stage3_not_computed_table.csv": "08_atlas/missed",
    "stage3_laptop_gpu_check.json": "08_atlas/gpu_check", "stage3_text_numbers.json": "08_atlas/summary",
    "stage3_relax_sensitivity.json": "09_relax_sensitivity/relax_sensitivity",
    "v040_verification.json": "consistency check of the SevenNet outputs (not part of the pipeline)",
}
DIR_INFO = {   # data directory -> produced_by
    "upstream/mp": "upstream", "upstream/exfor": "upstream", "upstream/endfb80_hdf5": "upstream",
    "upstream/ncrystal_datacmp": "upstream", "intermediate/nucdata": "07_figures_of_merit/nucdata_extract",
    "intermediate/mlip": "03_dynamics_mlip/mlip_benchmark_runs", "intermediate/mlip_sensitivity": "05_validation/pbbi_qha",
    "intermediate/mlip_atlas": "03_dynamics_mlip/mlip_atlas_runs", "intermediate/atlas_eval": "04_cross_sections/evaluate_sevennet",
    "intermediate/logs": "03_dynamics_mlip (run logs)",
    "intermediate/mlip_relax_sensitivity": "03_dynamics_mlip/mlip_atlas_runs --relax-sensitivity",
    "intermediate/relax_sensitivity_eval": "09_relax_sensitivity/relax_sensitivity",
}


def upstream_source(rel):
    meta = json.loads((DATA / "upstream" / "mp" / "summary_meta.json").read_text())
    if rel.startswith("upstream/mp"):
        return f"Materials Project database {meta.get('db_version')} (fetched {meta.get('fetched')}; phonon DOS / force constants from the public Delta tables)"
    if rel.startswith("upstream/exfor"):
        return f"EXFOR entry {Path(rel).stem} (IAEA exfor_master GitHub mirror, retrieved 2026-10-20/21)"
    if rel.startswith("upstream/endfb80_hdf5"):
        return "ENDF/B-VIII.0, NNDC HDF5 build via openmc_data_downloader 0.6.2 (retrieved 2026-10-22)"
    if rel.startswith("upstream/ncrystal_datacmp"):
        return "NCrystal published data-comparison plot (NCrystal 4.4.6 documentation)"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", default="")
    ap.add_argument("--previous", default="", help="previous MANIFEST.csv: unchanged files keep its source")
    a = ap.parse_args()
    prev = {}
    if a.previous:
        with open(a.previous, encoding="utf-8") as fh:
            prev = {r["path"]: r for r in csv.DictReader(fh)}
    reg = {v: k for k, v in PATHS["files"].items()}
    rows = []
    for cls in ("upstream", "intermediate", "derived", "curated"):
        for p in sorted((DATA / cls).rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(DATA).as_posix()
            name = reg.get(rel)
            if cls == "curated":
                prod, src = "curated", "hand-built; see data/curated/*README*"
            elif cls == "upstream":
                prod, src = "upstream", upstream_source(rel)
            else:
                prod = PRODUCER.get(name) or DIR_INFO.get("/".join(rel.split("/")[:2]), "")
                if rel.startswith("intermediate/logs/"):
                    prod = DIR_INFO["intermediate/logs"]
                src = a.commit or "computed by this pipeline"
            h = sha256(p)
            old = prev.get(f"data/{rel}")
            if cls in ("intermediate", "derived") and old and old["sha256"] == h:
                src = old["source"]
            rows.append({"path": f"data/{rel}", "bytes": p.stat().st_size, "sha256": h, "class": cls,
                         "produced_by": prod, "source": src, "notes": old["notes"] if old else ""})
    out = DATA / "MANIFEST.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    by = {}
    for r in rows:
        by.setdefault(r["class"], [0, 0])
        by[r["class"]][0] += 1
        by[r["class"]][1] += r["bytes"]
    print(out, {k: (n, f"{b / 1e6:.1f} MB") for k, (n, b) in by.items()})


if __name__ == "__main__":
    main()
