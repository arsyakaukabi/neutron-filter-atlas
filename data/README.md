# data/

**The data products of version 0.5.0 are distributed separately (four zips with SHA-256 checksums). They are available from the authors on request and will be archived under a DOI on acceptance of the accompanying article.** Unzip them into the repository root so that `data/upstream`, `data/intermediate` and `data/derived` match this tree, then run `make check-data`.

Only this file, `MANIFEST.csv` and `curated/` are tracked in git. `MANIFEST.csv` lists every data file with its path, size, SHA-256, class, the pipeline step that produced it and its source.

## Tree

```
data/
  README.md                this file (tracked)
  MANIFEST.csv             path, bytes, sha256, class, produced_by, source, notes (tracked)
  curated/                 hand-built tables (tracked): growable_crystals.csv (+README), hygroscopic_families.csv,
                           exfor_candidates_extra.csv, stage4_dossier.csv (+README)
  upstream/                downloaded from third parties
    mp/                    Materials Project: summary_exp_nearhull.parquet, summary_meta.json, structures.jsonl.gz,
                           phonon_dos/{pheasy,dfpt}/<id>.json.gz, force_constants/<id>.npz
    exfor/                 EXFOR entries <entry>.x4
    endfb80_hdf5/          ENDF/B-VIII.0 NNDC HDF5 (neutron + photon) and cross_sections.xml
    ncrystal_datacmp/      NCrystal's published sapphire data-comparison plot (qualitative check only)
  intermediate/            expensive computed products that later steps read
    mlip/                  phonon DOS of the learned-potential benchmark, per model
    mlip_atlas/            SevenNet phonon DOS for the atlas
    mlip_relax_sensitivity/, relax_sensitivity_eval/   sensitivity reruns with 2,000 relaxation steps and their
                           evaluation (kept separate from the atlas, which uses the 500-step rule)
    mlip_sensitivity/      Pb/Bi quasi-harmonic DOS (diagnostic)
    stage4_transport/, stage4_activation/              OpenMC inputs and tallies of the slab transport, activation
                           inputs and results
    atlas_eval/            per-compound cache of the SevenNet-DOS evaluation (regenerable by step 04)
    nucdata/               elemental ENDF/B-VIII.0 tables (endfb80_elements.npz, endfb80_meta.json)
    logs/                  the run logs the code reads (MLIP status and wall times, *.jsonl)
    *.parquet, *.csv       large per-compound tables (MP-DOS models, stability screens, proxy, Debye, curves, EXFOR table)
  derived/                 the atlas and the computed results (atlas.parquet, atlas_counts.json and the JSON, parquet
                           and CSV files of each pipeline step)
```

The data root can be moved: set `FA_DATA=/path/to/data` and every path in `configs/paths.yaml` resolves under it.

## Which files you need

| Use | Needed | Approx. size |
|---|---|---|
| **(a) Read the atlas** and the results | `derived/` | about 10 MB |
| **(b) Rerun atlas assembly, validation and uncertainty** without recomputing phonons (steps 04–08) | `derived/` + `intermediate/` + `upstream/endfb80_hdf5/` + `upstream/exfor/` + `upstream/mp/` | about 3.7 GB (upstream 3.2 GB, of which ENDF 2.4 GB; intermediate 0.52 GB) |
| **(c) Full reproduction** | everything, or re-fetch upstream with `make fetch` | about 3.7 GB + the learned-potential phonon computations |

Step 08 reads the MP phonon DOS (`upstream/mp/phonon_dos`) for the Bragg section and the controls, so (b) includes `upstream/mp/`.

`make fetch` needs a Materials Project API key in `.env` (see `.env.example`). **The Materials Project serves only its current release**, so a re-fetch may differ from the manifest; `make check-data` will say so. Use the distributed upstream files to reproduce the atlas exactly.

Check only what you need:

```
python scripts/check_data.py --class derived          # (a)
python scripts/check_data.py                           # everything (SHA-256; a few minutes)
python scripts/check_data.py --fast                    # sizes only
```

## Licences of the data

- **Materials Project data:** CC BY 4.0.
- **EXFOR:** free use; cite the original measurements.
- **ENDF/B-VIII.0:** public.
- **NCrystal plot:** NCrystal documentation.
- **Computed data (intermediate, derived):** produced by this repository's code (MIT). Cite them as in `CITATION.cff`, together with the upstream sources: the Materials Project, EXFOR, ENDF/B-VIII.0 and NCrystal.
