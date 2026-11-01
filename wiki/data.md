# Data

What this page is for: the data release, what each archive holds, where to put it, how to verify it and what runs
without it.
Who should read it: everyone; the atlas itself is in the data release, not in the git repository.

## The data release

**Neutron filter atlas data release, version 0.5.0**, produced by the pipeline of the same version. It is available
from the authors on request (the e-mail addresses are in `CITATION.cff`) and will be archived openly under a
persistent identifier; the citation is in [Citing](citing.md). Only `data/README.md`, `data/MANIFEST.csv` and `data/curated/` are tracked in git.

The release is a folder with four zip archives, a checksum list and the manifest:

| Archive | Holds | Packed | Files | Needed for |
|---|---|---|---|---|
| `derived/nfa-data-derived-v0.5.0.zip`, the **results archive** | the atlas and every computed result | about 6 MB | 65 | reading the atlas (start here) |
| `intermediate/nfa-data-intermediate-v0.5.0.zip` | learned-potential phonon DOS, per-crystal evaluation cache, run logs, transport and activation tallies, large per-crystal tables, the elemental nuclear-data tables | about 420 MB | 4,942 | recomputing steps 04-10, computing your own crystals |
| `upstream/nfa-data-upstream-mp-v0.5.0.zip` | Materials Project structures, phonon DOS and force constants as retrieved for this release | about 830 MB | 9,025 | recomputing from the phonons (steps 01-10) |
| `upstream/nfa-data-upstream-exfor-ncrystal-v0.5.0.zip` | the EXFOR entries of the validation and one NCrystal comparison plot | about 4 MB | 74 | recomputing the validation (step 05) |
| `MANIFEST.csv` | every data file of the version: path, size, SHA-256, class, producing step, source | 3.7 MB | | verifying |
| `SHA256SUMS.txt` | the SHA-256 of each archive and of the manifest | | | verifying |

Not in the archives: the ENDF/B-VIII.0 library (public, fetched separately, see
[below](#the-endfb-viii0-library)) and the hand-made tables, which are in the repository under `data/curated/`
([Curated tables](curated-tables.md)). The manifest lists 14,416 files: 14,106 in the archives, 304 of the ENDF
library and 6 curated files.

**Licences.** See [Licences](#licences) at the end of this page.

## Where the files go

Unzip every archive into the **repository root**; each archive unpacks into a folder named `data`:

```
unzip path/to/release/derived/nfa-data-derived-v0.5.0.zip -d .
unzip path/to/release/intermediate/nfa-data-intermediate-v0.5.0.zip -d .
unzip path/to/release/upstream/nfa-data-upstream-mp-v0.5.0.zip -d .
unzip path/to/release/upstream/nfa-data-upstream-exfor-ncrystal-v0.5.0.zip -d .
```

The result is the tree below. Every path the code uses is registered in `configs/paths.yaml` and resolved by
`filteratlas.common.paths` (`D("name")` for a directory, `F("name")` for a file). To keep the data elsewhere, set the
environment variable `FA_DATA` to the data folder (the folder that contains `derived/`).

```
data/
  README.md, MANIFEST.csv      tracked in git
  curated/                     hand-made tables (tracked in git)
  upstream/
    mp/                        Materials Project: summary table, structures, phonon_dos/{pheasy,dfpt}/, force_constants/
    exfor/                     EXFOR entries (.x4)
    endfb80_hdf5/              ENDF/B-VIII.0 HDF5 files (fetched separately)
    ncrystal_datacmp/          NCrystal's published sapphire comparison plot
  intermediate/                computed products that later steps read (phonon DOS, caches, logs, large tables)
  derived/                     the atlas and the results
```

## Verify the data

First the archives, inside the release folder (Linux, macOS):

```
sha256sum -c SHA256SUMS.txt
```

On Windows, `certutil -hashfile <archive> SHA256` prints the checksum to compare by eye.

Then the unpacked files against the manifest, from the repository root. Check only what you unpacked:

```
python scripts/check_data.py --class derived
```

A pass looks like this (exit code 0):

```
manifest files: 65; missing: 0; mismatched: 0; extra (not in manifest): 0
```

The full check, `make check-data` or `python scripts/check_data.py`, compares every file of the manifest (SHA-256;
a few minutes). With all four archives but without the ENDF library it fails, and says why (exit code 1):

```
manifest files: 14416; missing: 304; mismatched: 0; extra (not in manifest): 0
  MISSING data/upstream/endfb80_hdf5/cross_sections.xml
  MISSING data/upstream/endfb80_hdf5/download.log
  MISSING data/upstream/endfb80_hdf5/ENDFB-8.0-NNDC_Al.h5
  ...
```

Those 304 missing files are the ENDF library; everything else is present and identical. A **mismatched** file
differs from the release: an unpacking error, a file you changed, or (for `upstream/mp`) a re-download from the
current Materials Project, which serves only its newest database version. `--fast` compares sizes only. Classes for
`--class`: `upstream`, `intermediate`, `derived`, `curated` (comma-separated).

## The ENDF/B-VIII.0 library

The fast-neutron and photon cross sections come from ENDF/B-VIII.0 (NNDC HDF5 build, about 2.3 GB, 304 files). You
need the library only for the steps that run OpenMC: the extraction of the elemental tables (step 07), the transport
check (step 08), and the slab transport and activation of the shortlist (step 10). The elemental tables that every
other step uses are in the intermediate archive (`data/intermediate/nucdata/`).

The release used `openmc_data_downloader` 0.6.2 (installed by `environment.yml`) with the library `ENDFB-8.0-NNDC`,
neutron and photon data, into `data/upstream/endfb80_hdf5/`. The command has this form (it needs the OpenMC
environment and was not rerun for this page; `openmc_data_downloader --help` lists every option):

```
openmc_data_downloader -l ENDFB-8.0-NNDC -p neutron photon -e <elements> -d data/upstream/endfb80_hdf5
```

with the elements of the candidate set after `-e`. A new download reproduces the HDF5 files but not the
`download.log` of the release, so the full manifest check keeps reporting that one file.

## Disk space

| Content | Unpacked |
|---|---|
| `derived/` (results archive) | about 10 MB |
| `intermediate/` | about 510 MB |
| `upstream/mp/` and `upstream/exfor/`, `upstream/ncrystal_datacmp/` | about 830 MB |
| `upstream/endfb80_hdf5/` (fetched separately) | about 2.3 GB |
| Everything | about 3.7 GB |

## What runs without the data

| You have | You can |
|---|---|
| nothing (a fresh clone) | run the tests (data tests are skipped); compute NCrystal's own materials with the [Python interface](api-cross-sections.md) (`evaluate("stdlib::Si_sg227.ncmat", ...)`); read the configuration |
| the results archive | everything on the *Reading the atlas* pages: the atlas, the fronts, the validation, the shortlist |
| + the intermediate archive | compute L₀.₈ and the rejections of your own crystal from a phonon DOS ([Extending](extending.md)); rerun steps 04 and 06 and the atlas assembly |
| + both upstream archives | rerun every step except 00, the learned-potential runs and the OpenMC sub-steps, without network access |
| + the ENDF library and OpenMC | rerun the OpenMC steps (07, 08 transport check, 10 transport and activation) |
| a Materials Project API key | download upstream data again (step 00; the result can differ from the release) |
| a learned-potential environment | compute new learned-potential phonons (step 03) |

## What the results archive holds

Tables (columns in [Atlas columns](atlas-columns.md) and [Results tables](results-tables.md)):

| File | Contents |
|---|---|
| `atlas.parquet` | the atlas, one row per crystal and temperature |
| `stage0_candidates.parquet` | the candidate set with its Materials Project fields |
| `stage3_funnel.parquet`, `stage3_tiers.parquet` | which level of physics each crystal received; the commercially grown crystals |
| `stage3_fronts.parquet` | membership of every Pareto front |
| `stage3_not_computed_table.csv` | crystals without a full model, with their elastic-Debye and bound positions |
| `stage2_exfor_benchmark.parquet`, `stage2_exfor_binned.parquet` | validation against measured total cross-sections |
| `stage2_mlip_benchmark.parquet` | learned-potential benchmark |
| `stage1_controls.parquet` | control crystals |
| `stage4_shortlist.parquet` | the shortlist |

JSON files (structured results; each holds its own field names):

| Files | Contents |
|---|---|
| `atlas_counts.json`, `stage3_proxy.json`, `stage3_tiers.json`, `stage3_mlip_set.json`, `stage2_mlip_set.json` | counts of the funnel and the atlas; the proxy bound; the tier mapping; the learned-potential selections |
| `stage1_stats.json`, `stage1_controls.json`, `stage1_enantiomorphs.json`, `stage1_anchors.json` | triviality test (proxy against full model), controls, enantiomorph pairs, literature anchors |
| `stage1r_*.json` | dynamical-stability screens, the validation of the 2 meV floor, reproducibility checks |
| `stage2_exfor.json`, `stage2_exfor_sc_only.json`, `stage2_kill2_D11.json`, `stage2_d11_robustness.json`, `stage2_mgo_audit.json`, `stage2_lit_compare.json` | the acceptance criterion against measured cross-sections, the count restricted to single-crystal data and its robustness variants, the MgO audit, the sapphire literature comparison |
| `stage2_mlip.json`, `stage2_mlip_sensitivity.json`, `stage2_orb_diagnostic.json`, `stage2_pbbi_diag.json`, `stage2_pbbi_qha.json` | learned-potential benchmark and diagnostics |
| `stage2_residuals.json`, `stage2_uncertainty_model.json` | residual classes and the model uncertainty |
| `stage3_fronts.json`, `stage3_controls.json`, `stage3_bragg.json`, `stage3_coherent.json`, `stage3_anchors.json`, `stage3_missed.json`, `stage3_relax_sensitivity.json` | fronts and their sensitivity runs, atlas controls, the Bragg cut-off section, the coherent one-phonon check, sapphire anchors, what the funnel could have missed, the relaxation-limit sensitivity |
| `stage3_nucdata_check.json`, `stage3_openmc_input.json`, `stage3_openmc_check.json` | nuclear-data grid check; the transport check of the uncollided rejection |
| `stage3_laptop_gpu_check.json`, `v040_verification.json` | reproducibility checks of the learned-potential runs (GPU against CPU; consistency of the stored outputs) |
| `stage4_rank_probabilities.json`, `stage4_transport.json`, `stage4_dips.json`, `stage4_temperature.json`, `stage4_dossier.json`, `stage4_controls.json` | shortlist: Monte Carlo rank probabilities, slab transport, Bragg dips and orientation, temperature series, practicality, controls |
| `stage2_text_numbers.json`, `stage3_text_numbers.json`, `stage4_text_numbers.json` | summary numbers of the validation, the atlas and the shortlist |
| `cache_manifest.json` | SHA-256 of the upstream inputs as retrieved |

The file names are part of the release and its checksums and are kept as they are. Their prefixes group the files by
the part of the computation that produced them:

| Prefix | Produced by | Contents |
|---|---|---|
| `stage0_` | steps 00-02 | candidate set, Materials Project metadata, model control |
| `stage1_`, `stage1r_` | steps 02 and 05 | proxy against full model, controls, stability screens, floor and reproducibility checks |
| `stage2_` | steps 03, 05 and 06 | validation, learned-potential benchmark, model uncertainty |
| `stage3_` | steps 01, 03, 07, 08 and 09 | proxy bound, tiers, funnel, atlas checks and fronts |
| `stage4_` | step 10 | shortlist |

The column `produced_by` of `data/MANIFEST.csv` names the exact step and sub-step of every file.

## Licences

- Computed data (`intermediate/`, `derived/`): under the licence of the code (MIT).
- Materials Project data (`upstream/mp/`): CC BY 4.0.
- EXFOR (`upstream/exfor/`): free use; cite the original measurements.
- ENDF/B-VIII.0: public.
- NCrystal comparison plot: from the NCrystal documentation.

Cite the release together with the upstream sources ([Citing](citing.md)).
