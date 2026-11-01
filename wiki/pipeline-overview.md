# Pipeline overview

What this page is for: the steps from a crystal structure to an atlas entry, with the command, inputs, outputs, run
time and requirements of each.
Who should read it: users who want to recompute or extend the atlas.

## The flow

```
 Materials Project ──► 00 fetch ──► 01 candidates ──► 02 MP phonons ──► 03 learned-potential phonons
 EXFOR, ENDF/B-VIII.0     (network)    proxy bound,        VDOS models,        selection (funnel),
                                       tiers               stability screens   SevenNet runs
                                                                │                    │
                                                                ▼                    ▼
 05 validation ◄──────────────────────────────────── 04 cross sections of the SevenNet crystals
 (EXFOR benchmark,                                                │
  learned-potential benchmark)                                    │
        │                                                         │
        ▼                                                         ▼
 06 model uncertainty ──► (configs/uncertainty.yaml) ──► 08 atlas: assembly, fronts, controls, checks
                                                                  ▲            │
 07 nuclear data (ENDF/B-VIII.0 → elemental tables) ──────────────┘            ├──► 09 relaxation-limit sensitivity
                                                                               └──► 10 shortlist: rank probabilities,
                                                                                    transport, Bragg dips, temperature,
                                                                                    practicality, activation
```

Each step is one script `pipelines/NN_name.py` made of named sub-steps, and one Makefile target. Data flow between
steps goes only through files under `data/` (registered in `configs/paths.yaml`), so any step can be rerun on its own
when its inputs are in place.

## The commands

| Step | Script | Makefile target | Environment | Network |
|---|---|---|---|---|
| 00 | `pipelines/00_fetch_upstream.py --with network` | `make fetch` | main | yes (Materials Project API key, EXFOR) |
| 01 | `pipelines/01_candidates.py` | `make candidates` | main | no, if `data/upstream/mp` is complete |
| 02 | `pipelines/02_dynamics_mp.py` | `make dynamics-mp` | main | no |
| 03 | `pipelines/03_dynamics_mlip.py` | `make dynamics-mlip` | main; runs need a learned-potential environment (`--with mlip`) | no |
| 04 | `pipelines/04_cross_sections.py` | `make cross-sections` | main | no |
| 05 | `pipelines/05_validation.py` | `make validation` | main; two sub-steps need `--with mlip` | no |
| 06 | `pipelines/06_uncertainty.py` | `make uncertainty` | main | no |
| 07 | `pipelines/07_figures_of_merit.py` | `make figures-of-merit` | extraction needs OpenMC (`--with openmc`) | no |
| 08 | `pipelines/08_atlas.py` | `make atlas` | main; transport check needs `--with openmc` | no |
| 09 | `pipelines/09_relax_sensitivity.py` | `make relax-sensitivity` | main | no |
| 10 | `pipelines/10_shortlist.py` | `make shortlist` | main; transport and activation need `--with openmc` | no |

Run every command from the repository root as `python pipelines/NN_name.py [options]`. The Makefile adds
`PYTHONIOENCODING=utf-8` and `-W ignore` (warnings off); `make all` runs 01 to 10 in order. The options, the same for
every script:

| Option | Effect |
|---|---|
| `--list` | list the sub-steps with their environment and a one-line description |
| `--only a,b` | run only these sub-steps (whatever environment they need) |
| `--from a` | start at sub-step `a` and run the rest |
| `--with openmc,mlip,network` | also run the sub-steps that need these environments; `main` sub-steps always run |

Sub-steps marked `archive` need inputs that are not part of the data release and are skipped. The output of
`--list` for step 08:

```
$ python pipelines/08_atlas.py --list
assemble                     [main] atlas.parquet, stage3_curves.parquet, atlas_counts.json (atlas part)
fronts                       [main] Pareto fronts, fidelity ladder, who moves, sensitivity, front test preview
controls                     [main] positive, negative, model and drift controls
bragg                        [main] Bragg cut-off section (Tier A)
anchors_stage3               [main] sapphire literature anchors
coherent                     [main] coherent one-phonon check, sapphire and Si
openmc_input                 [main] input of the OpenMC check
openmc_check                 [openmc] uncollided vs transported A_fast (OpenMC 0.16)
missed                       [main] what the funnel could have missed; supplementary table
gpu_check                    [archive] GPU vs CPU reproducibility (archive-class inputs)
summary                      [main] summary numbers of the atlas
```

Some sub-step names are short labels: in step 05, `kill2_d11` is the acceptance criterion counted on single-crystal
data only and `d11_robustness` its robustness variants; in step 02, `d7a` is the stability screen of the crystals
that pass the proxy bound; `anchors_stage1` and `anchors_stage3` compare with literature values; `stage2_summary`
(step 06) writes the summary numbers of the validation. In the descriptions, "fidelity" means the comparison of the levels of phonon description, "who
moves" the crystals whose rank changes most between levels or temperatures, "drift controls" crystals computed with
both DFT and SevenNet phonons, and SIMEX a noise-correction method ([Glossary](glossary.md)).

## The steps

**00 Fetch upstream data.** Downloads the Materials Project summary table, structures, phonon DOS and force
constants, and the EXFOR entries, into `data/upstream/`; writes the candidate table, the EXFOR table
(`data/intermediate/exfor_benchmark.csv`) and the checksums of the inputs (`cache_manifest.json`). Cached: files
already present are not downloaded again. Needs `MP_API_KEY` in `.env`. The Materials Project serves only its current
database version, so a new download can differ from the release; reproduce from the released upstream archive
instead. The ENDF/B-VIII.0 library is fetched separately ([Data](data.md#the-endfb-viii0-library)).

**01 Candidates, proxy bound, tiers.** Applies the element filters and the hull cut (`candidates`), computes the
composition-only proxy and excludes crystals that cannot transmit half the band through 1 cm (`proxy_bound`), and
matches the vendor table to Materials Project entries (`tiers`). About 10 s.

**02 Materials Project phonons.** Converts every Materials Project phonon DOS to per-element vibrational DOS (2 meV
floor), screens dynamical stability with the DOS rule and evaluates the thermal quantities with NCrystal
(`mp_models`, about 1-2 h with 4 workers); then the statistics of proxy against full model (`fidelity`), the
validation of the 2 meV floor against force-constant DOS, and the force-constant stability screens (phonopy on about
2,500 crystals, hours).

**03 Learned-potential phonons.** Chooses the crystals that get learned-potential phonons: the benchmark set, the
elastic-Debye scores, the budget of 1,500 by priority (`funnel`). The runs themselves are separate processes in a
learned-potential environment: `python -m filteratlas.dynamics.mlip_run <model>` for the benchmark and
`python -m filteratlas.dynamics.scheduler` for the atlas (one process per crystal; machine settings in
`configs/compute.yaml`). Each SevenNet calculation takes seconds to about twenty minutes depending on the supercell;
the logged runs of the released atlas add up to about 80 hours of wall time (`wall_s` in
`data/intermediate/logs/mlip_atlas.jsonl`). The released phonon DOS are in the intermediate archive.

**04 Cross sections of the SevenNet crystals.** Passes every SevenNet DOS through the same conversion and NCrystal
code as the DFT DOS and caches one result per crystal in `data/intermediate/atlas_eval/`. About 25 minutes for the
released DOS with 4 workers.

**05 Validation.** Compares the models with measured total cross-sections from EXFOR and evaluates the acceptance
criterion (`exfor_benchmark`), with the count restricted to single-crystal data and its variants; compares the
learned potentials with DFT (`mlip_eval`). The sub-steps that need a learned-potential environment are skipped by
default.

**06 Model uncertainty.** Turns the validation residuals into the model uncertainty. The values the atlas uses are
copied into `configs/uncertainty.yaml`; `tests/test_configs.py` checks that the two agree. Seconds.

**07 Nuclear data.** Extracts elemental total cross sections (1 eV-10 MeV, 5,000 points per decade) and photon
attenuation at 1 and 6 MeV from ENDF/B-VIII.0 into `data/intermediate/nucdata/` (OpenMC needed), then checks the grid
against the exact union grid. The extracted tables are in the intermediate archive.

**08 Atlas.** Assembles the atlas table from the earlier outputs (`assemble`), builds the Pareto fronts and the
fidelity statistics (`fronts`), runs the controls and the Bragg cut-off, coherent and literature checks, and writes
the coverage table of the crystals without a full model (`missed`) and the summary numbers. `assemble` re-evaluates
the SevenNet DOS when the cache of step 04 is missing or older than the DOS files.

**09 Relaxation-limit sensitivity.** Reads the atlas and the learned-potential reruns with 2,000 relaxation steps
(kept apart from the atlas, which uses 500) and reports how the result would change. No atlas file depends on it.

**10 Shortlist.** Selects the shortlist, computes Monte Carlo rank probabilities (20,000 draws), prepares and runs
OpenMC slab transport (all transmitted particles), scans Bragg dips over 500 beam directions (up to 2 h per crystal),
computes the temperature series 293-77 K and the activation, and collects the practicality data.

## Expected run times

With 4 worker processes (`FA_WORKERS`) from the released inputs. "Measured" values come from a rerun on a
workstation for this documentation; the others are the estimates in the module notes.

| Step | Time | Source |
|---|---|---|
| 01 | about 10 s | measured |
| 02 | several hours (NCrystal models 1-2 h, force-constant screens most of the rest) | module notes |
| 03 | about 80 h of summed wall time for the atlas runs | run log |
| 04 | about 25 min | measured (as part of step 08) |
| 05 | about 4 min (main sub-steps) | measured |
| 06 | about 2 s | measured |
| 08 `assemble` | about 24 min when the SevenNet evaluation is redone, 45 s with the cache | measured |
| 08 `fronts` | about 15 s | measured |
| 10 | hours (Bragg-dip scans up to 2 h per crystal, and transport) | module notes |
