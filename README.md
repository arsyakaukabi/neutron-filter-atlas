# neutron-filter-atlas

Pipeline for a database-wide atlas of single-crystal thermal-neutron filters. For every candidate crystal from the Materials Project, it computes the off-Bragg thermal transmission with NCrystal, using Materials Project DFT phonons or phonons from a machine-learned interatomic potential (SevenNet), and from it the reference thickness L₀.₈ that transmits 80 % of a thermal beam at 293 K and 77 K. At that thickness it computes the rejection of fast neutrons and of photons from ENDF/B-VIII.0 data. A benchmark against measured EXFOR total cross-sections gives the model uncertainty, which is carried into the fronts of the atlas and into Monte Carlo rank probabilities. For a shortlist of grown crystals it adds energy-resolved slab transport with OpenMC, Bragg dips, a temperature series and practicality data.

## Install

```
micromamba create -f environment.yml && micromamba activate neutron-filter-atlas
pip install -e .            # optional; the pipelines add src/ to the path themselves
```

OpenMC (used by steps 07, 08 and 10) is available from conda-forge on Linux and macOS; on Windows run those steps in a Linux environment. The phonon computations with machine-learned potentials need separate environments: `requirements-mlip.txt` (CPU), `requirements-mlip-gpu.txt` (CUDA) and `requirements-mlip-benchmark.txt` (MACE-MP-0 and ORB v3). `requirements-lock.txt` pins the main environment exactly (Windows; install command in its header).

## Data

The data products of version 0.5.0 (about 3.7 GB; four zips with SHA-256 checksums) are available from the authors on request and will be archived under a DOI on acceptance of the accompanying article. Unzip them into the repository root (the paths inside start with `data/`) and verify them:

```
make check-data             # every file against data/MANIFEST.csv (size and SHA-256)
python scripts/check_data.py --class derived    # only the atlas and the computed results
```

The ENDF/B-VIII.0 library is public and is fetched with `openmc_data_downloader`. `data/README.md` describes the tree and which files each use needs. Upstream data can be re-fetched with `make fetch` (needs `MP_API_KEY` in `.env`; see `.env.example`), but the Materials Project serves only its current release, so a re-fetch may differ from the manifest.

## Run

Each step is one script in `pipelines/` and one Makefile target, in data-flow order:

```
make candidates dynamics-mp dynamics-mlip cross-sections validation uncertainty figures-of-merit atlas relax-sensitivity shortlist
python pipelines/08_atlas.py --list         # the sub-steps of a step
python pipelines/08_atlas.py --only fronts  # run selected sub-steps
```

Sub-steps that need OpenMC, a machine-learned-potential environment or the network are skipped unless enabled with `--with openmc`, `--with mlip` or `--with network`. Every tunable number is in `configs/`, and every data path in `configs/paths.yaml`; set `FA_DATA=/path/to/data` to move the data root.

## Tests

```
make test                   # pytest; data-dependent tests skip cleanly when data/ is empty
```

## Cite

See `CITATION.cff`. Please also cite the upstream sources: the Materials Project and its phonon database, the EXFOR entries used, ENDF/B-VIII.0, NCrystal and SevenNet.

## Licence

Code: MIT (`LICENSE`). Data: see `data/README.md`; upstream data keep their own licences.
