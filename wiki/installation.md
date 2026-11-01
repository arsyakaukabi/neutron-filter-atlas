# Installation

What this page is for: installing the package from a clone, the optional dependencies, and checking that it works.
Who should read it: everyone who wants to run code from the repository, including the examples of the other pages.

## Requirements

- **Python 3.12** (`pyproject.toml` requires 3.12 or newer; the released atlas was computed with 3.12).
- **Operating system:** the main environment runs on Linux, macOS and Windows. OpenMC, which only a few steps need
  (see below), is available from conda-forge for Linux and macOS only; on Windows run those steps in a Linux
  environment (for example a virtual machine or a container).
- **Disk:** about 1 GB for the environment; the data need their own space (see [Data](data.md#disk-space)).
- **Memory:** 16 GB is enough for every step of the main environment with the default four worker processes.

## Install from a clone

```
git clone https://github.com/arsyakaukabi/neutron-filter-atlas.git
cd neutron-filter-atlas
```

Then one of the two routes below.

**Route 1, conda or micromamba (Linux, macOS).** `environment.yml` is the portable recipe and includes OpenMC:

```
micromamba create -f environment.yml
micromamba activate neutron-filter-atlas
pip install -e .
```

On Windows, `micromamba create -f environment.yml` fails because conda-forge has no OpenMC for Windows; delete the
`openmc` line from a copy of the file first, or use route 2.

**Route 2, the exact lock (Windows, or anywhere with uv).** `requirements-lock.txt` pins every package of the
environment used for the released atlas (Windows 11, Python 3.12). It does not contain OpenMC.

```
uv venv --python 3.12 .venv
uv pip install --python .venv -r requirements-lock.txt --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match
uv pip install --python .venv -e .
```

Activate it with `.venv\Scripts\activate` (Windows) or `source .venv/bin/activate`. The extra index provides the CPU
build of PyTorch, which one of the pinned packages requires.

`pip install -e .` makes `import filteratlas` work from anywhere. The pipeline scripts in `pipelines/` do not need it:
they add `src/` to the import path themselves.

## Optional dependencies

The table lists where each heavy dependency comes from. Most are in the main environment; OpenMC is in the main
environment only on Linux and macOS (route 1); the learned potentials need environments of their own. The pipeline
skips the sub-steps that need OpenMC, a learned potential or the network unless you ask for them with `--with`.

| Needed for | Package | Where it comes from | Steps |
|---|---|---|---|
| Thermal cross sections | NCrystal 4.4.6 | main environment | 01-05, 08-10 |
| Force-constant phonons and stability screens | phonopy 4.8.3, spglib | main environment | 02, 08 (and the learned-potential runs of 03 and 05) |
| Structures, compositions | pymatgen, ASE | main environment | all |
| Materials Project download | mp-api, emmet-core (and an API key) | main environment | 00 (`--with network`) |
| ENDF/B-VIII.0 download | openmc_data_downloader 0.6.2 | `environment.yml` (pip section) | manual, see [Data](data.md#the-endfb-viii0-library) |
| Nuclear-data extraction, slab transport, activation, transport check | OpenMC 0.16.0 | `environment.yml`, Linux and macOS | 07, 08, 10 (`--with openmc`) |
| Learned-potential phonons (SevenNet) | SevenNet and PyTorch | `requirements-mlip.txt` (CPU) or `requirements-mlip-gpu.txt` (CUDA) | 03, 05 (`--with mlip`) |
| Learned-potential benchmark (MACE-MP-0, ORB v3) | MACE, ORB | `requirements-mlip-benchmark.txt` | 03, 05 |

A learned-potential environment is a separate virtual environment, for example:

```
uv venv --python 3.12 .venv-mlip-cpu
uv pip install --python .venv-mlip-cpu -r requirements-mlip.txt
```

The scheduler of step 03 finds it under the paths set in `configs/compute.yaml` (`scheduler.python_cpu`,
`scheduler.python_gpu`) or in the environment variables `FA_PY_CPU` and `FA_PY_GPU`.

## Check the installation

Run the tests. Without data, the tests that need data are skipped with a message that names the missing file:

```
python -m pytest -q tests
```

On a fresh clone without data (output of the released version):

```
s....................s....sssss                                          [100%]
24 passed, 7 skipped in 0.84s
```

With data, pytest also prints a count of warnings from the libraries (for example `59 warnings`); they are harmless.

With the data release unpacked (see [Data](data.md)), 30 or 31 tests pass. One test,
`test_check_data_complete_tree`, compares the whole data tree with `data/MANIFEST.csv` and therefore also needs the
ENDF/B-VIII.0 library, which is not part of the release archives; without it the result is `1 failed, 30 passed`
(see [Troubleshooting](troubleshooting.md#a-test-fails-test_check_data_complete_tree)).

`make test` runs the same command. `make` is optional: every Makefile target is one Python command, listed in
[Pipeline overview](pipeline-overview.md#the-commands).

## Environment variables

| Variable | Effect |
|---|---|
| `FA_DATA` | data root to use instead of `data/` in the repository |
| `FA_REPO` | repository whose `configs/` are read (default: the one the package is installed from) |
| `FA_ENV` | `.env` file with API keys (default `.env` in the repository root) |
| `FA_WORKERS` | number of worker processes (default `n_workers` in `configs/compute.yaml`, 4) |
| `MP_API_KEY` | Materials Project API key, usually set in `.env` (copy `.env.example`); needed only by step 00 |
| `FA_PY_CPU`, `FA_PY_GPU` | Python of the learned-potential environments (step 03) |
| `FA_GPU_MEM_FRACTION`, `FA_NO_ADAPTIVE` | GPU memory cap and the out-of-memory guard of the scheduler (step 03) |
