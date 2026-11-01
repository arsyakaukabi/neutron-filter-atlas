# Troubleshooting

What this page is for: error messages you can meet, what causes them and how to fix them.
Who should read it: anyone whose command failed; search this page for the last line of the error.

The messages below are real output of the released version; paths are shortened to `path/to/...`.

## Missing data

```
FileNotFoundError: [Errno 2] No such file or directory: 'path/to/data/derived/atlas.parquet'
```

A step or snippet needs a data file that is not there. Unzip the archive that holds it into the repository root
([Data](data.md#where-the-files-go)); `data/MANIFEST.csv` lists every file with its class, and
`python scripts/check_data.py --class derived` (or `intermediate`, `upstream`) shows what is missing. The table in
[Data](data.md#what-runs-without-the-data) says which archive each kind of work needs.

The same error for `curated/hygroscopic_families.csv` (or another file under `curated/`) appears when `FA_DATA`
points to a data folder without the curated tables:

```
FileNotFoundError: [Errno 2] No such file or directory: 'path/to/data/curated/hygroscopic_families.csv'
```

The curated tables are tracked in git under `data/curated/`; copy that folder into the folder `FA_DATA` names.

```
FileNotFoundError: [Errno 2] No such file or directory: 'path/to/data/intermediate/nucdata/endfb80_elements.npz'
```

`filteratlas.fom.rejection` needs the elemental nuclear-data tables of the intermediate archive. Unzip it (or run
step 07 with OpenMC to rebuild the tables).

## A data file is reported as mismatched

```
manifest files: 65; missing: 0; mismatched: 1; extra (not in manifest): 0
  MISMATCH data/derived/stage0_mp_meta.json
Mismatched upstream files usually mean a re-fetch from the current Materials Project release (MP serves only its current release); derived results built from them can differ from the published atlas.
```

The file differs from the release. Causes, most likely first: you reran a step that rewrote it (this example: step 01
rewrites `stage0_mp_meta.json`, which records the run time in `wall_s`, so it never matches after a rerun; the other
steps reproduce their files exactly, see [Reproducing](reproducing.md)); the archive was not unpacked completely; a
curated table was edited; for `upstream/mp/`, the files were downloaded again from the current Materials Project,
which serves only its newest database version (the note at the end of the output is about this case). Unzip the
archive again to restore the released file. `check_data.py` always checks `data/` of the repository it belongs to,
not `FA_DATA`.

## A test fails: test_check_data_complete_tree

```
FAILED tests/test_check_data.py::test_check_data_complete_tree - AssertionErr...
1 failed, 30 passed
```

This test compares the whole data tree with the manifest, including the ENDF/B-VIII.0 library, which is not part of
the release archives. With the four archives unpacked and no ENDF library, this one failure is expected and every
other test passes. It goes away once the library is in `data/upstream/endfb80_hdf5/`, except that a new download has
its own `download.log` ([Data](data.md#the-endfb-viii0-library)).

## Missing optional dependency

```
    import openmc.data as od
ModuleNotFoundError: No module named 'openmc'
```

You asked for a sub-step that needs OpenMC (step 07 `nucdata_extract`, step 08 `openmc_check`, step 10 `transport`
and `activation`), or named it in `--only`, which runs it whatever its environment. Run it in an environment with
OpenMC (`environment.yml`, Linux or macOS), or leave it out: without `--with openmc` it is skipped with

```
-- skip nucdata_extract (needs environment 'openmc'; use --with openmc)
```

The same pattern applies to `sevenn`, `mace` and `orb_models` (learned-potential environments, `--with mlip`) and to
`mp_api` (main environment; `pip install mp-api` if missing).

## Materials Project API key absent

```
mp_api.client.core.exceptions.MPRestError: REST query returned with error status code 401 on URL https://api.materialsproject.org/materials/summary/?... with message:
Response {
  "message":"No API key found in request"
}
```

Step 00 needs a Materials Project API key. Copy `.env.example` to `.env` in the repository root and fill in
`MP_API_KEY` (or set `FA_ENV` to another file). `.env` is never committed. Only step 00 needs the key; every other
step reads the downloaded files. If you have the released upstream archive you do not need step 00 at all.

## A data file name is not registered

```
KeyError: "data file 'my_table.parquet' is not registered in configs/paths.yaml"
```

`F("name")` only knows the files listed under `files` in `configs/paths.yaml`. Add your file there (name and path
under the data root) or build the path yourself.

## Out of memory

**Worker processes.** The parallel steps (02 `mp_models`, 04, 08 `assemble`) use `n_workers` processes (4 by default,
`configs/compute.yaml`). If the machine runs out of memory (the process is killed, or `MemoryError` or a
`BrokenProcessPool` error appears), lower the number: `FA_WORKERS=2 python pipelines/04_cross_sections.py`.

**Learned-potential phonons (step 03).** Memory grows with the number of atoms in the supercell. The scheduler sends
supercells of up to `scheduler.gpu_max_atoms` atoms (330) to the GPU and larger ones to the CPU, keeps
`scheduler.ram_reserve_gb` free, and estimates CPU memory as `cpu_memory_model`. When a GPU job runs out of memory,
the error text (for example `CUDA out of memory` or `OutOfMemoryError`) is written to the run log, and after
`adaptive_failures` such failures within `adaptive_window` tries larger crystals are no longer sent to the GPU; set
`FA_NO_ADAPTIVE=1` to try each once anyway. Lower `gpu_max_atoms` for a GPU with less than 8 GB, and
`FA_GPU_MEM_FRACTION` caps the share of GPU memory a job may take. Supercells above 480 atoms are never attempted
(`max_supercell_atoms` in `configs/funnel.yaml`).

## Step 08 takes much longer than expected

`assemble` re-evaluates every SevenNet DOS (about 25 minutes with 4 workers) when the per-crystal cache in
`data/intermediate/atlas_eval/` is missing or not newer than the DOS files. Right after unzipping the release the
cache files carry the same timestamps as the DOS files, so the first run re-evaluates them all and rewrites the
cache; later runs take under a minute. The re-evaluated values are the same as the cached ones ([Reproducing](reproducing.md)).

## A command works in one folder but not another

Run commands from the repository root. The pipeline scripts find `src/` themselves, but the snippets on these pages
use paths relative to the root (`data/derived/...`). If you installed the package with `pip install -e .`, you can
also use `filteratlas.common.paths.F("atlas.parquet")`, which works from anywhere.

## Garbled characters on the Windows console

Some outputs contain characters such as Å or ₀. Set `PYTHONIOENCODING=utf-8` (the Makefile does) or use a UTF-8
terminal.
