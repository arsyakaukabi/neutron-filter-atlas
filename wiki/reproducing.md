# Reproducing

What this page is for: how to recompute the released atlas from the released data, how to check the result against
the release, and which differences to expect.
Who should read it: pipeline users who want to confirm the release or start their own changes from a verified base.

## What can be reproduced, and from what

The pipeline of version 0.5.0 run on the data release of version 0.5.0 reproduces the release. How far back you start
decides what you need and how long it takes:

| Start from | You need | Steps to run | Time with 4 workers |
|---|---|---|---|
| the computed phonons (recommended) | all four archives | 01, 04 (optional), 05, 06, 08, 09, 10 | under an hour, plus the Bragg-dip scans and transport of step 10 |
| the Materials Project phonons | all four archives | 01, 02, then as above | plus several hours (step 02) |
| a new download | an API key, the learned-potential environments, OpenMC and the ENDF library | 00 to 10 | days; differs from the release (below) |

Use the released upstream archives, not a new download: the Materials Project serves only its current database
version, so a new download can change structures, phonons and even the set of entries.

## Recompute

Recompute in a clone of its own with all four archives unpacked: the steps overwrite their outputs in `data/`, and
the check below compares `data/` of that clone with the manifest. Unzipping the archives again restores the release.
From the repository root:

```
python pipelines/01_candidates.py
python pipelines/05_validation.py
python pipelines/06_uncertainty.py
python pipelines/08_atlas.py --only assemble,fronts
```

Measured on the released data (4 workers), steps 01, 06 and 08 (step 05 prints one such line per sub-step and took
3 min 50 s, most of it in `mlip_eval` and `exfor_benchmark`):

```
== candidates done in 1 s
== proxy_bound done in 6 s
== tiers done in 1 s
== residuals done in 1 s
== stage2_summary done in 0 s
evaluated 1507 SevenNet DOS in 1368s
== assemble done in 1409 s
== fronts done in 13 s
```

The 1,507 evaluated DOS are the SevenNet files of the atlas runs (the drift-control set is skipped); 1,279 of them
become full models, the others are rejected by the relaxation or stability rules (`funnel_step`, `mlip_unstable`).

`assemble` re-evaluated the SevenNet DOS because the per-crystal cache of the release is not newer than the DOS files
after unzipping ([Troubleshooting](troubleshooting.md#step-08-takes-much-longer-than-expected)); a second run uses
the cache and takes under a minute (43 s). Run further steps the same way ([Pipeline overview](pipeline-overview.md)); steps 02 and
10 take hours.

## Check the result against the release

The package's check is `scripts/check_data.py`, which compares every file with its size and SHA-256 in
`data/MANIFEST.csv` (it always checks `data/` of the repository it belongs to, not `FA_DATA`). A recomputation that
reproduces the release byte for byte passes:

```
python scripts/check_data.py --class derived
manifest files: 65; missing: 0; mismatched: 0; extra (not in manifest): 0
```

Reruns of steps 05, 06 and 08 (`assemble`, `fronts`) in a fresh clone with the released data, in the environment
that `requirements-lock.txt` pins, give exactly this result: every file they rewrite is identical to the release.
Step 01 rewrites one file that cannot match, because it records the run time of the step:

```
manifest files: 65; missing: 0; mismatched: 1; extra (not in manifest): 0
  MISMATCH data/derived/stage0_mp_meta.json
```

The only difference in that file is the field `wall_s` (seconds the download and selection took); everything else
in it, and every other file of step 01, is identical.

When bytes differ (another platform, other library versions), compare values instead. For the atlas:

```python
import numpy as np, pandas as pd
a = pd.read_parquet("path/to/release/data/derived/atlas.parquet").sort_values(["material_id", "T"]).reset_index(drop=True)
b = pd.read_parquet("data/derived/atlas.parquet").sort_values(["material_id", "T"]).reset_index(drop=True)
assert a[["material_id", "T"]].equals(b[["material_id", "T"]])
for c in a.columns:
    if a[c].dtype.kind == "f":
        x, y = a[c].to_numpy(), b[c].to_numpy()
        ok = ~np.isnan(x) & ~np.isnan(y)
        rel = np.max(np.abs(x[ok] - y[ok]) / np.maximum(np.abs(x[ok]), 1e-12)) if ok.any() else 0.0
        if rel > 0 or not np.array_equal(np.isnan(x), np.isnan(y)):
            print(c, f"largest relative difference {rel:.1e}")
    elif not a[c].equals(b[c]):
        print(c, "differs")
```

On the run above it prints nothing. As a guide: relative differences below about 10⁻⁶ are rounding and change
nothing; a difference larger than about 1 % of `sigma_model` in `L08` (that is, 0.001 in ln L₀.₈) means an input or a
setting differs.

Then run the tests, which include the controls of the atlas (known filters reproduce their validation values, strong
absorbers rank at the bottom, counts are consistent):

```
python -m pytest -q tests
```

## Deterministic and non-deterministic steps

| Part | Deterministic? | Expected difference |
|---|---|---|
| NCrystal evaluation (steps 01, 02, 04, 08), L₀.₈ and rejections | yes, for the same inputs and library versions | none on the same platform; rounding (≲10⁻⁶ relative) across platforms |
| Bootstrap intervals, Monte Carlo rank probabilities, random samples | yes: every random generator has a fixed seed in `configs/` (`seed`, `screen_seed`, `rank_mc.seed`, ...) | none with the same numpy version; a numpy release that changes its random streams changes intervals and probabilities within their statistical spread |
| Force-constant screens and DOS (phonopy, step 02) | yes | none |
| Learned-potential phonons (step 03) | no: GPU and CPU differ, and SevenNet runs in single precision | GPU against CPU on 18 benchmark crystals: median \|Δ ln L₀.₈\| 0.0004 at 293 K, largest 0.007 (`stage3_laptop_gpu_check.json`), small against σ_model = 0.138 |
| OpenMC transport and activation (step 10) | statistical, run to a target relative uncertainty of 2 % | within the stated uncertainties of `stage4_transport.json` |
| A new Materials Project download (step 00) | no | entries, structures and phonons can change with the database version |

The released atlas uses the released learned-potential DOS, so steps 04 and 08 reproduce it exactly; only a rerun
of step 03 brings in the differences of the learned potential.
