# Quickstart

What this page is for: five minutes with the atlas: load it, look up a crystal, rank, filter and export.
Who should read it: users of the atlas; you need Python with pandas and pyarrow, and the results archive.

Set-up: install as in [Installation](installation.md) (any Python with `pandas` and `pyarrow` is enough for this
page), unzip the results archive into the repository root ([Data](data.md#where-the-files-go)) and check it:

```
python scripts/check_data.py --class derived
manifest files: 65; missing: 0; mismatched: 0; extra (not in manifest): 0
```

Run the snippets below one after the other from the repository root, in one Python interpreter. The outputs are from the released data,
version 0.5.0.

## 1. Load the atlas

```python
import pandas as pd
pd.set_option("display.width", 120)
at = pd.read_parquet("data/derived/atlas.parquet")
print(at.shape)
print(at["T"].value_counts())
```

```
(36898, 63)
T
293.15    18449
77.00     18449
Name: count, dtype: int64
```

One row per crystal and temperature: 18,449 crystals at 293.15 K and at 77 K. The 63 columns are explained in
[Atlas columns](atlas-columns.md). With the package installed, `filteratlas.common.paths.F("atlas.parquet")` gives
the same path and follows `FA_DATA`.

## 2. Look up one crystal

Look crystals up by formula and space group (several entries can share a formula) or by `material_id`. Sapphire:

```python
cols = ["material_id", "formula", "sg_symbol", "T", "vdos_source", "L08", "sigma_model", "A_fast_1E", "A_gamma_1MeV"]
print(at[(at.formula == "Al2O3") & (at.sg_symbol == "R-3c")][cols].to_string(index=False))
```

```
material_id formula sg_symbol      T vdos_source       L08  sigma_model  A_fast_1E  A_gamma_1MeV
mp-aaaaabrz   Al2O3      R-3c 293.15      mp_dft  5.977090        0.138   2.008650      1.477801
mp-aaaaabrz   Al2O3      R-3c  77.00      mp_dft 11.125792        0.138   3.699925      2.750788
```

Read it as: at room temperature a 5.98 cm sapphire crystal passes 80 % of the thermal beam; at that thickness the
uncollided 1/E epithermal flux falls by e^2.01 (a factor of 7.5) and 1 MeV photons by e^1.48. Cooled to 77 K the
crystal can be 11.1 cm thick for the same thermal transmission, and rejects more. `sigma_model` = 0.138 means the
thickness is uncertain by a factor exp(±0.138), about ±15 %; see [Concepts](concepts.md#model-uncertainty).

A crystal without a full phonon model (`vdos_source == "none"`) has empty `L08` and rejections; its
composition-only values (`proxy_L08`) and, where available, its elastic-Debye value (`debye_L08`) are still there.

## 3. Rank on one figure of merit

The ten highest epithermal rejections at 293 K among the crystals that enter the main comparison (full model,
inside the validated domain):

```python
room = at[(at["T"] == 293.15) & at.main_front_eligible]
print(room.nlargest(10, "A_fast_1E")[["formula", "sg_symbol", "tier", "L08", "A_fast_1E", "A_fast_1E_lo", "A_fast_1E_hi"]].to_string(index=False))
```

```
formula sg_symbol tier      L08  A_fast_1E  A_fast_1E_lo  A_fast_1E_hi
   Be2C     Fm-3m    B 3.815357   3.203114      2.790621      3.676491
    BeO    P6_3mc    B 4.028688   2.880998      2.509883      3.306931
    SiC    P6_3mc    B 8.111888   2.613141      2.277541      2.998057
    SiC      P3m1    B 8.053801   2.594006      2.260859      2.976111
    SiC    P6_3mc    B 8.029336   2.585225      2.253204      2.966040
    SiC    P6_3mc    B 7.972509   2.568201      2.238362      2.946516
    SiC       R3m    B 8.010773   2.544083      2.200907      2.940629
    SiC       R3m    B 7.952451   2.526025      2.185280      2.919763
    SiC     F-43m    B 7.816399   2.517599      2.194247      2.888480
    SiC       R3m    B 7.922040   2.515119      2.175844      2.907163
```

`_lo` and `_hi` are the same rejection at L08·exp(∓sigma_model): rows whose bands overlap are not distinguishable
by the model. A ranking on one axis ignores the photon rejection; [Comparing crystals](comparing-crystals.md) shows
how to rank on two axes at once.

## 4. Filter on flags

Commercially grown crystals (`tier == "A"`) that are not water-sensitive and contain none of the flagged elements:

```python
flags = ["flag_Be", "flag_Tl", "flag_As", "flag_Pb", "flag_Bi"]
safe = room[(room.tier == "A") & ~room.hygroscopic & ~room[flags].any(axis=1)]
print(len(room), len(safe))
```

```
4438 65
```

The flags are listed in [Concepts](concepts.md#domain-and-practicality-flags). `main_front_eligible` already removes
the crystals outside the validated domain (`outside_domain`); to look at those, select `at.outside_domain`.

## 5. Export a subset

```python
safe[cols].to_csv("grown_safe_293K.csv", index=False)
```

The first lines of the file:

```
material_id,formula,sg_symbol,T,vdos_source,L08,sigma_model,A_fast_1E,A_gamma_1MeV
mp-aaaaaaac,Pd,Fm-3m,293.15,sevennet,0.5212434369124594,0.14546477236774547,0.280740615324515,0.3495266013345904
mp-aaaaaaan,Fe,Im-3m,293.15,sevennet,0.6197998847457001,0.14546477236774547,0.48968022300227965,0.2982298883199846
```

## Next

- What every number means: [Atlas columns](atlas-columns.md), [Concepts](concepts.md).
- Fair comparisons and Pareto fronts: [Comparing crystals](comparing-crystals.md).
- Your own crystal: [Extending](extending.md).
