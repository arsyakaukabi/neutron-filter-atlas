# FAQ

What this page is for: answers to the questions users of the atlas ask most often.
Who should read it: users of the atlas.

The snippets run from the repository root after this set-up:

```python
import pandas as pd
at = pd.read_parquet("data/derived/atlas.parquet")
```

## Why is my crystal missing?

The candidate set contains Materials Project entries that are experimentally known (not marked theoretical), lie
within 0.05 eV/atom of the convex hull, and contain no hydrogen, no radioactive element and no strong thermal
absorber (B, Cd, Gd, Sm, Eu, Li, Hg, In, Ir, Rh, Au, Hf, Dy, Er, Co, Ag). A crystal fails one of these, or the
Materials Project has no entry for it, or it is there under another formula or space group. Search loosely first:

```python
print(at[at.formula.str.contains("Al2O3")][["material_id", "formula", "sg_symbol", "T"]].drop_duplicates("material_id"))
```

The Materials Project writes reduced formulas in its own order (for example `Ge3(BiO3)4` for BGO); the space-group
symbol uses `_` for screw axes (`P6_3mc`). If the crystal is present but its `L08` is empty, it has no full phonon
model; `funnel_step` says why ([Concepts](concepts.md#coverage-which-crystals-have-a-full-phonon-model)), and
[Extending](extending.md) shows how to compute it yourself.

## Why does the proxy value differ from the full value?

`proxy_L08` uses only the free elements' absorption and incoherent scattering. It leaves out the inelastic
scattering by phonons, which is the main thermal loss of many good filter crystals at room temperature, so the proxy
attenuation is too small and `proxy_L08` too large: it is an upper bound on `L08`. The gap is largest where phonons
matter most: light elements, soft lattices, room temperature. Part of the difference is also density:
`proxy_L08` uses the Materials Project summary structure, `L08` the structure of the phonon calculation; compare
`L08` with `L08_proxy_vdos_n` (proxy at the full-model density) to see the phonon effect alone. How well the proxy
orders crystals compared with the full model is one of the results in `stage1_stats.json` and `stage3_fronts.json`
(`ladder`).

## Why do the identifiers look like `mp-aaaaabrz` and not `mp-1143`?

The atlas uses the identifiers in the form the Materials Project's current database returns (letters after `mp-`).
The Materials Project also knows the older numeric form of an entry, but the atlas does not store it. If you know a
numeric identifier, find the entry by formula and space group, as above, and compare the structure if several
entries match.

## Why are hydrogen-containing crystals excluded?

Hydrogen has an incoherent scattering cross section of about 80 barn per atom, far more than any other element. A
crystal with hydrogen loses a large part of a thermal beam within millimetres, so it cannot be a thermal filter. The
candidate rules exclude it before any computation, and a test checks that no hydrogen-bearing crystal is in the atlas.

## How do I cite one entry?

Cite the data release and the software as in [Citing](citing.md), and name the entry by `material_id`, temperature
and column, for example: "L₀.₈ = 5.98 cm (σ_model = 0.138 in ln L₀.₈) for sapphire, mp-aaaaabrz, R-3c, at 293.15 K
(neutron filter atlas data release 0.5.0)". Cite the Materials Project too.

## How does temperature change the thickness?

Cooling removes phonons, so inelastic scattering, and with it the thermal attenuation, falls; the crystal can be
thicker for the same 80 % transmission, and the rejections at that thickness rise. How much depends on the crystal:

```python
w = at.pivot_table(index=["material_id", "formula", "sg_symbol"], columns="T", values="L08").dropna()
w["ratio_77_over_293"] = w[77.0] / w[293.15]
print(w.loc["mp-aaaaabrz"])
```

```
T                       77.0   293.15  ratio_77_over_293
formula sg_symbol
Al2O3   R-3c       11.125792  5.97709           1.861406
```

A crystal whose attenuation is dominated by absorption (which does not depend on temperature) gains little from
cooling. Temperatures other than 293.15 K and 77 K need a rerun ([Extending](extending.md#3-change-the-temperature)).

## Why are there several rows for one formula?

Polymorphs and polytypes are separate Materials Project entries with their own structures and phonons (SiC, for
example, has several polytypes). Each has its own `material_id` and `sg_symbol`. The commercially grown crystals
(`tier == "A"`) are matched to one entry each, the one that best fits the vendor's stated phase (`tierA_match`).

## A measured cross-section disagrees with the atlas. What should I do?

First check that you compare like with like:

1. **Sample form.** The atlas describes a large single crystal away from Bragg reflections. A polycrystal or powder
   adds Bragg scattering and has a higher cross-section below its Bragg cut-off.
2. **Orientation.** A single crystal shows narrow Bragg dips in transmission at particular wavelengths and
   orientations; compare with the baseline between dips, not with the dips.
3. **Temperature.** Compare at the temperature of the measurement. The atlas has 293.15 K and 77 K; for other
   temperatures compute μ(λ) yourself (`filteratlas.xs.evaluate.mu_offbragg`, [Python interface](api-cross-sections.md)).
4. **Size of the disagreement.** The model uncertainty on L₀.₈ is about ±15 % (σ = 0.138 in ln L₀.₈); a difference
   in the thermal attenuation within that is consistent with the model.
5. **The data.** Impurities (especially absorbers), water, normalisation and the stated temperature of a measurement
   all matter. The validation found datasets that needed renormalisation or were superseded (column `audit` in
   `stage2_exfor_benchmark.parquet`).

If a careful single-crystal measurement disagrees by more than the model uncertainty, that is valuable: the
validation data are thin for many chemistries. Open an issue in the repository with the crystal, its temperature,
the sample form and the data or their EXFOR entry.

## Why is `melting_point_unknown` true everywhere?

The atlas holds no melting points. Melting points, growth methods and sizes are collected only for the shortlist, in
`data/curated/stage4_dossier.csv`.

## Can I use the values for a crystal at another wavelength band or moderator temperature?

Not directly from the table. For the commercially grown crystals with a full model, the attenuation curves on
0.3–4.5 Å are stored in `data/intermediate/stage3_curves.parquet`, and `filteratlas.xs.evaluate.l_target_band`
recomputes L₀.₈ for any band and moderator temperature from them; see
[Extending](extending.md#5-change-the-thermal-band-or-the-transmission-level).
