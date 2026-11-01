# Comparing crystals

What this page is for: how to compare crystals fairly with the atlas, how to build a Pareto front, and what the table
cannot tell you.
Who should read it: users who choose between filter crystals; read [Concepts](concepts.md) first.

## Rules for a fair comparison

1. **Same temperature.** Compare rows with the same `T`. Cooling changes L₀.₈ and every rejection, and not by the
   same factor for every crystal.
2. **Same figure of merit.** Epithermal (`A_fast_1E`), fission-spectrum (`A_fast_watt`) and photon (`A_gamma_1MeV`,
   `A_gamma_6MeV`) rejections order crystals differently. Choose the one that matches your background, or compare on
   two at once (Pareto front, below).
3. **Within the uncertainty band.** Each rejection has a band `[_lo, _hi]` from the model uncertainty. Overlapping
   bands mean the model cannot tell the two crystals apart; see the rule of thumb in
   [Concepts](concepts.md#model-uncertainty).
4. **Inside the validated domain.** Use `main_front_eligible` (a full model and not `outside_domain`). Crystals
   outside the domain (soft lattices, soft lattices with heavy elements, molecular crystals) are kept in the table
   but their values are extrapolations; compare them only among themselves (the fronts do this in a separate panel).
5. **Same level of physics.** Compare `L08` with `L08`, not with `proxy_L08` or `debye_L08`. The proxy is an upper
   bound on L₀.₈ and can be several times too large; the elastic-Debye value is closer but still an approximation.
6. **Availability.** `tier == "A"` restricts the comparison to commercially grown crystals. A tier-B crystal may not
   exist as a large single crystal.

Practicality flags (`flag_*`, `hygroscopic`) do not change any number; filter on them according to what your
facility accepts.

## A Pareto front from the table

A crystal is on the **Pareto front** of two figures of merit if no other crystal is better on both. The package
provides the two functions the atlas uses (`filteratlas.atlas.fronts`, see [Python interface](api.md)):

- `pareto(x, y)`: True for the points not dominated on both axes (higher is better);
- `within(xhi, yhi, fx, fy)`: True for the points whose +1σ value (`_hi` columns) is not dominated by any front
  member, that is, the crystals the model cannot rule out.

The main front of the commercially grown crystals at 293 K, epithermal against 1 MeV photon rejection:

```python
import pandas as pd
from filteratlas.atlas.fronts import pareto, within

at = pd.read_parquet("data/derived/atlas.parquet")
d = at[(at["T"] == 293.15) & (at.tier == "A") & at.main_front_eligible]
on = pareto(d.A_fast_1E, d.A_gamma_1MeV)
near = within(d.A_fast_1E_hi, d.A_gamma_1MeV_hi, d.A_fast_1E[on], d.A_gamma_1MeV[on])
print(d[on][["formula", "sg_symbol", "A_fast_1E", "A_gamma_1MeV", "L08"]].to_string(index=False))
print(f"{len(d)} eligible, {on.sum()} on the front, {near.sum()} within one model uncertainty")
```

```
formula sg_symbol  A_fast_1E  A_gamma_1MeV      L08
    SiC    P6_3mc   2.459533      1.560723 7.633516
82 eligible, 1 on the front, 3 within one model uncertainty
```

A front with one member says that no other grown crystal beats it on both axes at the central values. The three
crystals within one model uncertainty include the front member itself; the other two are the crystals the model
cannot rule out. Print `d[near]` to see them.

The released front table `data/derived/stage3_fronts.parquet` holds the same result for every combination: key
`"293K|A|main|main"` is the front above (same members, same three crystals within one uncertainty). The key has four
parts: temperature (`293K`, `77K`), tier (`A`, `B`), panel (`main` = `main_front_eligible`; `outside_domain`;
`no_heavy_no_4f` = eligible without elements of Z ≥ 80 or partly filled 4f shells) and axes (`main` = `A_fast_1E`
against `A_gamma_1MeV`; `secondary` = `A_fast_watt` against `A_gamma_6MeV`). Each front holds one tier only;
there is no released front of all candidates together, but `pareto` builds one from the table in one line. The front step also leaves out
crystals flagged `unstable_dft`, which have no full model anyway.

```python
fr = pd.read_parquet("data/derived/stage3_fronts.parquet")
members = fr[(fr.key == "77K|B|main|secondary") & fr.front].material_id
print(at[at.material_id.isin(members) & (at["T"] == 77.0)][["formula", "sg_symbol", "A_fast_watt", "A_gamma_6MeV"]])
```

For the probability that one crystal beats another on a front (rather than the central values), the pipeline runs a
Monte Carlo over the model errors for the shortlist; the results are in `stage4_rank_probabilities.json` and the
settings under `rank_mc` in [Configuration](configuration.md#shortlistyaml).

## What the table cannot tell

| Question | What the atlas assumes | Where the package helps |
|---|---|---|
| **Orientation and Bragg dips** | the crystal is oriented away from every Bragg reflection over the whole band (off-Bragg baseline) | for the commercially grown shortlist, step 10 scans 500 beam directions for Bragg dips and finds the best orientation (`stage4_dips.json`; settings `dips` in `configs/shortlist.yaml`) |
| **Mosaic spread** | not in the off-Bragg baseline | the dip scan uses a Gaussian mosaic of 0.3°, with 0.1° and 1.0° as sensitivity |
| **Beam divergence and in-scattering** | narrow beam: only uncollided particles count | slab transport with OpenMC for the shortlist counts every transmitted particle (`stage4_transport.json`); nothing for other crystals |
| **Polycrystalline (Bragg cut-off) filters** | single crystals only | the cold pass band of polycrystalline forms of the commercially grown crystals (`stage3_bragg.json`) |
| **Other temperatures** | 293.15 K and 77 K only | a series 293-77 K for the commercially grown shortlist (`stage4_temperature.json`); rerun for others ([Extending](extending.md#3-change-the-temperature)) |
| **Manufacturability, size, price** | not modelled | vendor offers for the commercially grown crystals (`data/curated/growable_crystals.csv`); melting point, growth method, size, water sensitivity and price for the shortlist (`data/curated/stage4_dossier.csv`, `stage4_dossier.json`) |
| **Activation** | not in the atlas | activity after irradiation for the shortlist (step 10, OpenMC depletion; in `stage4_dossier.json`) |
| **Impurities, defects, radiation damage** | ideal, pure crystals | not computed |
| **Crystals without a full model** | no rejections | elastic-Debye and upper-bound positions in `stage3_not_computed_table.csv`; compute the crystal yourself ([Extending](extending.md)) |

The atlas ranks candidates. Before buying a crystal, check the ranking with a measurement of its total
cross-section or its transmission, ideally on the actual crystal.
