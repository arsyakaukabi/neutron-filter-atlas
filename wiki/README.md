# neutron-filter-atlas documentation

What this page is for: the index of the user documentation.
Who should read it: everyone; start here.

The atlas lists, for every candidate crystal of the Materials Project, how thick a single-crystal filter must be to
pass 80 % of a thermal neutron beam (L₀.₈), and how strongly that filter then rejects epithermal neutrons, fission
neutrons and 1 MeV and 6 MeV photons, at 293 K and 77 K, each with a model uncertainty and flags for where the model
is not validated. If you only want to read the atlas, you need the results archive of the data release and the pages
under *Reading the atlas*; you never have to run the pipeline. If you want to recompute the atlas, change its settings
or add your own crystals, read *Running the pipeline* as well.

These pages are plain Markdown in the repository: they travel with a clone and match the code of the same version.

## Getting started

| Page | What it covers |
|---|---|
| [Installation](installation.md) | Python and operating system, installing from a clone, optional dependencies, checking the installation |
| [Data](data.md) | The data release, its four archives, where to put them, how to verify them, what runs without them |
| [Quickstart](quickstart.md) | Five minutes with the atlas: load, look up, rank, filter, export |

## Reading the atlas

| Page | What it covers |
|---|---|
| [Atlas columns](atlas-columns.md) | The column dictionary of the atlas table: meaning, unit, definition, range, missing values |
| [Results tables](results-tables.md) | The columns of the other released tables (candidates, validation, fronts, coverage, shortlist) |
| [Curated tables](curated-tables.md) | The hand-made tables kept in the repository (vendor offers, water-sensitive families, EXFOR entries, practicality) |
| [Concepts](concepts.md) | Definitions behind the numbers: thermal band, L₀.₈, rejections, levels of phonon description, model uncertainty, flags |
| [Comparing crystals](comparing-crystals.md) | How to compare two crystals fairly, build a Pareto front, and what the table cannot tell |
| [FAQ](faq.md) | Missing crystals, proxy against full values, hydrogen, citing one entry, temperature, disagreeing measurements |

## Running the pipeline

| Page | What it covers |
|---|---|
| [Pipeline overview](pipeline-overview.md) | The steps from structure to atlas entry, their commands, inputs, outputs and run times |
| [Configuration](configuration.md) | Every configuration key: value, meaning, which step reads it |
| [Reproducing](reproducing.md) | Recomputing the released atlas and checking the result against the release |
| [Extending](extending.md) | Recipes: your own crystal or phonon DOS, another temperature, spectrum, band or transmission level, one step only |
| [Python interface](api.md) | The public functions: configuration, data paths, phonon files, Pareto fronts |
| [Python interface: materials](api-materials.md) | Building NCrystal materials from a phonon DOS or a Debye temperature |
| [Python interface: cross sections](api-cross-sections.md) | Off-Bragg attenuation, L₀.₈, rejections |
| [Troubleshooting](troubleshooting.md) | Error messages you can meet, their cause and the fix |

## Reference

| Page | What it covers |
|---|---|
| [Glossary](glossary.md) | One-line definitions of the technical terms used here |
| [Citing](citing.md) | How to cite the software and the data release, and the upstream sources to cite with them |
| [Changelog](changelog.md) | Released versions and what changed for the user |

## Keeping the generated pages current

Four groups of pages are generated from the code and the data and must not be edited by hand: the column pages
(`tools/docs/gen_columns.py`, needs the released data), [Configuration](configuration.md) (`tools/docs/gen_config.py`)
and the three Python interface pages (`tools/docs/gen_api.py`). Each script has a `--check` option that fails when a
page is out of date. `tools/docs/check_docs.py` checks links, anchors, paths and page form for every page.
