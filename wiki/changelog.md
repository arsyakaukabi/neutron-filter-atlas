# Changelog

What this page is for: the released versions and what changed for the user in each.
Who should read it: users who compare results across versions or update a clone.

## 0.5.0 (2026-10-31): first public release

- The atlas: 18,449 candidate crystals of the Materials Project at 293.15 K and 77 K, with L₀.₈, the epithermal,
  fission-spectrum and 1 MeV and 6 MeV photon rejections, the model uncertainty and the domain and practicality
  flags ([Atlas columns](atlas-columns.md)).
- Full phonon models from Materials Project DFT phonons and from SevenNet phonons computed for the atlas; composition
  and elastic-Debye levels for every crystal that has the data.
- Validation against measured total cross-sections and a model uncertainty derived from it.
- Pareto fronts, Monte Carlo rank probabilities, slab transport, Bragg dips, temperature series, activation and
  practicality data for a shortlist of commercially grown crystals.
- The pipeline (steps 00-10), configuration, tests and the data manifest; the data release of the same version
  ([Data](data.md)).

The repository's `CHANGELOG.md` records the version number only; this page describes what it contains. The pages in
`wiki/` describe version 0.5.0.
