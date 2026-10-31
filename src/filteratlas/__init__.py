"""filteratlas: the pipeline of the single-crystal neutron filter atlas.

Packages:
- common:     configuration (configs/*.yaml), data paths (configs/paths.yaml), constants, hashing, step runner
- ingest:     Materials Project, EXFOR and ENDF/B-VIII.0 fetch and parsing
- dynamics:   phonon DOS stability screens, MLIP phonons, the SevenNet scheduler
- xs:         NCrystal materials from a VDOS, off-Bragg cross sections, proxy, Debye model, coherent check
- fom:        fast-neutron and photon rejection, nuclear-data grid check
- atlas:      tiers, funnel, atlas assembly, Pareto fronts, controls, fidelity statistics
- validation: EXFOR benchmark, MLIP benchmark, residuals and the uncertainty model, literature anchors
- transport:  OpenMC slab transport and activation
- shortlist:  shortlist selection, rank Monte Carlo, Bragg dips, temperature series, practicality
Pipeline entry points: pipelines/00_fetch_upstream.py ... pipelines/10_shortlist.py.
"""

__version__ = "0.5.0"
