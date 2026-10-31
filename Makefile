# neutron-filter-atlas: one target per pipeline step, plus data checks and tests.
# Run from the repository root in the activated main environment (environment.yml). PYTHON can be overridden,
# e.g. `make atlas PYTHON=.venv/Scripts/python.exe`. Targets are listed in data-flow order.

PYTHON ?= python
RUN = PYTHONIOENCODING=utf-8 $(PYTHON) -W ignore

.PHONY: help fetch candidates dynamics-mp dynamics-mlip cross-sections validation uncertainty figures-of-merit \
        atlas relax-sensitivity shortlist all check-data check-data-fast manifest test

help:
	@echo "fetch | candidates | dynamics-mp | dynamics-mlip | cross-sections | validation | uncertainty |"
	@echo "figures-of-merit | atlas | relax-sensitivity | shortlist | all | check-data | check-data-fast | manifest | test"

fetch:             ## 00: download upstream data (network, MP_API_KEY in .env)
	$(RUN) pipelines/00_fetch_upstream.py --with network

candidates:        ## 01: candidate set, proxy bound, tiers
	$(RUN) pipelines/01_candidates.py

dynamics-mp:       ## 02: MP phonons -> VDOS models, fidelity statistics, stability screens
	$(RUN) pipelines/02_dynamics_mp.py

dynamics-mlip:     ## 03: learned-potential selection steps (runs need an MLIP environment: add --with mlip)
	$(RUN) pipelines/03_dynamics_mlip.py

cross-sections:    ## 04: SevenNet DOS -> NCrystal
	$(RUN) pipelines/04_cross_sections.py

validation:        ## 05: EXFOR benchmark, learned-potential benchmark
	$(RUN) pipelines/05_validation.py

uncertainty:       ## 06: uncertainty model
	$(RUN) pipelines/06_uncertainty.py

figures-of-merit:  ## 07: nuclear-data grid check (extraction needs OpenMC: add --with openmc)
	$(RUN) pipelines/07_figures_of_merit.py

atlas:             ## 08: atlas, fronts, controls, checks
	$(RUN) pipelines/08_atlas.py

relax-sensitivity: ## 09: sensitivity of the 500-step relaxation limit (reads the atlas; not used by it)
	$(RUN) pipelines/09_relax_sensitivity.py

shortlist:         ## 10: shortlist, rank probabilities, dips, temperature (transport and activation need OpenMC: add --with openmc)
	$(RUN) pipelines/10_shortlist.py

all: candidates dynamics-mp dynamics-mlip cross-sections validation uncertainty figures-of-merit atlas relax-sensitivity shortlist

check-data:        ## verify data/ against data/MANIFEST.csv (sizes and SHA-256)
	$(RUN) scripts/check_data.py

check-data-fast:   ## sizes only
	$(RUN) scripts/check_data.py --fast

manifest:          ## regenerate data/MANIFEST.csv (before distributing data)
	$(RUN) scripts/build_manifest.py

test:
	$(RUN) -m pytest -q tests
