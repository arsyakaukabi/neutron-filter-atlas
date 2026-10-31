"""07 - Fast-neutron and photon data for the rejection figures of merit.

Extraction of the elemental ENDF/B-VIII.0 tables needs OpenMC (openmc.data); run it in the OpenMC environment:
    python pipelines/07_figures_of_merit.py --with openmc
The grid-convergence check then runs in the main environment.
"""
import _bootstrap  # noqa: F401

from filteratlas.common.paths import D
from filteratlas.common.pipeline import Step, run
from filteratlas.fom import nucdata_check


def _extract():
    from filteratlas.ingest import nucdata
    nucdata.main(D("endfb80_hdf5"), D("nucdata"))


STEPS = [
    Step("nucdata_extract", _extract, "openmc", "ENDF/B-VIII.0 HDF5 -> elemental tables (intermediate/nucdata)"),
    Step("nucdata_check", nucdata_check.main, "main", "5,000-points-per-decade grid vs exact union grid (< 0.1 %)"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
