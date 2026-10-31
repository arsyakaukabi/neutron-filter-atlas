"""10 - The shortlist.

Reads the atlas (step 08); never the 2,000-step relaxation reruns of step 09. The OpenMC steps (slab transport,
activation) need the OpenMC environment and are skipped by default; run them there as the module docstrings say, then the "main" steps
that read their outputs.
"""
import _bootstrap  # noqa: F401

from filteratlas.common.paths import D
from filteratlas.common.pipeline import Step, run
from filteratlas.shortlist import dips, rank_mc, select, temperature, transport_inputs, transport_results


def _transport():
    from filteratlas.transport import slab_transport
    slab_transport.main(D("stage4_transport"), D("stage4_transport") / "work")


def _activation_prepare():
    from filteratlas.transport import activation
    activation.prepare()


def _activation_run():
    from filteratlas.transport import activation
    activation.run(D("stage4_activation"), D("stage4_activation") / "work")


def _dossier():
    from filteratlas.shortlist import dossier
    dossier.main()


def _summary():
    from filteratlas.shortlist import summary4
    summary4.main()


def _controls():
    from filteratlas.shortlist import controls
    controls.main()


STEPS = [
    Step("select", select.main, "main", "the shortlist and the reason for each entry"),
    Step("rank_mc", rank_mc.main, "main", "rank probabilities, front test, sensitivities"),
    Step("transport_inputs", transport_inputs.main, "main", "NCMAT files and slab cases for OpenMC"),
    Step("transport", _transport, "openmc", "energy-resolved slab transport"),
    Step("transport_results", transport_results.main, "main", "in-scatter corrections, effect on fronts and statements, framing"),
    Step("dips", dips.main, "main", "Bragg-dip density and best orientation, Tier A"),
    Step("temperature", temperature.main, "main", "temperature series 293-77 K, Tier A"),
    Step("activation_prepare", _activation_prepare, "main", "activation inputs"),
    Step("activation", _activation_run, "openmc", "OpenMC depletion, 7 d at 1e8 n/cm2/s"),
    Step("dossier", _dossier, "main", "practicality dossier"),
    Step("controls", _controls, "main", "controls after transport and literature anchors"),
    Step("summary", _summary, "main", "summary numbers of the shortlist"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
