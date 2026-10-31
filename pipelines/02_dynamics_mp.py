"""02 - Materials Project phonons: VDOS models, triviality-test statistics and dynamical-stability screens.

The MP phonon DOS of every candidate is converted to per-element VDOS (2 meV floor), screened with the DOS rule
and evaluated with NCrystal (L0.8, mu, curves, proxy, Debye nulls). The usable set and the statistics of the
triviality test (proxy vs full model) follow, then the force-constant stability screens that need them. Heavy steps: mp_models (about 1-2 h with 4 workers) and the
force-constant screens (phonopy on about 2,500 compounds).
"""
import _bootstrap  # noqa: F401

from filteratlas.atlas import enantiomorphs, fidelity
from filteratlas.common.pipeline import Step, run
from filteratlas.dynamics import (dimensionality, fc_validation, floor_sensitivity, reproducibility, screen,
                                  stability_screen, stability_summary)
from filteratlas.xs import model_control, mp_models

STEPS = [
    Step("model_control", model_control.main, "main", "Debye-300 K fallback on sapphire vs EXFOR 31743 (model control)"),
    Step("mp_models", mp_models.main, "main", "MP DOS -> NCrystal models, fidelity quantities and controls"),
    Step("dimensionality", dimensionality.main, "main", "molecular-crystal flag input"),
    Step("fidelity", fidelity.main, "main", "usable set and fidelity statistics (proxy vs full, nulls, sensitivity)"),
    Step("enantiomorphs", enantiomorphs.main, "main", "enantiomorph pairs as a reproducibility check"),
    Step("fc_validation", fc_validation.main, "main", "2 meV floor vs force-constant tetrahedron DOS"),
    Step("floor_sensitivity", floor_sensitivity.main, "main", "fidelity statistics at 1 and 3 meV floors"),
    Step("reproducibility", reproducibility.main, "main", "pheasy vs DFPT, SIMEX"),
    Step("stability_screen_sample", stability_screen.run, "main", "dense-mesh stability screen, stratified sample"),
    Step("stability_screen_all", stability_screen.run, "main", "dense-mesh stability screen, all usable compounds with force constants",
         {"all_usable": True}),
    Step("stability_screen_mesh60", stability_screen.run, "main", "mesh-density sensitivity",
         {"all_usable": True, "mesh": 60.0}),
    Step("stability_summary", stability_summary.main, "main", "prevalence, restricted reproducibility, large cells"),
    Step("d7a", screen.d7a, "main", "stability screen of the funnel survivors with MP phonons"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
