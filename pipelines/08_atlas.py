"""08 - Atlas assembly, fronts and checks.

Reads the outputs of every earlier step. The OpenMC transport check needs the OpenMC environment; the GPU-vs-CPU check
needs its archive-class inputs; both are skipped by default.
"""
import _bootstrap  # noqa: F401

from filteratlas.atlas import assemble, bragg, controls, fronts, missed, summary
from filteratlas.common.paths import D, F
from filteratlas.common.pipeline import Step, run
from filteratlas.dynamics import gpu_check
from filteratlas.transport import openmc_input
from filteratlas.validation import anchors_stage3
from filteratlas.xs import coherent


def _openmc():
    from filteratlas.transport import openmc_check
    work = D("logs").parent / "openmc_work"
    openmc_check.main(F("stage3_openmc_input.json"), F("stage3_openmc_check.json"), work)


STEPS = [
    Step("assemble", assemble.main, "main", "atlas.parquet, stage3_curves.parquet, atlas_counts.json (atlas part)"),
    Step("fronts", fronts.main, "main", "Pareto fronts, fidelity ladder, who moves, sensitivity, front test preview"),
    Step("controls", controls.main, "main", "positive, negative, model and drift controls"),
    Step("bragg", bragg.main, "main", "Bragg cut-off section (Tier A)"),
    Step("anchors_stage3", anchors_stage3.main, "main", "sapphire literature anchors"),
    Step("coherent", coherent.main, "main", "coherent one-phonon check, sapphire and Si"),
    Step("openmc_input", openmc_input.main, "main", "input of the OpenMC check"),
    Step("openmc_check", _openmc, "openmc", "uncollided vs transported A_fast (OpenMC 0.16)"),
    Step("missed", missed.main, "main", "what the funnel could have missed; supplementary table"),
    Step("gpu_check", gpu_check.main, "archive", "GPU vs CPU reproducibility (archive-class inputs)"),
    Step("summary", summary.main, "main", "summary numbers of the atlas"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
