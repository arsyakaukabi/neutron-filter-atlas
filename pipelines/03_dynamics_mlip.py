"""03 - Machine-learned phonons: which compounds, and the runs.

MLIP benchmark: the benchmark set, then one run per model in its own environment (requirements-mlip*.txt):
    python -m filteratlas.dynamics.mlip_run <model> [--threads N] [--device cpu|cuda]
Atlas: Debye scores and supercell sizes for the survivors without MP phonons, the funnel that selects the
SevenNet set, then the SevenNet runs through the scheduler (one process per compound):
    python -m filteratlas.dynamics.scheduler        (machine settings: configs/compute.yaml, section scheduler)
The run steps are compute-heavy (hours to days).
"""
import _bootstrap  # noqa: F401

from filteratlas.atlas import funnel
from filteratlas.common.pipeline import Step, run
from filteratlas.dynamics import screen
from filteratlas.validation import mlip_select


def _benchmark_runs():
    raise SystemExit("run `python -m filteratlas.dynamics.mlip_run <model>` per model in its MLIP environment")


def _atlas_runs():
    from filteratlas.dynamics import scheduler
    scheduler.main()


STEPS = [
    Step("mlip_benchmark_set", mlip_select.main, "main", "learned-potential benchmark set (stratified random + EXFOR crystals)"),
    Step("mlip_benchmark_runs", _benchmark_runs, "mlip", "per-model MLIP phonons of the benchmark set"),
    Step("debye_scores", screen.debye, "main", "Debye L0.8 and SevenNet supercell size"),
    Step("funnel", funnel.main, "main", "SevenNet budget selection and drift set"),
    Step("mlip_atlas_runs", _atlas_runs, "mlip", "SevenNet atlas phonons via the scheduler"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
