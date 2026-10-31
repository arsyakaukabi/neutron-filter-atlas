"""06 - Model uncertainty (residual classes, band constructions A-D, MLIP term) and the validation summary numbers.

The atlas sigma values in configs/uncertainty.yaml were fixed from the uncertainty-model JSON written by the residuals
step; tests/test_configs.py checks that the two still agree.
"""
import _bootstrap  # noqa: F401

from filteratlas.common.pipeline import Step, run
from filteratlas.validation import residuals, summary

STEPS = [
    Step("residuals", residuals.main, "main", "residual classes and the uncertainty model"),
    Step("stage2_summary", summary.main, "main", "summary numbers of the validation and the uncertainty model"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
