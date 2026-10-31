"""09 - Relaxation-limit sensitivity (added after inspection of the results).

A sensitivity rerun of failed relaxations with 2,000 steps, kept separate from the atlas, which uses the 500-step rule.
Reads the main atlas (step 08) and the 2,000-step reruns in data/intermediate/mlip_relax_sensitivity/; writes only
the relaxation-sensitivity JSON in data/derived. Kept outside steps 04-08 so that no atlas code path reads the
sensitivity outputs.
"""
import _bootstrap  # noqa: F401

from filteratlas.atlas import relax_sensitivity
from filteratlas.common.pipeline import Step, run

STEPS = [Step("relax_sensitivity", relax_sensitivity.main, "main",
              "2,000-step reruns: convergence, dln L0.8 vs 500 steps, fronts and top-ranking tau with them added")]

if __name__ == "__main__":
    run(STEPS, __doc__)
