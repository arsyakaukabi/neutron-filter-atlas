"""01 - Candidate set, proxy bound and tiers.

Reads data/upstream/mp and data/curated/growable_crystals.csv; offline when data/upstream/mp is complete.
Writes the candidate table, the MP metadata, the proxy bound and the tier table (file names in configs/paths.yaml).
"""
import _bootstrap  # noqa: F401

from filteratlas.atlas import tiers
from filteratlas.common.pipeline import Step, run
from filteratlas.ingest import candidates
from filteratlas.xs import proxy_bound

STEPS = [
    Step("candidates", candidates.main, "main", "element filters, E_hull cut; offline if the MP cache is complete"),
    Step("proxy_bound", proxy_bound.main, "main", "composition-only bound; excluded_by_proxy_bound"),
    Step("tiers", tiers.main, "main", "Tier A = growable_crystals.csv matched to MP entries"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
