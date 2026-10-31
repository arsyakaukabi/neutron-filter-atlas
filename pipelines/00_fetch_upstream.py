"""00 - Fetch upstream data (network; Materials Project API key in .env as MP_API_KEY).

Every fetcher is cached: files already in data/upstream are not downloaded again, so this step is a no-op on a
complete data tree. Materials Project serves only its current release; a re-fetch may differ from data/MANIFEST.csv
(`make check-data` reports it). ENDF/B-VIII.0 is fetched with openmc_data_downloader.
Writes data/upstream/mp and data/upstream/exfor, plus the candidate table, the MP metadata, exfor_benchmark.csv
(EXFOR tidy table), the EXFOR summary and cache_manifest.json (file names in configs/paths.yaml).
"""
import _bootstrap  # noqa: F401

from filteratlas.common.pipeline import Step, run
from filteratlas.ingest import candidates, exfor_build, mp_manifest, structures_fc

STEPS = [
    Step("mp_candidates_and_dos", candidates.main, "network", "MP summary table, candidate set, phonon DOS (Delta table)"),
    Step("mp_structures_and_fc", structures_fc.main, "network", "MP structures (no-phonon survivors), force constants for the stability screen"),
    Step("exfor", exfor_build.main, "network", "EXFOR entries -> tidy table (intermediate/exfor_benchmark.csv)"),
    Step("mp_manifest", mp_manifest.main, "main", "SHA-256 of the frozen MP/EXFOR inputs -> derived/cache_manifest.json"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
