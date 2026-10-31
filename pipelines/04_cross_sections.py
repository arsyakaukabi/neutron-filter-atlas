"""04 - Cross sections of the SevenNet-phonon compounds.

Every SevenNet DOS goes through the identical conversion and NCrystal code path as the MP DFT DOS; results are cached
per compound in data/intermediate/atlas_eval (atlas assembly, step 08, reads the cache). The coherent one-phonon
check runs in step 08 because it uses the atlas rows.
"""
import _bootstrap  # noqa: F401

from multiprocessing import Pool

from filteratlas.common.pipeline import Step, run


def evaluate_sevennet():
    from filteratlas.atlas.assemble import MLIP_DIR, eval_mlip
    from filteratlas.common.constants import N_WORKERS
    files = sorted(str(p) for p in MLIP_DIR.glob("*.json.gz"))
    with Pool(N_WORKERS) as pool:
        recs = pool.map(eval_mlip, files, chunksize=4)
    print(f"evaluated {len(recs)} SevenNet DOS; errors: {sum(r['status'] != 'ok' for r in recs)}")


STEPS = [Step("evaluate_sevennet", evaluate_sevennet, "main", "SevenNet DOS -> NCrystal (cache: intermediate/atlas_eval)")]

if __name__ == "__main__":
    run(STEPS, __doc__)
