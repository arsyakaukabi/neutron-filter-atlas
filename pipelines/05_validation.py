"""05 - Validation against measured total cross-sections, and the MLIP benchmark.

Order matters: the benchmark JSON feeds the single-crystal-only and reclassification variants, the MgO audit and the
Pb/Bi diagnostic. The Pb/Bi quasi-harmonic run and the ORB diagnostic need an MLIP environment; the MLIP sensitivity
evaluation needs archive-class inputs.
"""
import _bootstrap  # noqa: F401

from filteratlas.common.pipeline import Step, run
from filteratlas.validation import (anchors_stage1, d11_robustness, exfor_benchmark, exfor_sc_only, kill2_d11,
                                    lit_compare, mgo_audit, mlip_eval, mlip_sensitivity, orb_diagnostic, pbbi_diag,
                                    pbbi_qha)

STEPS = [
    Step("anchors_stage1", anchors_stage1.main, "main", "literature anchors (Freund / Barker forms)"),
    Step("mlip_eval", mlip_eval.main, "main", "MLIP vs DFT L0.8 -> stage2_mlip.json"),
    Step("mlip_sensitivity", mlip_sensitivity.main, "archive", "MLIP sensitivity runs added after inspection of the results (archive-class inputs)"),
    Step("orb_diagnostic", orb_diagnostic.main, "mlip", "ORB failure diagnostic"),
    Step("exfor_benchmark", exfor_benchmark.run, "main", "EXFOR benchmark and acceptance criterion (audited data)"),
    Step("kill2_d11", kill2_d11.main, "main", "count restricted to single-crystal data"),
    Step("exfor_sc_only", exfor_sc_only.main, "main", "single-crystal-only variant"),
    Step("d11_robustness", d11_robustness.main, "main", "envelope, temperature, leave-one-out and no-floor variants"),
    Step("mgo_audit", mgo_audit.main, "main", "MgO 13761004 vs 13761005 (reclassification rule)"),
    Step("pbbi_qha", pbbi_qha.main, "mlip", "Pb/Bi quasi-harmonic volumes (SevenNet)"),
    Step("pbbi_diag", pbbi_diag.main, "main", "Pb/Bi VDOS-scale and QHA diagnostic"),
    Step("lit_compare", lit_compare.main, "main", "sapphire literature comparison (Mildner, Tennant, Born)"),
]

if __name__ == "__main__":
    run(STEPS, __doc__)
