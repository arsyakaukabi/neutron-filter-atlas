"""Model, positive and negative controls of the atlas, run as tests."""
import json

import numpy as np
import pytest

from conftest import DATA, FILES, cfg, need, need_dir


def test_model_control_debye300_sapphire():
    """The Debye-300 K fallback on sapphire overestimates EXFOR 31743 by the reference factor 2.85."""
    need("exfor_benchmark.csv", "stage0_model_control.json")
    need_dir("mp")
    from filteratlas.atlas.controls import model_control
    r = model_control()
    mc = cfg("controls")["model_control"]
    assert r["median_ratio_now"] == pytest.approx(mc["phase0_reference_ratio"], rel=mc["tolerance"])
    assert r["pass"]


def test_proxy_is_a_bound(atlas):
    """mu_proxy <= mu_th within the tolerance of stage1_stats.json: L0.8,full / L0.8,proxy (same density) never exceeds it."""
    need("stage1_stats.json")
    tol = json.loads((DATA / FILES["stage1_stats.json"]).read_text())["proxy_lower_bound"]
    full = atlas[atlas.vdos_source.ne("none")]
    for T, tag in ((293.15, "293K"), (77.0, "77K")):
        d = full[full["T"] == T]
        r = d.L08 / d.L08_proxy_vdos_n
        assert r.max() <= tol[f"L08_{tag}"]["max_ratio_full_over_proxy"] + 5e-5, (tag, d.loc[r.idxmax(), "formula"])


def test_positive_controls_reproduce_benchmark(atlas):
    """Sapphire, Si, MgO and quartz through the atlas pipeline reproduce their benchmark single-crystal R_c to 3 decimals."""
    need("exfor_benchmark.csv", "stage2_exfor_benchmark.parquet", "stage2_d11_robustness.json", "stage2_text_numbers.json")
    from filteratlas.atlas.controls import positive
    res = positive(atlas)
    for f in cfg("controls")["stage3"]["positive"]:
        assert res[f]["R_c_match_3dp"], (f, res[f]["R_c_sc_atlas"], res[f]["R_c_sc_stage2"])
    assert res["pass"]


def test_negative_controls_and_no_hydrogen(atlas):
    """CdS, natural-Li LiF and Gd2O3 rank below the 1st percentile of the atlas L0.8; no H-bearing compound."""
    need("stage1_controls.parquet")
    from filteratlas.atlas.controls import negative
    res = negative(atlas)
    assert res["n_H_bearing_in_atlas"] == 0
    assert res["pass"]


def test_atlas_counts_consistent(atlas, derived_json):
    c = derived_json("atlas_counts.json")["atlas"]
    one = atlas[atlas["T"] > 100]
    assert c["compounds"] == one.material_id.nunique() == len(one)
    assert c["by_funnel_step"] == one.funnel_step.value_counts().to_dict()
    sig = atlas[atlas.vdos_source.ne("none")].sigma_model
    assert np.all(sig >= cfg("uncertainty")["sigma_model_dft"])
