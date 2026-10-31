"""Config values equal the values fixed for the analysis, and parse as numbers."""
import re
from pathlib import Path

import pytest

from conftest import REPO, cfg, need


def walk(x, path=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from walk(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from walk(v, f"{path}[{i}]")
    else:
        yield path, x


@pytest.mark.parametrize("name", sorted(p.stem for p in (REPO / "configs").glob("*.yaml")))
def test_no_number_parsed_as_string(name):
    """A YAML number that loads as a string (e.g. 1e5 under YAML 1.1) would silently change behaviour."""
    num = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
    bad = [(p, v) for p, v in walk(cfg(name)) if isinstance(v, str) and num.match(v)
           and p not in ("/scheduler/gpu_memory_fraction", "/gpu_memory_fraction")]
    assert not bad, bad


def test_vdos_floor():
    """MP DOS bins below 2 meV are dropped."""
    assert cfg("vdos")["vdos_floor_meV"] == 2.0


def test_ncrystal_minimum_energy():
    assert cfg("vdos")["ncrystal_emin_eV"] == 1e-5


def test_stability_rule():
    s = cfg("stability")
    assert s["unstable_freq_THz"] == -0.5 and s["unstable_weight_max"] == 0.005   # DOS rule
    assert s["fc_unstable_THz"] == -0.5 and s["fc_mesh_length_A"] == 40.0       # force-constant rule


def test_single_crystal_forms():
    forms = set(cfg("benchmark")["single_crystal_forms"])
    assert forms == {"single_crystal", "single_crystal_assumed", "pyrolytic_graphite"}   # pyrolytic graphite counts


def test_reclassification_factor():
    assert cfg("benchmark")["d12_reclassification"]["factor"] == 2.0


def test_sigma_values():
    u = cfg("uncertainty")
    assert u["sigma_model_dft"] == 0.138
    assert u["sigma_mlip_add"] == {"293K": 0.046, "77K": 0.040}
    assert u["sigma_soft_floor"] == 0.144


def test_sigma_matches_uncertainty_model(derived_json):
    """configs/uncertainty.yaml was fixed from stage2_uncertainty_model.json; the two must agree."""
    m = derived_json("stage2_uncertainty_model.json")
    u = cfg("uncertainty")
    assert round(m["band_DFT_constructions"]["B_per_crystal_rms_then_rms"]["value"], 3) == u["sigma_model_dft"]
    assert round(m["band_DFT_constructions"]["D_within_band_bins"]["value"], 3) == u["sigma_construction_D"]
    assert round(m["heavy_soft_metal_abs_max"], 3) == u["sigma_soft_floor"]
    for t in ("293K", "77K"):
        assert round(m["mlip_vs_dft_dlnL08"][t]["p68_abs_dlog"], 3) == u["sigma_mlip_add"][t]


def test_relax_sensitivity_is_separate():
    """The extended relaxation is a sensitivity run in its own folder; the limit of the atlas stays 500 steps."""
    m = cfg("mlip")
    assert m["relax_max_steps"] == 500 and m["relax_sensitivity_max_steps"] == 2000
    dirs = cfg("paths")["dirs"]
    assert dirs["mlip_relax_sensitivity"] != dirs["mlip_atlas"]


def test_atlas_never_reads_the_sensitivity_run():
    """Only the MLIP runner (writes) and step 09 (reads) may touch the 2,000-step outputs; steps 04-08 do not."""
    users = {p.relative_to(REPO / "src").as_posix() for p in (REPO / "src" / "filteratlas").rglob("*.py")
             if "relax_sensitivity" in p.read_text(encoding="utf-8")}
    assert users == {"filteratlas/dynamics/mlip_atlas_run.py", "filteratlas/atlas/relax_sensitivity.py"}, users
    for n in ("04", "05", "06", "07", "08"):
        src = next((REPO / "pipelines").glob(f"{n}_*.py")).read_text(encoding="utf-8")
        assert "relax_sensitivity" not in src, n


def test_lambda_T():
    """lambda_T = 1.77889 A is the fixed constant, stated as the value of h / sqrt(2 m_n k_B 300 K). With CODATA 2018
    constants the formula gives 1.778856 A (relative difference 1.9e-5); the pipeline uses the fixed constant, and no
    result depends on the difference at this level."""
    from filteratlas.xs.evaluate import lambda_T
    val = cfg("figures_of_merit")["lambda_T_A"]
    assert val == 1.77889
    assert abs(lambda_T(300.0) - val) / val < 1e-4


def test_every_data_file_used_by_code_is_registered():
    names = set()
    for p in (REPO / "src" / "filteratlas").rglob("*.py"):
        names |= set(re.findall(r'\bF\("([^"{}]+)"\)', p.read_text(encoding="utf-8")))
    missing = names - set(cfg("paths")["files"])
    assert not missing, missing
