"""The fixed constants under the names the pipeline code has always used, read from configs/*.yaml.

This module only renames: every value is defined (with a comment on its origin) in the YAML file
named next to it. Lists become tuples and element lists become sets, exactly as in the original constants.
"""
import os

from dotenv import load_dotenv

from .config import cfg, tup
from .paths import env_file

load_dotenv(env_file())
MP_KEY = os.environ.get("MP_API_KEY", "")
os.environ.setdefault("MATERIALSPROJECT_USER_API_KEY", MP_KEY)  # NCrystal reads this name

_c, _v, _s, _f, _a = cfg("candidates"), cfg("vdos"), cfg("stability"), cfg("figures_of_merit"), cfg("analysis")

# candidates.yaml
E_HULL_MAX = _c["e_hull_max_eV_per_atom"]
RADIOACTIVE = set(_c["radioactive"])
STRONG_ABSORBERS = set(_c["strong_absorbers"])
FLAG_ELEMENTS = set(_c["flag_elements"])
# stability.yaml
UNSTABLE_FREQ_THZ = _s["unstable_freq_THz"]
UNSTABLE_WEIGHT_MAX = _s["unstable_weight_max"]
# vdos.yaml
THZ_TO_EV = _v["thz_to_eV"]
VDOS_NPTS = _v["vdos_npts"]
# figures_of_merit.yaml
TEMPS = tup(_f["temperatures_K"])
BAND = tup(_f["band_A"])
HOT_BAND = tup(_f["hot_band_A"])
MODERATOR_T = _f["moderator_T_K"]
LAMBDA_T = _f["lambda_T_A"]
N_LAMBDA = _f["n_lambda"]
T_TARGET = _f["T_target"]
PROBE_WL = tup(_f["probe_wavelengths_A"])
# analysis.yaml
USABLE_MU_CUT = _a["stage1"]["usable_mu_cut_per_cm"]
USABLE_SENSITIVITY = tup(_a["stage1"]["usable_mu_cut_sensitivity"])
KILL_RHO = _a["stage1"]["kill_rho"]
N_BOOT = _a["n_boot"]
N_WORKERS = int(os.environ.get("FA_WORKERS", cfg("compute")["n_workers"]))
