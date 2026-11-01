# Extending

What this page is for: recipes for computing your own crystals and changing the settings of the atlas.
Who should read it: pipeline users; you need the installation ([Installation](installation.md)) and, for most
recipes, the intermediate archive ([Data](data.md)).

Two ways to extend: compute **single crystals** with the Python interface (fast, no rerun of the pipeline; recipes
1-5), or change **the whole atlas** through `configs/` and rerun steps (recipe 6 and the notes in each recipe). Run
the snippets from the repository root with the package installed (`pip install -e .`). Every snippet was run on the
released data and the outputs shown are real, except the call to the SevenNet potential in recipe 1, which needs its
own environment (the same code was run with ASE's EMT potential instead). The snippets print some `UserWarning`s from
pymatgen and NCrystal; they are harmless, and `python -W ignore` hides them.

## 1. Add one crystal from a structure file

The phonons come from a learned potential, which runs in its own environment; the cross sections need NCrystal,
which is in the main environment. So the recipe has two halves, joined by a DOS file.

**In the learned-potential environment** (install it as in [Installation](installation.md#optional-dependencies),
then `pip install -e .` there too, so that `filteratlas` imports): relax the structure, compute the phonon DOS by
finite displacements with the functions of the atlas (`filteratlas.dynamics.mlip_run`, see
[Python interface: materials](api-materials.md)), check stability and save the DOS.

```python
import gzip, json
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from filteratlas.dynamics.mlip_run import calculator, phonons, relax
from filteratlas.xs.ncmat import stability

calc, _ = calculator("sevennet_mf_ompa")                 # the potential of the atlas
s = Structure.from_file("my_crystal.cif").get_primitive_structure()
atoms, converged, steps, fmax = relax(AseAtomsAdaptor.get_atoms(s), calc)
n = [max(1, int(np.ceil(12.0 / L))) for L in atoms.cell.lengths()]          # supercell of at least 12 Å per axis
dos, fmin_mesh, fmin_comm, ndisp = phonons(list(atoms.get_chemical_symbols()), atoms.cell[:],
                                           atoms.get_scaled_positions(), np.diag(n), np.eye(3), calc=calc)
print(f"relaxed: {converged} in {steps} steps; lowest frequency {fmin_comm:.3f} THz; {ndisp} displaced supercells")
print(stability(dos))                                    # (unstable DOS weight, passes the DOS rule)
with gzip.open("my_crystal_dos.json.gz", "wt", encoding="utf-8") as fh:
    json.dump(dos, fh)
```

Accept the DOS only if `converged` is True (the atlas treats a relaxation that does not converge within 500 steps as
failed), the lowest commensurate frequency is not below −0.5 THz and `stability(dos)[1]` is True: the
dynamical-stability rules of the atlas.

**In the main environment:** read the DOS and compute the crystal's row with the helper below, which computes every
value-bearing column of the atlas the way the atlas does (it reproduces the released rows exactly; for sapphire at
77 K the largest relative difference is 0).

```python
import numpy as np
from filteratlas.common.config import cfg
from filteratlas.fom import rejection as rj
from filteratlas.ingest.mp import load_dos
from filteratlas.xs import evaluate as xe, ncmat
from filteratlas.xs.mp_models import mean_phonon

U = cfg("uncertainty")


def atlas_row(dos, formula, name, T=293.15, source="sevennet"):
    """The atlas columns of one crystal at one temperature, computed as the atlas computes them."""
    tag = "293K" if T > 100 else "77K"
    r = xe.evaluate(ncmat.build_vdos_material(dos, name), (T,))
    n, L = r["numberdensity_A3"], r["L08_" + tag]
    mph = mean_phonon(dos)
    s = U["sigma_model_dft"] if source == "mp_dft" else float(np.hypot(U["sigma_model_dft"], U["sigma_mlip_add"][tag]))
    if mph < U["soft_lattice_mean_phonon_meV"]:
        s = max(s, U["sigma_soft_floor"])
    Ls = np.array([L, L * np.exp(-s), L * np.exp(s)])
    sig = rj.sigma_tot(formula)
    a1, aw = rj.a_fast(sig, n, Ls, "1E"), rj.a_fast(sig, n, Ls, "watt")
    g1, g6 = rj.mu_gamma(formula, n)
    row = {"T": T, "L08": L, "L08_hot": r["L08_hot_" + tag], "sigma_model": s, "mean_phonon_meV": mph,
           "numberdensity_A3": n, "mu_1.0": r["mu_1.0_" + tag], "mu_1.8": r["mu_1.8_" + tag],
           "L08_proxy_vdos_n": r["L08_proxy_" + tag], "mu_epi_1eV": rj.mu_epi(formula, n)}
    for k, a in (("A_fast_1E", a1), ("A_fast_watt", aw), ("A_gamma_1MeV", g1 * Ls), ("A_gamma_6MeV", g6 * Ls)):
        row.update({k: a[0], k + "_lo": a[1], k + "_hi": a[2]})
    row["nieman_ratio"] = row["mu_epi_1eV"] / row["mu_1.8"]
    return row


dos = load_dos("my_crystal_dos.json.gz")
print(atlas_row(dos, "Al", "my_crystal", 293.15))
```

`formula` must be the crystal's formula; `source="mp_dft"` gives the uncertainty of a DFT DOS. Set the flags
yourself from the row: `soft_lattice` if `mean_phonon_meV` < 10, `heavy_soft` and `untested_heavy` from the largest
atomic number ([Concepts](concepts.md#domain-and-practicality-flags)). `rejection` needs the intermediate archive
(nuclear-data tables).

**Trying it without a learned potential.** Any ASE calculator works in `relax` and `phonons`, so the whole recipe can
run in the main environment with ASE's built-in EMT potential, which knows a few metals (not accurate enough for a
real filter study): replace `calculator("sevennet_mf_ompa")` by `from ase.calculators.emt import EMT; calc = EMT()`.
For an fcc aluminium CIF (a = 4.05 Å) the first half prints

```
relaxed: True in 4 steps; lowest frequency 2.084 THz; 1 displaced supercells
(0.0, True)
```

and the row has `L08` = 4.889 cm at 293 K (8.724 cm at 77 K), `A_fast_1E` = 0.462 and `A_gamma_1MeV` = 0.843.

## 2. Supply your own phonon density of states

If you already have phonons (from DFT, another code, or a measurement), put them in the pipeline's DOS format: a dict
with `frequencies` (THz, one grid), `densities` (total DOS on that grid), `projected_densities` (one row per atom of
`structure`, on the same grid) and `structure` (pymatgen `Structure.as_dict()`). Units of the densities do not matter;
each element's DOS is normalised. From phonopy force constants, `phonons(..., fc=force_constants)` produces this dict
without any calculator.

The Materials Project DOS files in the release have this format; sapphire through the same path reproduces its atlas
row exactly:

```python
import pandas as pd
from filteratlas.ingest.mp import load_dos
from filteratlas.xs import ncmat, evaluate as xe
from filteratlas.fom import rejection as rj

dos = load_dos("data/upstream/mp/phonon_dos/pheasy/aaaegyio.json.gz")   # or your own dict
name = ncmat.build_vdos_material(dos, "my_sapphire")
res = xe.evaluate(name, (293.15, 77.0))
n, L = res["numberdensity_A3"], res["L08_293K"]
sig = rj.sigma_tot("Al2O3")
print(f"L08(293 K) = {L:.4f} cm, L08(77 K) = {res['L08_77K']:.4f} cm")
print(f"A_fast_1E = {rj.a_fast(sig, n, L, '1E')[0]:.4f}, A_gamma_1MeV = {rj.mu_gamma('Al2O3', n)[0] * L:.4f}")
at = pd.read_parquet("data/derived/atlas.parquet")
row = at[(at.material_id == "mp-aaaaabrz") & (at["T"] > 100)].iloc[0]
print(f"atlas:  L08 = {row.L08:.4f} cm, A_fast_1E = {row.A_fast_1E:.4f}, A_gamma_1MeV = {row.A_gamma_1MeV:.4f}")
```

```
L08(293 K) = 5.9771 cm, L08(77 K) = 11.1258 cm
A_fast_1E = 2.0086, A_gamma_1MeV = 1.4778
atlas:  L08 = 5.9771 cm, A_fast_1E = 2.0086, A_gamma_1MeV = 1.4778
```

The DOS below 2 meV is dropped (`vdos_floor_meV` in `configs/vdos.yaml`); `build_vdos_material(dos, name,
floor_meV=...)` changes it for one crystal.

## 3. Change the temperature

For single crystals, pass any temperature to `evaluate`, one at a time. Its result keys are labelled `293K` for any
temperature above 100 K and `77K` otherwise:

```python
dos = load_dos("data/upstream/mp/phonon_dos/pheasy/aaaegyio.json.gz")
name = ncmat.build_vdos_material(dos, "sapphire")
for T in (200.0, 150.0):
    r = xe.evaluate(name, (T,))
    tag = "293K" if T > 100 else "77K"
    print(f"T = {T:.0f} K: L08 = {r['L08_' + tag]:.3f} cm")
```

```
T = 200 K: L08 = 8.027 cm
T = 150 K: L08 = 9.398 cm
```

For the whole atlas, `temperatures_K` in `configs/figures_of_merit.yaml` sets the two temperatures of steps 02, 04 and
08. **Limitation:** the atlas assembly supports exactly two temperatures, one above and one below 100 K, and writes
the column `T` as 293.15 and 77.0 whatever the configured values are. Change the two values only if you also accept
that label, and rerun steps 02, 04 and 08. For a series of temperatures, step 10 computes the commercially grown
shortlist at the temperatures of `temperature_series_K` in `configs/shortlist.yaml`.

## 4. Add a neutron or photon spectrum

The atlas has two neutron spectra (1/E over 1 eV-100 keV, Watt over 100 keV-10 MeV; `spectrum_1E_eV`,
`spectrum_watt_eV`, `watt_a_eV`, `watt_b_per_eV` in `configs/figures_of_merit.yaml`) and two photon energies
(`photon_energies_eV`). Changing the ranges or the Watt parameters and rerunning step 08 changes `A_fast_1E` and
`A_fast_watt`. Any other spectrum, for a single crystal, is a weighted average on the nuclear-data grid (1 eV-10 MeV):

```python
import numpy as np
t = rj.tables()
E = t["E"]
m = (E >= 1e3) & (E <= 1e4)                              # 1/E over 1-10 keV
w = np.zeros(m.sum()); de = np.diff(E[m]); s = 1.0 / E[m]   # trapezoid weights of S(E) = 1/E
w[:-1] += 0.5 * de * s[:-1]; w[1:] += 0.5 * de * s[1:]; w /= w.sum()
mu, info = xe.mu_offbragg(name, 293.15, 1.8)
n, L = info.numberdensity, 5.977
print(f"rejection: {-np.log(np.dot(w, np.exp(-n * rj.sigma_tot('Al2O3')[m] * L))):.3f}")
```

```
rejection: 1.949
```

Other photon energies need new elemental tables: set `photon_energies_eV` and rerun step 07 (OpenMC and the ENDF
library). The code keeps two photon energies (`mu_gamma` returns the first two) and the column names say 1 MeV and
6 MeV whatever the values are.

## 5. Change the thermal band or the transmission level

For a single crystal, compute μ(λ) on your band and solve for your transmission:

```python
lam = np.linspace(0.5, 4.0, 200)
mu, info = xe.mu_offbragg(name, 293.15, lam)
print(f"L0.9 over 0.5-4.0 A: {xe.l_target(mu, lam, target=0.90):.3f} cm")
```

```
L0.9 over 0.5-4.0 A: 2.794 cm
```

`l_target` weights the band with the 300 K Maxwellian; for a cold moderator use `l_target_band` with
`lambda_T(T_moderator)` on a curve from `curves` (0.3-4.5 Å). For the commercially grown crystals the curves are
stored in `data/intermediate/stage3_curves.parquet` (columns `curve_full_293K`, `curve_full_77K`, ...), so no
NCrystal run is needed.

For the whole atlas: `band_A`, `hot_band_A`, `moderator_T_K`/`lambda_T_A` and `T_target` in
`configs/figures_of_merit.yaml`, then rerun steps 02, 04 and 08 (and 05, 06 if the validation should follow).
`lambda_T_A` is stored as a number; keep it equal to h/√(2mₙk_B·T_moderator) when you change the moderator.

## 6. Run only one step, or part of one

```
python pipelines/08_atlas.py --list                 # the sub-steps
python pipelines/08_atlas.py --only fronts          # one sub-step
python pipelines/08_atlas.py --from controls        # from a sub-step to the end
python pipelines/10_shortlist.py --with openmc      # include the OpenMC sub-steps
```

A step reads its inputs from `data/` and overwrites its outputs there. To keep the released files, work on a copy:
use a second clone with its own unpacked data (`scripts/check_data.py` verifies the `data/` of the clone it is in), or
copy `data/` elsewhere and set `FA_DATA` to the copy. After changing a configuration value, rerun the steps that read
it ([Configuration](configuration.md) lists them) and every later step that uses their outputs (the flow in
[Pipeline overview](pipeline-overview.md)). Then `python -m pytest -q tests` reruns the controls on the new atlas.

## Add or correct a commercially grown crystal

The vendor table `data/curated/growable_crystals.csv` ([Curated tables](curated-tables.md)) decides which entries
are tier A. Add a row with the formula, the stated crystal system or phase (with the space-group symbol if the
vendor gives it), the vendor and the page; rerun step 01 (`--only tiers`) and then steps 03 and 08. The matching
rules are in `tier_a_matching` of `configs/funnel.yaml`. Changing a curated file changes its checksum, so
`check_data.py --class curated` reports it as mismatched; rebuild the manifest with `make manifest` before
distributing your data.
