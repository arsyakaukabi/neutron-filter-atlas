"""Build NCrystal materials from MP structures with per-element VDOS."""
import numpy as np
import NCrystal as NC
from ase import Atoms
from pymatgen.core import Structure

from filteratlas.common.config import cfg
from filteratlas.common.constants import THZ_TO_EV, UNSTABLE_FREQ_THZ, UNSTABLE_WEIGHT_MAX, VDOS_NPTS

H_PLANCK = 6.62607015e-34
K_B = 1.380649e-23
NCRYSTAL_EMIN_EV = cfg("vdos")["ncrystal_emin_eV"]   # NCrystal's minimum VDOS energy
VDOS_FLOOR_MEV = cfg("vdos")["vdos_floor_meV"]       # 2 meV floor: validated against force-constant tetrahedron DOS (stage1r_fc_validation.json)
N_A = 6.02214076e23


def to_ase(structure):
    """ASE Atoms from a pymatgen Structure (or its dict), keeping only species, positions and cell.

    Parameters:
        structure: pymatgen Structure or Structure.as_dict().
    Returns:
        ase.Atoms with periodic boundary conditions.
    Example:
        >>> from filteratlas.xs.ncmat import to_ase
        >>> to_ase(dos["structure"]).get_chemical_formula()   # dos from filteratlas.ingest.mp.load_dos
    """
    s = Structure.from_dict(structure) if isinstance(structure, dict) else structure
    return Atoms(symbols=[str(site.specie) for site in s], scaled_positions=s.frac_coords,
                 cell=s.lattice.matrix, pbc=True)


def composer(structure):
    """NCrystal NCMATComposer for a structure, and the map from its atom labels to elements.

    Parameters:
        structure: pymatgen Structure or its dict.
    Returns:
        (composer, labels): the NCMATComposer and a dict {label: element symbol}.
    Example:
        >>> from filteratlas.xs.ncmat import composer
        >>> c, labels = composer(dos["structure"])
    """
    c = NC.NCMATComposer.from_ase(to_ase(structure), quiet=True)
    labels = {lab: comp[0][1] for lab, comp in c.to_dict()["compositions"].items()}
    return c, labels


def stability(dos):
    """DOS rule for dynamical stability: the fraction of the DOS weight below `unstable_freq_THz` (−0.5 THz).

    Parameters:
        dos: phonon DOS dict with "frequencies" (THz) and "densities".
    Returns:
        (weight, stable): the unstable fraction and whether it is at most `unstable_weight_max` (0.5 %).
    Example:
        >>> from filteratlas.xs.ncmat import stability
        >>> stability(dos)
    """
    f = np.asarray(dos["frequencies"], float)
    d = np.clip(np.asarray(dos["densities"], float), 0, None)
    tot = np.trapezoid(d, f)
    neg = np.trapezoid(np.where(f < UNSTABLE_FREQ_THZ, d, 0.0), f)
    w = float(neg / tot) if tot > 0 else 1.0
    return w, w <= UNSTABLE_WEIGHT_MAX


def element_vdos(dos, floor_meV=None):
    """Per-element vibrational DOS from the site-projected phonon DOS, on an energy grid in eV.

    Bins below the floor are dropped and NCrystal extrapolates g(E) ~ E^2 below the first kept point (the DOS near zero
    frequency carries a smearing tail of the acoustic modes).

    Parameters:
        dos: phonon DOS dict with "frequencies" (THz), "projected_densities" (one row per site) and "structure".
        floor_meV: lowest energy kept in meV (None = `vdos_floor_meV` of `configs/vdos.yaml`, 2 meV).
    Returns:
        dict {element: (energy grid in eV, density)}.
    Raises:
        ValueError if the projected DOS does not match the structure or an element has no positive-frequency DOS.
    Example:
        >>> from filteratlas.xs.ncmat import element_vdos
        >>> sorted(element_vdos(dos))
    """
    floor_meV = VDOS_FLOOR_MEV if floor_meV is None else floor_meV
    s = Structure.from_dict(dos["structure"])
    f = np.asarray(dos["frequencies"], float)
    proj = np.asarray(dos["projected_densities"], float)
    if proj.shape != (len(s), len(f)):
        raise ValueError(f"projected DOS shape {proj.shape} != sites x freq ({len(s)}, {len(f)})")
    out = {}
    for el in sorted({str(x.specie) for x in s}):
        dens = np.clip(proj[[i for i, x in enumerate(s) if str(x.specie) == el]].sum(axis=0), 0, None)
        m = f * THZ_TO_EV >= max(NCRYSTAL_EMIN_EV, floor_meV * 1e-3)   # NCrystal minimum and the VDOS floor
        if m.sum() < 5 or dens[m].max() <= 0:
            raise ValueError(f"no positive-frequency DOS for {el}")
        fp, dp = f[m], dens[m]
        grid = np.linspace(fp[0], fp[-1], VDOS_NPTS)
        out[el] = (grid * THZ_TO_EV, np.interp(grid, fp, dp))
    return out


def build_vdos_material(dos, name, floor_meV=None, vdos=None):
    """Register an in-memory NCrystal material from a structure and per-element vibrational DOS.

    Parameters:
        dos: phonon DOS dict with "structure" (and "frequencies", "projected_densities" unless vdos is given).
        name: name for the material; ".ncmat" is appended.
        floor_meV: VDOS floor in meV (None = 2 meV).
        vdos: optional precomputed {element: (energy grid in eV, density)}; otherwise taken from the DOS.
    Returns:
        the registered material name, usable wherever NCrystal takes a configuration string.
    Raises:
        RuntimeError if NCrystal would fall back to a Debye model for any atom.
    Example:
        >>> from filteratlas.xs.ncmat import build_vdos_material
        >>> name = build_vdos_material(dos, "my_crystal")   # then filteratlas.xs.evaluate.evaluate(name, (293.15, 77.0))
    """
    c, labels = composer(dos["structure"])
    vd = vdos if vdos is not None else element_vdos(dos, floor_meV)
    for lab, el in labels.items():
        eg, dd = vd[el]
        c.set_dyninfo_vdos(lab, eg, dd)
    text = c.create_ncmat()
    if "vdosdebye" in text or "debye_temp" in text:  # every label must carry its own VDOS
        raise RuntimeError(f"{name}: fallback Debye dynamics present in generated NCMAT")
    fname = f"{name}.ncmat"
    c.register_as(fname)
    return fname


def debye_temperature(K_gpa, G_gpa, density_gcc, nsites, volume_A3):
    """Debye temperature from the Voigt-Reuss-Hill moduli (Anderson 1963).

    Parameters:
        K_gpa: bulk modulus in GPa.
        G_gpa: shear modulus in GPa.
        density_gcc: mass density in g/cm³.
        nsites: atoms in the cell.
        volume_A3: cell volume in Å³.
    Returns:
        Debye temperature in K, or nan if a modulus is missing or not positive.
    Example:
        >>> from filteratlas.xs.ncmat import debye_temperature
        >>> debye_temperature(165.0, 66.0, 2.33, 8, 160.2)   # silicon
    """
    if not (K_gpa and G_gpa) or K_gpa <= 0 or G_gpa <= 0:
        return float("nan")
    rho = density_gcc * 1e3
    vs = np.sqrt(G_gpa * 1e9 / rho)
    vl = np.sqrt((K_gpa + 4 * G_gpa / 3) * 1e9 / rho)
    vm = ((2 / vs ** 3 + 1 / vl ** 3) / 3) ** (-1 / 3)
    n_atoms_per_m3 = nsites / (volume_A3 * 1e-30)
    return float(H_PLANCK / K_B * (3 * n_atoms_per_m3 / (4 * np.pi)) ** (1 / 3) * vm)


def build_debye_material(structure, theta_d, name):
    """Register an NCrystal material in which every atom has the same Debye temperature.

    Parameters:
        structure: pymatgen Structure or its dict.
        theta_d: Debye temperature in K.
        name: name for the material; ".ncmat" is appended.
    Returns:
        the registered material name.
    Example:
        >>> from filteratlas.xs.ncmat import build_debye_material
        >>> name = build_debye_material(dos["structure"], 1000.0, "my_debye")
    """
    c, labels = composer(structure)
    for lab in labels:
        c.set_dyninfo_debyetemp(lab, float(theta_d))
    fname = f"{name}.ncmat"
    c.register_as(fname)
    return fname
