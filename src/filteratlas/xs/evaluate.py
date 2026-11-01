"""Cross sections, attenuation and figures of merit (definitions fixed before the computation).

mu in cm^-1: NCrystal numberdensity is atoms/A^3 and sigma is barn, so n[A^-3]*sigma[b] = mu[cm^-1].
"""
import numpy as np
import NCrystal as NC
from scipy.optimize import brentq

from filteratlas.common.config import cfg
from filteratlas.common.constants import BAND, HOT_BAND, LAMBDA_T, N_LAMBDA, PROBE_WL, T_TARGET

H = 6.62607015e-34
M_N = 1.67492749804e-27
K_B = 1.380649e-23


def lambda_T(T=300.0):
    """Characteristic wavelength of a Maxwellian of temperature T, h / sqrt(2 m_n k_B T).

    Parameters:
        T: temperature in K (default 300).
    Returns:
        wavelength in Å.
    Example:
        >>> from filteratlas.xs.evaluate import lambda_T
        >>> round(lambda_T(300.0), 4)
    """
    return H / np.sqrt(2 * M_N * K_B * T) * 1e10


def maxwellian_flux(lam, lam_t=LAMBDA_T):
    """Maxwellian flux per unit wavelength, phi(lambda) ~ lambda^-5 exp(-lambda_T^2 / lambda^2), unnormalised.

    Parameters:
        lam: wavelengths in Å (array).
        lam_t: characteristic wavelength in Å (default: `lambda_T_A` of `configs/figures_of_merit.yaml`, 300 K).
    Returns:
        array of weights, same shape as lam.
    Example:
        >>> from filteratlas.xs.evaluate import LAM, maxwellian_flux
        >>> maxwellian_flux(LAM).argmax()
    """
    return lam ** -5 * np.exp(-(lam_t / lam) ** 2)


LAM = np.linspace(*BAND, N_LAMBDA)
LAM_HOT = np.linspace(*HOT_BAND, N_LAMBDA)
_W = cfg("figures_of_merit")["wide_curve_A"]
LAM_WIDE = np.linspace(_W["min"], _W["max"], _W["n"])   # stored curves for band/moderator sensitivity


def curves(cfg, T):
    """Full-model and composition-only attenuation on the wide wavelength grid `LAM_WIDE` (0.3-4.5 Å).

    Used to recompute L0.8 for another band or moderator temperature without rerunning NCrystal (see l_target_band).

    Parameters:
        cfg: NCrystal configuration string or registered material name.
        T: crystal temperature in K.
    Returns:
        (mu_full, mu_proxy): two lists of float32 values in 1/cm on `LAM_WIDE`.
    Example:
        >>> from filteratlas.xs.evaluate import curves
        >>> mu_full, mu_proxy = curves("stdlib::Al2O3_sg167_Corundum.ncmat", 293.15)
    """
    mu, info = mu_offbragg(cfg, T, LAM_WIDE)
    return mu.astype("float32").tolist(), mu_proxy(info, LAM_WIDE).astype("float32").tolist()


def l_target_band(mu_wide, band, lam_t=LAMBDA_T, n=N_LAMBDA):
    """L0.8 for any band and moderator temperature, from a curve stored on `LAM_WIDE`.

    Parameters:
        mu_wide: attenuation in 1/cm on `LAM_WIDE` (from curves).
        band: (lowest, highest) wavelength in Å, inside 0.3-4.5 Å.
        lam_t: characteristic wavelength of the moderator in Å (lambda_T(T_moderator)).
        n: grid points in the band.
    Returns:
        thickness in cm at which the band transmission is T_target (0.80); inf if the attenuation is too small.
    Example:
        >>> from filteratlas.xs.evaluate import curves, l_target_band, lambda_T
        >>> l_target_band(curves("stdlib::Al2O3_sg167_Corundum.ncmat", 293.15)[0], (0.5, 4.5), lambda_T(40.0))
    """
    lam = np.linspace(*band, n)
    mu = np.interp(lam, LAM_WIDE, mu_wide)
    w = maxwellian_flux(lam, lam_t)
    f = lambda L: float(np.trapezoid(w * np.exp(-mu * L), lam) / np.trapezoid(w, lam)) - T_TARGET  # noqa: E731
    hi = 1.0
    while f(hi) > 0:
        hi *= 4
        if hi > 1e6:
            return float("inf")
    return float(brentq(f, 0.0, hi, xtol=1e-6, rtol=1e-8))


def mu_offbragg(cfg, T, wl):
    """Off-Bragg single-crystal attenuation coefficient: inelastic plus incoherent elastic scattering plus absorption.

    Coherent elastic (Bragg) scattering is left out, as for a large crystal oriented away from Bragg reflections.

    Parameters:
        cfg: NCrystal configuration string or registered material name (without the temperature).
        T: crystal temperature in K.
        wl: wavelengths in Å (scalar or array).
    Returns:
        (mu, info): mu in 1/cm at the wavelengths, and the NCrystal Info object.
    Example:
        >>> from filteratlas.xs.evaluate import mu_offbragg
        >>> mu, info = mu_offbragg("stdlib::Al2O3_sg167_Corundum.ncmat", 293.15, [1.0, 1.8])
    """
    c = f"{cfg};temp={T}K"
    info = NC.createInfo(c)
    sc = NC.createScatter(f"{c};comp=inelas,incoh_elas")
    ab = NC.createAbsorption(c)
    wl = np.atleast_1d(np.asarray(wl, float))
    return info.numberdensity * (sc.xsect(wl=wl) + ab.xsect(wl=wl)), info


def mu_proxy(info, wl):
    """Composition-only attenuation: n * sum_i x_i (sigma_abs,i(lambda) + sigma_inc,i), free elements, no phonons.

    Parameters:
        info: NCrystal Info object (from mu_offbragg), which gives the composition and the atom density.
        wl: wavelengths in Å (scalar or array).
    Returns:
        mu in 1/cm at the wavelengths.
    Example:
        >>> from filteratlas.xs.evaluate import mu_offbragg, mu_proxy
        >>> mu_proxy(mu_offbragg("stdlib::Al2O3_sg167_Corundum.ncmat", 293.15, 1.8)[1], [1.0, 1.8])
    """
    wl = np.atleast_1d(np.asarray(wl, float))
    sabs = sum(x * ad.captureXS() for x, ad in info.composition)
    sinc = sum(x * ad.incoherentXS() for x, ad in info.composition)
    return info.numberdensity * (sabs * wl / 1.798 + sinc)


def t_band(mu, L, lam=LAM):
    """Transmission of a slab of thickness L averaged over the thermal band with Maxwellian weights.

    Parameters:
        mu: attenuation in 1/cm on the wavelength grid lam.
        L: thickness in cm.
        lam: wavelength grid in Å (default `LAM`, the thermal band 0.5-3 Å).
    Returns:
        transmission between 0 and 1.
    Example:
        >>> from filteratlas.xs.evaluate import LAM, mu_offbragg, t_band
        >>> t_band(mu_offbragg("stdlib::Al2O3_sg167_Corundum.ncmat", 293.15, LAM)[0], 5.0)
    """
    w = maxwellian_flux(lam)
    return float(np.trapezoid(w * np.exp(-mu * L), lam) / np.trapezoid(w, lam))


def l_target(mu, lam=LAM, target=T_TARGET):
    """Thickness at which the band transmission equals the target (L0.8 for the default target 0.80).

    Parameters:
        mu: attenuation in 1/cm on the grid lam.
        lam: wavelength grid in Å (default `LAM`).
        target: transmission (default `T_target`, 0.80).
    Returns:
        thickness in cm (Brent root); inf if the attenuation is zero or too small.
    Example:
        >>> from filteratlas.xs.evaluate import LAM, l_target, mu_offbragg
        >>> l_target(mu_offbragg("stdlib::Si_sg227.ncmat", 293.15, LAM)[0])
    """
    if np.all(mu <= 0):
        return float("inf")
    f = lambda L: t_band(mu, L, lam) - target  # noqa: E731
    hi = 1.0
    while f(hi) > 0:
        hi *= 4
        if hi > 1e6:
            return float("inf")
    return float(brentq(f, 0.0, hi, xtol=1e-6, rtol=1e-8))


def evaluate(cfg, temps):
    """All per-material thermal quantities of the atlas, for each temperature: full model and composition-only proxy.

    Parameters:
        cfg: NCrystal configuration string or registered material name.
        temps: temperatures in K; a temperature above 100 K is labelled 293K, otherwise 77K.
    Returns:
        dict with mu_1.0_<tag>, mu_1.8_<tag> (1/cm), proxy_1.0_<tag>, proxy_1.8_<tag>, L08_<tag>, L08_proxy_<tag>,
        L08_hot_<tag>, L08_hot_proxy_<tag> (cm) and numberdensity_A3 (atoms/Å³).
    Example:
        >>> from filteratlas.xs.evaluate import evaluate
        >>> evaluate("stdlib::Al2O3_sg167_Corundum.ncmat", (293.15, 77.0))["L08_77K"]
    """
    out = {}
    for T in temps:
        mu, info = mu_offbragg(cfg, T, LAM)
        mu_hot, _ = mu_offbragg(cfg, T, LAM_HOT)
        mu_probe, _ = mu_offbragg(cfg, T, PROBE_WL)
        mp_ = mu_proxy(info, LAM)
        tag = "293K" if T > 100 else "77K"
        out.update({
            f"mu_1.0_{tag}": float(mu_probe[0]), f"mu_1.8_{tag}": float(mu_probe[1]),
            f"proxy_1.0_{tag}": float(mu_proxy(info, PROBE_WL)[0]),
            f"proxy_1.8_{tag}": float(mu_proxy(info, PROBE_WL)[1]),
            f"L08_{tag}": l_target(mu), f"L08_proxy_{tag}": l_target(mp_),
            f"L08_hot_{tag}": l_target(mu_hot, LAM_HOT),
            f"L08_hot_proxy_{tag}": l_target(mu_proxy(info, LAM_HOT), LAM_HOT),
        })
    out["numberdensity_A3"] = float(info.numberdensity)
    return out
