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
    return H / np.sqrt(2 * M_N * K_B * T) * 1e10


def maxwellian_flux(lam, lam_t=LAMBDA_T):
    """Maxwellian flux per unit wavelength, phi(lambda) ~ lambda^-5 exp(-lambda_T^2/lambda^2)."""
    return lam ** -5 * np.exp(-(lam_t / lam) ** 2)


LAM = np.linspace(*BAND, N_LAMBDA)
LAM_HOT = np.linspace(*HOT_BAND, N_LAMBDA)
_W = cfg("figures_of_merit")["wide_curve_A"]
LAM_WIDE = np.linspace(_W["min"], _W["max"], _W["n"])   # stored curves for band/moderator sensitivity


def curves(cfg, T):
    """mu_th and mu_proxy on LAM_WIDE [cm^-1] (float32 lists, for storage)."""
    mu, info = mu_offbragg(cfg, T, LAM_WIDE)
    return mu.astype("float32").tolist(), mu_proxy(info, LAM_WIDE).astype("float32").tolist()


def l_target_band(mu_wide, band, lam_t=LAMBDA_T, n=N_LAMBDA):
    """L0.8 for an arbitrary band / moderator temperature from a stored LAM_WIDE curve."""
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
    """Off-Bragg single-crystal baseline mu(lambda) [cm^-1]: inelastic + incoherent elastic + absorption."""
    c = f"{cfg};temp={T}K"
    info = NC.createInfo(c)
    sc = NC.createScatter(f"{c};comp=inelas,incoh_elas")
    ab = NC.createAbsorption(c)
    wl = np.atleast_1d(np.asarray(wl, float))
    return info.numberdensity * (sc.xsect(wl=wl) + ab.xsect(wl=wl)), info


def mu_proxy(info, wl):
    """Composition-only (elemental) proxy n*sum x_i(sigma_abs,i(lambda) + sigma_inc,i) [cm^-1]."""
    wl = np.atleast_1d(np.asarray(wl, float))
    sabs = sum(x * ad.captureXS() for x, ad in info.composition)
    sinc = sum(x * ad.incoherentXS() for x, ad in info.composition)
    return info.numberdensity * (sabs * wl / 1.798 + sinc)


def t_band(mu, L, lam=LAM):
    w = maxwellian_flux(lam)
    return float(np.trapezoid(w * np.exp(-mu * L), lam) / np.trapezoid(w, lam))


def l_target(mu, lam=LAM, target=T_TARGET):
    """Thickness [cm] at which band transmission equals target (Brent); inf if mu ~ 0."""
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
    """All per-material quantities (mu at the probe wavelengths, L0.8 in the thermal and hot bands): per temperature,
    full model and proxy."""
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
