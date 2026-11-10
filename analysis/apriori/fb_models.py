"""Film-boiling dataset and candidate closures (all gravity-free, low-g)."""
from __future__ import annotations

import functools

import numpy as np
import pandas as pd

from common import CASES, D, STATIONS, Z, load, sat, vap
from prep import film_mask, station_data

TAU_G = 2.5   # s, global flow-meter delay (a02_lag.py pooled optimum; paper: ~3 s for case 10)
X_IN = 0.0    # inlet equilibrium quality assumed (saturated liquid; inlet state not measured except case 7)


def g_testsection(case: int, tau: float = TAU_G) -> np.ndarray:
    d = load(case)
    t = d.t_s.values
    return np.interp(t + tau, t, d.G_kgm2s.values)


def x_eq(q_up, q_dn, G, hfg, st):
    """quasi-steady equilibrium quality from an energy balance inlet -> station using the measured
    flux at the two stations: q = q_up on [0, z_up], linear between z_up and z_dn."""
    zu, zd = Z["up"], Z["dn"]
    I = q_up * zu if st == "up" else q_up * zu + 0.5 * (q_up + q_dn) * (zd - zu)
    return X_IN + 4.0 * I / (np.clip(G, 5, None) * D * hfg)


@functools.lru_cache(maxsize=None)
def dataset(tau: float = TAU_G, all_samples: bool = False) -> pd.DataFrame:
    """film-boiling samples with properties.  all_samples=True returns the whole valid record
    (for the 0-D surrogate)."""
    rows = []
    for c in CASES:
        Gl = g_testsection(c, tau)
        Graw = load(c).G_kgm2s.values
        su, sd = station_data(c, "up"), station_data(c, "dn")
        n = min(len(su), len(sd))
        for st in STATIONS:
            s = station_data(c, st)
            m = np.ones(len(s), bool) if all_samples else film_mask(c, st)
            m[n:] = False
            df = pd.DataFrame({"case": c, "st": st, "t": s.t[m], "Ti": s.Ti[m], "Tsat": s.Tsat[m],
                               "P": s.P[m], "dT": s.dT[m], "h_ec": s.h_ec[m], "h_raw": s.h[m],
                               "q_ec": s.q_ec[m], "G": np.clip(Gl[:len(s)][m], 5, None),
                               "G_raw": np.clip(Graw[:len(s)][m], 5, None), "z": Z[st],
                               "q_up": su.q_ec.values[:len(s)][m], "q_dn": sd.q_ec.values[:len(s)][m]})
            rows.append(df)
    df = pd.concat(rows, ignore_index=True)
    sp = sat(df.P.values)
    for k, v in sp.items():
        df[k] = v
    Tf = 0.5 * (df.Ti.values + df.Tsat.values)
    vf = vap(df.P.values, Tf)
    for k, v in vf.items():
        df[f"{k}_vf"] = v
    vw = vap(df.P.values, df.Ti.values)
    df["Pr_vw"] = vw["cp"] * vw["mu"] / vw["k"]
    df["x"] = x_eq(df.q_up.values, df.q_dn.values, df.G.values, df.h_fg.values, "up")
    dn = df.st.values == "dn"
    df.loc[dn, "x"] = x_eq(df.q_up.values[dn], df.q_dn.values[dn], df.G.values[dn], df.h_fg.values[dn], "dn")
    df["x_raw"] = x_eq(df.q_up.values, df.q_dn.values, df.G_raw.values, df.h_fg.values, "up")
    df.loc[dn, "x_raw"] = x_eq(df.q_up.values[dn], df.q_dn.values[dn], df.G_raw.values[dn], df.h_fg.values[dn], "dn")
    return df


# ------------------------------------------------------------------ closures
# every closure: f(df, G) -> h [W/m2K]; G column passed explicitly (lag-corrected or raw)

def hfg_mod(df, beta=0.5):
    return df.h_fg + beta * df.cp_v * df.dT


def bromley_unit(df, G, props="film", beta=0.5):
    """Bromley, LeRoy & Robbers (1953) forced-convection film boiling without the C:
    sqrt(k_v rho_v h_fg' u_l / (D dT)), u_l = G/rho_l, h_fg' = h_fg + beta cp_v dT."""
    k, rho = (df.k_vf, df.rho_vf) if props == "film" else (df.k_v, df.rho_v)
    u = G / df.rho_l
    return np.sqrt(k * rho * hfg_mod(df, beta) * u / (D * df.dT))


def darr2016(df, G, x=None):
    """Darr et al. 2016 npj eqs 1-4 with low-g constants (c7 = 0 -> no buoyancy term).
    sat. properties at local Tsat; L = distance from the tube inlet (no quench front upstream)."""
    c1, c2, c3, c4, c5, c6 = 0.150, -0.76, 0.52, 5e-6, 0.10, 2.0
    x = df.x if x is None else x
    L = df.z
    Re_v = G * D / df.mu_v
    Pr_v = df.cp_v * df.mu_v / df.k_v
    htc = c1 * (1 + (D / (L + 0.1 * D)) ** 0.7) * (1 - np.clip(x, 0, 0.99)) ** c2 * Re_v**c3 * Pr_v**0.4 * df.k_v / D
    We = G**2 * D / (df.rho_l * df.sigma)
    Ja = df.cp_l * df.dT / df.h_fg
    hdw = c4 * df.k_l / D * np.exp(-c5 * L / D) * We**c6 / Ja
    return htc + hdw


def dittus_boelter_v(df, G, props="sat"):
    if props == "sat":
        mu, k, cp = df.mu_v, df.k_v, df.cp_v
    else:
        mu, k, cp = df.mu_vf, df.k_vf, df.cp_vf
    Re = G * D / mu
    return 0.023 * Re**0.8 * (cp * mu / k) ** 0.4 * k / D


def miropolskii(df, G, x):
    x = np.clip(x, 0, 1)
    Re = G * D / df.mu_v * (x + df.rho_v / df.rho_l * (1 - x))
    Y = 1 - 0.1 * (df.rho_l / df.rho_v - 1) ** 0.4 * (1 - x) ** 0.4
    return 0.023 * Re**0.8 * df.Pr_vw**0.8 * Y * df.k_v / D


def groeneveld(df, G, x):
    """Groeneveld (1973) tube correlation (a=1.09e-3, b=0.989, c=1.41, d=-1.15)."""
    x = np.clip(x, 0, 1)
    Re = G * D / df.mu_v * (x + df.rho_v / df.rho_l * (1 - x))
    Y = 1 - 0.1 * (df.rho_l / df.rho_v - 1) ** 0.4 * (1 - x) ** 0.4
    return 1.09e-3 * Re**0.989 * df.Pr_vw**1.41 * Y**-1.15 * df.k_v / D


def couette(df, G, ell, beta=0.5):
    """laminar Couette vapour film between the wall (u=0) and a liquid core moving at u_l,
    growing from a leading edge: vapour flow per unit perimeter = rho_v u_l delta/2,
    d/dz(rho_v u_l delta/2) = k_v dT /(delta h_fg')  ->  delta^2 = 4 k_v dT ell/(rho_v u_l h_fg'),
    h = k_v/delta = sqrt(rho_v u_l h_fg' k_v /(4 dT ell)).  Film-T vapour properties."""
    u = G / df.rho_l
    return np.sqrt(df.rho_vf * u * hfg_mod(df, beta) * df.k_vf / (4 * df.dT * ell))


def lambda_kh(df, G, urel_factor=1.0):
    """most-dangerous Kelvin-Helmholtz wavelength for a thin vapour layer sheared by the liquid
    (gravity-free): lambda = 2 pi sigma (rho_l + rho_v)/(rho_l rho_v u_rel^2) ~ 2 pi sigma/(rho_v u_rel^2),
    with u_rel = urel_factor * u_l."""
    u = urel_factor * G / df.rho_l
    return 2 * np.pi * df.sigma * (df.rho_l + df.rho_vf) / (df.rho_l * df.rho_vf * u**2)


def penetration_unit(df, G, ell=D):
    """sensible-heat (vapour-superheat) limited film / surface-renewal model:
    vapour elements in contact with the wall over a length ell while convected at u_l absorb heat by
    transient conduction: h = (2/sqrt(pi)) k_v sqrt(u_l/(alpha_v ell)) (gravity-free; no dT
    dependence except via properties).  Film-T properties.  Returned without a constant
    (C=1 corresponds to the 2/sqrt(pi) penetration value)."""
    u = G / df.rho_l
    alpha = df.k_vf / (df.rho_vf * df.cp_vf)
    return 2 / np.sqrt(np.pi) * df.k_vf * np.sqrt(u / (alpha * ell))
