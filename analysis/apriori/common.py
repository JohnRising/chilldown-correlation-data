"""Shared helpers for the a-priori analysis of the Darr et al. (2016) low-g LN2 chilldown data.

Geometry, data loading, nitrogen properties (CoolProp, cached), 304SS properties (NIST
cryogenic material properties database fits), and small numerical utilities.
"""
from __future__ import annotations

import functools
import os

import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from scipy.signal import savgol_filter
from scipy.special import erfc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "uf2016")
FIG = os.path.join(HERE, "figs")
OUT = os.path.join(HERE, "out")
os.makedirs(FIG, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- geometry
D = 11.68e-3          # inner diameter, m
DO = 12.70e-3         # outer diameter, m
LTUBE = 0.572         # m
Z = {"up": 0.149, "dn": 0.401}   # station distance from test-section inlet, m
RHO_SS = 7900.0       # kg/m3 (304SS)
# wall mass per unit INNER area: rho_w * (Do^2 - Di^2) / (4 Di)  [kg/m2]
MW_A = RHO_SS * (DO**2 - D**2) / (4 * D)
DT = 1.0 / 15.0
CASES = list(range(1, 11))
STATIONS = ("up", "dn")
TC_OF = {"up": ["TC1_K", "TC2_K", "TC3_K"], "dn": ["TC4_K", "TC5_K", "TC6_K"]}

# Paper Table 1 (mean over low-g period): Re_in, P1 kPa, G
TABLE1 = {1: (2838, 169.9, 31.3), 2: (6145, 242.0, 67.8), 3: (8801, 100.3, 116.2),
          4: (12133, 135.3, 143.8), 5: (13318, 126.0, 166.2), 6: (13881, 288.2, 130.4),
          7: (28704, 317.6, 262.3), 8: (35291, 365.5, 304.4), 9: (46245, 441.7, 375.3),
          10: (74177, 711.9, 503.1)}

FLUID = "Nitrogen"
TC_N2 = CP.PropsSI("Tcrit", FLUID)
PC_N2 = CP.PropsSI("pcrit", FLUID)


def load(case: int) -> pd.DataFrame:
    d = pd.read_csv(os.path.join(DATA, f"case{case}.csv"))
    return d


def station_P(d: pd.DataFrame, st: str) -> np.ndarray:
    """Local pressure [Pa] at a station.  P1 is at the upstream station; P2 at ~0.37 m.
    For the downstream station (0.401 m) we use P2 (closest)."""
    return (d["P1_kPa"].values if st == "up" else d["P2_kPa"].values) * 1e3


def smooth(y, win=9, order=2, deriv=0):
    win = min(win, len(y) - (1 - len(y) % 2))
    if win % 2 == 0:
        win -= 1
    return savgol_filter(y, win, order, deriv=deriv, delta=DT)


# ---------------------------------------------------------------- 304 SS (NIST)
# NIST Cryogenic Material Properties database, 304 Stainless (UNS S30400):
# log10(y) = sum a_i (log10 T)^i ; k: 1-300 K (2 %), cp: 4-300 K (5 %)
_SS_K = [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199]
_SS_CP = [22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212, 24.7593, -2.239153, 0.0]


def _nistpoly(c, T):
    x = np.log10(np.clip(np.asarray(T, float), 4.0, 300.0))
    return 10 ** sum(ci * x**i for i, ci in enumerate(c))


def ss_k(T):
    return _nistpoly(_SS_K, T)


def ss_cp(T):
    return _nistpoly(_SS_CP, T)


def ss_e(T):
    """thermal effusivity sqrt(k rho c) of 304SS."""
    return np.sqrt(ss_k(T) * RHO_SS * ss_cp(T))


# ---------------------------------------------------------------- N2 properties
@functools.lru_cache(maxsize=200000)
def _sat(Pr: float):
    P = Pr * 1e2  # Pa (Pr in hPa resolution)
    g = lambda k, q: CP.PropsSI(k, "P", P, "Q", q, FLUID)
    Ts = g("T", 0)
    return dict(
        Tsat=Ts, rho_l=g("D", 0), rho_v=g("D", 1), h_fg=g("H", 1) - g("H", 0),
        mu_l=g("V", 0), mu_v=g("V", 1), k_l=g("L", 0), k_v=g("L", 1),
        cp_l=g("C", 0), cp_v=g("C", 1), sigma=g("I", 0))


def sat(P):
    """saturated properties dict of arrays for pressure array P [Pa]."""
    P = np.atleast_1d(np.asarray(P, float))
    rows = [_sat(round(p / 1e2)) for p in P]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


@functools.lru_cache(maxsize=400000)
def _vap(Pr: float, Tr: float):
    P, T = Pr * 1e2, Tr * 0.1
    Ts = CP.PropsSI("T", "P", P, "Q", 1, FLUID)
    T = max(T, Ts + 0.05)
    g = lambda k: CP.PropsSI(k, "P", P, "T", T, FLUID)
    return dict(rho=g("D"), mu=g("V"), k=g("L"), cp=g("C"))


def vap(P, T):
    """superheated vapour properties at (P [Pa], T [K]) arrays."""
    P = np.atleast_1d(np.asarray(P, float))
    T = np.broadcast_to(np.asarray(T, float), P.shape)
    rows = [_vap(round(p / 1e2), round(t * 10)) for p, t in zip(P, T)]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------- superheat limits
def T_hn_lienhard(Tsat):
    """Lienhard (1976) homogeneous-nucleation (superheat-limit) correlation:
    (T_hn - Tsat)/Tc = 0.905 - Tr_sat + 0.095 Tr_sat^8  ->  T_hn = Tc(0.905 + 0.095 Tr^8)."""
    Tr = np.asarray(Tsat) / TC_N2
    return TC_N2 * (0.905 + 0.095 * Tr**8)


def T_spinodal_vdw(P):
    """liquid spinodal temperature of a van der Waals fluid at reduced pressure P/Pc."""
    from scipy.optimize import brentq
    Pr = np.atleast_1d(np.asarray(P, float)) / PC_N2
    out = []
    for p in Pr:
        # on the spinodal: Tr = (3v-1)^2/(4v^3), Pr = (3v-2)/v^3 ; liquid branch v in (2/3, 1)
        v = brentq(lambda v: (3 * v - 2) / v**3 - p, 0.6667 + 1e-9, 1.0)
        out.append(TC_N2 * (3 * v - 1) ** 2 / (4 * v**3))
    return np.array(out)


def T_spinodal_eos(P, Tlo=112.0, Thi=124.0):
    """liquid spinodal from the CoolProp (Span et al. 2000) Helmholtz EOS evaluated in the
    metastable region (imposed liquid phase): T at which the spinodal pressure equals P.
    NB: multiparameter EOS are not constrained in the metastable region -> indicative only."""
    AS = CP.AbstractState("HEOS", FLUID)
    AS.specify_phase(CP.iphase_liquid)

    def psp(T):
        rl = CP.PropsSI("D", "T", T, "Q", 0, FLUID)
        for rho in np.arange(rl, rl * 0.5, -0.1):
            AS.update(CP.DmassT_INPUTS, rho, T)
            if AS.first_partial_deriv(CP.iP, CP.iDmass, CP.iT) <= 0:
                return AS.p()
        return np.nan
    Ts = np.arange(Tlo, Thi + 1e-9, 0.25)
    ps = np.array([psp(t) for t in Ts])
    ok = np.isfinite(ps)
    out = np.interp(np.atleast_1d(P), ps[ok], Ts[ok])
    return out


def contact_wall_T(T_target, T_l, P):
    """wall temperature for which the semi-infinite contact temperature with liquid at T_l
    equals T_target.  Solves T_target = (e_w(Tw) Tw + e_l T_l)/(e_w + e_l)."""
    s = sat(P)
    e_l = np.sqrt(s["k_l"] * s["rho_l"] * s["cp_l"])
    Tw = np.array(T_target, float).copy()
    for _ in range(30):
        e_w = ss_e(Tw)
        Tw = T_target + e_l / e_w * (T_target - T_l)
    return Tw, e_l / ss_e(Tw)


def baumeister_B(Tw):
    """Baumeister & Simon (1973) transient-contact factor B = exp(3.06e6 beta) erfc(1751.5 sqrt(beta)),
    beta = 1/(k rho c)_w (SI)."""
    beta = 1.0 / (ss_k(Tw) * RHO_SS * ss_cp(Tw))
    return np.exp(3.06e6 * beta) * erfc(1751.5 * np.sqrt(beta))


# ---------------------------------------------------------------- statistics
def err_stats(pred, meas):
    pred, meas = np.asarray(pred, float), np.asarray(meas, float)
    ok = np.isfinite(pred) & np.isfinite(meas) & (pred > 0) & (meas > 0)
    r = np.log(pred[ok] / meas[ok])
    rat = pred[ok] / meas[ok]
    return dict(n=int(ok.sum()), bias=float(np.exp(r.mean()) - 1),
                rms_log=float(np.sqrt(np.mean(r**2))),
                sd_log=float(r.std()),
                w20=float(np.mean(np.abs(rat - 1) <= 0.2)),
                w30=float(np.mean(np.abs(rat - 1) <= 0.3)),
                mape=float(np.mean(np.abs(rat - 1))))


def md_table(df: pd.DataFrame, floatfmt=3) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("" if not np.isfinite(v) else (f"{v:.{floatfmt}g}" if abs(v) < 1e5 else f"{v:.3e}"))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)
