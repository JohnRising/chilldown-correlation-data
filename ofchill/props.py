"""Closed-form property fits (generated from CoolProp) usable both in numpy and
as OpenFLUME expression strings (via ofchill.expr).

Saturation properties f(P):  y = poly(s),  s = (ln P - m)/h in [-1, 1]
    Tsat: y = 1/Tsat ; all others: y = ln(f)
Superheated-vapour properties f(T, P) on T in [Tsat(P), Tmax]:
    y = sum_ij c_ij xa^i s^j,  xa = normalised ln(T / Tsat(P))
    k, mu, cp: y = ln(f) ;  rho: y = ln(P / (rho R T))^-1 = ln(rho R T / P)
P is clipped to [Pmin, Pmax] and T to [Tsat(P), Tmax] inside the fit
(polynomials must not extrapolate).

Fits are generated on first use and cached in ofchill/_fits/<fluid>.json.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path

import numpy as np

from . import expr as E

FIT_DIR = Path(__file__).resolve().parent / "_fits"
FIT_VERSION = 4

# Default ranges (P in Pa, Tmax in K).
FLUID_RANGES = {
    "Nitrogen": dict(Pmin=80e3, Pmax=1500e3, Tmax=320.0),
    "Oxygen": dict(Pmin=50e3, Pmax=2000e3, Tmax=320.0),
    "Hydrogen": dict(Pmin=50e3, Pmax=1100e3, Tmax=320.0),
}
SAT_NAMES = ["Tsat", "rho_l", "rho_v", "h_fg", "sigma", "mu_l", "mu_v", "k_l", "k_v", "cp_l", "cp_v"]
VAP_NAMES = ["k", "mu", "cp", "rho"]
SAT_TARGET = {"Tsat": 1e-5, "rho_l": 1e-5, "rho_v": 2e-5}  # tighter: used for phase tests
SAT_TARGET_DEFAULT = 1e-4
VAP_TARGET = 2e-3


# ----------------------------------------------------------------------------
# CoolProp reference values
# ----------------------------------------------------------------------------
def _cp():
    import CoolProp.CoolProp as CP

    return CP


def coolprop_sat(fluid: str, P: np.ndarray) -> dict[str, np.ndarray]:
    CP = _cp()
    q = {
        "Tsat": (0, "T"), "rho_l": (0, "D"), "rho_v": (1, "D"), "mu_l": (0, "V"), "mu_v": (1, "V"),
        "k_l": (0, "L"), "k_v": (1, "L"), "cp_l": (0, "C"), "cp_v": (1, "C"), "sigma": (0, "I"),
    }
    out = {k: np.array([CP.PropsSI(p, "P", float(x), "Q", Q, fluid) for x in P]) for k, (Q, p) in q.items()}
    out["h_fg"] = np.array(
        [CP.PropsSI("H", "P", float(x), "Q", 1, fluid) - CP.PropsSI("H", "P", float(x), "Q", 0, fluid) for x in P]
    )
    return out


def coolprop_vap(fluid: str, T: np.ndarray, P: np.ndarray) -> dict[str, np.ndarray]:
    """Vapour properties; at T == Tsat(P) uses the saturated-vapour state."""
    CP = _cp()
    out = {k: np.empty(len(T)) for k in VAP_NAMES}
    code = dict(k="L", mu="V", cp="C", rho="D")
    for i, (t, p) in enumerate(zip(T, P)):
        ts = CP.PropsSI("T", "P", float(p), "Q", 1, fluid)
        for k, c in code.items():
            if t <= ts * (1 + 1e-9):
                out[k][i] = CP.PropsSI(c, "P", float(p), "Q", 1, fluid)
            else:
                out[k][i] = CP.PropsSI(c, "P", float(p), "T", float(t), fluid)
    return out


def gas_constant(fluid: str) -> float:
    CP = _cp()
    return CP.PropsSI("gas_constant", fluid) / CP.PropsSI("molar_mass", fluid)


# ----------------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------------
def _fit_sat(fluid, Pmin, Pmax):
    P = np.geomspace(Pmin, Pmax, 500)
    ref = coolprop_sat(fluid, P)
    m = 0.5 * (np.log(Pmax) + np.log(Pmin))
    h = 0.5 * (np.log(Pmax) - np.log(Pmin))
    s = (np.log(P) - m) / h
    fits = {}
    for k in SAT_NAMES:
        v = ref[k]
        y = 1 / v if k == "Tsat" else np.log(v)
        tgt = SAT_TARGET.get(k, SAT_TARGET_DEFAULT)
        for deg in range(2, 15):
            c = np.polynomial.polynomial.polyfit(s, y, deg)
            yf = np.polynomial.polynomial.polyval(s, c)
            vf = 1 / yf if k == "Tsat" else np.exp(yf)
            err = float(np.max(np.abs(vf / v - 1)))
            if err < tgt:
                break
        fits[k] = dict(coefs=c.tolist(), deg=deg, fit_err=err)
    return dict(m=m, h=h, fits=fits)


# Vapour-fit segmentation in a = ln(T/Tsat(P)) (fractions of a_max).  H2 needs
# segments (near-critical sat-vapour line + rotational cp hump).
VAP_SEGMENTS = {"Hydrogen": [0.0, 0.03, 0.12, 0.4, 1.0]}
VAP_DEGREES = [(4, 2), (6, 3), (8, 4), (8, 6), (10, 5), (10, 6), (12, 6), (12, 7), (14, 7), (14, 8), (12, 10), (14, 12)]


def _fit_vap(fluid, Pmin, Pmax, Tmax, sat):
    CP = _cp()
    Ps = np.geomspace(Pmin, Pmax, 60)
    T, P = [], []
    for p in Ps:
        ts = CP.PropsSI("T", "P", float(p), "Q", 1, fluid)
        a = np.linspace(0, 1, 120) ** 2 * np.log(Tmax / ts)
        T += list(ts * np.exp(a))
        P += [p] * len(a)
    T, P = np.array(T), np.array(P)
    ref = coolprop_vap(fluid, T, P)
    R = gas_constant(fluid)
    s = (np.log(P) - sat["m"]) / sat["h"]
    Ts = np.array([CP.PropsSI("T", "P", float(p), "Q", 1, fluid) for p in P])
    a = np.log(T / Ts)
    a_max = float(np.log(Tmax / CP.PropsSI("T", "P", Pmin, "Q", 1, fluid)))
    bounds = [b * a_max for b in VAP_SEGMENTS.get(fluid, [0.0, 1.0])]
    bounds[-1] = a_max * (1 + 1e-9)
    fits = {}
    for k in VAP_NAMES:
        v = ref[k]
        y = np.log(v * R * T / P) if k == "rho" else np.log(v)
        segs = []
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            # fit on a slightly widened window for continuity at the joins
            w = 0.1 * (hi - lo)
            msk = (a >= lo - w) & (a <= hi + w)
            xa = (2 * a[msk] - (lo + hi)) / (hi - lo)
            for dx, ds in VAP_DEGREES:
                A = np.stack([xa**i * s[msk] ** j for i in range(dx + 1) for j in range(ds + 1)], 1)
                c, *_ = np.linalg.lstsq(A, y[msk], rcond=None)
                err = float(np.max(np.abs(np.exp(A @ c - y[msk]) - 1)))
                if err < VAP_TARGET:
                    break
            segs.append(dict(lo=lo, hi=hi, coefs=np.asarray(c).reshape(dx + 1, ds + 1).tolist(), dx=dx, ds=ds,
                             fit_err=err))
        fits[k] = dict(segments=segs, fit_err=max(sg["fit_err"] for sg in segs))
    return dict(a_max=a_max, R=R, fits=fits)


@functools.lru_cache(maxsize=None)
def fluid_props(fluid: str = "Nitrogen") -> "FluidProps":
    rng = FLUID_RANGES[fluid]
    FIT_DIR.mkdir(exist_ok=True)
    f = FIT_DIR / f"{fluid}.json"
    data = None
    if f.exists():
        data = json.loads(f.read_text())
        if data.get("version") != FIT_VERSION or any(data.get(k) != v for k, v in rng.items()):
            data = None
    if data is None:
        sat = _fit_sat(fluid, rng["Pmin"], rng["Pmax"])
        vap = _fit_vap(fluid, rng["Pmin"], rng["Pmax"], rng["Tmax"], sat)
        data = dict(version=FIT_VERSION, fluid=fluid, **rng, sat=sat, vap=vap)
        f.write_text(json.dumps(data, indent=1))
    return FluidProps(data)


class FluidProps:
    """Dual-mode property fits.  All methods accept numpy arrays or Exprs.

    >>> n2 = fluid_props("Nitrogen")
    >>> n2.Tsat(3e5); n2.sat("h_fg", P); n2.vap("k", T, P)
    """

    def __init__(self, data: dict):
        self.data = data
        self.fluid = data["fluid"]
        self.Pmin, self.Pmax, self.Tmax = data["Pmin"], data["Pmax"], data["Tmax"]
        self.R = data["vap"]["R"]

    # -- helpers --------------------------------------------------------------
    def _s(self, P):
        Pc = E.clip(P, self.Pmin, self.Pmax)
        return (E.log(Pc) - self.data["sat"]["m"]) / self.data["sat"]["h"]

    def sat(self, name: str, P):
        f = self.data["sat"]["fits"][name]
        y = E.poly(f["coefs"], self._s(P))
        return 1.0 / y if name == "Tsat" else E.exp(y)

    def Tsat(self, P):
        return self.sat("Tsat", P)

    def vap(self, name: str, T, P):
        """Superheated/saturated vapour property at (T, P); T clipped to
        [Tsat(P), Tmax], P to [Pmin, Pmax].  Piecewise in a = ln(T/Tsat)
        (nested ternary; only the selected segment is evaluated)."""
        f = self.data["vap"]["fits"][name]
        Pc = E.clip(P, self.Pmin, self.Pmax)
        Ts = self.Tsat(Pc)
        Tc = E.clip(T, Ts, self.Tmax)
        a = E.log(Tc / Ts)
        s = self._s(Pc)

        def seg_val(sg):
            xa = (2 * a - (sg["lo"] + sg["hi"])) / (sg["hi"] - sg["lo"])
            return E.poly([E.poly(row, s) for row in sg["coefs"]], xa)

        segs = f["segments"]
        y = seg_val(segs[-1])
        for sg in reversed(segs[:-1]):
            y = E.where(a < sg["hi"], seg_val(sg), y)
        if name == "rho":
            return Pc / (self.R * Tc) * E.exp(y)
        return E.exp(y)

    def __getattr__(self, name):
        # n2.rho_l(P), n2.h_fg(P), ...
        if name in SAT_NAMES:
            return lambda P: self.sat(name, P)
        raise AttributeError(name)


def verify_fits(fluid: str = "Nitrogen", n_sat: int = 997, n_vap: int = 4000, seed: int = 1) -> dict:
    """Max relative error of every fit vs CoolProp on an independent grid
    (random P log-uniform; random T uniform in ln T between Tsat and Tmax)."""
    fp = fluid_props(fluid)
    rng = np.random.default_rng(seed)
    P = np.exp(rng.uniform(np.log(fp.Pmin), np.log(fp.Pmax), n_sat))
    ref = coolprop_sat(fluid, P)
    out = {"fluid": fluid, "P_range_kPa": [fp.Pmin / 1e3, fp.Pmax / 1e3], "Tmax": fp.Tmax, "sat": {}, "vap": {}}
    for k in SAT_NAMES:
        out["sat"][k] = float(np.max(np.abs(fp.sat(k, P) / ref[k] - 1)))
    CP = _cp()
    Pv = np.exp(rng.uniform(np.log(fp.Pmin), np.log(fp.Pmax), n_vap))
    Ts = np.array([CP.PropsSI("T", "P", float(p), "Q", 1, fluid) for p in Pv])
    u = rng.uniform(0, 1, n_vap)
    u[: n_vap // 10] = 0.0  # include the saturated-vapour line itself
    Tv = Ts * np.exp(u * np.log(fp.Tmax / Ts))
    refv = coolprop_vap(fluid, Tv, Pv)
    for k in VAP_NAMES:
        out["vap"][k] = float(np.max(np.abs(fp.vap(k, Tv, Pv) / refv[k] - 1)))
    out["max_sat"] = max(out["sat"].values())
    out["max_vap"] = max(out["vap"].values())
    return out
