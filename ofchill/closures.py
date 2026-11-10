"""Library of closures + registry used by the CLI."""

from __future__ import annotations

from . import expr as E
from .closure import Closure, CustomClosure, NamedClosure, closure, NAMED_MODELS


def dittus_boelter(G, D, mu, k, cp, C=0.023, m=0.8, n=0.4):
    """OpenFLUME dittusBoelterH (v0.3.1): Nu = 3.66 for Re < 2000, C Re^m Pr^n
    for Re >= 4000, linear blend in between; h = Nu k / D."""
    Re = G * D / mu
    nu_t = C * E.maximum(Re, 1e-6) ** m * (cp * mu / k) ** n
    # == 3.66 (Re<2000), nu_t (Re>=4000), linear blend between (single nu_t copy)
    nu = 3.66 + (nu_t - 3.66) * E.clip((Re - 2000.0) / 2000.0, 0.0, 1.0)
    return nu * k / D


@closure(params={"C": 0.023, "m": 0.8, "n": 0.4}, label="db_vapor")
def db_vapor(s):
    """Reproduces the prior agent's PATCHED dispatcher (UF_DB_TWOPHASE=vapor):
    inside the dome Dittus-Boelter at full G with saturated-VAPOUR properties
    at P; outside the dome the shipped DB with the node's own properties."""
    C, m, n = s.param("C"), s.param("m"), s.param("n")
    P = s.P
    h_dome = dittus_boelter(s.G, s.D, s.fp.mu_v(P), s.fp.k_v(P), s.fp.cp_v(P), C, m, n)
    h_single = dittus_boelter(s.G, s.D, s.mu, s.k, s.cp, C, m, n)  # cp only valid outside the dome
    return E.where(s.in_dome(), h_dome, h_single)


@closure(
    params={"C": 0.023, "m": 0.8, "n": 0.4, "yC": 0.1, "yD": 0.4, "yQ": 0.4},
    label="miropolskii_expr",
)
def miropolskii_expr(s):
    """Expression re-implementation of OpenFLUME's named 'miropolskii' (v0.3.1):
    dome -> Miropolskii with sat-vapour/liquid props and x clipped to
    [0.01, 0.99]; single-phase -> DB with node properties.  Used to validate
    the closure compiler end-to-end against the named model."""
    P = s.P
    fp = s.fp
    x = E.clip(s.x_eq(), 0.01, 0.99)
    rl, rv, mug, kg, cpg = fp.rho_l(P), fp.rho_v(P), fp.mu_v(P), fp.k_v(P), fp.cp_v(P)
    Re_g = s.G * s.D / mug
    Pr_g = cpg * mug / kg
    Y = 1 - s.param("yC") * (rl / rv - 1) ** s.param("yD") * (1 - x) ** s.param("yQ")
    nu = s.param("C") * (Re_g * (x + rv / rl * (1 - x))) ** s.param("m") * Pr_g ** s.param("n") * Y
    h_dome = nu * kg / s.D
    h_single = dittus_boelter(s.G, s.D, s.mu, s.k, s.cp, s.param("C"), s.param("m"), s.param("n"))
    return E.where(s.in_dome(), h_dome, h_single)


@closure(params={"h0": 50.0, "a": 2.0}, label="probe_all")
def probe_all(s):
    """Test closure touching Tw, Tf, P, G, D, rho, mu, k and the fits (never
    cp/quality, so it is valid in every phase)."""
    return s.param("h0") + s.param("a") * (s.Tw - s.Tf) + 1e-3 * s.G * s.D / s.mu * s.k / s.D + 0.01 * s.rho \
        + 1e-4 * (s.P - 1e5) / 1e3 + 10 * s.fp.vap("k", 0.5 * (s.Tw + s.fp.Tsat(s.P)), s.P) / s.k


@closure(params={"h_film": 150.0, "h_nb": 5000.0, "dT_rewet": 40.0, "width": 0.0}, label="boil_demo")
def boil_demo(s):
    """Toy boiling curve for solver probing (NOT a physical model):
    h = h_film for Tw - Tsat > dT_rewet + width/2, h_nb below dT_rewet - width/2,
    smoothstep in between (width = 0 -> a sharp step, 33x jump)."""
    dT = s.Tw - s.fp.Tsat(s.P)
    w = s.param("width")
    f = E.smoothstep(s.param("dT_rewet") - 0.5 * w, s.param("dT_rewet") + 0.5 * w, dT)
    return s.param("h_nb") + (s.param("h_film") - s.param("h_nb")) * f


REGISTRY: dict[str, Closure] = {
    "boil_demo": boil_demo,
    "miropolskii": NamedClosure("miropolskii"),
    "dittusBoelter": NamedClosure("dittusBoelter"),
    "db_vapor": db_vapor,
    "miropolskii_expr": miropolskii_expr,
}


def get_closure(name: str, **params) -> Closure:
    if name in NAMED_MODELS:
        return NamedClosure(name, extra=params)
    c = REGISTRY[name]
    if params and isinstance(c, CustomClosure):
        c = c.with_params(**params)
    return c


# =============================================================================
# Phase 2: parametric unified boiling-curve closure
# =============================================================================
import json as _json  # noqa: E402
import math as _math  # noqa: E402
from pathlib import Path as _Path  # noqa: E402

from .closure import CustomClosure as _CC  # noqa: E402

_FB_CONST_FILE = _Path(__file__).resolve().parent / "fb_constants.json"
N2_TC = 126.192  # K   (CoolProp Nitrogen)
N2_PC = 3.3958e6  # Pa
N2_M = 28.0134  # g/mol
FB_OPTIONS = ("pe_l", "renewal_l", "renewal_v", "darr2016", "bromley", "pe_lz")
TMFB_OPTIONS = ("const", "contact_Tc", "contact_Tl", "Thn_plus", "R5")
NB_OPTIONS = ("power", "cooper", "cooper_add")
CHF_OPTIONS = ("const", "katto")


def fb_constants() -> dict:
    """A-priori film-boiling constants (ofchill/fit_fb.py -> fb_constants.json)."""
    if _FB_CONST_FILE.exists():
        return _json.loads(_FB_CONST_FILE.read_text())
    return {}


def _default_fb_params(fb: str, props: str) -> dict:
    c = fb_constants()
    if fb == "pe_l":
        e = c.get(f"pe_l_{props}", {"C": 7.78 if props == "film" else 12.93, "n": 0.345 if props == "film" else 0.293})
        return {"C_fb": e["C"], "n_fb": e["n"]}
    if fb in ("renewal_l", "renewal_v", "bromley"):
        dflt = {"renewal_l": 2.65, "renewal_v": 0.226, "bromley": 1.89}[fb]
        return {"C_fb": c.get(f"{fb}_{props}", {"C": dflt})["C"]}
    if fb == "darr2016":
        return {"C_fb": 1.0}  # published constants (x 1.19 fits the energy-consistent h)
    if fb == "pe_lz":  # CLOSURE_SPEC §10 (developing Leveque layer, saturating at Lambda)
        return {"C_fb": 15.8, "Lam": 0.101, "x_ul_max": 0.9}
    raise ValueError(fb)


def ss304_effusivity(T: float, rho: float = 8027.0) -> float:
    """sqrt(k rho cp) of 304SS from the NIST cryogenic fits (as in the solver presets)."""
    x = _math.log10(T)
    pk = [-1.4087, 1.3982, 0.2543, -0.626, 0.2334, 0.4256, -0.4658, 0.165, -0.0199]
    pc = [22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212, 24.7593, -2.239153, 0]
    k = 10 ** sum(c * x**i for i, c in enumerate(pk))
    cp = 10 ** sum(c * x**i for i, c in enumerate(pc))
    return _math.sqrt(k * rho * cp)


def contact_beta(P_mean: float, T_wall: float = 130.0, fluid: str = "Nitrogen") -> float:
    """beta = e_l / e_w with e = sqrt(k rho cp): sat. liquid at P_mean, 304SS at T_wall."""
    import CoolProp.CoolProp as CP

    e_l = _math.sqrt(CP.PropsSI("L", "P", P_mean, "Q", 0, fluid) * CP.PropsSI("D", "P", P_mean, "Q", 0, fluid)
                     * CP.PropsSI("C", "P", P_mean, "Q", 0, fluid))
    return e_l / ss304_effusivity(T_wall)


pi_ = _math.pi


def make_boiling_closure(
    fb: str = "pe_l",
    props: str = "sat",
    tmfb: str = "const",
    nb: str = "power",
    gravity: bool | str = False,
    vapor_film: bool = False,
    chf: str = "const",
    label: str | None = None,
    fluid: str = "Nitrogen",
    strat: bool = False,
    no_dbx: bool = False,
    **params,
) -> _CC:
    """Unified, q-continuous boiling curve as ONE custom expression.

    Structural (Python) options -> expression shape:
      fb    film-boiling law  'pe_l' | 'renewal_l' | 'renewal_v' | 'darr2016' | 'bromley'
      props vapour props for fb: 'sat' (sat. vapour at P; small expression) | 'film' ((Tw+Tsat)/2 fits)
      tmfb  rewet temperature 'const' (T_MFB) | 'contact_Tc' (Tc + beta (Tc - Tf)) | 'Thn_plus' (Lienhard T_hn(Tsat) + dThn)
      nb    'power' q = q_CHF (dT/dT_CHF)^2.5 | 'cooper' (Cooper 1984 pool NB, h-form; dT_CHF where q = q_CHF)
      gravity  add a Bromley natural-convection film term, cube-sum (h_f^3 + h_b^3)^(1/3), param g [m/s^2]:
               True/'horizontal': horizontal tube, h_b = 0.62 [k^3 rho (rl-rv) g h_fg'/(D mu dT)]^1/4
               'vertical': vertical wall (vapour film rises along the wall), 0.943 [...]^1/4 with the
               length L = max(z, D) instead of D (z = distance from the tube inlet)
      fluid  'Nitrogen' | 'Hydrogen' | 'Oxygen' (property fits; Tc, Pc, M defaults from CoolProp)
      strat  (needs gravity='horizontal') 1-g stratification: only the wetted bottom fraction f_w of the
             perimeter is in film boiling.  eps_l = 1 - alpha_Zivi(x), f_w = delta/pi with Biberg's (1999)
             wetted half-angle delta(eps_l); blended to full wetting by the superficial-liquid Froude number
             Fr = G(1-x)/(rho_l sqrt(g D)): f = s + (1-s) f_w, s = Fr^2/(1+Fr^2).  No fitted constant.
      no_dbx drop the vapour-DB(G x) floor in two-phase nodes (HEM forces that heat into evaporation;
             physically it superheats the vapour) - the film law then applies up to x_dry
      vapor_film  superheated-vapour nodes use max(DB vapour, film law if Tw > T_MFB) (default: DB vapour only)
      chf   'const' (q_CHF param) | 'katto' (Katto-Ohno 1984, sat. inlet, L = max(z, z_min); q_CHF >= q_MFB)
      tmfb 'contact_Tl': CLOSURE_SPEC form Tc + beta (Tc - min(Tf, Tsat))
      nb   'cooper_add': CLOSURE_SPEC form q = h_DBl dT + (A dT)^3, dT_CHF from the cubic (Cardano)
      Nu_floor (param, default 0): h_film >= Nu_floor k_v/D (CLOSURE_SPEC: 4.36)
    Numeric params (correlation.params; override by keyword): C_fb, n_fb, T_MFB, beta, dThn, Tc,
      dT_CHF, q_CHF, m, x_dry, g, Pc, M, z (axial position, set per conductor by the builder).

    Regimes (dT = Tw - Tsat(P), all fluxes on the Tsat basis, h = q / (Tw - Tf)):
      vapour node (rho < rho_v(P)):  h = DB(G, node mu, k, cp)
      otherwise (two-phase / subcooled liquid; x = x_eq from rho):
        h_FB(Tw)  = max(h_film(Tw) * (1 - smoothstep(x_dry, 1, x)), DB(G x, sat-vapour props))
        Tw >= T_MFB*:          q = h_FB(Tw) dT
        T_CHF < Tw < T_MFB*:   q = q_MFB + (q_CHF - q_MFB) theta^m,  q_MFB = h_FB(T_MFB*) (T_MFB* - Tsat)
        Tw <= T_CHF:           q = max(q_NB(dT), DB(G (1-x), sat-liquid props) (Tw - Tf))
        dT <= 0.01 K:          h = DB liquid
      T_CHF = min(Tsat + dT_CHF, T_MFB - 2 K (1 K for katto)); dT_CHF = param (power), q_CHF^0.33/A (cooper),
      cubic root (cooper_add).  (T_MFB* below = T_MFB.)
    q is continuous at T_MFB* and at T_CHF by construction (h therefore too).
    """
    if (fb not in FB_OPTIONS or tmfb not in TMFB_OPTIONS or nb not in NB_OPTIONS or props not in ("sat", "film")
            or chf not in CHF_OPTIONS):
        raise ValueError((fb, props, tmfb, nb, chf))
    if gravity is True:
        gravity = "horizontal"
    if gravity not in (False, "horizontal", "vertical"):
        raise ValueError(gravity)
    Tc0, Pc0, M0 = fluid_crit(fluid)
    p = {"T_MFB": 134.0, "beta": 0.09, "dThn": 18.8, "Tc": Tc0, "dT_CHF": 24.0, "q_CHF": 1.0e5, "m": 0.5,
         "x_dry": 0.95, "Pc": Pc0, "M": M0}
    p.update(_default_fb_params(fb, props))
    if fb in ("darr2016", "pe_lz") or chf == "katto" or gravity == "vertical":
        p["z"] = 0.0
    if chf == "katto":
        p["z_min"] = 2 * 0.01168
    if gravity:
        p["g"] = 0.0
    if tmfb == "R5":
        p["e_hn"] = 0.122
    p.update(params)

    def fn(s):
        fp, P, Tw, Tf, G, D = s.fp, s.P, s.Tw, s.Tf, s.G, s.D
        par = s.param
        Tsat = fp.Tsat(P)
        rl, rv = fp.rho_l(P), fp.rho_v(P)
        hfg = fp.h_fg(P)
        x = s.x_eq()
        mug, kg, cpg = fp.mu_v(P), fp.k_v(P), fp.cp_v(P)

        def h_film(Twx):
            dT = E.maximum(Twx - Tsat, 0.5)
            if props == "film":
                Tfi = 0.5 * (Twx + Tsat)
                kv, rhov, cpv, muv = fp.vap("k", Tfi, P), fp.vap("rho", Tfi, P), fp.vap("cp", Tfi, P), mug
            else:
                kv, rhov, cpv, muv = kg, rv, cpg, mug
            ul = G / rl
            C = par("C_fb")
            if fb == "pe_lz":
                # h = C k^(2/3) (u_l rho cp)^(1/3) D^(-1/3) (1/z + 1/Lam)^(1/3), u_l = G (1 - x)/rho_l
                ulx = G * (1 - E.clip(x, 0.0, par("x_ul_max"))) / rl
                zz = E.maximum(par("z"), D)
                h = C * kv ** (2.0 / 3.0) * (ulx * rhov * cpv) ** (1.0 / 3.0) * D ** (-1.0 / 3.0) \
                    * (1.0 / zz + 1.0 / par("Lam")) ** (1.0 / 3.0)
            elif fb == "pe_l":
                n_ = par("n_fb")
                h = C * kv ** (1 - n_) * (ul * D * rhov * cpv) ** n_ / D  # == C kv/D Pe^n, kv once
            elif fb == "renewal_l":
                h = C * E.sqrt(kv * rhov * cpv * ul / D)
            elif fb == "renewal_v":
                h = C * E.sqrt(kv * cpv * G / D)
            elif fb == "bromley":
                h = C * E.sqrt(kv * rhov * (hfg + 0.5 * cpv * dT) * ul / (D * dT))
            else:  # darr2016 low-g (Darr et al. 2016 eqs 1-4; c7 = 0), sat. props, L = z
                L = par("z")
                xc = E.clip(x, 0.0, 0.99)
                Re_v = G * D / mug
                Pr_v = cpg * mug / kg
                htc = 0.150 * (1 + (D / (L + 0.1 * D)) ** 0.7) * (1 - xc) ** (-0.76) * Re_v**0.52 * Pr_v**0.4 * kg / D
                We = G * G * D / (rl * fp.sigma(P))
                Ja = fp.cp_l(P) * dT / hfg
                hdw = 5e-6 * fp.k_l(P) / D * E.exp(-0.10 * L / D) * We**2.0 / Ja
                h = C * (htc + hdw)
            if "Nu_floor" in p:
                h = E.maximum(h, par("Nu_floor") * kv / D)
            if gravity:
                cg, Lg = (0.62, D) if gravity == "horizontal" else (0.943, E.maximum(par("z"), D))
                hb = cg * (kv**3 * rhov * (rl - rhov) * par("g") * (hfg + 0.4 * cpv * dT) / (Lg * muv * dT)) ** 0.25
                h = (h**3 + hb**3) ** (1.0 / 3.0)
            return h

        h_dbv = dittus_boelter(G * x, D, mug, kg, cpg)
        taper = 1 - E.smoothstep(par("x_dry"), 1.0, x)
        if strat:
            if gravity != "horizontal":
                raise ValueError("strat needs gravity='horizontal'")
            xz = E.clip(x, 1e-6, 1.0)
            alpha = 1.0 / (1.0 + (1.0 - xz) / xz * (rv / rl) ** (2.0 / 3.0))
            el = E.clip(1.0 - alpha, 0.0, 1.0)
            delta = pi_ * el + (1.5 * pi_) ** (1.0 / 3.0) * (1.0 - 2.0 * el + el ** (1.0 / 3.0) - (1.0 - el) ** (1.0 / 3.0))
            fw = E.clip(delta / pi_, 0.0, 1.0)
            Fr = G * (1.0 - E.clip(x, 0.0, 1.0)) / (rl * E.sqrt(E.maximum(par("g"), 1e-9) * D))
            sF = Fr * Fr / (1.0 + Fr * Fr)
            fwet = sF + (1.0 - sF) * fw
        else:
            fwet = 1.0

        def h_FB(Twx):
            hf = h_film(Twx) * taper * fwet
            return hf if no_dbx else E.maximum(hf, h_dbv)

        # rewet / CHF temperatures
        if tmfb == "const":
            Tm = par("T_MFB")
        elif tmfb == "contact_Tc":
            Tm = par("Tc") + par("beta") * (par("Tc") - Tf)
        elif tmfb == "contact_Tl":
            Tm = par("Tc") + par("beta") * (par("Tc") - E.minimum(Tf, Tsat))
        elif tmfb == "R5":  # CLOSURE_SPEC §11: T* = T_hn(Lienhard) + e Tc ; T_MFB = T* + beta (T* - T_l)
            Trs = Tsat / par("Tc")
            Tstar = par("Tc") * (0.905 + 0.095 * Trs**8) + par("e_hn") * par("Tc")
            Tm = Tstar + par("beta") * (Tstar - E.minimum(Tf, Tsat))
        else:
            Tr = Tsat / par("Tc")
            Tm = par("Tc") * (0.905 + 0.095 * Tr**8) + par("dThn")
        dTw = Tw - Tsat
        rl_mu, kl, cpl = fp.mu_l(P), fp.k_l(P), fp.cp_l(P)
        h_dbl = dittus_boelter(G * (1 - x), D, rl_mu, kl, cpl)
        dTf = Tw - Tf
        if chf == "const":
            qchf = par("q_CHF")
        else:  # Katto-Ohno (1984), R < 0.15 branch, saturated inlet, heated length L
            L = E.maximum(par("z"), par("z_min"))
            Gc = E.maximum(G, 5.0)
            R = rv / rl
            l_ = L / D
            W = fp.sigma(P) * rl / (Gc * Gc * L)
            CK = E.clip(0.25 + 0.0009 * (l_ - 50.0), 0.25, 0.34)
            q1 = CK * W**0.043 / l_
            q2 = 0.10 * R**0.133 * W ** (1.0 / 3.0) / (1 + 0.0031 * l_)
            q3 = 0.098 * R**0.133 * W**0.433 * l_**0.27 / (1 + 0.0031 * l_)
            qchf = Gc * hfg * E.where(q1 < q2, q1, E.minimum(q2, q3))
        pr = P / par("Pc")
        A = 55.0 * pr**0.12 * (-E.log(pr) / _math.log(10.0)) ** (-0.55) * par("M") ** (-0.5)
        dTp = E.clip(dTw, 0.0, 1e3)
        if nb == "power":
            dTchf = par("dT_CHF")
            q_nb = qchf * (dTp / dTchf) ** 2.5
        elif nb == "cooper":
            dTchf = qchf**0.33 / A
            q_nb = A ** (1 / 0.33) * dTp ** (1 / 0.33)
        else:  # cooper_add: q = h_DBl dT + (A dT)^3 ; dT_CHF = real root of A^3 dT^3 + h dT - q_CHF (Cardano)
            A3 = A * A * A
            a_ = 0.5 * qchf / A3
            sq = E.sqrt(a_ * a_ + (h_dbl / (3.0 * A3)) ** 3)
            dTchf = (a_ + sq) ** (1.0 / 3.0) - E.maximum(sq - a_, 0.0) ** (1.0 / 3.0)
            q_nb = h_dbl * dTp + (A * dTp) ** 3
        # ordering guard T_CHF <= T_MFB - gap applied to T_CHF (keeps the large T_CHF expression out of
        # every T_MFB occurrence; T_MFB itself is never modified)
        Tchf = E.minimum(Tsat + dTchf, Tm - (1.0 if chf == "katto" else 2.0))
        q_mfb = h_FB(Tm) * (Tm - Tsat)
        theta = E.clip((Tm - Tw) / (Tm - Tchf), 0.0, 1.0)
        thm = theta ** par("m")
        q_tr = qchf * thm + q_mfb * (1 - thm)  # == q_MFB + (q_CHF - q_MFB) theta^m
        q = E.where(Tw >= Tm, h_FB(Tw) * dTw, E.where(Tw > Tchf, q_tr, q_nb))
        h_two = E.where(dTw > 0.01, q / E.maximum(dTf, 0.01), h_dbl)
        # superheated vapour node (no liquid): single-phase DB with the node's own props
        h_vap = dittus_boelter(G, D, s.mu, s.k, s.cp)
        if vapor_film:
            h_vap = E.maximum(h_vap, E.where(Tw > Tm, h_film(Tw), 0.0))
        is_vap = s.rho < rv * (1 - 2e-4)
        return E.where(is_vap, h_vap, h_two)

    lab = label or f"{fb}-{props}_{tmfb}_{nb}" + ("_g" if gravity else "")
    return _CC(fn, p, fluid, lab)


def fluid_crit(fluid: str) -> tuple[float, float, float]:
    """(Tc [K], Pc [Pa], M [kg/kmol]) from CoolProp."""
    if fluid == "Nitrogen":
        return N2_TC, N2_PC, N2_M
    import CoolProp.CoolProp as CP

    return CP.PropsSI("Tcrit", fluid), CP.PropsSI("pcrit", fluid), 1e3 * CP.PropsSI("molar_mass", fluid)


def copper_effusivity(T: float, rho: float = 8960.0, rrr: int = 100) -> float:
    """sqrt(k rho cp) of OFHC copper, NIST fits (k: RRR 100 rational form; cp: log-poly), as the
    OpenFLUME 'ofhc-copper' preset."""
    a = [2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281, 0.003207]
    t = T
    num = a[0] + a[2] * t**0.5 + a[4] * t + a[6] * t**1.5 + a[8] * t**2
    den = 1 + a[1] * t**0.5 + a[3] * t + a[5] * t**1.5 + a[7] * t**2
    k = 10 ** (num / den)
    pc = [-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0]
    x = _math.log10(T)
    cp = 10 ** sum(c * x**i for i, c in enumerate(pc))
    return _math.sqrt(k * rho * cp)


def contact_beta_spec(P: float, fluid: str = "Nitrogen", wall: str = "ss304", T_l: float | None = None) -> tuple[float, float]:
    """CLOSURE_SPEC §2: beta = e_l / e_w(T_MFB), T_MFB = Tc + beta (Tc - T_l), fixed point.
    Returns (beta, T_MFB).  T_l defaults to Tsat(P)."""
    import CoolProp.CoolProp as CP

    Tc = fluid_crit(fluid)[0]
    if T_l is None:
        T_l = CP.PropsSI("T", "P", P, "Q", 0, fluid)
    g = lambda k: CP.PropsSI(k, "P", P, "Q", 0, fluid)  # noqa: E731
    e_l = _math.sqrt(g("L") * g("D") * g("C"))
    ew = ss304_effusivity if wall == "ss304" else copper_effusivity
    Tm = Tc + 0.1 * (Tc - T_l)
    for _ in range(6):
        b = e_l / ew(Tm)
        Tm = Tc + b * (Tc - T_l)
    return b, Tm


def boiling_variant(name: str, **params) -> _CC:
    """Named phase-2 variants (see README)."""
    V = {
        "a": dict(fb="pe_l", props="film", tmfb="const", nb="power"),
        "a_sat": dict(fb="pe_l", props="sat", tmfb="const", nb="power"),
        "b": dict(fb="pe_l", props="film", tmfb="contact_Tc", nb="power"),
        "c_rl": dict(fb="renewal_l", props="sat", tmfb="const", nb="power"),
        "c_rv": dict(fb="renewal_v", props="sat", tmfb="const", nb="power"),
        "d": dict(fb="darr2016", props="sat", tmfb="const", nb="power"),
        "d119": dict(fb="darr2016", props="sat", tmfb="const", nb="power", C_fb=1.1864),
        "e_brom": dict(fb="bromley", props="film", tmfb="const", nb="power"),
        "a_cooper": dict(fb="pe_l", props="film", tmfb="const", nb="cooper"),
        "a_thn": dict(fb="pe_l", props="film", tmfb="Thn_plus", nb="power"),
        # analysis/apriori/CLOSURE_SPEC.md (recommended): Pe^(1/3), C = 8.55, film props, laminar floor,
        # contact-T_c rewet with T_l = min(Tf, Tsat), Katto-Ohno CHF, Cooper(additive) NB, n = 0.5
        "spec": dict(fb="pe_l", props="film", tmfb="contact_Tl", nb="cooper_add", chf="katto",
                     C_fb=8.55, n_fb=1.0 / 3.0, Nu_floor=4.36),
        # CLOSURE_SPEC §10: position-dependent film law (C 15.8, Lambda 0.101 m) + spec §2-6
        "spec10": dict(fb="pe_lz", props="film", tmfb="contact_Tl", nb="cooper_add", chf="katto", Nu_floor=4.36),
        # §10 film law + §11 rewet rule R5 (final default from Phase 6)
        "spec11": dict(fb="pe_lz", props="film", tmfb="R5", nb="cooper_add", chf="katto", Nu_floor=4.36),
        "spec11s": dict(fb="pe_lz", props="film", tmfb="R5", nb="cooper_add", chf="katto", Nu_floor=4.36, strat=True),
        "spec11sn": dict(fb="pe_lz", props="film", tmfb="R5", nb="cooper_add", chf="katto", Nu_floor=4.36, strat=True,
                         no_dbx=True),
        "spec11n": dict(fb="pe_lz", props="film", tmfb="R5", nb="cooper_add", chf="katto", Nu_floor=4.36, no_dbx=True),
        # Phase 8 rows with the final R5 rewet (film law is the only difference to spec11)
        "specz11": dict(fb="pe_l", props="film", tmfb="R5", nb="cooper_add", chf="katto",
                        C_fb=8.55, n_fb=1.0 / 3.0, Nu_floor=4.36),
        "ds11": dict(fb="darr2016", props="sat", tmfb="R5", nb="cooper_add", chf="katto", C_fb=1.0),
        "d119s11": dict(fb="darr2016", props="sat", tmfb="R5", nb="cooper_add", chf="katto", C_fb=1.1864),
        # Darr 2016 x 1.19 film boiling with the spec's rewet / transition / NB / CHF (fluid/wall general)
        "d119s": dict(fb="darr2016", props="sat", tmfb="contact_Tl", nb="cooper_add", chf="katto", C_fb=1.1864),
    }
    kw = dict(V[name])
    kw.update(params)
    lab = f"v{name}"
    return make_boiling_closure(label=lab, **kw)
