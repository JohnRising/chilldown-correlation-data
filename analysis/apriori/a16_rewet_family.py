"""Phase 5: family of pressure-dependent, corresponding-states rewet rules; one constant each fitted on
the 11 UF low-g stations (leave-one-TEST-out); pressure slopes; predictions for NBS 9264 copper lines.

Every rule defines a liquid-side contact temperature T* and the wall temperature follows from the
semi-infinite contact relation  T_MFB = T* + beta (T* - T_l),  beta = e_l/e_w(T_MFB) (iterated),
T_l = Tsat (saturated inlet).  T_hn = Tc (0.905 + 0.095 Tr^8), Tr = Tsat/Tc (Lienhard 1976, reduced).
  R0  T* = Tc                                    (Phase-2 rule, 0 constants)
  R1a T* = T_hn + a (Tc - T_hn)                  (a = 1 -> R0)
  R1b T* = Tsat + c (T_hn - Tsat)
  R2  T* = Tsat + k (Tc - Tsat)                  (reduced linear)
  R3  Henry (1974) on T_B = T_hn:  T* = T_B + A (T_B - T_l) X^0.6, X = (e_l/e_w) h_fg/(c_w (T_B - Tsat))  (A fitted; Henry 0.42)
  R4  T* = Tsat + d Tc                           (constant reduced superheat; slope 1)
  R5  T* = T_hn + e Tc                           (constant reduced offset above T_hn)
"""
import json

import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar

from common import OUT, md_table, ss_cp, ss_e
from a09_rewet_theory import cu_e, cu_cp, liq

FLU = {"Nitrogen": {}, "Hydrogen": {}}
for f in FLU:
    FLU[f] = dict(Tc=CP.PropsSI("Tcrit", f), Pc=CP.PropsSI("pcrit", f))


def sat_props(P, fluid):
    g = lambda k, q=0: CP.PropsSI(k, "P", P, "Q", q, fluid)
    return dict(Tsat=g("T"), e_l=np.sqrt(g("L") * g("D") * g("C")), h_fg=g("H", 1) - g("H", 0))


def T_hn(Tsat, Tc):
    return Tc * (0.905 + 0.095 * (Tsat / Tc) ** 8)


def rule(name, p, P, fluid="Nitrogen", e_wall=ss_e, c_wall=ss_cp):
    Tc = FLU[fluid]["Tc"]
    out, betas = [], []
    for PP in np.atleast_1d(P):
        s = sat_props(PP, fluid)
        Ts, Th = s["Tsat"], T_hn(s["Tsat"], Tc)
        Tw = Tc
        for _ in range(60):
            ew = e_wall(Tw) if callable(e_wall) else e_wall
            b = s["e_l"] / ew
            if name == "R0": Ts_ = Tc
            elif name == "R1a": Ts_ = Th + p * (Tc - Th)
            elif name == "R1b": Ts_ = Ts + p * (Th - Ts)
            elif name == "R2": Ts_ = Ts + p * (Tc - Ts)
            elif name == "R3":
                X = b * s["h_fg"] / (c_wall(Tw) * (Th - Ts))
                Ts_ = Th + p * (Th - Ts) * X**0.6
            elif name == "R4": Ts_ = Ts + p * Tc
            elif name == "R5": Ts_ = Th + p * Tc
            Tn = Ts_ + b * (Ts_ - Ts)
            if abs(Tn - Tw) < 1e-4:
                break
            Tw = 0.5 * (Tw + Tn)
        out.append(Tw); betas.append(b)
    return np.array(out), np.array(betas)


RULES = {"R0": None, "R1a": (-1.0, 3.0), "R1b": (0.5, 4.0), "R2": (0.3, 3.0), "R3": (0.0, 3.0), "R4": (0.05, 0.8), "R5": (-0.2, 0.5)}


def fit(name, P, T):
    if RULES[name] is None:
        return np.nan
    return minimize_scalar(lambda p: np.mean((rule(name, p, P)[0] - T) ** 2), bounds=RULES[name], method="bounded").x


def slope(Tsat, y, dn):
    X = np.c_[np.ones(len(y)), Tsat, dn]
    return np.linalg.lstsq(X, y, rcond=None)[0][1]


def main():
    rt = pd.read_csv(f"{OUT}/rewet_theory_table.csv")
    P, T, Ts, dn = rt.P_kPa.values * 1e3, rt.T_mfb.values, rt.Tsat.values, (rt.st == "dn").values.astype(float)
    # data slopes: all, and P > 370 kPa subset (tests 8-10)
    rng = np.random.default_rng(5)
    def boot(mask):
        cs = rt.case[mask].unique(); bs = []
        for _ in range(2000):
            pick = rng.choice(cs, len(cs))
            ii = np.concatenate([np.where((rt.case.values == c) & mask)[0] for c in pick])
            if np.ptp(Ts[ii]) > 1 and len(np.unique(dn[ii])) == 2:
                bs.append(slope(Ts[ii], T[ii], dn[ii]))
        return np.percentile(bs, [2.5, 97.5]).tolist()
    hi = P > 365e3
    dslope = dict(all=slope(Ts, T, dn), all_ci=boot(np.ones(len(T), bool)),
                  hiP=slope(Ts[hi], T[hi], dn[hi]), hiP_ci=boot(hi))
    rows, params = [], {}
    for name in RULES:
        p = fit(name, P, T)
        pred = rule(name, p, P)[0]
        # LOTO
        cv = np.empty(len(T))
        for c in rt.case.unique():
            m = rt.case.values == c
            pc = fit(name, P[~m], T[~m])
            cv[m] = rule(name, pc, P[m])[0]
        Pg = np.array([1e5, 3e5, 7.5e5])
        Tg = rule(name, p, Pg)[0]
        e = pred - T
        rows.append(dict(rule=name, const=p, mean_err=e.mean(), MAE=np.abs(e).mean(), rmse=np.sqrt(np.mean(e**2)),
                         loto_rmse=np.sqrt(np.mean((cv - T) ** 2)), up_bias=e[dn == 0].mean(), dn_bias=e[dn == 1].mean(),
                         slope_vs_Tsat=slope(Ts, pred, dn), T_100kPa=Tg[0], T_300kPa=Tg[1], T_750kPa=Tg[2]))
        params[name] = p
    tab = pd.DataFrame(rows)
    # NBS predictions, copper wall
    nbs = []
    for fluid, Ps in (("Nitrogen", [1, 2.5, 3.4, 4.2, 5.9]), ("Hydrogen", [1, 2.5, 4.2, 5.1, 7.6, 11.0])):
        for Pa in Ps:
            row = dict(fluid=fluid, P_atm=Pa, Tsat=sat_props(Pa * 101325, fluid)["Tsat"])
            for name in RULES:
                Tw, b = rule(name, params[name], Pa * 101325, fluid, cu_e, cu_cp)
                row[name] = Tw[0]
            row["beta_Cu_R0"] = rule("R0", np.nan, Pa * 101325, fluid, cu_e, cu_cp)[1][0]
            nbs.append(row)
    nbs = pd.DataFrame(nbs)
    pd.set_option("display.width", 250)
    print(json.dumps(dslope)); print(tab.round(3).to_string()); print(nbs.round(2).to_string())
    json.dump({k: (None if np.isnan(v) else float(v)) for k, v in params.items()}, open(f"{OUT}/rewet_family_params.json", "w"), indent=1)
    with open(f"{OUT}/t16_rewet_family.md", "w") as f:
        f.write(__doc__ + "\ndata slope dT_MFB/dTsat (station dummy; test bootstrap): " + json.dumps(dslope)
                + "\n\n" + md_table(tab.round(3), 4) + "\n\n## predictions, copper wall (NBS 9264 drive pressures)\n\n"
                + md_table(nbs.round(2), 4) + "\n")


if __name__ == "__main__":
    main()
