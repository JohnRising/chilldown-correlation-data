"""Phase 3 / Task 1: position-dependent (developing vapour layer) film-boiling laws on the consistent
dataset (a07/a08).  All laws use film-temperature vapour properties and the lag-corrected G.

  h = C k_v (u_l rho_v cp_v / (k_v L_eff))^(1/3)        (local Leveque-type, D-free)
  L0  : L_eff = D                                (Phase-2 law, C = 8.55)
  L1  : L_eff = z            (inlet = leading edge, 1 constant)
  L2  : L_eff = z + z0       (z0 fitted, 2 constants)
  L3  : L_eff = z Lam/(z + Lam) (saturation at renewal length Lam, z0 = 0; 2 constants)
  L4  : L_eff = (z+z0) Lam/(z+z0+Lam) (3 constants)
  L5  : D-law x Darr entrance factor [1 + (D/(z+0.1D))^0.7]   (1 constant)
  Xm  : D-law x (1 - x_eq)^m  (liquid-core depletion; 2 constants)  [x built from measured q: endogenous]
"""
import json

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import D, OUT, md_table
from a08_film_laws import load_fb


def base(df, L):
    """Leveque-type with shear length D and development length L:  (k/D) Pe_l^(1/3) (D/L)^(1/3)."""
    ul = df.G / df.rho_l
    return df.k_f / D * (ul * D * df.rho_f * df.cp_f / df.k_f) ** (1 / 3) * (D / L) ** (1 / 3)


def pen(df, L):
    """D-free penetration (renewal) form sqrt(k rho cp u_l / L)."""
    return np.sqrt(df.k_f * df.rho_f * df.cp_f * df.G / df.rho_l / L)


E = np.exp
LAWS = {
    "L0 L=D": (lambda d, p: p[0] * base(d, D), [8.5]),
    "L1 L=z": (lambda d, p: p[0] * base(d, d.z), [8.5]),
    "L2 L=z+z0": (lambda d, p: p[0] * base(d, d.z + E(p[1])), [8.5, np.log(0.3)]),
    "L3 L=z*Lam/(z+Lam)": (lambda d, p: p[0] * base(d, d.z * E(p[1]) / (d.z + E(p[1]))), [8.5, np.log(0.3)]),
    "Lz power (D/z)^a": (lambda d, p: p[0] * base(d, D) * (D / d.z) ** p[1] * 30**p[1], [8.5, 0.1]),
    "L5 D-law x Darr entrance": (lambda d, p: p[0] * base(d, D) * (1 + (D / (d.z + 0.1 * D)) ** 0.7), [7.6]),
    "Xm D-law x (1-x)^m": (lambda d, p: p[0] * base(d, D) * (1 - np.clip(d.x, 0, 0.9)) ** p[1], [8.5, 1.0]),
    "L3x L3 with u_l=G(1-x)/rho_l [m=1/3, 0 extra]": (lambda d, p: p[0] * base(d, d.z * E(p[1]) / (d.z + E(p[1]))) * (1 - np.clip(d.x, 0, 0.9)) ** (1 / 3), [16.0, np.log(0.13)]),
    "L3m L3 x (1-x)^m": (lambda d, p: p[0] * base(d, d.z * E(p[1]) / (d.z + E(p[1]))) * (1 - np.clip(d.x, 0, 0.9)) ** p[2], [16.0, np.log(0.13), 0.5]),
    "P1 D-free sqrt(k rho cp u_l/(z+z0))": (lambda d, p: p[0] * pen(d, d.z + E(p[1])), [1.0, np.log(0.3)]),
    "P0 D-free sqrt(k rho cp u_l/z)": (lambda d, p: p[0] * pen(d, d.z), [1.0]),
}


def fit(df, f, p0):
    w = np.sqrt(df.w.values); y = np.log(df.h_E.values)
    return least_squares(lambda p: w * (np.log(np.clip(f(df, p), 1e-6, None)) - y), p0).x


def stats(df, r):
    w = df.w.values
    cs = pd.Series(r).groupby([df.case.values, df.st.values]).mean()
    updn = [cs[(c, "up")] - cs[(c, "dn")] for c in sorted(df.case.unique())]
    return dict(rms=float(np.sqrt(np.sum(w * r**2) / w.sum())), sd_cs=float(cs.std()),
                updn_mean=float(np.mean(updn)), updn_sd=float(np.std(updn, ddof=1)),
                updn_by_case=[round(v, 3) for v in updn])


def main():
    df = load_fb()
    rows, params = [], {}
    for name, (f, p0) in LAWS.items():
        p = fit(df, f, p0)
        r = np.log(f(df, p) / df.h_E.values)
        st = stats(df, r)
        rcv = np.empty(len(df)); ps = []
        for c in df.case.unique():
            m = df.case.values != c
            pc = fit(df[m].reset_index(drop=True), f, p); ps.append(pc)
            rcv[~m] = np.log(f(df[~m].reset_index(drop=True), pc) / df.h_E.values[~m])
        w = df.w.values
        ps = np.array(ps); pa = np.array(p)
        rows.append(dict(law=name, n_fit=len(p0), params=",".join(f"{v:.4g}" for v in pa),
                         loco_range=";".join(f"{a:.3g}-{b:.3g}" for a, b in zip(ps.min(0), ps.max(0))),
                         rms=st["rms"], loco=float(np.sqrt(np.sum(w * rcv**2) / w.sum())), sd_cs=st["sd_cs"],
                         updn_resid_mean=st["updn_mean"], updn_resid_sd=st["updn_sd"],
                         updn_by_case=str(st["updn_by_case"])))
        params[name] = pa.tolist()
    tab = pd.DataFrame(rows)
    # measured up/dn h ratio per case at matched times (both stations in film boiling)
    rat = []
    for c in sorted(df.case.unique()):
        u = df[(df.case == c) & (df.st == "up")].set_index("t"); d_ = df[(df.case == c) & (df.st == "dn")].set_index("t")
        tt = u.index.intersection(d_.index)
        rat.append(dict(case=c, G=u.G.loc[tt].median(), h_up_over_dn=float(np.exp(np.mean(np.log(u.h_E.loc[tt] / d_.h_E.loc[tt])))),
                        x_up=u.x.loc[tt].median(), x_dn=d_.x.loc[tt].median(),
                        dT_up_minus_dn=float((u.dT.loc[tt] - d_.dT.loc[tt]).median())))
    rat = pd.DataFrame(rat)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    print(tab.drop(columns="updn_by_case").round(3).to_string()); print(tab[["law", "updn_by_case"]].to_string())
    print(rat.round(3).to_string())
    json.dump(params, open(f"{OUT}/position_law_params.json", "w"), indent=1)
    with open(f"{OUT}/t11_position_laws.md", "w") as f:
        f.write(__doc__ + "\nupdn_resid = mean ln(pred/meas) upstream minus downstream (0 = no station bias; "
                "negative = law under-predicts upstream relative to downstream)\n\n")
        f.write(md_table(tab.round(3), 4) + "\n\n## measured h_up/h_dn at matched times\n\n" + md_table(rat.round(3), 3) + "\n")


if __name__ == "__main__":
    main()
