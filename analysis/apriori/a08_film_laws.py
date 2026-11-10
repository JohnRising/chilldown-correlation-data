"""Phase 2 / Task 2: film-boiling law discrimination on the consistent wall-energy-balance dataset
(out/q_consistent.csv from a07).  Film-boiling samples: phase-1 film mask (t>0.5 s, before the
earliest MFB onset - 0.5 s, Tw > T_MFB + 15 K).  Equal weight per case-station.  LOCO = leave one
CASE (both stations) out.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import D, FIG, OUT, md_table, sat, vap

COL = plt.cm.viridis(np.linspace(0, 0.95, 10))


def load_fb():
    df = pd.read_csv(f"{OUT}/q_consistent.csv")
    df = df[df.fb & (df.t > 0.5) & (df.q_E > 0) & (df.dT > 5)].reset_index(drop=True)
    s = sat(df.P.values)
    for k, v in s.items():
        df[k] = v
    vf = vap(df.P.values, 0.5 * (df.Tw.values + df.Tsat.values))
    for k, v in vf.items():
        df[f"{k}_f"] = v
    df["w"] = 1.0 / df.groupby(["case", "st"]).t.transform("size")
    return df


def props(df, which):
    if which == "sat":
        return df.k_v, df.rho_v, df.cp_v, df.mu_v
    return df.k_f, df.rho_f, df.cp_f, df.mu_f


def unit_forms(df, G):
    """1-constant forms: name -> unit (h = C * unit)."""
    U = {}
    ul = G / df.rho_l
    for pw in ("sat", "film"):
        k, rho, cp, mu = props(df, pw)
        U[f"F1 sqrt(k rho cp u_l/D) [{pw}]"] = np.sqrt(k * rho * cp * ul / D)
        U[f"F2 sqrt(k cp G/D) [{pw}]"] = np.sqrt(k * cp * G / D)
        U[f"F3 (k/D) Pe_l^1/3 [{pw}]"] = k / D * (ul * D * rho * cp / k) ** (1 / 3)
    U["F5 Darr2016 exp[-(L/D)^c5] (x C)"] = darr(df, G, corrected=True)
    U["F5' Darr2016 exp(-c5 L/D) (x C)"] = darr(df, G, corrected=False)
    U["Bromley sqrt(k rho hfg' u_l/(D dT)) [film]"] = np.sqrt(df.k_f * df.rho_f * (df.h_fg + 0.5 * df.cp_f * df.dT) * ul / (D * df.dT))
    return U


def darr(df, G, corrected=True):
    c1, c2, c3, c4, c5, c6 = 0.150, -0.76, 0.52, 5e-6, 0.10, 2.0
    L = df.z
    x = np.clip(df.x, 0, 0.99)
    Re = G * D / df.mu_v
    Pr = df.cp_v * df.mu_v / df.k_v
    htc = c1 * (1 + (D / (L + 0.1 * D)) ** 0.7) * (1 - x) ** c2 * Re**c3 * Pr**0.4 * df.k_v / D
    We = G**2 * D / (df.rho_l * df.sigma)
    Ja = df.cp_l * df.dT / df.h_fg
    ex = np.exp(-(L / D) ** c5) if corrected else np.exp(-c5 * L / D)
    return htc + c4 * df.k_l / D * ex * We**c6 / Ja


def trends(df, r):
    X = np.c_[np.ones(len(df)), np.log(df.Gu / 150), np.log(df.P / 3e5), np.log(df.dT / 100),
              (df.st == "dn").astype(float), df.x]
    sw = np.sqrt(df.w.values)
    b = np.linalg.lstsq(X * sw[:, None], r * sw, rcond=None)[0]
    key = df.case.astype(str) + df.st
    rd = r - pd.Series(r).groupby(key.values).transform("mean").values
    td = df.t.values - df.groupby(key.values).t.transform("mean").values
    return dict(tr_lnG=b[1], tr_lnP=b[2], tr_lndT=b[3], tr_dn=b[4], tr_x=b[5],
                tr_per10s=float(np.sum(rd * td) / np.sum(td * td) * 10))


def wstats(df, r):
    w = df.w.values
    cs = df.groupby(["case", "st"]).apply(lambda g: np.exp(np.mean(r[g.index])) - 1, include_groups=False)
    return dict(bias=float(np.exp(np.sum(w * r) / w.sum()) - 1), rms=float(np.sqrt(np.sum(w * r**2) / w.sum())),
                sd_cs_bias=float(np.std(np.log1p(cs.values), ddof=1)),
                w30=float(np.sum(w * (np.abs(np.exp(r) - 1) <= 0.3)) / w.sum()))


def eval_one_const(df, unit, h):
    lr = np.log(h / unit)
    w = df.w.values
    C = np.exp(np.sum(w * lr) / w.sum())
    r = np.log(C * unit / h)
    rcv = np.empty(len(df)); Cs = []
    for c in df.case.unique():
        tr = df.case.values != c
        Cc = np.exp(np.sum(w[tr] * lr[tr]) / w[tr].sum()); Cs.append(Cc)
        rcv[~tr] = np.log(Cc * unit[~tr] / h[~tr])
    out = dict(C=C, C_loco_sd=float(np.std(Cs, ddof=1) * np.sqrt((len(Cs) - 1) ** 2 / len(Cs))),
               C_loco_min=min(Cs), C_loco_max=max(Cs), **wstats(df, r),
               loco_rms=float(np.sqrt(np.sum(w * rcv**2) / w.sum())))
    out.update(trends(df, r))
    return out, C * unit


def eval_free(df, f, p0, h, names):
    w = df.w.values; sw = np.sqrt(w)
    fit = lambda m, p0: least_squares(lambda p: sw[m] * (np.log(f(df[m], p)) - np.log(h[m])), p0).x
    allm = np.ones(len(df), bool)
    p = fit(allm, p0)
    r = np.log(f(df, p) / h)
    rcv = np.empty(len(df)); ps = []
    for c in df.case.unique():
        m = df.case.values != c
        pc = fit(m, p); ps.append(pc)
        rcv[~m] = np.log(f(df[~m], pc) / h[~m])
    ps = np.array(ps)
    out = dict(C=",".join(f"{n}={v:.4g}" for n, v in zip(names, p)),
               C_loco_min=";".join(f"{v:.3g}" for v in ps.min(0)), C_loco_max=";".join(f"{v:.3g}" for v in ps.max(0)),
               **wstats(df, r), loco_rms=float(np.sqrt(np.sum(w * rcv**2) / w.sum())))
    out.update(trends(df, r))
    return out


def p_profile(df, h, Gcol):
    """profile of a pressure exponent p at fixed G exponent 0.5 on sat-vapour F2 and on the bare
    dimensional law h = A G^0.5 P^p; case-bootstrap CI; leverage from cases 3/4/5 vs 6."""
    w = df.w.values; sw = np.sqrt(w)
    res = {}
    for lab, base in (("h/G^0.5", np.log(h) - 0.5 * np.log(df[Gcol])),
                      ("h/F2sat-unit", np.log(h / np.sqrt(df.k_v * df.cp_v * df[Gcol] / D)))):
        ps = np.linspace(-0.4, 1.0, 141)
        rss = []
        for p in ps:
            y = base - p * np.log(df.P / 3e5)
            X = np.c_[np.ones(len(df)), (df.st == "dn")]
            b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
            rss.append(np.sum(w * (y - X @ b) ** 2))
        rss = np.array(rss)
        # effective number of independent units = case-stations (20); profile CI via F-like criterion
        ncs = df.groupby(["case", "st"]).ngroups
        dev = ncs * np.log(rss / rss.min())
        ci = ps[dev <= 3.84]
        rng = np.random.default_rng(1)
        cases = df.case.unique()
        boot = []
        for _ in range(400):
            pick = rng.choice(cases, len(cases))
            idx = np.concatenate([np.where(df.case.values == c)[0] for c in pick])
            X = np.c_[np.ones(len(idx)), np.log(df.P.values[idx] / 3e5), (df.st.values[idx] == "dn")]
            b = np.linalg.lstsq(X * sw[idx, None], base.values[idx] * sw[idx], rcond=None)[0]
            boot.append(b[1])
        res[lab] = dict(p_hat=float(ps[np.argmin(rss)]), ci_profile=[float(ci.min()), float(ci.max())],
                        ci_boot=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                        ps=ps.tolist(), dev=dev.tolist())
        # leverage pairs
        lev = {}
        m = lambda c: np.exp(np.average(base[df.case == c], weights=w[df.case == c]))
        Pm = lambda c: np.average(df.P[df.case == c], weights=w[df.case == c])
        for grp in ([3], [4], [5], [3, 4, 5]):
            num = np.log(m(6) / np.exp(np.mean([np.log(m(c)) for c in grp])))
            den = np.log(Pm(6) / np.exp(np.mean([np.log(Pm(c)) for c in grp])))
            lev["6 vs " + "/".join(map(str, grp))] = float(num / den)
        res[lab]["leverage_p"] = lev
    return res


def implied_P_exponent(df, G):
    """d ln(unit)/d ln P at fixed G, T_w for each form, median over data (finite difference)."""
    out = {}
    d2 = df.copy(); d2["P"] = df.P * 1.1
    s = sat(d2.P.values)
    for k, v in s.items():
        d2[k] = v
    d2["Tsat"] = s["Tsat"]; d2["dT"] = df.Tw - d2.Tsat
    vf = vap(d2.P.values, 0.5 * (df.Tw.values + d2.Tsat.values))
    for k, v in vf.items():
        d2[f"{k}_f"] = v
    U1, U2 = unit_forms(df, G), unit_forms(d2, G)
    for k in U1:
        out[k] = float(np.median(np.log(U2[k] / U1[k]) / np.log(1.1)))
    return out


def main():
    df = load_fb()
    print("n samples", len(df), "case-stations", df.groupby(["case", "st"]).ngroups)
    h = df.h_E.values
    allrows, store = [], {}
    for Gcol, lab in (("G", "tau=2.5 s"), ("G_raw", "tau=0")):
        df["Gu"] = df[Gcol]
        G = df[Gcol].values
        U = unit_forms(df, G)
        for name, unit in U.items():
            o, pred = eval_one_const(df, unit.values if hasattr(unit, "values") else unit, h)
            allrows.append(dict(G=lab, law=name, n_fit=1, **o))
            store[(lab, name)] = pred
        # free exponent forms
        for pw in ("sat", "film"):
            kk = "k_v" if pw == "sat" else "k_f"
            rr, cc = ("rho_v", "cp_v") if pw == "sat" else ("rho_f", "cp_f")
            mu = "mu_v" if pw == "sat" else "mu_f"
            f3 = lambda d, p, kk=kk, rr=rr, cc=cc: p[0] * d[kk] / D * (d.Gu / d.rho_l * D * d[rr] * d[cc] / d[kk]) ** p[1]
            allrows.append(dict(G=lab, law=f"F3n C (k/D) Pe_l^n [{pw}]", n_fit=2, **eval_free(df, f3, [1.0, 0.4], h, ["C", "n"])))
            f4 = lambda d, p, kk=kk, cc=cc, mu=mu: p[0] * (d.Gu * D / d[mu]) ** p[1] * (d[cc] * d[mu] / d[kk]) ** p[2] * d[kk] / D
            allrows.append(dict(G=lab, law=f"F4 C Re_v^n Pr^m [{pw}]", n_fit=3, **eval_free(df, f4, [0.2, 0.5, 0.5], h, ["C", "n", "m"])))
    tab = pd.DataFrame(allrows)
    tab.to_csv(f"{OUT}/film_laws_phase2.csv", index=False)
    df["Gu"] = df.G
    prof = {"tau=2.5": p_profile(df, h, "G"), "tau=0": p_profile(df, h, "G_raw")}
    impl = implied_P_exponent(df, df.G.values)
    # entrance factor and quality checks on the chosen law (F2 sat)
    CH = "F3 (k/D) Pe_l^1/3 [film]"
    unit = (df.k_f / D * (df.G / df.rho_l * D * df.rho_f * df.cp_f / df.k_f) ** (1 / 3)).values
    ent = 1 + (D / (df.z + 0.1 * D)) ** 0.7
    o_ent, _ = eval_one_const(df, unit * ent.values, h)
    lowx = df.x.values <= 0.05
    o_lowx, _ = eval_one_const(df[lowx].reset_index(drop=True), unit[lowx], h[lowx])
    # station bias of F2 by case
    C2 = tab[(tab.G == "tau=2.5 s") & (tab.law == CH)].C.values[0]
    r = np.log(C2 * unit / h)
    sb = []
    for c in sorted(df.case.unique()):
        ru = r[(df.case == c) & (df.st == "up")].mean(); rd = r[(df.case == c) & (df.st == "dn")].mean()
        sb.append(dict(case=c, up_minus_dn_resid=ru - rd))
    sb = pd.DataFrame(sb)
    # x-dependence identifiability: x_eq is built from the measured q (x ~ q/G) -> endogenous
    extra0 = dict(corr_x_with_lnq_within_runs=float(np.corrcoef(
        df.x - df.groupby(["case", "st"]).x.transform("mean"),
        np.log(df.q_E) - np.log(df.q_E).groupby([df.case, df.st]).transform("mean"))[0, 1]))
    extra = dict(chosen=CH, **extra0, entrance_factor=dict(C=o_ent["C"], rms=o_ent["rms"], loco=o_ent["loco_rms"], tr_dn=o_ent["tr_dn"]),
                 lowx_only=dict(n=int(lowx.sum()), C=o_lowx["C"], rms=o_lowx["rms"], tr_x=o_lowx["tr_x"]),
                 station_resid_up_minus_dn=sb.set_index("case").up_minus_dn_resid.round(3).to_dict(),
                 implied_P_exponent_at_fixed_G=impl)
    with open(f"{OUT}/film_laws_phase2.json", "w") as f:
        json.dump(dict(profile={k: {kk: {a: b for a, b in vv.items() if a not in ("ps", "dev")} for kk, vv in v.items()}
                                for k, v in prof.items()}, **extra), f, indent=1)
    cols = ["G", "law", "n_fit", "C", "C_loco_min", "C_loco_max", "bias", "rms", "sd_cs_bias", "w30", "loco_rms",
            "tr_lnG", "tr_lnP", "tr_lndT", "tr_dn", "tr_x", "tr_per10s"]
    with open(f"{OUT}/t08_film_laws.md", "w") as f:
        f.write(__doc__ + "\nrms = case-station-weighted RMS of ln(pred/meas); sd_cs_bias = sd over the 20 case-stations "
                "of their mean ln ratio; tr_* = slopes of ln(pred/meas) (multivariate) on ln G, ln P, ln dT, station dn, x_eq; "
                "tr_per10s = within-run time slope.\n\n")
        f.write(md_table(tab[cols].round(3), 4) + "\n\n```\n")
        f.write(json.dumps(dict(profile={k: {kk: {a: b for a, b in vv.items() if a not in ("ps", "dev")} for kk, vv in v.items()}
                                         for k, v in prof.items()}, **extra), indent=1) + "\n```\n")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 42)
    print(tab[cols].round(3).to_string())
    print(json.dumps({k: {kk: {a: b for a, b in vv.items() if a not in ("ps", "dev")} for kk, vv in v.items()} for k, v in prof.items()}, indent=1))
    print(json.dumps(extra, indent=1))
    # figures: parity of chosen law + residual trends + profile
    pred = store[("tau=2.5 s", CH)]
    fig, axs = plt.subplots(1, 4, figsize=(21, 5))
    for i, c in enumerate(sorted(df.case.unique())):
        for st, mk in (("up", "."), ("dn", "x")):
            m = (df.case == c).values & (df.st == st).values
            axs[0].loglog(h[m], pred[m], mk, ms=2, color=COL[i], label=f"{c}" if st == "up" else None)
            axs[1].plot(df.G[m], np.log(pred[m] / h[m]), mk, ms=2, color=COL[i])
            axs[2].plot(df.P[m] / 1e3, np.log(pred[m] / h[m]), mk, ms=2, color=COL[i])
    x = np.array([15, 600]); axs[0].plot(x, x, "k"); axs[0].plot(x, 1.3 * x, "k:"); axs[0].plot(x, x / 1.3, "k:")
    row = tab[(tab.G == "tau=2.5 s") & (tab.law == CH)].iloc[0]
    axs[0].set_title(f"h = {row.C:.2f} (k_v/D) Pe_l^1/3 [film props]; rms(ln) {row.rms:.2f}, LOCO {row.loco_rms:.2f}", fontsize=9)
    axs[0].set_xlabel("h_E measured (wall energy balance) [W/m2K]"); axs[0].set_ylabel("predicted"); axs[0].legend(fontsize=6, ncol=2)
    for a, lab in ((axs[1], "G lag-corr [kg/m2s]"), (axs[2], "P [kPa]")):
        a.axhline(0, color="k"); a.set_xscale("log"); a.set_ylim(-1, 1); a.set_xlabel(lab); a.set_ylabel("ln(pred/meas)")
    for k, v in prof["tau=2.5"].items():
        axs[3].plot(v["ps"], v["dev"], label=f"{k}: p={v['p_hat']:.2f}")
    axs[3].axhline(3.84, color="k", ls=":"); axs[3].set_ylim(0, 20)
    axs[3].set_xlabel("pressure exponent p (G exponent fixed 0.5)"); axs[3].set_ylabel("profile deviance (20 case-stations)")
    axs[3].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig14_film_law_phase2.png", dpi=100); plt.close(fig)


if __name__ == "__main__":
    main()
