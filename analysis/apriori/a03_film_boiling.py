"""Goal 3: film-boiling scaling and a-priori closure evaluation."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import D, FIG, OUT, err_stats, md_table
import fb_models as fm

COL = plt.cm.viridis(np.linspace(0, 0.95, 10))


def weights(df):
    n = df.groupby(["case", "st"]).t.transform("size")
    return 1.0 / n.values


# ------------------------------------------------------------- candidate registry
def candidates():
    """name -> (n_fit, p0, f(df, G, p) -> h)."""
    C = {}
    C["Bromley C=2.7 (film props)"] = (0, [], lambda d, G, p: 2.7 * fm.bromley_unit(d, G, "film"))
    C["Bromley C=2.7 (sat props)"] = (0, [], lambda d, G, p: 2.7 * fm.bromley_unit(d, G, "sat"))
    C["Darr 2016 low-g (x_eq)"] = (0, [], lambda d, G, p: fm.darr2016(d, G))
    C["Darr 2016 low-g (x=0)"] = (0, [], lambda d, G, p: fm.darr2016(d, G, x=0 * d.x))
    C["Dittus-Boelter vapour, Re=GD/mu_v (sat)"] = (0, [], lambda d, G, p: fm.dittus_boelter_v(d, G, "sat"))
    C["Dittus-Boelter vapour (film props)"] = (0, [], lambda d, G, p: fm.dittus_boelter_v(d, G, "film"))
    C["Miropolskii (x_eq)"] = (0, [], lambda d, G, p: fm.miropolskii(d, G, d.x))
    C["Miropolskii (x=0)"] = (0, [], lambda d, G, p: fm.miropolskii(d, G, 0 * d.x))
    C["Groeneveld 1973 tube (x_eq)"] = (0, [], lambda d, G, p: fm.groeneveld(d, G, d.x))
    C["Couette film, ell=z (from inlet)"] = (0, [], lambda d, G, p: fm.couette(d, G, d.z))
    C["Couette film, ell=D"] = (0, [], lambda d, G, p: fm.couette(d, G, D))
    C["Couette film, ell=lambda_KH(u_l)"] = (0, [], lambda d, G, p: fm.couette(d, G, fm.lambda_kh(d, G)))
    C["Penetration (sensible), ell=D, C=1"] = (0, [], lambda d, G, p: fm.penetration_unit(d, G))
    # one fitted constant
    C["Bromley, C fit (film props)"] = (1, [2.7], lambda d, G, p: p[0] * fm.bromley_unit(d, G, "film"))
    C["Couette, ell=C*lambda_KH fit"] = (1, [1.0], lambda d, G, p: fm.couette(d, G, p[0] * fm.lambda_kh(d, G)))
    C["Penetration, C fit"] = (1, [1.0], lambda d, G, p: p[0] * fm.penetration_unit(d, G))
    C["Dittus-Boelter vapour (sat), C fit"] = (1, [0.5], lambda d, G, p: p[0] * fm.dittus_boelter_v(d, G, "sat"))
    C["Darr 2016 x C fit"] = (1, [1.0], lambda d, G, p: p[0] * fm.darr2016(d, G))
    # two fitted constants
    C["Bromley, C & beta fit (hfg'=hfg+beta cp dT)"] = (
        2, [2.7, 0.5], lambda d, G, p: p[0] * fm.bromley_unit(d, G, "film", beta=abs(p[1])))
    C["Pe_l power law: C (k_v/D) Pe_l^n"] = (
        2, [0.5, 0.5], lambda d, G, p: p[0] * d.k_vf / D * (G / d.rho_l * D * d.rho_vf * d.cp_vf / d.k_vf) ** p[1])
    C["Pe_l power law x Darr entrance factor"] = (
        2, [0.5, 0.5], lambda d, G, p: p[0] * (1 + (D / (d.z + 0.1 * D)) ** 0.7) * d.k_vf / D
        * (G / d.rho_l * D * d.rho_vf * d.cp_vf / d.k_vf) ** p[1])
    C["Re_v power law: C Re_v^n Pr^0.4 k_v/D (sat)"] = (
        2, [0.02, 0.8], lambda d, G, p: p[0] * (G * D / d.mu_v) ** p[1] * (d.cp_v * d.mu_v / d.k_v) ** 0.4 * d.k_v / D)
    return C


def fit(df, G, h, f, p0, w):
    if len(p0) == 0:
        return np.array([])
    sw = np.sqrt(w)

    def res(p):
        return sw * (np.log(np.clip(f(df, G, p), 1e-9, None)) - np.log(h))
    r = least_squares(res, p0, x_scale="jac")
    return r.x


def residual_trends(df, r):
    """(i) multivariate regression of the log-residual on ln G, ln P, ln dT and station;
    (ii) within-run time trend: slope of the case/station-demeaned residual vs demeaned t (per 10 s)."""
    X = np.c_[np.ones(len(df)), np.log(df.G / 150), np.log(df.P / 3e5), np.log(df.dT / 100),
              (df.st == "dn").astype(float)]
    b = np.linalg.lstsq(X, r, rcond=None)[0]
    out = dict(zip(["r0", "dlnG", "dlnP", "dlndT", "dn"], b))
    key = df.case.astype(str) + df.st
    rd = r - pd.Series(r).groupby(key.values).transform("mean").values
    td = df.t.values - df.groupby(key.values).t.transform("mean").values
    out["per10s"] = float(np.sum(rd * td) / np.sum(td * td) * 10)
    return out


def evaluate(df, Gcol, hcol, label, loco=True):
    G = df[Gcol].values
    h = df[hcol].values
    w = weights(df)
    rows, preds, params = [], {}, {}
    for name, (nf, p0, f) in candidates().items():
        p = fit(df, G, h, f, p0, w)
        hp = f(df, G, p)
        st = err_stats(hp, h)
        r = np.log(hp / h)
        st["rms_log_w"] = float(np.sqrt(np.sum(w * r**2) / np.sum(w)))
        # case-station mean of per-cs RMS
        row = dict(closure=name, n_fit=nf, params=",".join(f"{v:.4g}" for v in p), **st)
        if nf and loco:
            rcv = np.empty(len(df))
            pl = []
            for c in sorted(df.case.unique()):
                tr = df.case.values != c
                pc = fit(df[tr], G[tr], h[tr], f, p, w[tr])
                pl.append(pc)
                rcv[~tr] = np.log(f(df[~tr], G[~tr], pc) / h[~tr])
            row["loco_rms_log"] = float(np.sqrt(np.mean(rcv**2)))
            row["loco_bias_casemax"] = float(max(abs(np.exp(rcv[df.case.values == c].mean()) - 1)
                                                 for c in df.case.unique()))
            row["param_range_loco"] = ";".join(f"{np.min(a):.3g}-{np.max(a):.3g}" for a in np.array(pl).T)
        row.update({k: v for k, v in residual_trends(df, r).items() if k != "r0"})
        rows.append(row)
        preds[name] = hp
        params[name] = p.tolist()
    out = pd.DataFrame(rows)
    out.insert(0, "data", label)
    return out, preds, params


def free_regression(df, Gcol, hcol):
    """diagnostic power law: ln h = a + b ln G + c ln dT + d ln P + e*dn (equal case-station weights)
    and the dimensionless analogue: ln(Nu_vf) = a + b ln Re_l + c ln(dT/Tsat) + d ln(rho_v/rho_l) + e dn."""
    w = weights(df)
    sw = np.sqrt(w)
    out = {}
    X = np.c_[np.ones(len(df)), np.log(df[Gcol]), np.log(df.dT), np.log(df.P / 1e5), (df.st == "dn")]
    y = np.log(df[hcol])
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    r = y - X @ b
    # LOCO
    rcv = np.empty(len(df))
    for c in df.case.unique():
        tr = df.case.values != c
        bc = np.linalg.lstsq(X[tr] * sw[tr, None], y[tr] * sw[tr], rcond=None)[0]
        rcv[~tr] = y[~tr] - X[~tr] @ bc
    # bootstrap by case for CI
    rng = np.random.default_rng(0)
    cs = df.case.unique()
    bs = []
    for _ in range(400):
        pick = rng.choice(cs, len(cs))
        idx = np.concatenate([np.where(df.case.values == c)[0] for c in pick])
        bs.append(np.linalg.lstsq(X[idx] * sw[idx, None], y.values[idx] * sw[idx], rcond=None)[0])
    bs = np.array(bs)
    out["dimensional"] = dict(names=["a", "b_lnG", "c_lndT", "d_lnP", "e_dn"], coef=b.tolist(),
                              ci95_lo=np.percentile(bs, 2.5, 0).tolist(), ci95_hi=np.percentile(bs, 97.5, 0).tolist(),
                              rms=float(np.sqrt(np.mean(r**2))), loco_rms=float(np.sqrt(np.mean(rcv**2))))
    # G + P only (no dT)
    X2 = X[:, [0, 1, 3, 4]]
    b2 = np.linalg.lstsq(X2 * sw[:, None], y * sw, rcond=None)[0]
    out["no_dT"] = dict(coef=b2.tolist(), rms=float(np.sqrt(np.mean((y - X2 @ b2) ** 2))))
    # dimensionless
    Nu = df[hcol] * D / df.k_vf
    Re_l = df[Gcol] * D / df.mu_l
    Xd = np.c_[np.ones(len(df)), np.log(Re_l), np.log(df.dT / df.Tsat), np.log(df.rho_v / df.rho_l), (df.st == "dn")]
    bd = np.linalg.lstsq(Xd * sw[:, None], np.log(Nu) * sw, rcond=None)[0]
    out["dimensionless"] = dict(names=["a", "b_lnRe_l", "c_ln(dT/Tsat)", "d_ln(rho_v/rho_l)", "e_dn"],
                                coef=bd.tolist(), rms=float(np.sqrt(np.mean((np.log(Nu) - Xd @ bd) ** 2))))
    return out


def parity_plot(df, preds, h, names, fn, title):
    n = len(names)
    fig, axs = plt.subplots(int(np.ceil(n / 4)), 4, figsize=(18, 4.4 * np.ceil(n / 4)))
    axs = np.atleast_1d(axs).ravel()
    for a, nm in zip(axs, names):
        for i, c in enumerate(sorted(df.case.unique())):
            m = df.case.values == c
            a.loglog(h[m], preds[nm][m], ".", ms=2, color=COL[i], label=f"{c}")
        lo, hi = 15, 800
        x = np.array([lo, hi])
        a.plot(x, x, "k-", lw=0.8); a.plot(x, 1.3 * x, "k:", lw=0.6); a.plot(x, x / 1.3, "k:", lw=0.6)
        a.set_xlim(lo, hi); a.set_ylim(lo / 3, hi * 3)
        st = err_stats(preds[nm], h)
        a.set_title(f"{nm}\nbias {st['bias']*100:+.0f}%  rms(ln) {st['rms_log']:.2f}  ±30%: {st['w30']*100:.0f}%", fontsize=8)
        a.set_xlabel("h measured [W/m2K]"); a.set_ylabel("h predicted")
    axs[0].legend(fontsize=6, ncol=2, title="case", title_fontsize=6)
    for a in axs[n:]:
        a.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(fn, dpi=95)
    plt.close(fig)


def residual_plot(df, preds, h, names, fn):
    fig, axs = plt.subplots(len(names), 4, figsize=(17, 3.3 * len(names)))
    for i, nm in enumerate(names):
        r = np.log(preds[nm] / h)
        for j, (col, lab, lg) in enumerate([("G", "G (lag-corr.) [kg/m2s]", True), ("P", "P [Pa]", True),
                                             ("dT", "Ti-Tsat [K]", True), ("t", "t [s]", False)]):
            a = axs[i, j]
            for k, c in enumerate(sorted(df.case.unique())):
                for st, mk in (("up", "."), ("dn", "x")):
                    m = (df.case.values == c) & (df.st.values == st)
                    a.plot(df[col].values[m], r[m], mk, ms=2, color=COL[k])
            if lg:
                a.set_xscale("log")
            a.axhline(0, color="k", lw=0.8)
            a.set_ylim(-1.2, 1.2)
            a.set_xlabel(lab)
            if j == 0:
                a.set_ylabel(f"ln(pred/meas)\n{nm}", fontsize=8)
    fig.tight_layout()
    fig.savefig(fn, dpi=90)
    plt.close(fig)


def main():
    df = fm.dataset()
    df.to_csv(f"{OUT}/film_boiling_dataset.csv", index=False)
    print("film-boiling samples:", len(df), df.groupby(["case", "st"]).size().to_dict())
    allres = []
    store = {}
    for Gcol, hcol, lab in [("G", "h_ec", "h_ec, G lag-corr (PRIMARY)"), ("G_raw", "h_ec", "h_ec, G raw"),
                            ("G", "h_raw", "h_published, G lag-corr"), ("G_raw", "h_raw", "h_published, G raw")]:
        res, preds, params = evaluate(df, Gcol, hcol, lab, loco=True)
        allres.append(res)
        store[lab] = (preds, params, hcol)
    tab = pd.concat(allres, ignore_index=True)
    tab.to_csv(f"{OUT}/film_boiling_closures.csv", index=False)
    freg = {lab: free_regression(df, G, h) for G, h, lab in
            [("G", "h_ec", "h_ec, G lag-corr"), ("G_raw", "h_ec", "h_ec, G raw"), ("G", "h_raw", "h_published, G lag-corr")]}
    with open(f"{OUT}/film_boiling_free_regression.json", "w") as f:
        json.dump(freg, f, indent=1)
    # params for the surrogate
    with open(f"{OUT}/film_boiling_params.json", "w") as f:
        json.dump(store["h_ec, G lag-corr (PRIMARY)"][1], f, indent=1)
    cols = ["closure", "n_fit", "params", "n", "bias", "rms_log", "rms_log_w", "w20", "w30", "loco_rms_log",
            "loco_bias_casemax", "param_range_loco", "dlnG", "dlnP", "dlndT", "dn", "per10s"]
    with open(f"{OUT}/t03_film_boiling.md", "w") as f:
        for lab in tab.data.unique():
            f.write(f"### {lab}\n\n")
            f.write(md_table(tab[tab.data == lab][cols].round(3), 3) + "\n\n")
        f.write("### free power-law regressions\n\n```\n" + json.dumps(freg, indent=1) + "\n```\n")
    pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 45)
    for lab in tab.data.unique():
        print("\n==", lab)
        print(tab[tab.data == lab][["closure", "n_fit", "params", "bias", "rms_log", "w20", "w30", "loco_rms_log",
                                    "dlnG", "dlnP", "dlndT", "dn", "per10s"]].round(3).to_string())
    print(json.dumps(freg, indent=1))
    # irreducible scatter: residual of ln h about a per-case/station quadratic in t
    fl = []
    for (c, st), g in df.groupby(["case", "st"]):
        y = np.log(g.h_ec.values); t = g.t.values
        fl.append(y - np.polyval(np.polyfit(t, y, 2), t))
    floor = float(np.sqrt(np.mean(np.concatenate(fl) ** 2)))
    print("noise floor rms(ln h) about per-run quadratic trend:", floor)
    with open(f"{OUT}/t03_film_boiling.md", "a") as f:
        f.write(f"\nIrreducible scatter (rms of ln h about a per-case/station quadratic in t): {floor:.3f}\n")
    preds, params, hcol = store["h_ec, G lag-corr (PRIMARY)"]
    h = df[hcol].values
    names = list(candidates().keys())
    parity_plot(df, preds, h, names, f"{FIG}/fig08_fb_parity.png",
                "Film boiling a-priori parity (energy-consistent h, lag-corrected G, tau=2.5 s)")
    best = ["Bromley, C fit (film props)", "Penetration, C fit", "Pe_l power law: C (k_v/D) Pe_l^n",
            "Darr 2016 low-g (x_eq)", "Dittus-Boelter vapour (sat), C fit"]
    residual_plot(df, preds, h, best, f"{FIG}/fig09_fb_residuals.png")


if __name__ == "__main__":
    main()
