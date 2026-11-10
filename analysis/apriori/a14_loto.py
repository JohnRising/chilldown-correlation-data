"""Phase 4 / Tasks 1-3 (referee M3, M8, M9): nested leave-one-test-out (LOTO) with the meter lag tau
re-fitted in each fold; the same on the PUBLISHED flux; paired test-bootstrap of LOTO differences.

Final law (L3x):  h = C (k_vf/D) Pe_l^(1/3) (D/L_eff)^(1/3),  L_eff = z Lam/(z+Lam),  u_l = G(1-x_e)/rho_l
Inner fit per fold: grid over tau in {0, 0.25, ..., 5} s and Lam (log grid 0.01-3 m); C closed-form
(weighted mean log ratio).  x_e is recomputed for each tau from the station fluxes and G(t+tau).
Equal weight per test-station.  Errors: RMS of ln(pred/meas)."""
import json

import numpy as np
import pandas as pd

from common import D, OUT, Z, load, md_table
from a08_film_laws import load_fb, darr

TAUS = np.arange(0.0, 5.01, 0.25)
LAMS = np.exp(np.linspace(np.log(0.01), np.log(3.0), 61))


def prep(flux="E"):
    df = load_fb()
    q = pd.read_csv(f"{OUT}/q_consistent.csv")
    col = "q_E" if flux == "E" else "q_pub"
    qq = q.pivot_table(index=["case", "t"], columns="st", values=col).reset_index()
    qq.columns = ["case", "t", "q_dn", "q_up"]
    df = df.merge(qq, on=["case", "t"], how="left")
    df[["q_up", "q_dn"]] = df[["q_up", "q_dn"]].ffill().clip(lower=0)
    if flux == "pub":
        df["h_E"] = df.h_pub            # target replaced by published h (q_pub/(Ti_pub - Tsat))
        df = df[df.h_E > 0].reset_index(drop=True)
    I = np.where(df.st == "up", df.q_up * Z["up"], df.q_up * Z["up"] + 0.5 * (df.q_up + df.q_dn) * (Z["dn"] - Z["up"]))
    df["I"] = I
    meter = {c: (load(c).t_s.values, load(c).G_kgm2s.values) for c in range(1, 11)}
    G = {}
    for tau in TAUS:
        g = np.empty(len(df))
        for c in range(1, 11):
            m = df.case.values == c
            t, Gm = meter[c]
            g[m] = np.interp(df.t.values[m] + tau, t, Gm)
        G[round(tau, 2)] = np.clip(g, 5, None)
    return df, G


def units(df, G):
    """log-units (without C) of each law for one G array; L3x without the Lam factor."""
    x = np.clip(4 * df.I.values / (G * D * df.h_fg.values), 0, 0.9)
    ul = G * (1 - x) / df.rho_l.values
    pe = lambda u: (df.k_f / D * (u * D * df.rho_f * df.cp_f / df.k_f) ** (1 / 3)).values
    d2 = df.copy(); d2["x"] = x
    return dict(L3x=np.log(pe(ul)), L0=np.log(pe(G / df.rho_l.values)),
                Darr=np.log(darr(d2, G).values), F2=np.log(np.sqrt(df.k_v * df.cp_v * G / D).values))


def lamfac(df, lam):
    return (1 / 3) * np.log(D * (df.z.values + lam) / (df.z.values * lam))


class Fitter:
    def __init__(self, df, G):
        self.df, self.y, self.w = df, np.log(df.h_E.values), df.w.values
        self.U = {tau: units(df, g) for tau, g in G.items()}
        self.fac = {lam: lamfac(df, lam) for lam in LAMS}

    def fit(self, law, m, tau=None):
        """returns (tau, lam, lnC, rss) minimising weighted RSS on mask m."""
        best = (None, None, None, np.inf)
        w = self.w[m]
        for t in (self.U if tau is None else [tau]):
            base = self.y[m] - self.U[t][law][m]
            for lam in (LAMS if law == "L3x" else [None]):
                r = base - (self.fac[lam][m] if lam is not None else 0)
                lc = np.sum(w * r) / w.sum()
                rss = np.sum(w * (r - lc) ** 2)
                if rss < best[3]:
                    best = (t, lam, lc, rss)
        return best

    def pred(self, law, t, lam, lc):
        return lc + self.U[t][law] + (self.fac[lam] if lam is not None else 0)


def loto(F, law, tau=None):
    df = F.df
    r = np.empty(len(df)); consts = {}
    for k in range(1, 11):
        tr = df.case.values != k
        t, lam, lc, _ = F.fit(law, tr, tau)
        r[~tr] = (F.pred(law, t, lam, lc) - F.y)[~tr]
        consts[k] = dict(tau=t, C=float(np.exp(lc)), Lambda=lam)
    t, lam, lc, _ = F.fit(law, np.ones(len(df), bool), tau)
    rin = F.pred(law, t, lam, lc) - F.y
    per = pd.DataFrame(dict(case=df.case.values, w=F.w, r2=F.w * r**2)).groupby("case").sum()
    return dict(all=dict(tau=t, C=float(np.exp(lc)), Lambda=lam), folds=consts,
                rms_in=float(np.sqrt(np.sum(F.w * rin**2) / F.w.sum())),
                loto=float(np.sqrt(np.sum(F.w * r**2) / F.w.sum())), per=per)


def boot_diff(res, pairs, nb=4000, seed=0):
    rng = np.random.default_rng(seed)
    cases = np.arange(1, 11)
    out = []
    idx = rng.integers(0, 10, (nb, 10))
    for a, b in pairs:
        pa, pb = res[a]["per"], res[b]["per"]
        ra = np.sqrt(pa.r2.values[idx].sum(1) / pa.w.values[idx].sum(1))
        rb = np.sqrt(pb.r2.values[idx].sum(1) / pb.w.values[idx].sum(1))
        d = rb - ra
        out.append(dict(comparison=f"{b} - {a}", diff=res[b]["loto"] - res[a]["loto"], se=d.std(),
                        ci_lo=np.percentile(d, 2.5), ci_hi=np.percentile(d, 97.5), p_b_better=float(np.mean(d < 0))))
    return pd.DataFrame(out)


def free_law(df, G, nb=400):
    w = df.w.values; sw = np.sqrt(w)
    X = np.c_[np.ones(len(df)), np.log(G), np.log(df.P), (df.st == "dn")]
    y = np.log(df.h_E.values)
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    rng = np.random.default_rng(2); bs = []
    for _ in range(nb):
        pick = rng.choice(np.arange(1, 11), 10)
        ii = np.concatenate([np.where(df.case.values == c)[0] for c in pick])
        bs.append(np.linalg.lstsq(X[ii] * sw[ii, None], y[ii] * sw[ii], rcond=None)[0])
    bs = np.array(bs)
    return dict(bG=b[1], bG_ci=[np.percentile(bs[:, 1], 2.5), np.percentile(bs[:, 1], 97.5)],
                bP=b[2], bP_ci=[np.percentile(bs[:, 2], 2.5), np.percentile(bs[:, 2], 97.5)])


def main():
    out, tables = {}, []
    for flux in ("E", "pub"):
        df, G = prep(flux)
        F = Fitter(df, G)
        res = {}
        for law in ("L3x", "L0", "Darr", "F2"):
            res[law] = loto(F, law)                  # nested tau
            res[law + " tau=2.5"] = loto(F, law, 2.5)
            res[law + " tau=0"] = loto(F, law, 0.0)
        for k, v in res.items():
            sp = pd.DataFrame(v["folds"]).T
            tables.append(dict(flux=flux, law=k, tau_all=v["all"]["tau"], C_all=v["all"]["C"], Lam_all=v["all"]["Lambda"],
                               tau_folds=f"{sp.tau.min():.2f}-{sp.tau.max():.2f}", C_folds=f"{sp.C.min():.3g}-{sp.C.max():.3g}",
                               Lam_folds=(f"{sp.Lambda.min():.3g}-{sp.Lambda.max():.3g}" if sp.Lambda.notna().all() else "-"),
                               rms_in=v["rms_in"], loto=v["loto"]))
        fl = {f"tau={t}": free_law(df, G[t]) for t in (0.0, 2.5)}
        bd = boot_diff(res, [("L3x", "L0"), ("L3x", "Darr"), ("L3x", "F2"), ("L0", "Darr"), ("L0", "F2"), ("Darr", "F2")])
        out[flux] = dict(res=res, free=fl, boot=bd, df=df)
    tab = pd.DataFrame(tables)
    # constants file (consistent flux, final law, nested tau)
    r = out["E"]["res"]["L3x"]
    js = {str(k): {kk: (None if vv is None else float(vv)) for kk, vv in v.items()} for k, v in r["folds"].items()}
    js["all"] = {k: float(v) for k, v in r["all"].items()}
    js["note"] = "nested LOTO: tau, C, Lambda re-fitted on the 9 training tests (consistent flux q_E); Lambda in m, tau in s"
    json.dump(js, open(f"{OUT}/loto_constants.json", "w"), indent=1)
    # matched-G comparison test 4 vs 3/5/6
    q = pd.read_csv(f"{OUT}/q_consistent.csv")
    q = q[q.fb & (q.t > 1)]
    mg = []
    for c in (3, 4, 5, 6):
        s = q[q.case == c]
        mg.append(dict(case=c, corrected=c != 4, G=s.G.median(), P_kPa=s.P.median() / 1e3, dT=s.dT.median(),
                       h_E=s.h_E.median(), h_pub=s.h_pub.median()))
    mg = pd.DataFrame(mg)
    h4E, h4p = mg[mg.case == 4].h_E.values[0], mg[mg.case == 4].h_pub.values[0]
    mg["h_E/h_E(4)"] = mg.h_E / h4E; mg["h_pub/h_pub(4)"] = mg.h_pub / h4p
    ho = open(f"{OUT}/t07_consistent_q.md").read().split("## hold-out arbiter")[1]
    pd.set_option("display.width", 250)
    print(tab.round(3).to_string())
    for f in ("E", "pub"):
        print(f, json.dumps(out[f]["free"], default=float)); print(out[f]["boot"].round(3).to_string())
    print(mg.round(3).to_string())
    print(json.dumps(js, indent=0)[:900])
    with open(f"{OUT}/t14_loto.md", "w") as fh:
        fh.write(__doc__ + "\n## constants and errors (E = energy-consistent q, pub = published q)\n\n" + md_table(tab.round(3), 4))
        for f in ("E", "pub"):
            fh.write(f"\n\n## free law ln h = a + b ln G + c ln P + d dn ({f}; test-bootstrap 95 % CI)\n\n```\n"
                     + json.dumps(out[f]["free"], default=float, indent=1) + "\n```\n\n## paired test-bootstrap of LOTO differences "
                     f"(nested tau; {f})\n\n" + md_table(out[f]["boot"].round(3), 3))
        fh.write("\n\n## matched-G/P comparison, test 4 (uncorrected) vs 3/5/6 (corrected)\n\n" + md_table(mg.round(3), 3))
        fh.write("\n\n## hold-out arbiter (from t07)" + ho)


if __name__ == "__main__":
    main()
