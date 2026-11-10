"""Phase 6: parsimony study.  Ablation ladder M0..M5 with leave-one-test-out (LOTO) re-fitting of every
fitted constant; a-priori film-h log error, rewet error (UF and NBS LN2/Cu knees), and 0-D surrogate
wall-T metrics (a10 machinery, fold constants applied to the held-out test only).
Meter lag tau = 2.5 s is held fixed as a rig boundary condition (Phase 4: nesting it changes the final
law's LOTO error 0.240 -> 0.247).

Ladder (film law / rewet law; k = number of constants fitted to UF data):
  M0 k=0  Bromley forced-convection C=2.7 (film props)      / contact with T* = T_hn (Lienhard)
  M1 k=1  C (k_vf/D) Pe_l^(1/3), u_l=G(1-x)/rho_l          / T* = T_hn
  M2 k=2  as M1                                            / T* = T_hn + e Tc   (e fitted)
  M3 k=3  as M1 x (D/L_eff)^(1/3), L_eff = z Lam/(z+Lam)    / T* = T_hn + e Tc   (FINAL)
  M4 k=5  C (k_vf/D) Pe_l^n (P/300 kPa)^m x Lam factor      / T* = T_hn + e Tc
  M5 k=7  Darr 2016 eqs 2-3 with c1..c6 refitted            / T* = T_hn + e Tc
"""
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import D, OUT, md_table
from a08_film_laws import load_fb
from a16_rewet_family import rule
from a09_rewet_theory import cu_e, cu_cp


# ---------------------------------------------------------------- film laws (work on DataFrame or surrogate state)
def _pe(d, G, n=1 / 3):
    ul = G * (1 - np.clip(d.x, 0, 0.9)) / d.rho_l
    return d.k_vf / D * (ul * D * d.rho_vf * d.cp_vf / d.k_vf) ** n


def _lam(d, lam):
    return (D * (d.z + lam) / (d.z * lam)) ** (1 / 3)


def bromley(d, G, p):
    ul = G / d.rho_l
    return 2.7 * np.sqrt(d.k_vf * d.rho_vf * (d.h_fg + 0.5 * d.cp_vf * d.dT) * ul / (D * d.dT))


def m1(d, G, p):
    return p[0] * _pe(d, G)


def m3(d, G, p):
    return p[0] * _pe(d, G) * _lam(d, np.exp(p[1]))


def m4(d, G, p):
    return p[0] * _pe(d, G, p[2]) * (d.P / 3e5) ** p[3] * _lam(d, np.exp(p[1]))


def darr(d, G, p):
    c1, c2, c3, c4, c5, c6 = p
    L = d.z
    x = np.clip(d.x, 0, 0.99)
    htc = c1 * (1 + (D / (L + 0.1 * D)) ** 0.7) * (1 - x) ** c2 * (G * D / d.mu_v) ** c3 * (d.cp_v * d.mu_v / d.k_v) ** 0.4 * d.k_v / D
    We = G**2 * D / (d.rho_l * d.sigma)
    Ja = d.cp_l * d.dT / d.h_fg
    return htc + abs(c4) * 1e-6 * d.k_l / D * np.exp(-((L / D) ** abs(c5))) * We**c6 / Ja


FILM = {"M0": (bromley, []), "M1": (m1, [15.0]), "M2": (m1, [15.0]), "M3": (m3, [15.7, np.log(0.1)]),
        "M4": (m4, [15.7, np.log(0.1), 1 / 3, 0.0]), "M5": (darr, [0.15, -0.76, 0.52, 5.0, 0.1, 2.0])}
REWET_E = {"M0": False, "M1": False, "M2": True, "M3": True, "M4": True, "M5": True}
KFIT = {"M0": 0, "M1": 1, "M2": 2, "M3": 3, "M4": 5, "M5": 7}


def to_ns(df):
    d = df.rename(columns={"k_f": "k_vf", "rho_f": "rho_vf", "cp_f": "cp_vf", "mu_f": "mu_vf"})
    return d


def fit_film(f, p0, d, h, w):
    if not p0:
        return np.array([])
    sw = np.sqrt(w)
    return least_squares(lambda p: sw * (np.log(np.clip(f(d, d.G.values, p), 1e-3, None)) - np.log(h)), p0,
                         max_nfev=400).x


def film_ladder(df):
    d = to_ns(df)
    h, w = d.h_E.values, d.w.values
    res = {}
    for M, (f, p0) in FILM.items():
        p = fit_film(f, p0, d, h, w)
        rin = np.log(f(d, d.G.values, p) / h)
        rcv = np.empty(len(d)); folds = {}
        for k in range(1, 11):
            tr = d.case.values != k
            pk = fit_film(f, list(p) if len(p) else [], d[tr], h[tr], w[tr]) if len(p) else p
            folds[k] = pk
            rcv[~tr] = np.log(f(d[~tr], d.G.values[~tr], pk) / h[~tr])
        cs_in = pd.Series(rin).groupby([d.case.values, d.st.values]).mean().values
        cs_cv = pd.Series(rcv).groupby([d.case.values, d.st.values]).mean().values
        # lag-1 autocorrelation of within-run residual (for n_eff)
        rho = []
        for (c, s), g in pd.Series(rin).groupby([d.case.values, d.st.values]):
            v = g.values - g.values.mean()
            rho.append(np.sum(v[1:] * v[:-1]) / np.sum(v * v))
        res[M] = dict(p=p, folds=folds, rms_in=float(np.sqrt(np.sum(w * rin**2) / w.sum())),
                      loto=float(np.sqrt(np.sum(w * rcv**2) / w.sum())), cs_in=cs_in, cs_cv=cs_cv,
                      rho1=float(np.median(rho)), N=len(d))
    return res


def rewet_ladder():
    rt = pd.read_csv(f"{OUT}/rewet_theory_table.csv")
    P, T = rt.P_kPa.values * 1e3, rt.T_mfb.values
    nb = pd.read_csv(f"{OUT}/nbs_knees.csv")
    nb = nb[(nb.fluid == "LN2") & (nb.ratio >= 1.9) & (nb.dT_seg_at_knee <= 20) & (nb.z_m < 60)]
    from scipy.optimize import minimize_scalar
    fit = lambda m: minimize_scalar(lambda e: np.mean((rule("R5", e, P[m]) [0] - T[m]) ** 2), bounds=(-0.2, 0.5), method="bounded").x
    out = {}
    for lab, fitted in (("T_hn contact", False), ("T_hn + e Tc contact", True)):
        e_all = fit(np.ones(len(T), bool)) if fitted else 0.0
        cv = np.empty(len(T)); ef = {}
        for c in rt.case.unique():
            m = rt.case.values == c
            ef[int(c)] = fit(~m) if fitted else 0.0
            cv[m] = rule("R5", ef[int(c)], P[m])[0]
        pin = rule("R5", e_all, P)[0]
        pn = np.array([rule("R5", e_all, p * 1e3, "Nitrogen", cu_e, cu_cp)[0][0] for p in nb.P_local_kPa])
        out[lab] = dict(e=e_all, folds=ef, mae_in=np.abs(pin - T).mean(), loto_mae=np.abs(cv - T).mean(),
                        loto_rmse=np.sqrt(np.mean((cv - T) ** 2)), nbs_mae=np.abs(pn - nb.T_knee.values).mean(),
                        nbs_bias=(pn - nb.T_knee.values).mean(), n_nbs=len(nb))
    return out


def surrogate(film, rew):
    """LOTO 0-D surrogate: fold constants of test k applied to test k's two stations."""
    import a10_boiling_curve_sensitivity as a10
    recs = a10.inputs()
    a10.DT_STEP = 0.01
    rows = {}
    for M in FILM:
        f = FILM[M][0]
        mets = []
        for k in range(1, 11):
            pk = film[M]["folds"][k]
            ek = rew["T_hn + e Tc contact"]["folds"].get(k, rew["T_hn + e Tc contact"]["e"]) if REWET_E[M] else 0.0
            rr = [r for r in recs if r["case"] == k]
            Pg = np.linspace(min(r["P"].min() for r in rr) * 0.98, max(r["P"].max() for r in rr) * 1.02, 12)
            Tg = rule("R5", ek, Pg)[0]
            rr = [dict(r, Tmfb=np.interp(r["P"], Pg, Tg)) for r in rr]
            a10.h_film = (lambda d, G, f=f, pk=pk: f(d, G, pk))
            tt, sim = a10.simulate(rr, dict(a10.BASE))
            mets.append(a10.metrics(rr, tt, sim))
        m = pd.concat(mets, ignore_index=True)
        rows[M] = dict(rms_mean=m.rms.mean(), rms_med=m.rms.median(), t130_mae=m.dt130.abs().mean(),
                       t100_mae=m.dt100.abs().mean(), miss130=int(m.miss130.sum()), miss100=int(m.miss100.sum()),
                       t100_n=int(m.dt100.notna().sum()))
        print(M, rows[M], flush=True)
    return rows


def info(film):
    """AIC/BIC with two effective sample sizes: (a) n = 20 test-stations, residual = station-mean log error
    (within-run residuals are strongly autocorrelated and the dominant error is the run-level bias);
    (b) AR(1) n_eff = N (1 - rho1)/(1 + rho1) on the pointwise weighted RMS."""
    out = {}
    for M, r in film.items():
        k = len(r["p"])
        n = len(r["cs_in"]); s2 = np.mean(r["cs_in"] ** 2)
        aic = n * np.log(s2) + 2 * k; bic = n * np.log(s2) + k * np.log(n)
        ne = r["N"] * (1 - r["rho1"]) / (1 + r["rho1"])
        s2p = r["rms_in"] ** 2
        out[M] = dict(k_film=k, n_station=n, AIC_station=aic, BIC_station=bic, rho1=r["rho1"], n_eff_ar1=ne,
                      AIC_ar1=ne * np.log(s2p) + 2 * k, BIC_ar1=ne * np.log(s2p) + k * np.log(ne))
    return out


def main(do_sur=True):
    t0 = time.time()
    df = load_fb()
    film = film_ladder(df)
    rew = rewet_ladder()
    ic = info(film)
    print("film/rewet done", time.time() - t0, flush=True)
    sur = surrogate(film, rew) if do_sur else {}
    rows = []
    for M in FILM:
        rw = rew["T_hn + e Tc contact"] if REWET_E[M] else rew["T_hn contact"]
        r = dict(model=M, k_fit=KFIT[M], film_consts=",".join(f"{v:.4g}" for v in film[M]["p"]),
                 h_rms_in=film[M]["rms_in"], h_loto=film[M]["loto"],
                 station_bias_sd_loto=float(np.std(film[M]["cs_cv"], ddof=1)),
                 Tmfb_e=rw["e"], Tmfb_loto_mae=rw["loto_mae"], Tmfb_loto_rmse=rw["loto_rmse"], nbs_knee_mae=rw["nbs_mae"],
                 **{k: v for k, v in ic[M].items() if k != "k_film"}, **sur.get(M, {}))
        rows.append(r)
    tab = pd.DataFrame(rows)
    tab.to_csv(f"{OUT}/parsimony.csv", index=False)
    json.dump({M: {"p": film[M]["p"].tolist(), "folds": {k: np.asarray(v).tolist() for k, v in film[M]["folds"].items()}}
               for M in FILM} | {"rewet": {k: {"e": v["e"], "folds": v["folds"]} for k, v in rew.items()}},
              open(f"{OUT}/parsimony_constants.json", "w"), indent=1)
    pd.set_option("display.width", 250)
    print(tab.round(3).to_string())


if __name__ == "__main__":
    main(do_sur="nosur" not in sys.argv)
