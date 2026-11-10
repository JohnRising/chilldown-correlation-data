"""Phase 2 / Task 4: full gravity-free boiling curve + 0-D surrogate sensitivity (tornado).

Wall: m''_w cp(T) dT/dt = -q''(T) + q''_par(T),  q''_par = 98.6 + 2.2 (293 - T) W/m2 (a07),
compared with the station-mean TC history T_m (what a 1-D model is validated on).  Inputs at each
station: measured P(t), Tsat(t), lag-corrected G(t) (tau = 2.5 s).
Baseline boiling curve (all regimes as functions of local node variables):
  film       : q_FB = h_FB (T - Tsat), h_FB = C (k_vf/D) Pe_l^(1/3)  (a08, C = 8.55)
  rewet      : T_MFB = T_c + beta (T_c - Tsat), beta = e_l/e_w(T_MFB)       (a09, 0 constants)
  nucleate   : q_NB = h_DB,l dT + (A dT)^3,  A = 55 pr^0.12 (-log10 pr)^-0.55 M^-0.5  (Cooper 1984
               solved for q; Dittus-Boelter liquid floor h_DB,l = 0.023 Re_l^0.8 Pr_l^0.4 k_l/D, Re_l = G D/mu_l)
  CHF        : q_CHF = Katto-Ohno (1984) saturated, L = z; dT_CHF = root of q_NB(dT) = q_CHF
  transition : q_TB = q_FB(T_MFB) + (q_CHF - q_FB(T_MFB)) theta^n, theta = (T_MFB-T)/(T_MFB-T_CHF), n = 0.5
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CASES, D, FIG, MW_A, OUT, PC_N2, STATIONS, TC_N2, Z, md_table, sat, ss_cp, ss_e
import a06_surrogate as s6
from a05_chf_transition import katto_ohno

C_FB = 8.55
DT_STEP = 0.005
M_N2 = 28.013


def tc_rule_vec(P, Tsat):
    s = sat(P)
    el = np.sqrt(s["k_l"] * s["rho_l"] * s["cp_l"])
    Tw = np.full_like(P, TC_N2 + 4.0)
    for _ in range(20):
        Tw = TC_N2 + el / ss_e(Tw) * (TC_N2 - Tsat)
    return Tw


def inputs():
    q = pd.read_csv(f"{OUT}/q_consistent.csv")
    recs = []
    for c in CASES:
        for st in STATIONS:
            s = q[(q.case == c) & (q.st == st)].sort_values("t")
            P = s.P.values; G = s.G.values
            recs.append(dict(case=c, st=st, t=s.t.values, Tm=s.Tm.values, Tsat=s.Tsat.values, P=P, G=G,
                             x=np.clip(s.x.values, 0, 0.9), qE=s.q_E.values,
                             Tmfb=tc_rule_vec(P, s.Tsat.values),
                             qchf_ko=katto_ohno(P, G, np.full(len(P), Z[st]))[0]))
    return recs


def h_film(d, G):
    return C_FB * d.k_vf / D * (G / d.rho_l * D * d.rho_vf * d.cp_vf / d.k_vf) ** (1 / 3)


def boiling_q(T, d, P, Tmfb, qchf, opt):
    Tsat = d.Tsat
    dT = np.maximum(T - Tsat, 0.0)
    hfb = opt["f_fb"] * h_film(d, d.G)
    pr = P / PC_N2
    A = 55 * pr**0.12 * (-np.log10(pr)) ** -0.55 * M_N2**-0.5
    hdb = 0.023 * (d.G * D / d.mu_l) ** 0.8 * (d.cp_l * d.mu_l / d.k_l) ** 0.4 * d.k_l / D
    fnb = opt["f_nb"]
    qnb = fnb * (hdb * dT + (A * dT) ** 3)
    qchf = opt["f_chf"] * qchf
    # dT_CHF: root of fnb (hdb x + A^3 x^3) = qchf (Newton, vectorised)
    x = np.full_like(dT, 10.0)
    for _ in range(30):
        f = fnb * (hdb * x + (A * x) ** 3) - qchf
        x = np.clip(x - f / (fnb * (hdb + 3 * A**3 * x**2)), 0.1, 200)
    Tchf = Tsat + x
    Tm = Tmfb + opt["dTmfb"]
    Tchf = np.minimum(Tchf, Tm - 1.0)
    # film flux at T_MFB (h evaluated with film props at T_MFB -> approximate with current h, small error)
    qfb_mfb = hfb * (Tm - Tsat)
    theta = np.clip((Tm - T) / (Tm - Tchf), 0, 1)
    qtb = qfb_mfb + (qchf - qfb_mfb) * theta ** opt["n_tb"]
    return np.where(T >= Tm, hfb * dT, np.where(T > Tchf, qtb, qnb))


def simulate(recs, opt):
    n = len(recs)
    tmax = max(r["t"][-1] for r in recs)
    tt = np.arange(0, tmax + DT_STEP + 1e-9, DT_STEP)
    gi = {k: np.array([np.interp(tt, r["t"], r[k]) for r in recs]).T for k in ("P", "Tsat", "G", "x", "Tmfb", "qchf_ko")}
    z = np.array([Z[r["st"]] for r in recs])
    T = np.array([r["Tm"][0] for r in recs])
    out = np.empty((len(tt), n))

    def rhs(T, i):
        P, Ts = gi["P"][i], gi["Tsat"][i]
        d = s6.state(P, T, Ts, gi["G"][i], gi["x"][i], z)
        d.P = P
        q = boiling_q(T, d, P, gi["Tmfb"][i], gi["qchf_ko"][i], opt)
        return (-q + 98.6 + 2.2 * (293 - T)) / (MW_A * ss_cp(T))

    for i in range(len(tt) - 1):
        out[i] = T
        k1 = rhs(T, i)
        k2 = rhs(np.maximum(T + DT_STEP * k1, 65.0), i + 1)
        T = np.maximum(T + 0.5 * DT_STEP * (k1 + k2), 65.0)
    out[-1] = T
    return tt, out


def metrics(recs, tt, sim):
    rows = []
    for j, r in enumerate(recs):
        Ts = np.interp(r["t"], tt, sim[:, j])
        e = Ts - r["Tm"]
        l2 = max(100.0, np.median(r["Tsat"]) + 8)
        row = dict(case=r["case"], st=r["st"], rms=np.sqrt(np.mean(e**2)))
        for lab, lv in (("130", 130.0), ("100", l2)):
            tm, ts = s6.t_cross(r["t"], r["Tm"], lv), s6.t_cross(r["t"], Ts, lv)
            row[f"dt{lab}"] = ts - tm
            row[f"miss{lab}"] = int(np.isnan(tm) != np.isnan(ts))
        rows.append(row)
    return pd.DataFrame(rows)


def summ(m):
    return dict(rms_mean=m.rms.mean(), rms_med=m.rms.median(), t130_mae=m.dt130.abs().mean(), t130_bias=m.dt130.mean(),
                t100_mae=m.dt100.abs().mean(), t100_bias=m.dt100.mean(), miss130=m.miss130.sum(), miss100=m.miss100.sum())


BASE = dict(f_fb=1.0, dTmfb=0.0, f_chf=1.0, n_tb=0.5, f_nb=1.0)
VARIANTS = [("film h x0.7", dict(f_fb=0.7)), ("film h x1.3", dict(f_fb=1.3)),
            ("T_MFB -5 K", dict(dTmfb=-5.0)), ("T_MFB +5 K", dict(dTmfb=5.0)),
            ("CHF x0.5", dict(f_chf=0.5)), ("CHF x1.5", dict(f_chf=1.5)),
            ("transition n=0.25", dict(n_tb=0.25)), ("transition n=1", dict(n_tb=1.0)),
            ("NB x0.5", dict(f_nb=0.5)), ("NB x2", dict(f_nb=2.0)),
            ("CHF = 100 kW/m2 const (meas.)", dict(chf_const=100e3))]


def main():
    recs = inputs()
    ko = np.concatenate([r["qchf_ko"] for r in recs])
    print("Katto-Ohno q_CHF range over records [kW/m2]:", np.percentile(ko, [5, 50, 95]).round(0) / 1e3)
    rows, sims = [], {}
    for name, ch in [("baseline", {})] + VARIANTS:
        opt = dict(BASE, **ch)
        rr = recs
        if "chf_const" in ch:
            rr = [dict(r, qchf_ko=np.full(len(r["t"]), ch["chf_const"])) for r in recs]
        tt, sim = simulate(rr, opt)
        m = metrics(recs, tt, sim)
        rows.append(dict(variant=name, **summ(m)))
        sims[name] = (tt, sim, m)
        print(rows[-1], flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/sensitivity.csv", index=False)
    sims["baseline"][2].to_csv(f"{OUT}/surrogate_final_per_station.csv", index=False)
    with open(f"{OUT}/t10_sensitivity.md", "w") as f:
        f.write(__doc__ + "\n" + md_table(res.round(3), 3) + "\n\n## baseline per station\n\n"
                + md_table(sims["baseline"][2].round(2), 3) + "\n")
    # tornado: change vs baseline
    b = res.iloc[0]
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    pairs = [("film h x0.7", "film h x1.3"), ("T_MFB -5 K", "T_MFB +5 K"), ("CHF x0.5", "CHF x1.5"),
             ("transition n=0.25", "transition n=1"), ("NB x0.5", "NB x2")]
    for ax, met, lab in ((axs[0], "rms_mean", "mean RMS(T) [K]"), (axs[1], "t130_bias", "mean t130 error [s]"),
                         (axs[2], "t100_bias", "mean t100* error [s]")):
        for i, (lo, hi) in enumerate(pairs):
            vlo = res[res.variant == lo][met].values[0] - b[met]; vhi = res[res.variant == hi][met].values[0] - b[met]
            ax.barh(i, vlo, color="C0", alpha=0.7); ax.barh(i, vhi, color="C3", alpha=0.7)
            ax.text(0, i + 0.3, f"{lo}: {vlo:+.2f} | {hi}: {vhi:+.2f}", fontsize=7, ha="center")
        ax.set_yticks(range(len(pairs))); ax.set_yticklabels([p[0].split(" x")[0].split(" -")[0].split(" n=")[0] for p in pairs])
        ax.axvline(0, color="k"); ax.set_xlabel(f"change in {lab} vs baseline ({b[met]:.2f})")
    axs[0].set_title("blue: low setting  red: high setting")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig16_sensitivity_tornado.png", dpi=100); plt.close(fig)
    # baseline traces
    tt, sim, m = sims["baseline"]
    fig, axs = plt.subplots(5, 4, figsize=(18, 17))
    for j, (ax, r) in enumerate(zip(axs.ravel(), recs)):
        ax.plot(r["t"], r["Tm"], "k", lw=2, label="measured (station-mean TC)")
        msk = tt <= r["t"][-1]
        ax.plot(tt[msk], sim[msk, j], "C3", label="final closure (0-D)")
        ax.plot(r["t"], r["Tmfb"], "C0:", lw=0.8, label="T_MFB rule")
        ax.set_title(f"case {r['case']} {r['st']}  rms {m.rms[j]:.1f} K", fontsize=9); ax.set_ylim(70, 310)
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig17_surrogate_final.png", dpi=80); plt.close(fig)
    print(res.round(2).to_string())


if __name__ == "__main__":
    main()
