"""Goal 6: 0-D lumped-wall surrogate at each station.

m''_w cp(T) dT/dt = -q_conv(T, t) + q_par,     m''_w = rho_w (Do^2 - Di^2)/(4 Di) = {MW_A:.3f} kg/m2
driven by measured Tsat(t), P(t) at the station and the lag-corrected G(t) (tau = 2.5 s) or raw G(t).
Boiling curve:
  T > T_MFB          : film boiling, q = h_FB(closure) (T - Tsat)
  T_CHF < T < T_MFB  : transition, q = q_FB(T_MFB) + (q_CHF - q_FB(T_MFB)) theta^0.5,
                       theta = (T_MFB - T)/(T_MFB - T_CHF)       (theta^0.5 = best fit in a05)
  T < T_CHF          : nucleate, q = q_CHF ((T - Tsat)/dT_CHF)^2.5
  q_CHF = 100 kW/m2, dT_CHF = 24 K (data means, a05); q_par = +0.33 kW/m2 (heat leak, a01).
All 20 station records are integrated simultaneously (Heun, dt = 5 ms).
"""
import json
import os
import types

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from common import CASES, D, FIG, MW_A, OUT, STATIONS, TC_N2, T_hn_lienhard, Z, contact_wall_T, md_table, sat, ss_cp, vap
import fb_models as fm
from prep import station_data

Q_PAR = 330.0
Q_CHF = 100e3
DT_CHF = 24.0
N_NB = 2.5
DT_STEP = 0.005

# ------------------------------------------------------------- property tables
_PG = np.arange(60e3, 900e3 + 1, 5e3)
_SAT = sat(_PG)
_PV = np.arange(60e3, 900e3 + 1, 30e3)
_TV = np.arange(64.0, 312.0, 1.0)
_VT = {}
for key in ("rho", "mu", "k", "cp"):
    _VT[key] = np.empty((len(_PV), len(_TV)))
for i, p in enumerate(_PV):
    v = vap(np.full(len(_TV), p), _TV)
    for key in _VT:
        _VT[key][i] = v[key]
_VI = {k: RegularGridInterpolator((_PV, _TV), v, bounds_error=False, fill_value=None) for k, v in _VT.items()}


def state(P, T, Tsat, G, x, z):
    """namespace with the fields fb_models closures use."""
    d = types.SimpleNamespace()
    for k, v in _SAT.items():
        setattr(d, k, np.interp(P, _PG, v))
    d.Tsat = Tsat
    d.dT = np.maximum(T - Tsat, 0.5)
    Tf = 0.5 * (T + Tsat)
    pts = np.c_[P, np.maximum(Tf, Tsat + 0.1)]
    d.k_vf, d.rho_vf, d.cp_vf, d.mu_vf = (_VI[k](pts) for k in ("k", "rho", "cp", "mu"))
    pw = np.c_[P, np.maximum(T, Tsat + 0.1)]
    d.Pr_vw = _VI["cp"](pw) * _VI["mu"](pw) / _VI["k"](pw)
    d.x = x
    d.z = z
    d.G = G
    return d


def closures(params):
    pC = params
    return {
        "Bromley C=2.7": lambda d, G: 2.7 * fm.bromley_unit(d, G, "film"),
        f"Bromley C={pC['Bromley, C fit (film props)'][0]:.2f} (fit)":
            lambda d, G: pC["Bromley, C fit (film props)"][0] * fm.bromley_unit(d, G, "film"),
        f"Penetration C={pC['Penetration, C fit'][0]:.2f} (fit)":
            lambda d, G: pC["Penetration, C fit"][0] * fm.penetration_unit(d, G),
        "Pe_l power law (2 fit)": lambda d, G: pC["Pe_l power law: C (k_v/D) Pe_l^n"][0] * d.k_vf / D
            * (G / d.rho_l * D * d.rho_vf * d.cp_vf / d.k_vf) ** pC["Pe_l power law: C (k_v/D) Pe_l^n"][1],
        "Darr 2016 low-g (published)": lambda d, G: fm.darr2016(d, G),
        "Dittus-Boelter vapour (sat)": lambda d, G: fm.dittus_boelter_v(d, G, "sat"),
        "Miropolskii (x_eq)": lambda d, G: fm.miropolskii(d, G, d.x),
    }


def tmfb_rules(st_meas):
    return {
        "T_MFB measured (oracle)": lambda P, Tsat, G, key: st_meas.get(key, 134.0),
        # (stations that never quench get 134 K)
        "T_MFB = 134 K (low-g const.)": lambda P, Tsat, G, key: 134.0,
        "T_MFB contact-T_hn (~117 K)": lambda P, Tsat, G, key: contact_wall_T(T_hn_lienhard(Tsat), Tsat, P)[0],
        "T_MFB Darr-Hartwig We": lambda P, Tsat, G, key: 0.844 * TC_N2 * (1 + 0.060 * (G**2 * D / (np.interp(P, _PG, _SAT["rho_l"]) * np.interp(P, _PG, _SAT["sigma"]))) ** 0.208),
    }


def load_inputs(use_lag=True):
    """per-station time series on a common time base (held at the last value beyond the record)."""
    full = fm.dataset(all_samples=True)
    recs = []
    for c in CASES:
        for st in STATIONS:
            s = station_data(c, st)
            sub = full[(full.case == c) & (full.st == st)].sort_values("t")
            recs.append(dict(case=c, st=st, t=s.t.values, Ti=s.Ti.values, Tsat=s.Tsat.values, P=s.P.values,
                             q=s.q_ec.values,
                             G=np.interp(s.t.values, sub.t.values, sub.G.values if use_lag else sub.G_raw.values),
                             x=np.interp(s.t.values, sub.t.values, np.clip(sub.x.values, 0, 0.9))))
    return recs


def boiling_q(T, d, h_fb_fun, Tmfb):
    Tsat = d.Tsat
    qfb = h_fb_fun(d, d.G) * np.maximum(T - Tsat, 0.0)
    # film flux evaluated at T_MFB (for the transition anchor)
    dm = types.SimpleNamespace(**vars(d))
    dm.dT = np.maximum(Tmfb - Tsat, 0.5)
    Tf = 0.5 * (Tmfb + Tsat)
    pts = np.c_[d.P, Tf]
    dm.k_vf, dm.rho_vf, dm.cp_vf, dm.mu_vf = (_VI[k](pts) for k in ("k", "rho", "cp", "mu"))
    pw = np.c_[d.P, Tmfb]
    dm.Pr_vw = _VI["cp"](pw) * _VI["mu"](pw) / _VI["k"](pw)
    q_mfb = h_fb_fun(dm, d.G) * (Tmfb - Tsat)
    Tchf = Tsat + DT_CHF
    theta = np.clip((Tmfb - T) / np.maximum(Tmfb - Tchf, 1e-3), 0, 1)
    qtb = q_mfb + (Q_CHF - q_mfb) * np.sqrt(theta)
    qnb = Q_CHF * (np.maximum(T - Tsat, 0) / DT_CHF) ** N_NB
    return np.where(T >= Tmfb, qfb, np.where(T > Tchf, qtb, qnb))


def simulate(recs, h_fb_fun, tmfb_fun, oracle_q=False):
    n = len(recs)
    tmax = max(r["t"][-1] for r in recs)
    tt = np.arange(0, tmax + DT_STEP + 1e-9, DT_STEP)
    keys = [(r["case"], r["st"]) for r in recs]
    z = np.array([Z[r["st"]] for r in recs])
    grid = {k: np.array([np.interp(tt, r["t"], r[k]) for r in recs]).T for k in ("P", "Tsat", "G", "x", "q")}
    if not oracle_q:
        Tm = np.empty((len(tt), n))
        for j, r in enumerate(recs):
            v = tmfb_fun(r["P"], r["Tsat"], r["G"], keys[j])
            Tm[:, j] = np.interp(tt, r["t"], np.broadcast_to(np.asarray(v, float), r["t"].shape))
    T = np.array([r["Ti"][0] for r in recs], float)
    out = np.empty((len(tt), n))

    def rhs(T, i):
        if oracle_q:
            q = grid["q"][i]
        else:
            P, Ts = grid["P"][i], grid["Tsat"][i]
            d = state(P, T, Ts, grid["G"][i], grid["x"][i], z)
            d.P = P
            q = boiling_q(T, d, h_fb_fun, Tm[i])
        return (-q + Q_PAR) / (MW_A * ss_cp(T))

    for i in range(len(tt) - 1):
        out[i] = T
        k1 = rhs(T, i)
        k2 = rhs(np.maximum(T + DT_STEP * k1, 65.0), i + 1)
        T = np.maximum(T + 0.5 * DT_STEP * (k1 + k2), 65.0)
    out[-1] = T
    return tt, out


def t_cross(t, T, level):
    k = np.where(T < level)[0]
    if len(k) == 0:
        return np.nan
    k = k[0]
    if k == 0:
        return t[0]
    return t[k - 1] + (T[k - 1] - level) / (T[k - 1] - T[k]) * (t[k] - t[k - 1])


def metrics(recs, tt, sim):
    rows = []
    for j, r in enumerate(recs):
        Ts_ = np.interp(r["t"], tt, sim[:, j])
        e = Ts_ - r["Ti"]
        lvl2 = max(100.0, np.median(r["Tsat"]) + 8.0)
        row = dict(case=r["case"], st=r["st"], rms_K=np.sqrt(np.mean(e**2)), max_abs_K=np.abs(e).max(),
                   T_end_err=e[-1])
        for lab, lvl in (("130", 130.0), ("100*", lvl2)):
            tm = t_cross(r["t"], r["Ti"], lvl)
            ts = t_cross(r["t"], Ts_, lvl)
            row[f"t{lab}_meas"] = tm
            row[f"t{lab}_sim"] = ts
            row[f"dt{lab}"] = ts - tm
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(m):
    both130 = m.dropna(subset=["t130_meas", "t130_sim"])
    both100 = m.dropna(subset=["t100*_meas", "t100*_sim"])
    miss130 = int((m.t130_meas.isna() != m.t130_sim.isna()).sum())
    miss100 = int((m["t100*_meas"].isna() != m["t100*_sim"].isna()).sum())
    return dict(rms_K_mean=m.rms_K.mean(), rms_K_median=m.rms_K.median(), maxabs_K=m.max_abs_K.max(),
                mae_t130_s=np.abs(both130.dt130).mean() if len(both130) else np.nan,
                bias_t130_s=both130.dt130.mean() if len(both130) else np.nan, n_t130=len(both130),
                mismatch130=miss130,
                mae_t100_s=np.abs(both100["dt100*"]).mean() if len(both100) else np.nan,
                bias_t100_s=both100["dt100*"].mean() if len(both100) else np.nan, n_t100=len(both100),
                mismatch100=miss100)


def main():
    params = json.load(open(f"{OUT}/film_boiling_params.json"))
    stm = pd.read_csv(f"{OUT}/mfb_chf_station.csv")
    st_meas = {(int(r.case), r.st): r.T_mfb for _, r in stm.iterrows()}
    rows = []
    sims = {}
    for use_lag in (True, False):
        recs = load_inputs(use_lag)
        if use_lag:
            tt, so = simulate(recs, None, None, oracle_q=True)
            m = metrics(recs, tt, so)
            rows.append(dict(G="lag", closure="measured q'' (oracle check)", rewet="-", **summarize(m)))
        for cn, cf in closures(params).items():
            rules = tmfb_rules(st_meas)
            if not use_lag:
                rules = {k: v for k, v in rules.items() if k == "T_MFB = 134 K (low-g const.)"}
            for rn, rf in rules.items():
                tt, sim = simulate(recs, cf, rf)
                m = metrics(recs, tt, sim)
                fn = f"surrogate_{'lag' if use_lag else 'raw'}_{cn[:12]}_{rn[:12]}.csv".replace(" ", "_").replace("/", "_")
                os.makedirs(f"{OUT}/surrogate_runs", exist_ok=True)
                m.to_csv(f"{OUT}/surrogate_runs/{fn}", index=False)
                rows.append(dict(G="lag" if use_lag else "raw", closure=cn, rewet=rn, **summarize(m)))
                if use_lag and cn.startswith("Pe_l") and rn.startswith("T_MFB = 134"):
                    best_m = m
                sims[(use_lag, cn, rn)] = (tt, sim)
                print(rows[-1])
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/surrogate_summary.csv", index=False)
    with open(f"{OUT}/t06_surrogate.md", "w") as f:
        f.write(f"wall areal mass m''_w = {MW_A:.3f} kg/m2 (rho 7900, Do 12.70, Di 11.68 mm); "
                f"m''_w cp(290 K) = {MW_A * ss_cp(290):.0f} J/m2K, (120 K) = {MW_A * ss_cp(120):.0f} J/m2K\n\n")
        f.write(md_table(res.round(3), 3) + "\n")
        f.write("\n## per station: Pe_l power law + T_MFB = 134 K, lag-corrected G\n\n")
        f.write(md_table(best_m.round(2), 4) + "\n")
    # figure: best few closures with 134 K rule vs data
    recs = load_inputs(True)
    show = [k for k in closures(params)][:6]
    fig, axs = plt.subplots(5, 4, figsize=(20, 20))
    for j, (ax, r) in enumerate(zip(axs.ravel(), recs)):
        ax.plot(r["t"], r["Ti"], "k", lw=2, label="measured Ti")
        ax.plot(r["t"], r["Tsat"], "b:", lw=1)
        for i, cn in enumerate(show):
            tt, sim = sims[(True, cn, "T_MFB = 134 K (low-g const.)")]
            ax.plot(tt[tt <= r["t"][-1]], sim[tt <= r["t"][-1], j], lw=1, color=f"C{i}", label=cn)
        ax.set_title(f"case {r['case']} {r['st']}")
        ax.set_ylim(70, 310)
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("0-D lumped-wall surrogate, lag-corrected G, T_MFB = 134 K, theta^0.5 transition")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig12_surrogate.png", dpi=80)
    plt.close(fig)
    pd.set_option("display.width", 250)
    print(res.round(2).to_string())


if __name__ == "__main__":
    main()
