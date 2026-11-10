"""Phase 2 / Task 1: ONE consistent heat-flux dataset from a lumped wall energy balance on the
station-mean outer-wall TC history, for all 10 cases; plus tests of the 0.75 anomaly.

Wall energy balance per unit INNER area (thin wall, Bi < 0.02):
    q''_E(t) = -m''_w cp(T_m) dT_m/dt + q''_par(T_m)
    m''_w = rho_w (Do^2 - Di^2)/(4 Di) = 4.20 kg/m2, cp = NIST 304SS cp(T)
    T_m   = mean of the three station TCs
    dT/dt : Savitzky-Golay, 2nd order, 9 points (0.60 s) [film boiling]; 5 points (0.33 s) used for
            CHF peaks (a01). Sensitivity to 5/9/15 points reported.
    q''_par(T) = a + b (293 - T): parasitic heat leak into the wall (radiation from the 293 K vacuum
            chamber + residual-gas conduction), fitted to (q_inverse - q_lump) of the five
            self-consistent cases (4,7,8,9,10) -- i.e. the authors' own parasitic model; compared with
            the radiation-only estimate eps sigma (T_ch^4 - T^4) Do/Di, eps = 0.453 (Kassemi et al.).
Inner-surface temperature (quasi-steady radial profile with uniform storage sink, adiabatic outer
surface):  T_i = T_m - q''_E r_i/k(T) [ r_o^2 ln(r_o/r_i)/(r_o^2 - r_i^2) - 1/2 ]  (= q'' * 0.254 mm / k,
0.2 K in film boiling, ~2.5 K at the CHF peak).  Radial diffusion time delta^2/alpha ~ 0.07 s
(one sample) -> no time-shift correction applied.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (CASES, D, DO, FIG, MW_A, OUT, RHO_SS, STATIONS, Z, load, md_table, smooth, ss_cp, ss_k)
from prep import film_mask, inverse_scale, quench_info, station_data, t_valid_end
import fb_models as fm

RI, RO = D / 2, DO / 2
GEO = RI * (RO**2 * np.log(RO / RI) / (RO**2 - RI**2) - 0.5)   # [m] ~ 2.54e-4
GOOD = [4, 7, 8, 9, 10]
BAD = [1, 2, 3, 5, 6]
SIGMA = 5.670374e-8


def dTdt(T, win):
    return smooth(np.asarray(T, float), win, 2, deriv=1)


def parasitic_fit(win=9):
    """fit q_inv - q_lump(Ti) = a + b (293 - Ti) on the consistent cases, film boiling."""
    X, Y = [], []
    for c in GOOD:
        for st in STATIONS:
            s = station_data(c, st)
            m = film_mask(c, st) & (s.t.values > 1.0)
            ql = -MW_A * ss_cp(s.Ti.values) * dTdt(s.Ti.values, win)
            X.append(293.0 - s.Ti.values[m]); Y.append(s.q.values[m] - ql[m])
    X = np.concatenate(X); Y = np.concatenate(Y)
    b, a = np.polyfit(X, Y, 1)
    return a, b


def build(win=9, a=None, b=None, tau=fm.TAU_G):
    rows = []
    for c in CASES:
        d = load(c)
        Gl = fm.g_testsection(c, tau)
        for st in STATIONS:
            s = station_data(c, st)
            n = len(s)
            Tm = s[["TC0", "TC1", "TC2"]].mean(axis=1).values
            qs = -MW_A * ss_cp(Tm) * dTdt(Tm, win)
            qp = a + b * (293.0 - Tm)
            qE = qs + qp
            Ti = Tm - qE * GEO / ss_k(Tm)
            df = pd.DataFrame(dict(case=c, st=st, t=s.t.values, Tm=Tm, Tw=Ti, Tsat=s.Tsat.values, P=s.P.values,
                                   q_E=qE, q_store=qs, q_par=qp, q_pub=s.q.values, Ti_pub=s.Ti.values,
                                   G=np.clip(Gl[:n], 5, None), G_raw=np.clip(d.G_kgm2s.values[:n], 5, None),
                                   z=Z[st]))
            df["dT"] = df.Tw - df.Tsat
            df["h_E"] = df.q_E / df.dT
            df["h_pub"] = df.q_pub / (df.Ti_pub - df.Tsat)
            df["fb"] = film_mask(c, st)
            rows.append(df)
    df = pd.concat(rows, ignore_index=True)
    # equilibrium quality from the consistent fluxes (same quasi-steady energy balance as phase 1)
    out = []
    for c in CASES:
        su = df[(df.case == c) & (df.st == "up")].reset_index(drop=True)
        sd = df[(df.case == c) & (df.st == "dn")].reset_index(drop=True)
        n = min(len(su), len(sd))
        for st, s in (("up", su), ("dn", sd)):
            from common import sat
            hfg = sat(s.P.values)["h_fg"]
            qu = np.r_[su.q_E.values[:n], np.full(len(s) - n, su.q_E.values[n - 1])][:len(s)]
            qd = np.r_[sd.q_E.values[:n], np.full(len(s) - n, sd.q_E.values[n - 1])][:len(s)]
            s = s.copy()
            s["x"] = fm.x_eq(np.clip(qu, 0, None), np.clip(qd, 0, None), s.G.values, hfg, st)
            out.append(s)
    return pd.concat(out, ignore_index=True)


def anomaly_tests(df, a, b):
    """Is the 0.75 factor a scaling of the storage term only (wall heat capacity in the inverse),
    of the total q, or something physical?"""
    rows = []
    for c in CASES:
        for st in STATIONS:
            s = df[(df.case == c) & (df.st == st) & df.fb & (df.t > 1)]
            # q_pub = alpha * q_store + beta * q_par_model
            X = np.c_[s.q_store.values, s.q_par.values]
            coef = np.linalg.lstsq(X, s.q_pub.values, rcond=None)[0]
            # q_pub = k * q_E  (total scaling)
            k_tot = np.sum(s.q_pub * s.q_E) / np.sum(s.q_E**2)
            r_store = np.sqrt(np.mean((s.q_pub - X @ coef) ** 2))
            r_tot = np.sqrt(np.mean((s.q_pub - k_tot * s.q_E) ** 2))
            rows.append(dict(case=c, st=st, alpha_store=coef[0], beta_par=coef[1], rms_storefit=r_store,
                             k_total=k_tot, rms_totalfit=r_tot,
                             Ti_pub_minus_Tm=np.median(s.Ti_pub - s.Tm)))
    an = pd.DataFrame(rows)
    # implied wall thickness if the storage factor were geometric (same rho, cp, OD)
    alpha_bad = an[an.case.isin(BAD)].alpha_store.mean()
    alpha_good = an[an.case.isin(GOOD)].alpha_store.mean()
    ratio = alpha_bad / alpha_good
    from scipy.optimize import brentq
    m0 = (DO**2 - D**2) / (4 * D)
    di = brentq(lambda di: (DO**2 - di**2) / (4 * di) / m0 - ratio, D, DO - 1e-6)
    thick = (DO - di) / 2
    return an, ratio, thick


def holdout_test(df):
    """Arbiter: fit ln h = a + b ln G + c ln P on the 5 self-consistent cases only, then predict the
    5 anomalous cases with (i) h_E (0.51 mm wall) and (ii) the published h.  If the anomalous runs
    really had a thinner (0.38 mm) wall, the published h would be the correct one and match."""
    fb = df[df.fb]
    w = 1.0 / fb.groupby(["case", "st"]).t.transform("size").values
    g = fb.case.isin(GOOD).values
    X = np.c_[np.ones(len(fb)), np.log(fb.G), np.log(fb.P), (fb.st == "dn")]
    rows = []
    for lab, form in (("free lnG,lnP,dn", X), ("F2 form sqrt(k_v cp_v G/D)", None)):
        if form is None:
            from common import sat
            s = sat(fb.P.values)
            unit = np.sqrt(s["k_v"] * s["cp_v"] * fb.G.values / D)
            C = np.exp(np.sum(w[g] * np.log(fb.h_E.values[g] / unit[g])) / np.sum(w[g]))
            pred = C * unit
        else:
            bcoef = np.linalg.lstsq(form[g] * np.sqrt(w[g])[:, None], np.log(fb.h_E.values[g]) * np.sqrt(w[g]), rcond=None)[0]
            pred = np.exp(form @ bcoef)
        for c in BAD:
            m = fb.case.values == c
            rows.append(dict(model=lab, case=c, G=fb.G[m].median(), P_kPa=fb.P[m].median() / 1e3,
                             ratio_hE_over_pred=np.exp(np.mean(np.log(fb.h_E.values[m] / pred[m]))),
                             ratio_hpub_over_pred=np.exp(np.mean(np.log(fb.h_pub.values[m] / pred[m])))))
    return pd.DataFrame(rows)


def main():
    a, b = parasitic_fit()
    print(f"parasitic fit: q_par = {a:.1f} + {b:.3f} (293 - T)  W/m2")
    rad = lambda T: 0.453 * SIGMA * (293.0**4 - T**4) * DO / D
    print("radiation-only at 250/180/130 K:", rad(250), rad(180), rad(130), " fit:", a + b * 43, a + b * 113, a + b * 163)
    df = build(9, a, b)
    df.to_csv(f"{OUT}/q_consistent.csv", index=False)
    # smoothing sensitivity
    sens = {}
    for win in (5, 9, 15):
        dw = build(win, a, b)
        fb = dw[dw.fb & (dw.t > 1)]
        sens[win] = fb.h_E.values
    ref = sens[9]
    sens_txt = "; ".join(f"win {w}: median h ratio to win9 {np.median(v/ref):.4f}, rms ln ratio {np.sqrt(np.mean(np.log(np.clip(v,1,None)/np.clip(ref,1,None))**2)):.3f}"
                         for w, v in sens.items())
    print(sens_txt)
    an, ratio, thick = anomaly_tests(df, a, b)
    ho = holdout_test(df)
    # per-case comparison table q_E vs published
    comp = []
    for c in CASES:
        for st in STATIONS:
            s = df[(df.case == c) & (df.st == st) & df.fb & (df.t > 1)]
            comp.append(dict(case=c, st=st, n=len(s), hE_med=s.h_E.median(), hpub_med=s.h_pub.median(),
                             ratio_pub_E=np.sum(s.q_pub) / np.sum(s.q_E), Tw_minus_Tipub_med=np.median(s.Tw - s.Ti_pub),
                             phase1_k=inverse_scale(c)))
    comp = pd.DataFrame(comp)
    pd.set_option("display.width", 220)
    print(comp.round(3).to_string()); print(an.round(3).to_string()); print(ho.round(3).to_string())
    # radial-conduction length implied by the published Ti: L = (T_m - Ti_pub) k / q_pub
    fbq = df[df.fb & (df.t > 1) & (df.q_pub > 3000)]
    Lr = ((fbq.Tm - fbq.Ti_pub) * ss_k(fbq.Tm) / fbq.q_pub * 1e6).groupby(fbq.case).median()
    print("implied radial length [um] per case:", Lr.round(1).to_dict(), " 0.51-mm geometry:", GEO * 1e6,
          " 0.388-mm geometry would give ~", 0.388 / 0.51 * GEO * 1e6)
    with open(f"{OUT}/t07_consistent_q.md", "a") as f:
        f.write(f"\nradial-conduction length implied by published Ti (um): {Lr.round(1).to_dict()}; "
                f"0.51-mm wall geometry gives {GEO*1e6:.0f} um, a 0.388-mm wall would give ~{0.388/0.51*GEO*1e6:.0f} um\n")
    print(f"storage-term factor bad/good = {ratio:.4f}; equivalent wall thickness {thick*1e3:.3f} mm (same OD)")
    with open(f"{OUT}/t07_consistent_q.md", "w") as f:
        f.write(__doc__ + "\n")
        f.write(f"parasitic fit: q_par = {a:.1f} + {b:.3f} (293 - T) W/m2 "
                f"(radiation-only eps=0.453: {rad(250):.0f}/{rad(180):.0f}/{rad(130):.0f} W/m2 at 250/180/130 K; "
                f"fit {a+b*43:.0f}/{a+b*113:.0f}/{a+b*163:.0f})\n\n")
        f.write("smoothing sensitivity (film boiling h_E): " + sens_txt + "\n\n")
        f.write("## consistent vs published (film boiling, t>1 s)\n\n" + md_table(comp.round(3), 4) + "\n\n")
        f.write("## anomaly: decomposition q_pub = alpha q_store + beta q_par\n\n" + md_table(an.round(3), 4) + "\n\n")
        f.write(f"storage factor bad/good = {ratio:.4f} -> equivalent wall thickness {thick*1e3:.3f} mm at the same OD\n\n")
        f.write("## hold-out arbiter (law fitted on cases 4,7-10, predicting 1,2,3,5,6)\n\n" + md_table(ho.round(3), 4) + "\n")
    # figure
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    for c in CASES:
        s = df[(df.case == c) & (df.st == "up") & (df.t > 0.5)]
        col = "C3" if c in BAD else "C0"
        axs[0].plot(s.Tm, s.q_pub / s.q_E.clip(lower=100), color=col, lw=0.6, alpha=0.7)
    axs[0].set_ylim(0.5, 1.3); axs[0].invert_xaxis()
    axs[0].set_xlabel("station-mean TC T [K]"); axs[0].set_ylabel("q_published / q_E (upstream)")
    axs[0].set_title("red: cases 1,2,3,5,6  blue: 4,7-10")
    for c in BAD:
        for lab, mk in (("ratio_hE_over_pred", "o"), ("ratio_hpub_over_pred", "x")):
            sub = ho[(ho.case == c)]
            axs[1].plot(sub.G, sub[lab], mk, color="C3" if mk == "x" else "C0")
    axs[1].axhline(1, color="k"); axs[1].set_xscale("log")
    axs[1].set_xlabel("G [kg/m2s]"); axs[1].set_ylabel("measured / predicted (law fitted on 4,7-10)")
    axs[1].set_title("o: h_E (0.51 mm wall)   x: published h")
    T = np.linspace(100, 290)
    axs[2].plot(T, a + b * (293 - T), label="fitted parasitic (inverse model, good cases)")
    axs[2].plot(T, rad(T), label="radiation only, eps=0.453")
    axs[2].set_xlabel("T [K]"); axs[2].set_ylabel("q_par [W/m2]"); axs[2].legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/fig13_consistent_q.png", dpi=100); plt.close(fig)


if __name__ == "__main__":
    main()
