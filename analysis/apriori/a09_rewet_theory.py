"""Phase 2 / Task 3: zero-constant contact-temperature rewet rule
    rewet when the liquid-wall contact temperature falls to T_c:
    T_w,MFB = T_c + beta (T_c - T_l),  beta = e_l / e_w(T_w,MFB)  (iterated),  T_l = Tsat(P_local)
e_l = sqrt(k rho cp) of saturated liquid at P; e_w from NIST fits (304SS; OFHC copper RRR=50).
Compared with measured T_MFB, Henry(1974), T_hn+19 K, constant 134 K; upstream-offset diagnosis;
predictions for Cu/LN2, Cu/LH2, FEP-coated SS (Hartwig 2022 / Chung 2021), 1-g Darr tests,
and checks on digitized external data (Darr test 31/40 1 g, Hartwig F2 0.05 g, NBS 9264 Cu LN2).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

from common import (FIG, OUT, ROOT, TC_N2, T_hn_lienhard, md_table, sat, smooth, ss_e)
from prep import detect_quench, station_data
from a04_rewet import henry
from common import contact_wall_T

RHO_CU = 8960.0


def cu_cp(T):
    """NIST cryogenic database, OFHC copper specific heat, log10 polynomial (4-300 K)."""
    c = [-1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0.0]
    x = np.log10(np.clip(T, 4, 300))
    return 10 ** sum(ci * x**i for i, ci in enumerate(c))


def cu_k(T, ):
    """NIST OFHC copper thermal conductivity, RRR = 50 (rational fit in T^0.5)."""
    a, b, c, d, e, f, g, h, i = 1.8743, -0.41538, -0.6018, 0.13294, 0.26426, -0.0219, -0.051276, 0.0014871, 0.003723
    T = np.clip(T, 4, 300)
    s = np.sqrt(T)
    return 10 ** ((a + c * s + e * T + g * T * s + i * T**2) / (1 + b * s + d * T + f * T * s + h * T**2))


def cu_e(T):
    return np.sqrt(cu_k(T) * RHO_CU * cu_cp(T))


FEP = dict(k=0.195, rho=2150.0, cp=700.0)   # DuPont k (Chung 2021 ref 32); cp at ~130 K assumed (+-30 %)
E_FEP = np.sqrt(FEP["k"] * FEP["rho"] * FEP["cp"])


def liq(P, fluid="Nitrogen"):
    g = lambda k: CP.PropsSI(k, "P", P, "Q", 0, fluid)
    return dict(Tsat=g("T"), e=np.sqrt(g("L") * g("D") * g("C")), Tc=CP.PropsSI("Tcrit", fluid))


def tc_rule(P, fluid="Nitrogen", e_wall=ss_e, Tl=None):
    P = np.atleast_1d(np.asarray(P, float))
    out, betas = [], []
    for p in P:
        L = liq(p, fluid)
        T_l = L["Tsat"] if Tl is None else Tl
        Tw = L["Tc"] + 5
        for _ in range(50):
            beta = L["e"] / (e_wall(Tw) if callable(e_wall) else e_wall)
            Tw = L["Tc"] + beta * (L["Tc"] - T_l)
        out.append(Tw); betas.append(beta)
    return np.array(out), np.array(betas)


def onset_generic(t, T, Tsat, base=(3.0, 1.0), win=5):
    """same 2x-cooling-rate onset as prep.detect_quench but with a scalable baseline window."""
    dt = np.median(np.diff(t))
    r = -smooth(T, win, 2, deriv=1) / 1.0 * (1 / 15) / dt   # smooth() assumes dt=1/15
    k = np.argmax(np.where(t > t[0] + base[0], r, -np.inf))
    m = (t >= t[k] - base[0]) & (t <= t[k] - base[1])
    rfb = max(np.median(r[m]), 1e-3)
    j = k
    while j > 0 and r[j] > 2 * rfb:
        j -= 1
    a = np.clip((2 * rfb - r[j]) / max(r[j + 1] - r[j], 1e-9), 0, 1)
    return T[j] + a * (T[j + 1] - T[j]), t[k], r[k] / rfb


def external_checks():
    rows = []
    d = ROOT + "/data/"
    # Darr 2016 IJHMT 1 g upflow tests 31/40 (station means)
    for test, P in ((31, 314.6e3), (40, 325.9e3)):
        df = pd.read_csv(d + f"darr2016ijhmt_1g_vertical/test{test}_LN2_1g_upflow.csv", comment="#")
        for col, st in ((df.columns[1], "up"), (df.columns[2], "dn")):
            T = df[col].values; t = df.t_s.values
            Tm, tq, ratio = onset_generic(t, T, None, base=(3.0, 1.0), win=7)
            pr, b = tc_rule(P)
            rows.append(dict(data=f"Darr 1g test {test} {st}", fluid="LN2", wall="304SS", P_kPa=P / 1e3,
                             T_mfb_meas=Tm, peak_rate_ratio=ratio, Tc_rule=pr[0], beta=b[0]))
    # Hartwig 2022 F2 (bare SS, 0.05 g, exit station mean)
    df = pd.read_csv(d + "hartwig2022_F2/F2_bare_continuous_0p05g.csv", comment="#")
    t = df.t_s.values; T = df.Twall_exit_avg_K.values
    Tm, tq, ratio = onset_generic(t, T, None, base=(3.0, 1.0), win=5)
    P = np.median(df.P_exit_kPa.values[(t > 3) & (t < tq)]) * 1e3
    pr, b = tc_rule(P)
    rows.append(dict(data="Hartwig 2022 F2 exit (0.05 g)", fluid="LN2", wall="304SS", P_kPa=P / 1e3, T_mfb_meas=Tm,
                     peak_rate_ratio=ratio, Tc_rule=pr[0], beta=b[0]))
    # NBS 9264 copper, LN2 saturated runs (hand-clicked markers -> knee only indicative)
    for fig, Pd in ((10, 2.5), (11, 3.4), (12, 5.9)):
        df = pd.read_csv(d + f"nbs9264/nbs9264_fig{fig}_ln2_{str(Pd).replace('.', 'p')}atm.csv", comment="#")
        for zst in (6.1, 24.4):
            s = df[df.station_m == zst]
            Tm, tq, ratio = onset_generic(s.time_s.values, s.T_wall_K.values, None, base=(30.0, 8.0), win=5)
            # local pressure: linear between P_drive and 0.82 atm exit (61 m) -- rough
            P = (Pd + (0.82 - Pd) * zst / 61.0) * 101325
            pr, b = tc_rule(P, e_wall=cu_e)
            rows.append(dict(data=f"NBS9264 fig{fig} z={zst} m", fluid="LN2", wall="Cu", P_kPa=P / 1e3, T_mfb_meas=Tm,
                             peak_rate_ratio=ratio, Tc_rule=pr[0], beta=b[0]))
    return pd.DataFrame(rows)


def main():
    rt = pd.read_csv(f"{OUT}/rewet_table.csv")
    tr = pd.read_csv(f"{OUT}/mfb_chf_per_trace.csv")
    P = rt.P_kPa.values * 1e3
    rt["Tc_rule"], rt["beta"] = tc_rule(P)
    Ts = rt.Tsat.values
    rt["Henry"] = henry(contact_wall_T(T_hn_lienhard(Ts), Ts, P)[0], P)
    rt["Thn+19"] = T_hn_lienhard(Ts) + 18.8
    rt["const134"] = 134.0
    # alternative measured definitions (for the upstream-offset question)
    defs = {}
    for (c, st), g in tr[tr.trace != "Ti"].groupby(["case", "st"]):
        defs[(c, st)] = dict(median_TC=g.T_mfb.median(), first_TC=g.loc[g.t_mfb.idxmin(), "T_mfb"],
                             max_TC=g.T_mfb.max(), min_TC=g.T_mfb.min(), thr15=g.T_mfb_15.median(),
                             thr3=g.T_mfb_3.median(), curv=g.T_mfb_curv.median())
    for (c, st) in defs:
        s = station_data(c, st)
        Tm = s[["TC0", "TC1", "TC2"]].mean(axis=1).values
        q = detect_quench(s.t.values, Tm, s.Tsat.values, s.attrs["t_end"])
        defs[(c, st)]["station_mean_TC"] = q["T_mfb"] if q else np.nan
        ti = tr[(tr.case == c) & (tr.st == st) & (tr.trace == "Ti")]
        defs[(c, st)]["Ti_inverse"] = ti.T_mfb.values[0] if len(ti) else np.nan
    dd = pd.DataFrame(defs).T
    dd.index.names = ["case", "st"]
    paired = [c for c in rt.case.unique() if (rt.case == c).sum() == 2]
    off = {k: float(np.mean([dd.loc[(c, "up"), k] - dd.loc[(c, "dn"), k] for c in paired])) for k in dd.columns}
    # other station differences that could explain it
    diag = []
    for c in paired:
        u, d_ = rt[(rt.case == c) & (rt.st == "up")].iloc[0], rt[(rt.case == c) & (rt.st == "dn")].iloc[0]
        diag.append(dict(case=c, dT_up_minus_dn=u.T_mfb - d_.T_mfb, dP_kPa=u.P_kPa - d_.P_kPa,
                         dTsat=u.Tsat - d_.Tsat, dG=u.G_lag - d_.G_lag, t_mfb_up=u.t_mfb, t_mfb_dn=d_.t_mfb,
                         rule_up_minus_dn=u.Tc_rule - d_.Tc_rule))
    diag = pd.DataFrame(diag)
    # performance
    rows = []
    for k in ("Tc_rule", "Henry", "Thn+19", "const134"):
        e = rt[k] - rt.T_mfb
        rows.append(dict(model=k, mean_err=e.mean(), MAE=e.abs().mean(), rmse=np.sqrt((e**2).mean()),
                         mean_err_up=e[rt.st == "up"].mean(), mean_err_dn=e[rt.st == "dn"].mean(),
                         MAE_up=e[rt.st == "up"].abs().mean(), MAE_dn=e[rt.st == "dn"].abs().mean()))
    perf = pd.DataFrame(rows)
    # predictions for other systems
    preds = []
    for P in (101325, 300e3, 600e3):
        for lab, fl, ew in (("LN2 / 304SS", "Nitrogen", ss_e), ("LN2 / Cu (OFHC)", "Nitrogen", cu_e),
                            ("LH2 / Cu (OFHC)", "Hydrogen", cu_e), ("LH2 / 304SS", "Hydrogen", ss_e),
                            ("LN2 / FEP surface (semi-inf.)", "Nitrogen", E_FEP)):
            T, b = tc_rule(P, fl, ew)
            preds.append(dict(system=lab, P_kPa=P / 1e3, Tsat=liq(P, fl)["Tsat"], Tc=liq(P, fl)["Tc"], beta=b[0],
                              T_wMFB=T[0], dT_MFB=T[0] - liq(P, fl)["Tsat"]))
    preds = pd.DataFrame(preds)
    # coating: surface vs metal temperature at rewet, 1L..4L, q_MFB ~ 8 kW/m2 ; contact time for semi-inf validity
    alpha_c = FEP["k"] / (FEP["rho"] * FEP["cp"])
    coat = []
    Tsurf = tc_rule(550e3, e_wall=E_FEP)[0][0]
    for lay, dl in (("1L", 15.3e-6), ("2L", 25.8e-6), ("3L", 45.3e-6), ("4L", 64.8e-6)):
        coat.append(dict(coating=lay, thickness_um=dl * 1e6, t_semiinf_ms=dl**2 / alpha_c * 1e3,
                         T_surface_MFB_550kPa=Tsurf, T_metal_at_rewet=Tsurf + 8e3 * dl / FEP["k"],
                         bare_SS_550kPa=tc_rule(550e3)[0][0]))
    coat = pd.DataFrame(coat)
    ext = external_checks()
    # Cu property sanity print
    sanity = f"Cu NIST fits at 30/77/130/293 K: k = {cu_k(np.array([30,77,130,293.])).round(0)}, cp = {cu_cp(np.array([30,77,130,293.])).round(1)}"
    pd.set_option("display.width", 250)
    print(sanity)
    print(rt[["case", "st", "P_kPa", "Tsat", "T_mfb", "Tc_rule", "beta", "Henry", "Thn+19"]].round(2).to_string())
    print(perf.round(2).to_string()); print(dd.round(1).to_string()); print("mean up-dn offset by definition:", {k: round(v, 1) for k, v in off.items()})
    print(diag.round(2).to_string()); print(preds.round(2).to_string()); print(coat.round(2).to_string()); print(ext.round(1).to_string())
    with open(f"{OUT}/t09_rewet_theory.md", "w") as f:
        f.write(__doc__ + "\n" + sanity + "\n\n## per station\n\n")
        f.write(md_table(rt[["case", "st", "P_kPa", "Tsat", "G_lag", "T_mfb", "T_mfb_TCmin", "T_mfb_TCmax", "Tc_rule", "beta",
                             "Henry", "Thn+19"]].round(2), 4))
        f.write("\n\n## performance [K]\n\n" + md_table(perf.round(2), 3))
        f.write("\n\n## measured T_MFB under alternative definitions\n\n" + md_table(dd.reset_index().round(1), 4))
        f.write("\n\nmean up-minus-dn offset by definition [K]: " + str({k: round(v, 1) for k, v in off.items()}))
        f.write("\n\n## station differences (paired cases)\n\n" + md_table(diag.round(2), 3))
        f.write("\n\n## predictions of the T_c-contact rule\n\n" + md_table(preds.round(2), 4))
        f.write("\n\n## FEP-coated SS (Chung 2021 thicknesses; k=0.195, rho=2150, cp~700 assumed)\n\n" + md_table(coat.round(2), 4))
        f.write("\n\n## external digitized data (indicative; same 2x-rate onset definition)\n\n" + md_table(ext.round(2), 4) + "\n")
    rt.to_csv(f"{OUT}/rewet_theory_table.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for st, mk in (("up", "o"), ("dn", "s")):
        m = rt.st == st
        axs[0].errorbar(rt.P_kPa[m], rt.T_mfb[m], yerr=[rt.T_mfb[m] - rt.T_mfb_TCmin[m], rt.T_mfb_TCmax[m] - rt.T_mfb[m]],
                        fmt=mk + "k", label=f"measured low-g {st}")
    Pg = np.linspace(90e3, 800e3, 30)
    axs[0].plot(Pg / 1e3, tc_rule(Pg)[0], "C0-", label="T_c + beta(T_c - Tsat) [0 constants]")
    axs[0].plot(Pg / 1e3, T_hn_lienhard(sat(Pg)["Tsat"]) + 18.8, "C1--", label="T_hn + 19 K")
    axs[0].axhline(134, color="C2", ls=":", label="134 K")
    axs[0].axhline(TC_N2, color="gray", lw=0.6)
    for _, r in ext[ext.wall == "304SS"].iterrows():
        axs[0].plot(r.P_kPa, r.T_mfb_meas, "^", color="C3", label=r.data)
    axs[0].set_xlabel("P [kPa]"); axs[0].set_ylabel("T_w at MFB [K]"); axs[0].legend(fontsize=6); axs[0].set_ylim(110, 160)
    sysl = preds[preds.P_kPa == 300]
    axs[1].barh(sysl.system, sysl.T_wMFB - sysl.Tc, color="C0")
    axs[1].set_xlabel("T_w,MFB - T_c [K] at 300 kPa (rule)")
    for i, (_, r) in enumerate(sysl.iterrows()):
        axs[1].text(r.T_wMFB - r.Tc + 0.3, i, f"{r.T_wMFB:.1f} K (beta={r.beta:.3f})", va="center", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig15_rewet_theory.png", dpi=100); plt.close(fig)


if __name__ == "__main__":
    main()
