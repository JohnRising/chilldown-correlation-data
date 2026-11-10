"""Goal 4: measured rewet (MFB) temperature vs thermodynamic-limit theories and published
correlations."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (D, FIG, OUT, TC_N2, T_hn_lienhard, T_spinodal_eos, T_spinodal_vdw, baumeister_B,
                    contact_wall_T, md_table, sat, ss_cp, ss_e, vap)
import fb_models as fm
from prep import station_data

G_EARTH = 9.81


def berenson_T(P, g):
    """Berenson (1961) minimum film boiling superheat (pool, horizontal) -- gravity dependent,
    vapour props at film T (iterated once)."""
    s = sat(P)
    Ts = s["Tsat"]
    dT = 30.0 * np.ones_like(Ts)
    for _ in range(20):
        v = vap(P, Ts + dT / 2)
        dr = s["rho_l"] - s["rho_v"]
        dT = (0.127 * v["rho"] * s["h_fg"] / v["k"] * (g * dr / (s["rho_l"] + s["rho_v"])) ** (2 / 3)
              * (s["sigma"] / (g * dr)) ** 0.5 * (v["mu"] / (g * dr)) ** (1 / 3))
    return Ts + dT


def henry(T_base, P):
    """Henry (1974): (T_MFB - T_B)/(T_B - T_l) = 0.42 [sqrt((k rho c)_l/(k rho c)_w) h_fg/(c_w (T_B - Tsat))]^0.6
    (T_l = Tsat, saturated liquid; wall properties at T_B)."""
    s = sat(P)
    Tl = s["Tsat"]
    e_l = np.sqrt(s["k_l"] * s["rho_l"] * s["cp_l"])
    X = (e_l / ss_e(T_base)) * s["h_fg"] / (ss_cp(T_base) * (T_base - Tl))
    return T_base + 0.42 * (T_base - Tl) * X**0.6


def predictions(P, G):
    s = sat(P)
    Ts = s["Tsat"]
    We = G**2 * D / (s["rho_l"] * s["sigma"])
    out = {}
    Thn = T_hn_lienhard(Ts)
    out["T_hn Lienhard"] = Thn
    out["contact: T_c = T_hn(Lienhard)"] = contact_wall_T(Thn, Ts, P)[0]
    Tvdw = T_spinodal_vdw(P)
    out["T_spinodal vdW(P)"] = Tvdw
    out["contact: T_c = T_sp,vdW"] = contact_wall_T(Tvdw, Ts, P)[0]
    Teos = T_spinodal_eos(P)
    out["T_spinodal EOS(P) [indicative]"] = Teos
    out["contact: T_c = T_sp,EOS"] = contact_wall_T(Teos, Ts, P)[0]
    out["27/32 Tc (Spiegler)"] = 27 / 32 * TC_N2 * np.ones_like(Ts)
    out["Henry(1974) on T_hn-contact base"] = henry(out["contact: T_c = T_hn(Lienhard)"], P)
    out["Berenson 1g"] = berenson_T(P, G_EARTH)
    out["Henry on Berenson 1g"] = henry(out["Berenson 1g"], P)
    out["Berenson 0.01g"] = berenson_T(P, 0.01 * G_EARTH)
    # Darr-Hartwig: simplified (LH2 2019/2020 paper eq. 18)
    out["Darr-Hartwig simpl.: 0.844Tc(1+0.060We^0.208)"] = 0.844 * TC_N2 * (1 + 0.060 * We**0.208)
    # Darr et al. 2019 (rewet paper eq. 28, modified Baumeister-Simon, surface-energy factor -> 1)
    Tguess = out["Darr-Hartwig simpl.: 0.844Tc(1+0.060We^0.208)"]
    B = baumeister_B(Tguess)
    out["Darr 2019 full (B-factor, S=1)"] = Ts + (0.844 * TC_N2 - Ts) / B * (1 + 0.060 * We**0.208)
    # De Salve & Panella (as quoted in Darr 2019 eq. 17)
    out["De Salve-Panella"] = Ts + 0.29 * (0.844 * TC_N2 - Ts) / B * (1 + 0.279 * G**0.49)
    return out, We, B


def main():
    st_df = pd.read_csv(f"{OUT}/mfb_chf_station.csv")
    tc_df = pd.read_csv(f"{OUT}/mfb_chf_per_trace.csv")
    # lag-corrected G at the MFB time
    Gl = []
    xl = []
    fbd = fm.dataset(all_samples=True)
    for _, r in st_df.iterrows():
        g = fm.g_testsection(int(r.case))
        s = station_data(int(r.case), r.st)
        j = np.argmin(np.abs(s.t.values - r.t_mfb))
        Gl.append(g[j])
        sub = fbd[(fbd.case == r.case) & (fbd.st == r.st)]
        k = np.argmin(np.abs(sub.t.values - r.t_mfb))
        xl.append(sub.x.values[max(k - 15, 0)])
    st_df["G_lag"] = Gl
    st_df["x_eq"] = xl
    P = st_df.P_kPa.values * 1e3
    pred, We, B = predictions(P, st_df.G_lag.values)
    st_df["We_D"] = We
    st_df["B_BS"] = B
    _, e_ratio = contact_wall_T(T_hn_lienhard(sat(P)["Tsat"]), sat(P)["Tsat"], P)
    st_df["e_l/e_w"] = e_ratio
    for k, v in pred.items():
        st_df[k] = v
    meas = st_df.T_mfb.values
    rows = []
    for k, v in pred.items():
        e = v - meas
        rows.append(dict(model=k, mean_err_K=e.mean(), mae_K=np.abs(e).mean(), max_abs_K=np.abs(e).max(),
                         mape_dT=np.mean(np.abs(e) / (meas - st_df.Tsat.values)) * 100,
                         corr_with_meas=np.corrcoef(v, meas)[0, 1] if np.std(v) > 0 else np.nan))
    perf = pd.DataFrame(rows)
    # trends of the data
    X = np.c_[np.ones(len(st_df)), np.log(st_df.G_lag), np.log(st_df.P_kPa), (st_df.st == "dn")]
    b = np.linalg.lstsq(X, st_df.dT_mfb.values, rcond=None)[0]
    X2 = np.c_[np.ones(len(st_df)), st_df.Tsat, (st_df.st == "dn")]
    b2 = np.linalg.lstsq(X2, st_df.T_mfb.values, rcond=None)[0]
    r2 = 1 - np.sum((st_df.T_mfb.values - X2 @ b2) ** 2) / np.sum((meas - meas.mean()) ** 2)
    up = st_df[st_df.st == "up"].set_index("case").T_mfb
    dn = st_df[st_df.st == "dn"].set_index("case").T_mfb
    both = up.index.intersection(dn.index)
    trend_txt = (f"dT_MFB = {b[0]:.1f} + {b[1]:.2f} ln G + {b[2]:.2f} ln P[kPa] + {b[3]:.1f}*dn\n"
                 f"T_MFB = {b2[0]:.1f} + {b2[1]:.2f} Tsat + {b2[2]:.1f}*dn  (R2={r2:.2f})\n"
                 f"up - dn (cases {list(both)}): {(up[both]-dn[both]).round(1).tolist()} K, mean {(up[both]-dn[both]).mean():.1f} K\n"
                 f"all-trace T_MFB (2x): mean {tc_df[tc_df.trace!='Ti'].T_mfb.mean():.1f} K, sd {tc_df[tc_df.trace!='Ti'].T_mfb.std():.1f} K;"
                 f" station medians: mean {meas.mean():.1f} K, sd {meas.std(ddof=1):.1f} K\n")
    # candidate one-constant low-g rewet rules: which normalisation is most invariant?
    s_ = sat(P)
    e_l = np.sqrt(s_["k_l"] * s_["rho_l"] * s_["cp_l"])
    Tcont = (ss_e(meas) * meas + e_l * st_df.Tsat.values) / (ss_e(meas) + e_l)
    st_df["T_contact_meas"] = Tcont
    norms = {"T_MFB [K]": meas, "T_MFB - Tsat [K]": meas - st_df.Tsat.values,
             "T_MFB - T_hn [K]": meas - pred["T_hn Lienhard"], "T_contact at MFB [K]": Tcont,
             "T_MFB / Tc": meas / TC_N2, "(T_MFB-Tsat)/(T_hn-Tsat)": (meas - st_df.Tsat.values) / (pred["T_hn Lienhard"] - st_df.Tsat.values)}
    nrows = []
    for k, v in norms.items():
        # LOCO error in T_MFB when the rule's constant is the mean of the other cases
        errs = []
        for c in st_df.case.unique():
            tr = st_df.case.values != c
            cst = v[tr].mean()
            # back out T_MFB from the rule
            if k == "T_MFB [K]": Tp = cst * np.ones((~tr).sum())
            elif k == "T_MFB - Tsat [K]": Tp = st_df.Tsat.values[~tr] + cst
            elif k == "T_MFB - T_hn [K]": Tp = pred["T_hn Lienhard"][~tr] + cst
            elif k == "T_contact at MFB [K]": Tp = contact_wall_T(cst * np.ones((~tr).sum()), st_df.Tsat.values[~tr], P[~tr])[0]
            elif k == "T_MFB / Tc": Tp = cst * TC_N2 * np.ones((~tr).sum())
            else: Tp = st_df.Tsat.values[~tr] + cst * (pred["T_hn Lienhard"][~tr] - st_df.Tsat.values[~tr])
            errs.append(Tp - meas[~tr])
        errs = np.concatenate(errs)
        nrows.append(dict(rule=k, mean=v.mean(), sd=v.std(ddof=1), loco_rmse_K=np.sqrt(np.mean(errs**2))))
    norm_df = pd.DataFrame(nrows)
    print(norm_df.round(3).to_string())
    # quench-front velocity between stations
    qf = []
    for c in st_df.case.unique():
        u = st_df[(st_df.case == c) & (st_df.st == "up")]
        d_ = st_df[(st_df.case == c) & (st_df.st == "dn")]
        if len(u) and len(d_):
            dt_ = d_.t_mfb.values[0] - u.t_mfb.values[0]
            qf.append(dict(case=c, dt_s=dt_, U_qf_mm_s=(0.401 - 0.149) / dt_ * 1e3,
                           u_l_mm_s=u.G_lag.values[0] / sat(u.P_kPa.values * 1e3)["rho_l"][0] * 1e3))
    qf = pd.DataFrame(qf)
    print(qf.round(2).to_string())
    trend_txt += "\none-constant rewet rules (LOCO by case):\n" + norm_df.round(3).to_string() + "\n"
    trend_txt += "\nquench-front velocity between stations (2x-onset times):\n" + qf.round(2).to_string() + "\n"
    print(trend_txt)
    cols = ["case", "st", "P_kPa", "Tsat", "G_lag", "We_D", "x_eq", "T_mfb", "T_mfb_TCmin", "T_mfb_TCmax",
            "T_mfb_unc", "dT_mfb"] + list(pred.keys())
    with open(f"{OUT}/t04_rewet.md", "w") as f:
        f.write("## measured vs predicted T_MFB [K] per station\n\n")
        f.write(md_table(st_df[cols].round(2), 4) + "\n\n## model performance (K)\n\n")
        f.write(md_table(perf.round(2), 3) + "\n\n## data trends\n\n```\n" + trend_txt + "```\n")
        f.write(f"\ncontact-temperature amplification e_l/e_w (LN2 on 304SS at T_hn): "
                f"{st_df['e_l/e_w'].min():.3f}-{st_df['e_l/e_w'].max():.3f}; Baumeister B: {B.min():.3f}-{B.max():.3f}\n")
    st_df.to_csv(f"{OUT}/rewet_table.csv", index=False)
    pd.set_option("display.width", 250)
    print(st_df[cols].round(1).to_string())
    print(perf.round(2).to_string())

    # figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.5))
    show = ["contact: T_c = T_hn(Lienhard)", "contact: T_c = T_sp,EOS", "Henry(1974) on T_hn-contact base",
            "Darr-Hartwig simpl.: 0.844Tc(1+0.060We^0.208)", "Darr 2019 full (B-factor, S=1)", "De Salve-Panella",
            "Berenson 1g"]
    mk = dict(up="o", dn="s")
    for k, m in zip(show, ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]):
        for st in ("up", "dn"):
            sel = st_df.st == st
            axs[0].plot(meas[sel], st_df[k][sel], mk[st], color=m, label=k if st == "up" else None)
    axs[0].errorbar(meas, meas, xerr=[meas - st_df.T_mfb_TCmin, st_df.T_mfb_TCmax - meas], fmt="none", ecolor="gray")
    axs[0].plot([100, 175], [100, 175], "k-", lw=0.8)
    axs[0].set_xlabel("measured T_MFB [K] (median of TCs; bars = TC range)")
    axs[0].set_ylabel("predicted [K]"); axs[0].legend(fontsize=7)
    axs[0].set_xlim(105, 150); axs[0].set_ylim(100, 175)
    Pg = np.linspace(90e3, 800e3, 60)
    pr, _, _ = predictions(Pg, 300 * np.ones_like(Pg))
    for k, m in zip(show[:5], ["C0", "C1", "C2", "C3", "C4"]):
        axs[1].plot(Pg / 1e3, pr[k], "-", color=m, label=k + (" (G=300)" if "We" in k or "Darr" in k else ""))
    axs[1].plot(Pg / 1e3, sat(Pg)["Tsat"], "b:", label="Tsat")
    for st in ("up", "dn"):
        sel = st_df.st == st
        axs[1].errorbar(st_df.P_kPa[sel], meas[sel], yerr=[meas[sel] - st_df.T_mfb_TCmin[sel],
                        st_df.T_mfb_TCmax[sel] - meas[sel]], fmt=mk[st] + "k", label=f"measured {st}")
    tcs = tc_df[tc_df.trace != "Ti"].merge(st_df[["case", "st", "P_kPa", "G_lag"]], on=["case", "st"])
    axs[1].plot(tcs.P_kPa, tcs.T_mfb, "x", color="gray", ms=4, label="per TC")
    axs[1].set_xlabel("P [kPa]"); axs[1].set_ylabel("T_MFB [K]"); axs[1].legend(fontsize=6)
    for st in ("up", "dn"):
        sel = st_df.st == st
        axs[2].plot(st_df.G_lag[sel], st_df.dT_mfb[sel], mk[st] + "k", label=f"measured {st}")
        for _, r in st_df[sel].iterrows():
            axs[2].annotate(str(int(r.case)), (r.G_lag, r.dT_mfb), fontsize=7)
    axs[2].plot(tcs.G_lag, tcs.T_mfb - tcs.Tsat, "x", color="gray", ms=4, label="per TC")
    axs[2].set_xlabel("G lag-corrected at MFB [kg/m2s]"); axs[2].set_ylabel("T_MFB - Tsat [K]")
    axs[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig10_rewet.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
