"""Goal 5: CHF and transition-boiling curve shape."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import D, FIG, OUT, MW_A, TC_OF, Z, md_table, sat, ss_cp
from prep import quench_info, station_data


def zuber(P, g, C=0.131):
    s = sat(P)
    return C * np.sqrt(s["rho_v"]) * s["h_fg"] * (s["sigma"] * g * (s["rho_l"] - s["rho_v"])) ** 0.25


def katto_ohno(P, G, L):
    """Katto & Ohno (1984) saturated-inlet CHF in uniformly heated tubes (rho_v/rho_l < 0.15),
    heated length L (here the distance from the tube inlet), no inlet subcooling."""
    s = sat(P)
    R = s["rho_v"] / s["rho_l"]
    l = L / D
    We = s["sigma"] * s["rho_l"] / (G**2 * L)        # Katto's (inverse) Weber number
    C = np.where(l < 50, 0.25, np.where(l < 150, 0.25 + 0.0009 * (l - 50), 0.34))
    q1 = C * We**0.043 / l
    q2 = 0.10 * R**0.133 * We ** (1 / 3) / (1 + 0.0031 * l)
    q3 = 0.098 * R**0.133 * We**0.433 * l**0.27 / (1 + 0.0031 * l)
    q4 = 0.0384 * R**0.6 * We**0.173 / (1 + 0.28 * We**0.233 * l)
    qco = np.where(q1 < q2, q1, np.where(q2 < q3, q2, q3))
    # (q4 -- the N-regime -- only enters the R > 0.15 branch of Katto-Ohno; reported for reference)
    return qco * G * s["h_fg"], dict(q1=q1 * G * s["h_fg"], q2=q2 * G * s["h_fg"], q3=q3 * G * s["h_fg"],
                                     q4=q4 * G * s["h_fg"])


def transition_shapes():
    """per-TC: compare interpolation rules between (dT_MFB, q_MFB) and (dT_CHF, q_CHF) with the
    measured lumped q'' in between (samples strictly between onset and peak)."""
    rows = []
    allpts = []
    for c in range(1, 11):
        for st in ("up", "dn"):
            s = station_data(c, st)
            qi = quench_info(c, st)
            for k in range(3):
                v = qi[f"TC{k}"]
                if not v:
                    continue
                T = s[f"TC{k}"].values; q = s[f"qTC{k}"].values; t = s.t.values
                Ts = s.Tsat.values
                j0 = np.argmin(np.abs(t - v["t_mfb"])); j1 = np.argmin(np.abs(t - v["t_q"]))
                if j1 - j0 < 2:
                    continue
                T0, T1 = T[j0], T[j1]
                q0 = MW_A * ss_cp(T0) * v["r_fb"]          # film-boiling flux at onset
                q1 = q[j1]
                idx = np.arange(j0 + 1, j1)
                dT = T[idx] - Ts[idx]
                dT0, dT1 = T0 - Ts[j0], T1 - Ts[j1]
                n = np.log(q1 / q0) / np.log(dT1 / dT0)
                q_ll = q0 * (dT / dT0) ** n
                q_lin = q0 + (q1 - q0) * (T[idx] - T0) / (T1 - T0)
                theta = (T[idx] - T0) / (T1 - T0)
                allpts.append(pd.DataFrame(dict(case=c, st=st, tc=k, theta=theta, qn=(q[idx] - q0) / (q1 - q0),
                                                q=q[idx], q_ll=q_ll, q_lin=q_lin, q0=q0, q1=q1)))
                rows.append(dict(case=c, st=st, trace=TC_OF[st][k], dT_mfb=dT0, q_mfb=q0 / 1e3, dT_chf=dT1,
                                 q_chf=q1 / 1e3, n_endpoints=n, npts=len(idx),
                                 rms_loglog=np.sqrt(np.mean(np.log(q_ll / q[idx]) ** 2)),
                                 rms_linear=np.sqrt(np.mean(np.log(q_lin / q[idx]) ** 2))))
    return pd.DataFrame(rows), pd.concat(allpts, ignore_index=True)


def main():
    st = pd.read_csv(f"{OUT}/rewet_table.csv")
    P = st.P_kPa.values * 1e3
    G = st.G_lag.values
    L = np.array([Z[s] for s in st.st])
    st["q_chf_meas"] = st.q_chf_TC
    st["Zuber_1g"] = zuber(P, 9.81) / 1e3
    st["Zuber_0.01g"] = zuber(P, 0.0981) / 1e3
    st["Zuber_1e-4g"] = zuber(P, 9.81e-4) / 1e3
    st["Katto-Ohno"] = katto_ohno(P, G, L)[0] / 1e3
    st["KO/meas"] = st["Katto-Ohno"] / st.q_chf_meas
    st["Zuber0.01g/meas"] = st["Zuber_0.01g"] / st.q_chf_meas
    # power-law fit q_CHF ~ G^m
    m, a = np.polyfit(np.log(G), np.log(st.q_chf_meas), 1)
    mP, aP = np.polyfit(np.log(P), np.log(st.q_chf_meas), 1)
    tr, pts = transition_shapes()
    # theta^c fit (universal-poster style q = q_FB + (q_CHF - q_FB) theta^c)
    cs = np.linspace(0.2, 4, 77)
    rr = [np.sqrt(np.mean((np.log(np.clip(pts.q0 + (pts.q1 - pts.q0) * pts.theta**c, 1, None)) - np.log(np.clip(pts.q, 1, None))) ** 2)) for c in cs]
    c_best = cs[int(np.argmin(rr))]
    txt = (f"q_CHF (per-TC lumped, station median) = {np.exp(a):.3g} G^{m:.2f}  [kW/m2]; vs P: exponent {mP:.2f}\n"
           f"q_CHF range {st.q_chf_meas.min():.0f}-{st.q_chf_meas.max():.0f} kW/m2, dT_CHF (TC) {st.dT_chf_TC.min():.1f}-{st.dT_chf_TC.max():.1f} K "
           f"(mean {st.dT_chf_TC.mean():.1f}); inverse-method dT_CHF {st.dT_chf_inv.min():.1f}-{st.dT_chf_inv.max():.1f} K\n"
           f"Katto-Ohno / measured: mean {st['KO/meas'].mean():.2f} (range {st['KO/meas'].min():.2f}-{st['KO/meas'].max():.2f})\n"
           f"Zuber(0.01 g)/measured: mean {st['Zuber0.01g/meas'].mean():.2f}; Zuber(1 g) {st.Zuber_1g.mean():.0f} kW/m2 mean\n"
           f"transition exponent n (q ~ dT^n between MFB and CHF, per TC): median {tr.n_endpoints.median():.2f}, "
           f"IQR {tr.n_endpoints.quantile(.25):.2f}..{tr.n_endpoints.quantile(.75):.2f}\n"
           f"interpolation RMS(ln q) inside transition: log-log {np.sqrt(np.mean(np.log(np.clip(pts.q_ll,1,None)/np.clip(pts.q,1,None))**2)):.3f}, "
           f"linear-in-T {np.sqrt(np.mean(np.log(np.clip(pts.q_lin,1,None)/np.clip(pts.q,1,None))**2)):.3f}, "
           f"theta^c best c={c_best:.2f}: {min(rr):.3f}\n")
    print(txt)
    cols = ["case", "st", "P_kPa", "G_lag", "q_chf_meas", "q_chf_TCmax", "q_chf_inv", "dT_chf_TC", "dT_chf_inv",
            "Katto-Ohno", "Zuber_1g", "Zuber_0.01g", "Zuber_1e-4g", "KO/meas"]
    with open(f"{OUT}/t05_chf_transition.md", "w") as f:
        f.write(md_table(st[cols].round(2), 4) + "\n\n```\n" + txt + "```\n\n## per-TC transition\n\n")
        f.write(md_table(tr.round(3), 3) + "\n")
    print(st[cols].round(1).to_string())
    print(tr.round(2).to_string())
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    for s_, mk in (("up", "o"), ("dn", "s")):
        sel = st.st == s_
        axs[0].plot(G[sel], st.q_chf_meas[sel], mk + "k", label=f"measured (TC lumped) {s_}")
        axs[0].plot(G[sel], st["Katto-Ohno"][sel], mk, color="C0", mfc="none", label=f"Katto-Ohno {s_}")
        axs[0].plot(G[sel], st["Zuber_0.01g"][sel], mk, color="C1", mfc="none", label=f"Zuber 0.01 g {s_}")
        axs[0].plot(G[sel], st["Zuber_1g"][sel], mk, color="C2", mfc="none", label=f"Zuber 1 g {s_}")
    axs[0].set_xlabel("G lag-corr. [kg/m2s]"); axs[0].set_ylabel("q_CHF [kW/m2]"); axs[0].legend(fontsize=6)
    axs[0].set_yscale("log")
    for _, r in tr.iterrows():
        pass
    sc = axs[1].scatter(pts.theta, pts.qn, c=pts.case, s=6, cmap="viridis")
    th = np.linspace(0, 1, 50)
    axs[1].plot(th, th, "k-", label="linear in T (Darr-Hartwig)")
    axs[1].plot(th, th**c_best, "r-", label=f"theta^{c_best:.2f}")
    axs[1].set_xlabel("theta = (T - T_MFB)/(T_CHF - T_MFB)"); axs[1].set_ylabel("(q - q_MFB)/(q_CHF - q_MFB)")
    axs[1].legend(fontsize=7); axs[1].set_ylim(-0.3, 1.3)
    for c in sorted(tr.case.unique()):
        sub = tr[tr.case == c]
        axs[2].loglog(np.r_[sub.dT_chf.values], sub.q_chf.values, "o", color=plt.cm.viridis((c - 1) / 9))
        for _, r in sub.iterrows():
            axs[2].plot([r.dT_chf, r.dT_mfb], [r.q_chf, r.q_mfb], "-", color=plt.cm.viridis((c - 1) / 9), lw=0.8)
            axs[2].plot(r.dT_mfb, r.q_mfb, "v", color=plt.cm.viridis((c - 1) / 9))
    axs[2].set_xlabel("Tw - Tsat [K]"); axs[2].set_ylabel("q'' [kW/m2]")
    axs[2].set_title(f"per-TC transition segments (CHF o -> MFB v); median n = {tr.n_endpoints.median():.2f}")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig11_chf_transition.png", dpi=100)
    plt.close(fig)
    tr.to_csv(f"{OUT}/transition_per_tc.csv", index=False)
    st.to_csv(f"{OUT}/chf_table.csv", index=False)


if __name__ == "__main__":
    main()
