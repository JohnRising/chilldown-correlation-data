"""Goal 1: data characterisation -- regimes, boiling curves, MFB/CHF table, per-TC quench,
inverse-vs-lumped heat-flux consistency, case-7 inlet probe."""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (CASES, FIG, MW_A, OUT, STATIONS, TABLE1, TC_OF, Z, load, md_table, smooth,
                    ss_cp)
from prep import film_mask, quench_info, station_data

COL = plt.cm.viridis(np.linspace(0, 0.95, 10))
TCCOL = ["tab:blue", "tab:orange", "tab:green"]


def overview():
    fig, axs = plt.subplots(10, 4, figsize=(17, 30))
    for i, c in enumerate(CASES):
        for j, st in enumerate(STATIONS):
            s = station_data(c, st)
            qi = quench_info(c, st)
            ax = axs[i, 2 * j]
            for k in range(3):
                ax.plot(s.t, s[f"TC{k}"], color=TCCOL[k], lw=0.9, label=TC_OF[st][k])
                if qi[f"TC{k}"]:
                    v = qi[f"TC{k}"]
                    ax.plot(v["t_mfb"], v["T_mfb"], "v", color=TCCOL[k], mec="k")
                    ax.plot(v["t_q"], v["T_chf"], "o", color=TCCOL[k], mec="k", ms=4)
            ax.plot(s.t, s.Ti, "k--", lw=0.8, label="Ti (inverse)")
            ax.plot(s.t, s.Tsat, "b:", lw=1, label="Tsat")
            ax.set_ylabel("T [K]")
            ax.set_title(f"case {c} {st} (z={Z[st]} m)")
            if i == 0:
                ax.legend(fontsize=7)
            ax = axs[i, 2 * j + 1]
            ax.plot(s.t, s.q / 1e3, "r", lw=1, label="q'' inverse")
            for k in range(3):
                ax.plot(s.t, s[f"qTC{k}"] / 1e3, color=TCCOL[k], lw=0.5, alpha=0.7)
            ax.set_ylabel("q'' [kW/m2]")
            ax.set_ylim(-5, max(20, np.nanmax(s[[f"qTC{k}" for k in range(3)]].values) / 1e3 * 1.05))
            a2 = ax.twinx()
            a2.plot(s.t, s.G, "g", lw=0.8)
            a2.plot(s.t, s.P / 1e3 / 2, "m", lw=0.6)
            a2.set_ylabel("G [kg/m2s] (green), P/2 [kPa] (mag.)", fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("t [s]")
    fig.suptitle("Darr et al. 2016 low-g LN2 chilldown: wall T per TC (markers: v = MFB onset, o = CHF), "
                 "inverse q'' (red) and per-TC lumped q'' (thin)", y=1.0)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig01_overview.png", dpi=80)
    plt.close(fig)


def boiling_curves():
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    for j, st in enumerate(STATIONS):
        ax = axs[j]
        for i, c in enumerate(CASES):
            s = station_data(c, st)
            m = s.t > 0.5
            ax.loglog(s.dT[m], np.clip(s.q_ec[m], 50, None) / 1e3, "-", color=COL[i], lw=1, label=f"case {c}")
            qi = quench_info(c, st)
            for k in range(3):
                v = qi[f"TC{k}"]
                if v:
                    w = (s.t > v["t_q"] - 2.5) & (s.t < v["t_q"] + 2.0)
                    ax.loglog((s[f"TC{k}"] - s.Tsat)[w], np.clip(s[f"qTC{k}"][w], 50, None) / 1e3, ":",
                              color=COL[i], lw=0.9)
        ax.set_xlabel("Ti - Tsat [K]")
        ax.set_ylabel("q'' [kW/m2]")
        ax.set_title(f"{st}: solid = energy-consistent inverse q''(Ti); dotted = per-TC lumped q'' (quench)", fontsize=9)
        ax.set_ylim(0.3, 200)
        ax.set_xlim(1, 250)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig02_boiling_curves.png", dpi=110)
    plt.close(fig)


def htc_plots():
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    for j, st in enumerate(STATIONS):
        for i, c in enumerate(CASES):
            s = station_data(c, st)
            m = film_mask(c, st)
            axs[0, j].plot(s.t[m], s.h_ec[m], color=COL[i], lw=1, label=f"case {c}")
            axs[1, j].loglog(s.dT[m], s.h_ec[m], ".", ms=2, color=COL[i], label=f"case {c}")
        axs[0, j].set_xlabel("t [s]"); axs[0, j].set_ylabel("h [W/m2K]")
        axs[0, j].set_title(f"{st}: film-boiling h (energy-consistent) vs time")
        axs[1, j].set_xlabel("Ti - Tsat [K]"); axs[1, j].set_ylabel("h [W/m2K]")
        axs[1, j].set_title(f"{st}: film-boiling h vs superheat")
        axs[1, j].grid(True, which="both", alpha=0.3)
        axs[0, j].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig03_htc_film.png", dpi=110)
    plt.close(fig)


def case7_inlet():
    d = load(7)
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax[0].plot(d.t_s, d.T_inlet_K, label="T_inlet probe")
    ax[0].plot(d.t_s, d.Tsat_inlet_K, label="Tsat(P_inlet)")
    ax[0].plot(d.t_s, d.Tsat_up_K, label="Tsat(P1) upstream station")
    ax[0].set_ylabel("T [K]"); ax[0].legend()
    ax[1].plot(d.t_s, d.T_inlet_K - d.Tsat_inlet_K, label="T_inlet - Tsat(P_inlet)")
    ax[1].plot(d.t_s, d.T_inlet_K - d.Tsat_up_K, label="T_inlet - Tsat(P1)")
    ax[1].plot(d.t_s, (d.P_inlet_kPa - d.P1_kPa) / 10, label="(P_inlet - P1)/10 [kPa]")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_xlabel("t [s]"); ax[1].legend()
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig04_case7_inlet.png", dpi=110)
    plt.close(fig)
    t = d.t_s.values
    rows = []
    for a, b in [(0, 2), (2, 10), (10, 20), (20, 28.3)]:
        m = (t >= a) & (t < b)
        rows.append(dict(window=f"{a}-{b} s", T_inlet=d.T_inlet_K[m].mean(),
                         Tsat_inlet=d.Tsat_inlet_K[m].mean(), sup_inlet=(d.T_inlet_K - d.Tsat_inlet_K)[m].mean(),
                         sup_vs_P1=(d.T_inlet_K - d.Tsat_up_K)[m].mean(),
                         Pin_minus_P1_kPa=(d.P_inlet_kPa - d.P1_kPa)[m].mean()))
    return pd.DataFrame(rows)


def consistency():
    """inverse q'' vs lumped wall-storage q'' in film boiling."""
    rows = []
    for c in CASES:
        for st in STATIONS:
            s = station_data(c, st)
            m = film_mask(c, st) & (s.t > 3)
            dq = s.q[m] - s.q_lump[m]
            rows.append(dict(case=c, st=st, q_inv_med=s.q[m].median() / 1e3, q_lump_med=s.q_lump[m].median() / 1e3,
                             ratio=(s.q[m].sum() / s.q_lump[m].sum()), offset_kW=dq.median() / 1e3,
                             k_inv=s.attrs["k_inv"], ratio_ec=(s.q_ec[m].sum() / s.q_lump[m].sum())))
    return pd.DataFrame(rows)


def mfb_chf_table():
    rows, tcrows = [], []
    for c in CASES:
        for st in STATIONS:
            s = station_data(c, st)
            qi = quench_info(c, st)
            tcs = {k: v for k, v in qi.items() if k.startswith("TC") and v}
            if not tcs:
                continue
            for k, v in qi.items():
                if not v:
                    continue
                j = np.argmin(np.abs(s.t.values - v["t_mfb"]))
                jq = np.argmin(np.abs(s.t.values - v["t_q"]))
                Tts = s.Tsat.values[j]
                if k.startswith("TC"):
                    qchf = MW_A * ss_cp(v["T_chf"]) * v["r_max"]
                    qfb = MW_A * ss_cp(v["T_mfb"]) * v["r_fb"]
                else:
                    qchf = np.nan; qfb = np.nan
                # transition-curve log-log slope between onset and peak (per TC)
                slope = np.nan
                if k.startswith("TC"):
                    w = (s.t >= v["t_mfb"]) & (s.t <= v["t_q"])
                    x = np.log((s[k] - s.Tsat)[w]); y = np.log(np.clip(s[f"q{k}"][w], 1, None))
                    if w.sum() >= 3:
                        slope = np.polyfit(x, y, 1)[0]
                tcrows.append(dict(case=c, st=st, trace=k if k == "Ti" else TC_OF[st][int(k[2])],
                                   t_mfb=v["t_mfb"], T_mfb=v["T_mfb"], T_mfb_15=v["T_mfb_1.5x"],
                                   T_mfb_3=v["T_mfb_3x"], T_mfb_curv=v["T_mfb_curv"], Tsat=Tts,
                                   t_chf=v["t_q"], T_chf=v["T_chf"], dT_chf=v["T_chf"] - s.Tsat.values[jq],
                                   q_chf_lump=qchf / 1e3, q_fb_lump=qfb / 1e3, r_fb=v["r_fb"], r_max=v["r_max"],
                                   tb_slope=slope))
            T = np.array([v["T_mfb"] for v in tcs.values()])
            Tall = np.array([[v[f] for f in ("T_mfb_1.5x", "T_mfb_2x", "T_mfb_3x", "T_mfb_curv")] for v in tcs.values()])
            t0 = min(v["t_mfb"] for v in tcs.values())
            t1 = max(v["t_q"] for v in tcs.values())
            j = np.argmin(np.abs(s.t.values - np.median([v["t_mfb"] for v in tcs.values()])))
            pre = (s.t >= t0 - 1.0) & (s.t < t0)
            win = (s.t >= t0 - 0.5) & (s.t <= t1 + 1.5)
            jq = np.where(win)[0][np.argmax(s.q_ec.values[win])]
            qchf_tc = [MW_A * ss_cp(v["T_chf"]) * v["r_max"] for v in tcs.values()]
            dTchf_tc = [v["T_chf"] - s.Tsat.values[np.argmin(np.abs(s.t.values - v["t_q"]))] for v in tcs.values()]
            # nucleate tail: after the last TC quench + 1 s
            nb = (s.t > t1 + 1.0) & (s.dT > 1) & (s.dT < 15)
            rows.append(dict(
                case=c, st=st, nTC_quench=len(tcs), t_mfb=np.median([v["t_mfb"] for v in tcs.values()]),
                P_kPa=s.P.values[j] / 1e3, G_meas=s.G.values[j], G_tab=TABLE1[c][2], Tsat=s.Tsat.values[j],
                T_mfb=np.median(T), T_mfb_TCmin=T.min(), T_mfb_TCmax=T.max(),
                T_mfb_unc=0.5 * (Tall.max() - Tall.min()),
                dT_mfb=np.median(T) - s.Tsat.values[j],
                q_mfb=s.q_ec[pre].median() / 1e3,
                q_chf_inv=s.q_ec.values[jq] / 1e3, q_chf_inv_raw=s.q.values[jq] / 1e3, dT_chf_inv=s.dT.values[jq],
                q_chf_TC=np.median(qchf_tc) / 1e3, q_chf_TCmax=np.max(qchf_tc) / 1e3,
                dT_chf_TC=np.median(dTchf_tc),
                h_nb=(s.q_ec[nb] / s.dT[nb]).median() if nb.sum() > 5 else np.nan,
                dT_nb=s.dT[nb].median() if nb.sum() > 5 else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(tcrows)


def quench_zoom(tcdf):
    sel = tcdf[["case", "st"]].drop_duplicates().values
    n = len(sel)
    fig, axs = plt.subplots(int(np.ceil(n / 3)), 3, figsize=(16, 3.2 * np.ceil(n / 3)))
    axs = axs.ravel()
    for a, (c, st) in zip(axs, sel):
        s = station_data(c, st)
        qi = quench_info(c, st)
        tq = [v["t_q"] for v in qi.values() if v]
        w = (s.t > min(tq) - 3) & (s.t < max(tq) + 3)
        for k in range(3):
            a.plot(s.t[w], s[f"TC{k}"][w], color=TCCOL[k], lw=1)
            v = qi[f"TC{k}"]
            if v:
                for f, mk in (("T_mfb_1.5x", "^"), ("T_mfb_2x", "v"), ("T_mfb_curv", "s")):
                    tt = v["t_mfb_curv"] if f == "T_mfb_curv" else v[f"t_mfb_{f.split('_')[-1]}"]
                    a.plot(tt, v[f], mk, color=TCCOL[k], mec="k", ms=5)
        a.plot(s.t[w], s.Ti[w], "k--", lw=0.8)
        a.plot(s.t[w], s.Tsat[w], "b:")
        a2 = a.twinx()
        a2.plot(s.t[w], s.q[w] / 1e3, "r", lw=0.8)
        for k in range(3):
            a2.plot(s.t[w], s[f"qTC{k}"][w] / 1e3, color=TCCOL[k], lw=0.5, ls="--")
        a2.set_ylabel("q'' kW/m2", fontsize=7)
        a.set_title(f"case {c} {st}: ^1.5x v2x scurv", fontsize=9)
    for a in axs[n:]:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig05_quench_zoom.png", dpi=90)
    plt.close(fig)


if __name__ == "__main__":
    overview()
    boiling_curves()
    htc_plots()
    inl = case7_inlet()
    con = consistency()
    st_df, tc_df = mfb_chf_table()
    quench_zoom(tc_df)
    st_df.to_csv(f"{OUT}/mfb_chf_station.csv", index=False)
    tc_df.to_csv(f"{OUT}/mfb_chf_per_trace.csv", index=False)
    con.to_csv(f"{OUT}/q_consistency.csv", index=False)
    with open(f"{OUT}/t01_characterization.md", "w") as f:
        f.write("## MFB / CHF per case and station (per-TC medians; T_MFB = 2x-cooling-rate onset)\n\n")
        f.write(md_table(st_df.round(3), 4) + "\n\n")
        f.write("## Per trace (TC and inverse Ti)\n\n")
        f.write(md_table(tc_df.round(3), 4) + "\n\n")
        f.write("## Inverse vs lumped q'' (film boiling, t>3 s)\n\n")
        f.write(md_table(con.round(3), 3) + "\n\n")
        f.write("## Case 7 inlet probe\n\n")
        f.write(md_table(inl.round(3), 4) + "\n")
    pd.set_option("display.width", 250)
    print(st_df.round(2).to_string())
    print(tc_df.round(2).to_string())
    print(con.round(3).to_string())
    print(inl.round(2).to_string())
