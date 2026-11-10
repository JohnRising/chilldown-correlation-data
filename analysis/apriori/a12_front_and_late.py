"""Phase 3 / Tasks 2-3: (a) quench-front reset test on the downstream film h when the upstream
station rewets; (b) front-mechanism budget; (c) late-time / low-G behaviour (cases 1, 5).
Law reference: Phase-2 L0, h = 8.55 (k_vf/D) Pe_l^(1/3) (accounts for G, P, dT trends)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import D, DO, FIG, MW_A, OUT, RHO_SS, md_table, sat, smooth, ss_cp, ss_k, vap
from prep import station_data

Q = pd.read_csv(f"{OUT}/q_consistent.csv")
ST = pd.read_csv(f"{OUT}/mfb_chf_station.csv")
TR = pd.read_csv(f"{OUT}/mfb_chf_per_trace.csv")


def law(df):
    s = sat(df.P.values)
    v = vap(df.P.values, 0.5 * (df.Tw.values + s["Tsat"]))
    ul = df.G.values / s["rho_l"]
    return 8.55 * v["k"] / D * (ul * D * v["rho"] * v["cp"] / v["k"]) ** (1 / 3)


def ratio_series(c, st):
    s = Q[(Q.case == c) & (Q.st == st)].sort_values("t").reset_index(drop=True)
    s["r"] = np.log(s.h_E / law(s))
    # per-TC lumped h (with the same parasitic model)
    raw = station_data(c, st)
    for k in range(3):
        T = raw[f"TC{k}"].values[:len(s)]
        q = -MW_A * ss_cp(T) * smooth(T, 9, 2, deriv=1) + 98.6 + 2.2 * (293 - T)
        s[f"rTC{k}"] = np.log(np.clip(q / (T - s.Tsat), 1, None) / law(s.assign(Tw=T)))
    return s


def window_jump(s, t0, pre=(2.0, 0.2), post=(0.3, 2.3), tmax=None, col="r"):
    a = s[(s.t >= t0 - pre[0]) & (s.t <= t0 - pre[1])][col]
    hi = t0 + post[1] if tmax is None else min(t0 + post[1], tmax)
    b = s[(s.t >= t0 + post[0]) & (s.t <= hi)][col]
    if len(a) < 10 or len(b) < 10:
        return np.nan
    return b.mean() - a.mean()


def reset_test():
    rows, null = [], []
    for c in (3, 6, 7, 8, 9, 10):
        up = TR[(TR.case == c) & (TR.st == "up") & (TR.trace != "Ti")]
        t_up = up.t_mfb.min()                      # first upstream TC leaves film boiling
        dn = TR[(TR.case == c) & (TR.st == "dn") & (TR.trace != "Ti")]
        t_dn = dn.t_mfb.min() if len(dn) else 1e9
        s = ratio_series(c, "dn")
        tmax = min(t_dn - 0.5, s.t.max())
        row = dict(case=c, t_up_onset=t_up, t_dn_onset=t_dn if t_dn < 1e8 else np.nan,
                   jump_station=window_jump(s, t_up, tmax=tmax))
        for k in range(3):
            row[f"jump_TC{k}"] = window_jump(s, t_up, tmax=tmax, col=f"rTC{k}")
        row["G_pre"] = s[(s.t > t_up - 2) & (s.t < t_up)].G.mean()
        row["G_post"] = s[(s.t > t_up + 0.3) & (s.t < min(t_up + 2.3, tmax))].G.mean()
        rows.append(row)
        # null: same statistic at placebo times in clear film boiling (t from 4 s to t_up - 3 s)
        for tp in np.arange(4.0, t_up - 3.0, 0.5):
            null.append(window_jump(s, tp))
    df = pd.DataFrame(rows)
    null = np.array([v for v in null if np.isfinite(v)])
    return df, null


def front_budget():
    """mechanism budget per quenching case (velocities in mm/s)."""
    rows = []
    mw_per_len = RHO_SS * np.pi / 4 * (DO**2 - D**2)
    A = np.pi * D**2 / 4
    for c in (6, 7, 8, 9, 10):
        u, d_ = ST[(ST.case == c) & (ST.st == "up")].iloc[0], ST[(ST.case == c) & (ST.st == "dn")].iloc[0]
        dt = d_.t_mfb - u.t_mfb
        U = 0.252 / dt * 1e3
        s = sat(np.array([u.P_kPa * 1e3]))
        G = Q[(Q.case == c) & (Q.st == "up") & (abs(Q.t - u.t_mfb) < 1)].G.mean()
        x = Q[(Q.case == c) & (Q.st == "up") & (abs(Q.t - u.t_mfb) < 1)].x.mean()
        Tg = np.linspace(u.Tsat, u.T_mfb, 50)
        e_cool = np.trapezoid(ss_cp(Tg), Tg) * mw_per_len            # J/m to quench T_MFB -> Tsat
        U_supply = G * A * (1 - x) * s["h_fg"][0] / e_cool * 1e3
        # conduction-controlled (Yamanouchi) with wet-side h 5 kW/m2K, dry wall at T_d = T_MFB+20 K
        k = ss_k(120.0); rc = RHO_SS * ss_cp(120.0); dw = (DO - D) / 2
        th = 20.0 / (u.T_mfb - u.Tsat)
        U_cond = 1 / rc * np.sqrt(5e3 * k / dw) / np.sqrt(th * (1 + th)) * 1e3
        # kinematic (dispersed) front: dt_q = (T_dn - T_up) at t_up_mfb / cooling rate dn + dT_MFB/rate
        su = Q[(Q.case == c) & (Q.st == "up")]; sd = Q[(Q.case == c) & (Q.st == "dn")]
        j = np.argmin(abs(sd.t.values - u.t_mfb))
        Tdn = sd.Tm.values[j]
        rate = -np.polyfit(sd.t.values[j - 15:j + 1], sd.Tm.values[j - 15:j + 1], 1)[0]
        dt_kin = (Tdn - d_.T_mfb) / rate
        rows.append(dict(case=c, U_meas=U, dt_meas=dt, U_supply_limited=U_supply, U_conduction=U_cond,
                         Tdn_at_up_rewet=Tdn, T_mfb_dn=d_.T_mfb, dn_rate_K_s=rate, dt_kinematic=dt_kin,
                         part_from_hgrad=(Tdn - u.T_mfb) / rate, part_from_Tmfb_offset=(u.T_mfb - d_.T_mfb) / rate))
    return pd.DataFrame(rows)


def late_time():
    rows = []
    for c in range(1, 11):
        for st in ("up", "dn"):
            s = ratio_series(c, st)
            fb = s[s.fb & (s.t > 0.5)]
            if len(fb) < 30:
                continue
            tmid = fb.t.max() / 2
            e, l = fb[fb.t < tmid], fb[fb.t >= tmid]
            rows.append(dict(case=c, st=st, G_early=e.G.mean(), G_late=l.G.mean(),
                             Gcv_late=l.G.std() / l.G.mean(), P_late_kPa=l.P.mean() / 1e3, x_late=l.x.mean(),
                             meas_over_law_early=np.exp(e.r.mean()), meas_over_law_late=np.exp(l.r.mean()),
                             late_minus_early=l.r.mean() - e.r.mean()))
    return pd.DataFrame(rows)


def main():
    rt, null = reset_test()
    fb = front_budget()
    lt = late_time()
    nm, nsd = null.mean(), null.std()
    pd.set_option("display.width", 250)
    print(rt.round(3).to_string()); print(f"null placebo jumps: n={len(null)} mean {nm:.3f} sd {nsd:.3f}")
    print(fb.round(2).to_string()); print(lt.round(3).to_string())
    exp = dict(pure_z=np.log((0.401 / 0.252) ** (1 / 3)), L2=np.log((0.861 / 0.712) ** (1 / 3)),
               L3=np.log((0.0976 / (0.252 * 0.129 / 0.381)) ** (1 / 3)))
    print("expected ln-jumps if L resets to distance from front:", {k: round(v, 3) for k, v in exp.items()})
    with open(f"{OUT}/t12_front_late.md", "w") as f:
        f.write(__doc__ + "\n## (a) reset test: change of ln(h_meas/h_law) at the downstream station across the first "
                "upstream rewet (windows [-2,-0.2] s vs [+0.3,+2.3] s, truncated 0.5 s before downstream onset)\n\n")
        f.write(md_table(rt.round(3), 3) + f"\n\nplacebo (same statistic at times in clear film boiling): n={len(null)}, "
                f"mean {nm:.3f}, sd {nsd:.3f}\n\nexpected ln-jump if L restarts at the front: "
                + str({k: round(v, 3) for k, v in exp.items()}) + "\n\n## (b) front budget (mm/s, s)\n\n"
                + md_table(fb.round(2), 3) + "\n\n## (c) early vs late film boiling (meas/law)\n\n" + md_table(lt.round(3), 3) + "\n")
    # figure
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    for a, c in zip(axs.ravel(), (3, 6, 7, 8, 9, 10)):
        s = ratio_series(c, "dn")
        r = rt[rt.case == c].iloc[0]
        w = (s.t > r.t_up_onset - 6) & (s.t < (r.t_dn_onset if np.isfinite(r.t_dn_onset) else s.t.max()))
        a.plot(s.t[w], np.exp(s.r[w]), "k", label="station mean")
        for k in range(3):
            a.plot(s.t[w], np.exp(s[f"rTC{k}"][w]), lw=0.6)
        a.axvline(r.t_up_onset, color="r", ls="--", label="upstream rewet onset")
        a.set_ylim(0.4, 2); a.set_title(f"case {c} downstream h/law, jump {r.jump_station:+.2f}"); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig18_front_reset.png", dpi=90); plt.close(fig)
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for c, col in ((1, "C0"), (5, "C1"), (4, "C2")):
        for st, ls in (("up", "-"), ("dn", "--")):
            s = ratio_series(c, st); s = s[s.t > 0.5]
            axs[0].plot(s.t, smooth(np.exp(s.r.values), 31, 1), col, ls=ls, label=f"case {c} {st}")
        axs[1].plot(s.t, s.G, col, label=f"case {c} G (lag-corr.)")
    axs[0].axhline(1, color="k"); axs[0].set_ylabel("h_meas / h_law (2-s smoothed)"); axs[0].legend(fontsize=7)
    axs[1].set_ylabel("G"); axs[1].set_xlabel("t [s]"); axs[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig19_late_time.png", dpi=90); plt.close(fig)


if __name__ == "__main__":
    main()
