"""Phase 5: rewet (knee) temperatures from the digitized NBS 9264 copper-line histories.

The digitized curves are piecewise linear between hand-clicked markers (resampled at 0.6 s), so the
analysis works on the linear SEGMENTS between kinks.  Because Cu cp(T) falls steeply at low T, a jump in
dT/dt is not a jump in heat flux; we use the lumped wall flux per inner area
    q = m''_w cp_Cu(T) (-dT/dt),  m''_w = rho (Do^2 - Di^2)/(4 Di) = 15.25 kg/m2 (OFHC, 19.0/15.9 mm)
Knee (MFB onset) = upper temperature of the first segment, going back from the max-q segment, whose q
exceeds 2x the film-boiling flux q_fb (median q of segments 15-60 K (N2) / 8-40 K (H2) above the peak
segment).  Resolution: the knee lies between T_knee and T_prev (upper end of the preceding segment).
Local P ~ linear between P_drive and 0.82 atm along 61 m (rough); TC accuracy is poor at low T (report
sec. 3.1; Cu-constantan with LN2 reference -> LH2 readings are least reliable)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

from common import FIG, OUT, ROOT, md_table
from a09_rewet_theory import cu_cp

MW = 8960 * (0.0190**2 - 0.0159**2) / (4 * 0.0159)


def segments(t, T, tol=0.05):
    d = np.diff(T) / np.diff(t)
    brk = np.where(np.abs(np.diff(d)) > tol)[0] + 1
    idx = np.r_[0, brk, len(T) - 1]
    segs = []
    for a, b in zip(idx[:-1], idx[1:]):
        if b - a < 1:
            continue
        sl = (T[b] - T[a]) / (t[b] - t[a])
        Tm = 0.5 * (T[a] + T[b])
        segs.append(dict(t0=t[a], t1=t[b], T0=T[a], T1=T[b], q=MW * cu_cp(Tm) * (-sl)))
    return pd.DataFrame(segs)


def knee(t, T, fluid):
    S = segments(t, T)
    lim, band = (170.0, (15, 60)) if fluid == "LN2" else (90.0, (8, 40))
    cand = S[(S[["T0", "T1"]].mean(axis=1) < lim) & (S.q > 0)]
    if len(cand) < 3:
        return None
    k = cand.q.idxmax()
    Tp = 0.5 * (S.T0[k] + S.T1[k])
    fbm = (S.T1 > Tp + band[0]) & (S.T0 < Tp + band[1]) & (S.index < k)
    if fbm.sum() < 2:
        return None
    qfb = S.q[fbm].median()
    j = k
    while j - 1 >= 0 and S.q[j - 1] > 2 * qfb:
        j -= 1
    return dict(T_knee=S.T0[j], T_prev=S.T0[j - 1] if j > 0 else np.nan, T_peak_seg=Tp, q_fb=qfb, q_peak=S.q[k],
                ratio=S.q[k] / qfb, T_end=S.T1[k], n_seg=len(S), dT_seg_at_knee=S.T0[j] - S.T1[j])


def main():
    meta = pd.read_csv(ROOT + "/data/nbs9264/nbs9264_runs_metadata.csv", comment="#")
    rows = []
    fig, axs = plt.subplots(3, 4, figsize=(20, 12))
    for a, (_, m) in zip(axs.ravel(), meta.iterrows()):
        fl = m.fluid
        f = ROOT + f"/data/nbs9264/nbs9264_fig{int(m.figure):02d}_{fl.lower()}_{str(m.P_drive_atm).replace('.', 'p')}atm.csv"
        d = pd.read_csv(f, comment="#")
        cf = "Nitrogen" if fl == "LN2" else "Hydrogen"
        for z, g in d.groupby("station_m"):
            t, T = g.time_s.values, g.T_wall_K.values
            S = segments(t, T)
            Tm = 0.5 * (S.T0 + S.T1)
            a.semilogy(Tm, S.q.clip(lower=10), ".-", ms=3, label=f"z={z}")
            r = knee(t, T, fl)
            Pz = (m.P_drive_atm + (0.82 - m.P_drive_atm) * z / 61.0) * 101325
            Ts = CP.PropsSI("T", "P", Pz, "Q", 0, cf)
            row = dict(fig=int(m.figure), fluid=fl, state=m.liquid_state, P_drive_atm=m.P_drive_atm, z_m=z,
                       P_local_kPa=Pz / 1e3, Tsat_local=Ts, T_plateau=np.median(T[-5:]))
            if r:
                row.update(r)
                a.axvline(r["T_knee"], color=a.lines[-1].get_color(), ls=":")
            rows.append(row)
        a.set_title(f"fig {int(m.figure)} {fl} {m.liquid_state} {m.P_drive_atm} atm", fontsize=9)
        a.set_xlabel("T_wall [K]"); a.set_ylabel("q lumped [W/m2]"); a.legend(fontsize=6)
        a.set_xlim((70, 250) if fl == "LN2" else (15, 150))
    fig.tight_layout(); fig.savefig(f"{FIG}/fig20_nbs_knees.png", dpi=70); plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/nbs_knees.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.round(2).to_string())
    with open(f"{OUT}/t17_nbs_knees.md", "w") as fh:
        fh.write(__doc__ + "\n" + md_table(df.round(2), 4) + "\n")


if __name__ == "__main__":
    main()


def compare_rules():
    """LN2 knees with a clear flux jump (ratio >= 1.9, segment resolution <= 20 K, not the exit station)
    against the rewet-rule family (copper wall, local P); LH2: end-of-quench corner temperatures."""
    import json
    from a16_rewet_family import rule, RULES
    from a09_rewet_theory import cu_e
    from a15_renewal_rewet_gravity import kalinin, simon_simoneau
    df = pd.read_csv(f"{OUT}/nbs_knees.csv")
    par = json.load(open(f"{OUT}/rewet_family_params.json"))
    ln = df[(df.fluid == "LN2") & (df.ratio >= 1.9) & (df.dT_seg_at_knee <= 20) & (df.z_m < 60)].copy()
    rows = []
    for nm in RULES:
        p = np.nan if par[nm] is None else par[nm]
        ln[nm] = [rule(nm, p, P * 1e3, "Nitrogen", cu_e, cu_cp)[0][0] for P in ln.P_local_kPa]
    ln["Kalinin"] = [kalinin(P * 1e3, "Nitrogen", cu_e)[0] for P in ln.P_local_kPa]
    ln["Simon-Simoneau"] = simon_simoneau(ln.P_local_kPa.values * 1e3)
    for nm in list(RULES) + ["Kalinin", "Simon-Simoneau"]:
        e = ln[nm] - ln.T_knee
        rows.append(dict(rule=nm, mean_err=e.mean(), MAE=e.abs().mean(), n=len(e)))
    X = np.c_[np.ones(len(ln)), ln.Tsat_local]
    b = np.linalg.lstsq(X, ln.T_knee.values, rcond=None)[0]
    perf = pd.DataFrame(rows)
    # LH2: corner = lower end of the steepest segment below 120 K (end of quench), plateau reading
    lh = []
    meta = pd.read_csv(ROOT + "/data/nbs9264/nbs9264_runs_metadata.csv", comment="#")
    for _, m in meta[meta.fluid == "LH2"].iterrows():
        f = ROOT + f"/data/nbs9264/nbs9264_fig{int(m.figure):02d}_lh2_{str(m.P_drive_atm).replace('.', 'p')}atm.csv"
        d = pd.read_csv(f, comment="#")
        for z, g in d.groupby("station_m"):
            S = segments(g.time_s.values, g.T_wall_K.values)
            S["rate"] = -(S.T1 - S.T0) / (S.t1 - S.t0)
            c = S[(S.T0 < 120)]
            if len(c) < 3:
                continue
            k = c.rate.idxmax()
            lh.append(dict(fig=int(m.figure), z_m=z, corner_T=S.T1[k], steep_rate_K_s=S.rate[k],
                           seg_dt_s=S.t1[k] - S.t0[k], plateau=np.median(g.T_wall_K.values[-5:])))
    lh = pd.DataFrame(lh)
    print(ln[["fig", "z_m", "P_local_kPa", "Tsat_local", "T_knee", "T_prev", "ratio", "R0", "R5", "R4", "Kalinin"]].round(1).to_string())
    print(perf.round(2).to_string()); print("LN2 knee slope vs Tsat:", b[1]); print(lh.round(1).to_string())
    with open(f"{OUT}/t17_nbs_knees.md", "a") as fh:
        fh.write("\n## LN2 knees (ratio >= 1.9, resolution <= 20 K, z < 60 m) vs rules (Cu wall, local P)\n\n"
                 + md_table(ln.round(2), 4) + "\n\n" + md_table(perf.round(2), 3)
                 + f"\n\nLN2 knee slope vs local Tsat: {b[1]:.2f} K/K (n={len(ln)}; local P estimate rough)\n\n"
                 "## LH2: end-of-quench corner (lower end of the steepest segment below 120 K) and plateau readings\n\n"
                 + md_table(lh.round(2), 3) + "\n")


if __name__ == "__main__":
    compare_rules()
