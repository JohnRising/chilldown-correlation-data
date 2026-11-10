"""Phase 4 / Tasks 4-6 (referee M1e, M4, M6): surface-renewal scales; rewet criteria comparison
(T_c rule, 134 K, Kalinin et al. 1969, Simon-Simoneau 1969, Spiegler 27/32 Tc); gravity estimates."""
import json

import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

from common import D, OUT, ROOT, md_table, sat, ss_e
from a08_film_laws import load_fb
from a09_rewet_theory import tc_rule, cu_e, E_FEP, liq, onset_generic
from a04_rewet import henry
from common import contact_wall_T, T_hn_lienhard

C3X, LAM = 15.69, 0.098


def law(df):
    ul = df.G * (1 - df.x.clip(0, 0.9)) / df.rho_l
    return C3X * df.k_f / D * (ul * D * df.rho_f * df.cp_f / df.k_f) ** (1 / 3) * (D * (df.z + LAM) / (df.z * LAM)) ** (1 / 3)


def wreg(df, y):
    w = np.sqrt(df.w.values)
    X = np.c_[np.ones(len(df)), np.log(df.G), np.log(df.P), (df.st == "dn")]
    b = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]
    return b


# ------------------------------------------------------------------ 4. renewal
def renewal(df):
    x = df.x.clip(0, 0.9)
    ul = df.G * (1 - x) / df.rho_l
    al = df.k_f / (df.rho_f * df.cp_f)
    e2 = df.k_f * df.rho_f * df.cp_f
    h = df.h_E
    dl = 0.807**3 * df.k_f**3 * ul / (al * h**3)          # delta * l_r  [m2] (Leveque restarted every l_r)
    tr = 4 * e2 / (np.pi * h**2)                          # Higbie renewal time [s]
    lr = ul * tr
    out = {}
    per = pd.DataFrame(dict(case=df.case, st=df.st, G=df.G, P=df.P, dl=dl, tr=tr, lr=lr, kh=df.k_f / h))
    agg = per.groupby(["case", "st"]).median().reset_index()
    for nm in ("dl", "tr", "lr"):
        b = wreg(df, np.log(per[nm].values))
        out[nm] = dict(median=float(np.median(per[nm])), p10=float(np.percentile(per[nm], 10)), p90=float(np.percentile(per[nm], 90)),
                       exp_G=float(b[1]), exp_P=float(b[2]))
    # candidate renewal-length closures: Higbie with t_r = l_r/u_l
    rho_l, rho_v, sig = df.rho_l, df.rho_v, df.sigma
    jv = df.G * x / rho_v
    cands = {
        "l_r = D": D + 0 * df.G,
        "l_r = lambda_KH(du = u_l)": 2 * np.pi * sig * (rho_l + rho_v) / (rho_l * rho_v * ul**2),
        "l_r = lambda_KH(du = j_v) [x endogenous]": 2 * np.pi * sig * (rho_l + rho_v) / (rho_l * rho_v * np.clip(jv, 0.05, None) ** 2),
        "l_r = lambda_KH(du = j_v - u_l)": 2 * np.pi * sig * (rho_l + rho_v) / (rho_l * rho_v * np.clip(np.abs(jv - ul), 0.05, None) ** 2),
        "l_r = capillary length-free: sigma/(rho_l u_l^2) (liquid We=1)": sig / (rho_l * ul**2),
    }
    rows = []
    for nm, l in cands.items():
        hp = 2 * np.sqrt(e2 * ul / (np.pi * l))
        lm = np.log(h / hp)
        M = float(np.exp(np.sum(df.w * lm) / df.w.sum()))
        r = np.log(M * hp / h)
        b = wreg(df, np.log(hp.values))
        rows.append(dict(renewal_length=nm, multiplier_M=M, rms=float(np.sqrt(np.sum(df.w * r**2) / df.w.sum())),
                         exp_G=b[1], exp_P=b[2], median_l_mm=float(np.median(l) * 1e3)))
    bd = wreg(df, np.log(h.values))
    rows.append(dict(renewal_length="DATA (free law)", multiplier_M=np.nan, rms=np.nan, exp_G=bd[1], exp_P=bd[2], median_l_mm=np.nan))
    return out, pd.DataFrame(rows), agg


# ------------------------------------------------------------------ 5. rewet criteria
def kalinin(P, fluid="Nitrogen", e_wall=ss_e, Tl=None):
    """Kalinin et al. 1969 (as transcribed by Darr et al. 2019): T_wet = T_sat + 1.65[0.1 + 1.5 sqrt(bw/bl)
    + 0.6 bw/bl](T_c - T_sat), beta = 1/(k rho c) -> bw/bl = (e_l/e_w)^2; e_w at T_wet (iterated)."""
    out = []
    for p in np.atleast_1d(P):
        L = liq(p, fluid)
        Ts = L["Tsat"] if Tl is None else Tl
        Tw = L["Tc"]
        for _ in range(40):
            r = (L["e"] / (e_wall(Tw) if callable(e_wall) else e_wall)) ** 2
            Tw = Ts + 1.65 * (0.1 + 1.5 * np.sqrt(r) + 0.6 * r) * (L["Tc"] - Ts)
        out.append(Tw)
    return np.array(out)


def simon_simoneau(P, fluid="Nitrogen"):
    Tc, Pc = CP.PropsSI("Tcrit", fluid), CP.PropsSI("pcrit", fluid)
    return Tc * (0.87 + 0.13 * np.atleast_1d(P) / Pc)


def rewet():
    rt = pd.read_csv(f"{OUT}/rewet_theory_table.csv")
    P = rt.P_kPa.values * 1e3
    rt["Kalinin"] = kalinin(P)
    rt["Simon-Simoneau"] = simon_simoneau(P)
    rt["Spiegler 27/32Tc"] = 27 / 32 * 126.192
    models = ["Tc_rule", "const134", "Kalinin", "Simon-Simoneau", "Spiegler 27/32Tc", "Henry", "Thn+19"]
    rows = []
    for m in models:
        e = rt[m] - rt.T_mfb
        rows.append(dict(model=m, mean_err=e.mean(), MAE=e.abs().mean(), up_mean=e[rt.st == "up"].mean(),
                         dn_mean=e[rt.st == "dn"].mean(), pred_range=f"{rt[m].min():.1f}-{rt[m].max():.1f}",
                         corr=np.corrcoef(rt[m], rt.T_mfb)[0, 1] if rt[m].std() > 1e-6 else np.nan))
    perf = pd.DataFrame(rows)
    # slope test: T_MFB = a + b Tsat + c dn, case bootstrap; compare to Henry slope
    def slope(df, col):
        X = np.c_[np.ones(len(df)), df.Tsat, (df.st == "dn")]
        return np.linalg.lstsq(X, df[col].values, rcond=None)[0][1]
    rng = np.random.default_rng(3); bs = []
    cases = rt.case.unique()
    for _ in range(2000):
        pick = rng.choice(cases, len(cases))
        sub = pd.concat([rt[rt.case == c] for c in pick])
        if sub.Tsat.std() > 0.5 and sub.st.nunique() == 2:
            bs.append(slope(sub, "T_mfb"))
    sl = dict(meas=slope(rt, "T_mfb"), ci=[np.percentile(bs, 2.5), np.percentile(bs, 97.5)],
              henry=slope(rt, "Henry"), tc_rule=slope(rt, "Tc_rule"), kalinin=slope(rt, "Kalinin"), simon=slope(rt, "Simon-Simoneau"))
    # discriminating predictions
    pr = []
    for lab, fl, ew, P in (("LN2/304SS", "Nitrogen", ss_e, 300e3), ("LN2/Cu", "Nitrogen", cu_e, 300e3),
                           ("LN2/FEP surface", "Nitrogen", E_FEP, 300e3), ("LH2/304SS", "Hydrogen", ss_e, 300e3),
                           ("LH2/Cu", "Hydrogen", cu_e, 300e3), ("LH2/Cu", "Hydrogen", cu_e, 150e3)):
        Tc = CP.PropsSI("Tcrit", fl)
        pr.append(dict(system=lab, P_kPa=P / 1e3, Tsat=liq(P, fl)["Tsat"], Tc_rule=tc_rule(P, fl, ew)[0][0],
                       const_1p061Tc=1.061 * Tc, Kalinin=kalinin(P, fl, ew)[0], Simon_Simoneau=simon_simoneau(P, fl)[0],
                       Spiegler=27 / 32 * Tc))
    pr = pd.DataFrame(pr)
    # NBS LH2 knees (indicative)
    nb = []
    meta = pd.read_csv(ROOT + "/data/nbs9264/nbs9264_runs_metadata.csv", comment="#")
    for _, m in meta[meta.fluid == "LH2"].iterrows():
        f = ROOT + f"/data/nbs9264/nbs9264_fig{int(m.figure):02d}_lh2_{str(m.P_drive_atm).replace('.', 'p')}atm.csv"
        try:
            d = pd.read_csv(f, comment="#")
        except FileNotFoundError:
            continue
        for zst in (6.1, 24.4):
            s = d[d.station_m == zst]
            if len(s) < 20 or s.T_wall_K.min() > 40:
                continue
            try:
                Tm, tq, ratio = onset_generic(s.time_s.values, s.T_wall_K.values, None, base=(10.0, 3.0), win=5)
            except Exception:
                continue
            Pz = (m.P_drive_atm + (0.82 - m.P_drive_atm) * zst / 61.0) * 101325
            nb.append(dict(fig=int(m.figure), state=m.liquid_state, z=zst, P_kPa=Pz / 1e3, T_onset=Tm, peak_ratio=ratio,
                           Tc_rule=tc_rule(Pz, "Hydrogen", cu_e)[0][0], Kalinin=kalinin(Pz, "Hydrogen", cu_e)[0],
                           Simon=simon_simoneau(Pz, "Hydrogen")[0]))
    return rt, perf, sl, pr, pd.DataFrame(nb)


# ------------------------------------------------------------------ 6. gravity
def gravity(df):
    def hb(g):
        hfg = df.h_fg + 0.4 * df.cp_f * df.dT
        return 0.62 * (df.k_f**3 * df.rho_f * (df.rho_l - df.rho_f) * g * hfg / (D * df.mu_f * df.dT)) ** 0.25
    hl = law(df)
    rows = []
    for c in range(1, 11):
        m = df.case.values == c
        row = dict(case=c, G=df.G[m].median(), h_meas=df.h_E[m].median(), h_law=hl[m].median())
        for g, lab in ((9.81, "1g"), (0.981, "0.1g"), (0.196, "0.02g")):
            b = hb(g)[m]
            row[f"hb_{lab}"] = b.median()
            row[f"cube_ratio_{lab}"] = np.median(((hl[m] ** 3 + b**3) ** (1 / 3)) / hl[m])
        rows.append(row)
    gt = pd.DataFrame(rows)
    # sensitivity: C refit excluding t < 10 s (the +-0.1 g period); Lam fixed
    r = np.log(df.h_E / (hl / C3X))
    C_all = float(np.exp(np.sum(df.w * r) / df.w.sum()))
    m = df.t.values >= 10
    w = df.w.values[m]
    C_late = float(np.exp(np.sum(w * r[m]) / w.sum()))
    per = {}
    for c in (1, 2, 3):
        mm = (df.case.values == c)
        per[c] = dict(early=float(np.exp(np.mean(r[mm & ~m]))), late=float(np.exp(np.mean(r[mm & m]))))
    return gt, dict(C_all=C_all, C_t_ge_10s=C_late, n_late=int(m.sum()), C_by_period_tests123=per)


def main():
    df = load_fb()
    ren, cand, agg = renewal(df)
    rt, perf, sl, pr, nb = rewet()
    gt, gs = gravity(df)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(json.dumps(ren, indent=0)); print(cand.round(3).to_string())
    print(perf.round(2).to_string()); print(json.dumps(sl, default=float)); print(pr.round(2).to_string()); print(nb.round(2).to_string())
    print(gt.round(2).to_string()); print(json.dumps(gs, indent=0))
    with open(f"{OUT}/t15_renewal_rewet_gravity.md", "w") as f:
        f.write(__doc__ + "\n## 4. renewal scales implied by the data (Leveque restarted every l_r: delta*l_r; Higbie: t_r, "
                "l_r = u_l t_r; u_l = G(1-x)/rho_l, film props). exp_G/exp_P = slopes of ln(quantity) on ln G, ln P\n\n```\n"
                + json.dumps(ren, indent=1) + "\n```\n\n### Higbie with a physical renewal length (h = M * 2 sqrt(k rho cp u_l/(pi l_r)))\n\n"
                + md_table(cand.round(3), 3) + "\n\n## 5. rewet criteria (11 stations)\n\n" + md_table(perf.round(2), 3)
                + "\n\nslope dT_MFB/dTsat (with station dummy): " + json.dumps(sl, default=float)
                + "\n\n### predictions that discriminate the criteria\n\n" + md_table(pr.round(2), 4)
                + "\n\n### NBS 9264 LH2 copper-line knees (hand-clicked markers; indicative only)\n\n" + md_table(nb.round(2), 3)
                + "\n\n## 6. gravity: Bromley horizontal term and cube sum vs the final law\n\n" + md_table(gt.round(2), 3)
                + "\n\n```\n" + json.dumps(gs, indent=1) + "\n```\n")


if __name__ == "__main__":
    main()
