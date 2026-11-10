"""Regenerate the paper figures from results/ and analysis/apriori/out/.

  cd chilldown-lowg && .venv/bin/python analysis/make_paper_figs.py [F1 F2 ...]

Writes paper/figures/<name>.pdf and .png (serif fonts; widths 3.5 in single / 7 in double column).
Inputs:
  analysis/apriori/out/film_boiling_dataset.csv, rewet_table.csv, q_consistent.csv   (a-priori data)
  results/final11/{uf,sens,num,val} (ofchill.final11), results/final11/chilldown_stations.csv (ofchill.chilldown),
  results/final/uf (old-rewet comparison rows: Miropolskii, DB-vapor), results/nbs6, nbs7 (NBS runs)
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ofchill.data import UF_STATION_TCS, load_uf_case  # noqa: E402
from ofchill.metrics import cross  # noqa: E402
from ofchill.runner import load_run  # noqa: E402

OUT = ROOT / "paper" / "figures"
APR = ROOT / "analysis" / "apriori" / "out"
RES = ROOT / "results"
D = 0.01168
W1, W2 = 3.5, 7.0
COL = {"data": "k", "spec11": "#0072B2", "spec10": "#0072B2", "spec": "#56B4E9", "miropolskii": "#D55E00", "db_vapor": "#009E73",
       "ds": "#CC79A7", "d119s": "#E69F00", "g": "#0072B2", "specz11": "#56B4E9", "ds11": "#CC79A7",
       "d119s11": "#E69F00"}
LAB = {"spec11": "this closure", "spec10": "this closure", "spec": r"film law with $L_{eff}=D$",
       "miropolskii": "Miropolskii", "db_vapor": "vapor Dittus–Boelter, full $G$", "ds": "Darr et al. (2016)",
       "d119s": "Darr et al. (2016) ×1.19", "specz11": r"film law with $L_{eff}=D$",
       "ds11": "Darr et al. (2016) film law", "d119s11": "Darr et al. (2016) film law ×1.19"}
F11 = "final11"


def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman", "Times"], "mathtext.fontset": "dejavuserif",
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "lines.linewidth": 1.0, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.direction": "in", "ytick.direction": "in", "legend.frameon": False,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "figure.dpi": 150})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", name)


def run(tag, rid):
    try:
        m, a = load_run(tag, rid)
        return a if m.get("status") == "ok" else None
    except (FileNotFoundError, KeyError):
        return None


def uf_run(c, v):
    if v == "miropolskii":
        a = run("final/uf", f"c{c}__miropolskii")
        return a if a is not None else run("final/uf", f"c{c}__miropolskii_N10dt0.1")
    if v in ("spec11", "specz11", "ds11", "d119s11"):
        return run(f"{F11}/uf", f"c{c}__{v}")
    return run("final/uf", f"c{c}__{v}")


# ------------------------------------------------------------------ F1 film-law parity
def film_laws(df):
    G, rl = df.G, df.rho_l
    u = G / rl
    z = np.maximum(df.z, D)
    ux = G * (1 - df.x.clip(0, 0.9)) / rl
    L = {}
    L[r"this film law ($C$, $\Lambda$ fitted)"] = 15.8 * df.k_vf ** (2 / 3) * (ux * df.rho_vf * df.cp_vf) ** (1 / 3) \
        * D ** (-1 / 3) * (1 / z + 1 / 0.101) ** (1 / 3)
    L[r"$L_{eff}=D$, $C$ = 8.55"] = 8.55 * df.k_vf / D * (u * D * df.rho_vf * df.cp_vf / df.k_vf) ** (1 / 3)
    # Darr et al. 2016 low-g, published constants
    x = df.x.clip(0, 0.99)
    Re = G * D / df.mu_v
    Pr = df.cp_v * df.mu_v / df.k_v
    htc = 0.150 * (1 + (D / (df.z + 0.1 * D)) ** 0.7) * (1 - x) ** -0.76 * Re**0.52 * Pr**0.4 * df.k_v / D
    We = G**2 * D / (rl * df.sigma)
    Ja = df.cp_l * df.dT / df.h_fg
    L["Darr et al. (2016), published"] = htc + 5e-6 * df.k_l / D * np.exp(-0.10 * df.z / D) * We**2 / Ja
    L["Bromley et al. (1953), C = 2.7"] = 2.7 * np.sqrt(df.k_vf * df.rho_vf * (df.h_fg + 0.5 * df.cp_vf * df.dT) * u / (D * df.dT))
    L["vapor Dittus–Boelter, full G"] = 0.023 * Re**0.8 * Pr**0.4 * df.k_v / D
    xm = df.x.clip(0, 1)
    Rem = Re * (xm + df.rho_v / rl * (1 - xm))
    Y = 1 - 0.1 * (rl / df.rho_v - 1) ** 0.4 * (1 - xm) ** 0.4
    L["Miropolskii"] = (0.023 * Rem**0.8 * df.Pr_vw**0.8 * Y * df.k_v / D).clip(lower=1e-3)
    return L


def F1():
    df = pd.read_csv(APR / "film_boiling_dataset.csv")
    w = 1.0 / df.groupby(["case", "st"]).t.transform("size")
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.5), sharex=True, sharey=True)
    cases = sorted(df.case.unique())
    cm = plt.get_cmap("viridis", len(cases))
    for ax, (name, hp) in zip(axs.flat, film_laws(df).items()):
        for i, c in enumerate(cases):
            s = df.case == c
            ax.scatter(df.h_ec[s], hp[s], s=1.5, color=cm(i), alpha=0.35, lw=0, rasterized=True,
                       label=f"test {c}" if ax is axs.flat[0] else None)
        r = np.log(hp / df.h_ec)
        rms = np.sqrt(np.sum(w * r**2) / np.sum(w))
        bias = np.exp(np.sum(w * r) / np.sum(w)) - 1
        lim = [10, 3000]
        ax.plot(lim, lim, "k-", lw=0.6)
        for f in (0.7, 1.3):
            ax.plot(lim, [f * lim[0], f * lim[1]], "k:", lw=0.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_title(name, fontsize=7.5)
        ax.text(0.04, 0.96, f"RMS(ln) = {rms:.2f}\nbias = {bias:+.0%}", transform=ax.transAxes, va="top", fontsize=6.5)
    for ax in axs[1]:
        ax.set_xlabel(r"measured $h$ [W m$^{-2}$ K$^{-1}$]")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"predicted $h$ [W m$^{-2}$ K$^{-1}$]")
    h, l = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, l, markerscale=5, ncol=len(l), fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, -0.005),
               handletextpad=0.1, columnspacing=1.0)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "F1_film_parity")


# ------------------------------------------------------------------ F2 rewet
def F2():
    import CoolProp.CoolProp as CP

    r = pd.read_csv(APR / "rewet_table.csv")
    k = pd.read_csv(APR / "nbs_knees.csv")
    k = k[(k.fluid == "LN2") & (k.ratio >= 1.9) & (k.dT_seg_at_knee <= 20) & (k.z_m < 60)]
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    P = np.geomspace(80e3, 900e3, 40)
    Ts = np.array([CP.PropsSI("T", "P", p, "Q", 0, "Nitrogen") for p in P])
    Tc = 126.192
    Thn = Tc * (0.905 + 0.095 * (Ts / Tc) ** 8)
    Ts_ = Thn + 0.122 * Tc
    for beta, ls, lab in ((0.085, "-", "rule, 304SS wall"), (0.012, "-.", "rule, copper wall")):
        ax.plot(P / 1e3, Ts_ + beta * (Ts_ - Ts), ls, color=COL["spec11"], label=lab)
    ax.plot(P / 1e3, Tc + 0.085 * (Tc - Ts), ":", color=COL["spec11"], lw=0.8, label=r"$T^*=T_c$, 304SS")
    ax.plot(P / 1e3, Thn, "--", color="0.4", label=r"$T_\mathrm{hn}$ (Lienhard 1976)")
    ax.plot(P / 1e3, Ts, "-", color="0.75", lw=0.8, label=r"$T_\mathrm{sat}$")
    for st, mk in (("up", "o"), ("dn", "s")):
        s = r.st == st
        ax.errorbar(r.P_kPa[s], r.T_mfb[s], yerr=[r.T_mfb[s] - r.T_mfb_TCmin[s], r.T_mfb_TCmax[s] - r.T_mfb[s]],
                    fmt=mk, ms=3.5, mfc="w", mec="k", ecolor="k", elinewidth=0.6, capsize=1.5,
                    label=f"low g, 304SS, z = {'0.149' if st == 'up' else '0.401'} m")
    ax.errorbar(k.P_local_kPa, k.T_knee, yerr=[np.zeros(len(k)), k.T_prev - k.T_knee], fmt="^", ms=3.5,
                color="#D55E00", elinewidth=0.6, capsize=1.5, label="NBS 9264, 1 g, copper")
    ax.set_xscale("log")
    ax.set_xlabel("pressure [kPa]")
    ax.set_ylabel(r"rewet temperature $T_\mathrm{MFB}$ [K]")
    ax.set_ylim(70, 152)
    ax.legend(fontsize=5.6, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    save(fig, "F2_rewet")


# ------------------------------------------------------------------ F3 UF histories
def F3(cases=(1, 4, 6, 7, 8, 10)):
    vs = ("miropolskii", "db_vapor", "ds11", "spec11")
    fig, axs = plt.subplots(3, 4, figsize=(W2, 6.0), sharey=True)
    for k, c in enumerate(cases):
        d = load_uf_case(c)
        runs = {v: uf_run(c, v) for v in vs}
        for j, st in enumerate(("up", "dn")):
            ax = axs[k // 2, 2 * (k % 2) + j]
            lo, hi = d.station_band(st)
            ax.fill_between(d.t, lo, hi, color="0.8", lw=0)
            ax.plot(d.t, d.station_mean(st), "k-", lw=1.3)
            for v in vs:
                a = runs[v]
                if a is not None:
                    ax.plot(a["t"], a[f"T_{st}"], color=COL[v], lw=0.9)
            ax.set_title(f"test {c}, z = {0.149 if st == 'up' else 0.401} m", fontsize=7.5)
            ax.set_xlim(0, d.t[-1])
            if k // 2 == 2:
                ax.set_xlabel("t [s]")
            if j == 0 and k % 2 == 0:
                ax.set_ylabel("wall T [K]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    h = [Patch(color="0.8", label="data (TC range)"), Line2D([], [], color="k", lw=1.3, label="data (TC mean)")] + \
        [Line2D([], [], color=COL[v], label=LAB[v]) for v in vs]
    fig.legend(handles=h, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.01), fontsize=6.5)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    save(fig, "F3_uf_histories")


# ------------------------------------------------------------------ F4 quench front
def F4(cases=(6, 7, 8, 9, 10)):
    fig, axs = plt.subplots(1, len(cases), figsize=(W2, 2.2), sharey=False)
    for ax, c in zip(axs, cases):
        d = load_uf_case(c)
        for v, ls in (("specz11", "--"), ("spec11", "-")):
            a = uf_run(c, v)
            if a is None:
                continue
            tq = [cross(a["t"], T, 130.0) for T in a["wall_T"]]
            x = a["wall_x"]
            ok = [i for i, q in enumerate(tq) if q is not None]
            ax.plot(x[ok], [tq[i] for i in ok], ls, color=COL[v], label=LAB[v])
        for st, z in (("up", 0.149), ("dn", 0.401)):
            tcs = [cross(d.t, d.col(tc), 130.0) for tc in UF_STATION_TCS[st]]
            tcs = [q for q in tcs if q is not None]
            ax.plot([z] * len(tcs), tcs, "k.", ms=3)
            tm = cross(d.t, d.station_mean(st), 130.0)
            if tm is not None:
                ax.plot(z, tm, "o", mfc="w", mec="k", ms=4)
        ax.set_title(f"test {c}", fontsize=7.5)
        ax.set_xlabel("z [m]")
        ax.set_xlim(0, 0.572)
    axs[0].set_ylabel("time to 130 K [s]")
    axs[0].plot([], [], "o", mfc="w", mec="k", ms=4, label="data, TC mean")
    axs[0].plot([], [], "k.", ms=3, label="data, single TC")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), fontsize=6.5)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, "F4_quench_front")


# ------------------------------------------------------------------ F5 independent validation
def F5():
    from ofchill.validate4 import darr_data, f2_data

    tag = f"{F11}/val"
    c = COL["spec11"]
    fig, axs = plt.subplots(1, 5, figsize=(W2, 2.2), sharey=True)
    d = f2_data()
    ax = axs[0]
    ax.plot(d.t, d.Twall_exit_avg_K, "k-", lw=1.3)
    for rid, ls in (("f2__spec11", "-"),):
        a = run(tag, rid)
        if a is not None:
            ax.plot(a["t"], a["T_exit"], ls, color=c)
    ax.set_title("F2 (0.05 g)\nz = 0.610 m", fontsize=7)
    ax.set_ylabel("wall T [K]")
    k = 1
    for test in (31, 40):
        dd = darr_data(test)
        for st, colname, z in (("up", "Twall_TC1-3_avg_K", 0.149), ("dn", "Twall_TC4-6_avg_K", 0.401)):
            ax = axs[k]
            k += 1
            ax.plot(dd.t, dd[colname], "k-", lw=1.3)
            for rid, ls in ((f"darr{test}__spec11", "-"), (f"darr{test}__spec11_nog", ":")):
                a = run(tag, rid)
                if a is not None:
                    ax.plot(a["t"], a[f"T_{st}"], ls, color=c)
            ax.set_title(f"test {test} (1 g, upflow)\nz = {z} m", fontsize=7)
    for ax in axs:
        ax.set_xlabel("t [s]")
    axs[0].plot([], [], "k-", lw=1.3, label="data")
    axs[0].plot([], [], "-", color=c, label="closure + buoyancy term")
    axs[0].plot([], [], ":", color=c, label="no buoyancy term")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02), fontsize=6.5)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save(fig, "F5_validation")


# ------------------------------------------------------------------ F6 NBS LN2
def _station_T(a, x):
    return np.array([np.interp(x, a["wall_x"], a["wall_T"][:, i]) for i in range(len(a["t"]))])


def F6():
    from ofchill.validate4 import NBS_STATIONS, nbs_data

    fig, axs = plt.subplots(2, 4, figsize=(W2, 4.2), sharey="row")
    rows = [(10, [("nbs6", "10:spec11_gh:Pout=83082,end=400", COL["spec11"], "-", "this closure"),
                  ("nbs7", "10:spec11_gh:Pout=83082,end=400,f=0.0085", COL["spec11"], "--", "this closure, f = 0.0085"),
                  ("nbs7", "10:spec11_gh:Pout=83082,end=400,f=0.0043", COL["spec11"], ":", "this closure, f = 0.0043"),
                  ("nbs6", "10:miropolskii:Pout=83082,end=400", COL["miropolskii"], "-", "Miropolskii")], "LN$_2$, 2.5 atm"),
            (2, [("nbs6", "2:spec11_gh:Pout=83082,end=128", COL["spec11"], "-", "this closure"),
                 ("nbs6", "2:miropolskii:Pout=83082,end=128", COL["miropolskii"], "--", "Miropolskii")], "LH$_2$, 5.1 atm")]
    for r, (fig_no, curves, title) in enumerate(rows):
        d = nbs_data(fig_no)
        runs = [(run(t, rid), col, ls, lab) for t, rid, col, ls, lab in curves]
        tmax = 0
        for j, (st, x) in enumerate(NBS_STATIONS.items()):
            ax = axs[r, j]
            dd = d[np.isclose(d.station_m, x, atol=0.3)]
            if dd.empty:
                dd = d[np.isclose(d.station_m, 60.4, atol=0.3)]
            dd = dd.sort_values("time_s")
            tmax = max(tmax, dd.time_s.max())
            ax.plot(dd.time_s, dd.T_wall_K, "k-", lw=1.3, label="data")
            for a, col, ls, lab in runs:
                if a is not None:
                    ax.plot(a["t"], _station_T(a, x), ls, color=col, label=lab)
            ax.set_title(f"{title}, z = {x:.1f} m", fontsize=7)
            if r == 1:
                ax.set_xlabel("t [s]")
        for ax in axs[r]:
            ax.set_xlim(0, 1.15 * tmax)
        axs[r, 0].set_ylabel("wall T [K]")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.subplots_adjust(hspace=0.62)
    from matplotlib.lines import Line2D

    for r, (fig_no, curves, title) in enumerate(rows):
        hs = [Line2D([], [], color="k", lw=1.3, label="data")] + \
             [Line2D([], [], color=col, ls=ls, label=lab) for t, rid, col, ls, lab in curves]
        y0 = axs[r, 0].get_position().y0
        fig.legend(handles=hs, loc="upper center", ncol=len(hs), fontsize=6.5, columnspacing=1.5,
                   bbox_to_anchor=(0.5, y0 - (0.045 if r == 0 else 0.095)))
    save(fig, "F6_nbs")


# ------------------------------------------------------------------ F7 tornado
SENS_LABEL = {"lag0": r"meter delay $\tau$ = 0 (2.5 s)", "lag3.5": r"meter delay $\tau$ = 3.5 s (2.5 s)",
              "x0.02": r"inlet quality $x_e$ = 0.02 (0)", "x0.05": r"inlet quality $x_e$ = 0.05 (0)",
              "Pex0.95": "exit pressure ×0.95", "Pex1.05": "exit pressure ×1.05",
              "TMFB-5": r"$T_\mathrm{MFB}$ − 5 K", "TMFB+5": r"$T_\mathrm{MFB}$ + 5 K", "C14.9": r"$C$ = 14.9 (15.8)",
              "C16.4": r"$C$ = 16.4 (15.8)", "Lam0.081": r"$\Lambda$ = 0.081 m (0.101)",
              "Lam0.118": r"$\Lambda$ = 0.118 m (0.101)",
              "parasitic": r"parasitic heat $q''_\mathrm{par}$ included", "no_axial": "axial wall conduction off"}


def F7():
    s = pd.read_csv(RES / F11 / "sens_aggregate.csv")
    u = pd.read_csv(RES / F11 / "uf_aggregate.csv")
    base = u[u.run == "spec11"].iloc[0]
    s["key"] = s.run.str.replace("spec11_", "", regex=False)
    s = s[s.key.isin(SENS_LABEL)].copy()
    s["dRMS"] = s.RMS - base.RMS
    s["dT130"] = s.t130_MAE - base.t130_MAE
    s = s.reindex(s.dRMS.abs().sort_values().index)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.8), sharey=True)
    y = np.arange(len(s))
    for ax, col, lab in ((axs[0], "dRMS", r"$\Delta$ mean RMS [K]"), (axs[1], "dT130", r"$\Delta$ MAE $t_{130}$ [s]")):
        ax.barh(y, s[col], color=[("#0072B2" if v < 0 else "#D55E00") for v in s[col]], height=0.6)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel(lab)
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([SENS_LABEL[k] for k in s.key])
    axs[0].set_title(f"reference: mean station RMS {base.RMS:.1f} K", fontsize=7.5)
    axs[1].set_title(f"reference: MAE of $t_{{130}}$ {base.t130_MAE:.2f} s", fontsize=7.5)
    fig.tight_layout()
    save(fig, "F7_sensitivity")


# ------------------------------------------------------------------ F8 boiling curve
def F8():
    from ofchill.closures import boiling_variant, contact_beta_spec
    from ofchill.props import fluid_props

    q = pd.read_csv(APR / "q_consistent.csv")
    fp = fluid_props("Nitrogen")
    fig, ax = plt.subplots(figsize=(W1, 2.7))
    sel = [(4, "up", "#009E73"), (8, "up", "#0072B2"), (10, "up", "#D55E00")]
    for c, st, col in sel:
        s = (q.case == c) & (q.st == st) & (q.dT > 0.5)
        ax.scatter(q.dT[s], q.q_E[s] / 1e3, s=1.5, color=col, alpha=0.4, lw=0, rasterized=True)
        film = s & q.fb.astype(bool)
        P, G = q.P[film].median(), q.G[film].median()
        b = contact_beta_spec(P, "Nitrogen", "ss304")[0]
        cl = boiling_variant("spec11", beta=b)
        Ts = float(fp.Tsat(P))
        Tw = Ts + np.geomspace(0.3, 210, 600)
        x = 0.03
        rho = 1 / (x / fp.rho_v(P) + (1 - x) / fp.rho_l(P))
        h = cl.evaluate({"z": 0.149}, Tw=Tw, Tf=Ts, P=P, G=G, D=D, rho=rho, mu=1e-4, k=0.1, cp=2000.0)
        ax.plot(Tw - Ts, h * (Tw - Ts) / 1e3, "-", color=col, lw=1.1,
                label=f"test {c}: P = {P / 1e3:.0f} kPa, G = {G:.0f} kg m$^{{-2}}$s$^{{-1}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 250)
    ax.set_ylim(0.5, 400)
    ax.set_xlabel(r"$T_w - T_\mathrm{sat}$ [K]")
    ax.set_ylabel(r"$q''$ [kW m$^{-2}$]")
    ax.legend(fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.2), title="closure (lines), data (dots), z = 0.149 m",
              title_fontsize=5.8)
    save(fig, "F8_boiling_curve")


# ------------------------------------------------------------------ F9 chilldown time
F9_ORDER = ["spec11", "specz11", "d119s11", "ds11", "d119s", "ds", "miropolskii", "db_vapor"]
F9_LAB = {"spec11": "this closure", "specz11": r"film law with $L_{eff}=D$", "d119s11": "Darr (2016) ×1.19",
          "ds11": "Darr (2016)", "d119s": "Darr ×1.19, other rewet rule", "ds": "Darr, other rewet rule",
          "miropolskii": "Miropolskii", "db_vapor": "vapor Dittus–Boelter"}
F9_COL = {"spec11": "#0072B2", "specz11": "#56B4E9", "d119s11": "#E69F00", "ds11": "#CC79A7", "d119s": "#E69F00",
          "ds": "#CC79A7", "miropolskii": "#D55E00", "db_vapor": "#009E73"}
F9_MK = {"UF": "o", "F2 (0.05 g)": "D", "Darr 1 g": "^", "NBS LN2": "s", "NBS LH2": "v"}


def F9():
    df = pd.read_csv(RES / F11 / "chilldown_stations.csv")
    df = df[(df.outcome != "no run") & df.closure.isin(["spec11", "specz11", "d119s11", "ds11", "miropolskii", "db_vapor"])]
    fig = plt.figure(figsize=(W2, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25])
    ax = fig.add_subplot(gs[0])
    lo, hi = 0.5, 400.0
    for cl in ["db_vapor", "miropolskii", "ds11", "d119s11", "specz11", "spec11"]:
        g = df[df.closure == cl]
        for ds, mk in F9_MK.items():
            h = g[(g.dataset == ds) & (g.outcome == "hit")]
            ax.scatter(h.tcd_d, h.tcd_m, marker=mk, s=10 if cl == "spec11" else 6, facecolor=F9_COL[cl] if cl == "spec11" else "none",
                       edgecolor=F9_COL[cl], lw=0.6, zorder=3 if cl == "spec11" else 2)
            m = g[(g.dataset == ds) & (g.outcome == "miss")]
            ax.scatter(m.tcd_d, np.full(len(m), hi * 0.93), marker=mk, s=8, facecolor="none", edgecolor=F9_COL[cl],
                       lw=0.6, clip_on=False)
    ax.plot([lo, hi], [lo, hi], "k-", lw=0.6)
    ax.fill_between([lo, hi], [lo * 0.9, hi * 0.9], [lo * 1.1, hi * 1.1], color="0.88", lw=0, zorder=0)
    ax.axhline(hi * 0.86, color="0.5", lw=0.4, ls=":")
    ax.text(lo * 1.15, hi * 0.82, "not chilled within the test", fontsize=5.5, va="top", color="0.3")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("measured time to chill down [s]")
    ax.set_ylabel("predicted time to chill down [s]")
    from matplotlib.lines import Line2D

    hds = [Line2D([], [], marker=mk, ls="", mfc="none", mec="k", ms=4, label=ds) for ds, mk in F9_MK.items()]
    fig.legend(handles=hds, fontsize=6, loc="lower center", ncol=len(hds), bbox_to_anchor=(0.5, -0.005), columnspacing=1.5)
    hcl = [Line2D([], [], marker="o", ls="", mfc=F9_COL[c], mec=F9_COL[c], ms=3.5, label=F9_LAB[c])
           for c in ["spec11", "specz11", "d119s11", "ds11", "miropolskii", "db_vapor"]]
    fig.legend(handles=hcl, fontsize=6, loc="lower center", ncol=len(hcl), bbox_to_anchor=(0.5, 0.06), columnspacing=1.2)
    # right: share of stations predicted within 10 % / 2 s, and missed
    ax2 = fig.add_subplot(gs[1])
    order = ["spec11_all", "spec11", "specz11", "d119s11", "ds11", "miropolskii", "db_vapor"]
    common = df.dataset.isin(["UF", "F2 (0.05 g)", "Darr 1 g"])
    stats = []
    for cl in order:
        if cl == "spec11_all":
            g = df[(df.closure == "spec11") & df.tcd_d.notna()]
        else:
            g = df[(df.closure == cl) & df.tcd_d.notna() & common]
        n = len(g)
        hit = g[g.outcome == "hit"]
        rel = (hit.err / hit.tcd_d).abs()
        stats.append((n, (rel <= 0.1).sum() / n, ((rel <= 0.25) & (rel > 0.1)).sum() / n,
                      (rel > 0.25).sum() / n, (g.outcome == "miss").sum() / n))
    stats = np.array(stats)
    y = np.arange(len(order))[::-1]
    left = np.zeros(len(order))
    for k, (lab, col) in enumerate((("within 10 %", "#0072B2"), ("10–25 %", "#56B4E9"), ("> 25 %", "#E69F00"),
                                    ("not chilled in model", "#D55E00"))):
        ax2.barh(y, stats[:, k + 1], left=left, color=col, height=0.6, label=lab)
        left += stats[:, k + 1]
    ax2.set_yticks(y)
    names = {**F9_LAB, "spec11_all": "this closure, incl. NBS lines"}
    ax2.set_yticklabels([f"{names[c]} (n = {int(n)})" for c, n in zip(order, stats[:, 0])])
    ax2.axhline(y[0] - 0.5, color="0.5", lw=0.5)
    ax2.set_xlim(0, 1)
    ax2.set_xlabel("share of stations that chill down in the test\n(UF, F2 and 1-g tests; top bar: all five data sets)")
    ax2.legend(fontsize=5.5, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=4, handlelength=1.0)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    save(fig, "F9_chilldown_time")


# ------------------------------------------------------------------ F10 parsimony ladder
def F10():
    """Held-out error against the number of fitted constants (data: analysis/apriori/out/parsimony.csv).
    Redrawn at the paper's figure width (the a19 script draws 15 in wide, which prints at a quarter of its size)."""
    from matplotlib.lines import Line2D

    tab = pd.read_csv(APR / "parsimony.csv")
    k = tab.k_fit.values
    final = (tab.model == "M3").values
    blue, orange, red = "#0072B2", "#D55E00", "#CC3311"
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.6))

    def ladder(ax, y, **kw):
        ax.plot(k, y, "-", color="k", lw=0.8, zorder=1)
        ax.scatter(k[~final], y[~final], s=16, color=blue, zorder=3)
        ax.scatter(k[final], y[final], s=22, color=red, zorder=4)
        for m, x_, y_ in zip(tab.model, k, y):
            ax.annotate(m, (x_, y_), textcoords="offset points", xytext=kw.get("off", {}).get(m, (4, 4)), fontsize=6.5)

    ax = axs[0]
    ladder(ax, tab.h_loto.values)
    ax.plot(k, tab.h_rms_in, "--", color="0.5", lw=0.9)
    ax.axhline(0.174, color="#117733", ls=":", lw=1.0)
    ax.set_yscale("log")
    ax.set_ylabel(r"RMS of $\ln(h_\mathrm{pred}/h_\mathrm{meas})$")
    ax.legend(handles=[Line2D([], [], color="k", lw=0.8, marker="o", ms=3, mfc=blue, mec=blue, label="withheld test"),
                       Line2D([], [], color="0.5", ls="--", lw=0.9, label="fitted on all ten tests"),
                       Line2D([], [], color="#117733", ls=":", lw=1.0, label="within-run noise floor")],
              loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.5, ncol=1, handlelength=2.0)

    ax = axs[1]
    ladder(ax, tab.Tmfb_loto_mae.values, off={"M0": (5, -12), "M1": (5, -12)})
    ax.plot(k, tab.nbs_knee_mae, "s--", color=orange, ms=3.5, lw=0.9, zorder=2)
    ax.set_ylabel(r"rewet temperature MAE [K]")
    ax.legend(handles=[Line2D([], [], color="k", lw=0.8, marker="o", ms=3, mfc=blue, mec=blue, label="calibration stations (withheld)"),
                       Line2D([], [], color=orange, ls="--", marker="s", ms=3.5, lw=0.9, label=r"NBS LN$_2$ knees on copper (not fitted)")],
              loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.5, ncol=1, handlelength=2.0)

    ax = axs[2]
    ladder(ax, tab.t130_mae.values)
    ax.plot(k, tab.t100_mae, "^--", color=orange, ms=3.5, lw=0.9, zorder=2)
    ax.set_ylabel(r"quench-time MAE [s]")
    ax.legend(handles=[Line2D([], [], color="k", lw=0.8, marker="o", ms=3, mfc=blue, mec=blue, label=r"time to 130 K"),
                       Line2D([], [], color=orange, ls="--", marker="^", ms=3.5, lw=0.9, label=r"time to $T_\mathrm{sat}$+8 K / 100 K")],
              loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.5, ncol=1, handlelength=2.0)

    for ax in axs:
        ax.set_xlabel("fitted constants")
        ax.set_xticks(range(0, 8))
        ax.grid(alpha=0.3, lw=0.4)
    fig.tight_layout(w_pad=1.2)
    save(fig, "F10_parsimony")


ALL = {"F1": F1, "F2": F2, "F3": F3, "F4": F4, "F5": F5, "F6": F6, "F7": F7, "F8": F8, "F9": F9, "F10": F10}

if __name__ == "__main__":
    style()
    for k in (sys.argv[1:] or ALL):
        try:
            ALL[k]()
        except Exception as e:  # keep going; report
            print("FAILED", k, type(e).__name__, e)
