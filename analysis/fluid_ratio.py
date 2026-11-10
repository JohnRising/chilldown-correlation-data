"""Vapour-sensible-to-latent ratio S and liquid-supply ratio Phi for cryogenic line chilldown.

  cd chilldown-lowg && .venv/bin/python analysis/fluid_ratio.py      -> analysis/out_fluid/{S_table.md, phi_table.md,
                                                                         fluid_ratio.json}; paper/figures/F11_fluid_regimes.*

S   = [h_v(T_w0, P) - h_v,sat(P)] / h_fg(P) = int_{Tsat}^{T_w0} cp_v dT / h_fg   (CoolProp, isobaric)
      enthalpy a unit mass of saturated vapour can still absorb by warming to the initial wall temperature,
      relative to the latent heat.
Phi = m_liq,delivered(t_cd) * h_fg(P) / E_wall
      latent capacity of all liquid delivered up to the measured chilldown time, relative to the enthalpy
      removed from the wall (initial wall T -> Tsat(P)).  Phi < 1 means that part of the wall energy must be
      carried away as vapour superheat.  m from the simulated inlet mass flow (final closure runs);
      E_wall from the NIST wall cp (OpenFLUME presets) integrated node by node.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ofchill.runner import load_run, wall_cp  # noqa: E402

OUT = ROOT / "analysis" / "out_fluid"
FIG = ROOT / "paper" / "figures"
ATM = 101325.0
FLUIDS = {"Hydrogen": "H$_2$", "Helium": "He", "Neon": "Ne", "Nitrogen": "N$_2$", "Argon": "Ar",
          "Oxygen": "O$_2$", "Methane": "CH$_4$"}


def S_ratio(fluid, P, Tw0):
    if P >= CP.PropsSI("Pcrit", fluid) * 0.95:
        return np.nan, np.nan, np.nan
    hl = CP.PropsSI("H", "P", P, "Q", 0, fluid)
    hv = CP.PropsSI("H", "P", P, "Q", 1, fluid)
    hT = CP.PropsSI("H", "P", P, "T", Tw0, fluid)
    return (hT - hv) / (hv - hl), hv - hl, hT - hv


def wall_enthalpy(a, material, rho, D, OD, L, Tsat):
    x = np.asarray(a["wall_x"])
    T0 = np.asarray(a["wall_T"])[:, 0]
    N = len(x) - 1
    m_node = rho * np.pi / 4 * (OD**2 - D**2) * L / N * np.where((np.arange(N + 1) == 0) | (np.arange(N + 1) == N), 0.5, 1)
    E = 0.0
    for m, T in zip(m_node, T0):
        Tg = np.linspace(Tsat, T, 400)
        E += m * np.trapezoid(wall_cp(material, Tg), Tg)
    return E, float(np.mean(T0))


def delivered(a, t_end):
    t = np.asarray(a["t"])
    md = np.asarray(a["mdot_in"])
    k = t <= t_end
    return float(np.trapezoid(md[k], t[k]))


def phi_rows():
    rows = []
    nbs = pd.read_csv(ROOT / "results" / "nbs6" / "final_runs.csv")
    nbs = nbs[nbs.variant == "spec11_gh"]
    from ofchill.validate4 import nbs_meta
    meta = nbs_meta().set_index("figure")
    for _, r in nbs.iterrows():
        fig = int(r.fig)
        fluid = "Nitrogen" if r.fluid == "LN2" else "Hydrogen"
        P = float(meta.loc[fig, "P_drive_Pa"])
        _, a = load_run("nbs6", r.run)
        Ts = CP.PropsSI("T", "P", P, "Q", 0, fluid)
        E, Tw0 = wall_enthalpy(a, "ofhc-copper", 8960.0, 0.015875, 0.01905, 60.96, Ts)
        m = delivered(a, float(r.chill_d))
        S, hfg, _ = S_ratio(fluid, P, Tw0)
        rows.append(dict(set=f"NBS {r.fluid}", case=fig, inlet=r.inlet, P_kPa=P / 1e3, Tw0=Tw0, t_cd=float(r.chill_d),
                         m_kg=m, latent_MJ=m * hfg / 1e6, E_wall_MJ=E / 1e6, Phi=m * hfg / E, S=S,
                         local_ok=(fluid == "Hydrogen")))
    meta_uf = json.loads((ROOT / "data" / "uf2016" / "meta.json").read_text())["cases"]
    st = pd.read_csv(ROOT / "results" / "final11" / "chilldown_stations.csv")
    for c in range(1, 11):
        try:
            _, a = load_run("final11/uf", f"c{c}__spec11")
        except Exception:
            continue
        P = meta_uf[str(c)]["P1_mean"] * 1e3
        Ts = CP.PropsSI("T", "P", P, "Q", 0, "Nitrogen")
        E, Tw0 = wall_enthalpy(a, "stainless-steel-304", 8027.0, 0.01168, 0.0127, 0.572, Ts)
        s = st[(st.dataset == "UF") & (st.test.astype(str) == str(c)) & (st.closure == "spec11")]
        tcd = s.tcd_d.max() if s.tcd_d.notna().all() and len(s) else np.nan
        t_end = tcd if np.isfinite(tcd) else float(a["t"][-1])
        m = delivered(a, t_end)
        S, hfg, _ = S_ratio("Nitrogen", P, Tw0)
        rows.append(dict(set="UF low g", case=c, inlet="saturated", P_kPa=P / 1e3, Tw0=Tw0, t_cd=tcd, m_kg=m,
                         latent_MJ=m * hfg / 1e6, E_wall_MJ=E / 1e6, Phi=m * hfg / E, S=S, local_ok=True,
                         completed=bool(np.isfinite(tcd))))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    S_rows = []
    for f in FLUIDS:
        for Tw0 in (300.0, 200.0):
            for p in (1, 2, 3, 5):
                S, hfg, dh = S_ratio(f, p * ATM, Tw0)
                S_rows.append(dict(fluid=f, P_atm=p, Tw0=Tw0, S=S, h_fg_kJ=hfg / 1e3, dh_vap_kJ=dh / 1e3))
    Sdf = pd.DataFrame(S_rows)
    piv = Sdf.pivot_table(index="fluid", columns=["Tw0", "P_atm"], values="S").round(2)
    (OUT / "S_table.md").write_text("S = (h_v(T_w0,P) - h_v,sat)/h_fg\n\n" + piv.to_markdown() + "\n")
    ph = phi_rows()
    (OUT / "phi_table.md").write_text(ph.round(3).to_markdown(index=False) + "\n")
    summ = {"S": Sdf.round(3).to_dict(orient="records"), "phi": ph.round(4).to_dict(orient="records")}
    (OUT / "fluid_ratio.json").write_text(json.dumps(summ, indent=1, default=float))
    print(piv.to_string())
    print(ph.round(3).to_string())
    figure(Sdf, ph)


def figure(Sdf, ph):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
                         "font.size": 8, "legend.fontsize": 6.5, "legend.frameon": False, "xtick.direction": "in",
                         "ytick.direction": "in", "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.25]})
    order = ["Hydrogen", "Helium", "Neon", "Methane", "Nitrogen", "Oxygen", "Argon"]
    y = np.arange(len(order))
    for Tw0, mk, col in ((300.0, "o", "#0072B2"), (200.0, "s", "#D55E00")):
        for i, f in enumerate(order):
            s = Sdf[(Sdf.fluid == f) & (Sdf.Tw0 == Tw0)].dropna()
            if s.empty:
                continue
            a1.plot([s.S.min(), s.S.max()], [i, i], "-", color=col, lw=2, alpha=0.6)
            a1.plot(s[s.P_atm == 1].S, [i], mk, color=col, ms=4,
                    label=(f"$T_{{w0}}$ = {Tw0:.0f} K" if i == 0 else None))
    a1.set_yticks(y)
    a1.set_yticklabels([FLUIDS[f] for f in order])
    a1.set_xscale("log")
    a1.axvline(1, color="0.5", lw=0.6, ls=":")
    a1.set_xlabel(r"$S=\int_{T_\mathrm{sat}}^{T_{w0}} c_{p,v}\,\mathrm{d}T\,/\,h_{fg}$")
    a1.set_title("bars: 1–5 atm (He: 1–2 atm); symbols: 1 atm", fontsize=7)
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2)
    for (setname, mk), g in ph.groupby(["set", ph.set.map({"UF low g": "o", "NBS LN2": "s", "NBS LH2": "^"})]):
        ok = g.local_ok.values
        if ok.any():
            a2.scatter(g.S[ok], g.Phi[ok], marker=mk, s=22, color="#0072B2", label=f"{setname}: local knees reproduced")
        if (~ok).any():
            a2.scatter(g.S[~ok], g.Phi[~ok], marker=mk, s=22, facecolor="w", edgecolor="#D55E00",
                       label=f"{setname}: local knees not reproduced")
    a2.axhline(1, color="0.5", lw=0.6, ls=":")
    a2.axvline(1, color="0.5", lw=0.6, ls=":")
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlabel(r"$S$ at the run conditions")
    a2.set_ylabel(r"$\Phi$ = latent capacity delivered / wall enthalpy")
    a2.text(0.03, 0.04, "UF: 0.57 m, low g; NBS: 61 m, 1 g", transform=a2.transAxes, fontsize=6)
    a2.legend(loc="upper right", fontsize=5.8)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "F11_fluid_regimes.pdf")
    fig.savefig(FIG / "F11_fluid_regimes.png", dpi=300)


if __name__ == "__main__":
    main()
