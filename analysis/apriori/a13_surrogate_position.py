"""Phase 3 / Task 4: 0-D surrogate (a10 boiling curve, station-mean TC comparison) with
position-dependent film laws and an optional inlet (entrance) rewet increment.
Front metric: dt_q = t130(dn) - t130(up), model vs measured -> apparent front speed 0.252 m / dt_q."""
import numpy as np
import pandas as pd

from common import D, OUT, md_table
import a10_boiling_curve_sensitivity as a10
import a06_surrogate as s6

LAM3, C3 = np.exp(-2.045), 16.38      # L3
LAM3X, C3X = np.exp(-2.294), 15.80    # L3x


def pe13(d, G, ul):
    return d.k_vf / D * (ul * D * d.rho_vf * d.cp_vf / d.k_vf) ** (1 / 3)


H = {
    "L0 (Phase 2)": lambda d, G: 8.55 * pe13(d, G, G / d.rho_l),
    "L3 Lam=0.129 m": lambda d, G: C3 * pe13(d, G, G / d.rho_l) * (D * (d.z + LAM3) / (d.z * LAM3)) ** (1 / 3),
    "L3x Lam=0.101 m, u_l=G(1-x)/rho_l": lambda d, G: C3X * pe13(d, G, G * (1 - np.clip(d.x, 0, 0.9)) / d.rho_l)
        * (D * (d.z + LAM3X) / (d.z * LAM3X)) ** (1 / 3),
    "L0 x Darr entrance": lambda d, G: 7.60 * pe13(d, G, G / d.rho_l) * (1 + (D / (d.z + 0.1 * D)) ** 0.7),
}
A_ENT = 8.4 / np.exp(-0.149 / LAM3X)   # K, inlet rewet increment amplitude, decays with Lam


def run(recs, hname, ent):
    a10.h_film = H[hname]
    rr = [dict(r, Tmfb=r["Tmfb"] + (A_ENT * np.exp(-(0.149 if r["st"] == "up" else 0.401) / LAM3X) if ent else 0.0))
          for r in recs]
    tt, sim = a10.simulate(rr, dict(a10.BASE))
    m = a10.metrics(recs, tt, sim)
    # front: t130 crossing up vs dn
    fr = []
    for c in sorted(m.case.unique()):
        ru = [r for r in recs if r["case"] == c and r["st"] == "up"][0]
        rd = [r for r in recs if r["case"] == c and r["st"] == "dn"][0]
        ju = [j for j, r in enumerate(recs) if r is ru][0]; jd = [j for j, r in enumerate(recs) if r is rd][0]
        tmu, tmd = s6.t_cross(ru["t"], ru["Tm"], 130), s6.t_cross(rd["t"], rd["Tm"], 130)
        tsu = s6.t_cross(tt, sim[:, ju], 130); tsd = s6.t_cross(tt, sim[:, jd], 130)
        fr.append(dict(case=c, dtq_meas=tmd - tmu, dtq_sim=tsd - tsu))
    fr = pd.DataFrame(fr).dropna(subset=["dtq_meas"])
    up, dn = m[m.st == "up"], m[m.st == "dn"]
    out = dict(film=hname, entrance_dTmfb=ent, **a10.summ(m),
               t130_bias_up=up.dt130.mean(), t130_bias_dn=dn.dt130.mean(),
               rms_up=up.rms.mean(), rms_dn=dn.rms.mean(),
               dtq_meas=";".join(f"{v:.1f}" for v in fr.dtq_meas), dtq_sim=";".join(f"{v:.1f}" for v in fr.dtq_sim),
               U_meas_mm_s=float(np.nanmedian(252 / fr.dtq_meas)), U_sim_mm_s=float(np.nanmedian(252 / fr.dtq_sim.clip(lower=0.05))))
    return out, m


def main():
    recs = a10.inputs()
    rows, per = [], {}
    for hn in H:
        for ent in (False, True):
            if ent and hn not in ("L0 (Phase 2)", "L3x Lam=0.101 m, u_l=G(1-x)/rho_l", "L3 Lam=0.129 m"):
                continue
            o, m = run(recs, hn, ent)
            rows.append(o); per[(hn, ent)] = m
            print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in o.items()}, flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/surrogate_position.csv", index=False)
    with open(f"{OUT}/t13_surrogate_position.md", "w") as f:
        f.write(__doc__ + f"\nentrance rewet increment: dT_MFB,ent = {A_ENT:.1f} K exp(-z/{LAM3X:.3f} m) (+8.4 K at z=0.149, "
                f"+{A_ENT*np.exp(-0.401/LAM3X):.1f} K at 0.401)\n\n" + md_table(res.round(2), 3) + "\n")
        for k in [("L3x Lam=0.101 m, u_l=G(1-x)/rho_l", True), ("L3x Lam=0.101 m, u_l=G(1-x)/rho_l", False)]:
            f.write(f"\n## per station: {k}\n\n" + md_table(per[k].round(2), 3) + "\n")


if __name__ == "__main__":
    main()
