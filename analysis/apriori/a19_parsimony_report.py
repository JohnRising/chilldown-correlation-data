"""Phase 6: constant inventory, comparison with published correlation sets, parsimony figure and table
(reads out/parsimony.csv from a18)."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIG, OUT, md_table

INV = [
    # (item, value, class, transfer note)
    ("Film-law prefactor C", "15.7 (LOTO 14.9-16.4)", "FITTED (UF)", "renewal intensity; re-calibrate per system (first knob)"),
    ("Development length Lambda", "0.098 m = 8.4 D (LOTO 0.081-0.118)", "FITTED (UF)", "entrance/flow-development length; expect ~ D-scaled, re-check per geometry"),
    ("Rewet offset e", "0.122 (e Tc = 15.4 K N2)", "FITTED (UF)", "apparent-onset offset above T_hn; keep, check with coated/LH2 data"),
    ("Exponent 1/3 on Pe_l and on D/L_eff", "1/3", "THEORY (Leveque scaling; free fit 0.34)", "keep"),
    ("u_l = G(1-x_e)/rho_l (superficial liquid velocity)", "-", "MODELLING CHOICE (x dependence unidentified)", "keep; suspect for x_e > 0.2"),
    ("Property states (film-T vapour; sat. liquid)", "-", "MODELLING CHOICE (selected on data)", "keep"),
    ("T_hn = Tc(0.905 + 0.095 Tr^8)", "0.905, 0.095", "LITERATURE (Lienhard 1976)", "keep"),
    ("Contact factor beta = e_l/e_w", "0.085 SS, 0.012 Cu", "THEORY (semi-inf. contact) + NIST props", "keep"),
    ("Transition exponent n", "0.5", "FITTED to per-TC transition shape (insensitive)", "keep (n = 0.25-1 changes <0.5 s)"),
    ("Cooper NB 55, 0.12, -0.55, -0.5; R_p = 1 um", "-", "LITERATURE (Cooper 1984)", "keep (insensitive)"),
    ("Dittus-Boelter liquid floor 0.023, 0.8, 0.4", "-", "LITERATURE", "keep"),
    ("Katto-Ohno CHF regime constants", "-", "LITERATURE (Katto-Ohno 1984)", "keep (insensitive)"),
    ("Vapour DB limit for x_e >= 0.15 / laminar floor 4.36", "-", "LITERATURE (guard, untested)", "keep"),
    ("Gravity term (Bromley 0.62/0.943, cube sum)", "-", "LITERATURE (optional; over-predicts g effect ~2x at low G)", "1 g only"),
    ("Meter delay tau", "2.5 s (LOTO 2.0-2.75)", "RIG BC (fitted)", "not part of closure"),
    ("Parasitic heat q_par", "98.6 + 2.2(293-T) W/m2", "RIG BC (fitted to authors' inverse model)", "not part of closure"),
    ("Energy-consistent flux correction", "1/0.735 for tests 1,2,3,5,6", "DATA REDUCTION", "not part of closure"),
]
OTHER = [
    ("This work (final)", "film 2 + rewet 1 (+ n = 0.5, insensitive)", 3, "LN2, 1 tube, low g"),
    ("Darr et al. 2016 npj, low-g FB", "c1-c7 (c7 = 0 in low g); rewet/TB/NB from 1-g sets", 7, "LN2, same tube"),
    ("Darr et al. 2019 JHT 'UF universal'", "FB c1-c9, Twet c1-c4, TB c0-c1, CHF c1, NB c1-c5", 21, "LN2+LH2, one tube size"),
    ("Darr & Hartwig 2020 (LH2)", "FB c1-c10 + p (+ rewet 2 and NB 3 adopted)", 11, "LH2 GRC"),
    ("Darr et al. 2019b rewet (eq. 28)", "c1, c2 (+ Baumeister constants)", 2, "LN2 1 g"),
    ("Miropolskii 1963", "0 fitted here (water constants)", 0, "water, high P"),
    ("Groeneveld 1973 tube", "a, b, c, d (literature)", 0, "water"),
]


def main():
    tab = pd.read_csv(f"{OUT}/parsimony.csv")
    cons = json.load(open(f"{OUT}/parsimony_constants.json"))
    tab["d_loto_per_const"] = tab.h_loto.diff() / tab.k_fit.diff()
    # figure
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    lab = {"M0": "M0 Bromley + T_hn", "M1": "M1 C", "M2": "M2 C + e", "M3": "M3 C + Λ + e (final)",
           "M4": "M4 + free n, m", "M5": "M5 Darr refit"}
    for a, (y, yl) in zip(axs, (("h_loto", "held-out RMS ln(h_pred/h_meas)"), ("Tmfb_loto_mae", "rewet MAE [K]"),
                                ("t130_mae", "held-out quench-time MAE, lumped model [s]"))):
        a.plot(tab.k_fit, tab[y], "k-", lw=0.8, zorder=1)
        for _, r in tab.iterrows():
            c = "C3" if r.model == "M3" else "C0"
            a.scatter(r.k_fit, r[y], s=50, color=c, zorder=2)
            if y != "Tmfb_loto_mae" or r.model in ("M1", "M2"):
                a.annotate(lab[r.model] if y != "Tmfb_loto_mae" else ("T* = T_hn (0 fitted)" if r.model == "M1" else "T* = T_hn + e Tc"),
                           (r.k_fit, r[y]), textcoords="offset points", xytext=(4, 4), fontsize=7)
        a.set_xlabel("number of constants fitted to the calibration data"); a.set_ylabel(yl); a.grid(alpha=0.3)
        if y == "h_loto":
            a.plot(tab.k_fit, tab.h_rms_in, "C7--", lw=0.8, label="fitted on all ten tests")
            a.axhline(0.174, color="C2", ls=":", lw=1, label="within-run noise floor")
            a.legend(fontsize=7); a.set_yscale("log")
        if y == "Tmfb_loto_mae":
            a.plot(tab.k_fit, tab.nbs_knee_mae, "s--", color="C1", ms=4, label="NBS LN2 knees on copper (not fitted)")
            a.legend(fontsize=7)
        if y == "t130_mae":
            a.plot(tab.k_fit, tab.t100_mae, "^--", color="C1", ms=4, label="time to Tsat+8 K / 100 K")
            a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig21_parsimony.png", dpi=150); fig.savefig(f"{FIG}/fig21_parsimony.pdf"); plt.close(fig)
    inv = pd.DataFrame(INV, columns=["item", "value", "class", "transferability"])
    oth = pd.DataFrame(OTHER, columns=["closure", "fitted constants", "n_fitted", "calibration data"])
    cols = ["model", "k_fit", "film_consts", "h_rms_in", "h_loto", "d_loto_per_const", "station_bias_sd_loto", "Tmfb_e",
            "Tmfb_loto_mae", "nbs_knee_mae", "AIC_station", "BIC_station", "n_eff_ar1", "AIC_ar1", "BIC_ar1",
            "rms_mean", "t130_mae", "t100_mae", "miss130", "miss100"]
    with open(f"{OUT}/t18_parsimony.md", "w") as f:
        f.write("# Parsimony study (Phase 6)\n\n## 1. Inventory of numbers in the final closure\n\n" + md_table(inv) +
                "\n\n## 2. Fitted-constant counts of related closures\n\n" + md_table(oth) +
                "\n\n## 3. Ablation ladder (LOTO by test; tau = 2.5 s fixed as rig BC)\n\n"
                "Models: M0 Bromley C=2.7 + contact at T_hn (0 fitted); M1 C (k_vf/D)Pe_l^(1/3), u_l=G(1-x)/rho_l + T_hn contact (1); "
                "M2 = M1 + rewet offset e (2); M3 = M2 + Lambda (3, final); M4 = M3 + free Pe and P exponents (5); "
                "M5 Darr 2016 c1-c6 refit + e (7).\n\n"
                "AIC/BIC (film part only, k = film constants): *station* = n 20 test-stations, residual = station-mean ln error "
                "(within-run residual lag-1 autocorrelation 0.97, so pointwise samples are not independent and the error is "
                "dominated by run-level bias); *ar1* = n_eff = N(1-rho)/(1+rho) ~ 70-77.\n\n"
                + md_table(tab[cols].round(3), 4) + "\n\nfold constants: out/parsimony_constants.json\n")
    print(tab[cols].round(3).to_string())


if __name__ == "__main__":
    main()
