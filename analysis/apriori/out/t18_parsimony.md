# Parsimony study (Phase 6)

## 1. Inventory of numbers in the final closure

| item | value | class | transferability |
|---|---|---|---|
| Film-law prefactor C | 15.7 (LOTO 14.9-16.4) | FITTED (UF) | renewal intensity; re-calibrate per system (first knob) |
| Development length Lambda | 0.098 m = 8.4 D (LOTO 0.081-0.118) | FITTED (UF) | entrance/flow-development length; expect ~ D-scaled, re-check per geometry |
| Rewet offset e | 0.122 (e Tc = 15.4 K N2) | FITTED (UF) | apparent-onset offset above T_hn; keep, check with coated/LH2 data |
| Exponent 1/3 on Pe_l and on D/L_eff | 1/3 | THEORY (Leveque scaling; free fit 0.34) | keep |
| u_l = G(1-x_e)/rho_l (superficial liquid velocity) | - | MODELLING CHOICE (x dependence unidentified) | keep; suspect for x_e > 0.2 |
| Property states (film-T vapour; sat. liquid) | - | MODELLING CHOICE (selected on data) | keep |
| T_hn = Tc(0.905 + 0.095 Tr^8) | 0.905, 0.095 | LITERATURE (Lienhard 1976) | keep |
| Contact factor beta = e_l/e_w | 0.085 SS, 0.012 Cu | THEORY (semi-inf. contact) + NIST props | keep |
| Transition exponent n | 0.5 | FITTED to per-TC transition shape (insensitive) | keep (n = 0.25-1 changes <0.5 s) |
| Cooper NB 55, 0.12, -0.55, -0.5; R_p = 1 um | - | LITERATURE (Cooper 1984) | keep (insensitive) |
| Dittus-Boelter liquid floor 0.023, 0.8, 0.4 | - | LITERATURE | keep |
| Katto-Ohno CHF regime constants | - | LITERATURE (Katto-Ohno 1984) | keep (insensitive) |
| Vapour DB limit for x_e >= 0.15 / laminar floor 4.36 | - | LITERATURE (guard, untested) | keep |
| Gravity term (Bromley 0.62/0.943, cube sum) | - | LITERATURE (optional; over-predicts g effect ~2x at low G) | 1 g only |
| Meter delay tau | 2.5 s (LOTO 2.0-2.75) | RIG BC (fitted) | not part of closure |
| Parasitic heat q_par | 98.6 + 2.2(293-T) W/m2 | RIG BC (fitted to authors' inverse model) | not part of closure |
| Energy-consistent flux correction | 1/0.735 for tests 1,2,3,5,6 | DATA REDUCTION | not part of closure |

## 2. Fitted-constant counts of related closures

| closure | fitted constants | n_fitted | calibration data |
|---|---|---|---|
| This work (final) | film 2 + rewet 1 (+ n = 0.5, insensitive) | 3 | LN2, 1 tube, low g |
| Darr et al. 2016 npj, low-g FB | c1-c7 (c7 = 0 in low g); rewet/TB/NB from 1-g sets | 7 | LN2, same tube |
| Darr et al. 2019 JHT 'UF universal' | FB c1-c9, Twet c1-c4, TB c0-c1, CHF c1, NB c1-c5 | 21 | LN2+LH2, one tube size |
| Darr & Hartwig 2020 (LH2) | FB c1-c10 + p (+ rewet 2 and NB 3 adopted) | 11 | LH2 GRC |
| Darr et al. 2019b rewet (eq. 28) | c1, c2 (+ Baumeister constants) | 2 | LN2 1 g |
| Miropolskii 1963 | 0 fitted here (water constants) | 0 | water, high P |
| Groeneveld 1973 tube | a, b, c, d (literature) | 0 | water |

## 3. Ablation ladder (LOTO by test; tau = 2.5 s fixed as rig BC)

Models: M0 Bromley C=2.7 + contact at T_hn (0 fitted); M1 C (k_vf/D)Pe_l^(1/3), u_l=G(1-x)/rho_l + T_hn contact (1); M2 = M1 + rewet offset e (2); M3 = M2 + Lambda (3, final); M4 = M3 + free Pe and P exponents (5); M5 Darr 2016 c1-c6 refit + e (7).

AIC/BIC (film part only, k = film constants): *station* = n 20 test-stations, residual = station-mean ln error (within-run residual lag-1 autocorrelation 0.97, so pointwise samples are not independent and the error is dominated by run-level bias); *ar1* = n_eff = N(1-rho)/(1+rho) ~ 70-77.

| model | k_fit | film_consts | h_rms_in | h_loto | d_loto_per_const | station_bias_sd_loto | Tmfb_e | Tmfb_loto_mae | nbs_knee_mae | AIC_station | BIC_station | n_eff_ar1 | AIC_ar1 | BIC_ar1 | rms_mean | t130_mae | t100_mae | miss130 | miss100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M0 | 0 |  | 0.512 | 0.512 |  | 0.303 | 0 | 16.52 | 17.47 | -31.46 | -31.46 | 66.09 | -88.41 | -88.41 | 31.03 | 7.944 | 8.511 | 3 | 5 |
| M1 | 1 | 8.711 | 0.235 | 0.245 | -0.268 | 0.174 | 0 | 16.52 | 17.47 | -72.6 | -71.61 | 76.56 | -219.9 | -217.5 | 13.39 | 1.635 | 2.499 | 1 | 2 |
| M2 | 2 | 8.711 | 0.235 | 0.245 | 0 | 0.174 | 0.122 | 5.388 | 4.482 | -72.6 | -71.61 | 76.56 | -219.9 | -217.5 | 13.66 | 1.685 | 2.614 | 1 | 2 |
| M3 | 3 | 15.8,-2.294 | 0.23 | 0.24 | -0.004 | 0.168 | 0.122 | 5.388 | 4.482 | -72.69 | -70.7 | 76.56 | -221.2 | -216.6 | 13.03 | 1.379 | 2.428 | 1 | 2 |
| M4 | 5 | 14.33,-2.283,0.3464,-0.04108 | 0.229 | 0.296 | 0.028 | 0.242 | 0.122 | 5.388 | 4.482 | -69.52 | -65.54 | 74.94 | -212.8 | -203.5 | 15.29 | 1.547 | 2.95 | 2 | 4 |
| M5 | 7 | 0.5271,-0.7952,0.4337,7.374e-09,0.1709,4.923 | 0.229 | 0.689 | 0.197 | 0.62 | 0.122 | 5.388 | 4.482 | -61.57 | -55.59 | 69.02 | -191.6 | -178.2 | 21.26 | 2.945 | 3.74 | 3 | 6 |

fold constants: out/parsimony_constants.json
