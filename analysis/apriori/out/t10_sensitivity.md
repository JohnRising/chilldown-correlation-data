Phase 2 / Task 4: full gravity-free boiling curve + 0-D surrogate sensitivity (tornado).

Wall: m''_w cp(T) dT/dt = -q''(T) + q''_par(T),  q''_par = 98.6 + 2.2 (293 - T) W/m2 (a07),
compared with the station-mean TC history T_m (what a 1-D model is validated on).  Inputs at each
station: measured P(t), Tsat(t), lag-corrected G(t) (tau = 2.5 s).
Baseline boiling curve (all regimes as functions of local node variables):
  film       : q_FB = h_FB (T - Tsat), h_FB = C (k_vf/D) Pe_l^(1/3)  (a08, C = 8.55)
  rewet      : T_MFB = T_c + beta (T_c - Tsat), beta = e_l/e_w(T_MFB)       (a09, 0 constants)
  nucleate   : q_NB = h_DB,l dT + (A dT)^3,  A = 55 pr^0.12 (-log10 pr)^-0.55 M^-0.5  (Cooper 1984
               solved for q; Dittus-Boelter liquid floor h_DB,l = 0.023 Re_l^0.8 Pr_l^0.4 k_l/D, Re_l = G D/mu_l)
  CHF        : q_CHF = Katto-Ohno (1984) saturated, L = z; dT_CHF = root of q_NB(dT) = q_CHF
  transition : q_TB = q_FB(T_MFB) + (q_CHF - q_FB(T_MFB)) theta^n, theta = (T_MFB-T)/(T_MFB-T_CHF), n = 0.5

| variant | rms_mean | rms_med | t130_mae | t130_bias | t100_mae | t100_bias | miss130 | miss100 |
|---|---|---|---|---|---|---|---|---|
| baseline | 13.1 | 9.3 | 1.78 | 0.542 | 1.62 | -0.283 | 1 | 2 |
| film h x0.7 | 19.8 | 20.3 |  |  |  |  | 11 | 8 |
| film h x1.3 | 28.1 | 26.5 | 4.12 | -3.99 | 4.92 | -4.68 | 3 | 6 |
| T_MFB -5 K | 13 | 9.76 | 1.78 | 0.542 | 1.19 | 0.319 | 1 | 4 |
| T_MFB +5 K | 13.8 | 11 | 1.69 | -0.516 | 2.4 | -1.57 | 1 | 2 |
| CHF x0.5 | 12.9 | 9.53 | 1.78 | 0.542 | 1.66 | 0.14 | 1 | 2 |
| CHF x1.5 | 13.1 | 9.22 | 1.78 | 0.542 | 1.68 | -0.431 | 1 | 2 |
| transition n=0.25 | 13.1 | 9.2 | 1.78 | 0.542 | 1.67 | -0.426 | 1 | 2 |
| transition n=1 | 13 | 9.62 | 1.78 | 0.542 | 1.66 | 0.206 | 1 | 2 |
| NB x0.5 | 13 | 9.29 | 1.78 | 0.542 | 1.63 | -0.294 | 1 | 2 |
| NB x2 | 13.1 | 9.31 | 1.78 | 0.542 | 1.62 | -0.274 | 1 | 2 |
| CHF = 100 kW/m2 const (meas.) | 13 | 9.49 | 1.78 | 0.542 | 1.65 | -0.05 | 1 | 2 |

## baseline per station

| case | st | rms | dt130 | miss130 | dt100 | miss100 |
|---|---|---|---|---|---|---|
| 1 | up | 20.2 |  | 0 |  | 0 |
| 1 | dn | 29.6 |  | 0 |  | 0 |
| 2 | up | 9.04 |  | 0 |  | 0 |
| 2 | dn | 19.1 |  | 0 |  | 0 |
| 3 | up | 8.28 |  | 1 |  | 0 |
| 3 | dn | 3.83 |  | 0 |  | 0 |
| 4 | up | 3.44 |  | 0 |  | 0 |
| 4 | dn | 3.54 |  | 0 |  | 0 |
| 5 | up | 23.9 |  | 0 |  | 0 |
| 5 | dn | 28.6 |  | 0 |  | 0 |
| 6 | up | 13.4 | 1.7 | 0 |  | 1 |
| 6 | dn | 22.8 | -3.54 | 0 |  | 1 |
| 7 | up | 6.42 | 0.66 | 0 | -2.05 | 0 |
| 7 | dn | 12.2 | -1.74 | 0 | -3.82 | 0 |
| 8 | up | 7.32 | 1.56 | 0 | 0 | 0 |
| 8 | dn | 8.25 | -0.91 | 0 | -1.63 | 0 |
| 9 | up | 9.57 | 2.54 | 0 | 0.43 | 0 |
| 9 | dn | 2.66 | 0.4 | 0 | -0.13 | 0 |
| 10 | up | 22.6 | 4.21 | 0 | 4.54 | 0 |
| 10 | dn | 6.47 | 0.54 | 0 | 0.4 | 0 |
