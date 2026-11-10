Phase 2 / Task 1: ONE consistent heat-flux dataset from a lumped wall energy balance on the
station-mean outer-wall TC history, for all 10 cases; plus tests of the 0.75 anomaly.

Wall energy balance per unit INNER area (thin wall, Bi < 0.02):
    q''_E(t) = -m''_w cp(T_m) dT_m/dt + q''_par(T_m)
    m''_w = rho_w (Do^2 - Di^2)/(4 Di) = 4.20 kg/m2, cp = NIST 304SS cp(T)
    T_m   = mean of the three station TCs
    dT/dt : Savitzky-Golay, 2nd order, 9 points (0.60 s) [film boiling]; 5 points (0.33 s) used for
            CHF peaks (a01). Sensitivity to 5/9/15 points reported.
    q''_par(T) = a + b (293 - T): parasitic heat leak into the wall (radiation from the 293 K vacuum
            chamber + residual-gas conduction), fitted to (q_inverse - q_lump) of the five
            self-consistent cases (4,7,8,9,10) -- i.e. the authors' own parasitic model; compared with
            the radiation-only estimate eps sigma (T_ch^4 - T^4) Do/Di, eps = 0.453 (Kassemi et al.).
Inner-surface temperature (quasi-steady radial profile with uniform storage sink, adiabatic outer
surface):  T_i = T_m - q''_E r_i/k(T) [ r_o^2 ln(r_o/r_i)/(r_o^2 - r_i^2) - 1/2 ]  (= q'' * 0.254 mm / k,
0.2 K in film boiling, ~2.5 K at the CHF peak).  Radial diffusion time delta^2/alpha ~ 0.07 s
(one sample) -> no time-shift correction applied.

parasitic fit: q_par = 98.6 + 2.200 (293 - T) W/m2 (radiation-only eps=0.453: 97/177/198 W/m2 at 250/180/130 K; fit 193/347/457)

smoothing sensitivity (film boiling h_E): win 5: median h ratio to win9 0.9995, rms ln ratio 0.045; win 9: median h ratio to win9 1.0000, rms ln ratio 0.000; win 15: median h ratio to win9 1.0007, rms ln ratio 0.054

## consistent vs published (film boiling, t>1 s)

| case | st | n | hE_med | hpub_med | ratio_pub_E | Tw_minus_Tipub_med | phase1_k |
|---|---|---|---|---|---|---|---|
| 1 | up | 249 | 48.38 | 34.39 | 0.731 | -0.042 | 0.735 |
| 1 | dn | 249 | 35.1 | 26.28 | 0.728 | -0.033 | 0.735 |
| 2 | up | 194 | 85.66 | 63.33 | 0.732 | -0.067 | 0.732 |
| 2 | dn | 194 | 69.67 | 50.85 | 0.732 | -0.056 | 0.732 |
| 3 | up | 233 | 91.77 | 67.6 | 0.736 | -0.064 | 0.746 |
| 3 | dn | 312 | 79.01 | 57.87 | 0.737 | -0.054 | 0.746 |
| 4 | up | 260 | 80.99 | 81.19 | 0.995 | -0.008 | 1 |
| 4 | dn | 260 | 72.73 | 71.86 | 0.995 | -0.006 | 1 |
| 5 | up | 262 | 68.91 | 49.17 | 0.737 | -0.054 | 0.741 |
| 5 | dn | 262 | 62.95 | 44.87 | 0.737 | -0.048 | 0.741 |
| 6 | up | 195 | 112.9 | 82.92 | 0.735 | -0.076 | 0.739 |
| 6 | dn | 247 | 106 | 77.06 | 0.735 | -0.072 | 0.739 |
| 7 | up | 197 | 130.4 | 131.6 | 0.997 | -0.008 | 1 |
| 7 | dn | 236 | 123.4 | 123.8 | 0.996 | -0.007 | 1 |
| 8 | up | 186 | 157.5 | 156.6 | 0.996 | -0.009 | 1 |
| 8 | dn | 215 | 142.3 | 142.6 | 0.996 | -0.007 | 1 |
| 9 | up | 152 | 167.5 | 170.4 | 0.996 | -0.007 | 1 |
| 9 | dn | 189 | 173.5 | 174.2 | 0.996 | -0.007 | 1 |
| 10 | up | 72 | 272.5 | 269.7 | 0.996 | -0.008 | 1 |
| 10 | dn | 118 | 205.7 | 206.7 | 0.996 | -0.01 | 1 |

## anomaly: decomposition q_pub = alpha q_store + beta q_par

| case | st | alpha_store | beta_par | rms_storefit | k_total | rms_totalfit | Ti_pub_minus_Tm |
|---|---|---|---|---|---|---|---|
| 1 | up | 0.751 | -0.039 | 582.4 | 0.735 | 586.1 | -0.099 |
| 1 | dn | 0.74 | 0.335 | 425.6 | 0.731 | 426.5 | -0.079 |
| 2 | up | 0.728 | 1 | 483.3 | 0.731 | 484.7 | -0.165 |
| 2 | dn | 0.727 | 0.982 | 394.4 | 0.731 | 395.6 | -0.14 |
| 3 | up | 0.731 | 0.928 | 414.9 | 0.736 | 415.6 | -0.154 |
| 3 | dn | 0.735 | 0.805 | 420.5 | 0.737 | 420.6 | -0.12 |
| 4 | up | 1.001 | 0.802 | 547.3 | 0.996 | 548 | -0.181 |
| 4 | dn | 1 | 0.831 | 543.6 | 0.996 | 543.9 | -0.171 |
| 5 | up | 0.787 | -0.972 | 892.6 | 0.746 | 914.5 | -0.127 |
| 5 | dn | 0.813 | -1.833 | 856.1 | 0.748 | 894.3 | -0.114 |
| 6 | up | 0.747 | 0.242 | 772.2 | 0.738 | 775.7 | -0.196 |
| 6 | dn | 0.749 | 0.271 | 777.2 | 0.739 | 781 | -0.179 |
| 7 | up | 1.001 | 0.807 | 489.2 | 0.998 | 490.6 | -0.258 |
| 7 | dn | 0.999 | 0.882 | 493.7 | 0.996 | 494.3 | -0.24 |
| 8 | up | 0.995 | 1.049 | 528.2 | 0.996 | 528.2 | -0.338 |
| 8 | dn | 0.996 | 0.994 | 542.6 | 0.996 | 542.6 | -0.303 |
| 9 | up | 0.997 | 0.97 | 737.1 | 0.997 | 737.1 | -0.344 |
| 9 | dn | 0.997 | 0.936 | 741.8 | 0.996 | 741.9 | -0.326 |
| 10 | up | 1.002 | 0.458 | 688.1 | 0.998 | 696.7 | -0.462 |
| 10 | dn | 0.999 | 0.78 | 528.9 | 0.997 | 531.2 | -0.334 |

storage factor bad/good = 0.7518 -> equivalent wall thickness 0.388 mm at the same OD

## hold-out arbiter (law fitted on cases 4,7-10, predicting 1,2,3,5,6)

| model | case | G | P_kPa | ratio_hE_over_pred | ratio_hpub_over_pred |
|---|---|---|---|---|---|
| free lnG,lnP,dn | 1 | 22.44 | 176.8 | 0.556 | 0.402 |
| free lnG,lnP,dn | 2 | 68.5 | 233 | 0.817 | 0.598 |
| free lnG,lnP,dn | 3 | 122.7 | 92.04 | 1.435 | 1.055 |
| free lnG,lnP,dn | 5 | 176.5 | 121.7 | 0.92 | 0.663 |
| free lnG,lnP,dn | 6 | 143.7 | 290.1 | 0.926 | 0.672 |
| F2 form sqrt(k_v cp_v G/D) | 1 | 22.44 | 176.8 | 1.02 | 0.737 |
| F2 form sqrt(k_v cp_v G/D) | 2 | 68.5 | 233 | 1.332 | 0.974 |
| F2 form sqrt(k_v cp_v G/D) | 3 | 122.7 | 92.04 | 1.124 | 0.826 |
| F2 form sqrt(k_v cp_v G/D) | 5 | 176.5 | 121.7 | 0.731 | 0.527 |
| F2 form sqrt(k_v cp_v G/D) | 6 | 143.7 | 290.1 | 1.2 | 0.871 |
