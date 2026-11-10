Phase 3 / Task 1: position-dependent (developing vapour layer) film-boiling laws on the consistent
dataset (a07/a08).  All laws use film-temperature vapour properties and the lag-corrected G.

  h = C k_v (u_l rho_v cp_v / (k_v L_eff))^(1/3)        (local Leveque-type, D-free)
  L0  : L_eff = D                                (Phase-2 law, C = 8.55)
  L1  : L_eff = z            (inlet = leading edge, 1 constant)
  L2  : L_eff = z + z0       (z0 fitted, 2 constants)
  L3  : L_eff = z Lam/(z + Lam) (saturation at renewal length Lam, z0 = 0; 2 constants)
  L4  : L_eff = (z+z0) Lam/(z+z0+Lam) (3 constants)
  L5  : D-law x Darr entrance factor [1 + (D/(z+0.1D))^0.7]   (1 constant)
  Xm  : D-law x (1 - x_eq)^m  (liquid-core depletion; 2 constants)  [x built from measured q: endogenous]

updn_resid = mean ln(pred/meas) upstream minus downstream (0 = no station bias; negative = law under-predicts upstream relative to downstream)

| law | n_fit | params | loco_range | rms | loco | sd_cs | updn_resid_mean | updn_resid_sd | updn_by_case |
|---|---|---|---|---|---|---|---|---|---|
| L0 L=D | 1 | 8.546 | 8.27-8.8 | 0.235 | 0.246 | 0.167 | -0.115 | 0.095 | [np.float64(-0.273), np.float64(-0.216), np.float64(-0.059), np.float64(-0.111), np.float64(-0.085), np.float64(-0.025), np.float64(-0.068), np.float64(-0.092), np.float64(0.013), np.float64(-0.236)] |
| L1 L=z | 1 | 23.55 | 22.8-24.2 | 0.252 | 0.262 | 0.191 | 0.215 | 0.095 | [np.float64(0.057), np.float64(0.114), np.float64(0.271), np.float64(0.219), np.float64(0.245), np.float64(0.305), np.float64(0.262), np.float64(0.238), np.float64(0.343), np.float64(0.094)] |
| L2 L=z+z0 | 2 | 33.85,-0.7729 | 32.2-36.8;-0.962--0.526 | 0.228 | 0.24 | 0.156 | 0 | 0.095 | [np.float64(-0.157), np.float64(-0.101), np.float64(0.056), np.float64(0.004), np.float64(0.03), np.float64(0.09), np.float64(0.048), np.float64(0.024), np.float64(0.128), np.float64(-0.121)] |
| L3 L=z*Lam/(z+Lam) | 2 | 16.38,-2.045 | 15.6-17.1;-2.29--1.86 | 0.228 | 0.24 | 0.156 | 0 | 0.095 | [np.float64(-0.157), np.float64(-0.101), np.float64(0.056), np.float64(0.004), np.float64(0.03), np.float64(0.09), np.float64(0.048), np.float64(0.024), np.float64(0.128), np.float64(-0.121)] |
| Lz power (D/z)^a | 2 | 8.195,0.1163 | 7.91-8.49;0.0986-0.131 | 0.228 | 0.24 | 0.156 | -0 | 0.095 | [np.float64(-0.157), np.float64(-0.101), np.float64(0.056), np.float64(0.004), np.float64(0.03), np.float64(0.09), np.float64(0.048), np.float64(0.024), np.float64(0.128), np.float64(-0.121)] |
| L5 D-law x Darr entrance | 1 | 7.597 | 7.35-7.82 | 0.229 | 0.24 | 0.157 | -0.041 | 0.095 | [np.float64(-0.198), np.float64(-0.142), np.float64(0.015), np.float64(-0.037), np.float64(-0.011), np.float64(0.049), np.float64(0.007), np.float64(-0.018), np.float64(0.087), np.float64(-0.162)] |
| Xm D-law x (1-x)^m | 2 | 8.664,0.2387 | 8.37-8.92;-0.906-0.497 | 0.235 | 0.265 | 0.161 | -0.103 | 0.088 | [np.float64(-0.237), np.float64(-0.19), np.float64(-0.048), np.float64(-0.1), np.float64(-0.078), np.float64(-0.012), np.float64(-0.063), np.float64(-0.086), np.float64(0.018), np.float64(-0.23)] |
| L3x L3 with u_l=G(1-x)/rho_l [m=1/3, 0 extra] | 2 | 15.8,-2.294 | 14.9-16.5;-2.52--2.11 | 0.23 | 0.24 | 0.151 | 0 | 0.086 | [np.float64(-0.125), np.float64(-0.082), np.float64(0.054), np.float64(0.002), np.float64(0.022), np.float64(0.09), np.float64(0.037), np.float64(0.014), np.float64(0.118), np.float64(-0.13)] |
| L3m L3 x (1-x)^m | 3 | 16.55,-1.969,-0.1066 | 15.7-17.7;-2.21--1.38;-1.66-0.189 | 0.228 | 0.272 | 0.158 | -0 | 0.098 | [np.float64(-0.168), np.float64(-0.107), np.float64(0.057), np.float64(0.004), np.float64(0.033), np.float64(0.09), np.float64(0.051), np.float64(0.027), np.float64(0.131), np.float64(-0.119)] |
| P1 D-free sqrt(k rho cp u_l/(z+z0)) | 2 | 27.57,-0.02885 | 25.3-30.8;-0.184-0.209 | 0.333 | 0.357 | 0.273 | 0 | 0.093 | [np.float64(-0.159), np.float64(-0.103), np.float64(0.061), np.float64(0.016), np.float64(0.037), np.float64(0.089), np.float64(0.045), np.float64(0.014), np.float64(0.116), np.float64(-0.116)] |
| P0 D-free sqrt(k rho cp u_l/z) | 1 | 12.24 | 11.6-12.8 | 0.386 | 0.407 | 0.339 | 0.394 | 0.093 | [np.float64(0.235), np.float64(0.291), np.float64(0.455), np.float64(0.41), np.float64(0.431), np.float64(0.482), np.float64(0.438), np.float64(0.407), np.float64(0.509), np.float64(0.278)] |

## measured h_up/h_dn at matched times

| case | G | h_up_over_dn | x_up | x_dn | dT_up_minus_dn |
|---|---|---|---|---|---|
| 1 | 22.4 | 1.3 | 0.068 | 0.172 | -11.6 |
| 2 | 68.5 | 1.23 | 0.06 | 0.156 | -10.5 |
| 3 | 122 | 1.18 | 0.024 | 0.063 | -2.37 |
| 4 | 151 | 1.15 | 0.018 | 0.047 | -9.35 |
| 5 | 176 | 1.12 | 0.015 | 0.039 | -6.87 |
| 6 | 148 | 1.05 | 0.027 | 0.072 | -4.93 |
| 7 | 310 | 1.08 | 0.011 | 0.031 | -4.97 |
| 8 | 326 | 1.11 | 0.015 | 0.04 | -6.72 |
| 9 | 438 | 1.01 | 0.012 | 0.034 | -6.07 |
| 10 | 646 | 1.18 | 0.011 | 0.031 | -12 |
