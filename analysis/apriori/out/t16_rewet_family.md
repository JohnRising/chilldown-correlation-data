Phase 5: family of pressure-dependent, corresponding-states rewet rules; one constant each fitted on
the 11 UF low-g stations (leave-one-TEST-out); pressure slopes; predictions for NBS 9264 copper lines.

Every rule defines a liquid-side contact temperature T* and the wall temperature follows from the
semi-infinite contact relation  T_MFB = T* + beta (T* - T_l),  beta = e_l/e_w(T_MFB) (iterated),
T_l = Tsat (saturated inlet).  T_hn = Tc (0.905 + 0.095 Tr^8), Tr = Tsat/Tc (Lienhard 1976, reduced).
  R0  T* = Tc                                    (Phase-2 rule, 0 constants)
  R1a T* = T_hn + a (Tc - T_hn)                  (a = 1 -> R0)
  R1b T* = Tsat + c (T_hn - Tsat)
  R2  T* = Tsat + k (Tc - Tsat)                  (reduced linear)
  R3  Henry (1974) on T_B = T_hn:  T* = T_B + A (T_B - T_l) X^0.6, X = (e_l/e_w) h_fg/(c_w (T_B - Tsat))  (A fitted; Henry 0.42)
  R4  T* = Tsat + d Tc                           (constant reduced superheat; slope 1)
  R5  T* = T_hn + e Tc                           (constant reduced offset above T_hn)

data slope dT_MFB/dTsat (station dummy; test bootstrap): {"all": 0.17010939723346122, "all_ci": [-0.8270105565976967, 0.8754592947097288], "hiP": 1.224817987041536, "hiP_ci": [0.6134897630921762, 3.250571049552359]}

| rule | const | mean_err | MAE | rmse | loto_rmse | up_bias | dn_bias | slope_vs_Tsat | T_100kPa | T_300kPa | T_750kPa |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 |  | -4.574 | 5.623 | 7.595 | 7.595 | -8.403 | 0.021 | -0.111 | 130.7 | 129.4 | 128.2 |
| R1a | 1.377 | -0.082 | 5.087 | 6.122 | 6.662 | -3.912 | 4.514 | -0.14 | 135.5 | 134 | 132.4 |
| R1b | 1.571 | -0.852 | 6.111 | 7.188 | 7.68 | -4.373 | 3.374 | -0.607 | 140.9 | 133.8 | 127.4 |
| R2 | 1.109 | -0.301 | 5.161 | 6.243 | 6.694 | -4.072 | 4.224 | -0.229 | 136.4 | 133.9 | 131.3 |
| R3 | 0.426 | -0.392 | 5.647 | 6.619 | 7.126 | -4.043 | 3.989 | -0.415 | 138.8 | 134 | 129.6 |
| R4 | 0.325 | 0.018 | 6.352 | 7.741 | 8.717 | -4.342 | 5.249 | 0.957 | 122.2 | 132.4 | 143.5 |
| R5 | 0.122 | -0.001 | 5.025 | 5.935 | 6.465 | -3.826 | 4.588 | -0.039 | 134.6 | 133.8 | 133.8 |

## predictions, copper wall (NBS 9264 drive pressures)

| fluid | P_atm | Tsat | R0 | R1a | R1b | R2 | R3 | R4 | R5 | beta_Cu_R0 |
|---|---|---|---|---|---|---|---|---|---|---|
| Nitrogen | 1 | 77.35 | 126.9 | 131.4 | 136.5 | 132.3 | 122.3 | 119 | 130.6 | 0.01 |
| Nitrogen | 2.5 | 86.07 | 126.7 | 131.1 | 131.8 | 131.2 | 121.2 | 127.7 | 130.8 | 0.01 |
| Nitrogen | 3.4 | 89.47 | 126.7 | 131 | 130.1 | 130.7 | 120.9 | 131 | 130.9 | 0.01 |
| Nitrogen | 4.2 | 91.98 | 126.6 | 130.8 | 128.9 | 130.4 | 120.7 | 133.5 | 131.1 | 0.01 |
| Nitrogen | 5.9 | 96.33 | 126.6 | 130.6 | 127 | 129.8 | 120.5 | 137.9 | 131.4 | 0.01 |
| Hydrogen | 1 | 20.37 | 33.31 | 34.48 | 35.78 | 34.71 | 43.39 | 31.29 | 34.28 | 0.01 |
| Hydrogen | 2.5 | 23.92 | 33.27 | 34.38 | 33.97 | 34.29 | 42.15 | 34.84 | 34.41 | 0.01 |
| Hydrogen | 4.2 | 26.39 | 33.24 | 34.25 | 32.95 | 33.99 | 40.95 | 37.31 | 34.67 | 0.01 |
| Hydrogen | 5.1 | 27.42 | 33.23 | 34.17 | 32.63 | 33.86 | 40.36 | 38.34 | 34.84 | 0.02 |
| Hydrogen | 7.6 | 29.73 | 33.2 | 33.9 | 32.26 | 33.58 | 38.82 | 40.66 | 35.45 | 0.02 |
| Hydrogen | 11 | 32.12 | 33.17 | 33.44 | 32.64 | 33.28 | 36.72 | 43.11 | 36.59 | 0.02 |
