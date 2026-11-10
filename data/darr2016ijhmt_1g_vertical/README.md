# Darr et al. 2016 IJHMT (terrestrial LN2 chilldown, Part I) — Tests 31 & 40, 1 g vertical upflow

**Data:** Darr, S.R., Hu, H., Glikin, N.G., Hartwig, J.W., Majumdar, A.K., Leclair, A.C., Chung, J.N., *Int. J. Heat Mass Transfer* 103 (2016) 1225–1242, doi:10.1016/j.ijheatmasstransfer.2016.05.019. Part II (flow direction) is doi:10.1016/j.ijheatmasstransfer.2016.08.044.
**Figures + BCs taken from:** Kassemi, Hylton, Kartuzova, Ahuja, Hosangadi, "Predictions of Line Chilldown Boiling Regime Transitions by a Coupled CFD-Subgrid Boiling Model Validated against 1G LN2 Experiments", NASA NTRS 20240013323 (2024 preprint), Figs 3a/b and 9a/b, plus the text of §5 (`lit/raw/ntrs/20240013323.pdf`).

## Data availability
No supplementary data or repository was found for Darr IJHMT 2016 I/II (Elsevier; Crossref has no dataset relations). The curves here are digitized from Kassemi et al.'s Excel replot of the measured data (embedded 600-px rasters, extracted without loss).

## Files
- `test31_LN2_1g_upflow.csv`, `test40_LN2_1g_upflow.csv`: `t_s, Twall_TC1-3_avg_K, Twall_TC4-6_avg_K`, 0.05 s grid.
- `qc_overlay.png`; script `analysis/digitize_kassemi2024_darr_tests.py`.

## Boundary conditions (Kassemi et al., time-averaged)
| | Test 31 | Test 40 |
|---|---|---|
| T_wall,0 | 289.9 K | 288.54 K |
| T_in | 86.96 K | 88.61 K |
| P_out | 314.6 kPa | 325.9 kPa |
| Tsat(P_in) | 88.44 K (quoted for both tests) | 88.44 K |
| mdot | 0.0334 kg/s (G = 312 kg/m²/s) | 0.0346 kg/s (G = 323) |
| Re_in | 33 940 | 35 626 |

- Geometry: the same UF tube as `data/uf2016`, i.e. 304SS, L = 0.572 m, OD 12.70 / ID 11.68 mm. TC1-3 at 0.149 m and TC4-6 at 0.401 m (outer wall, 3 TCs averaged); PT1 at 0.180 m, PT2 at 0.371 m.
- Vacuum chamber: 293 K, 1.33 Pa, ε = 0.453. About 2.54 cm of tube sticks out of each side of the chamber.

## Digitization error
- Event times read from the digitized curves were compared with Kassemi Tables 2–3 (12 events): mean +0.08 s, rms 0.11 s, max 0.18 s. The table values are rounded to 0.1 s.
- Pixel resolution: 0.46 K and 0.055 s per pixel. Film-boiling values are good to about ±1 K. The first ~0.3 s reads up to 3 K low because the plot markers are thick at the plot edge; use T_wall,0 from the table above.

## Missing
- **Only mean mdot and P_out** are available, with no time histories (the original paper has them in figures). The inlet pressure is not given separately.
- No other Darr IJHMT tests have been digitized. The original paper's figures (for example its mass-flux sweep) would need the article PDF, which is paywalled and not in `lit/raw`.
