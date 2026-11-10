# Hartwig et al. 2022 (npj Microgravity 8:33) — LN2 chilldown in microgravity, test F2 (bare tube, continuous flow)

**Citation:** Hartwig, J., Chung, J.N., Dong, J., Han, B., Wang, H., Darr, S., Taliaferro, M., Jain, S., Doherty, M. (journal author list), "Nitrogen flow boiling and chilldown experiments in microgravity using pulse flow and low-thermally conductive coatings", *npj Microgravity* 8, 33 (2022). doi:10.1038/s41526-022-00220-9. CC-BY 4.0.

## Data availability (checked 2026-11-06)
The Data Availability statement says only "Supplementary information accompanies the paper". The single ESM file is the Nature **Reporting Summary** (`lit/raw/hartwig2022_npj_MOESM1_ESM.pdf`); it contains **no data**. Zenodo, figshare and data.nasa.gov were searched with no hits. So there is **no public numeric time series**, unlike Darr 2016 npj (MOESM26 xls), and the curves had to be digitized.

## Files
- `F2_bare_continuous_0p05g.csv` – **main file**, 0.1 s grid, 0.1–19.8 s: average exit outer-wall T (Fig 2a, plus an independent Fig 1a copy as a check), exit pressure (Fig 2b), mass flux G (Fig 2c).
- `F2_fig*_raw_columns.csv` – one row per image column (≈0.033 s), with the min/max of the curve pixels in each column.
- `hartwig2022_flight_test_matrix_ln2_microgravity.csv` – Table 1 (all 28 flight tests), copied from OpenFLUME.
- `qc_overlay.png` – digitized curves overlaid on the source raster.
- Script: `analysis/digitize_hartwig2022.py` (`python analysis/digitize_hartwig2022.py write`).

## Test F2 conditions
- LN2, 0.05 g (ZeroG parabolic flight), horizontal, bare SS304 tube: OD 12.70 mm, wall 0.51 mm, ID 11.68 mm. This is the same tube cross-section as Darr 2016.
- Inlet pressure 550 kPa (Table 1, time-averaged); the inlet was subcooled by a precooler, but its temperature is not reported. Steady-state G 533 kg/m²/s, Re 67 324; chilldown time 14.6 s; LN2 consumed 0.872 kg.
- Wall TCs: stations at 4/6/8/16/24 in (0.102/0.152/0.203/0.406/0.610 m) from flange A (Fig 5d), 3 TCs per station. The "exit" curve is the mean of TC5/10/15 at 0.610 m.
- Vacuum jacket at ~1 Pa. Initial wall T ≈ 288 K.
- Geometry inconsistency: the text gives 0.914 m test sections, while Fig 5 shows 26.6 in (0.676 m) inside the vacuum and labels it a "0.61 m test section".

## Digitization error
- Fig 1a and Fig 2a show the same data. The two digitized copies agree to 0.69 K rms (max 2.5 K) outside the quench knee and 4.4 K rms (max 8.7 K) across the 13–14 s knee, which corresponds to a time shift of about 0.1 s.
- Pixel resolution: 0.62 K, 1.7 kPa and 2.6 kg/m²/s per pixel; 0.033 s per pixel.
- Quoted uncertainty: T ±1.5 K outside the knee (instrument ±2–3 K), knee timing ±0.15 s, P ±5 kPa, G ±5 kg/m²/s (instrument 2.8%).

## Missing / caveats
- **Only the exit-station average** was digitized. Fig 1c has all 15 TCs (15 overlapping coloured curves; could be digitized for quench-front timing at each station).
- No inlet temperature history and **no inlet pressure history** (only exit P(t)). G(t) comes from gas flow meters downstream of the vaporizers, so it is lagged and smoothed relative to the test-section flow. The very high initial G (~900–1000) is likely a startup/meter artefact.
- Other bare-tube continuous-flow tests (F1 at 0.38 g, F5/F8/F10 at low Re, F12) appear only in Table 1, with no curves. The pulse-flow tests F3/F4 (Fig 3) could be digitized the same way, but they are not continuous flow.
