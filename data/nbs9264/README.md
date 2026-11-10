# NBS Report 9264 (Brennan et al. 1966) — long copper transfer-line cooldown, LH2 & LN2

Copied verbatim (headers intact) from `OpenFLUME/validation/data/digitized/chilldown/` (read-only source) on 2026-11-06.

**Citation:** Brennan, J.A., Brentari, E.G., Smith, R.V., Steward, W.G., *Cooldown of Cryogenic Transfer Lines – An Experimental Report*, NBS Report 9264 (NASA-CR-81338), 7 Nov 1966. NTRS 19670007291.

## Files
- `nbs9264_runs_metadata.csv` – one row per digitized run (figure, fluid, liquid state, drive pressure, inlet T, geometry, TC stations, time span).
- `nbs9264_figNN_<fluid>_<P>atm.csv` – long format `time_s, station_m, T_wall_K` (hand-clicked markers, linearly resampled; per-file header gives method/uncertainty).
- `brennan1966_nbs9264_inlet_restriction_results.csv` – Table 1 (peak line pressure with inlet restrictors); no transients.

## Boundary conditions (all runs)
| item | value |
|---|---|
| line | 61.0 m (200 ft) horizontal, vacuum-jacketed copper, ID 15.9 mm, OD 19.0 mm (wall 1.55 mm) |
| initial wall T | ~ambient 270–290 K (first digitized point per station; e.g. Fig 2: 278.8 K, Fig 10: 269.3 K) |
| inlet | 300 L supply dewar pressurized to P_drive (LH2: 2.5–11.0 atm; LN2: 2.5–5.9 atm abs.); "saturated" = equilibrated at P_drive (T_in = Tsat(P_drive): LH2 27.4 K @5.1 atm; LN2 86.1/89.5/96.3 K); "subcooled" = saturated at ~0.82 atm then pressurized (T_in ≈ 19.5 K H2 / 76 K N2) |
| outlet | discharge to atmosphere, ~0.82 atm (Boulder) |
| wall TCs | Cu-constantan on outer wall at 6.1 / 24.4 / 43.0 / 60.4 m from dewar (20/80/141/198 ft) per metadata |
| runs | LH2: Figs 2 (sat, 5.1 atm), 3–7 (subcooled 2.5/4.2/5.9/7.6/11.0 atm); LN2: Figs 10–12 (sat 2.5/3.4/5.9 atm), 13–14 (subcooled 4.2/5.9 atm) |
| durations | LH2 28–150 s, LN2 90–250 s |

## Missing / caveats for 1-D validation
- **Mass flow rate history not digitized** (report has flow/pressure oscillograph figures among Figs 2–51; these would need digitizing). Only P_drive and outlet ≈ atmospheric known → simulation must predict the flow (pressure-driven BC), which is actually a stringent test.
- Inlet valve/restriction geometry and the dewar-to-line plumbing losses are not fully specified; valve type "did not affect T histories" (report sec. 3).
- Station-location inconsistency between files: `runs_metadata` uses 20/80/141/198 ft, while the header of `brennan1966_nbs9264_inlet_restriction_results.csv` says ~20/60/100/140 ft. Check report Fig. 1 before use.
- Low-T TC readings unreliable (curve crossings near the cold end are measurement error, report sec. 3.1). Digitization ±0.5–1 K, ±0.5 s; steep knees ±2–5 K.
- 1 g only, horizontal; large L/D (~3800) → good for quench-front propagation and pressure-driven flow, not for gravity effects.
