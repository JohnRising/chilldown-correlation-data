# A-priori analysis of the Darr et al. (2016) low-g LN2 chilldown data

Data: `data/uf2016/case{1..10}.csv`. These are the 304SS tube results from the 2016 npj Microgravity paper (2:16033): ID 11.68 mm, OD 12.70 mm, TC stations at z = 0.149 m and z = 0.401 m.

Run the scripts in order, from `analysis/apriori/`, with `../../.venv/bin/python`. The total run time is about 4 min.

| script | content |
|---|---|
| `common.py` | geometry, N2 properties (CoolProp, cached), NIST 304SS k(T), cp(T), superheat limits, statistics |
| `prep.py` | valid time window, per-TC lumped q'', quench/MFB detection, film-boiling mask, energy-consistent q'' |
| `a01_characterize.py` | regimes, boiling curves, MFB/CHF table, inverse-vs-energy-balance check, case-7 inlet probe |
| `a02_lag.py` | flow-meter lag |
| `fb_models.py`, `a03_film_boiling.py` | film-boiling dataset and candidate closures, fits, LOCO-CV, residual trends |
| `a04_rewet.py` | T_MFB compared with homogeneous-nucleation/contact theory and published correlations |
| `a05_chf_transition.py` | CHF (Katto-Ohno, Zuber) and transition-curve shape |
| `a06_surrogate.py` | 0-D lumped-wall surrogate for every closure family |

Figures are in `figs/fig01..fig12`. Tables are in `out/t01..t06*.md`, with CSV versions in `out/`.

---------------------------------------------------------------------------------------------------

## 0. Data-reduction issue: the published q'' is inconsistent in 5 of the 10 cases

The wall is thin, with Biot number below 0.02, so q'' can be checked against the wall energy balance:

    q''_lump = -m''_w cp(T) dTi/dt,   m''_w = ρ_w (Do² - Di²)/(4 Di) = 4.20 kg/m²

Here cp is the NIST 304SS value, so m''_w·cp is 2.0 kJ/m²K at 290 K and 1.1 kJ/m²K at 120 K.

- **Cases 4 and 7-10:** inverse q''/q''_lump = 1.016-1.027 at all times. The remaining +0.3 kW/m² is a plausible parasitic heat leak.
- **Cases 1, 2, 3, 5 and 6:** the ratio is 0.746-0.761 (`out/q_consistency.csv`). It stays at that value over the whole run, from 290 K to 150 K and at both stations.

A constant ratio like this is a scaling factor in the data reduction for those runs, such as a wall-property or geometry factor. It is not a physical heat loss.

**Fix used:** for those 5 cases I divide the tabulated q'' and h by k = 0.732-0.746 (`prep.inverse_scale`). I call the result the *energy-consistent* q''_ec and h_ec, and all primary results use it. Results with the published h are also reported. The Darr 2016 correlation was fitted to the published h, so this issue goes into its constants.

## 1. Characterization (fig01-fig05, `out/t01_characterization.md`)

**End of the usable record.** It is cut where P1 drops below 70 % of its median, which marks valve closure or venting: case 2 after 18.7 s, case 4 after 18.4 s, case 10 after 13.8 s.

**Film boiling covers most of each run.** h_ec ranges from 40 W/m²K (case 1) to 280 W/m²K (case 10). h is almost independent of superheat (fig03). Pressure and flow oscillations with a 1.3-2 s period (cases 1, 5, 6) cause ±30 % swings in h.

**Quenching stations.** Rewet is reached at both stations in cases 6-10, and only upstream in case 3.

**Strong asynchrony between TCs.** The 3rd TC of each station (TC3/TC6) is always the last to quench, 2-3 s after the others (cases 7, 8, 9 upstream; case 7 downstream). In cases 3 and 6 it does not quench at all within the record. The ±0.02 g residual gravity or the ~10° skew of the TCs are likely causes.

The inverse-method Ti and q'' average over the three TCs, which gives double or smeared q'' peaks (fig05). For that reason, MFB and CHF are taken from per-TC lumped q'' = -m''_w cp dT_TC/dt.

**T_MFB definition.** r = -dT/dt, from a Savitzky-Golay filter with 5 points, 2nd order:
- Quench: the maximum of r, required to exceed 4× the film-boiling baseline r_fb, to exceed 10 K/s, and to come with a ≥15 K drop within ±1 s. r_fb is the median of r over [t_q-3, t_q-1] s.
- T_MFB: the last T before the quench at which r ≤ 2·r_fb, i.e. the point where q'' has doubled relative to the film-boiling trend.
- Sensitivity: thresholds of 1.5× and 3×, and the point of maximum curvature. The 1.5×/2×/3× thresholds agree within 1-3 K. The curvature point sits 3-25 K lower because it falls inside the drop.
- Station value: the median over the TCs that quench. The spread between TCs is 1-12 K.

**Measured T_MFB, compared with predictions** (`out/t04_rewet.md`; columns from "contact T_hn" onward are predictions):

| case | st | P kPa | Tsat | G_lag | T_MFB | TC range | ±meth | dT_MFB | q_MFB kW/m² | contact T_hn | contact T_sp,EOS | Henry/T_hn | Darr-H We | Darr19 full |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | up | 93 | 76.7 | 127 | 136.3 | 136.0-136.7 | 14.4 | 59.7 | 7.7 | 118.1 | 123.8 | 140.4 | 119.1 | 124.2 |
| 6 | up | 265 | 86.5 | 95 | 144.4 | 143.1-145.7 | 4.4 | 57.8 | 8.9 | 117.4 | 122.7 | 135.4 | 118.5 | 118.2 |
| 6 | dn | 220 | 84.6 | 101 | 130.7 | 126.9-134.5 | 5.9 | 46.1 | 7.2 | 117.5 | 122.9 | 136.4 | 118.5 | 119.4 |
| 7 | up | 293 | 87.6 | 209 | 134.3 | 129.3-136.3 | 10.9 | 46.7 | 6.2 | 117.3 | 122.6 | 134.9 | 123.2 | 118.6 |
| 7 | dn | 289 | 87.5 | 191 | 128.5 | 125.0-130.7 | 12.6 | 41.0 | 6.6 | 117.3 | 122.6 | 135.0 | 122.6 | 118.5 |
| 8 | up | 382 | 90.7 | 293 | 129.0 | 124.4-131.6 | 6.1 | 38.3 | 6.2 | 117.2 | 122.3 | 133.5 | 126.3 | 117.1 |
| 8 | dn | 370 | 90.3 | 269 | 126.3 | 125.9-129.0 | 3.0 | 35.9 | 5.9 | 117.2 | 122.3 | 133.7 | 125.5 | 117.2 |
| 9 | up | 450 | 92.7 | 298 | 140.1 | 132.6-142.1 | 12.4 | 47.4 | 9.1 | 117.2 | 122.1 | 132.6 | 126.7 | 115.9 |
| 9 | dn | 432 | 92.2 | 294 | 127.3 | 125.8-129.5 | 4.1 | 35.1 | 6.6 | 117.2 | 122.2 | 132.8 | 126.5 | 116.2 |
| 10 | up | 762 | 99.7 | 411 | 142.5 | 141.5-143.7 | 3.0 | 42.8 | 12.5 | 117.3 | 121.6 | 129.7 | 131.4 | 111.5 |
| 10 | dn | 573 | 95.8 | 238 | 133.3 | 132.5-134.5 | 2.2 | 37.5 | 8.7 | 117.2 | 121.9 | 131.3 | 125.5 | 113.6 |

Notes on the columns:
- "±meth" is half the spread across the four definitions. It is large where the curvature point falls deep inside a slow quench.
- q_MFB is the q''_ec median over the 1 s before the earliest TC onset, i.e. the minimum film-boiling flux. It is 6-12.5 kW/m².

**CHF.** The per-TC lumped peak q'' is 75-123 kW/m² (station median), at ΔT_CHF = 19-37 K (mean 24 K). The inverse method reports only 34-116 kW/m², because averaging and the inverse solution smear the peak.

**Nucleate tail.** After the quench, h settles at 80-550 W/m²K at ΔT of 2.5-5.5 K.

**Case-7 inlet probe (fig04).** The probe reads 2-3 K above Tsat(P_inlet) for the whole run: +3.0 K for 0-2 s, +2.6 K for 2-10 s, +2.0 K after 20 s. P_inlet is about 1-3 kPa below P1, which is not physically possible for flow in that direction, so there is a pressure-transducer offset of order 2 kPa. That offset accounts for only 0.2 K of Tsat. Either the probe sees superheated vapour in a two-phase inflow, or it has a ~2.5 K bias. In either case the inflow is not subcooled; x_in = 0 is assumed below.

## 2. Mass-flux lag (fig06, fig07, `out/t02_lag.md`)

- **Fast fluctuations (high-pass xcorr).** The meter's G fluctuations are in phase with P1 at the test section (τ_P = -0.3 to +0.6 s in all cases). So the fluctuating part of the meter signal is pressure-driven and is not a transported signal. Cross-correlating h with G is ambiguous because the oscillations are periodic, which aliases the peak at multiples of the 1.3-2 s period.
- **Slow (bulk) lag.** I fitted a pooled regression ln h = a + b ln G(t+τ) + c ln ΔT + d·dn over all film-boiling samples. RMS(ln h) falls from 0.336 at τ = 0 to a clear minimum of 0.305 at **τ = 2.5 s** for a pure delay. A first-order-lag model is worse.
- **Per case.** With the pooled coefficients held fixed, τ = 4.5, 1.75, 4.75, 1.5, 0.25, 3.0, ≥6, 4.0, 3.25, 2.5 s for cases 1-10 (median 3.1 s). These are noisy because within a single case the error surface is flat.
- **Case 10, as in the paper:** the early h peak is at 1.53 s and the G peak at 4.53 s, so **τ = 2.9 s**, which matches the paper's "≈3 s".
- **Choice:** a global delay, G_ts(t) = G_meas(t + 2.5 s), held at the last value at the end of the record. Results with uncorrected G are also given. The lag correction lowers the free-regression RMS from 0.283 to 0.256 and the LOCO RMS from 0.376 to 0.332.

## 3. Film-boiling scaling (fig08, fig09, `out/t03_film_boiling.md`)

**Dataset.** 4422 samples (20 case-stations). A sample is included if t > 0.5 s, it is before the earliest MFB onset at that station minus 0.5 s, and Ti > max(T_MFB) + 15 K. Case-stations are weighted equally in the fits.

**Free power law** (h_ec, lag-corrected G; 95 % case-bootstrap CI):

    h ∝ G^0.35 [0.08, 0.48] · ΔT^-0.01 [-0.22, 0.28] · P^0.37 [0.15, 0.67] · (0.89 [0.85, 0.94] downstream)
    RMS(ln h) = 0.256, LOCO 0.332; irreducible scatter about per-run trends = 0.174

G and P are correlated across cases (r = 0.65 in log), so their separate exponents are only loosely identified. Their combined effect is well constrained.

**The data rule out** a ΔT^-0.25..-0.5 dependence with constant h_fg, and they rule out G^0.8 turbulent vapour convection.

Closure results (h_ec, lag-corrected G). bias is the geometric-mean ratio minus 1. Trend columns are the slopes of ln(pred/meas) against each variable; "per10s" is the within-run time slope.

| closure (props) | n fit | constants | bias | RMS ln | ±20 % | ±30 % | LOCO RMS | trend lnG | lnP | lnΔT | dn |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Bromley FC, C = 2.7 (film T) | 0 | – | +40 % | 0.51 | 31 % | 42 % | – | +0.16 | +0.19 | -0.30 | +0.10 |
| Bromley FC, C fit (film T) | 1 | C = 1.85 (LOCO 1.76-1.94) | -4 % | 0.38 | 36 % | 55 % | 0.41 | +0.16 | +0.19 | -0.30 | +0.10 |
| Bromley, C and β fit (h_fg' = h_fg + β cp ΔT) | 2 | β → 2900 (sensible-heat limit) | -4 % | 0.37 | 39 % | 57 % | 0.39 | | | | |
| Darr 2016 low-g, as published (sat, x_eq) | 0 (7 pre-fitted) | – | -15 % | 0.32 | 52 % | 73 % | – | +0.15 | -0.18 | +0.11 | +0.07 |
| Darr 2016 × C | 1 | 1.19 | +1 % | 0.27 | 62 % | 78 % | 0.28 | | | | |
| Dittus-Boelter vapour, Re_v = GD/μ_v (sat) | 0 | – | +271 % | 1.38 | 1 % | 1 % | – | +0.47 | | | |
| Dittus-Boelter vapour × C | 1 | 0.27 | -1 % | 0.44 | 33 % | 49 % | 0.46 | +0.47 | -0.26 | | |
| Miropolskii (x_eq) / (x = 0) | 0 | – | -88 % / -97 % | 2.2 / 3.6 | 0 | 0 | – | | | | |
| Groeneveld 1973 tube (x_eq) | 0 | – | -78 % | 1.61 | 0 | 1 % | – | | | | |
| Couette laminar film, ℓ = z | 0 | – | -94 % | 2.9 | 0 | 0 | – | | | | -0.39 |
| Couette, ℓ = D (≡ Bromley with C = 0.5) | 0 | – | -74 % | 1.40 | 0 | 0 | – | | | | |
| Couette, ℓ = λ_KH = 2πσ/(ρ_v u_l²) | 0 / 1 | – | -95 % | 3.3 / 1.5 | 0 | 0 | 1.67 | +1.18 | +0.87 | | |
| Penetration (sensible-heat) model, C fit | 1 | C = 2.35 (LOCO 2.22-2.44) | -4 % | 0.35 | 40 % | 58 % | 0.37 | +0.17 | +0.21 | +0.05 | +0.10 |
| **Pe-law: h = C (k_v/D) Pe_l^n** | **2** | **C = 7.78, n = 0.345** (LOCO 6.6-10.5, 0.31-0.36) | 0 % | **0.27** | 61 % | **78 %** | **0.30** | +0.01 | +0.07 | +0.18 | +0.10 |
| Pe-law × Darr entrance factor [1 + (D/(z+0.1D))^0.7] | 2 | C = 6.95, n = 0.345 | 0 % | 0.26 | 62 % | 79 % | 0.29 | +0.01 | +0.07 | +0.18 | +0.03 |
| Re_v power law, C Re_v^n Pr^0.4 k_v/D (sat) | 2 | C = 0.67, n = 0.43 | +3 % | 0.28 | 61 % | 77 % | 0.31 | | | | |

Definitions used in the table:
- Pe_l = u_l D/α_v,f, with u_l = G/ρ_l. Vapour properties are at the film temperature (Ti+Tsat)/2.
- Penetration model: h = C·(2/√π) k_v √(u_l/(α_v D)). It assumes heat goes into superheating vapour elements that touch the wall over a length D while convected at u_l. It is gravity-free and has no explicit ΔT dependence.

### Interpretation

- **The laminar Couette film model.** For a vapour film growing from the inlet between the wall and a liquid core moving at u_l, the derivation gives δ² = 4 k_v ΔT z/(ρ_v u_l h_fg'). This model under-predicts h by a factor of about 18, and it predicts a z^-1/2 drop from upstream to downstream (-39 %) where the data show -11 %.
  - The measured film thickness, k_v/h, is 50-150 µm. It does not grow with z, so the film is limited locally by mixing or interfacial dynamics, not by vapour accumulation.
  - A Kelvin-Helmholtz length cutoff gives G^1.5 dependence, which the data reject.
- **Bromley forced-convection film boiling.** Its √(u_l) = G^0.5 form is about right. C = 2.7 over-predicts by 40 %; with vapour properties at the film temperature the best C is 1.85.
  - Its residual still decreases with ΔT (-0.30): the data have almost no ΔT dependence.
  - When β in h_fg' is left free, the fit pushes β → ∞. The latent term then becomes negligible and Bromley collapses onto the sensible-heat penetration form.
  - So in low g the wall heat seems to go mainly into superheating vapour that is constantly renewed by the moving core, not into evaporation across a stagnant film.
- **Best scaling with the fewest constants: h ∝ G^≈0.35 ΔT^0 P^≈0.35.** The Pe-law (2 constants) captures it with no residual trend in G or P. Its RMS(ln) of 0.27 sits against a scatter floor of 0.17.
  - In dimensional form: h = 7.78 k_v^0.655 (ρ_v cp_v u_l)^0.345 / D^0.655, with film-temperature properties.
  - The 1-constant Darr × 1.19 has a similar CV error (0.28). However, its structure carries 7 constants fitted to these same data plus the 1-g data.
- **Remaining trends.**
  - Every closure without a position term over-predicts the downstream station by about 10 %. Darr's entrance factor removes most of that.
  - In case 10, the start-up burst at 0.5-3 s is under-predicted even with the lag-corrected G (fig09).
  - Within-run time trends are small (|slope| < 0.03 per 10 s) for all forms except Bromley and the vapour-convection forms.
- **Uncorrected G or published h:**
  - With uncorrected G, all RMS values rise by 0.01-0.05.
  - With the published h, the Darr correlation is unbiased (+1 %, RMS 0.30), as expected since it was fitted to that h. The Pe-law gives 0.29 and Bromley needs C = 1.59.

## 4. Rewet temperature (fig10, `out/t04_rewet.md`)

**Theory inputs:**
- Homogeneous nucleation: Lienhard (1976), T_hn = Tc(0.905 + 0.095 Tr^8), giving 114.4-116.0 K.
- Liquid spinodals:
  - van der Waals at the local P: 107-110 K.
  - The CoolProp (Span 2000) EOS evaluated with liquid phase imposed: about 120 K. This is indicative only, since the EOS is not constrained in the metastable region and is non-monotonic below 116 K.
- Contact temperature: e = √(kρc). For 304SS, k and cp come from the NIST cryogenic material-properties database fits (log10 polynomials, `common.ss_k/ss_cp`). For LN2, saturated-liquid properties. e_l/e_w = 0.09-0.10, so the contact criterion adds only +1.5-3 K to T_hn.
- Darr-Hartwig: T_wet = 0.844 Tc (1 + 0.060 We_D^0.208), with We_D = G²D/(ρ_l σ) (Darr & Hartwig, LH2 paper, eq. 18). The full Darr 2019 form uses the Baumeister factor B = exp(3.06e6 β_w) erfc(1751.5 β_w^0.5), which is 0.69-0.71 for 304SS, with the surface-energy factor set to 1.

| model | mean error K | MAE K | corr. with measured |
|---|---|---|---|
| T_hn (Lienhard) | -18.8 | 18.8 | 0.23 |
| contact T = T_hn | -16.5 | 16.5 | 0.22 |
| contact T = T_sp(EOS) | -11.4 | 11.4 | -0.05 |
| contact T = T_sp(vdW) | -23.9 | 23.9 | 0.40 |
| 27/32 Tc | -27.4 | 27.4 | – |
| Henry (1974) correction on the contact-T_hn base | +0.3 | 6.0 | -0.07 |
| Darr-Hartwig 0.844Tc(1+0.06We^0.208) | -9.9 | 9.9 | -0.02 |
| Darr 2019 full (B-factor) | -16.6 | 16.6 | -0.17 |
| De Salve-Panella | -11.6 | 11.6 | -0.55 |
| Berenson 1 g / Henry-Berenson 1 g | +32.9 / +53.0 | | (g-dependent; ∝ g^-1/6, diverges as g → 0) |

**What the measurements show:**
- T_MFB = **133.9 ± 6.3 K** (11 stations; 133.3 ± 6.3 K over 30 TCs). That is 1.06 ± 0.05 Tc.
- The contact temperature at rewet is 130.5 ± 6.0 K. That is about 15 K above T_hn (115 K) and about 10 K above the indicative EOS spinodal (120 K). The thermodynamic-limit theory is therefore biased low by 11-17 K.
- **Pressure.** There is no pressure trend in T_MFB from 93 to 762 kPa: the regression gives T_MFB = 122.6 + 0.17·Tsat - 8.8·dn. So ΔT_MFB falls with P, from 60 K at 93 kPa to 37-43 K at 570-760 kPa.
  - Case 10, at the highest P, is not special: 142.5 K upstream and 133.3 K downstream.
  - Theories based on T_hn predict a nearly constant T_MFB (117 K), which is consistent with the absence of a trend but too low by ~17 K.
  - Henry's correction predicts T_MFB decreasing with P (140 → 130 K), which the data do not show.
- **Mass flux.** No G trend is resolved: the ΔT_MFB slope with ln G is negative and not significant. The +12 K rise from G = 100 to 400 kg/m²s predicted by the We term of Darr-Hartwig is not seen.
- **Station effect.** Systematically, the upstream station rewets 8.9 K hotter than downstream (+3 to +14 K in all 5 paired cases).
- **The best one-constant rules** (leave-one-case-out RMSE):
  - T_MFB = 134 K: 6.5 K
  - T_MFB = T_hn + 18.8 K: 6.4 K
  - T_contact = 130.5 K: 6.6 K
  - ΔT_MFB = 44 K: 8.9 K
- **Quench-front velocity** between the stations: 44-89 mm/s, i.e. 0.11-0.35 times u_l.

## 5. CHF and transition (fig11, `out/t05_chf_transition.md`)

- **Measured quench CHF** (per-TC lumped, station median): 75-123 kW/m². It depends weakly on G (∝ G^0.17) and on P (∝ P^0.21). ΔT_CHF is 19-37 K (mean 24 K), far larger than the 2 K assumed in the Darr-Hartwig LH2 set.
- **Katto-Ohno (1984)**, saturated inlet, L = z: over-predicts by a factor of 1.66 on average (range 1.2-2.5). It is the steady heated-tube CHF and is gravity-free, but the transient quench peak is lower.
- **Zuber pool CHF:** 230 kW/m² at 1 g and 78 kW/m² at 0.01 g (0.78 of measured), falling as g^1/4 toward 0. Agreement at residual g = 0.01 is a coincidence. The flow-driven CHF does not vanish in low g; Zuber is not a usable low-g closure.
- **Transition curve shape.**
  - Between (ΔT_MFB, q_MFB ≈ 6-12 kW/m²) and (ΔT_CHF, q_CHF), the endpoint log-log slope is very steep: n = -4.3 (IQR -4.9 to -3.6).
  - The curve is not log-log linear: interpolating that way gives RMS(ln q) of 1.18. Linear-in-T interpolation (Darr-Hartwig) gives 0.60.
  - The flux jumps right after MFB. The best form is q = q_MFB + (q_CHF - q_MFB) θ^0.5, with θ = (T_MFB - T)/(T_MFB - T_CHF), giving RMS 0.19.

## 6. 0-D lumped-wall surrogate (fig12, `out/t06_surrogate.md`)

**Model.** At each station:

    m''_w cp(T) dT/dt = -q''(T) + 0.33 kW/m²

with Tsat(t) and P(t) measured at the station, G(t) lag-corrected, and x from the energy balance. The boiling curve is: film boiling (closure) → transition θ^0.5 → nucleate q = q_CHF (ΔT/24 K)^2.5. q_CHF = 100 kW/m². The 20 stations are integrated with Heun steps of 5 ms.

**Oracle check.** Driving the model with the measured q''_ec gives a median RMS of 6 K, so the energy bookkeeping is consistent.

Mean of the per-station RMS(T_sim - Ti) with T_MFB = 134 K. The t130 and t100* columns are errors in the time to reach 130 K and max(100 K, Tsat + 8 K), averaged over stations where both measurement and simulation cross. "mismatch" is the number of stations that quench in only one of the two.

| film-boiling closure | RMS K (mean / median) | t130 MAE (bias) s | t100* MAE (bias) s | mismatch 130/100 |
|---|---|---|---|---|
| Pe-law (2 fit) | **12.4 / 10.0** | **1.7 (-0.5)** | **2.1 (-1.3)** | 1/2 |
| Penetration, C = 2.34 | 15.3 / 14.6 | 3.6 (-3.6) | 4.8 (-4.8) | 1/2 |
| Bromley, C = 1.85 | 16.7 / 16.6 | 4.6 (-4.6) | 5.6 (-5.6) | 1/2 |
| Bromley, C = 2.7 | 33.9 / 37.3 | 8.4 (-8.4) | 9.2 (-9.2) | 3/6 |
| Darr 2016 (published constants) | 13.4 / 11.8 | 3.3 (+3.3) | 2.0 (+1.8) | 5/2 |
| Dittus-Boelter vapour (sat) | 77.8 / 73.1 | 14.4 (-14.4) | 13.5 | 7/10 |
| Miropolskii (x_eq) | 71.6 / 65.3 | never quenches | – | 11/8 |

How each closure fails:
- **Darr 2016:** fails to quench at 5 stations because its h is about 15 % low against the energy-consistent flux.
- **Bromley and penetration:** cool too fast and quench 4-8 s early.
- **Pe-law:** its largest errors are
  - case 10 upstream, quenching 3 s late because the start-up burst is missed;
  - the downstream stations of cases 6 and 7, which quench 3-5 s early because the 0-D model has no quench-front propagation and only a 0-10 % downstream h penalty;
  - cases 1 and 5, which run 25-40 K too cold by the end of the run.

**Choice of T_MFB rule.** With the Pe-law, the choice changes the RMS by less than 1.2 K and the t130 MAE by less than 0.3 s:

| T_MFB rule | RMS K | t130 MAE s | t100* MAE s | mismatch 100 |
|---|---|---|---|---|
| measured | 11.8 | 1.4 | 1.6 | |
| 134 K | 12.4 | 1.7 | 2.1 | |
| contact-T_hn (117 K) | 12.8 | 1.7 | 3.0 | 3 |
| Darr-Hartwig | 11.7 | 1.7 | 1.6 | |

In low g the film-boiling closure dominates the wall-temperature history. Within the transition, T_MFB mostly shifts how fast the final drop happens.

With uncorrected G the Pe-law gives nearly identical results: 12.8 K RMS and 1.6 s t130 MAE.

**Caveats.**
- The constants of the fitted closures are in-sample for the surrogate. Their LOCO error in h is only 3-10 % larger than in-sample.
- x for the Miropolskii and Darr closures comes from the measured q'', so those runs are semi-a-priori.

## Recommendations for the 1-D OpenFLUME runs

1. **Film boiling:** h = 7.78 (k_v/D) Pe_l^0.345, Pe_l = (G/ρ_l) D ρ_v cp_v / k_v, with vapour properties at the film temperature. Optionally multiply by [1 + (D/(z+0.1D))^0.7], which brings C to 6.95. Alternatives are Bromley with C = 1.85, or Darr 2016 × 1.19.
2. **Mass flux:** use the lag-corrected G (2.5 s), and remember that the published q'' and h of cases 1, 2, 3, 5 and 6 are 26 % low.
3. **Rewet:** T_MFB = 134 ± 6 K (pressure-independent in these data), or T_hn(Lienhard) + 19 K. Upstream is about 9 K higher than downstream.
4. **Transition and nucleate boiling:** θ^0.5 transition up to q_CHF ≈ 100 kW/m² at ΔT_CHF ≈ 24 K; nucleate q ∝ ΔT^2.5.

> **These Phase-1 recommendations are superseded by Phase 2 below and by `CLOSURE_SPEC.md`.**

---------------------------------------------------------------------------------------------------

# Phase 2: consistent dataset and final closure

New scripts: `a07_consistent_q.py`, `a08_film_laws.py`, `a09_rewet_theory.py`, `a10_boiling_curve_sensitivity.py`. Tables are `out/t07..t10*.md` and figures are `fig13`–`fig17`. The final closure is in **`CLOSURE_SPEC.md`**.

## P2.1 One consistent heat-flux dataset (fig13, `out/t07_consistent_q.md`)

**Definition.** For all 10 cases:

    q''_E = −m''_w cp(T_m) dT_m/dt + q''_par(T_m)

- T_m is the mean of the three station TCs.
- m''_w = 4.20 kg/m², and cp is NIST 304SS.
- dT/dt comes from a Savitzky–Golay filter, 9 points (0.6 s), 2nd order. Using 5 or 15 points changes the median h by < 0.1 % (RMS 0.05 in ln h).
- q''_par = 98.6 + 2.2 (293 − T) W/m². It is fitted to (q_inverse − q_lump) on the self-consistent cases 4 and 7–10, i.e. it is the authors' parasitic model. That is about 2× a radiation-only estimate with ε = 0.453, which is plausible once gas conduction and the tube ends are included.
- Radial correction to the inner surface: T_i = T_m − q'' · 0.254 mm / k. This is 0.2 K in film boiling and about 2.5 K at the CHF peak.
- The radial diffusion time is 0.07 s, about one sample, so no time-shift correction is applied.

**Check against the published flux.** q_pub/q_E is 0.995–0.997 for cases 4 and 7–10, and 0.728–0.737 for cases 1, 2, 3, 5 and 6.

**The 0.75 factor is a data-reduction scaling, not a different tube.**
1. Splitting the published flux as q_pub = α q_store + β q_par gives α = 0.73–0.81 in the anomalous cases and 1.00 elsewhere, with β unchanged. Only the wall-storage term (the m''cp the inverse model used) is scaled.
2. A storage factor of 0.752 would correspond to a 0.388 mm wall at the same OD. But the published T_i − T_m implies a radial conduction length (T_m − T_i)k/q_pub of 241–247 µm in **all** cases. That matches the 0.51 mm geometry (251 µm); a 0.388 mm wall would give 191 µm. So the inverse model used the 0.51 mm wall for conduction and a 0.75× heat capacity for storage, which is internally inconsistent. The paper also describes a single test section.
3. Hold-out test: the √(k_v cp_v G/D) law fitted only on cases 4 and 7–10 predicts the anomalous cases with a geometric-mean measured/predicted ratio of 1.06 using h_E (per case 1.02, 1.33, 1.12, 0.73, 1.20). With the published h the ratio is 0.77.

## P2.2 Film-law discrimination (fig14, `out/t08_film_laws.md`)

Setup:
- Dataset: 4422 samples, 20 case-stations, equal weights; τ = 2.5 s lag.
- LOCO is by case.
- Trend columns are multivariate slopes of ln(pred/meas).

| law (one constant) | C (LOCO range) | RMS ln | LOCO | σ case-station | trend lnG | lnP | lnΔT | RMS at τ=0 |
|---|---|---|---|---|---|---|---|---|
| F1 √(k ρ cp u_l/D), sat | 2.38 (2.23–2.52) | 0.373 | 0.402 | 0.32 | +0.07 | **+0.35** | +0.06 | 0.388 |
| F1, film | 2.68 (2.53–2.79) | 0.336 | 0.360 | 0.28 | +0.06 | +0.27 | +0.02 | 0.361 |
| F2 √(k cp G/D), sat | 0.277 (0.269–0.288) | 0.259 | 0.270 | 0.18 | +0.07 | −0.15 | +0.06 | 0.306 |
| F2, film | 0.229 | 0.283 | 0.299 | 0.22 | +0.05 | −0.24 | +0.23 | 0.321 |
| F3 (k/D) Pe_l^(1/3), sat | 9.39 (9.02–9.58) | 0.239 | 0.250 | 0.17 | −0.10 | +0.20 | +0.06 | 0.271 |
| **F3, film** | **8.55 (8.27–8.80)** | **0.235** | **0.246** | **0.17** | −0.11 | +0.11 | +0.16 | **0.267** |
| F5 Darr 2016 × C, corrected exp[−(L/D)^c5] | 1.19 | 0.243 | 0.253 | 0.17 | +0.10 | −0.15 | +0.05 | 0.306 |
| Bromley √(k ρ h_fg' u_l/(D ΔT)), film | 1.91 | 0.377 | 0.402 | 0.30 | +0.05 | +0.23 | **−0.34** | 0.418 |

Forms with more fitted parameters:

| form | fitted values | RMS ln | LOCO |
|---|---|---|---|
| F3 with free exponent, film | n = 0.343 (LOCO 0.31–0.36), C = 8.0 | – | 0.263 |
| F3 with free exponent, sat | n = 0.29 | – | 0.261 |
| F4 Re_v^n Pr^m, sat (3 constants) | n = 0.33, m = 2.9 | 0.227 | 0.280 |
| F4 Re_v^n Pr^m, film (3 constants) | – | – | 0.303 |

F4 over-fits: its LOCO is worse than the one-constant F3. Darr's droplet term is negligible here; the "corrected" and original exponentials give the same RMS, 0.243.

**F1 Higbie renewal length.** The implied renewal length is λ/D = 4/(πC²) = 0.18–0.22, i.e. about 2.3 mm. F1 is rejected because its implied pressure exponent at fixed G (0.55–0.63) lies outside the data's interval.

**Pressure identifiability.** Profile of a pressure exponent p with the G exponent fixed at 0.5:
- h/G^0.5 ∝ P^p gives p = 0.22, with 95 % interval [0.04, 0.40] from the profile and [0.06, 0.39] from a case bootstrap. With τ = 0, p = 0.26 [0.05, 0.47].
- Leverage pairs at similar G: case 6 against cases 3, 4 and 5 gives p = 0.17, 0.41 and 0.68 (pooled 0.40). Case 4 is uncorrected and cases 3, 5 and 6 are corrected, and the 6-vs-4 pair still gives 0.41, so the 0.735 correction does not create the P effect.
- Implied exponents at fixed G for each form: F1 0.55–0.63 (rejected), F2-sat 0.13, F3-film 0.38 combined with G^(1/3).
- The free fit (G^0.35 P^0.37) is matched almost exactly by F3-film. G and P are only jointly identified, so F2 (G^0.5 P^0.13) remains admissible but scores worse (LOCO 0.270 vs 0.246, and less robust to τ).

**Choice: F3 with film properties, C = 8.55 ± 0.27 (LOCO range).** Physical reading: a Lévêque thermal layer in the linear shear of the vapour gap. The vapour velocity is 0 at the wall and about u_l at the interface of the inverted-annular liquid core, and the layer is renewed by interfacial waves and slugs.
- Lévêque gives h = 0.807 (k/λ)(γλ²/α)^(1/3). With γ = u_l/δ this becomes C = 0.807 (D²/(δλ))^(1/3).
- C = 8.55 then implies δ ≈ λ ≈ 0.34 mm, which is consistent with the measured equivalent film thickness k_v/h of 0.05–0.3 mm.
- The 1/3 exponent is supported by the free-exponent film fit (0.343, LOCO 0.31–0.36).

**Entrance and station factor.** Upstream h is measured higher than downstream in 9 of 10 cases, by +12 % on average (ln 0.115; largest in cases 1, 2 and 10).
- Darr's entrance factor gives C = 7.60, RMS 0.229, LOCO 0.240.
- That is a marginal gain (−0.006) against a case-to-case σ of 0.17. It is optional and not required.

**Quality.** No x dependence is needed for x_e ≤ 0.05: on that subset C = 8.50 and RMS 0.21.
- The apparent x trends are not identifiable, because x_e is built from the measured q (within-run correlation with ln q is 0.22). Their sign also flips between τ = 2.5 and τ = 0.
- Darr's (1−x)^−0.76 factor is < 4 % in this range.

## P2.3 Rewet: contact temperature equals T_c (fig15, `out/t09_rewet_theory.md`)

**Rule:** T_MFB = T_c + β(T_c − T_sat), with β = e_l/e_w(T_MFB) = 0.08–0.09, giving 128.2–130.8 K.

| rule | constants | mean error K | MAE K | upstream mean / MAE | downstream mean / MAE |
|---|---|---|---|---|---|
| **T_c contact** | 0 | −4.6 | 5.6 | −8.4 / 8.4 | **+0.0 / 2.2** |
| Henry 1974 (on T_hn contact base) | 0 | +0.3 | 6.0 | −3.3 / 6.4 | +4.6 / 5.4 |
| T_hn + 19 K | 1 | 0.0 | 5.1 | −3.9 / 5.5 | +4.7 / 4.7 |
| 134 K | 1 | +0.1 | 5.1 | −3.8 / 5.4 | +4.8 / 4.8 |

**External checks** (same onset definition, digitized data):

| data | measured T_MFB | rule |
|---|---|---|
| Darr 1 g tests 31/40, both stations | 128.5–129.9 K | 129.3–129.4 K |
| Hartwig 2022 F2 exit (0.05 g, 556 kPa) | 135.7 K | 128.6 K |
| NBS 9264 copper LN2, 6 knees (hand-clicked, low quality) | 121–141 K, mean 130.3 | 126.7 K |

For tests 31/40, Kassemi's quoted 135 K uses an earlier knee definition.

**The upstream +9 K offset.** It is robust to the onset definition: 6.2–13 K across 9 definitions (median/first/min/max TC, 1.5×/3× thresholds, curvature, station-mean TC, inverse T_i). The candidate explanations rule out as follows:
- **Pressure or T_sat:** the rule predicts only −0.02 to −0.4 K upstream minus downstream.
- **Mass flux:** case 6 has ΔG = −5 kg/m²s but a +13.7 K offset.
- **Start-up burst:** rewet happens at 8–23 s, well after the burst.
- **Axial conduction:** the quench front moves at 44–89 mm/s, so the dry-side precursory length k/(ρcU) is about 0.1 mm.
- **TC skew:** the ~10° skew applies to both stations.
- **Quality:** the downstream station has higher x, but the Iloeje-type effect is < 1 K.

The remaining candidate is an **entrance effect** at z/D = 12.8: coherent liquid-core impingement from the inlet, consistent with the +12 % film h at the same station. The optional fix is a fitted +8 K for z/D ≲ 15.

**Predictions of the rule:**

| system | T_MFB | β |
|---|---|---|
| LN2 / Cu (OFHC, NIST fits) | 126.6–126.9 K | 0.01 |
| LH2 / Cu | 33.2–33.3 K | 0.01–0.02 |
| LH2 / 304SS | 34.5–35.9 K | 0.22–0.28 |
| LN2 / 304SS, 1 atm / 300 / 600 kPa | 130.7 / 129.4 / 128.5 K | 0.08–0.09 |
| LN2 / FEP surface (Teflon 959G-203: k = 0.195, ρ = 2150, cp ≈ 700 assumed) | 170 / 158 / 149 K at 1 atm / 300 / 600 kPa | 0.76–0.90 |

For the FEP coating at 550 kPa, the surface rewets at 150 K, which puts the metal at 150.7–152.7 K for 1–4 layers, against 128.7 K on bare SS. That means about 15 K earlier rewet in metal temperature. The semi-infinite coating assumption holds for contacts shorter than 1.8 ms (1 layer) to 32 ms (4 layers).

**1 g.** Darr's 1-g tests at 315–326 kPa are predicted to within 1 K under the same onset definition.

## P2.4 Boiling curve and sensitivity (fig16 tornado, fig17, `out/t10_sensitivity.md`)

Boiling curve used: Cooper nucleate boiling plus a Dittus–Boelter liquid floor; Katto–Ohno CHF (50–195 kW/m² over the records); ΔT_CHF from the intersection of the NB curve with q_CHF; θ^0.5 transition; T_c-contact rewet; F3 film boiling.

The 0-D surrogate against station-mean TCs gives a mean RMS of **13.1 K (median 9.3 K)**, t130 MAE **1.8 s** (bias +0.5 s), and t100* MAE **1.6 s** (bias −0.3 s). Mismatches (stations that quench in only one of simulation and data): 1 at 130 K, 2 at 100 K.

| change | ΔRMS K | Δt130 bias s | Δt100* bias s |
|---|---|---|---|
| film h ×0.7 / ×1.3 | +6.7 / +15.0 | no quench / −4.5 | no quench / −4.4 |
| T_MFB −5 / +5 K | −0.1 / +0.8 | 0.0 / −1.1 | +0.6 / −1.3 |
| CHF ×0.5 / ×1.5 | −0.1 / +0.1 | 0 | +0.4 / −0.2 |
| transition n 0.25 / 1 | +0.1 / −0.1 | 0 | −0.1 / +0.5 |
| NB ×0.5 / ×2 | ±0.04 | 0 | ±0.01 |
| CHF = 100 kW/m² constant | −0.02 | 0 | +0.2 |

**So the film-boiling C matters most.** Its ±3 % LOCO uncertainty is equivalent to ≈ ±0.4 s on t130. T_MFB is second, and CHF, transition and NB can be crude.


---------------------------------------------------------------------------------------------------

# Phase 3: position dependence, the quench front, and late-time over-cooling

Scripts: `a11_position_laws.py`, `a12_front_and_late.py`, `a13_surrogate_position.py`. Tables are `out/t11`–`t13*.md`. Figures: `fig18_front_reset.png` (downstream h/law across the upstream rewet) and `fig19_late_time.png` (cases 1, 4, 5). The final recommendation is in **CLOSURE_SPEC.md §10**.

## P3.1 Developing-layer laws (`t11`)

All laws use the form h = C (k_vf/D) Pe_l^(1/3) (D/L_eff)^(1/3), with film-temperature properties, the lag-corrected G, the consistent q, and LOCO by case. The "up − down" column is the mean upstream-minus-downstream station residual (0 means no station bias).

| L_eff | fitted constants | LOCO | σ case-station | up − down |
|---|---|---|---|---|
| D (Phase 2) | C = 8.55 | 0.246 | 0.167 | −0.115 |
| z | C | 0.262 | 0.191 | +0.215 (over-corrects) |
| z + z0 | C, z0 = 0.46 m (0.38–0.59) | 0.240 | 0.156 | 0 |
| zΛ/(z+Λ) | C = 16.4, Λ = 0.129 m (0.10–0.16) | 0.240 | 0.156 | 0 |
| zΛ/(z+Λ) with u_l = G(1−x)/ρ_l | C = 15.8, Λ = 0.101 m (0.08–0.12) | 0.240 | **0.151** | 0 |
| D × Darr entrance factor | C = 7.60 | 0.240 | 0.157 | −0.04 |
| D-free penetration √(kρcp u_l/(z+z0)) | 2 | 0.357 | 0.273 | 0 |

**Notes:**
- The measured h_up/h_dn ratio at matched times is 1.01–1.31 (mean 1.12). It is largest in cases 1 and 2, which also have the largest downstream quality.
- A free (1−x)^m exponent is not identifiable, because x_e is built from the measured q.
- A dimensionally consistent local Lévêque law needs a shear length, taken here as D. The fully D-free form fails.

## P3.2 Front-reset test and front mechanism (`t12`, fig18)

**Reset test.** I measured the change in downstream ln(h/law) across the first upstream rewet, comparing the windows [−2, −0.2] s and [+0.3, +2.3] s.
- Per case: −0.13, −0.16, +0.29, +0.07, +0.10, +0.03 for cases 3, 6, 7, 8, 9, 10. The mean is +0.03 ± 0.045.
- The per-TC values agree in sign.
- The placebo spread at other times is σ = 0.11.
- A front-restarted pure-z law predicts +0.155 and is rejected at about 2.8σ. The saturating law predicts +0.045, which cannot be resolved.

**Front budget (mm/s).** The measured speed is 44–89 mm/s at the 2×-onset points; at the 130 K crossings of station-mean TCs it is 47–106 mm/s.

| mechanism | speed or delay | sets the front? |
|---|---|---|
| conduction-controlled (Yamanouchi; h_wet 5 kW/m²K, dry wall T_MFB + 20 K) | 4–6 mm/s | no, too slow |
| liquid supply ṁ(1−x)h_fg divided by the wall energy to quench | 680–3300 mm/s | no, not limiting |
| kinematic/dispersed: ΔT_wall/rate + ΔT_MFB/rate | 2.8–3.8 s (cases 7–10) against 2.8–3.9 s measured; case 6: 8.0 against 5.7 s | **yes** |

In the kinematic estimate, the axial h gradient accounts for 0.5–3.3 s of the delay and the T_MFB offset for 0.8–4.8 s.

**Conclusion.** The front is a dispersed, kinematic consequence of an axial gradient in film h (and T_MFB). With a z-independent closure, every node reaches T_MFB together, so the model front runs at about 1–2 m/s. That is exactly the OpenFLUME symptom.

## P3.3 Late-time over-cooling, cases 1 and 5 (`t12` (c), fig19)

Measured h relative to the law, early half → late half of the film-boiling period:

| case | upstream | downstream |
|---|---|---|
| 1 | 0.86 → 0.90 | 0.67 → 0.68 |
| 5 | 0.84 → 0.90 | 0.71 → 0.89 |

- The law over-predicts from the start; it does not drift late. The cold end is error accumulated over 24 s.
- Case 5 has steady G (CV 0.07), low x and no pressure collapse. It is about 25 % below case 4 at the same G and P, which is run-to-run scatter.
- Case 1 (G 24–37, oscillating with CV 0.3; downstream x up to 0.19) is partly explained by liquid depletion: the (1−x)^(1/3) term reduces the downstream error by about 5 K.
- Only cases 7, 9 and 10 show the measured/law ratio falling late (−0.1 to −0.5). That is linked to the start-up burst and the meter lag of G, not to quality.
- There is no evidence that the G^(1/3) law breaks down at G ≈ 30: case 1 upstream sits at 0.86–0.90 of the law, and case 2 at 1.16–1.20.

## P3.4 0-D surrogate with position laws (`t13`)

| film law (+ T_c rewet) | RMS mean K | t130 MAE s | t130 bias up / dn s | t100* MAE s | delay dn−up at 130 K, s (measured 5.4, 2.6, 2.6, 2.4, 3.8) | U mm/s (measured 96) |
|---|---|---|---|---|---|---|
| Phase 2, L = D | 13.1 | 1.78 | +2.1 / −1.1 | 1.62 | 0.1–0.2 | 1360 |
| Darr entrance factor | 12.3 | 1.29 | +1.4 / −0.2 | 1.74 | 1.4–2.3 | 153 |
| L3 (Λ = 0.129 m) | 12.0 | 1.16 | +1.0 / +0.3 | 1.78 | 2.1–3.5 | 104 |
| **L3x (Λ = 0.101 m, (1−x))** | **11.8** | **1.10** | **+0.9 / −0.0** | 1.87 | 1.9–3.5 | **114** |
| L3x + entrance rewet +8 K | 12.8 | 1.32 | −0.6 / −0.0 | 2.50 | 3.1–6.1 | 69 (over-corrected) |

**Recommendation:** L3x, with no dynamic front tracking and no rewet entrance increment (CLOSURE_SPEC §10). Cases 1 and 5 remain about 25 K too cold at the end of the run; that is within the ±17 % run-to-run h scatter.

---------------------------------------------------------------------------------------------------

# Phase 4: referee-driven checks (`paper/review_physics.md`)

Scripts: `a14_loto.py` → `out/t14_loto.md` and `out/loto_constants.json`; `a15_renewal_rewet_gravity.py` → `out/t15_renewal_rewet_gravity.md`.

## P4.1 Nested leave-one-test-out (M8, M9)

Within each fold, τ (grid 0–5 s), C and Λ are all re-fitted on the 9 training tests. The training error is minimised over (τ, Λ) on a grid, with C in closed form. x_e is recomputed for every τ.

**Final law (L3x), consistent flux:**
- All data: τ = 2.50 s, C = 15.69, Λ = 0.098 m. (The Phase-3 values 15.8 / 0.101 differ only through grid resolution.)
- Fold spread: τ 2.00–2.75 s, C 14.9–16.4, Λ 0.081–0.118 m.
- Nested LOTO error is 0.247. With τ fixed at 2.5 s it is 0.240; with τ fixed at 0 it is 0.274 (C = 16.6, Λ = 0.118 m).
- Per-fold constants are in `out/loto_constants.json`.

| law (nested τ) | consistent flux q_E: C, LOTO | published flux q_pub: C, LOTO |
|---|---|---|
| L3x (final) | 15.69, Λ = 0.098 m; 0.247 | 13.56, Λ = 0.108 m, τ = 2.25 s; 0.323 |
| L0, (k/D)Pe^(1/3) | 8.55; 0.249 | 7.26; 0.326 |
| Darr 2016 × C | 1.19; 0.261 | 1.02; 0.310 |
| √(k_v cp_v G/D), saturated properties | 0.277; 0.279 | 0.235; 0.332 |

**Paired test-bootstrap of LOTO differences** (4000 resamples over tests, q_E). Each row is the second law's LOTO minus the first's:

| comparison | difference | SE | 95 % CI |
|---|---|---|---|
| L3x vs L0 | +0.003 | 0.005 | [−0.007, 0.014] |
| L3x vs Darr | +0.015 | 0.030 | [−0.045, 0.069] |
| L3x vs √-law | +0.032 | 0.031 | [−0.030, 0.088] |
| Darr vs √-law | +0.018 | 0.004 | [0.009, 0.025] |

The final law, L0 and Darr × C are **statistically indistinguishable**. Only the √-law is resolvably worse, and only relative to Darr. With q_pub, Darr × C is nominally best, and again none of the differences are significant except √ vs Darr.

## P4.2 The 0.735 correction and the exponents (M3)

**Free-law G exponent depends on the correction.**
- With q_E: b_G = 0.35 [0.12, 0.46] and b_P = 0.36 [0.16, 0.60] at τ = 2.5 s. At τ = 0: b_G = 0.31 [−0.12, 0.42].
- With q_pub: b_G = 0.45 [0.17, 0.55] and b_P = 0.42 [0.23, 0.73].
- So the correction moves the G exponent from about 1/2 toward about 1/3, as the referee suspected. The confidence intervals overlap.

**Matched-G/P comparison** (film boiling, medians; test 4 uncorrected):

| test | G | P kPa | h_E/h_E(4) | h_pub/h_pub(4) | expectation with no exponent assumed |
|---|---|---|---|---|---|
| 3 | 123 | 92 | 1.10 | 0.81 | < 1 (lower G and P) → favours q_pub |
| 5 | 177 | 122 | 0.86 | 0.63 | ≳ 1 (higher G, similar P) → favours q_E |
| 6 | 145 | 289 | 1.41 | 1.05 | > 1 (higher P) → both plausible, q_E closer to the P^0.36–0.42 expectation |

This exponent-free test is **mixed**. The strongest arguments for q_E remain internal to the inverse model:
- The scaling is confined to the storage term.
- The radial-conduction length implied by the published T_i is about 245 µm in all 10 tests, which matches a 0.51 mm wall and not 0.388 mm.

**Hold-out ratios** (laws fitted on tests 4 and 7–10, measured/predicted for tests 1, 2, 3, 5, 6):

| law | with q_E | with q_pub |
|---|---|---|
| free law | 0.56, 0.82, 1.44, 0.92, 0.93 (geometric mean 0.89) | 0.40, 0.60, 1.06, 0.66, 0.67 (0.65) |
| √-law | 1.02, 1.33, 1.12, 0.73, 1.20 (1.06) | 0.74, 0.97, 0.83, 0.53, 0.87 (0.77) |

These extrapolate from G ≥ 132 down to G = 22 and carry about ±30 % per-test scatter.

## P4.3 Surface-renewal reading (M1e)

All quantities use u_l = G(1−x)/ρ_l and film properties, over 4422 samples. The exponents are the weighted regression slopes on ln G and ln P.

| quantity implied by measured h | median (p10–p90) | exponent on G | exponent on P |
|---|---|---|---|
| δ·ℓ_r, Lévêque restarted every ℓ_r: 0.807³ k³ u_l/(α h³) | 1.2×10⁻⁷ m² (0.4–2.4×10⁻⁷) | −0.03 | +0.07 |
| Higbie renewal time t_r = 4kρcp/(πh²) | 9 ms (5–28) | −0.69 | +0.33 |
| renewal length u_l t_r | 1.9 mm (0.7–4.1) | +0.38 | +0.35 |

**Empirical result:** δ·ℓ_r is constant to within the scatter over G = 12–700 and P = 90–770 kPa. Examples: a 100 µm layer renewed every 1.2 mm, or 300 µm every 0.4 mm. This constancy is what makes C constant.

**Physical renewal lengths tried.** Each is inserted in Higbie as h = M·2√(kρcp u_l/(π ℓ_r)) with a fitted multiplier M. The data's free law is G^0.35 P^0.37.

| renewal length ℓ_r | exponent on G | exponent on P | RMS | M |
|---|---|---|---|---|
| D | 0.54 | 0.54 | 0.35 | 2.4 |
| λ_KH with Δu = u_l | 1.62 | 1.20 | 1.59 | 7.3 |
| λ_KH with Δu = j_v | 0.74 | 0.64 | 0.76 | 2.6 |
| λ_KH with Δu = j_v − u_l | 0.55 | 0.70 | 0.93 | 4.0 |
| σ/(ρ_l u_l²) | 1.62 | 0.70 | 1.35 | 0.34 |

**None reproduces G^~1/3 P^~0.37, and none fixes the magnitude (M ≠ 1).** C is an empirical renewal intensity. The fitted renewal length (about 2 mm) is of the order of a Kelvin–Helmholtz wavelength with Δu ≈ 1 m/s, but no tested velocity scale gives its G and P dependence.

## P4.4 Rewet criteria (M4)

Evaluated on 11 stations. Kalinin: T_sat + 1.65[0.1 + 1.5√(β_w/β_l) + 0.6 β_w/β_l](T_c − T_sat), with β = 1/(kρc) as transcribed by Darr 2019. Simon–Simoneau: T_c(0.87 + 0.13 P/P_c).

| criterion | mean error K | MAE K | predicted range K |
|---|---|---|---|
| T_c-contact rule | −4.6 | 5.6 | 128.2–130.8 |
| constant 134 K | +0.1 | 5.1 | 134 |
| T_hn + 19 K | 0.0 | 5.1 | 133–135 |
| Henry 1974 | +0.3 | 6.0 | 130–140 |
| Simon–Simoneau 1969 | −22.3 | 22.3 | 110–114 |
| Kalinin et al. 1969 | −29.3 | 29.3 | 99–110 |
| Spiegler 27/32 T_c | −27.4 | 27.4 | 106.5 |

**On these data:**
- The T_c rule, 134 K, T_hn + 19 K and Henry cannot be distinguished. No criterion correlates with the station scatter (|r| ≤ 0.23).
- The measured slope dT_MFB/dT_sat is 0.17, with a test-bootstrap 95 % CI of [−0.85, 0.93]. Henry's slope of −0.46 **cannot be excluded**.
- Kalinin and Simon–Simoneau, both critical-temperature-referenced, are 22–29 K low.

**What would discriminate them:**

| system | T_c rule | 1.061 T_c (scaled constant) | Kalinin | Simon–Simoneau |
|---|---|---|---|---|
| LH₂ on Cu (150–300 kPa) | 33.3 | 35.2 | 24–26 | 29.3–29.8 |
| LH₂ on 304SS | 35.2 | 35.2 | 32.6 | 29.8 |
| LN₂ on Cu | 126.7 | 133.9 | 95.6 | 111.2 |
| LN₂ on an FEP surface | 157.7 | 133.9 | 197.7 | 111.2 |

So copper/LH₂ separates Kalinin and Simon–Simoneau from the T_c rule, and coated tubes separate everything.

**Data that could do this:**
- NBS 9264 LH₂ copper line, Figs 2–7. Our knee detection on the hand-clicked digitization fails (peak-rate ratios < 2 in 10 of 12 traces), so a re-digitization at the knees is needed.
- NBS LN₂ copper (indicative: 121–141 K, mean 130 K).
- Chung 2021 / Hartwig 2022 coated-tube curves, 1L–4L.
- The Darr–Hartwig LH₂ (GRC) database.

## P4.5 Gravity (M6)

The cube sum uses the horizontal Bromley term (Eq. 9, film properties) on the final law:

| test | G | h_law | h_b at 1 g | (h³ + h_b³)^(1/3)/h at 1 g | at 0.1 g | at 0.02 g |
|---|---|---|---|---|---|---|
| 1 | 22 | 49 | 116 | **2.42** | 1.50 | 1.19 |
| 2 | 69 | 78 | 126 | 1.75 | 1.21 | 1.07 |
| 3 | 123 | 62 | 102 | 1.75 | 1.21 | 1.07 |
| 7 | 303 | 144 | 142 | 1.25 | 1.05 | 1.02 |
| 9 | 431 | 179 | 152 | 1.17 | 1.03 | 1.01 |
| 10 | 629 | 254 | 179 | 1.10 | 1.02 | 1.01 |

**Comparison with Darr 2016.** Darr reports that low-g h is up to 25 % below 1 g at Re ≈ 2800 (test 1), i.e. a 1-g/low-g ratio of ≤ 1.33, and no difference at high Re. The cube sum predicts 2.4 for test 1 and 1.17 for test 9. **It over-predicts the gravity effect at low G by about 2×**, so its 1-g success on the Darr tests at G ≈ 312 (predicted ratio about 1.2) is not a strong test.

**Residual g during the fitted window.**
- If Bromley held at 0.1 g (the first 5–10 s), test 1 would be lifted 50 % and tests 2–3 by 21 %. At 0.02 g the lifts are 19 % and 7 %.
- The data show **no** such early elevation. The per-test C for t < 10 s against t ≥ 10 s is 12.7/12.4 (test 1), 16.8/17.2 (test 2) and 20.5/21.4 (test 3).
- Re-fitting C on t ≥ 10 s alone gives 15.97 against 15.69 for all data (+1.8 %).
- So the buoyancy contribution at ≤ 0.1 g is below about 5 % even at G = 22. Bromley over-predicts it there by ≳ 3×, and C is not contaminated by residual gravity.

---------------------------------------------------------------------------------------------------

# Phase 5: pressure-dependent rewet rules and the NBS knees

Scripts: `a16_rewet_family.py` → `out/t16_rewet_family.md`; `a17_nbs_knees.py` → `out/t17_nbs_knees.md`, `out/nbs_knees.csv`, `fig20_nbs_knees.png`. The final rule is in **CLOSURE_SPEC §11**.

**Rule family.** Each rule has one constant, is written in corresponding-states form with a contact correction, and is fitted with LOTO on the 11 UF stations.

| rule | constant | UF LOTO RMSE K | NBS LN2/Cu MAE K | predicted LN2/Cu 1–5.9 atm | predicted LH2/Cu 1–11 atm |
|---|---|---|---|---|---|
| R0 T* = T_c | – | 7.6 | 6.6 | 126.9 → 126.6 | 33.3 → 33.2 |
| **R5 T* = T_hn + e T_c** | e = 0.122 | **6.5** | 4.5 | 130.6 → 131.4 | 34.3 → 36.6 |
| R1a T* = T_hn + a(T_c − T_hn) | a = 1.38 | 6.7 | 4.4 | 131.4 → 130.6 | 34.5 → 33.4 |
| R2 T* = Tsat + k(T_c − Tsat) | k = 1.11 | 6.7 | 4.4 | 132.3 → 129.9 | 34.7 → 33.3 |
| R3 Henry on T_hn | A = 0.43 | 7.1 | 11.3 | 122.3 → 120.5 | 43.4 → 36.7 |
| R4 T* = Tsat + d T_c | d = 0.325 | 8.7 | 5.7 | 119.0 → 137.9 | 31.3 → 43.1 |
| Kalinin / Simon–Simoneau | – | MAE 29 / 22 | 39 / 22 | 95.6 / 111.2 | 24–32 / 29–31 |

**Pressure slopes dT_MFB/dTsat:**
- UF, all stations: 0.17 [−0.83, 0.88].
- UF, P > 370 kPa: 1.22 [0.61, 3.25]. These are only tests 8–10.
- NBS LN2 knees: 0.16.

A rule following T_hn rises only weakly below T_r ≈ 0.8. The steep high-P rise would need a Tsat-anchored rule (R4), which the low-P UF stations and the NBS knees disfavour.

**NBS knees.** I took the lumped copper-wall flux m''cp_Cu(T)(−dT/dt) on the linear segments between digitized markers. The knee is where q first exceeds 2× the film-boiling q before its peak.
- LN2: 6 well-resolved traces give 124.5–141.6 K (median 133 K).
- LH2: not resolvable. The copper cp is tiny below 45 K, so the transition takes < 1 s, comparable to the digitization step. The plateau readings (13–38 K) also show TC errors of about ±10 K at LH2 temperatures.

**Literature.** Darr & Hartwig 2020 also could not extract the LH2 rewet temperature (their correlation gives 31–33 K).

**Recommendation:** R5. It ties for best on UF and NBS-LN2, its contact temperature rises with P, it transfers to other fluids through T_c and T_hn(T_r), and it has one constant. It replaces the Phase-2 T_c rule.

---------------------------------------------------------------------------------------------------

# Phase 6: parsimony (is the closure over-fit?)

Scripts: `a18_parsimony.py` (ladder, LOTO, surrogate) and `a19_parsimony_report.py`. Outputs: `out/t18_parsimony.md` (full inventory and ladder table), `out/parsimony_constants.json` (fold constants), and `figs/fig21_parsimony.png/.pdf`.

**Inventory.** The final closure contains:
- **3 constants fitted to the UF data:** film C = 15.7 (renewal intensity), Λ = 0.098 m ≈ 8.4 D (development length), and the rewet offset e = 0.122.
- **1 shape constant fitted to the per-TC transition data:** n = 0.5. It is insensitive.
- **Fixed by scaling theory:** the 1/3 exponents. The free fit gives 0.34.
- **Modelling choices:** u_l = G(1−x)/ρ_l and the property states.
- **Literature, unchanged:** Lienhard T_hn, contact β with NIST properties, Cooper, Dittus–Boelter, Katto–Ohno, Bromley gravity term, and the vapour Dittus–Boelter guard.
- **Rig boundary conditions, not closure:** τ = 2.5 s, q_par, and the 0.735 data correction.

For comparison, the fitted-constant counts of related closures are: Darr 2016 low g, 7 (6 active); Darr 2019 UF universal, 21; Darr–Hartwig 2020 LH2 film boiling, 11 (plus adopted rewet and nucleate-boiling sets); Darr 2019b rewet, 2.

**Ablation ladder.** Every fitted constant is re-fitted in each leave-one-test-out fold, and the 0-D surrogate applies the fold constants to the held-out test.

| model | constants fitted to UF | h held-out RMS ln (Δ per constant) | station-bias σ | rewet held-out MAE K (NBS LN2/Cu knees) | surrogate RMS K / t130 MAE s / t100* MAE s / stations missing a quench (130 K, 100*) |
|---|---|---|---|---|---|
| M0 Bromley 2.7 + T_hn contact | 0 | 0.512 | 0.30 | 16.5 (17.5) | 31.0 / 7.9 / 8.5 / 3, 5 |
| M1 C·Pe^(1/3) + T_hn contact | 1 | 0.245 (−0.27) | 0.17 | 16.5 (17.5) | 13.4 / 1.63 / 2.50 / 1, 2 |
| M2 + rewet offset e | 2 | 0.245 | 0.17 | **5.4 (4.5)** | 13.7 / 1.68 / 2.61 / 1, 2 |
| **M3 + Λ (final)** | 3 | **0.240** (−0.004) | 0.17 | 5.4 (4.5) | **13.0 / 1.38 / 2.43** / 1, 2 |
| M4 + free Pe and P exponents | 5 | 0.296 (+0.03) | 0.24 | 5.4 (4.5) | 15.3 / 1.55 / 2.95 / 2, 4 |
| M5 Darr 2016, c1–c6 refit | 7 | 0.689 (+0.20) | 0.62 | 5.4 (4.5) | 21.3 / 2.94 / 3.74 / 3, 6 |

For reference, Darr 2016 with its published constants (0 fitted) gives 0.315, and Darr × C (1 fitted) gives 0.253.

**Information criteria.** These cover the film part only. I use n = 20 test-stations with station-mean residuals: the within-run residual has lag-1 autocorrelation of 0.97, and the error is dominated by run-level bias. The AR(1) effective n of about 70–77 gives the same ranking.
- AIC: M1 −72.6, M3 −72.7, M4 −69.5, M5 −61.6.
- BIC: M1 −71.6, M3 −70.7, M4 −65.5, M5 −55.6.

**Where the returns go.**
- **The first constant (C) buys almost everything:** the held-out h error drops 2×, and the quench-time error drops from 7.9 s to 1.6 s. M1 is already within 0.07 of the within-run noise floor (0.17); the rest is run-to-run scatter (σ ≈ 0.17).
- **The rewet offset e** is required for the rewet temperature: 16.5 → 5.4 K on UF, and 17.5 → 4.5 K on the independent NBS LN2/Cu knees. Its fold values are 0.115–0.132. It barely changes wall-T metrics on thin steel walls.
- **Λ** is cheap and stable (fold range 0.081–0.118 m). It is neutral under BIC, gives −0.004 held-out log error, and −0.25 s in quench-time MAE. Its role is kinematic: it gives the front speed (Phase 3), which pointwise statistics do not reward.
- **Beyond 3 constants, held-out performance gets worse.**
  - M4: the P exponent is collinear with G, and fold C ranges 11–28.
  - M5: the droplet term collapses (c4 → 0) and the We exponent jumps to about 5. A 7-constant form is over-parameterised for 10 tests.

**Transferability.**
- C is a renewal intensity. It is the first knob to re-calibrate for a new fluid, geometry or gravity level, preferably from one wall-temperature trace.
- Λ ≈ 8.4 D is a development length. Expect it to scale with D. Re-check it only if station-resolved data exist; otherwise keep 8.4 D, or drop it to get M2, which costs about 0.3 s in quench time and the front speed.
- e (15.4 K for N2) is an apparent-onset offset above T_hn. Keep it; it is confirmed on Cu/LN2 without refitting.
- The literature constants and exponents should not be re-fitted. The engineering recipe is the M3 structure with C re-calibrated per system.
