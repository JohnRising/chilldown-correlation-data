# Final wall-to-fluid closure for 1-D low-g cryogenic chilldown (OpenFLUME)

This closure was derived a priori from Darr et al. 2016 (npj Microgravity 2:16033): low-g LN2 flowing through a 304SS tube with ID 11.68 mm. Supporting scripts are `a07`–`a10` and the evidence is in `REPORT.md` §Phase 2.

It is **gravity-free**. Only **one constant is fitted** (the film-boiling C); every other piece is theory or a literature correlation used without refitting.

## 0. Node variables used

| symbol | meaning |
|---|---|
| P | local pressure |
| T_f | local fluid temperature |
| x_e | equilibrium quality of the node |
| G | local mass flux, magnitude, in kg/m²s |
| D | inner diameter |
| z | distance from the line inlet (heated length) |
| T_w | wall inner-surface temperature; for a thin wall, the wall node temperature |

Material: wall effusivity e_w(T) = √(k_w ρ_w c_w).

Saturation properties, all at P via CoolProp (subscript l = saturated liquid, v = saturated vapour): T_sat, ρ_l, μ_l, k_l, cp_l, h_fg, σ.

Film-temperature vapour properties are those of superheated vapour at (P, T_film), with T_film = (T_w + T_sat)/2. They are denoted ρ_vf, k_vf, cp_vf.

Other quantities:
- ΔT = T_w − T_sat
- T_l = min(T_f, T_sat)
- p_r = P/P_c
- M = molar mass in kg/kmol (28.013 for N2)

## 1. Film boiling (T_w ≥ T_MFB): Lévêque-type renewal law, one fitted constant

    Pe_l = (G/ρ_l) · D · ρ_vf · cp_vf / k_vf
    h_FB = C · (k_vf / D) · Pe_l^(1/3),      C = 8.55
    q_FB = h_FB · ΔT

**Where C comes from.**
- It is fitted to the consistent wall-energy-balance heat-flux dataset: 20 case-stations, 4422 film-boiling samples.
- The test-section mass flux is taken as the metered value 2.5 s later: G_ts(t) = G_meter(t + 2.5 s).
- Leave-one-case-out (LOCO) cross-validation gives C between 8.27 and 8.80 (±3 %).

**Accuracy.**
- RMS of ln(pred/meas): 0.235 in-sample, 0.246 under LOCO.
- Case-to-case scatter of the bias: σ_ln = 0.17. This is the realistic ±1σ for a single new run.
- 83 % of points fall within ±30 %.
- With an unlagged meter signal (τ = 0): C = 8.63, RMS 0.267.

**Equivalent dimensional form:** h_FB = 8.55 · k_vf^(2/3) · (ρ_vf cp_vf G/ρ_l)^(1/3) · D^(−2/3). This gives h ∝ G^(1/3) P^≈0.38 and no explicit ΔT dependence.

**Optional entrance factor.** Multiply by [1 + (D/(z+0.1D))^0.7] (Darr 2016) and use C = 7.60.
- LOCO improves only from 0.246 to 0.240.
- It removes about 8 % of the measured ≈12 % upstream/downstream difference.
- It is not required, and it needs a distance-from-front definition.

**Quality.** No quality factor is needed for x_e ≤ 0.15. That covers the data range, which is mostly ≤ 0.05; on the x_e ≤ 0.05 subset C = 8.50.

**Guards outside the data range (untested).**
- h_FB ≥ 4.36 k_vf/D (a laminar-vapour floor as G → 0).
- For x_e ≥ 0.15, use h = max(h_FB, 0.023 Re_vx^0.8 Pr_v^0.4 k_v/D), with Re_vx = G x_e D/μ_v.
- For x_e ≥ 1, use single-phase vapour Dittus–Boelter with bulk properties.

**Data range covered:**

| quantity | range |
|---|---|
| G | 12–700 kg/m²s |
| P | 90–770 kPa |
| ΔT | 40–220 K |
| z/D | 12.8 and 34.3 |
| gravity | 0 g (±0.02 g) |

For 1 g, add a buoyancy (Bromley-type) term. It is not part of this closure.

## 2. Rewet / minimum film boiling: contact temperature equals T_c (zero constants)

    β     = e_l / e_w(T_MFB),   e_l = √(k_l ρ_l cp_l) of saturated liquid at P
    T_MFB = T_c + β (T_c − T_l)

- The equation is implicit only through e_w(T_MFB). Three fixed-point iterations from T_MFB = T_c + 0.1(T_c − T_l) converge to < 0.01 K.
- Valid for P < P_c.
- Physical basis: on brief liquid–wall contact, the contact temperature (e_w T_w + e_l T_l)/(e_w + e_l) must fall below the thermodynamic critical temperature for liquid to wet the wall.
- For 304SS with LN2, β = 0.08–0.09 and T_MFB = 128.2–130.8 K.

**Comparison with the measured T_MFB** (11 stations; onset defined as the point where the cooling rate doubles):

| | mean error | MAE |
|---|---|---|
| all stations | −4.6 K | 5.6 K |
| downstream station (z/D = 34) | +0.0 K | 2.2 K |
| upstream station (z/D = 12.8) | −8.4 K | – |

The upstream offset of +6 to +13 K is robust to how rewet onset is defined; see REPORT. On the 1-g Darr tests 31 and 40, with the same onset definition, the rule gives 129.3–129.4 K against 128.5–129.9 K measured.

**Optional:** add ΔT_ent ≈ +8 K for z/D ≲ 15. This is a single fitted constant and is not cross-validated.

**Coated walls.** Apply the rule at the coating surface node with e_w = e_coating, where the contact time is shorter than δ_c²/α_c. For FEP this is 2 ms for one layer and 30 ms for four layers.

## 3. Nucleate boiling (T_w ≤ T_CHF)

Cooper (1984) pool boiling, solved explicitly for q, plus a single-phase liquid floor:

    A     = 55 · p_r^0.12 · (−log10 p_r)^(−0.55) · M^(−0.5)        [SI; R_p = 1 µm]
    h_DBl = 0.023 · Re_l^0.8 · Pr_l^0.4 · k_l/D,   Re_l = G(1−x_e) D/μ_l,  Pr_l = cp_l μ_l/k_l
    q_NB  = h_DBl · ΔT + (A · ΔT)^3          (ΔT ≥ 0;  for T_w < T_sat: q = h_DBl (T_w − T_f))

## 4. CHF: Katto–Ohno (1984), saturated inlet, heated length L = z

    R = ρ_v/ρ_l,  l = z/D,  W = σ ρ_l /(G² z)
    C_K = 0.25 (l<50); 0.25+0.0009(l−50) (50≤l≤150); 0.34 (l>150)
    q1 = C_K W^0.043 / l
    q2 = 0.10 R^0.133 W^(1/3) / (1+0.0031 l)
    q3 = 0.098 R^0.133 W^0.433 l^0.27 / (1+0.0031 l)
    q_CHF/(G h_fg) = q1 if q1<q2;  else (q2 if q2<q3 else q3)        (R < 0.15 branch)

**Superheat at CHF.** ΔT_CHF is the root of q_NB(ΔT) = q_CHF; Newton iteration converges in a few steps.

**Ordering constraint.** Enforce T_CHF ≤ T_MFB − 1 K and q_CHF ≥ q_FB(T_MFB).

**Accuracy.** Katto–Ohno over-predicts the transient quench peak by 1.66×. Because the sensitivity to CHF is negligible (§7), it is kept un-fitted.

## 5. Transition boiling (T_CHF < T_w < T_MFB)

    θ    = (T_MFB − T_w)/(T_MFB − T_CHF)
    q_TB = q_FB(T_MFB) + [q_CHF − q_FB(T_MFB)] · θ^n,    n = 0.5

- q_FB(T_MFB) is h_FB evaluated with film properties at (T_MFB + T_sat)/2, multiplied by (T_MFB − T_sat).
- n = 0.5 is the best fit to the per-thermocouple transition data (RMS in ln q of 0.19, against 0.60 for linear interpolation).
- n = 1 may be used if the infinite slope at θ → 0 causes solver difficulty. The effect on wall-temperature metrics is ≤ 0.1 K and ≤ 0.5 s.

## 6. Single expression for the wall heat flux

    q''(T_w) = 1[T_w ≥ T_MFB]        · h_FB(T_w) (T_w − T_sat)
             + 1[T_CHF < T_w < T_MFB] · { q_FB(T_MFB) + [q_CHF − q_FB(T_MFB)] θ^n }
             + 1[T_w ≤ T_CHF]        · q_NB(T_w)

- q'' is continuous at T_MFB and at T_CHF by construction.
- Optional smoothing for implicit solvers: replace 1[·] by logistic functions of width 1 K.
- Regime switching must be done per node, on the node's own T_w.

## 7. What matters (0-D surrogate, 20 station records, `a10`)

Baseline: the full closure above against station-mean TC histories.

| metric | baseline |
|---|---|
| RMS(T), mean / median | 13.1 / 9.3 K |
| time-to-130 K error, MAE (bias) | 1.8 s (+0.5) |
| time to T_sat+8 K (or 100 K) error, MAE (bias) | 1.6 s (−0.3) |

One-at-a-time sensitivities:

| sub-model change | effect |
|---|---|
| film h ×0.7 | no station quenches; RMS +6.7 K |
| film h ×1.3 | quench 4.5 s early; RMS +15 K |
| T_MFB ±5 K | Δt130 ≤ 1.1 s, Δt100 ≤ 1.3 s, ΔRMS ≤ 0.8 K |
| CHF ×0.5–1.5, transition n 0.25–1, NB ×0.5–2, CHF = 100 kW/m² constant | ΔRMS ≤ 0.12 K, Δt130 = 0, Δt100 ≤ 0.5 s |

**Conclusion:** the film-boiling law dominates, T_MFB is second, and CHF, transition and nucleate boiling can be crude.

## 8. Variants

- **Recommended:** §1–6 as written. One fitted constant, C = 8.55.
- **Literature alternative for film boiling:** h = 0.277 √(k_v cp_v G/D), with saturated-vapour properties. LOCO C range 0.269–0.288; RMS 0.259 and LOCO 0.270 at τ = 2.5 s; RMS 0.306 at τ = 0.
  - The literature value C = 0.218 was fitted to the published heat flux, which is 26 % low in 5 of the 10 cases.
- **Theory-only (zero fitted constants):** acceptable for §2–5.
  - It is **not** acceptable for film boiling. Bromley forced convection with C = 2.7 is 40 % high, which makes the surrogate quench about 8 s early. No physically derived film law came within 40 % without a fitted constant.
  - A zero-fit closure would therefore have to use Bromley C = 2.7; it is not recommended.

## 9. Rig-specific boundary conditions (Darr 2016 rig only, not closure)

- **Mass flux:** G_test-section(t) = G_meter(t + 2.5 s).
- **Wall parasitic heat:** q_par = 98.6 + 2.2 (293 − T_w) W/m² per inner area. It is fitted to the authors' inverse model on the self-consistent cases.
- **Wall properties:** areal mass m''_w = 4.20 kg/m² with NIST 304SS cp(T).
- **TC comparison:** TCs sit on the outer wall. T_outer − T_inner = q'' · 0.254 mm / k, i.e. ≤ 0.3 K in film boiling and about 2.5 K at the CHF peak.
- **Published heat flux and h:** those of cases 1, 2, 3, 5 and 6 must be multiplied by 1/0.735.

## 10. Position / front dependence (Phase 3 — supersedes the film law of §1 as the recommended default)

**Recommendation.** Use a developing-layer factor based on a static axial position, plus liquid-core velocity. **No dynamic front tracking is needed.**

    u_l   = G (1 − x_e) / ρ_l                              (liquid-core velocity; x_e clipped to [0, 0.9])
    Pe_l  = u_l D ρ_vf cp_vf / k_vf
    L_eff = z Λ / (z + Λ),      z = max(z_node_centre, D)  (distance from the start of the cryogen-wetted line)
    h_FB  = C (k_vf / D) Pe_l^(1/3) (D / L_eff)^(1/3)
          = C k_vf^(2/3) (u_l ρ_vf cp_vf)^(1/3) D^(−1/3) (1/z + 1/Λ)^(1/3)
    C = 15.8 (LOCO 14.9–16.5),   Λ = 0.101 m ≈ 8.6 D (LOCO 0.080–0.121 m)

Both constants are fitted. The (1 − x_e) term adds no constant.

**Physical picture.** A Lévêque thermal layer in the vapour gap, sheared over a length of order D, develops from the inlet leading edge. Interfacial renewal caps its development length at Λ.

**How the factor behaves:**
- At the two stations, (D/L_eff)^(1/3) equals 1.17·(D/Λ)^(1/3) upstream (z/D = 12.8) and 1.07·(D/Λ)^(1/3) downstream (z/D = 34). The resulting up/down ratio is 1.10; measured is 1.12.
- For z ≫ Λ, i.e. long lines, the law saturates to 7.7 (k_vf/D) Pe_l^(1/3), within 10 % of the Phase-2 law. It therefore stays finite and transferable to long lines.

**Accuracy and comparison.**
- In-sample RMS of ln(pred/meas) is 0.230 and LOCO 0.240 (Phase 2: 0.246). The case-station bias σ falls from 0.167 to 0.151, and the mean up/down station bias falls from −0.115 to 0.
- Saturating at Λ, an inlet offset (L = z + z0, z0 = 0.46 m), a power law z^−0.12 and Darr's entrance factor all give the same LOCO (0.240). Two stations cannot discriminate their shapes. Only the saturating form stays finite for long lines.
- A pure L = z law over-corrects (up/down ratio 1.39; LOCO 0.262).
- A fully D-free penetration form √(k ρ cp u_l/(z+z0)) fits poorly (LOCO 0.357). The shear length D is kept.

**Effect in the 0-D surrogate** (a13, against station-mean TCs):

| | Phase-2 law | §10 law |
|---|---|---|
| apparent front speed | 1360 mm/s | 114 mm/s (measured 96) |
| quench delay downstream − upstream at 130 K | 0.1–0.2 s | 1.9–3.5 s (measured 2.4–5.4 s) |
| t130 MAE | 1.78 s | 1.10 s |
| t130 bias, upstream / downstream | +2.1 / −1.1 s | +0.9 / −0.0 s |
| RMS(T), mean | 13.1 K | 11.8 K |
| t100* MAE | 1.62 s | 1.87 s |

**Why there is no dynamic front.** The test was whether, when the upstream station rewets, the downstream film h jumps as if its development length had restarted at the front.
- Measured jump in ln(h/law): +0.03 ± 0.045 (6 cases). The placebo spread at other times is σ = 0.11.
- A front-restarted pure-z law predicts +0.155, which is rejected at about 2.8σ.
- The saturating law predicts +0.045, indistinguishable from inlet-referenced L.
- So L referenced to the static z is sufficient.

**What sets the front speed.**
- The measured "front" at 44–106 mm/s is **kinematic (dispersed)**: each node reaches its own T_MFB at a different time.
- The kinematic estimate of the delay, ΔT_wall/rate + ΔT_MFB/rate, gives 2.8–3.8 s against 2.8–3.9 s measured for cases 7–10. Case 6 is predicted 8.0 s against 5.7 s measured.
- Conduction-controlled (Yamanouchi) propagation would be only 4–6 mm/s. The liquid-supply limit is ≥ 680 mm/s. Neither sets the speed.
- A 1-D code reproduces the front only if h, and optionally T_MFB, vary with z. Keep axial wall conduction for numerical smoothness; it does not set the speed.

**Rewet entrance increment** (+8 K upstream): **do not use it together with §10.** It over-corrects the front (65–69 mm/s; upstream t130 −0.6 s; t100 MAE 2.5 s). Keep the zero-constant T_c rule of §2.

**Low-G cases 1 and 5** (end 25–40 K too cold) are **not a late-time mechanism**:
- The measured h is 0.67–0.90 of the law from the start of the run. The ratio of late-to-early h relative to the law changes by only −0.0 to +0.23.
- G is steady in case 5 (coefficient of variation 0.07, x_e ≈ 0.02–0.04). No flow reversal or pressure collapse was found.
- Case 5 has about 25 % lower h than case 4 at the same G and P. That is run-to-run scatter of about −1.5σ to −2σ of the 0.17 case σ, accumulated over 24 s.
- The (1 − x_e) term recovers about 5 K in the downstream stations of cases 1 and 2 (x_e up to 0.19), but it slightly worsens the upstream stations of cases 1 and 5.
- No further physically justified fix is supported. Expect ±17 % (1σ) run-to-run h uncertainty.

## 11. Rewet rule (revised; Phase 5, final) — supersedes §2 as the default

### Recommended rule: R5, a fixed reduced offset above the superheat limit (1 fitted constant)

    Tsat  = Tsat(P)                     ofchill.props: fluid_props("Nitrogen"|"Hydrogen").Tsat(P)   (no new fits needed)
    Tr    = Tsat / Tc
    T_hn  = Tc · (0.905 + 0.095 · Tr^8)          Lienhard (1976) superheat limit, corresponding-states form
    T*    = T_hn + e · Tc,      e = 0.122        (fitted on the 11 UF low-g stations)
    T_MFB = T* + β · (T* − T_l),  T_l = min(T_f, Tsat)

**Constants.**
- T_c(N2) = 126.19 K and T_c(H2, normal) = 33.145 K.
- The offset e·T_c is 15.4 K for N2 and 4.04 K for H2.
- β = √(kρc)_l,sat / √(kρc)_wall, both evaluated at about T_MFB. A constant per fluid/wall pair is adequate:

| fluid / wall | β |
|---|---|
| LN2 / 304SS | 0.085 |
| LN2 / Cu | 0.012 |
| LH2 / Cu | 0.014–0.02 |
| LH2 / 304SS | 0.25 |

  In a code, either take β from `ofchill.closures.contact_beta_spec(P, fluid, wall)[0]` (wall "ss304" or copper; it evaluates e_w at the §2 T_MFB, a difference of < 0.1 K) or use these constants.

**Equivalent explicit form for 1-D codes:** T_MFB = (1+β)(T_hn(P) + 0.122 T_c) − β · min(T_f, T_sat(P)).

### Predictions

| system | T_MFB |
|---|---|
| LN2 / 304SS, 1 / 3 / 7.5 atm | 134.6 / 133.8 / 133.8 K |
| LN2 / Cu, 1–5.9 atm | 130.6 → 131.4 K |
| LH2 / Cu, 1 → 11 atm | 34.3 → 36.6 K |

For LH2/Cu, NBS line pressures give 34.4–35.5 K.

### Why this rule

Six one-constant rules were compared (a16, a17). Errors are in K.

| rule | UF 11 stations, LOTO RMSE | NBS LN2/Cu knees (n=6), MAE (bias) | P-slope dT_MFB/dTsat (UF range) |
|---|---|---|---|
| R0 T* = T_c (§2, 0 constants) | 7.6 | 6.6 (−5.9) | −0.11 |
| **R5 T* = T_hn + 0.122 T_c** | **6.5** | 4.5 (−1.9) | −0.04 for SS (T* itself rises, +0.06…+0.15) |
| R1a T* = T_hn + 1.38(T_c − T_hn) | 6.7 | 4.4 | −0.14 |
| R2 T* = Tsat + 1.11(T_c − Tsat) | 6.7 | 4.4 | −0.23 |
| R3 Henry on T_hn (A = 0.43) | 7.1 | 11.3 (−11) | −0.42 |
| R4 T* = Tsat + 0.325 T_c (constant reduced superheat) | 8.7 | 5.7 (−5.7) | +0.96 |
| Kalinin 1969 / Simon–Simoneau 1969 (no fit) | MAE 29 / 22 | 39 / 22 | – |

**Measured pressure slopes:**
- UF, all stations: 0.17 [−0.83, 0.88].
- UF, P > 370 kPa only (tests 8–10, 6 stations): 1.22 [0.61, 3.25].
- NBS LN2 knees against local Tsat: 0.16.

**Interpretation:**
- R5, R1a and R2 are statistically equivalent on both datasets.
- R5 is the only member of that group whose contact temperature increases with P, following T_hn, as physics suggests.
- For N2 below T_r ≈ 0.8, that increase is small: dT_hn/dTsat ≈ 0.06–0.15.
- The steep rise seen in the three high-P UF tests could only be reproduced by a Tsat-anchored rule (R4). R4 is rejected by the low-P UF stations (tests 3 and 6) and under-predicts the NBS LN2 knees.

**What R5 is and is not:**
- It is an empirical contact-temperature criterion.
- T* exceeds T_c by about 4 K for N2, so the measured T_MFB is an apparent onset temperature (outer-wall TC, cooling-rate threshold). It is not the thermodynamic contact limit.
- Upstream/downstream bias on UF is −3.8 / +4.6 K.

### Discriminating data

**NBS 9264 LN2 copper line.**
- Knees from flux-jump detection on the digitized segments, using the copper-cp lumped flux (`out/t17_nbs_knees.md`, fig20): 124.5–141.6 K, median 133 K, at 6 well-resolved station traces (segment resolution 1.5–17 K).
- This rules out Kalinin (95 K) and Simon–Simoneau (111 K).
- R0 is 6 K low. R5, R1a and R2 are within about 2 K on average.

**NBS 9264 LH2.**
- **Not resolvable.** At T < 45 K copper cp is tiny (26 J/kg/K at 30 K), so the wall drops from rewet to the plateau in < 1 s. That is comparable to the 0.6 s resampling and the integer-second click rounding.
- No flux jump is detectable in any of the 24 traces.
- The post-quench plateau readings are 13–38 K, against T_l ≈ 19.5–29 K. This is consistent with the report's warning that the Cu-constantan/LN2-reference TCs are unreliable at LH2 temperatures, so errors of ±5–15 K are possible.
- The NBS LH2 data therefore cannot discriminate 33–36 K (R0/R5) from Kalinin's 24–30 K at local P.

**Literature.** Darr & Hartwig (2020) state that the LH2 rewet temperature "is difficult to determine from the LH2 chilldown data". They use Darr 2019, which gives 31–33 K.

**A decisive test** would need LH2 quench data on a stainless tube with resolved TC knees (GRC; Hartwig 2019/2020), where β ≈ 0.25. There R5 predicts 36.6–36.8 K, Kalinin 32.6 K and Simon–Simoneau 29.8 K. Coated-tube data (FEP surface; R5 gives 165 K at 300 kPa) would also discriminate.
