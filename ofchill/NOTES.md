# ofchill — solver facts (OpenFLUME pip 0.3.1 = git tag v0.3.1)

Everything here was read from `git show v0.3.1:<path>` in /Users/john/Code/OpenFLUME
and, where marked **[measured]**, checked numerically with the pip package.

## 1. When and how convection h is evaluated

Source: `src/core/solver/step.ts`, `src/core/solver/conductorH.ts`,
`src/core/correlations/dispatcher.ts`, `src/core/transient/resultRecorder.ts`.

One time step = one *attempt* (fixed stepping: no rejection) = an **outer Picard
loop** (up to `maxOuter` = 1000) around an inner Newton solve of the fluid
[P, ṁ, h] system, followed by the wall (solid) thermal solve.

| moment | what happens to h |
|---|---|
| attempt start | `computeConductorHMap(ctx, state, undefined, t)` — a **fresh, un-relaxed** evaluation at the warm-start state (= previous accepted state: previous fluid state **and previous wall T**). `t` = step END time. |
| inside inner Newton | h is **frozen** (the fluid energy rows use h·A·(Tw−Tf) with that frozen h). |
| after each outer's inner Newton | `computeConductorHMap(ctx, state, prevMap, t)`: h_raw evaluated at the NEW fluid state but the **wall T of the previous outer** (the h-map is refreshed *before* `solveThermalSubsystem`), floored, then **under-relaxed**: `h = h_prev + H_RELAX·(h_raw − h_prev)`, `H_RELAX = 0.5` (`correlations/types.ts`, not configurable). The wall solve of this outer uses the relaxed h. |
| wall solve | implicit (backward Euler) in Tw with h FIXED (enthalpy-form storage, exact cp curve). |
| recorded result `heatTransferCoeff` | a **separate fresh, un-relaxed** evaluation at the accepted end-of-step state (`resultRecorder.ts`), NOT the (relaxed) h actually used in the step. `heatRate`/`heatFlux` are recomputed with that recorded h and the end-of-step temperatures. So Σ heatRate·dt is not exactly the energy the solver exchanged (energy diagnostic in runner.py shows this as `wallHeat_over_wallEnergy_m1`, typically −0.4…−2 %). |

Consequences:
* h is relaxed **only across outer iterations within a step**; there is no
  relaxation or lag across steps (each step restarts from a fresh h at the
  previous accepted state).  If the outer loop converges, h is effectively
  implicit (h = H(state at end of step)).
* Outer convergence test: `maxDeltaT < fluidTol` where maxDeltaT = max(|Δh_fluid|
  [J/kg], |ΔT_fluid|·1e4, **|ΔT_wall| [K]**) and fluidTol = tol·1e6 = **10** when any
  node is two-phase (tol·1e7 = 100 otherwise; tol = settings.tolerance = 1e-5).
  The wall-T part of the test is therefore very loose (10 K / 100 K per outer):
  convergence is governed by the fluid enthalpy change and by the inner-Newton
  scaled-residual bar (< tol·1e3).  h itself is not tested.
* In practice **[measured, case 8, Miropolskii, FN_DEBUG_OUTER=1]** steps take
  4–25 outers, and some steps run to **maxOuter = 1000** (dome-edge fluid limit
  cycles; fluid enthalpy relaxation hRelax = min(relaxation, 0.1) whenever any
  node is two-phase).  These 1000-outer steps dominate run time and make it very
  sensitive to inputs: case 9 Miropolskii 11 s with the prior xls-precision
  schedules vs 185 s with the CSV schedules; named DB (h at floor in the dome)
  629 s (case 8) and 1034 s (case 9).  The steps still end "converged" (best
  outer iterate restored, scaled residual < bar) — `stepResidualsScaled` stays ~1e-9.
* With a step-function boiling curve, the fixed-point iteration
  h_{k+1} = ½h_k + ½H(Tw(h_k)) switches regime only once the end-of-step Tw with
  the OLD regime's h crosses the threshold, i.e. switching is effectively
  explicit (≤ 1 step late, 0.1 s here).  **[measured, case 8, toy `boil_demo`,
  33× jump]**: sharp step: converged, 0 stalled steps, 10 s, one regime switch
  per conductor, no chatter; smoothstep widths 10/40 K: 20/67 s, the quench
  time moves earlier (14.0 → 13.3 → 11.6 s) because the effective rewet
  temperature is the bottom of the smoothing band.  Note h in the recorded
  results is fresh at end-of-step, so a recorded jump happens at the step whose
  END state crossed the threshold.
* Because Tw in the h evaluation lags by one outer (and the wall solve uses the
  relaxed h), a closure where h depends strongly on Tw is solved by Picard with
  a 0.5 damping — a limit cycle is possible if |dH/dTw · dTw/dh| > 3; it would
  show up as many outers (run time) rather than as a flagged failure.

## 2. Floor and fallback

`evaluateConvectionH`: after the model, `floor = cond.type.h ?? FALLBACK_H_FLOOR
(5 W/m²K)`; `if (hRaw < floor) hRaw = floor` (counted by `recordHFloorClamp`),
then relaxation.  Any failure (named: missing property / non-realFluid; custom:
parse error, unknown identifier, wrong type, non-finite result) returns the
*fallback* = `cond.type.h ?? 5`.  So:
* a closure cannot return h < 5 W/m²K unless the conductor sets `type.h` lower
  (then `type.h` is BOTH floor and fallback);
* zero / negative / NaN / ±Inf h → 5 W/m²K silently (no error, no warning);
* named `dittusBoelter` inside the dome: statePH gives **no cp** → floor (5).
  Named `miropolskii`: dome → Miropolskii with x clipped to [0.01, 0.99]
  (sat-vapour μ, k, cp, sat ρ_l, ρ_v); single phase → DB with node props.

## 3. Custom expression scope (`correlations/customCorrelation.ts`)

| name | present when | value |
|---|---|---|
| t | always | step END time at refresh (0 at t=0 recording) |
| Tf, P | always | node T, P (state; fallbacks 300 K / 1e5 Pa never seen) |
| Tw | conductor has a solid endpoint | wall node T (**previous outer's** value during a step) |
| area | always | conductor area |
| D, flowArea | if `diameter` or `flowArea` configured (one derived from the other) | |
| G | if flowArea known | ½Σ\|ṁ\| over branches at the node / flowArea (uses current-iterate ṁ) |
| rho, mu | realFluid, property call OK | from statePH(P, h_node) (P,h clamped). **In the dome:** rho = HEM mixture, mu = McAdams mixture 1/(x/μg+(1−x)/μf) |
| k | realFluid | single phase: CoolProp; **dome: mixture 1/(x/kg+(1−x)/kf)** (the RESULT reports thermalConductivity = NaN in the dome, but the scope has k) |
| cp | realFluid **single-phase only** | **ABSENT in the dome** |
| Pr | needs cp, mu, k | **ABSENT in the dome** |
| Re | needs G, D, mu | dome: G·D/μ_mixture |
| quality | **dome only** (statePH phase twoPhase) | **ABSENT for liquid/vapour nodes** |
| param('n') / params.n | always | correlation.params |

**An absent identifier makes the whole expression fail → fallback 5 W/m²K**
**[measured: test_absent_identifier_falls_back_to_floor]**.  There is no
"is defined" test.  The ternary `c ? a : b` and `&&`/`||` short-circuit, so an
identifier may be referenced only in a branch that is never taken when it is
absent.  ofchill provides `s.in_dome()` / `s.x_eq()` (from rho, P and
saturation-density fits, error 1e-5, margin 2e-4) that never touch cp/quality.

Language: `+ - * / % ^` (`^` right-assoc, binds tighter than unary minus),
comparisons, `&& || !`, `?:`, builtins `min max abs sqrt exp log sin cos tanh
clamp smoothstep`, constant `pi`.  No `pow`, `log10`, `if`, `let`/variables:
shared subexpressions are duplicated in the emitted string (ofchill
expressions are 13–29 kB; evaluation cost is negligible vs CoolProp but see
the run-time note below).  Booleans are not numbers (`(a<b)*2` is an error →
fallback); use `?:`.  `smoothstep(e,e,x)` = step (x<e ? 0 : 1).
The expression is compiled once per solve.

Expression cost **[measured]**: a film-temperature DB using three
superheated-vapour fits (k, μ, cp at T_film) emits 124 kB and costs ~0.3 ms per
evaluation (tiny network: +0.5 s per 400 steps for one conductor).  With
11 conductors × (4–25, sometimes 1000) outers per step this can matter; the
saturation-property fits are cheap (db_vapor: 13.5 kB, no measurable cost).

## 4. Other solver facts used by the builder

* Wall material preset `{material:'stainless-steel-304'}`: NIST log-log fits
  (4–300 K, spliced to ANL above 300 K) sampled into a piecewise-linear table;
  storage uses the exact enthalpy integral.  Conduction k evaluated at the
  mean of the two endpoint temperatures.
* `TransientResult` fields used: `times, nodes[*].{pressure, temperature,
  density, enthalpy, internalEnergy, quality (None outside dome), phase,
  viscosity, thermalConductivity (NaN in dome), specificHeat (NaN in dome)},
  branches[*].mdot, solidNodes[*].temperature, conductors[*].{heatRate,
  heatTransferCoeff, heatFlux}, stepResiduals, stepResidualsScaled
  (≥ 0.01 = stalled step), converged`.  `stats` is empty for fixed stepping.
* darrHartwig / ttWf need `axialPosition` (and ttWf `segmentLength`) per
  conductor and a solid wall endpoint; the builder fills them (x_j, L/N·frac).

## 5. Python API quirks (openflume 0.3.1)

* No solve timeout.  ofchill arms a `threading.Timer` that calls
  `session._kill()` (private); the blocked request then raises `BridgeError`
  and the session restarts on next use.
* One Node process per `Session`; the default session is per process, so
  `ProcessPoolExecutor` workers each get one (max_workers = 6 here).
* Node console output (e.g. `FN_DEBUG_OUTER=1`, set in `os.environ` before
  the session starts) arrives on the `openflume.node` logger at INFO.
* Result lists use `None` for missing values (quality) and NaN for others;
  `TransientResult.series()` returns object arrays when None is present —
  ofchill converts with `_arr`.

## 6. Data / reproduction notes

* **Reproduction of the prior Miropolskii N=10 dt=0.1 runs is bit-identical**
  (0.0 K wall-T / h / ṁ difference, cases 7, 8, 10; `python -m ofchill.repro_prior`,
  results/repro/repro_vs_prior.json) when the prior agent's input JSON (xls
  precision) is used.  Case 9 differs by ≤ 0.22 K: the prior branch
  (OpenFLUME-uf = v0.3.1 + 93 commits) contains 436fdc5 "floor the P-h clamp at
  the compressed-liquid minimum" (clampToValidPH); case 9 has subcooled-liquid
  node states.  Re-running the prior TS script on an exported copy of the
  branch with 436fdc5 reverted reproduces the pip result exactly (0.0 K).
* `db_vapor` (custom expression) vs the prior patched-dispatcher DB runs:
  ≤ 0.04 K wall-T difference, RMS identical to 0.01 K (cases 7–10).
* `miropolskii_expr` (expression re-implementation of the named model, fits
  instead of CoolProp sat props, density-based dome test) vs named: ≤ 0.97 K,
  RMS within 0.5 K, but 416 s instead of 56 s (case 8): tiny h perturbations
  change the outer-iteration path (more 1000-outer steps) — run time, not
  accuracy, is the fragile quantity.
* With the repo CSVs the same model differs by up to 11 K (case 8): the CSV
  times are rounded (0.2 vs 0.19999999 s in the xls) and the 0.2-s boxcar
  windows have inclusive edges that fall exactly on 15-Hz samples, so window
  membership (3 vs 4 samples) flips; G knots change by up to 5 kg/m²s.
  ofchill uses `window_tol = 1e-6` (insensitive to that rounding).
  The model is sensitive to such ~3 % schedule perturbations (and run time even
  more so).
* Case 4: P1/P2 collapse to 6–30 kPa at t ≈ 18.8 s (end of test, below the N2
  triple point; CoolProp fails at the boundary) → default end time 18.7 s.
* Energy diagnostic `fluidGain_over_wallEnergy_m1` is relative to the wall
  energy removed, which is small when the model barely cools the wall
  (cases 1, 5, 6: +23…+43 %), so it is only meaningful for cooled cases.

## 7. Phase 2: unified boiling-curve closure (custom expression)

Implementation: `closures.make_boiling_closure(fb, props, tmfb, nb, chf, gravity, vapor_film, **params)`
and the named variants `boiling_variant('a'|'b'|'c_rl'|'c_rv'|'d'|'d119'|'a_sat'|'a_cooper'|'a_thn'|'e_brom'|'spec')`.
The docstring is the formal spec. Notes:

* **Phase test without absent identifiers.** Vapour node ⇔ ρ < ρ_v,sat(P)(1−2e−4) → DB with the node's
  own μ, k, cp (cp exists there). Otherwise x = x_eq from ρ and the saturation-density fits; cp/quality are
  never read. Subcooled liquid gives x = 0, so the same boiling curve applies with Tf < Tsat
  (h = q/(Tw − Tf)).
* **q-continuity.** q is continuous at T_MFB (q_MFB = h_FB(T_MFB)(T_MFB − Tsat)) and at T_CHF
  (q_CHF = q_NB(ΔT_CHF) for power/cooper/cooper_add). With θ^0.5 the slope dq/dTw is infinite just
  below T_MFB. That is harmless here: no iteration chatter at rewet was seen in any run.
* **T_CHF ordering guard.** T_CHF,eff = min(Tsat + ΔT_CHF, T_MFB − gap). Putting a guard on T_MFB
  instead (max(T_MFB, T_CHF + gap)) copies the whole T_CHF expression (Katto + Cardano) into every
  T_MFB occurrence inside the film-property fits, and the `spec` expression grew to 6 MB.
* **Expression size.** Each superheated-vapour fit is about 20 kB because the language has no
  variables, so every copy of T and P is re-emitted. Variants with film properties are about 0.2 MB;
  `spec` is 0.7 MB. Only the taken ternary branch is evaluated, and the solve cost stayed similar to
  the sat-property variants.
* **Per-conductor axial position.** Use the reserved param `z`. `CustomClosure.correlation()` sets it
  per conductor to the node x; Darr 2016 and Katto–Ohno read it.
* **CHF at the inlet node.** Katto–Ohno needs L > 0, so L = max(z, z_min = 2D).
* **Darr 2016.** Uses x_eq (the HEM node quality) in (1 − x)^−0.76, and L = z.
* **Verification.** expression == numpy to ≤ 1.5e−13 relative for every variant: tiny network, three
  inlet states × four wall temperatures (film, transition and nucleate regimes), plus short UF runs.

### 1-D numerics findings (phase 2)
* **Run time and convergence.** N = 20, dt = 0.05: 25–80 s per case, and every run converged with 0
  stalled steps. The exception is case 6 (125–330 s), which spends more steps in two-phase
  fluid-iteration grinds.
* **Rewet.** No regime chatter. Every conductor crosses h = 1500 W/m²K once on the way up and once on
  the way down (on entering the nucleate tail). Two artefacts were seen:
  - Post-quench tail (Tw − Tsat < 2 K). h = q/(Tw − Tf) is ill-conditioned there and jumps by 10×,
    but q ≈ 0, so the wall is unaffected.
  - Quench-induced dry-out flash (case 6, node 12). The near-simultaneous quench of many nodes boils
    the local HEM node to x → 1 for 1–2 steps, with P rising from 225 to 280 kPa. h drops to vapour
    DB (≈20 W/m²K) and then recovers. This is physical within HEM rather than numerical.
* **No quench-front propagation.** All nodes rewet within 0.2–1.5 s of each other. The model front
  speed is 0.7–2.4 m/s against 66–106 mm/s in the data. Axial conduction in the 0.5 mm 304SS wall is
  negligible in film boiling (|q_axial|/|q_conv| ≈ 0). Rewet is purely local, set by
  T_w = T_MFB, so the 1-D model behaves like the 0-D surrogate along the tube.
* **Resolution (cases 7 and 10, variants a and spec; N = 10/20/40 at dt = 0.05 s, dt = 0.1/0.05/0.025 s
  at N = 20).** Converged in N and dt:
  - t130 varies by ≤ 0.13 s.
  - Station RMS varies by ≤ 0.5 K, except spec, case 7 dn, dt = 0.1, at +1 K.
  - Quench spread along the tube varies by ≤ 0.1 s.
  - Run time is NOT monotonic in dt. Case 7 at dt = 0.1 took 249–263 s, against 27–46 s at dt = 0.05,
    because of outer-iteration grinds on larger steps. **dt = 0.05 s is the sweet spot.**
* **Sensitivity of mathematically identical expressions.** Two forms of the same closure (variant a
  before and after a refactor of the expression layout) differ by ≤ 2.5 K (RMS 0.5–1 K), and only in
  the post-quench tail at Tw ≈ Tsat. There, h = q/(Tw − Tf) is ill-conditioned and the outer-iteration
  path depends on rounding. t130 and quench timing are identical.
* **Agreement with the 0-D surrogate.** The 1-D model matches the analysis 0-D surrogate:
  - Variant a: mean station RMS 12.1 K (1-D) against 12.4 K (0-D Pe-law, 134 K); t130 MAE 1.6 s
    against 1.7 s.
  - Model Tsat at the stations matches the measured Tsat (median −0.03 K), because the inlet and exit
    pressures are imposed from P1/P2.
  - Station quality stays small (median x ≈ 0.025 upstream, 0.055 downstream), and below x_dry the
    closure does not depend on x.
  - Axial conduction is negligible.
  - The 1-D model therefore adds nothing over 0-D except that it computes x and P itself. Both lack a
    propagating quench front.

## 8. Phase 4 (independent validation) — solver/runner facts

* **CoolProp-wasm abort poisons the Node module (v0.3.1).** Subcooled inlets (NBS figures 13, 14 and
  3–7) make CoolProp throw during compressed-liquid P–h flashes. Each throw leaks a C++ exception
  through WASM, and after enough of them the module fails. The symptoms are "memory access out of
  bounds", "WebAssembly.Table.get(): invalid address", or a valid boundary flash failing
  (`applyBoundaryConditions ... f0 at t=26 (P=425565, T=76)`).
  - This happens 2–13 minutes into the solve.
  - It is fixed upstream after 0.3.1 (OpenFLUME-uf 436fdc5, compressed-liquid floor in clampToValidPH).
  - Once poisoned, **every later solve in the same session fails instantly**. `runner._worker` now
    kills and restarts the session after any core or bridge error. Before that change, one crash
    cascaded through all later jobs of the worker.
  - Workaround used for the subcooled runs: a saturated inlet at P_drive (`force_saturated`,
    labelled "sat. inlet (fallback)"). This changes the inlet enthalpy, so it is a deviation from the
    experiment.
* **Darr 2016 film law on long or vapour-filled lines.** The (1 − x)^−0.76 factor (×33 at the x = 0.99
  clip) makes h jump by about 30× when liquid first reaches nodes that are near-dry. The outer
  iteration grinds: NBS fig 10 d119s reached only t ≈ 18 s in 90 s of wall time.
* **Runner.** Jobs record `last_sim_t` (from the progress callback), so a timeout says how far the
  solve got. NBS b/d119s runs are capped at 600 s, other NBS runs at 1500 s, and UF-tube
  Miropolskii/subcooled runs at 1200 s.
* **NBS time step.** dt = 2 s (LN2) and 0.5 s (LH2), against OpenFLUME's NBS validation values of
  5 s and 1.25 s. With dt = 1 s and 0.25 s a run took 5–10 min and the batch was too slow. The dt = 1 s
  results for figures 10–12 are kept in tag `p4` as a time-step reference.
* **NBS dt check (spec, fig 10 / fig 12).** dt = 1 s and dt = 2 s agree: station knee times within
  1.4 s, RMS within 1 K.
* **Choice of gravity term** (implemented as an option before any 1-g run). The Bromley natural-convection
  film term is cube-summed with the forced film law: h = (h_f³ + h_b³)^(1/3).
  - Horizontal line (NBS): horizontal-tube form, 0.62 [k³ρ_v(ρ_l−ρ_v) g h_fg'/(D μ ΔT)]^¼, with h_fg'
    = h_fg + 0.4 cp ΔT.
  - Vertical upflow (Darr IJHMT): vertical-wall form, 0.943 [...]^¼, with length L = max(z, D) from the
    tube inlet. The vapour film rises along the wall, so the length is along the flow rather than D.
  - Effect on the Darr 1-g tests: mean RMS 17.3 → 7.3 K; quench timing error 3.1 → 0.7 s.
  - Effect on NBS (horizontal, large G): small, ≤ 2 K RMS. On F2 (0.05 g, horizontal form): t130
    error 0.50 → 0.18 s.
* **Darr IJHMT boundary conditions.** Only mean ṁ and P_out are known, plus T_in.
  - ṁ is applied as a constant flowSource from t = 0 (a step start).
  - The inlet boundary is at max(Psat(88.44 K), P_out + 5 kPa); it only sets the inflow state.
  - Hydrostatic head (0.57 m) is ignored because ṁ is imposed.
  - Radiation from the 293 K chamber (about 0.2 kW/m², 2 % of the film-boiling flux) is not modelled.
* **F2 tube length.** 0.676 m is used: Fig. 5 shows 26.6 in inside the vacuum, and the "exit" TC
  station at 24 in = 0.610 m must lie inside the tube. Sensitivity: L = 0.61 m (station at the exit
  node) and L = 0.914 m (text value) change RMS by ≤ 0.05 K and t130 by ≤ 0.03 s. The station's thermal
  history is local in low g.

## 9. Phase 5a (production) — facts
* **Miropolskii on the UF tube with G(t + 2.5 s) is unreliable.** N = 20, dt = 0.05 had 4 failures in
  10 cases: timeouts at t = 7–12 s in cases 1, 6 and 9, and a CoolProp-wasm crash in case 5. N = 10,
  dt = 0.1 had 1 failure (case 6, crash at t = 26 s). Tables report the N = 20 run where it finished
  and the N = 10 run otherwise (`make_paper_figs.uf_run`); case 6 has no Miropolskii result.
* **h-continuity scan of spec10** (`results/final/h_continuity.csv`): 81 combinations of P, G, z and x.
  The maximum |Δq| between neighbouring scan points falls by about √10 when the scan step falls by
  10×. That is the θ^½ transition onset at T_MFB, which has infinite slope but no finite jump.
* **NBS inlet K.** K = 2 is OpenFLUME's NBS validation value (`NBS_INLET_RESTRICTION_K`, with
  sensitivity 1 and 5; OpenFLUME-uf `src/validation/nbsChilldown.ts`). K = 10 is about a fully open
  globe valve plus entrance (Crane TP-410). Over K = 1–10, mean RMS changes by < 0.5 K: the 61 m line's
  two-phase friction dominates the flow split. The documented 83 kPa outlet pressure lowers RMS by
  about 2 K.

## 10. Phase 6 — NBS 9264 boundary conditions and long-line physics (ofchill/nbs6.py)
* **The NBS report** is now in `lit/raw/nbs9264.pdf`/`.txt` (NTRS 19670007291). Facts used:
  - Inlet valve: a ¾-in ball valve inside the dewar, or a 1-in globe/gate valve outside it. Valve
    opening was varied by the actuator pressure.
  - Flow was measured by a turbine meter in the dewar.
  - Discharge to atmosphere at 0.82 atm (Boulder), so **P_out = 83 kPa is the documented value**.
  - Wall T < 70 K was not reliably recorded (§3.3). LH2 RMS is therefore also reported on data ≥ 70 K
    (RMS70), and LH2 timing on the 100 K crossing and the knee.
  - Fig. 41 (subcooled LN2, 5.9 atm): initial inflow ≈ 4 L/s, then reversals with surges to 21 atm.
    The steady-flow rate is not given; "experimental uncertainties in determining the mass flow rate"
    are noted.
  - Figs. 45–51: experimental film h of 200–600 W/m²K for LN2.
* **The inlet loss K is not identifiable for LN2.**
  - Fig 10, K = 0.01 … 100: mean RMS 47.7 → 48.7 K; mean ṁ 0.045 → 0.041 kg/s.
  - The flow is limited by friction of the warm vapour (and later two-phase) flow along the 61 m line:
    about 150 kPa of line ΔP at ṁ ≈ 0.02 kg/s early on.
  - Doubling or halving the film-law C changes RMS by < 0.7 K.
  - A static initial line at P_out with a 1 s valve-opening ramp gives the same result to 0.02 K.
  - So the LN2 error is not driven by K, P_out or the initial condition.
  - LH2 does respond to K (fig 2: RMS70 41 → 29 K at K = 100). K = 100 would need about 300 kPa of
    inlet loss at the LH2 flow rate, which is implausible for an open ball valve. K = 2 keeps knee
    errors within 1–4 s.
* **Closure sensitivity on NBS.**
  - Miropolskii's low two-phase h gives the lowest LN2 RMS (25 K) and good 150 K arrivals, but it never
    quenches station 1 (knee +117 to +214 s).
  - The present closure quenches on time upstream. Its station-4 chilldown times are +9 to +31 s
    (+7 to +13 %), comparable to GFSSP's published +7 to +12 %.
  - Its 150 K arrivals at stations 2–3 are 17–23 s late in the low-pressure LN2 runs. The data show
    slow precooling at 24 m from t ≈ 20 s, while station 1 is still at 220 K. That needs about 3× the
    model's early liquid supply, or liquid running ahead of the front: 1-g stratification and thermal
    non-equilibrium, which a 1-D homogeneous-equilibrium model cannot carry.
* **Tested and rejected (no fitted constant).**
  - `strat`: wetted-fraction 1-g stratification, from the Zivi void fraction, Biberg wetted angle and
    Froude blend.
  - `no_dbx`: drop the vapour-DB(Gx) floor in two-phase nodes.
  - Both lower LN2 RMS (fig 10: 45 → 36 K with both) but delay the station-1 quench and worsen every
    knee (fig 10 knee MAE 9 → 44 s).
  - They are kept as options and are not used.
* **Rewet R5 (CLOSURE_SPEC §11)** is implemented as `tmfb='R5'` (param `e_hn` = 0.122). Variant `spec11`
  is the §10 film law plus R5; the gravity term (`_gh` horizontal, `_gv` vertical) is used for g > 0.
* **Subcooled-inlet crash workaround that works:** `runner.solve_segmented`. The run is split into
  15 s segments, each solved in a fresh Node session and restarted from the previous end state
  (internal nodes take `quality` when two-phase, otherwise T; solid T; schedules shifted). A failed
  segment is retried at ¼ length.
  - Fresh sessions reset CoolProp-wasm's leaked-exception count.
  - NBS fig 3 (LH2 subcooled, 150 s) ran in 10 segments with 5 stalled steps.
  - Caveats: branch flows restart from the solver's initial guess at each restart (a small transient),
    and the per-solve timeout is not enforced inside segmented runs.
  - Inlet quality > 0 or a saturated inlet at P_drive also avoid the crash, but change the inlet
    enthalpy (the saturated fallback gives LN2 RMS ≈ 75 K).
* **Phase-6 final results** (`results/nbs6/aggregate.md`, `final_runs.csv`, `figs/`). Settings: P_out =
  83 kPa (documented), K = 2, N = 30, `spec11_gh`, horizons 1.6× the record.
  - **LH2, all 6 runs**, the 5 subcooled ones via segmented solves: station-4 chilldown time within
    −6 … +5.5 % (MAPE 5.0 %, bias −0.1 %); knee MAE 2.4 s; 100 K crossing MAE 2.4 s.
  - **LN2 saturated, figs 10–12**: chilldown time +8 … +13 % (MAPE 10 %); knee MAE 18 s; 150 K arrival
    MAE 16 s; RMS 42 K.
  - Miropolskii on the same LN2 runs: RMS 25 K, but its knees are 57 s late on average (station 1 never
    quenches) and chilldown time is +12 %.
  - **LN2 subcooled, figs 13–14**: did not complete. The segmented solves do not crash, but one segment
    ground for more than 35 min in the outer iteration; they were killed.
* **Inverse K.** LN2 RMS is minimised at K → 0, and the optimum is flat (ΔRMS < 1 K for K = 0–2, rising
  with K: fig 12 26.9 → 33.1 K at K = 100). LH2 fig 2 prefers K ≈ 100 (RMS70 44 → 32 K, knee 2.3 → 1.3 s),
  but that needs a 130 kPa inlet loss at 0.08 kg/s. That is implausible for the ¾-in ball valve and
  in-dewar turbine meter. The subcooled LH2 runs fit at K = 2. One K = 2 is consistent across both
  fluids.
* **UF / F2 / 1-g with `spec11`**:
  - UF mean RMS 11.9 K (spec10: 11.2); t130 MAE 1.17 s.
  - F2: RMS 9.0 K, t130 error +0.46 s (spec10: 11.0 K, +1.08 s).
  - Darr 1-g with `_gv`: RMS 8.1 K, quench MAE 0.83 s (spec10_gv: 8.7 K, 0.99 s).

## 11. Phase 7 — is two-phase friction the lever for NBS LN2? (results/nbs7)
* **OpenFLUME 0.3.1 pipe friction** (`src/core/components/pipe.ts`, `solver/kernel.ts`) is pure HEM:
  - ΔP_f = f(Re) (L/D) ṁ|ṁ| / (2 ρ A²), with Churchill / Swamee–Jain f.
  - ρ is the upstream node's HEM mixture density. With kineticEnergy on, it is the harmonic mean of the
    upstream and downstream densities.
  - μ is the McAdams mixture viscosity, 1/(x/μ_g + (1−x)/μ_f) (realFluid.domePH). Re uses that μ.
  - There is **no two-phase multiplier option** (no Friedel, Lockhart–Martinelli or MSH).
  - The only knobs are `roughness` and a fixed `frictionFactor`. A `userComponent` (JS) could implement
    Friedel (args give ρ, μ, P, T and a fluid accessor), but that was not needed, see below.
  - Acceleration is included through momentumFlux.
  - Roughness is 1.5 µm (drawn copper/brass tubing, standard tables). At ε/D ≈ 1e−4 and Re ≈ 1e5 the
    line is essentially hydraulically smooth (< 2 % in f).
* **Line ΔP budget, offline at mid-run (fig 11, spec11, t = 20–150 s)**:
  - Total line ΔP ≈ 250 kPa.
  - The two-phase cells carry only 1–21 kPa (HEM). Friedel gives 1.4×, MSH 1.6× that.
  - The rest is friction of the **warm superheated vapour** being expelled ahead of the front: smooth-pipe
    single-phase gas friction, which is not uncertain at the ±50 % level.
  - A Friedel or MSH multiplier would therefore *reduce* the flow, by ≤ 5 %. Two-phase friction is not
    the lever.
* **Direct test of the "more early flow" hypothesis**: fixed f = 0.017 (≈ baseline), 0.0085 and 0.0043.
  - More flow moves every station earlier together. On fig 10 at f/2, the station-2 onset error goes
    +68 → +38 s while station 4 goes +34 → −32 s, and RMS gets worse (43 → 49 K).
  - The data's spatial signature cannot be produced by any flow level: early, *gradual* cooling at 24 m
    (≈ 3× smaller slope than any model over 50–150 s) together with a sharp, on-time knee at 60 m.
  - Lower h ahead of the front (Miropolskii, no-DB(Gx), stratified wetted fraction) gives intermediate
    slopes but delays upstream quench.
* **Verdict.** The LN2 misfit is the precursory, gradual cooling far downstream of the main front. That is
  liquid carried ahead in a separated, non-equilibrium form (inrush slug or stratified rivulet in film
  boiling, droplet-laden vapour) while the main front stays on time. A homogeneous-equilibrium 1-D model
  cannot represent it, whatever the friction or the wall h.
  - NBS describes the first-inflow slug (Fig. 41: ≈ 4 L/s for ~1.5 s, then reversals) and links it to
    the surges.
  - The NBS chilldown times are still reproduced (+8 … +13 %).
* **Other closures on NBS LN2 (figs 10–12).**
  - z-independent spec (`spec_gh`): RMS 43.9 K, knee MAE 19.1 s, onset MAE 34.7 s. That is the same as
    spec11 (41.7 / 18.3 / 33.7): the film law has little effect, because the line is limited by liquid
    supply.
  - Miropolskii: lowest RMS (25.3 K) and onset MAE (19.9 s), but knees are 57 s late because station 1
    never quenches.
  - Darr ×1.19 (`d119s_gh`) grinds in the outer iteration on the long line, as in Phase 4: no result
    after more than 45 min per run. It was left running into `results/nbs7`.
  - Bug fixed: before Phase 7, `closure_for('d119s_gh')` silently built the spec closure. Phase-4 runs
    are unaffected; they used plain `d119s`.
* **No change was adopted.** Friction options don't help. The heat-transfer options (no-DB(Gx),
  stratified wetted fraction, Miropolskii-level h) trade onset error for knee error. Because the closure
  is unchanged, LH2, UF, F2 and 1-g need no rerun.

## 12. Phase 8 — consistency re-runs with spec11, and chilldown-completion metrics
* `ofchill/final11.py` runs everything into `results/final11/{uf,sens,num,val}`, in the same table formats
  as Phase 5a.
  - New closure variants:
    - `specz11`: z-independent Pe^(1/3) law (C 8.55) + R5;
    - `ds11` / `d119s11`: Darr 2016 film law (published / ×1.19) + R5, so the film law is the only
      difference from spec11.
  - Old-rewet rows (Darr published and ×1.19 with the contact-T_c rewet, DB-vapour, Miropolskii) are
    merged from `results/final/uf`.
  - The T_MFB ±5 K sensitivity uses Δe_hn = ±5/((1+β)T_c), which shifts the R5 T_MFB by exactly ±5 K.
* `ofchill/chilldown.py` computes chilldown completion (wall below T_sat + 10 K) and quench (maximum
  cooling rate) for every closure and station, with hit/miss/false outcomes and the whole-section
  (last-station) time. Output: `results/final11/chilldown_time.md`, `chilldown_stations.csv`,
  `chilldown_summary.csv`; figure F9.
* Miropolskii (N = 10, dt = 0.1) timed out on Darr test 31 (1800 s), so that row has no result.
