# ofchill — OpenFLUME pipe-chilldown simulation / validation pipeline

Runs OpenFLUME (pip `openflume` 0.3.1) on the UF 2016 low-g LN2 chilldown
cases (Darr et al. 2016, `data/uf2016/`) and compares with the measured
outer-wall temperatures.  Always use `.venv/bin/python` from `chilldown-lowg/`.
Solver facts (h evaluation/relaxation, expression scope, quirks): `NOTES.md`.

## CLI

```bash
.venv/bin/python -m ofchill.run --cases 7 8 9 10 --closure miropolskii --N 10 --dt 0.1 --tag baseline
.venv/bin/python -m ofchill.run --cases 7 8 9 10 --closure dittusBoelter db_vapor --tag baseline
.venv/bin/python -m ofchill.run --cases 8 --closure boil_demo --param width=10 --param dT_rewet=60 --tag probe
.venv/bin/python -m ofchill.run --tag baseline --report-only        # rebuild tables/figures
# options: --g-shift S  --inlet-quality X  --exit-p-scale F  --end-time T  --no-inertia
#          --workers (<=6)  --timeout SEC  --force  --retry-failed
.venv/bin/python -m ofchill.repro_prior                              # reproduce prior TS runs
.venv/bin/python -m ofchill.tests.test_closure_expr                  # tests
```

Long batches: `nohup .venv/bin/python -m ofchill.run ... > results/x.log 2>&1 &`.

Outputs: `results/_cache/<hash>.{npz,json}` (one per distinct NetworkConfig;
re-runs are skipped) and `results/<tag>/`: `index.json` (run id → hash, spec),
`metrics.md/.csv`, `runs.md/.csv` (status, stalled steps, solve time,
energy/mass diagnostics), `figs/case<N>.png`, `figs/summary.png`.

npz arrays: `t, wall_x, wall_T[node,t], T_up, T_dn` (station = linear
interpolation between wall nodes), `h, q, Q[conductor,t]` (h is the fresh
end-of-step value, see NOTES §1), `fluid_P, fluid_T, fluid_x (NaN outside the
dome), fluid_rho, phase_code (0 liq/1 two-phase/2 vap)`, `mdot_in, mdot_out,
step_res_scaled`.  json: status (`ok|failed|timeout`), converged, stalled steps
(scaled residual > 0.01), energy diagnostic, solve time, error text.

## Python API

```python
from ofchill.builder import build_uf_case, build_pipe_chilldown
from ofchill.closures import get_closure, db_vapor
from ofchill.closure import closure, NamedClosure
from ofchill import expr as E
from ofchill.props import fluid_props, verify_fits
from ofchill.runner import make_job, run_jobs, load_run
from ofchill.metrics import compare, data_metrics

m = build_uf_case(8, N=10, dt=0.1, closure=get_closure("miropolskii"),
                  g_shift_s=0.0, inlet_quality=0.0, exit_p_scale=1.0, end_time=None, inertia=True)
m.config            # NetworkConfig dict (openflume.Network.from_dict(m.config).solve())

# generic builder (NBS long copper line, LH2, ...)
m = build_pipe_chilldown(L=..., D=..., OD=..., wall_rho=8940, wall_material="ofhc-copper",
                         fluid="Hydrogen", closure=..., wall_T0=lambda x: 290.0,
                         p_in_schedule=[[t, Pa], ...], p_out_schedule=[...],
                         mdot_schedule=[[t, kg/s], ...]  # or inlet_component={...}
                         , N=20, dt=0.05, end_time=60.0)

# custom closure: one Python function -> numpy evaluation AND OpenFLUME expression
@closure(params={"C": 0.023}, fluid="Nitrogen")
def my_h(s):                      # s.Tw s.Tf s.P s.G s.D s.rho s.mu s.k (s.cp/s.quality: see NOTES §3)
    Tsat = s.fp.Tsat(s.P)
    Tfilm = 0.5 * (s.Tw + Tsat)
    kv, muv, cpv = (s.fp.vap(n, Tfilm, s.P) for n in ("k", "mu", "cp"))
    h_film = s.param("C") * (s.G * s.D / muv) ** 0.8 * (cpv * muv / kv) ** 0.4 * kv / s.D
    return E.where(s.in_dome(), h_film, 50.0)

my_h.expression                      # string for correlation.expression
my_h.evaluate(Tw=..., P=..., G=..., D=0.01168, rho=..., mu=..., k=...)   # numpy
build_uf_case(8, closure=my_h.with_params(C=0.03))

jobs = [make_job("c8_myh", m, {"case": 8, "label": "myh", "N": 10, "dt": 0.1})]
run_jobs(jobs, tag="mytag", max_workers=6)          # then: python -m ofchill.run --tag mytag --report-only
```

Closures: named `miropolskii`, `dittusBoelter`, `darrHartwig`, `ttWf`
(`NamedClosure(model, extra={...})`; axialPosition/segmentLength filled by the
builder); custom `db_vapor` (DB with sat-vapour props at full G in the dome =
the prior patched dispatcher), `miropolskii_expr` (expression copy of the
named model), `boil_demo` (toy step/smoothstep boiling curve), `probe_all` (test).

Property fits (`props.py`, cached in `_fits/`): saturation Tsat, ρl, ρv, h_fg,
σ, μl, μv, kl, kv, cpl, cpv vs P and superheated-vapour k, μ, cp, ρ vs (T, P)
for Nitrogen (80–1500 kPa), Oxygen (50–2000 kPa), Hydrogen (50–1100 kPa),
T ≤ 320 K.  `verify_fits(fluid)` checks against CoolProp on an independent
random grid: max error sat ≤ 7e-4, vapour ≤ 1.8e-3 (all fluids).

Metric definitions: docstring of `metrics.py` (RMS/bias over the record,
t at 200/150/130/100 K, slope 0–2 s, time/T of max cooling rate, T_end).
UF case 4 is truncated at 18.7 s (pressure collapse at the end of the record).

## Phase 2: unified boiling-curve closure

```bash
.venv/bin/python -m ofchill.fit_fb                                   # a-priori film-boiling constants -> fb_constants.json
.venv/bin/python -m ofchill.phase2 --variants a b c_rl c_rv d spec --tag p2       # all 10 cases, N=20 dt=0.05, G(t+2.5 s)
.venv/bin/python -m ofchill.phase2 --variants spec --parasitic --tag p2           # + CLOSURE_SPEC §9 rig heat leak
.venv/bin/python -m ofchill.phase2 --variants a --cases 7 10 --N 40 --dt 0.025 --tag p2num
.venv/bin/python -m ofchill.run --tag p2 --report-only               # metrics.md, aggregate.md, figs/
.venv/bin/python -m ofchill.diag_p2 --tag p2 --label va              # 1-D vs 0-D surrogate diagnostics
.venv/bin/python -m ofchill.chatter p2                               # regime-switch / re-dry counts
```

`make_boiling_closure(fb, props, tmfb, nb, chf, gravity, vapor_film, **params)` builds one
q-continuous boiling curve (film → θ^m transition → nucleate, with DB liquid/vapour limits) as a
single custom expression. Structural options select the expression shape and `params` hold the
numbers. See the docstring and NOTES §7.

| variant | film boiling | T_MFB | NB / CHF |
|---|---|---|---|
| a | Pe_l law, C 7.78, n 0.345, film props | 134 K | power, ΔT_CHF 24 K, q_CHF 100 kW/m² |
| a_sat | Pe_l law, C 12.93, n 0.293, sat props | 134 K | power |
| b | as a | T_c + β(T_c − T_f), β = e_l/e_w per case (≈0.08–0.09) | power |
| c_rl / c_rv | renewal: C 2.35 √(k ρ cp u_l/D) / C 0.274 √(k cp G/D), sat props | 134 K | power |
| d / d119 | Darr 2016 low-g (published) / ×1.19 | 134 K | power |
| a_cooper, a_thn, e_brom | a with Cooper NB / T_hn+18.8 K / Bromley C 1.89 | | |
| spec | CLOSURE_SPEC: Pe^(1/3), C 8.55, film props, floor Nu 4.36 | T_c + β(T_c − min(T_f, T_sat)), β = e_l/e_w(T_MFB) | Cooper + DB_l (additive), Katto–Ohno CHF |

The C values are least-squares fits in ln h to the a-priori film-boiling dataset (`fit_fb.py`,
energy-consistent h, G(t + 2.5 s)).

Results: `results/p2/` (all variants, N = 20, dt = 0.05; `aggregate.md`, `metrics.md`, `figs/`,
`figs_key/`: a, b, c_rv, d, spec, spec+parasitic; `diag_va.md` / `diag_vspec.md`: 1-D vs 0-D), `results/p2num/` (N/dt study, cases
7 & 10), `results/p2lag0/` (G without the 2.5 s lag).  Variant `spec_par` = spec + `--parasitic`.

## Phase 4: independent validation (no refitting)

```bash
.venv/bin/python -m ofchill.validate4 run --set f2 darr nbs_ln2 nbs_lh2 nbs_sens spec10 --tag p4
.venv/bin/python -m ofchill.validate4 report --tag p4      # metrics.md, aggregate.md, runs.md, figs/
```

Builders in `validate4.py` all use `build_pipe_chilldown`:
- `build_f2` (Hartwig 2022 F2; L = 0.676 m)
- `build_darr` (Darr 2016 IJHMT tests 31 and 40)
- `build_nbs` (NBS 9264 figures 2–7 and 10–14; K = 2 inlet `resistance`, OFHC copper, N = 30)

`closure_for(name, fluid, wall, ...)` handles:
- `spec`, and `spec_gh` / `spec_gv` (`spec` with the Bromley gravity term, horizontal or vertical form);
- `spec10` (CLOSURE_SPEC §10 position-dependent film law, `fb='pe_lz'`), with the same `_gh` / `_gv` suffixes;
- `b` (LN2 only);
- `d119s` (Darr 2016 ×1.19 film law with the spec rewet/NB/CHF);
- `miropolskii`.

β = e_l/e_w(T_MFB) is computed per dataset for the wall material (`contact_beta_spec`, SS304 or copper).

Phase-4 results:
- `results/p4v/` (final): `metrics.md`, `aggregate.md`, `runs.md` (status / `last_sim_t` / errors), `figs/`
- `results/p4/`: NBS runs at dt = 1 s (time-step reference)
- `results/p2/`: the UF rerun with `vspec10`; `figs_spec10/` compares spec, spec10 and b

## Phase 5a: production runs and paper figures

```bash
.venv/bin/python -m ofchill.final run --part uf num nbs sens    # -> results/final/<part>/ (cached)
.venv/bin/python -m ofchill.final report                        # -> results/final/*_aggregate.{csv,md}, *_per_station.csv
.venv/bin/python analysis/make_paper_figs.py [F1 ... F8]         # -> paper/figures/*.pdf|png
```

Tables in `results/final/`:
- `uf_aggregate.md` and `uf_per_station.md`/`.csv`: the UF closure comparison.
- `sens_aggregate.md`: spec10 sensitivities.
- `num_aggregate.md`: resolution study.
- `nbs_aggregate.md`: NBS boundary-condition study.

Figures:

| figure | content |
|---|---|
| F1 | film-law parity |
| F2 | rewet temperature |
| F3 | UF histories |
| F4 | quench front |
| F5 | F2 and Darr 1-g validation |
| F6 | NBS fig 10 |
| F7 | sensitivity tornado |
| F8 | closure boiling curve against the energy-consistent q |

## Phase 6: NBS 9264 boundary conditions and long lines

```bash
.venv/bin/python -m ofchill.nbs6 run --batch probe ic strat final_sat kfit final_sub   # results/nbs6
.venv/bin/python -m ofchill.nbs6 report                                               # per_run.md / per_run.csv
.venv/bin/python -m ofchill.nbs6_plot OUT.png FIG "nbs6|<run id>" ...                 # station plots
```

Run ids have the form `fig:closure:k=v,...`. Available keys:

| key | meaning |
|---|---|
| `K` | inlet loss coefficient |
| `Pout` | outlet pressure |
| `end` | end time |
| `seg` | segmented solve, the subcooled-crash workaround (`runner.solve_segmented`) |
| `static` / `ramp` | static initial line / valve-opening ramp |
| `x_in`, `Tin` | inlet quality / inlet temperature |
| `cp_<param>` | closure parameter override |

New closure variants:
- `spec11`: §10 film law plus the §11 rewet rule R5 (`tmfb='R5'`). This is the final default.
- `spec11s` / `spec11n` / `spec11sn`: tested 1-g stratification and no-vapour-DB(Gx) options; not adopted.

Phase 7 (`NBS_TAG=nbs7 .venv/bin/python -m ofchill.nbs6 run --batch p7_cmp p7_fric`; results in
`results/nbs7/aggregate.md` and `figs/`) tests friction as the lever on NBS LN2. New run keys: `f` (fixed
Darcy f for all pipe segments) and `rough` (roughness). The NBS per-run metrics now include `on20_*`, the
first-cooling (T₀ − 20 K) onset errors.
