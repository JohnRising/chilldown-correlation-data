"""Phase 4: validation of the boiling closures WITHOUT refitting on independent datasets.

  .venv/bin/python -m ofchill.validate4 run --set f2 darr nbs_ln2 nbs_lh2 [--workers 6] [--tag p4]
  .venv/bin/python -m ofchill.validate4 report [--tag p4]

Datasets (data/<name>/README.md):
  f2       Hartwig 2022 F2: low-g (0.05 g) LN2, UF-type 304SS tube, exit-station wall T (x = 0.610 m),
           P_exit(t), G(t) (gas meter).  Tube length: 0.676 m (Fig. 5: 26.6 in inside the vacuum; the TC
           "exit" station at 24 in = 0.610 m must lie inside the tube); sensitivity 0.61 / 0.914 m.
  darr     Darr 2016 IJHMT tests 31 / 40: 1 g vertical upflow LN2, UF tube, mean mdot / P_out / T_in only.
  nbs_*    NBS 9264: 61 m (200 ft) horizontal OFHC-copper line, 1 g, pressure driven: dewar at P_drive
           (saturated: quality 0; subcooled: T_in), K = 2 inlet resistance (OpenFLUME
           NBS_LINE_FLOW_MODEL: internal end layout, longitudinal inertia, OFHC copper), outlet 101325 Pa,
           N = 30.  Stations 6.1 / 24.4 / 43.0 / 60.35 m.
"""

from __future__ import annotations

import argparse
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .builder import build_pipe_chilldown
from .closure import NamedClosure
from .closures import boiling_variant, contact_beta_spec
from .data import DATA_DIR
from .metrics import cross
from .runner import RESULTS, load_index, load_run, make_job, run_jobs

UF_D, UF_OD, UF_RHO = 0.01168, 0.0127, 8027.0
NBS_D, NBS_OD, NBS_RHO, NBS_L = 0.015875, 0.01905, 8960.0, 60.96
NBS_STATIONS = {"s1": 6.1, "s2": 24.4, "s3": 43.0, "s4": 60.35}


def _csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(io.StringIO("".join(l for l in path.open() if not l.startswith("#"))))


def tsat(P, fluid="Nitrogen"):
    import CoolProp.CoolProp as CP

    return CP.PropsSI("T", "P", float(P), "Q", 0, fluid)


def psat(T, fluid="Nitrogen"):
    import CoolProp.CoolProp as CP

    return CP.PropsSI("P", "T", float(T), "Q", 0, fluid)


# ---------------------------------------------------------------------------- data
def f2_data() -> pd.DataFrame:
    d = _csv(DATA_DIR / "hartwig2022_F2" / "F2_bare_continuous_0p05g.csv")
    d["t"] = d.t_s - d.t_s.iloc[0]  # model t = 0 at the first digitized sample (0.1 s)
    return d


def darr_data(test: int) -> pd.DataFrame:
    d = _csv(DATA_DIR / "darr2016ijhmt_1g_vertical" / f"test{test}_LN2_1g_upflow.csv")
    d["t"] = d.t_s
    return d


DARR_BC = {31: dict(T0=289.9, T_in=86.96, P_out=314.6e3, mdot=0.0334),
           40: dict(T0=288.54, T_in=88.61, P_out=325.9e3, mdot=0.0346)}
DARR_TSAT_IN = 88.44


def nbs_meta() -> pd.DataFrame:
    return _csv(DATA_DIR / "nbs9264" / "nbs9264_runs_metadata.csv")


def nbs_data(fig: int) -> pd.DataFrame:
    f = next((DATA_DIR / "nbs9264").glob(f"nbs9264_fig{fig:02d}_*.csv"))
    return _csv(f)


# ---------------------------------------------------------------------------- closures
def closure_for(name: str, *, fluid: str, wall: str, P_ref: float, D: float, g: float = 0.0, cparams: dict | None = None):
    """name: spec | spec_g (gravity, orientation by dataset) | b | miropolskii | d119s."""
    if name == "miropolskii":
        return NamedClosure("miropolskii")
    if name == "db_vapor":
        from .closures import get_closure

        return get_closure("db_vapor")
    beta, _ = contact_beta_spec(P_ref, fluid, wall)
    kw = dict(fluid=fluid, beta=beta, **(cparams or {}))
    if name.startswith("spec") or name == "d119s":
        kw["z_min"] = 2 * D
    if name in ("ds11", "d119s11"):
        return boiling_variant(name, z_min=2 * D, **kw)
    if name.startswith("d119s"):
        if name.endswith("_gh") or name.endswith("_gv"):
            c = boiling_variant("d119s", gravity="vertical" if name.endswith("_gv") else "horizontal", g=g, **kw)
            return c.with_label(f"v{name}")
        return boiling_variant("d119s", **kw)
    base = name.split("_")[0] if name.split("_")[0] in ("spec", "spec10", "spec11", "spec11s", "spec11sn", "spec11n",
                                                         "spec11d", "spec11dn", "specz11") else "spec"
    if name == base:
        return boiling_variant(base, **kw)
    if name.endswith("_gv") or name.endswith("_gh"):
        orient = "vertical" if name.endswith("_gv") else "horizontal"
        c = boiling_variant(base, gravity=orient, g=g, **kw)
        return c.with_label(f"v{name}")
    if name == "b":
        if fluid != "Nitrogen":
            raise ValueError("variant b has LN2-specific NB/CHF constants")
        return boiling_variant("b", beta=beta)
    raise ValueError(name)


# ---------------------------------------------------------------------------- builders
def build_f2(closure_name: str, L: float = 0.676, g_lag: float = 0.0, inlet_T: float | None = None,
             dt: float = 0.05, N: int | None = None, g_level: float = 0.05 * 9.81):
    d = f2_data()
    t = d.t.to_numpy()
    tq = t + g_lag
    G = np.interp(tq, t, d.G_kgm2s.to_numpy())
    P = d.P_exit_kPa.to_numpy() * 1e3
    P[0] = P[1]  # first digitized P sample (287 kPa) is a start artefact
    A = math.pi / 4 * UF_D**2
    N = N or max(10, int(round(L / 0.0286)))  # ~ UF N=20 spacing
    sched = lambda y: [[float(a), float(b)] for a, b in zip(t, y)]  # noqa: E731
    P_in = 550e3
    cl = closure_for(closure_name, fluid="Nitrogen", wall="ss304", P_ref=float(np.median(P)), D=UF_D, g=g_level)
    return build_pipe_chilldown(
        L=L, D=UF_D, OD=UF_OD, wall_rho=UF_RHO, wall_material="stainless-steel-304", closure=cl,
        fluid="Nitrogen", wall_T0=288.0, p_in_schedule=[[0.0, P_in], [float(t[-1]), P_in]],
        p_out_schedule=sched(P), mdot_schedule=sched(G * A), N=N, dt=dt, end_time=float(t[-1]),
        inlet_quality=None if inlet_T else 0.0, inlet_temperature=inlet_T, stations={"exit": 0.610},
        name=f"Hartwig 2022 F2 L={L}")


def build_darr(test: int, closure_name: str, dt: float = 0.05, N: int = 20):
    bc = DARR_BC[test]
    d = darr_data(test)
    tend = float(d.t.iloc[-1])
    P_in = max(psat(DARR_TSAT_IN), bc["P_out"] + 5e3)
    sub = bc["T_in"] < tsat(P_in) - 0.05
    cl = closure_for(closure_name, fluid="Nitrogen", wall="ss304", P_ref=P_in, D=UF_D, g=9.81)
    return build_pipe_chilldown(
        L=0.572, D=UF_D, OD=UF_OD, wall_rho=UF_RHO, wall_material="stainless-steel-304", closure=cl,
        fluid="Nitrogen", wall_T0=bc["T0"], p_in_schedule=[[0.0, P_in], [tend, P_in]],
        p_out_schedule=[[0.0, bc["P_out"]], [tend, bc["P_out"]]],
        mdot_schedule=[[0.0, bc["mdot"]], [tend, bc["mdot"]]], N=N, dt=dt, end_time=tend,
        inlet_quality=None if sub else 0.0, inlet_temperature=bc["T_in"] if sub else None,
        stations={"up": 0.149, "dn": 0.401}, name=f"Darr 2016 IJHMT test {test}")


def build_nbs(fig: int, closure_name: str, N: int = 30, dt: float | None = None, K: float = 2.0,
              P_out: float = 101325.0, end_time: float | None = None, force_saturated: bool = False,
              cparams: dict | None = None, inlet_quality: float | None = None, inlet_T: float | None = None,
              static_init: bool = False, ramp_s: float = 0.0, friction_factor: float | None = None,
              roughness: float = 1.5e-6):
    m = nbs_meta().set_index("figure").loc[fig]
    fluid = "Nitrogen" if m.fluid == "LN2" else "Hydrogen"
    d = nbs_data(fig)
    tend = float(end_time or d.time_s.max())
    dt = dt or (2.0 if fluid == "Nitrogen" else 0.5)
    P_in = float(m.P_drive_Pa)
    sub = m.liquid_state == "subcooled" and not force_saturated
    # initial wall T: first digitized point per station, linear in x (held beyond the end stations)
    xs, T0 = [], []
    for st, x in NBS_STATIONS.items():
        dd = d[np.isclose(d.station_m, x, atol=0.3)] if not d[np.isclose(d.station_m, x, atol=0.3)].empty \
            else d[np.isclose(d.station_m, 60.4, atol=0.3)]
        xs.append(x)
        T0.append(float(dd.sort_values("time_s").T_wall_K.iloc[0]))
    wall_T0 = lambda x: float(np.interp(x, xs, T0))  # noqa: E731
    A = math.pi / 4 * NBS_D**2
    cl = closure_for(closure_name, fluid=fluid, wall="copper", P_ref=P_in, D=NBS_D, g=9.81, cparams=cparams)
    p_in_sched = [[0.0, P_in], [tend, P_in]]
    if ramp_s > 0:  # valve opening: dewar-side pressure ramps from the line pressure to P_drive
        p_in_sched = [[0.0, P_out], [ramp_s, P_in], [tend, P_in]]
    m_ = build_pipe_chilldown(
        L=NBS_L, D=NBS_D, OD=NBS_OD, wall_rho=NBS_RHO, wall_material="ofhc-copper", closure=cl, fluid=fluid,
        wall_T0=wall_T0, p_in_schedule=p_in_sched, p_out_schedule=[[0.0, P_out], [tend, P_out]],
        inlet_component={"type": "resistance", "k": K, "area": A}, N=N, dt=dt, end_time=tend,
        friction_factor=friction_factor, roughness=roughness,
        inlet_quality=(inlet_quality if inlet_quality is not None else (None if sub else 0.0)) if inlet_T is None else None,
        inlet_temperature=inlet_T if inlet_T is not None else (float(m.T_inlet_K) if sub and inlet_quality is None else None),
        stations=dict(NBS_STATIONS), name=f"NBS 9264 fig {fig}")
    if static_init:  # line initially static at the outlet (ambient) pressure, filled with warm vapour
        for n in m_.config["nodes"]:
            if n["type"] == "internal":
                n["pressure"] = P_out
    return m_


# ---------------------------------------------------------------------------- job sets
F2_VARIANTS = ["spec", "b", "miropolskii", "d119s", "spec_gh"]
DARR_VARIANTS = ["spec", "b", "miropolskii", "d119s", "spec_gv"]
NBS_LN2_VARIANTS = ["spec", "b", "miropolskii", "d119s", "spec_gh"]
NBS_LH2_VARIANTS = ["spec", "miropolskii", "d119s", "spec_gh"]


def jobs_for(sets: list[str], timeout: float) -> list[dict]:
    J = []

    def add(rid, model, spec):
        J.append(make_job(rid, model, spec, timeout=timeout))

    if "f2" in sets:
        for v in F2_VARIANTS:
            add(f"f2_L0.676_{v}", build_f2(v), {"dataset": "f2", "variant": v, "L": 0.676})
        for L in (0.61, 0.914):
            add(f"f2_L{L}_spec", build_f2("spec", L=L), {"dataset": "f2", "variant": "spec", "L": L, "sens": f"L={L}"})
        add("f2_L0.676_spec_lag2.5", build_f2("spec", g_lag=2.5), {"dataset": "f2", "variant": "spec", "sens": "G lag 2.5 s"})
        add("f2_L0.676_spec_Tin80", build_f2("spec", inlet_T=80.0), {"dataset": "f2", "variant": "spec", "sens": "T_in 80 K"})
    if "darr" in sets:
        for test in (31, 40):
            for v in DARR_VARIANTS:
                add(f"darr{test}_{v}", build_darr(test, v), {"dataset": "darr", "test": test, "variant": v})
    for key, figs, vs in (("nbs_ln2", (10, 11, 12, 13, 14), NBS_LN2_VARIANTS),
                          ("nbs_lh2", (2, 3, 4, 5, 6, 7), NBS_LH2_VARIANTS),
                          ("nbs_sat", (10, 11, 12), NBS_LN2_VARIANTS), ("nbs_sat", (2,), NBS_LH2_VARIANTS)):
        if key in sets:
            for fig in figs:
                for v in vs:
                    add(f"nbs{fig:02d}_{v}", build_nbs(fig, v), {"dataset": "nbs", "fig": fig, "variant": v})
    if "spec10" in sets:
        add("f2_L0.676_spec10", build_f2("spec10"), {"dataset": "f2", "variant": "spec10", "L": 0.676})
        add("f2_L0.676_spec10_gh", build_f2("spec10_gh"), {"dataset": "f2", "variant": "spec10_gh", "L": 0.676})
        for test in (31, 40):
            for v in ("spec10", "spec10_gv"):
                add(f"darr{test}_{v}", build_darr(test, v), {"dataset": "darr", "test": test, "variant": v})
        for fig in (10, 11, 12, 2):
            for v in ("spec10", "spec10_gh"):
                add(f"nbs{fig:02d}_{v}", build_nbs(fig, v), {"dataset": "nbs", "fig": fig, "variant": v})
    if "spec11" in sets:
        for v in ("spec11_gh", "spec11n_gh"):
            add(f"f2_L0.676_{v}", build_f2(v), {"dataset": "f2", "variant": v, "L": 0.676})
        for test in (31, 40):
            for v in ("spec11_gv", "spec11n_gv"):
                add(f"darr{test}_{v}", build_darr(test, v), {"dataset": "darr", "test": test, "variant": v})
    if "nbs_satin" in sets:  # v0.3.1 CoolProp-wasm crashes on subcooled inlets -> saturated-inlet fallback
        for fig in (13, 14, 3, 4, 5, 6, 7):
            for v in ("spec", "miropolskii"):
                add(f"nbs{fig:02d}_{v}_satin", build_nbs(fig, v, force_saturated=True),
                    {"dataset": "nbs", "fig": fig, "variant": v, "sens": "sat. inlet (fallback)"})
    if "nbs_sens" in sets:
        for fig in (10, 12):
            for v in ("spec",):
                add(f"nbs{fig:02d}_{v}_K5", build_nbs(fig, v, K=5.0), {"dataset": "nbs", "fig": fig, "variant": v, "sens": "K=5"})
                add(f"nbs{fig:02d}_{v}_Pout83", build_nbs(fig, v, P_out=83082.0),
                    {"dataset": "nbs", "fig": fig, "variant": v, "sens": "P_out 83 kPa"})
    # slow grinders (Miropolskii / subcooled inlet on the short UF tube) last, with a shorter timeout
    slow = [j for j in J if ("miropolskii" in j["run_id"] and not j["run_id"].startswith("nbs")) or "Tin80" in j["run_id"]]
    for j in slow:
        j["timeout"] = min(j["timeout"], 1200.0)
    for j in J:
        if j["run_id"].startswith("nbs"):
            j["timeout"] = min(j["timeout"], 600.0 if j["run_id"].endswith(("_b", "_d119s")) else 1500.0)
    late = [j for j in J if j["run_id"].endswith(("_b", "_d119s")) and j["run_id"].startswith("nbs")]
    return [j for j in J if j not in slow and j not in late] + late + slow


# ---------------------------------------------------------------------------- metrics
def _trace_metrics(tm, Tm, td, Td, thresholds: dict[str, float], onset_drop=None) -> dict:
    m = td <= tm[-1] + 1e-9
    Tmi = np.interp(td[m], tm, Tm)
    e = Tmi - Td[m]
    r = {"rms": float(np.sqrt(np.mean(e**2))), "bias": float(np.mean(e)), "T_end_err": float(Tm[-1] - Td[m][-1])}
    ok = Td[m] >= 70.0  # NBS 9264: wall temperatures < 70 K not reliably recorded (report sec. 3.3)
    r["rms70"] = float(np.sqrt(np.mean(e[ok] ** 2))) if ok.any() else np.nan
    for k, th in thresholds.items():
        a, b = cross(tm, Tm, th), cross(td[m], Td[m], th)
        r[f"{k}_m"], r[f"{k}_d"] = a, b
        r[f"{k}_err"] = None if (a is None or b is None) else a - b
    for drop in (20.0,):
        a = cross(tm, Tm, Tm[0] - drop)
        b = cross(td[m], Td[m], Td[m][0] - drop)
        r[f"on{int(drop)}_err"] = None if (a is None or b is None) else a - b
    if onset_drop is not None:
        a = cross(tm, Tm, Tm[0] - onset_drop)
        b = cross(td[m], Td[m], Td[m][0] - onset_drop)
        r["onset_m"], r["onset_d"] = a, b
        r["onset_err"] = None if (a is None or b is None) else a - b
    return r


def run_metrics(spec: dict, ar: dict) -> list[dict]:
    rows = []
    ds = spec["dataset"]
    tm = ar["t"]
    if ds == "f2":
        d = f2_data()
        Ts = tsat(np.median(d.P_exit_kPa) * 1e3)
        th = {"t130": 130.0, "t_q": max(100.0, Ts + 8)}
        rows.append({"station": "exit", **_trace_metrics(tm, ar["T_exit"], d.t.to_numpy(), d.Twall_exit_avg_K.to_numpy(), th)})
    elif ds == "darr":
        d = darr_data(spec["test"])
        Ts = tsat(DARR_BC[spec["test"]]["P_out"])
        th = {"t135": 135.0, "t126": 126.2, "t97": 97.0, "t_q": max(100.0, Ts + 8)}
        for st, col in (("up", "Twall_TC1-3_avg_K"), ("dn", "Twall_TC4-6_avg_K")):
            rows.append({"station": st, **_trace_metrics(tm, ar[f"T_{st}"], d.t.to_numpy(), d[col].to_numpy(), th)})
    elif ds == "nbs":
        m = nbs_meta().set_index("figure").loc[spec["fig"]]
        fluid = "Nitrogen" if m.fluid == "LN2" else "Hydrogen"
        d = nbs_data(spec["fig"])
        T_knee = tsat(m.P_drive_Pa, fluid) + 15.0  # OpenFLUME traceBaseline data-native knee threshold
        mid = 150.0 if fluid == "Nitrogen" else 100.0
        for st, x in NBS_STATIONS.items():
            dd = d[np.isclose(d.station_m, x, atol=0.3)]
            if dd.empty and st == "s4":
                dd = d[np.isclose(d.station_m, 60.4, atol=0.3)]
            if dd.empty:
                continue
            dd = dd.sort_values("time_s")
            rows.append({"station": st, **_trace_metrics(tm, ar[f"T_{st}"], dd.time_s.to_numpy(), dd.T_wall_K.to_numpy(),
                                                         {f"t{int(mid)}": mid, "t_knee": T_knee}, onset_drop=5.0)})
    return rows


def report(tag: str):
    idx = load_index(tag)
    out = RESULTS / tag
    rows, runs, curves = [], [], {}
    for rid, ent in sorted(idx.items()):
        spec = ent["spec"]
        try:
            meta, ar = load_run(tag, rid)
        except FileNotFoundError:
            continue
        runs.append(dict(run_id=rid, status=meta.get("status"), last_sim_t=meta.get("last_sim_t"), converged=meta.get("converged"),
                         n_stalled=meta.get("n_stalled"), solve_s=round(meta.get("solve_s", 0), 1),
                         error=(meta.get("error") or "")[:150]))
        if ar is None:
            continue
        key = (spec["dataset"], spec.get("test") or spec.get("fig") or "F2")
        if not spec.get("sens"):
            curves.setdefault(key, []).append((spec["variant"], ar))
        for r in run_metrics(spec, ar):
            rows.append(dict(run_id=rid, dataset=spec["dataset"], case=key[1], variant=spec["variant"],
                             sens=spec.get("sens", ""), solve_s=meta.get("solve_s"), **r))
    df = pd.DataFrame(rows)
    df.to_csv(out / "metrics.csv", index=False)
    pd.DataFrame(runs).to_csv(out / "runs.csv", index=False)
    (out / "runs.md").write_text(pd.DataFrame(runs).to_markdown(index=False) + "\n")
    errc = [c for c in df.columns if c.endswith("_err") and c != "T_end_err"]
    agg = df[df.sens == ""].groupby(["dataset", "variant"]).agg(
        n=("rms", "size"), RMS=("rms", "mean"), bias=("bias", "mean"),
        **{f"{c}_MAE": (c, lambda s: np.nanmean(np.abs(s.astype(float)))) for c in errc},
        **{f"{c}_n": (c, lambda s: int(s.notna().sum())) for c in errc}, solve_s=("solve_s", "mean")).reset_index()
    agg.to_csv(out / "aggregate.csv", index=False)
    (out / "aggregate.md").write_text(agg.to_markdown(index=False, floatfmt=".2f") + "\n")
    (out / "metrics.md").write_text(df.drop(columns=["run_id"]).to_markdown(index=False, floatfmt=".2f") + "\n")
    _figures(curves, out / "figs")
    return df, agg


def _figures(curves, figdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    order = ["spec", "spec_gh", "spec_gv", "b", "d119s", "miropolskii", "spec10", "spec10_gh", "spec10_gv"]
    for (ds, case), cur in curves.items():
        cur = sorted(cur, key=lambda c: order.index(c[0]) if c[0] in order else 9)
        if ds == "f2":
            d = f2_data()
            panels = [("exit (x = 0.610 m)", d.t, d.Twall_exit_avg_K, "T_exit")]
        elif ds == "darr":
            d = darr_data(case)
            panels = [("TC1-3 (0.149 m)", d.t, d["Twall_TC1-3_avg_K"], "T_up"), ("TC4-6 (0.401 m)", d.t, d["Twall_TC4-6_avg_K"], "T_dn")]
        else:
            d = nbs_data(case)
            panels = []
            for st, x in NBS_STATIONS.items():
                dd = d[np.isclose(d.station_m, x, atol=0.3)]
                if dd.empty:
                    dd = d[np.isclose(d.station_m, 60.4, atol=0.3)]
                dd = dd.sort_values("time_s")
                panels.append((f"{st} ({x} m)", dd.time_s, dd.T_wall_K, f"T_{st}"))
        fig, axs = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 3.8), sharey=True, squeeze=False)
        for ax, (title, td, Td, key) in zip(axs[0], panels):
            ax.plot(td, Td, "k-", lw=2, label="data")
            for i, (v, ar) in enumerate(cur):
                ax.plot(ar["t"], ar[key], lw=1.3, color=f"C{order.index(v) if v in order else 9}", label=v)
            ax.set_title(title, fontsize=9)
            ax.set_xlabel("t [s]")
            ax.grid(alpha=0.3)
        axs[0][0].set_ylabel("wall T [K]")
        axs[0][-1].legend(fontsize=7)
        name = {"f2": "Hartwig 2022 F2 (0.05 g)", "darr": f"Darr 2016 IJHMT test {case} (1 g up)",
                "nbs": f"NBS 9264 fig {case}"}[ds]
        fig.suptitle(name, fontsize=10)
        fig.tight_layout()
        fig.savefig(figdir / f"{ds}_{case}.png", dpi=120)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "report"])
    ap.add_argument("--set", nargs="+", default=["f2", "darr", "nbs_ln2", "nbs_lh2"])
    ap.add_argument("--tag", default="p4v")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=3600.0)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        run_jobs(jobs_for(a.set, a.timeout), a.tag, max_workers=min(a.workers, 6))
    df, agg = report(a.tag)
    print(agg.to_string())


if __name__ == "__main__":
    main()
