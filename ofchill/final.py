"""Phase 5a production runs (paper).  Results under results/final/<part>/.

  .venv/bin/python -m ofchill.final run --part uf sens num nbs
  .venv/bin/python -m ofchill.final report

uf    UF 10 cases, G(t + 2.5 s), N = 20, dt = 0.05: miropolskii (also N = 10, dt = 0.1), db_vapor,
      darr2016 published (ds) and x1.19 (d119s) [film law; spec rewet/transition/NB/CHF], spec, spec10
sens  spec10 boundary-condition / closure-constant sensitivities on all 10 cases
num   spec10 N = 10/20/40 (dt 0.05) and dt = 0.1/0.05/0.025 (N 20) on cases 3, 7, 10
nbs   NBS LN2 sat. figs 10, 12: spec10 / miropolskii x inlet K {1, 2, 5, 10} at 101.3 kPa, K = 2 at 83.1 kPa
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from .builder import build_uf_case
from .closures import N2_TC, boiling_variant, get_closure
from .metrics import compare, cross
from .phase2 import variant_for_case
from .runner import RESULTS, load_index, load_run, make_job, run_jobs

LAG = 2.5
PAR = (98.6, 2.2, 293.0)
CASES = list(range(1, 11))


def _beta(c):
    return variant_for_case("spec", c).params["beta"]


def uf_closure(name, c, **params):
    if name in ("spec", "spec10"):
        return variant_for_case(name, c, **params)
    if name == "ds":
        return boiling_variant("d119s", beta=_beta(c), C_fb=1.0, **params).with_label("vds")
    if name == "d119s":
        return boiling_variant("d119s", beta=_beta(c), **params)
    return get_closure(name)


def uf_job(rid, c, cl_name, N=20, dt=0.05, timeout=1500.0, cparams=None, **opts):
    opts.setdefault("g_shift_s", LAG)
    cl = uf_closure(cl_name, c, **(cparams or {}))
    m = build_uf_case(c, N=N, dt=dt, closure=cl, **opts)
    spec = {"case": c, "closure": cl_name, "label": rid.split("__")[1], "N": N, "dt": dt,
            "options": {k: (list(v) if isinstance(v, tuple) else v) for k, v in opts.items()}, "cparams": cparams}
    return make_job(rid, m, spec, timeout=timeout)


def jobs(part):
    J = []
    if part == "uf":
        for c in CASES:
            for v in ("spec10", "spec", "db_vapor", "ds", "d119s"):
                J.append(uf_job(f"c{c}__{v}", c, v))
            J.append(uf_job(f"c{c}__miropolskii_N10dt0.1", c, "miropolskii", N=10, dt=0.1, timeout=1200))
            J.append(uf_job(f"c{c}__miropolskii", c, "miropolskii", timeout=900))
    elif part == "sens":
        for c in CASES:
            b = _beta(c)
            S = {
                "lag0": dict(opts=dict(g_shift_s=0.0)), "lag3.5": dict(opts=dict(g_shift_s=3.5)),
                "x0.02": dict(opts=dict(inlet_quality=0.02)), "x0.05": dict(opts=dict(inlet_quality=0.05)),
                "Pex0.95": dict(opts=dict(exit_p_scale=0.95)), "Pex1.05": dict(opts=dict(exit_p_scale=1.05)),
                # T_MFB = Tc + beta (Tc - T_l): shifting the Tc param by d/(1+beta) moves T_MFB by d exactly
                "TMFB-5": dict(cp=dict(Tc=N2_TC - 5 / (1 + b))), "TMFB+5": dict(cp=dict(Tc=N2_TC + 5 / (1 + b))),
                "C14.9": dict(cp=dict(C_fb=14.9)), "C16.5": dict(cp=dict(C_fb=16.5)),
                "Lam0.080": dict(cp=dict(Lam=0.080)), "Lam0.121": dict(cp=dict(Lam=0.121)),
                "parasitic": dict(opts=dict(parasitic=PAR)), "no_axial": dict(opts=dict(axial_conduction=False)),
            }
            for k, v in S.items():
                J.append(uf_job(f"c{c}__spec10_{k}", c, "spec10", cparams=v.get("cp"), **v.get("opts", {})))
    elif part == "num":
        for c in (3, 7, 10):
            for N, dt in ((10, 0.05), (40, 0.05), (20, 0.1), (20, 0.025), (20, 0.05)):
                J.append(uf_job(f"c{c}__spec10_N{N}_dt{dt:g}", c, "spec10", N=N, dt=dt, timeout=2400))
    elif part == "nbs":
        from .validate4 import build_nbs

        for fig in (10, 12):
            for v in ("spec10", "miropolskii"):
                for K, Po in ((1.0, 101325.0), (2.0, 101325.0), (5.0, 101325.0), (10.0, 101325.0), (2.0, 83082.0)):
                    m = build_nbs(fig, v, K=K, P_out=Po)
                    rid = f"nbs{fig}__{v}_K{K:g}_Pout{Po / 1e3:.0f}"
                    J.append(make_job(rid, m, {"dataset": "nbs", "fig": fig, "variant": v, "K": K, "P_out": Po},
                                      timeout=1800))
    return J


# ------------------------------------------------------------------ reporting
def front_speed(ar, th=130.0):
    """model front speed [mm/s] between the two UF stations from station t130 crossings."""
    a, b = cross(ar["t"], ar["T_up"], th), cross(ar["t"], ar["T_dn"], th)
    return None if (a is None or b is None or b <= a) else 1e3 * 0.252 / (b - a)


def uf_rows(tag):
    rows = []
    for rid, ent in sorted(load_index(tag).items()):
        sp = ent["spec"]
        try:
            meta, ar = load_run(tag, rid)
        except FileNotFoundError:
            continue
        base = dict(run=rid.split("__")[1], case=sp["case"], status=meta.get("status"),
                    solve_s=meta.get("solve_s"), last_t=meta.get("last_sim_t"), n_stalled=meta.get("n_stalled"))
        if ar is None:
            rows += [{**base, "station": st} for st in ("up", "dn")]
            continue
        cm = compare(sp["case"], ar)
        fs = front_speed(ar)
        from .data import load_uf_case

        d = load_uf_case(sp["case"])
        fsd = None
        a, b = cross(d.t, d.station_mean("up"), 130), cross(d.t, d.station_mean("dn"), 130)
        if a is not None and b is not None and b > a:
            fsd = 1e3 * 0.252 / (b - a)
        for st in ("up", "dn"):
            r = cm[st]
            rows.append({**base, "station": st, "rms": r["rms"], "bias": r["bias"],
                         "t130_m": r["t_130_m"], "t130_d": r["t_130_d"], "t100s_m": r["t_100s_m"], "t100s_d": r["t_100s_d"],
                         "tq_m": r["t_quench_m"], "tq_d": r["t_quench_d"], "front_m": fs, "front_d": fsd})
    df = pd.DataFrame(rows)
    for k in ("t130", "t100s", "tq"):
        if f"{k}_m" in df:
            df[f"{k}_err"] = df[f"{k}_m"] - df[f"{k}_d"]
    return df


def aggregate(df):
    def mae(s):
        s = s.dropna().astype(float)
        return s.abs().mean() if len(s) else np.nan

    out = []
    for run, g in df.groupby("run"):
        ok = g[g.status == "ok"]
        both = ok[ok.t130_m.notna() & ok.t130_d.notna()]
        out.append(dict(run=run, n_ok=len(ok) // 2, RMS=ok.rms.mean(), RMS_med=ok.rms.median(), bias=ok.bias.mean(),
                        t130_MAE=mae(ok.t130_err), t130_bias=ok.t130_err.mean(),
                        t130_miss=int((ok.t130_m.notna() ^ ok.t130_d.notna()).sum()),
                        t100s_MAE=mae(ok.t100s_err), tq_MAE=mae(both.tq_err),
                        front_MAE=mae((ok.front_m - ok.front_d).drop_duplicates()),
                        solve_s=g.solve_s.mean(), fails=int((g.status != "ok").sum() // 2)))
    return pd.DataFrame(out).sort_values("RMS")


def report():
    out = RESULTS / "final"
    for part in ("uf", "sens", "num"):
        tag = f"final/{part}"
        if not (RESULTS / tag / "index.json").exists():
            continue
        df = uf_rows(tag)
        if df.empty or "run" not in df:
            continue
        df.to_csv(out / f"{part}_per_station.csv", index=False)
        ag = aggregate(df)
        ag.to_csv(out / f"{part}_aggregate.csv", index=False)
        (out / f"{part}_aggregate.md").write_text(ag.to_markdown(index=False, floatfmt=".2f") + "\n")
        if part == "uf":
            cols = ["case", "station", "run", "rms", "bias", "t130_err", "t100s_err", "tq_err", "front_m", "front_d"]
            (out / "uf_per_station.md").write_text(df[cols].to_markdown(index=False, floatfmt=".2f") + "\n")
        print(part)
        print(ag.round(2).to_string())
    tag = "final/nbs"
    if (RESULTS / tag / "index.json").exists():
        from .validate4 import run_metrics

        rows = []
        for rid, ent in sorted(load_index(tag).items()):
            try:
                meta, ar = load_run(tag, rid)
            except FileNotFoundError:
                continue
            sp = ent["spec"]
            base = dict(run=rid, fig=sp["fig"], variant=sp["variant"], K=sp["K"], P_out_kPa=sp["P_out"] / 1e3,
                        status=meta.get("status"), solve_s=meta.get("solve_s"), last_t=meta.get("last_sim_t"))
            if ar is None:
                rows.append(base)
                continue
            for r in run_metrics({"dataset": "nbs", "fig": sp["fig"]}, ar):
                rows.append({**base, **r})
        df = pd.DataFrame(rows)
        df.to_csv(out / "nbs_per_station.csv", index=False)
        ok = df[df.status == "ok"]
        ag = ok.groupby(["variant", "K", "P_out_kPa"]).agg(
            n=("rms", "size"), RMS=("rms", "mean"), bias=("bias", "mean"),
            knee_MAE=("t_knee_err", lambda s: s.abs().mean()), knee_bias=("t_knee_err", "mean"),
            knee_n=("t_knee_err", lambda s: int(s.notna().sum())), solve_s=("solve_s", "mean")).reset_index()
        fails = df[df.status != "ok"][["run", "status", "last_t"]].drop_duplicates()
        ag.to_csv(out / "nbs_aggregate.csv", index=False)
        (out / "nbs_aggregate.md").write_text(ag.to_markdown(index=False, floatfmt=".2f") + "\n\nfailures:\n"
                                              + fails.to_markdown(index=False) + "\n")
        print(ag.round(2).to_string())
        print(fails)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "report"])
    ap.add_argument("--part", nargs="+", default=["uf", "num", "nbs", "sens"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        for part in a.part:
            run_jobs(jobs(part), f"final/{part}", max_workers=a.workers)
    report()


if __name__ == "__main__":
    main()
