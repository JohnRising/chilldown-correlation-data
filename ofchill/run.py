"""CLI:  .venv/bin/python -m ofchill.run --cases 7 8 9 10 --closure miropolskii --N 10 --dt 0.1 --tag baseline

Runs every (case x closure) combination in parallel (cached), then writes
results/<tag>/{metrics.csv, metrics.md, runs.md, figs/case*.png, figs/summary.png}
covering ALL runs registered under the tag (repeated invocations accumulate).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .builder import build_uf_case
from .closures import get_closure
from .metrics import aggregate, compare, markdown_table
from .plotting import plot_case, plot_summary
from .runner import RESULTS, load_index, load_run, make_job, run_jobs


def spec_suffix(o: dict) -> str:
    s = ""
    if o.get("g_shift_s"):
        s += f"_gs{o['g_shift_s']:g}"
    if o.get("inlet_quality"):
        s += f"_x{o['inlet_quality']:g}"
    if o.get("exit_p_scale", 1.0) != 1.0:
        s += f"_pex{o['exit_p_scale']:g}"
    if o.get("end_time") is not None:
        s += f"_tend{o['end_time']:g}"
    if not o.get("inertia", True):
        s += "_noin"
    return s


def build_jobs(cases, closures, N, dt, options: dict, params: dict, timeout: float) -> list[dict]:
    jobs = []
    for cname in closures:
        cl = get_closure(cname, **params)
        for c in cases:
            m = build_uf_case(c, N=N, dt=dt, closure=cl, **options)
            label = cl.name + ("" if not params else "_" + "_".join(f"{k}{v:g}" for k, v in params.items()))
            rid = f"c{c}_{label}_N{N}_dt{dt:g}{spec_suffix(options)}"
            spec = {"case": c, "closure": cl.to_spec(), "label": label + spec_suffix(options), "N": N, "dt": dt,
                    "options": options}
            jobs.append(make_job(rid, m, spec, timeout=timeout))
    return jobs


def report(tag: str) -> list[dict]:
    idx = load_index(tag)
    tagdir = RESULTS / tag
    rows, runrows = [], []
    by_case: dict[int, list] = {}
    for rid, ent in sorted(idx.items()):
        spec = ent["spec"]
        try:
            meta, arrs = load_run(tag, rid)
        except FileNotFoundError:
            continue
        c = spec["case"]
        lab = f"{spec['label']} N{spec['N']} dt{spec['dt']:g}"
        diag = meta.get("diag") or {}
        runrows.append(dict(run_id=rid, status=meta.get("status"), converged=meta.get("converged"),
                            n_steps=meta.get("n_steps"), n_stalled=meta.get("n_stalled"),
                            solve_s=round(meta.get("solve_s", 0), 1), h_floor_frac=meta.get("h_at_floor_frac"),
                            E_fluid_m1=diag.get("fluidGain_over_wallEnergy_m1"),
                            E_wallQ_m1=diag.get("wallHeat_over_wallEnergy_m1"),
                            dM_kg=diag.get("mass_imbalance_kg"), error=meta.get("error", "")[:120]))
        base = dict(case=c, label=lab, run_id=rid, status=meta.get("status"), n_stalled=meta.get("n_stalled"),
                    solve_s=meta.get("solve_s"), energy=diag.get("fluidGain_over_wallEnergy_m1"))
        if meta.get("status") != "ok" or arrs is None:
            for st in ("up", "dn"):
                rows.append({**base, "station": st})
            continue
        cmpd = compare(c, arrs)
        for st in ("up", "dn"):
            rows.append({**base, "station": st, **cmpd[st]})
        by_case.setdefault(c, []).append((lab, arrs))
    df = pd.DataFrame(rows)
    df.to_csv(tagdir / "metrics.csv", index=False)
    rows_sorted = sorted(rows, key=lambda r: (r["case"], r["label"], r["station"] != "up"))
    (tagdir / "metrics.md").write_text(markdown_table(rows_sorted) + "\n")
    pd.DataFrame(runrows).to_csv(tagdir / "runs.csv", index=False)
    if _has_tabulate():
        (tagdir / "runs.md").write_text(pd.DataFrame(runrows).to_markdown(index=False) + "\n")
    for c, curves in sorted(by_case.items()):
        plot_case(c, curves, tagdir / "figs" / f"case{c}.png", title=f"UF 2016 case {c} ({tag})")
    if rows:
        plot_summary(rows, tagdir / "figs" / "summary.png")
        ag = aggregate(rows)
        ag.to_csv(tagdir / "aggregate.csv", index=False)
        if _has_tabulate():
            (tagdir / "aggregate.md").write_text(ag.to_markdown(index=False, floatfmt=".2f") + "\n")
    return rows


def _has_tabulate() -> bool:
    try:
        import tabulate  # noqa: F401

        return True
    except ImportError:
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ofchill.run", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--cases", type=int, nargs="+", default=[7, 8, 9, 10])
    ap.add_argument("--closure", nargs="+", default=["miropolskii"],
                    help="named: miropolskii dittusBoelter darrHartwig ttWf; custom: db_vapor miropolskii_expr")
    ap.add_argument("--param", action="append", default=[], help="closure param k=v (repeatable)")
    ap.add_argument("--N", type=int, default=10)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--g-shift", type=float, default=0.0)
    ap.add_argument("--inlet-quality", type=float, default=0.0)
    ap.add_argument("--exit-p-scale", type=float, default=1.0)
    ap.add_argument("--end-time", type=float, default=None)
    ap.add_argument("--no-inertia", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=1500.0, help="per-solve timeout [s]")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args(argv)
    if not a.report_only:
        options = dict(g_shift_s=a.g_shift, inlet_quality=a.inlet_quality, exit_p_scale=a.exit_p_scale,
                       end_time=a.end_time, inertia=not a.no_inertia)
        params = {k: float(v) for k, v in (p.split("=", 1) for p in a.param)}
        jobs = build_jobs(a.cases, a.closure, a.N, a.dt, options, params, a.timeout)
        run_jobs(jobs, a.tag, max_workers=min(a.workers, 6), force=a.force, retry_failed=a.retry_failed)
    rows = report(a.tag)
    print((RESULTS / a.tag / "metrics.md").read_text())
    print(f"-> {RESULTS / a.tag}")


if __name__ == "__main__":
    main()
