"""Reproduce the prior TypeScript runs (OpenFLUME-uf, uf-run/out/*.json) with
the pip package, using the prior agent's input JSON (full-precision xls data)
so that only the solver/closure differs.

  .venv/bin/python -m ofchill.repro_prior [--uf-run DIR] [--tag repro]

Pairs (our closure -> prior file):
  miropolskii       -> m<c>_N10.json   (named model; expect bit-identical)
  miropolskii_expr  -> m8_N10.json     (custom-expression re-implementation)
  db_vapor          -> d<c>_N10.json   (custom expression vs the PATCHED dispatcher)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .builder import build_uf_case
from .closures import get_closure
from .metrics import compare
from .runner import RESULTS, load_run, make_job, run_jobs

DEFAULT_UF_RUN = "/var/folders/sj/sc7cl8bx2f18cmk66lplkzyc0000gn/T/opencode/uf-run"


def diff_vs_prior(arrs: dict, prior: dict) -> dict:
    tp = np.asarray(prior["timesS"])
    n = min(len(tp), len(arrs["t"]))
    Wp = np.asarray(prior["wallTracesK"], float)[:, :n]
    W = arrs["wall_T"][:, :n]
    hp = np.asarray(prior["hWm2K"], float)[:, :n]
    return {
        "n_steps": [int(len(arrs["t"]) - 1), int(prior["nSteps"])],
        "max_abs_dTwall_K": float(np.max(np.abs(W - Wp))),
        "max_abs_dTwall_end_K": float(np.max(np.abs(W[:, -1] - Wp[:, -1]))),
        "max_rel_dh": float(np.nanmax(np.abs(arrs["h"][:, :n] - hp) / np.maximum(hp, 1.0))),
        "max_abs_dmdot_out": float(np.max(np.abs(arrs["mdot_out"][:n] - np.asarray(prior["mdotOut"])[:n]))),
        "prior_wall_s": prior.get("wallS"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--uf-run", default=DEFAULT_UF_RUN)
    ap.add_argument("--tag", default="repro")
    ap.add_argument("--cases", type=int, nargs="+", default=[7, 8, 9, 10])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    U = Path(a.uf_run)
    pairs = []
    for c in a.cases:
        pairs.append(("miropolskii", c, U / "out" / f"m{c}_N10.json"))
        pairs.append(("db_vapor", c, U / "out" / f"d{c}_N10.json"))
    pairs.append(("miropolskii_expr", 8, U / "out" / "m8_N10.json"))
    jobs, meta = [], {}
    for cname, c, pf in pairs:
        if not pf.exists():
            continue
        m = build_uf_case(c, N=10, dt=0.1, closure=get_closure(cname), data_json=str(U / f"uf_case{c}.json"))
        rid = f"c{c}_{cname}_priorjson"
        jobs.append(make_job(rid, m, {"case": c, "closure": cname, "label": f"{cname}[prior-json]", "N": 10,
                                      "dt": 0.1, "options": {"data_json": True}}))
        meta[rid] = pf
    run_jobs(jobs, a.tag, max_workers=a.workers)
    out = {}
    for rid, pf in meta.items():
        m, arrs = load_run(a.tag, rid)
        if arrs is None:
            out[rid] = {"status": m.get("status"), "error": m.get("error")}
            continue
        prior = json.loads(pf.read_text())
        d = diff_vs_prior(arrs, prior)
        d["solve_s"] = m.get("solve_s")
        d["prior_file"] = pf.name
        cm = compare(int(rid.split("_")[0][1:]), arrs)
        d["rms_up_dn"] = [round(cm["up"]["rms"], 2), round(cm["dn"]["rms"], 2)]
        out[rid] = d
        print(rid, json.dumps(d))
    (RESULTS / a.tag / "repro_vs_prior.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
