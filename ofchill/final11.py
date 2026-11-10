"""Phase 8: consistency re-runs with the final closure spec11 (§10 film law + §11 R5 rewet).
Results under results/final11/<part>/ (same table formats as ofchill.final).

  .venv/bin/python -m ofchill.final11 run --part uf sens num val
  .venv/bin/python -m ofchill.final11 report

uf    UF 10 tests, G(t+2.5 s), N=20, dt=0.05: spec11, specz11 (z-independent Pe^1/3 law + R5),
      ds11 / d119s11 (Darr 2016 published / x1.19 film law + R5).  Old-rewet rows (db_vapor, miropolskii,
      ds, d119s) are in results/final/uf and are merged into the tables by the report.
sens  spec11 sensitivities (Phase-5a list; C LOTO 14.9-16.4, Lambda 0.081-0.118)
num   spec11 N = 10/20/40, dt = 0.1/0.05/0.025 on tests 3, 7, 10
val   F2 (0.05 g) and Darr 1-g tests 31/40 with spec11 (+gravity), specz11 (+gravity), ds11, d119s11,
      db_vapor, miropolskii (N10 dt0.1)
"""
from __future__ import annotations

import argparse
import math

import numpy as np
import pandas as pd

from . import final as F
from .builder import build_uf_case
from .closures import N2_TC, boiling_variant, get_closure
from .phase2 import variant_for_case
from .runner import RESULTS, make_job, run_jobs

TAG = "final11"
LAG = 2.5
CASES = list(range(1, 11))


def _beta(c):
    return variant_for_case("spec", c).params["beta"]


def uf_closure(name, c, **cp):
    if name in ("spec11", "specz11", "spec10", "spec"):
        return variant_for_case(name, c, **cp)
    if name in ("ds11", "d119s11"):
        return boiling_variant(name, beta=_beta(c), **cp)
    return get_closure(name)


def uf_job(rid, c, cl, N=20, dt=0.05, timeout=1500.0, cparams=None, **opts):
    opts.setdefault("g_shift_s", LAG)
    m = build_uf_case(c, N=N, dt=dt, closure=uf_closure(cl, c, **(cparams or {})), **opts)
    spec = {"case": c, "closure": cl, "label": rid.split("__")[1], "N": N, "dt": dt,
            "options": {k: (list(v) if isinstance(v, tuple) else v) for k, v in opts.items()}, "cparams": cparams}
    return make_job(rid, m, spec, timeout=timeout)


def jobs(part):
    J = []
    if part == "uf":
        for c in CASES:
            for v in ("spec11", "specz11", "ds11", "d119s11"):
                J.append(uf_job(f"c{c}__{v}", c, v))
    elif part == "sens":
        for c in CASES:
            b = _beta(c)
            de = 5.0 / ((1 + b) * N2_TC)  # R5: dT_MFB = (1+beta) Tc de
            S = {"lag0": dict(o=dict(g_shift_s=0.0)), "lag3.5": dict(o=dict(g_shift_s=3.5)),
                 "x0.02": dict(o=dict(inlet_quality=0.02)), "x0.05": dict(o=dict(inlet_quality=0.05)),
                 "Pex0.95": dict(o=dict(exit_p_scale=0.95)), "Pex1.05": dict(o=dict(exit_p_scale=1.05)),
                 "TMFB-5": dict(cp=dict(e_hn=0.122 - de)), "TMFB+5": dict(cp=dict(e_hn=0.122 + de)),
                 "C14.9": dict(cp=dict(C_fb=14.9)), "C16.4": dict(cp=dict(C_fb=16.4)),
                 "Lam0.081": dict(cp=dict(Lam=0.081)), "Lam0.118": dict(cp=dict(Lam=0.118)),
                 "parasitic": dict(o=dict(parasitic=F.PAR)), "no_axial": dict(o=dict(axial_conduction=False))}
            for k, v in S.items():
                J.append(uf_job(f"c{c}__spec11_{k}", c, "spec11", cparams=v.get("cp"), **v.get("o", {})))
    elif part == "num":
        for c in (3, 7, 10):
            for N, dt in ((10, 0.05), (40, 0.05), (20, 0.1), (20, 0.025), (20, 0.05)):
                J.append(uf_job(f"c{c}__spec11_N{N}_dt{dt:g}", c, "spec11", N=N, dt=dt, timeout=2400))
    elif part == "val":
        from .validate4 import build_darr, build_f2

        for v, vg in (("spec11", "spec11_gh"), ("specz11", "specz11_gh"), ("ds11", "ds11"), ("d119s11", "d119s11"),
                      ("db_vapor", "db_vapor"), ("miropolskii", "miropolskii")):
            kw = dict(N=10, dt=0.1) if v == "miropolskii" else {}
            J.append(make_job(f"f2__{v}", build_f2(vg, **kw), {"dataset": "f2", "variant": v}, timeout=1800))
            for t in (31, 40):
                vgv = vg.replace("_gh", "_gv")
                J.append(make_job(f"darr{t}__{v}", build_darr(t, vgv, **kw),
                                  {"dataset": "darr", "test": t, "variant": v}, timeout=1800))
        for t in (31, 40):  # 1-g reference without the buoyancy term
            J.append(make_job(f"darr{t}__spec11_nog", build_darr(t, "spec11"),
                              {"dataset": "darr", "test": t, "variant": "spec11_nog"}, timeout=1800))
    return J


def report():
    out = RESULTS / TAG
    out.mkdir(parents=True, exist_ok=True)
    for part in ("uf", "sens", "num"):
        tag = f"{TAG}/{part}"
        if not (RESULTS / tag / "index.json").exists():
            continue
        df = F.uf_rows(tag)
        if part == "uf" and (RESULTS / "final/uf/index.json").exists():
            old = F.uf_rows("final/uf")
            old = old[old.run.isin(["db_vapor", "miropolskii", "miropolskii_N10dt0.1", "ds", "d119s"])]
            old["run"] = old["run"].replace({"ds": "ds (contact-Tc rewet)", "d119s": "d119s (contact-Tc rewet)"})
            df = pd.concat([df, old], ignore_index=True)
        if df.empty:
            continue
        df.to_csv(out / f"{part}_per_station.csv", index=False)
        ag = F.aggregate(df)
        ag.to_csv(out / f"{part}_aggregate.csv", index=False)
        (out / f"{part}_aggregate.md").write_text(ag.to_markdown(index=False, floatfmt=".2f") + "\n")
        if part == "uf":
            cols = ["case", "station", "run", "rms", "bias", "t130_err", "t100s_err", "tq_err", "front_m", "front_d"]
            (out / "uf_per_station.md").write_text(df[cols].to_markdown(index=False, floatfmt=".2f") + "\n")
        print(part)
        print(ag.round(2).to_string())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "report"])
    ap.add_argument("--part", nargs="+", default=["uf", "val", "num", "sens"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        for part in a.part:
            run_jobs(jobs(part), f"{TAG}/{part}", max_workers=a.workers)
    report()


if __name__ == "__main__":
    main()
