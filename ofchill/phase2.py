"""Phase 2: unified boiling-curve closure variants on the UF cases (1-D).

  .venv/bin/python -m ofchill.phase2 --variants a b c_rl c_rv d --cases 1 2 ... 10 --tag p2 [--N 20 --dt 0.05]
  .venv/bin/python -m ofchill.run --tag p2 --report-only      # tables / figures

Inlet G = meter G(t + 2.5 s) (held at the end of the record), sat. liquid at P1(t), exit P2(t).
Variant 'b' (contact_Tc) gets beta = e_l/e_w per case (sat. liquid at mean P1, 304SS at 130 K).
"""

from __future__ import annotations

import argparse

import numpy as np

from .builder import build_uf_case
from .closures import boiling_variant, contact_beta
from .data import load_uf_case
from .runner import make_job, run_jobs

G_LAG = 2.5


SPEC_PARASITIC = (98.6, 2.2, 293.0)


def variant_for_case(v: str, case: int, **params):
    if v in ("b", "spec") or v.startswith("b_") or v.startswith("spec"):
        P = float(np.mean(load_uf_case(case).col("P1_kPa"))) * 1e3
        if v.startswith("spec"):  # beta = e_l / e_w(T_MFB), fixed point (CLOSURE_SPEC §2)
            import CoolProp.CoolProp as CP

            Tc, Ts = 126.192, CP.PropsSI("T", "P", P, "Q", 0, "Nitrogen")
            Tm = Tc + 0.1 * (Tc - Ts)
            for _ in range(4):
                b = contact_beta(P, T_wall=Tm)
                Tm = Tc + b * (Tc - Ts)
            params.setdefault("beta", b)
        else:
            params.setdefault("beta", contact_beta(P))
    return boiling_variant(v, **params)


def jobs_for(variants, cases, N, dt, g_lag=G_LAG, timeout=5400.0, tag_extra="", parasitic=False, **params):
    jobs = []
    for v in variants:
        for c in cases:
            cl = variant_for_case(v, c, **params)
            kw = {"parasitic": SPEC_PARASITIC} if parasitic else {}
            m = build_uf_case(c, N=N, dt=dt, closure=cl, g_shift_s=g_lag, **kw)
            lab = f"v{v}{tag_extra}" + ("_par" if parasitic else "")
            rid = f"c{c}_{lab}_N{N}_dt{dt:g}_gs{g_lag:g}"
            spec = {"case": c, "closure": cl.to_spec(), "label": lab, "N": N, "dt": dt,
                    "options": {"g_shift_s": g_lag, **params}}
            jobs.append(make_job(rid, m, spec, timeout=timeout))
    return jobs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="+", default=["a", "b", "c_rl", "c_rv", "d"])
    ap.add_argument("--cases", type=int, nargs="+", default=list(range(1, 11)))
    ap.add_argument("--N", type=int, default=20)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--g-lag", type=float, default=G_LAG)
    ap.add_argument("--tag", default="p2")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=5400.0)
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--parasitic", action="store_true", help="CLOSURE_SPEC §9 rig wall heat leak")
    a = ap.parse_args(argv)
    params = {k: float(v) for k, v in (p.split("=", 1) for p in a.param)}
    extra = "" if not params else "_" + "_".join(f"{k}{v:g}" for k, v in params.items())
    jobs = jobs_for(a.variants, a.cases, a.N, a.dt, a.g_lag, a.timeout, extra, a.parasitic, **params)
    run_jobs(jobs, a.tag, max_workers=min(a.workers, 6))


if __name__ == "__main__":
    main()
