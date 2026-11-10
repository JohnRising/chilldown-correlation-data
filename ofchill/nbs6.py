"""Phase 6: NBS 9264 inverse boundary-condition tests and closure changes for long 1-g lines.

  .venv/bin/python -m ofchill.nbs6 run --batch probe            # job batches defined in BATCHES
  .venv/bin/python -m ofchill.nbs6 report                       # -> results/nbs6/{per_run.md, aggregate.md, figs/}

Job spec strings:  fig:closure[:k=v,...]   with keys K, Pout, N, dt, satin (1 = saturated-inlet fallback),
x_in (inlet quality), Tin (inlet temperature) and any closure param (cp_<name>=value).
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from .runner import RESULTS, load_index, load_run, make_job, run_jobs
from .validate4 import NBS_STATIONS, build_nbs, nbs_data, nbs_meta, run_metrics

import os

TAG = os.environ.get("NBS_TAG", "nbs6")
A_NBS = np.pi / 4 * 0.015875**2


def parse(spec: str):
    parts = spec.split(":")
    fig, cl = int(parts[0]), parts[1]
    kw = {}
    if len(parts) > 2 and parts[2]:
        for item in parts[2].split(","):
            k, v = item.split("=")
            kw[k] = v
    return fig, cl, kw


def build(spec: str):
    fig, cl, kw = parse(spec)
    cparams = {k[3:]: float(v) for k, v in kw.items() if k.startswith("cp_")}
    args = dict(K=float(kw.get("K", 2.0)), P_out=float(kw.get("Pout", 101325.0)))
    if "N" in kw:
        args["N"] = int(kw["N"])
    if "dt" in kw:
        args["dt"] = float(kw["dt"])
    if kw.get("satin") == "1":
        args["force_saturated"] = True
    if "x_in" in kw:
        args["inlet_quality"] = float(kw["x_in"])
    if "Tin" in kw:
        args["inlet_T"] = float(kw["Tin"])
    if kw.get("static") == "1":
        args["static_init"] = True
    if "ramp" in kw:
        args["ramp_s"] = float(kw["ramp"])
    if "f" in kw:
        args["friction_factor"] = float(kw["f"])
    if "rough" in kw:
        args["roughness"] = float(kw["rough"])
    if "end" in kw:
        args["end_time"] = float(kw["end"])
    return build_nbs(fig, cl, cparams=cparams or None, **args)


def job(spec: str, timeout: float = 2400.0):
    fig, cl, kw = parse(spec)
    seg = float(kw["seg"]) if "seg" in kw else None
    return make_job(spec, build(spec), {"dataset": "nbs", "fig": fig, "variant": cl, "kw": kw}, timeout=timeout,
                    segment_s=seg)


K_GRID = (0.01, 2.0, 20.0, 100.0)
BATCHES = {
    "probe": [f"10:spec10_gh:K={k}" for k in K_GRID] + ["10:spec10_gh:K=0.01,Pout=83082",
                                                        "10:spec10_gh:cp_C_fb=7.9", "10:spec10_gh:cp_C_fb=31.6"]
             + [f"2:spec10_gh:K={k}" for k in K_GRID],
    "ic": [f"{f}:{c}:Pout=83082,static=1,ramp=1" for f in (10, 2) for c in ("spec10_gh", "spec11_gh", "miropolskii")]
          + [f"{f}:spec11_gh:Pout=83082" for f in (10, 2)],
    # final: extended horizon (1.6x record) so late model knees are measured
    "final_sat": [f"{f}:{c}:Pout=83082,end={e}" for f, e in ((10, 400), (11, 384), (12, 208), (2, 128))
                  for c in ("spec11_gh", "spec11n_gh", "miropolskii")],
    "kfit": [f"{f}:spec11_gh:Pout=83082,K={k},end={e}" for f, e in ((11, 384), (12, 208), (2, 128))
             for k in (0.01, 20.0, 100.0)],
    # phase 7
    "p7_cmp": [f"{f}:{c}:Pout=83082,end={e}" for f, e in ((10, 400), (11, 384), (12, 208))
               for c in ("d119s_gh", "spec_gh")],
    "p7_fric": [f"{f}:spec11_gh:Pout=83082,end={e},f={ff}" for f, e in ((10, 400), (11, 384))
                for ff in (0.017, 0.0085, 0.0043)],
    "final_sub": [f"{f}:spec11_gh:Pout=83082,seg=15,end={e}" for f, e in ((13, 211), (14, 144), (4, 120), (5, 104),
                                                                          (6, 70), (7, 45))]
                 + [f"{f}:miropolskii:Pout=83082,seg=15,end={e}" for f, e in ((13, 211), (14, 144))],
    "strat": [f"{f}:{c}:Pout=83082" for f in (10, 12, 2) for c in ("spec11s_gh", "spec11sn_gh", "spec11n_gh")]
             + ["12:spec11_gh:Pout=83082"],
}


def per_run(tag=TAG) -> pd.DataFrame:
    rows = []
    for rid, ent in sorted(load_index(tag).items()):
        try:
            meta, ar = load_run(tag, rid)
        except FileNotFoundError:
            continue
        sp = ent["spec"]
        base = dict(run=rid, fig=sp["fig"], variant=sp["variant"], status=meta.get("status"),
                    last_t=meta.get("last_sim_t"), solve_s=meta.get("solve_s"))
        if ar is None:
            rows.append(base)
            continue
        t = ar["t"]
        m_in = ar["mdot_in"]
        base.update(mdot_mean=float(np.trapezoid(m_in, t) / (t[-1] - t[0])), mdot_max=float(m_in[1:].max()),
                    mdot_end=float(np.mean(m_in[-5:])))
        ms = run_metrics({"dataset": "nbs", "fig": sp["fig"]}, ar)
        base["RMS"] = float(np.mean([r["rms"] for r in ms]))
        base["RMS70"] = float(np.mean([r["rms70"] for r in ms]))
        t_key = "t100_err" if sp["fig"] < 10 else "t150_err"
        te = [r.get(t_key) for r in ms if r.get(t_key) is not None]
        base["tmid_MAE"] = float(np.mean(np.abs(te))) if te else np.nan
        base["tmid_n"] = len(te)
        base["bias"] = float(np.mean([r["bias"] for r in ms]))
        for r in ms:
            base[f"rms_{r['station']}"] = r["rms"]
            base[f"knee_err_{r['station']}"] = r.get("t_knee_err")
            base[f"onset_err_{r['station']}"] = r.get("onset_err")
        oe = [r.get("on20_err") for r in ms[1:] if r.get("on20_err") is not None]  # stations 2-4
        base["on20_MAE_s234"] = float(np.mean(np.abs(oe))) if oe else np.nan
        base["on20_bias_s234"] = float(np.mean(oe)) if oe else np.nan
        for r in ms:
            base[f"on20_err_{r['station']}"] = r.get("on20_err")
        ke = [r.get("t_knee_err") for r in ms if r.get("t_knee_err") is not None]
        base["knee_MAE"] = float(np.mean(np.abs(ke))) if ke else np.nan
        base["knee_n"] = len(ke)
        rows.append(base)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "report"])
    ap.add_argument("--batch", nargs="+", default=[])
    ap.add_argument("--spec", nargs="+", default=[])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=2400.0)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        specs = [s for b in a.batch for s in BATCHES[b]] + a.spec
        run_jobs([job(s, a.timeout) for s in specs], TAG, max_workers=min(6, a.workers))
    df = per_run()
    (RESULTS / TAG).mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / TAG / "per_run.csv", index=False)
    cols = [c for c in ["run", "status", "last_t", "RMS", "RMS70", "bias", "knee_MAE", "knee_n", "on20_MAE_s234",
                        "on20_bias_s234", "tmid_MAE", "tmid_n", "mdot_mean", "mdot_max",
                        "mdot_end", "rms_s1", "rms_s2", "rms_s3", "rms_s4", "knee_err_s1", "knee_err_s2",
                        "knee_err_s3", "knee_err_s4", "solve_s"] if c in df]
    (RESULTS / TAG / "per_run.md").write_text(df[cols].to_markdown(index=False, floatfmt=".3g") + "\n")
    print(df[cols].round(3).to_string())


if __name__ == "__main__":
    main()
