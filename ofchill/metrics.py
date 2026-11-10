"""Model-vs-data metrics (definitions follow the prior uf_metrics.py / uf_post.py).

Per case and station (data = mean of the 3 TCs; model = wall T linearly
interpolated between wall nodes at the station x):
  rms, bias        over data samples with t <= model end time (model linearly
                   interpolated to the data times); bias = mean(model - data)
  t_<T>            first downward crossing of threshold T (linear interp.)
  slope_0_2        least-squares slope of T(t) over 0 <= t <= 2 s [K/s]
  t_quench         time of the maximum cooling rate (data: 5-sample centred
                   rolling mean then gradient; model: gradient of the 0.1-s
                   trace smoothed with a 3-point centred mean), T_quench = T there,
                   rate_max = max cooling rate [K/s]
  T_end            T at the end of the (model) record
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .data import load_uf_case

THRESHOLDS = (200.0, 150.0, 130.0, 100.0)


def cross(t, T, th):
    t, T = np.asarray(t, float), np.asarray(T, float)
    i = np.where(T < th)[0]
    if len(i) == 0:
        return None
    i = i[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (T[i - 1] - th) / (T[i - 1] - T[i]) * (t[i] - t[i - 1]))


def slope(t, T, t0=0.0, t1=2.0):
    t, T = np.asarray(t), np.asarray(T)
    m = (t >= t0) & (t <= t1)
    if m.sum() < 2:
        return None
    return float(np.polyfit(t[m], T[m], 1)[0])


def quench(t, T, window: int):
    Ts = pd.Series(T).rolling(window, center=True, min_periods=1).mean().to_numpy()
    dT = np.gradient(Ts, t)
    i = int(np.argmin(dT))
    return float(t[i]), float(T[i]), float(-dT[i])


def trace_metrics(t, T, *, is_data: bool, thresholds=THRESHOLDS) -> dict:
    tq, Tq, rate = quench(t, T, 5 if is_data else 3)
    r = {f"t_{int(th)}": cross(t, T, th) for th in thresholds}
    r.update(slope_0_2=slope(t, T), t_quench=tq, T_quench=Tq, rate_max=rate, T_end=float(T[-1]), T0=float(T[0]))
    return r


def data_metrics(case_no: int, thresholds=THRESHOLDS) -> dict:
    d = load_uf_case(case_no)
    return {st: trace_metrics(d.t, d.station_mean(st), is_data=True, thresholds=thresholds) for st in ("up", "dn")}


def compare(case_no: int, arrs: dict, thresholds=THRESHOLDS) -> dict:
    """Metrics for one model run (arrays from runner.load_run) vs the data."""
    d = load_uf_case(case_no)
    td = d.t
    tm = arrs["t"]
    out = {}
    for st in ("up", "dn"):
        Td = d.station_mean(st)
        Tm = arrs[f"T_{st}"]
        m = td <= tm[-1] + 1e-9
        Tmi = np.interp(td[m], tm, Tm)
        e = Tmi - Td[m]
        dm = trace_metrics(td[m], Td[m], is_data=True, thresholds=thresholds)
        mm = trace_metrics(tm, Tm, is_data=False, thresholds=thresholds)
        r = {"rms": float(np.sqrt(np.mean(e**2))), "bias": float(np.mean(e)), "t_end": float(tm[-1])}
        # t100*: time to max(100 K, Tsat + 8 K) (Tsat = median measured station Tsat)
        th = max(100.0, float(np.nanmedian(d.col(f"Tsat_{st}_K"))) + 8.0)
        r["thr100s"] = th
        r["t_100s_m"] = cross(tm, Tm, th)
        r["t_100s_d"] = cross(td[m], Td[m], th)
        for k in mm:
            r[f"{k}_m"] = mm[k]
            r[f"{k}_d"] = dm[k]
        out[st] = r
    return out


def metrics_table(rows: list[dict], thresholds=(130.0, 100.0)) -> pd.DataFrame:
    """rows: dicts with keys case, label, station + compare() station dict."""
    df = pd.DataFrame(rows)
    return df


def _f(v, fmt="{:.1f}"):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    return fmt.format(v)


def markdown_table(rows: list[dict]) -> str:
    """rows: {case, label, N, dt, station, status, n_stalled, solve_s, energy, **compare-station}."""
    hdr = ("| case | closure | stn | RMS [K] | bias [K] | t130 m/d [s] | t100* m/d [s] | slope0-2 m/d [K/s] | "
           "t_quench m/d [s] | T_quench m/d [K] | T_end m/d [K] | E-bal | stalled | solve [s] |")
    lines = [hdr, "|" + "|".join(["---"] * (hdr.count("|") - 1)) + "|"]
    for r in rows:
        if r.get("status") != "ok":
            lines.append(f"| {r['case']} | {r['label']} | {r.get('station', '')} | {r.get('status')} | "
                         + " | ".join([""] * 10) + " |")
            continue
        eb = r.get("energy")
        lines.append(
            f"| {r['case']} | {r['label']} | {r['station']} | {_f(r['rms'])} | {_f(r['bias'], '{:+.1f}')} | "
            f"{_f(r['t_130_m'], '{:.2f}')}/{_f(r['t_130_d'], '{:.2f}')} | "
            f"{_f(r.get('t_100s_m'), '{:.2f}')}/{_f(r.get('t_100s_d'), '{:.2f}')} | "
            f"{_f(r['slope_0_2_m'])}/{_f(r['slope_0_2_d'])} | "
            f"{_f(r['t_quench_m'], '{:.2f}')}/{_f(r['t_quench_d'], '{:.2f}')} | "
            f"{_f(r['T_quench_m'])}/{_f(r['T_quench_d'])} | {_f(r['T_end_m'])}/{_f(r['T_end_d'])} | "
            f"{_f(eb * 100 if eb is not None else None, '{:+.1f}%')} | {r.get('n_stalled')} | {_f(r.get('solve_s'), '{:.0f}')} |"
        )
    return "\n".join(lines)


def aggregate(rows: list[dict]) -> "pd.DataFrame":
    """Per-label aggregate over case-stations: mean RMS, mean bias, MAE (bias) of t130 / t100*
    errors over stations where both cross, mismatches (crossed in only one), quench-time
    (max cooling rate) MAE over stations where both reach 130 K, mean slope and T_end errors."""
    df = pd.DataFrame([r for r in rows if r.get("status") == "ok"])
    out = []
    for lab, g in df.groupby("label"):
        rec = {"label": lab, "n": len(g), "RMS_mean": g.rms.mean(), "RMS_median": g.rms.median(),
               "bias_mean": g.bias.mean()}
        for k in ("t_130", "t_100s"):
            mm, dd = g[f"{k}_m"], g[f"{k}_d"]
            both = mm.notna() & dd.notna()
            err = (mm[both] - dd[both]).astype(float)
            rec[f"{k}_MAE"] = err.abs().mean() if len(err) else np.nan
            rec[f"{k}_bias"] = err.mean() if len(err) else np.nan
            rec[f"{k}_n"] = int(both.sum())
            rec[f"{k}_mismatch"] = int((mm.notna() ^ dd.notna()).sum())
        q = g[g.t_130_d.notna() & g.t_130_m.notna()]  # both quench
        rec["tq_MAE"] = (q.t_quench_m - q.t_quench_d).abs().mean() if len(q) else np.nan
        rec["Tend_err_mean"] = (g.T_end_m - g.T_end_d).mean()
        rec["slope_err_mean"] = (g.slope_0_2_m - g.slope_0_2_d).mean()
        rec["solve_s_mean"] = g.solve_s.mean() / 1.0
        rec["stalled"] = int(g.n_stalled.fillna(0).sum() / 2)
        out.append(rec)
    return pd.DataFrame(out).sort_values("RMS_mean")
