"""Per-case / per-station preprocessing shared by all analysis scripts.

valid-time window, smoothing, per-TC lumped heat flux, quench (rewet) detection per trace.
"""
from __future__ import annotations

import functools

import numpy as np
import pandas as pd

from common import (CASES, DT, MW_A, STATIONS, TC_OF, load, smooth, ss_cp, station_P)


def t_valid_end(d: pd.DataFrame) -> float:
    """end of the usable record: first time (t>3 s) at which the smoothed inlet pressure P1
    falls below 70 % of its median (flow shut-off / vent; cases 2, 4, 10) minus 0.3 s."""
    t = d.t_s.values
    P = smooth(d.P1_kPa.values, 15, 1)
    m = np.median(P[t > 3])
    bad = np.where((t > 3) & (P < 0.7 * m))[0]
    return (t[bad[0]] - 0.3) if len(bad) else t[-1]


def lumped_q(T, win=5):
    """wall-storage heat flux per unit INNER area from a wall temperature trace [W/m2]
    (thin-wall lumped: q = -m''_w cp(T) dT/dt).  Does not include parasitic heat."""
    dTdt = smooth(np.asarray(T, float), win, 2, deriv=1)
    return -MW_A * ss_cp(T) * dTdt


def detect_quench(t, T, Tsat, tmax, thr_factors=(1.5, 2.0, 3.0)):
    """Detect the rewet/quench in one wall-temperature trace.

    Definitions (documented in REPORT.md):
      r(t) = -dT/dt (Savitzky-Golay, 5 pts, 2nd order).
      quench  : global max of r for t<=tmax with r_max > max(4 r_fb, 10 K/s) and a >=15 K drop
                within +-1 s, where r_fb = median r over [t_q-3, t_q-1] s (film-boiling baseline).
      T_MFB   : wall T at the last sample before t_q where r <= 2 r_fb  (the cooling rate --
                and hence q'' -- has doubled w.r.t. the film-boiling trend), linearly
                interpolated between samples.  Sensitivity: thresholds 1.5x and 3x, and the
                curvature definition (min of d2T/dt2 in [t_q-1.5 s, t_q]).
      CHF     : T and q at t_q (max cooling rate)
    Returns dict or None.
    """
    t = np.asarray(t); T = np.asarray(T, float)
    ok = t <= tmax
    r = -smooth(T, 5, 2, deriv=1)
    r2 = smooth(T, 7, 2, deriv=2)
    idx = np.where(ok & (t > 1.5))[0]
    if len(idx) < 30:
        return None
    k = idx[np.argmax(r[idx])]
    base = (t >= t[k] - 3) & (t <= t[k] - 1)
    if base.sum() < 5:
        return None
    r_fb = max(np.median(r[base]), 0.05)
    win = (t >= t[k] - 1) & (t <= t[k] + 1)
    drop = T[win].max() - T[win].min()
    if not (r[k] > max(4 * r_fb, 10.0) and drop >= 15.0 and T[k] > Tsat[k] + 2):
        return None
    out = dict(t_q=t[k], T_chf=T[k], r_max=r[k], r_fb=r_fb, drop=drop)
    for f in thr_factors:
        j = k
        while j > 0 and r[j] > f * r_fb:
            j -= 1
        # interpolate between j (below thr) and j+1 (above thr)
        a = (f * r_fb - r[j]) / max(r[j + 1] - r[j], 1e-9)
        a = np.clip(a, 0, 1)
        out[f"T_mfb_{f:g}x"] = T[j] + a * (T[j + 1] - T[j])
        out[f"t_mfb_{f:g}x"] = t[j] + a * DT
    cw = (t >= t[k] - 1.5) & (t <= t[k])
    jc = np.where(cw)[0][np.argmin(r2[cw])]
    out["T_mfb_curv"] = T[jc]
    out["t_mfb_curv"] = t[jc]
    out["T_mfb"] = out["T_mfb_2x"]
    out["t_mfb"] = out["t_mfb_2x"]
    return out


@functools.lru_cache(maxsize=None)
def _raw_ratio(case: int) -> float:
    """sum(q_inverse)/sum(q_lumped(Ti)) over clear film boiling (t>1 s, Ti-Tsat>60 K), both stations."""
    d = load(case)
    tend = t_valid_end(d)
    num = den = 0.0
    for st in STATIONS:
        Ti = d[f"Ti_{st}_K"].values
        m = (d.t_s.values > 1) & (d.t_s.values < tend) & (Ti - d[f"Tsat_{st}_K"].values > 60)
        num += d[f"qi_{st}_kWm2"].values[m].sum()
        den += lumped_q(Ti)[m].sum() / 1e3
    return num / den


@functools.lru_cache(maxsize=None)
def inverse_scale(case: int) -> float:
    """The published inverse q'' of cases 1,2,3,5,6 is ~0.74x the wall energy balance
    (-m''cp dT/dt with NIST cp) at all times/temperatures, whereas cases 4,7-10 give ~1.02x
    (the +2 % being the parasitic heat leak).  We define an energy-consistent flux
    q_ec = q_inv / k with k = ratio_case / ratio_ref (ratio_ref = mean of the cases with
    ratio > 0.9); k = 1 for the consistent cases."""
    rs = {c: _raw_ratio(c) for c in CASES}
    ref = np.mean([r for r in rs.values() if r > 0.9])
    r = rs[case]
    return r / ref if r < 0.9 else 1.0


@functools.lru_cache(maxsize=None)
def station_data(case: int, st: str) -> pd.DataFrame:
    """Tidy per-station time series (valid window only)."""
    d = load(case)
    tend = t_valid_end(d)
    d = d[d.t_s <= tend].reset_index(drop=True)
    s = pd.DataFrame({"t": d.t_s.values})
    s["Ti"] = d[f"Ti_{st}_K"].values
    s["Tsat"] = d[f"Tsat_{st}_K"].values
    s["P"] = station_P(d, st)
    s["q"] = d[f"qi_{st}_kWm2"].values * 1e3
    s["dT"] = s.Ti - s.Tsat
    s["h"] = s.q / s.dT
    s["G"] = d.G_kgm2s.values
    s["P1"] = d.P1_kPa.values * 1e3
    for i, c in enumerate(TC_OF[st]):
        s[f"TC{i}"] = d[c].values
        s[f"qTC{i}"] = lumped_q(d[c].values)
    s["q_lump"] = lumped_q(s.Ti.values)
    k = inverse_scale(case)
    s["q_ec"] = s.q / k
    s["h_ec"] = s.q_ec / s.dT
    s.attrs["k_inv"] = k
    s.attrs["t_end"] = tend
    s.attrs["case"] = case
    s.attrs["st"] = st
    return s


@functools.lru_cache(maxsize=None)
def quench_info(case: int, st: str):
    """quench detection on the inverse-method Ti and on each TC."""
    s = station_data(case, st)
    tmax = s.attrs["t_end"]
    res = {"Ti": detect_quench(s.t.values, s.Ti.values, s.Tsat.values, tmax)}
    for i in range(3):
        res[f"TC{i}"] = detect_quench(s.t.values, s[f"TC{i}"].values, s.Tsat.values, tmax)
    return res


def film_mask(case: int, st: str, t_min=0.5, margin_K=15.0):
    """film-boiling sample mask: t>t_min and before the earliest detected MFB onset of ANY TC
    or of Ti at this station (minus 0.5 s) and Ti > T_MFB + margin; q>0."""
    s = station_data(case, st)
    qi = quench_info(case, st)
    m = (s.t.values > t_min) & (s.q.values > 0) & (s.dT.values > 5)
    ts = [v["t_mfb_1.5x"] for v in qi.values() if v]
    Ts = [v["T_mfb_1.5x"] for v in qi.values() if v]
    if ts:
        m &= s.t.values < min(ts) - 0.5
        m &= s.Ti.values > max(Ts) + margin_K
    return m
