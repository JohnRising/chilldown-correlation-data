"""'When is the line chilled down?' — chilldown-completion metrics for every closure and station.

  .venv/bin/python -m ofchill.chilldown     -> results/final11/chilldown_{stations.csv, summary.csv, time.md}

Definitions (same for model and data):
  t_cd     first time the wall temperature falls below T_sat + 10 K.  T_sat is the measured local saturation
           temperature for UF (station Tsat(t)) and F2 (Tsat(P_exit(t))); Tsat(P_out) for the Darr 1-g tests;
           Tsat(P_drive) for NBS (local line pressure not recorded).
  t_q      quench time = time of the maximum cooling rate (data: 5-sample running mean; model: 3-point), only
           counted when the trace drops by > 30 K within +-1.5 s of it (a real quench).
  outcome  hit (both complete within the data record), miss (data complete, model not), false (model
           completes, data not), none (neither).  The model is integrated at least to the end of the data
           record; model completion is judged within the model horizon (t_cd of a 'miss' is censored).
  whole section: last-station t_cd per test (all stations of the test must complete in the data).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .data import load_uf_case
from .metrics import cross
from .runner import RESULTS, load_index, load_run

OUT = RESULTS / "final11"

LABELS = {
    "spec11": "present closure",
    "specz11": "z-independent film law + R5",
    "ds11": "Darr 2016 film law + R5",
    "d119s11": "Darr 2016 x1.19 + R5",
    "ds": "Darr 2016 (contact-Tc rewet)",
    "d119s": "Darr 2016 x1.19 (contact-Tc rewet)",
    "db_vapor": "Dittus-Boelter (sat. vapour)",
    "miropolskii": "Miropolskii",
    "spec_old": "z-independent film law (contact-Tc rewet)",
}


def quench_time(t, T, w):
    s = pd.Series(T).rolling(w, center=True, min_periods=1).mean().to_numpy()
    r = -np.gradient(s, t)
    i = int(np.argmax(r))
    m = (t >= t[i] - 1.5) & (t <= t[i] + 1.5)
    return float(t[i]) if (T[m].max() - T[m].min()) > 30.0 else None


def tcd(t, T, Tsat):
    """first crossing of T below Tsat(t) + 10 K (Tsat array on t or scalar)."""
    thr = np.broadcast_to(np.asarray(Tsat, float), np.shape(T)) + 10.0
    d = np.asarray(T) - thr
    i = np.where(d < 0)[0]
    if len(i) == 0:
        return None
    i = i[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + d[i - 1] / (d[i - 1] - d[i]) * (t[i] - t[i - 1]))


def _run(tag, rid):
    try:
        m, a = load_run(tag, rid)
    except (FileNotFoundError, KeyError):
        return None
    return a if m.get("status") == "ok" else None


def station_rows():
    rows = []

    def add(dataset, test, station, closure, td, Td, Tsd, tm, Tm, Tsm, nd):
        tend = float(td[-1])
        a_d, a_m = tcd(td, Td, Tsd), (tcd(tm, Tm, Tsm) if tm is not None else None)
        qd = quench_time(td, Td, nd)
        qm = quench_time(tm, Tm, 3) if tm is not None else None
        if tm is None:
            outcome = "no run"
        elif a_d is not None and a_m is not None:
            outcome = "hit"
        elif a_d is not None:
            outcome = "miss"
        elif a_m is not None and a_m <= tend:
            outcome = "false"
        else:
            outcome = "none"
        rows.append(dict(dataset=dataset, test=test, station=station, closure=closure, t_end_data=tend,
                         tcd_d=a_d, tcd_m=a_m, tq_d=qd, tq_m=qm, outcome=outcome))

    # ---- UF
    for c in range(1, 11):
        d = load_uf_case(c)
        for cl, (tag, rid) in {
            "spec11": ("final11/uf", f"c{c}__spec11"), "specz11": ("final11/uf", f"c{c}__specz11"),
            "ds11": ("final11/uf", f"c{c}__ds11"), "d119s11": ("final11/uf", f"c{c}__d119s11"),
            "ds": ("final/uf", f"c{c}__ds"), "d119s": ("final/uf", f"c{c}__d119s"),
            "db_vapor": ("final/uf", f"c{c}__db_vapor"), "miropolskii": ("final/uf", f"c{c}__miropolskii_N10dt0.1"),
        }.items():
            a = _run(tag, rid)
            for st in ("up", "dn"):
                Tsd = d.col(f"Tsat_{st}_K")
                if a is None:
                    add("UF", c, st, cl, d.t, d.station_mean(st), Tsd, None, None, None, 5)
                else:
                    add("UF", c, st, cl, d.t, d.station_mean(st), Tsd, a["t"], a[f"T_{st}"],
                        np.interp(a["t"], d.t, Tsd), 5)
    # ---- F2 / Darr 1-g
    from .validate4 import DARR_BC, darr_data, f2_data, tsat

    f2 = f2_data()
    Ts_f2 = np.array([tsat(p * 1e3) for p in f2.P_exit_kPa])
    for cl in ("spec11", "specz11", "ds11", "d119s11", "db_vapor", "miropolskii"):
        a = _run("final11/val", f"f2__{cl}")
        add("F2 (0.05 g)", "F2", "exit", cl, f2.t.to_numpy(), f2.Twall_exit_avg_K.to_numpy(), Ts_f2,
            None if a is None else a["t"], None if a is None else a["T_exit"],
            None if a is None else np.interp(a["t"], f2.t, Ts_f2), 5)
        for test in (31, 40):
            dd = darr_data(test)
            Ts = tsat(DARR_BC[test]["P_out"])
            a = _run("final11/val", f"darr{test}__{cl}")
            for st, col in (("up", "Twall_TC1-3_avg_K"), ("dn", "Twall_TC4-6_avg_K")):
                add("Darr 1 g", test, st, cl, dd.t.to_numpy(), dd[col].to_numpy(), Ts,
                    None if a is None else a["t"], None if a is None else a[f"T_{st}"], Ts, 5)
    # ---- NBS
    from .validate4 import NBS_STATIONS, nbs_data, nbs_meta

    meta = nbs_meta().set_index("figure")
    nbs_runs = {
        10: {"spec11": "10:spec11_gh:Pout=83082,end=400", "miropolskii": "10:miropolskii:Pout=83082,end=400"},
        11: {"spec11": "11:spec11_gh:Pout=83082,end=384", "miropolskii": "11:miropolskii:Pout=83082,end=384"},
        12: {"spec11": "12:spec11_gh:Pout=83082,end=208", "miropolskii": "12:miropolskii:Pout=83082,end=208"},
        2: {"spec11": "2:spec11_gh:Pout=83082,end=128", "miropolskii": "2:miropolskii:Pout=83082,end=128"},
        3: {"spec11": "3:spec11_gh:Pout=83082,seg=15"}, 4: {"spec11": "4:spec11_gh:Pout=83082,seg=15,end=120"},
        5: {"spec11": "5:spec11_gh:Pout=83082,seg=15,end=104"}, 6: {"spec11": "6:spec11_gh:Pout=83082,seg=15,end=70"},
        7: {"spec11": "7:spec11_gh:Pout=83082,seg=15,end=45"},
    }
    nbs7 = {10: "10:spec_gh:Pout=83082,end=400", 11: "11:spec_gh:Pout=83082,end=384",
            12: "12:spec_gh:Pout=83082,end=208"}
    nbs7d = {11: "11:d119s_gh:Pout=83082,end=384", 12: "12:d119s_gh:Pout=83082,end=208"}
    for fig, runs in nbs_runs.items():
        fl = "Nitrogen" if meta.loc[fig].fluid == "LN2" else "Hydrogen"
        Ts = tsat(meta.loc[fig].P_drive_Pa, fl)
        dset = f"NBS {meta.loc[fig].fluid}"
        d = nbs_data(fig)
        allruns = {cl: ("nbs6", rid) for cl, rid in runs.items()}
        if fig in nbs7:
            allruns["spec_old"] = ("nbs7", nbs7[fig])
        if fig in nbs7d:
            allruns["d119s"] = ("nbs7", nbs7d[fig])
        for cl, (tag, rid) in allruns.items():
            a = _run(tag, rid)
            for st, x in NBS_STATIONS.items():
                dd = d[np.isclose(d.station_m, x, atol=0.3)]
                if dd.empty:
                    dd = d[np.isclose(d.station_m, 60.4, atol=0.3)]
                if dd.empty:
                    continue
                dd = dd.sort_values("time_s")
                add(dset, fig, st, cl, dd.time_s.to_numpy(), dd.T_wall_K.to_numpy(), Ts,
                    None if a is None else a["t"], None if a is None else a[f"T_{st}"], Ts, 3)
    df = pd.DataFrame(rows)
    df["err"] = df.tcd_m - df.tcd_d
    df["rel_err"] = df.err / df.tcd_d
    df["q_err"] = df.tq_m - df.tq_d
    return df


def summarize(df):
    out = []
    for (ds, cl), g in df[df.outcome != "no run"].groupby(["dataset", "closure"]):
        comp = g[g.tcd_d.notna()]
        hits = comp[comp.outcome == "hit"]
        n = len(comp)
        e = hits.err.abs()
        rec = dict(dataset=ds, closure=LABELS.get(cl, cl), n_stations=len(g), n_data_complete=n,
                   hit=int((g.outcome == "hit").sum()), miss=int((g.outcome == "miss").sum()),
                   false=int((g.outcome == "false").sum()),
                   MAE_s=e.mean() if len(e) else np.nan, bias_s=hits.err.mean() if len(e) else np.nan,
                   MAPE_pct=100 * hits.rel_err.abs().mean() if len(e) else np.nan,
                   within_1s=(e <= 1).sum() / n if n else np.nan, within_2s=(e <= 2).sum() / n if n else np.nan,
                   within_10pct=(hits.rel_err.abs() <= 0.1).sum() / n if n else np.nan)
        q = g[g.tq_d.notna() & g.tq_m.notna() & g.tcd_d.notna()]
        rec["quench_MAE_s"] = q.q_err.abs().mean() if len(q) else np.nan
        # whole section: last station per test
        ws = []
        for _, gt in g.groupby("test"):
            if gt.tcd_d.notna().all():
                td = gt.tcd_d.max()
                tm = gt.tcd_m.max() if gt.tcd_m.notna().all() else np.nan
                ws.append((tm - td) / td if np.isfinite(tm) else np.nan)
        ws = np.array(ws, float)
        rec["whole_n"] = len(ws)
        rec["whole_missed"] = int(np.isnan(ws).sum())
        rec["whole_MAPE_pct"] = 100 * np.nanmean(np.abs(ws)) if np.isfinite(ws).any() else np.nan
        rec["whole_bias_pct"] = 100 * np.nanmean(ws) if np.isfinite(ws).any() else np.nan
        out.append(rec)
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = station_rows()
    df.to_csv(OUT / "chilldown_stations.csv", index=False)
    sm = summarize(df)
    sm.to_csv(OUT / "chilldown_summary.csv", index=False)
    txt = ["# When is the line chilled down? (t_cd = wall reaches T_sat + 10 K)", "", __doc__, "",
           "## Summary per dataset and closure", "",
           sm.to_markdown(index=False, floatfmt=".2f"), "",
           "## Per station", "",
           df.assign(closure=df.closure.map(lambda c: LABELS.get(c, c))).to_markdown(index=False, floatfmt=".2f")]
    (OUT / "chilldown_time.md").write_text("\n".join(txt) + "\n")
    print(sm.round(2).to_string())


if __name__ == "__main__":
    main()
