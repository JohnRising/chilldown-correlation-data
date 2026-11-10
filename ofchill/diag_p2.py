"""1-D vs 0-D surrogate diagnostics for a phase-2 tag.

  .venv/bin/python -m ofchill.diag_p2 --tag p2 --label va  -> results/<tag>/diag_<label>.md

Per case-station: RMS 1-D vs 0-D surrogate (analysis/apriori/out/t06_surrogate.md, Pe-law + 134 K),
model Tsat at the station minus measured Tsat (median over the film period), model equilibrium
quality at the station (median, film period), axial-conduction flux relative to the convective
flux (median |q_ax|/|q_conv| in film boiling), and quench times (t130) up/dn and front speed.
"""

from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd

from .data import ROOT, UF_GEOM, load_uf_case
from .metrics import cross
from .props import fluid_props
from .runner import RESULTS, load_index, load_run


def ss_k(T):
    x = np.log10(T)
    c = [-1.4087, 1.3982, 0.2543, -0.626, 0.2334, 0.4256, -0.4658, 0.165, -0.0199]
    return 10 ** np.polynomial.polynomial.polyval(x, c)


def surrogate_table() -> pd.DataFrame:
    f = ROOT / "analysis" / "apriori" / "out" / "t06_surrogate.md"
    txt = f.read_text().split("## per station")[1]
    rows = []
    for line in txt.splitlines():
        m = re.match(r"\|\s*(\d+)\s*\|\s*(up|dn)\s*\|\s*([\d.]+)\s*\|.*", line)
        if m:
            cells = [c.strip() for c in line.strip("|").split("|")]
            rows.append(dict(case=int(cells[0]), station=cells[1], rms0d=float(cells[2]),
                             t130_0d=float(cells[7]) if cells[7] else np.nan))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="p2")
    ap.add_argument("--label", default="va")
    a = ap.parse_args(argv)
    fp = fluid_props("Nitrogen")
    idx = load_index(a.tag)
    g = UF_GEOM
    A_m = np.pi / 4 * (g["OD"] ** 2 - g["D"] ** 2)
    rows = []
    for rid, ent in sorted(idx.items(), key=lambda kv: kv[1]["spec"]["case"]):
        if ent["spec"]["label"] != a.label:
            continue
        c = ent["spec"]["case"]
        try:
            meta, ar = load_run(a.tag, rid)
        except FileNotFoundError:
            continue
        if ar is None:
            continue
        d = load_uf_case(c)
        t, wx, W = ar["t"], ar["wall_x"], ar["wall_T"]
        dx = wx[1] - wx[0]
        # axial conduction into each interior wall node per unit inner area
        k = ss_k(np.clip(W, 20, 400))
        qax = np.zeros_like(W)
        qax[1:-1] = k[1:-1] * A_m * (W[2:] - 2 * W[1:-1] + W[:-2]) / dx**2 / (np.pi * g["D"])
        tq = {}
        for st in ("up", "dn"):
            xs = g["stations"][st]
            j = int(np.argmin(np.abs(wx - xs)))
            Tst = ar[f"T_{st}"]
            Td = d.station_mean(st)
            film = (Tst > 150) & (t > 1.0)
            Tsat_m = fp.Tsat(ar["fluid_P"][j])
            Tsat_d = np.interp(t, d.t, d.col(f"Tsat_{st}_K"))
            xq = np.clip(np.nan_to_num(ar["fluid_x"][j], nan=np.where(ar["phase_code"][j] == 2, 1.0, 0.0)), 0, 1)
            qc = np.abs(ar["q"][j])
            ratio = np.abs(qax[j]) / np.maximum(qc, 1.0)
            tq[st] = (cross(t, Tst, 130), cross(d.t, Td, 130))
            rows.append(dict(case=c, station=st, rms1d=float(np.sqrt(np.mean((np.interp(d.t[d.t <= t[-1]], t, Tst)
                                                                         - Td[d.t <= t[-1]]) ** 2))),
                             dTsat_model_minus_data=float(np.median((Tsat_m - Tsat_d)[film])) if film.any() else np.nan,
                             x_station_med=float(np.median(xq[film])) if film.any() else np.nan,
                             x_station_max=float(np.max(xq)),
                             qax_over_qconv_film=float(np.median(ratio[film])) if film.any() else np.nan,
                             t130_m=tq[st][0], t130_d=tq[st][1]))
        dz = g["stations"]["dn"] - g["stations"]["up"]
        for r in rows[-2:]:
            if all(v is not None for v in (tq["up"][0], tq["dn"][0])):
                r["front_mm_s_m"] = 1e3 * dz / max(tq["dn"][0] - tq["up"][0], 1e-3)
            if all(v is not None for v in (tq["up"][1], tq["dn"][1])):
                r["front_mm_s_d"] = 1e3 * dz / max(tq["dn"][1] - tq["up"][1], 1e-3)
    df = pd.DataFrame(rows).merge(surrogate_table(), on=["case", "station"], how="left")
    out = RESULTS / a.tag / f"diag_{a.label}.md"
    out.write_text(df.to_markdown(index=False, floatfmt=".2f") + "\n")
    print(df.to_markdown(index=False, floatfmt=".2f"))


if __name__ == "__main__":
    main()
