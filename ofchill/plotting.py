"""Figures: per-case wall T vs time (data mean + TC min/max band + model
curves) and a summary figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .data import UF_GEOM, load_uf_case  # noqa: E402

STN_LABEL = {"up": "upstream (x = 0.149 m, TC1-3)", "dn": "downstream (x = 0.401 m, TC4-6)"}


def plot_case(case_no: int, curves: list[tuple[str, dict]], out: Path, title: str | None = None) -> Path:
    """curves: [(label, arrays from runner.load_run)]."""
    d = load_uf_case(case_no)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, st in zip(axs, ("up", "dn")):
        lo, hi = d.station_band(st)
        ax.fill_between(d.t, lo, hi, color="0.75", alpha=0.6, lw=0, label="data TC min/max")
        ax.plot(d.t, d.station_mean(st), "k-", lw=1.8, label="data mean")
        Ts = d.col(f"Tsat_{st}_K")
        ax.plot(d.t, Ts, "k:", lw=0.8, label="T_sat (data)")
        for i, (lab, a) in enumerate(curves):
            if a is None or f"T_{st}" not in a:
                continue
            ax.plot(a["t"], a[f"T_{st}"], lw=1.4, color=f"C{i}", label=lab)
        ax.set_title(STN_LABEL[st], fontsize=10)
        ax.set_xlabel("t [s]")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("outer wall T [K]")
    axs[1].legend(fontsize=8, loc="upper right")
    fig.suptitle(title or f"UF 2016 case {case_no}", fontsize=11)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_summary(rows: list[dict], out: Path) -> Path:
    """rows: metric rows (case, label, station, rms, t_130_m, t_130_d, ...)."""
    rows = [r for r in rows if r.get("status") == "ok"]
    labels = sorted({r["label"] for r in rows})
    cases = sorted({r["case"] for r in rows})
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    w = 0.8 / max(len(labels), 1)
    for i, lab in enumerate(labels):
        for k, st in enumerate(("up", "dn")):
            rr = {r["case"]: r for r in rows if r["label"] == lab and r["station"] == st}
            y = [rr[c]["rms"] if c in rr else np.nan for c in cases]
            axs[0].bar(np.arange(len(cases)) + i * w + (k - 0.5) * w / 2, y, w / 2, color=f"C{i}",
                       alpha=1.0 if st == "up" else 0.55, label=f"{lab} ({st})")
    axs[0].set_xticks(np.arange(len(cases)) + 0.4 - w / 2)
    axs[0].set_xticklabels([str(c) for c in cases])
    axs[0].set_xlabel("case")
    axs[0].set_ylabel("RMS wall-T error [K]")
    axs[0].legend(fontsize=7)
    for ax, key, name in ((axs[1], "t_130", "t(130 K)"), (axs[2], "t_quench", "t(max cooling rate)")):
        mx = 0
        for i, lab in enumerate(labels):
            for st, mk in (("up", "o"), ("dn", "s")):
                rr = [r for r in rows if r["label"] == lab and r["station"] == st]
                xd = [r[f"{key}_d"] for r in rr]
                ym = [r[f"{key}_m"] if r[f"{key}_m"] is not None else np.nan for r in rr]
                xd = [np.nan if v is None else v for v in xd]
                ax.plot(xd, ym, mk, color=f"C{i}", ms=6, label=f"{lab} {st}")
                mx = max(mx, np.nanmax(xd + ym + [0]))
        ax.plot([0, mx * 1.05], [0, mx * 1.05], "k--", lw=0.8)
        ax.set_xlabel(f"data {name} [s]")
        ax.set_ylabel(f"model {name} [s] (missing = never reached)")
        ax.grid(alpha=0.3)
    axs[1].legend(fontsize=7)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out
