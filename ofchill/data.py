"""Experimental data loading (UF 2016 low-g LN2 chilldown, Darr et al. 2016)."""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
UF_DIR = DATA_DIR / "uf2016"

# Darr et al. 2016 (Methods) test-section geometry.
UF_GEOM = dict(
    L=0.572,
    D=0.01168,
    OD=0.0127,
    rho_wall=8027.0,
    roughness=1.5e-6,  # drawn tube (friction only), as in the reference builder
    stations={"up": 0.149, "dn": 0.401},
)
UF_STATION_TCS = {"up": ["TC1_K", "TC2_K", "TC3_K"], "dn": ["TC4_K", "TC5_K", "TC6_K"]}


@dataclass
class CaseData:
    case_no: int
    df: pd.DataFrame

    @property
    def t(self) -> np.ndarray:
        return self.df["t_s"].to_numpy()

    def col(self, name: str) -> np.ndarray:
        return self.df[name].to_numpy()

    def station_mean(self, st: str) -> np.ndarray:
        return self.df[UF_STATION_TCS[st]].mean(axis=1).to_numpy()

    def station_band(self, st: str) -> tuple[np.ndarray, np.ndarray]:
        d = self.df[UF_STATION_TCS[st]]
        return d.min(axis=1).to_numpy(), d.max(axis=1).to_numpy()

    @property
    def t_end(self) -> float:
        return float(self.t[-1])


@functools.lru_cache(maxsize=None)
def load_uf_case(case_no: int) -> CaseData:
    df = pd.read_csv(UF_DIR / f"case{case_no}.csv")
    df = df.apply(pd.to_numeric, errors="coerce").dropna(subset=["t_s"]).reset_index(drop=True)
    return CaseData(case_no, df)


def uf_meta() -> dict:
    return json.loads((UF_DIR / "meta.json").read_text())


def boxcar(t: np.ndarray, y: np.ndarray, tau: float, half: float, tol: float = 0.0) -> float:
    """Mean of samples with tc - half - tol <= t <= tc + half + tol, tc = tau
    clipped to the record (tol = 0: identical to the reference TypeScript
    builder); nearest sample if empty.

    NOTE: at 15 Hz the window edges tc +- 0.1 s fall exactly on sample times for
    every other knot, so membership depends on the last digit of t.  The CSV
    times are rounded (0.2) while the original xls has 0.19999999; tol = 1e-6
    makes the membership insensitive to that rounding."""
    tc = min(max(tau, t[0]), t[-1])
    m = (t >= tc - half - tol) & (t <= tc + half + tol)
    if m.any():
        return float(y[m].mean())
    return float(y[np.argmin(np.abs(t - tc))])


def boxcar_schedule(
    t: np.ndarray, y: np.ndarray, t_end: float, knot: float = 0.1, half: float = 0.1,
    shift: float = 0.0, scale: float = 1.0, tol: float = 0.0,
) -> list[list[float]]:
    """[[tk, scale * boxcar(y, tk + shift)], ...] at tk = 0, knot, ... t_end."""
    nk = int(round(t_end / knot)) + 1
    return [[k * knot, boxcar(t, y, k * knot + shift, half, tol) * scale] for k in range(nk)]
