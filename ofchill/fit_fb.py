"""A-priori constants for the film-boiling options of the unified boiling closure.

Uses the analysis agent's vetted film-boiling dataset (analysis/apriori/fb_models.dataset:
energy-consistent h_ec, G lag-corrected by 2.5 s, 4422 samples / 20 case-stations) and the
same fit (least squares in ln h, every case-station weighted equally).

  .venv/bin/python -m ofchill.fit_fb       -> ofchill/fb_constants.json
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
APRIORI = HERE.parent / "analysis" / "apriori"
OUT = HERE / "fb_constants.json"


def _units(d, G, D):
    """unit (C = 1) forms of each family, both property sets."""
    u = G / d.rho_l
    out = {}
    for props in ("film", "sat"):
        if props == "film":
            k, rho, cp = d.k_vf, d.rho_vf, d.cp_vf
        else:
            k, rho, cp = d.k_v, d.rho_v, d.cp_v
        out[f"renewal_l_{props}"] = np.sqrt(k * rho * cp * u / D)
        out[f"renewal_v_{props}"] = np.sqrt(k * cp * G / D)
        hfgp = d.h_fg + 0.5 * cp * d.dT
        out[f"bromley_{props}"] = np.sqrt(k * rho * hfgp * u / (D * d.dT))
        out[f"pe_l_{props}"] = (k / D, u * D * rho * cp / k)  # (prefactor, Pe)
    return out


def main():
    sys.path.insert(0, str(APRIORI))
    cwd = os.getcwd()
    os.chdir(APRIORI)
    try:
        import fb_models as fm  # noqa: E402
        from common import D  # noqa: E402

        df = fm.dataset()
    finally:
        os.chdir(cwd)
    G = df.G.values
    h = df.h_ec.values
    n = df.groupby(["case", "st"]).t.transform("size").values
    sw = np.sqrt(1.0 / n)
    U = _units(df, G, D)
    res = {"source": "analysis/apriori fb_models.dataset (h_ec, G(t+2.5 s)), ln-LSQ, equal case-station weights",
           "n_samples": int(len(df))}

    def stats(pred):
        r = np.log(pred) - np.log(h)
        w = 1 / n
        return {"bias_gm": float(np.exp(np.sum(w * r) / np.sum(w)) - 1), "rms_ln": float(np.sqrt(np.sum(w * r**2) / np.sum(w)))}

    for key, val in U.items():
        if key.startswith("pe_l"):
            pre, Pe = val
            f = lambda p: sw * (np.log(p[0] * pre * Pe ** p[1]) - np.log(h))  # noqa: E731
            p = least_squares(f, [5.0, 0.35], x_scale="jac").x
            res[key] = {"C": float(p[0]), "n": float(p[1]), **stats(p[0] * pre * Pe ** p[1])}
        else:
            C = float(np.exp(np.sum(sw**2 * (np.log(h) - np.log(val))) / np.sum(sw**2)))
            res[key] = {"C": C, **stats(C * val)}
    hd = fm.darr2016(df, G).values
    C = float(np.exp(np.sum(sw**2 * (np.log(h) - np.log(hd))) / np.sum(sw**2)))
    res["darr2016_scale"] = {"C": C, **stats(C * hd), "published_stats": stats(hd)}
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
