"""Chatter / numerics diagnostics per run: regime switches (h crossing 1500 W/m2K) per conductor,
re-dry events, max wall-T rise between steps, quench time spread along the tube."""
import sys
import numpy as np
from .runner import load_index, load_run


def chatter(ar, thr=1500.0):
    from .props import fluid_props

    h = ar["h"][:, 1:]
    hi = h > thr
    dT = (ar["wall_T"] - fluid_props("Nitrogen").Tsat(ar["fluid_P"]))[:, 1:]
    chg = np.abs(np.diff(hi.astype(int), axis=1)) * (dT[:, 1:] > 10.0)  # switches with Tw - Tsat > 10 K
    sw = chg.sum(axis=1)
    rise = float(np.max(np.diff(ar["wall_T"], axis=1)))
    t = ar["t"][1:]
    tq = [t[np.argmax(r)] if r.any() else np.nan for r in hi]
    return dict(max_switches=int(sw.max()), redry_conductors=int(np.sum(sw > 1)), max_wallT_rise=rise,
                tq_first=float(np.nanmin(tq)), tq_last=float(np.nanmax(tq)))


if __name__ == "__main__":
    tag = sys.argv[1]
    for rid, ent in sorted(load_index(tag).items()):
        try:
            m, ar = load_run(tag, rid)
        except FileNotFoundError:
            continue
        if ar is None:
            print(rid, m.get("status")); continue
        print(rid, f"solve={m['solve_s']:.0f}s", chatter(ar))
