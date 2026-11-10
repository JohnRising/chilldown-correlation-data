"""Quick NBS comparison plots: python -m ofchill.nbs6_plot OUT.png fig run1 run2 ..."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from .runner import load_run
from .validate4 import NBS_STATIONS, nbs_data

out, fig_no, runs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
d = nbs_data(fig_no)
f, axs = plt.subplots(1, 4, figsize=(16, 3.6), sharey=True)
for ax, (st, x) in zip(axs, NBS_STATIONS.items()):
    dd = d[np.isclose(d.station_m, x, atol=0.3)]
    if dd.empty:
        dd = d[np.isclose(d.station_m, 60.4, atol=0.3)]
    dd = dd.sort_values("time_s")
    ax.plot(dd.time_s, dd.T_wall_K, "k-", lw=2, label="data")
    for i, r in enumerate(runs):
        tag, rid = r.split("|")
        try:
            m, a = load_run(tag, rid)
            ax.plot(a["t"], a[f"T_{st}"], color=f"C{i}", label=rid[:40])
        except Exception as e:
            print(r, e)
    ax.set_title(f"{st} {x} m")
axs[-1].legend(fontsize=6)
f.tight_layout()
f.savefig(out, dpi=90)
