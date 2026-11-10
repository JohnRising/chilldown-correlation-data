"""Per-case figures with a chosen subset of labels: python -m ofchill.keyfigs TAG OUTDIR label1 label2 ..."""
import sys
from .plotting import plot_case
from .runner import RESULTS, load_index, load_run

if __name__ == "__main__":
    tag, outdir, labels = sys.argv[1], sys.argv[2], sys.argv[3:]
    by = {}
    for rid, e in load_index(tag).items():
        lab = e["spec"]["label"]
        if lab in labels:
            try:
                _, a = load_run(tag, rid)
            except FileNotFoundError:
                continue
            by.setdefault(e["spec"]["case"], []).append((labels.index(lab), lab, a))
    for c, cur in by.items():
        cur.sort(key=lambda r: r[0])
        plot_case(c, [(l, a) for _, l, a in cur], RESULTS / tag / outdir / f"case{c}.png", title=f"UF case {c} ({tag})")
    print(sorted(by))
