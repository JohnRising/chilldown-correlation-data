"""Assemble the data/code repository for the chilldown-lowg paper (published through GitHub -> Zenodo).

    cd chilldown-lowg && .venv/bin/python analysis/build_zenodo_package.py

Writes ``zenodo-package/`` (deleted and rebuilt on each run) with the same directory layout as this
repository (``ofchill/``, ``analysis/``, ``data/``, ``results/``, ``paper/figures/``), so the scripts run
unchanged from inside the package. Only the run caches referenced by the result sets the paper uses are
copied. ``lit/`` (copyrighted source PDFs) is never included.
"""
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "zenodo-package"
RESULT_TAGS = ["final11", "final", "nbs6", "nbs7", "p4v"]
IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store", ".venv")


def copy_tree(src, dst, ignore=IGNORE):
    shutil.copytree(ROOT / src, PKG / dst, ignore=ignore)


def referenced_hashes(tag):
    hashes = set()
    for idx in (ROOT / "results" / tag).rglob("index.json"):
        def walk(x):
            if isinstance(x, dict):
                for v in x.values():
                    walk(v)
            elif isinstance(x, str) and len(x) >= 12 and (ROOT / "results" / "_cache" / f"{x}.json").exists():
                hashes.add(x)  # keys are a hex hash, sometimes with a suffix such as _seg15
        walk(json.load(open(idx)))
    return hashes


def main():
    if PKG.exists():
        shutil.rmtree(PKG)
    PKG.mkdir()
    copy_tree("ofchill", "ofchill")
    copy_tree("analysis", "analysis")
    copy_tree("data", "data")
    (PKG / "paper" / "figures").mkdir(parents=True)
    for f in sorted((ROOT / "paper" / "figures").glob("F*.pdf")):
        shutil.copy2(f, PKG / "paper" / "figures" / f.name)

    # results: aggregate tables, logs and run indexes of the sets the paper uses; caches only if referenced
    cache_dst = PKG / "results" / "_cache"
    cache_dst.mkdir(parents=True)
    hashes = set()
    for tag in RESULT_TAGS:
        shutil.copytree(ROOT / "results" / tag, PKG / "results" / tag, ignore=IGNORE)
        hashes |= referenced_hashes(tag)
    missing = []
    for h in sorted(hashes):
        for ext in ("npz", "json"):
            s = ROOT / "results" / "_cache" / f"{h}.{ext}"
            if s.exists():
                shutil.copy2(s, cache_dst / s.name)
            elif ext == "npz":
                missing.append(h)
    shutil.copy2(ROOT / "results" / "README.md", PKG / "results" / "README.md")
    print(f"result sets {RESULT_TAGS}: {len(hashes)} cached runs, {len(missing)} missing")

    # top-level documents are authored in zenodo-meta/ and copied verbatim
    for f in (ROOT / "zenodo-meta").iterdir():
        shutil.copy2(f, PKG / f.name)

    # checksum manifest (excluding itself)
    lines = []
    for f in sorted(PKG.rglob("*")):
        if f.is_file() and f.name != "MANIFEST.sha256":
            lines.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(PKG).as_posix()}")
    (PKG / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    size = sum(f.stat().st_size for f in PKG.rglob("*") if f.is_file())
    print(f"{len(lines)} files, {size / 1e6:.1f} MB -> {PKG}")


if __name__ == "__main__":
    main()
