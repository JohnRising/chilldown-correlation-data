# chilldown-lowg-data

Simulation results for the paper "A film-boiling and rewetting correlation
for one-dimensional simulation of cryogenic line chilldown" (John Rising).

The code that produced them is in the private repository
OpenFLUME/OpenFLUME-papers, under `chilldown-lowg/` (`ofchill/` and
`analysis/`). `ofchill/runner.py` expects this repository to be cloned at
`chilldown-lowg/results/`:

    cd OpenFLUME-papers/chilldown-lowg
    git clone https://github.com/OpenFLUME/chilldown-lowg-data.git results

## Layout

- `<tag>/index.json`: maps each run id to the hash of its full OpenFLUME
  network configuration and to its spec (case, closure, N, dt, options).
- `_cache/<hash>.npz`, `_cache/<hash>.json`: the result of each run
  (histories as compressed NumPy arrays, and run metadata). The hash covers
  the configuration, the result format version and the OpenFLUME version
  (0.3.1).
- `<tag>/*.csv`, `<tag>/*.md`: aggregate and per-station tables. The
  numbers in the paper are taken from `final11/`, `final/`, `nbs6/`,
  `nbs7/` and `p4v/` (see `chilldown-lowg/paper/AUDIT.md`).
- `*.log`, `<tag>/*.log`: run logs.

The paper's data statement offers these "from the author on request".
