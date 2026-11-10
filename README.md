# Data and code for "A film-boiling and rewetting correlation for one-dimensional simulation of cryogenic line chilldown"

John Rising, Independent researcher, Culver City, CA, United States.
Paper submitted to the *International Journal of Heat and Mass Transfer*.
Contact: johnrising@alum.mit.edu

This archive holds what is needed to check and rerun the numbers and figures in the paper: the simulation
and analysis code, the experimental time series used (some digitized from published figures; see
"Third-party data" below), and the simulation results from which the paper's tables are taken.

## Contents

The layout is the same as the working directory of the study, so the scripts run unchanged from the root of this
repository.

| Path | Content |
|---|---|
| `ofchill/` | Pipeline that builds OpenFLUME networks for each data set, applies the correlation (`closures.py`, `closure.py`), runs them, and computes the metrics. `NOTES.md` lists solver behavior that matters for the results; `tests/` tests the closure expression. |
| `analysis/make_paper_figs.py` | Generates the paper figures F1-F10 from `analysis/apriori/out/` and `results/`. |
| `analysis/fluid_ratio.py` | Vapor-sensible-to-latent ratio `S`, liquid-supply ratio, and figure F11; writes `analysis/out_fluid/`. |
| `analysis/apriori/` | Analysis of the low-gravity data before simulation: energy-consistent heat flux, film-boiling dataset, candidate closures, rewet-temperature theory, cross-validation, parsimony ladder (`a01`-`a19`). `REPORT.md` describes each script; `out/` holds its tables.  `a19_parsimony_report.py` writes the parsimony table (`out/parsimony.csv`) that figure F10 is drawn from. |
| `analysis/digitize_*.py` | Scripts used to digitize the curves in `data/hartwig2022_F2/` and `data/darr2016ijhmt_1g_vertical/`. |
| `data/` | Measured wall-temperature time series and boundary conditions (see below). |
| `results/` | Simulation results. `final11/`, `final/`, `nbs6/`, `nbs7/`, `p4v/` are the result sets the paper's numbers come from. `<set>/index.json` maps each run to a file in `results/_cache/` (histories as `.npz`, metadata as `.json`). `results/README.md` describes the layout. Only the runs referenced by those five sets are included. A few run ids in the `index.json` files have no result file (runs that failed, timed out or were not kept); they are missing from the working directory as well. |
| `paper/figures/` | The figures as they appear in the paper (PDF). |

## Software

Python 3.13; OpenFLUME 0.3.1 (`pip install openflume==0.3.1`, MIT license, <https://github.com/OpenFLUME/OpenFLUME>),
CoolProp 8.0.0, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.1, Matplotlib 3.11.2.
Every simulation in the paper used the unmodified OpenFLUME package; the correlation is supplied to it as an
expression (`ofchill/closures.py`).

## Reproducing

From the archive root, with those packages installed:

```bash
# figures from the stored results (about a minute; figures are written to paper/figures/)
python analysis/make_paper_figs.py F1 F2 F3 F4 F5 F6 F7 F8 F9 F10
python analysis/fluid_ratio.py                 # F11 and analysis/out_fluid/

# the a-priori analysis (about 4 minutes; see analysis/apriori/REPORT.md)
cd analysis/apriori && python a01_characterize.py && python a03_film_boiling.py   # etc.

# rerun the simulations of the final closure (each run is cached in results/_cache/, so identical
# configurations are not repeated; use --force to recompute)
python -m ofchill.final11 run --part uf sens num val
python -m ofchill.final11 report
python -m ofchill.final run --part uf sens num nbs
python -m ofchill.final report
```

The cache key of a run is a hash of its network configuration, the result-format version and the OpenFLUME
version, so results obtained with a different OpenFLUME version are recomputed rather than mixed with these.
Floating-point differences between platforms can move individual numbers in the last digits.

## Third-party data

The experiments were not done by the author. The files below are the measurements of the original authors,
either as distributed with their paper or digitized by the author from published figures. They are included so
that the paper can be checked, and they remain the work of the original authors: cite the original source when
using them. The license in `LICENSE-DATA.md` covers only the author's own files.

| Path | Source | How obtained |
|---|---|---|
| `data/uf2016/` | Darr, S. R. et al., *npj Microgravity* 2, 16033 (2016) (ten low-gravity LN2 tests, 304 stainless-steel tube) | Supplementary data file of the paper (MOESM26), columns renamed; `meta.json` gives the mapping |
| `data/hartwig2022_F2/` | Hartwig, J. et al., *npj Microgravity* 8, 33 (2022), doi:10.1038/s41526-022-00220-9 (flight test F2) | Digitized from Figs. 1a, 2a-c; Table 1 transcribed. Estimated digitization error is in the folder README |
| `data/darr2016ijhmt_1g_vertical/` | Darr, S. R. et al., *Int. J. Heat Mass Transf.* 103, 1225-1242 (2016), tests 31 and 40 (1 g upflow) | Digitized from the replot of the data in Kassemi, M. et al., NASA NTRS 20240013323 (2024); boundary conditions from the same report |
| `data/nbs9264/` | Brennan, J. A. et al., *Cooldown of Cryogenic Transfer Lines*, NBS Report 9264 (1966), NASA-CR-81338 | Digitized from the report's wall-temperature figures; uncertainty in the file headers (unreliable below 70 K) |

## Citation

Cite the paper (once published) and this archive (DOI assigned by Zenodo). `CITATION.cff` has the metadata.
