[README.md](https://github.com/user-attachments/files/29797054/README.md)
# Network Concentration Index (NCI) — Replication Package

Replication code for:

> Riso, L. and Zoia, M.G. (2026). *Structural Concentration in Weighted
> Networks: A Class of Topology-Aware Indices.*

All results reported in the paper — simulation study and the three
empirical applications — are fully reproducible from the scripts in this
repository.

## Repository structure

| Script | Paper section | Reproduces |
|---|---|---|
| `NCI_simulation_study.R` | Section 4, Appendix C | Table 2, Figures 1–4; Monte Carlo validation of Proposition 3.4; scalability analysis (Table 6, Appendix C.5) |
| `production_network.R` | Section 5.1 | Table 3, Figures 5–7; sector-level WIOD production network (θ = 0.01, δ = 0.0461, ψ/δ ≈ 1.69) |
| `trade_network.R` | Section 5.2 | Table 4, Figures 8–10; country-level WIOD trade network (θ = 0.005, δ = 0.0782, ψ/δ ≈ 3.39) |
| `equity_network.R` | Section 6 | Table 5, Figures 11–13; S&P 500 top-20 MST dependence network, rolling-window NCI with bootstrap bands, median-threshold robustness check (ψ = 0.6854) |

## Overview

`NCI_simulation_study.R` implements the baseline Network Concentration
Index and its six variants (density-adjusted, null-model,
degree-constrained, weighted, transformed-data, multi-layer) and
reproduces all simulation results reported in the paper:

- **Deterministic scenarios** — three N = 10 networks (core-periphery,
  peripheral, random) with identical weight vector and comparable density,
  isolating the effect of topology on concentration (Table 2, Figures 1–2).
- **Monte Carlo experiments** — R = 5,000 replications per
  network-generating mechanism with fixed weights, and R = 800 joint draws
  with random weights from the simplex (Figures 3–4).
- **Verification of theoretical propositions** — numerical confirmation of
  Proposition 3.4 (E[ψ] = p under Erdős–Rényi) and of the null-model and
  density-adjusted benchmarks over the grid p ∈ {0.05, …, 0.95}.
- **Scalability analysis** — mean NCI and standard deviation for
  N ∈ {10, 50, 100} across all three mechanisms (Appendix C.5). In the
  Monte Carlo experiments, the degree-constrained benchmark is evaluated
  through the expected-adjacency approximation of the configuration model,
  B̄ᵢⱼ = dᵢdⱼ/(2m) (see Appendix C.4).

The three application scripts construct the empirical networks, compute
the NCI against the HHI and Gini benchmarks, produce all figures, and
include the threshold-sensitivity analyses (Figures 7 and 10) and the
robustness checks reported in the paper (median-threshold filtering and
minimum pairwise correlation for the equity network).

## Data availability

No raw data are redistributed in this repository.

- **WIOD (production and trade networks).** The 2016-release input–output
  tables (R format, `WIOT2014_October16_ROW.RData`) are freely available
  at <https://www.rug.nl/ggdc/valuechain/wiod/wiod-2016-release>.
  Download the "WIOT tables in R format" and set `WIOD_PATH` at the top
  of `production_network.R` / `trade_network.R` to the local file
  location (instructions in the script headers).
- **Equity data.** Daily adjusted closing prices for the S&P 500 top-20
  constituents (1 January 2015 – 31 December 2025) are downloaded
  automatically from Yahoo Finance via the `quantmod` package; no manual
  download is required.

## Requirements

R ≥ 4.2 with: `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `reshape2`,
`igraph`, `scales`, `data.table`, `ggrepel`, `quantmod`.

## Reproducibility

All replications are seeded (`set.seed(42)`). Numerical errors are
handled via `tryCatch` and excluded listwise. Running the full
simulation script reproduces every figure and table of the simulation
section, including the LaTeX source of the scalability table, printed to
console. The application scripts print the concentration indices, the
network densities, and the assortative-connectivity ratios reported in
Sections 5–6.

## Citation

If you use this code, please cite:

```bibtex
@article{riso2026structural,
  title   = {Structural Concentration in Weighted Networks:
             A Class of Topology-Aware Indices},
  author  = {Riso, Luigi and Zoia, Maria Grazia},
  journal = {arXiv preprint arXiv:2603.21918},
  year    = {2026}
}
```

## License

[choose: MIT / GPL-3 / CC-BY-4.0]

## Contact

Luigi Riso — luigi.riso@unicatt.it
Department of Economic Policy, Università Cattolica del Sacro Cuore, Milan
