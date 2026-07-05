# Network Concentration Index (NCI) — Simulation Study

Replication code for the simulation study (Section 4 and Appendix B) of:

> Riso, L. and Zoia, M.G. (2026). *Structural Concentration in Weighted 
> Networks: A Class of Topology-Aware Indices.* 

## Overview

The script `NCI_simulation_study.R` implements the baseline Network 
Concentration Index and its six variants (density-adjusted, null-model, 
degree-constrained, weighted, transformed-data, multi-layer) and reproduces 
all simulation results reported in the paper:

- **Deterministic scenarios** — three N = 10 networks (core-periphery, 
  peripheral, random) with identical weight vector and comparable density, 
  isolating the effect of topology on concentration (Table 2, Figures 1–2).
- **Monte Carlo experiments** — R = 5,000 replications per network-generating 
  mechanism with fixed weights, and R = 800 joint draws with random weights 
  from the simplex (Figures 3–4).
- **Verification of theoretical propositions** — numerical confirmation of 
  Proposition 3.4 (E[ψ] = p under Erdős–Rényi) and of the null-model and 
  density-adjusted benchmarks over the grid p ∈ {0.05, …, 0.95}.
- **Scalability analysis** — mean NCI and standard deviation for 
  N ∈ {10, 50, 100} across all three mechanisms (Appendix B.5).

## Requirements

R ≥ 4.2 with: `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `reshape2`, 
`igraph`, `scales`.

## Reproducibility

All replications are seeded (`set.seed(42)`). Numerical errors are handled 
via `tryCatch` and excluded listwise. Running the full script reproduces 
every figure and table of the simulation section, including the LaTeX 
source of Table 1 and of the scalability table, printed to console.

## Citation

If you use this code, please cite: 
@article{riso2026structural,
  title={Structural Concentration in Weighted Networks: A Class of Topology-Aware Indices},
  author={Riso, Luigi and Zoia, Maria Grazia},
  journal={arXiv preprint arXiv:2603.21918},
  year={2026}
}

## License

[choose: MIT / GPL-3 / CC-BY-4.0]

## Contact

Luigi Riso — luigi.riso@unicatt.it  
Department of Economic Policy, Università Cattolica del Sacro Cuore, Milan
