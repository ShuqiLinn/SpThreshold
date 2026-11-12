# Paper replication code

Code and results accompanying:

> Lin, S. and Warren, J. L. (2026+). On the Need for Spatial Random Effects in Bayesian Regression Models for Multilevel Areal Data.

## Contents

- `validation/empirical_validation.R` — numerical checks of the Gaussian full conditional precision expressions and approximate replication threshold in Section 3, using direct matrix calculations.
- `figures/simulation_results.rds` — aggregated simulation results used to produce the figures.
- `figures/plotting_code.R` — plotting code for the simulation figures in the main text and Supplementary Material.

## Reproducing the figures

From the `paper/` directory, run:

```r
setwd("figures")
source("plotting_code.R")
```

The script reads `simulation_results.rds` and saves the generated figures in the working directory.