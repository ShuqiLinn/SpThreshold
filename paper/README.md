# Paper replication code

Code accompanying:

> Lin, S. and Warren, J. L. (2026+). On the Need for Spatial Random Effects
> in Bayesian Regression Models for Multilevel Areal Data.

## Contents

`validation/empirical_validation.R` — numerical verification that the
closed-form precision expressions (Section 3.1) and the asymptotic threshold
(Section 3.2) agree with direct matrix-algebra computations of the same
quantities. Run line-by-line in R.

`figures/plotting_code.R` — script that generates the simulation figures in
Sections 4 and the supplement. Reads `figures/simulation_results.rds`
(aggregated simulation output) and writes six PNG files into the working
directory. Run from inside `figures/`:

```r
setwd("figures")
source("plotting_code.R")
```

The methods themselves are implemented in the parent package. See
`?SpThreshold::m_star` for the closed-form bound and `?SpThreshold::spfit`
for fitting the spatial and nonspatial multilevel models via MCMC.
