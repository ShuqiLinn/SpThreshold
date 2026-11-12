# SpThreshold: Replication Threshold for Bayesian Multilevel Models of Areal Data

This package fits Bayesian multilevel models for areal data with Leroux spatial or independent Gaussian random intercepts. It supports Gaussian, Bernoulli, binomial, and negative binomial outcomes, with unequal replication and multiple covariates.

For Gaussian outcomes, the package compares full conditional coefficient variances under spatial and nonspatial models and computes approximate and exact replication thresholds for a specified tolerance. A two-stage procedure applies the criterion when covariance parameters are unknown, using a preliminary nonspatial fit.

## Installation

```r
devtools::install_github("ShuqiLinn/SpThreshold")
library(SpThreshold)
```

## Main functions

- `m_star()` — compute approximate (`method = "approx"`) or exact (`method = "exact"`) replication thresholds.
- `two_stage_estimate()` — estimate covariance parameters and compare conditional coefficient variances using the two-stage procedure.
- `variance_comparison()` — compare Gaussian full conditional coefficient variances at an observed design, including unequal replication and multiple covariates.
- `variance_calibration()` — match the average marginal random-intercept variance between the spatial and nonspatial models.
- `spfit()` — fit spatial (`model_indicator = 1`) or nonspatial (`model_indicator = 0`) models using MCMC.
- `posterior_draws()` and `posterior_summary()` — extract posterior draws and summarize means, medians, variances, and credible intervals.

## Examples and model details

Worked examples are provided in [SpThreshold_Example](SpThreshold_Example/). The model specifications and replication criterion are described in [SpThreshold_Model_Details](SpThreshold_Model_Details/).

## Reference

Lin, S. and Warren, J. L. (2026+). *On the Need for Spatial Random Effects in Bayesian Regression Models for Multilevel Areal Data.*
