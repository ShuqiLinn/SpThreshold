# SpThreshold: Replication Threshold for Bayesian Multilevel Models of Areal Data

This package fits Bayesian multilevel models for repeatedly measured areal spatial data using Leroux conditional autoregressive (CAR) or independent Gaussian random effects. It supports Gaussian, Bernoulli, binomial and negative binomial outcomes.

For Gaussian outcomes, the package evaluates how within-region replication changes the conditional posterior variances of regression coefficients. It provides a closed-form leading-order approximation to the replication threshold $m^*$ and a numerical threshold based on the exact conditional variance expression. Posterior means, variances and intervals are calculated directly from MCMC draws.

## Functions

- `m_star()` - replication threshold using the leading-order approximation (`method = "approx"`) or the exact conditional variance curve (`method = "exact"`)
- `variance_curve()` - exact and leading-order Gaussian relative variance differences across replication counts
- `variance_calibration()` - conversion between spatial and iid variance scales
- `variance_comparison()` - Gaussian conditional variance comparisons for a supplied design, including unequal replication and multiple covariates
- `spfit()` - MCMC posterior sampling for the spatial (`model_indicator = 1`) or nonspatial (`model_indicator = 0`) model
- `spfit_glmm()` - direct interface to the Bernoulli, binomial and negative binomial samplers
- `posterior_draws()`, `posterior_summary()` - extraction and summarization of posterior draws after removing the initial state and burn-in
- `create_random_W()`, `create_W_from_shapefile()` - adjacency matrix construction

The Gaussian variance calculations condition on specified covariance parameters and the other regression coefficients, with the random effects integrated out. They can be compared with marginal posterior variances obtained from MCMC.

Binary and count models use Pólya–Gamma augmentation through `pgdraw`. Negative binomial models require a fixed positive integer `nb_size`. The replication threshold formulas apply to Gaussian outcomes.

See the `SpThreshold_Example` folder for worked examples and `SpThreshold_Model_Details` for the underlying statistical model.

## Installation

```r
devtools::install_github("ShuqiLinn/SpThreshold")
```

## Reference

Lin, S. and Warren, J. L. (2026+). On the Need for Spatial Random Effects in Bayesian Regression Models for Multilevel Areal Data.