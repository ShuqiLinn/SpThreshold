# SpThreshold: Replication Threshold for Bayesian Multilevel Models of Areal Data

This package supports the analysis of repeatedly measured areal spatial data using a Bayesian hierarchical multilevel model. It computes a closed-form replication threshold $m^*$ beyond which a spatially structured model yields effectively equivalent inference to a simpler nonspatial alternative, and provides Markov chain Monte Carlo fitting samplers for both the spatial (Leroux conditional autoregressive) and nonspatial versions of the model.

## Functions

- `m_star()` - closed-form replication threshold derived in the manuscript
- `spfit()` - MCMC posterior sampling for the spatial (`model_indicator = 1`) or nonspatial (`model_indicator = 0`) model
- `create_random_W()`, `create_W_from_shapefile()` - adjacency matrix construction

See the `SpThreshold_Example` folder for worked examples and `SpThreshold_Model_Details` for the underlying statistical model.

## Installation
```r
devtools::install_github("ShuqiLinn/SpThreshold")
```

## Reference
Lin, S. and Warren, J. L. (2026+). On the Need for Spatial Random Effects in Bayesian Regression Models for Multilevel Areal Data.
