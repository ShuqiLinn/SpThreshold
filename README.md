# SpThreshold: Replication Threshold for Bayesian Multilevel Models of Areal Data

This package implements a hierarchical Bayesian multilevel model for areal spatial data and computes a closed-form replication threshold beyond which a spatially structured model yields effectively equivalent inference to a simpler non-spatial alternative. The method relies on the conditional autoregressive model and is fit using Markov chain Monte Carlo sampling techniques. Please see the "SpThreshold_Model_Details" and "SpThreshold_Example" folders for more specific information regarding the statistical model and package use details, respectively.

## Installation
```r
devtools::install_github("ShuqiLinn/SpThreshold")
```

## Reference
Lin, S. and Warren, J. L. (2026+). On the Need for Spatial Random Effects in Bayesian Regression Models for Multilevel Areal Data.
