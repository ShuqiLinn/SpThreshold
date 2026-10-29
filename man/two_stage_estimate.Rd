\name{two_stage_estimate}
\alias{two_stage_estimate}
\title{Two-Stage Gaussian Covariance Estimation}
\usage{
two_stage_estimate(y, X, loc, W, rho_cap = 0.9999)
}
\arguments{
\item{y}{Finite numeric response vector.}
\item{X}{Full-rank numeric design matrix with a first column of ones and at
least one predictor. No centering or standardization is performed.}
\item{loc}{Integer area indices aligned with the rows of \code{W}, one per
observation. Every area must be observed; counts may be unequal.}
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal and a
connected graph.}
\item{rho_cap}{Upper search limit for the spatial correlation, strictly
between zero and one.}
}
\description{
Fits a nonspatial random-intercept model by REML, then estimates spatial
correlation using a Gaussian working likelihood for the predicted random
intercepts, accounting for shrinkage, observation noise and fixed-effect
estimation. Includes conditional variance comparisons at the observed design.
}
\details{
The spatial scale is calibrated to match the average marginal random-effect
variance to the nonspatial variance. The second stage uses nonsingular
contrasts and a grid search followed by local optimization.
There must be more observations than areas and fixed-effect coefficients.
If REML gives zero random-effect variance, \code{rho} is \code{NA}, both
random-effect variances are zero, and both conditional variances coincide.
Small positive variances are not automatically set to zero.
The conditional comparisons integrate out random effects while holding
covariance parameters and all other regression coefficients fixed.
}
\value{
A list containing \code{estimates} (\code{rho}, \code{tau0_2} for the
nonspatial variance, \code{tau1_2} for the spatial scale, \code{sigma2},
and \code{upper_boundary}); \code{conditional_variances} for the predictors;
\code{iid_fit}, the fitted \code{lmer} model; and \code{stage2}, including
the contrast rank and estimation status.
}
\seealso{\code{\link{variance_comparison}}, \code{\link{m_star}}}
