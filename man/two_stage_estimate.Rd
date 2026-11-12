\name{two_stage_estimate}
\alias{two_stage_estimate}
\title{Two-Stage Gaussian Covariance Estimation}

\usage{
two_stage_estimate(y, X, loc, W, rho_cap = 0.9999)
}

\arguments{
\item{y}{Finite numeric response vector.}
\item{X}{Full-rank numeric design matrix with a first column of ones
and at least one predictor. No centering or standardization is performed.}
\item{loc}{Integer location indices from \code{1} to \code{nrow(W)},
one per observation. Every location must be represented.
Within-location sample sizes may be unequal.}
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal
representing a connected graph.}
\item{rho_cap}{Upper search limit for the spatial correlation parameter,
strictly between zero and one.}
}

\description{
Fits a nonspatial random-intercept model by REML, then estimates spatial
correlation using a Gaussian working likelihood for the predicted
random intercepts, accounting for shrinkage, observation noise,
and fixed-effect estimation. Also compares conditional coefficient
variances at the observed design.
}

\details{
The spatial covariance scale is calibrated so that its average marginal
random-intercept variance equals the estimated nonspatial variance.
The second stage forms the working likelihood using nonsingular
contrasts and estimates spatial correlation by numerical optimization.

The number of observations must exceed both the number of locations
and the number of regression coefficients.

If the estimated nonspatial random-intercept variance is zero,
\code{rho} is \code{NA}, both random-intercept variance parameters
are zero, and the conditional variances under the two models coincide.

Conditional comparisons use the estimated covariance parameters,
integrate out random intercepts, and condition on the remaining
regression coefficients.
}

\value{
A list containing:
\describe{
\item{\code{estimates}}{A data frame containing \code{rho},
\code{tau_ns_2}, \code{tau_sp_2}, \code{sigma2}, and
\code{upper_boundary}. The boundary indicator is \code{TRUE}
when the correlation estimate reaches \code{rho_cap}
within numerical tolerance.}
\item{\code{conditional_variances}}{A data frame with one row per
predictor, containing \code{coefficient}, \code{spatial_variance},
\code{iid_variance}, and \code{conditional_gap}.
The gap is the absolute relative variance difference using the
nonspatial variance as the denominator.}
\item{\code{iid_fit}}{The fitted nonspatial \code{lmer} model.}
\item{\code{stage2}}{A list containing \code{rho},
\code{upper_boundary}, \code{contrast_rank}, and \code{status}.
The contrast rank is the number of retained independent contrasts.
Status is \code{ok}, \code{rho_upper_cap}, or
\code{zero_variance_rho_unidentified}.}
}
}