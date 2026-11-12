\name{variance_comparison}
\alias{variance_comparison}
\title{Compare Conditional Variances for an Observed Design}

\usage{
variance_comparison(W, sigma2, rho, X, loc,
  tau_ns_2 = NULL, tau_sp_2 = NULL)
}

\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive residual variance shared by the two models.}
\item{rho}{Leroux spatial dependence parameter satisfying
\eqn{0 \leq \rho < 1}.}
\item{X}{Observation-level covariate vector or matrix, with one
covariate per column. No intercept is added.}
\item{loc}{Integer location indices from \code{1} to \code{nrow(W)},
one per row of \code{X}. Every location must be represented.
Within-location sample sizes may be unequal.}
\item{tau_ns_2,tau_sp_2}{Positive nonspatial variance and spatial
covariance scale. Supply either for variance matching, or both to
specify them separately. See \code{\link{variance_calibration}}.}
}

\description{
Compares Gaussian full conditional coefficient variances under
spatial and nonspatial models for an observed design, allowing
unequal replication and multiple covariates.
}

\details{
Covariance parameters and the remaining regression coefficients are
held fixed, with random intercepts integrated out and a flat prior
on the coefficient being compared. Each column of \code{X} is
evaluated separately. Covariates need not be centered or standardized.

The signed relative difference is the spatial variance divided by
the nonspatial variance, minus one. The absolute relative difference
is its absolute value.
}

\value{
A data frame with one row per coefficient, containing
\code{coefficient}, \code{variance_spatial}, \code{variance_iid},
\code{relative_signed}, and \code{relative_difference}.

The \code{calibration} attribute contains \code{g},
\code{tau_ns_2}, and \code{tau_sp_2}. The \code{counts} attribute
contains the observed sample size at each location.
}