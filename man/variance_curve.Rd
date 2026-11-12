\name{variance_curve}
\alias{variance_curve}
\title{Conditional Variance Comparisons across Replication Levels}

\usage{
variance_curve(W, sigma2, rho, m, xbar = NULL, X = NULL,
  loc = NULL, tau_ns_2 = NULL, tau_sp_2 = NULL,
  within_ss = NULL)
}

\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive residual variance shared by the two models.}
\item{rho}{Leroux spatial dependence parameter satisfying
\eqn{0 \leq \rho < 1}.}
\item{m}{Positive replication values at which to evaluate the curves.}
\item{xbar}{Location means for one covariate, or a matrix with one
covariate per column. Supply this or \code{X} with \code{loc}.}
\item{X}{Observation-level covariate vector or matrix for a balanced
design. No intercept is added.}
\item{loc}{Integer location indices from \code{1} to \code{nrow(W)},
one per row of \code{X}. Each location must have the same number
of observations in the input design.}
\item{tau_ns_2,tau_sp_2}{Positive nonspatial variance and spatial
covariance scale. Supply either for variance matching, or both to
specify them separately. See \code{\link{variance_calibration}}.}
\item{within_ss}{Within-location sum of squared covariate deviations
divided by the input design's within-location sample size, with one
value per covariate. Calculated automatically from \code{X} and
\code{loc}. See Details for the default when supplying only \code{xbar}.}
}

\description{
Evaluates Gaussian full conditional variances under spatial and
nonspatial models, their relative differences, and the corresponding
leading-order approximation and finite-sample bound as balanced
within-location replication varies.
}

\details{
Covariance parameters and the remaining regression coefficients are
held fixed, with random intercepts integrated out and a flat prior
on the coefficient being compared. Each covariate is evaluated
separately.

The location means and within-location variation per unit of
replication are held fixed as \code{m} varies. Covariates supplied
through \code{X} need not be centered or standardized.

When only \code{xbar} is supplied and \code{within_ss} is omitted,
centered covariates with population variance one are assumed.
For each covariate, this gives
\deqn{
\mathrm{within\_ss} = n - \sum_{i=1}^{n}\bar{x}_i^2.
}
For another scale, supply \code{within_ss} explicitly.

Relative differences use the nonspatial variance as the denominator.
The leading-order approximation and upper bound are proportional
to the inverse of \code{m}.

When within-location variation is zero, exact variances remain
available, while the leading approximation is \code{NA} and the
bound is \code{Inf}. If the two conditional variance expressions
are identical for the covariate, the differences and bound are zero.
}

\value{
A data frame with one row per replication value and coefficient,
containing \code{m}, \code{coefficient}, \code{variance_spatial},
\code{variance_iid}, \code{relative_signed},
\code{relative_difference}, \code{leading_signed},
\code{leading_difference}, and \code{tail_bound}.

The \code{relative_difference} and \code{leading_difference}
columns contain absolute relative differences.
The \code{calibration} attribute contains \code{g},
\code{tau_ns_2}, and \code{tau_sp_2}.
}