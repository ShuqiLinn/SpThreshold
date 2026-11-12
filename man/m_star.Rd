\name{m_star}
\alias{m_star}
\title{Approximate and Exact Replication Thresholds}

\usage{
m_star(W, sigma2, rho, xbar = NULL, X = NULL,
  loc = NULL, gamma = 0.05, tau_ns_2 = NULL, tau_sp_2 = NULL,
  method = c("approx", "exact"), within_ss = NULL,
  max_search = 1e6)
}

\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive residual variance shared by the two models.}
\item{rho}{Leroux spatial dependence parameter satisfying
\eqn{0 \leq \rho < 1}.}
\item{xbar}{Location means for one covariate, or a matrix with one
covariate per column. Supply this or \code{X} with \code{loc}.}
\item{X}{Observation-level covariate vector or matrix for a balanced
design. No intercept is added.}
\item{loc}{Integer location indices from \code{1} to \code{nrow(W)},
one per row of \code{X}. Each location must have the same number
of observations.}
\item{gamma}{Positive tolerance for the absolute relative variance
difference, using the nonspatial variance as the denominator.}
\item{tau_ns_2,tau_sp_2}{Positive nonspatial variance and spatial
covariance scale. Supply either for variance matching, or both to
specify them separately. See \code{\link{variance_calibration}}.}
\item{method}{\code{"approx"} uses the leading-order approximation.
\code{"exact"} searches for the smallest integer threshold at least 2
such that the tolerance holds at that and every larger integer
sample size.}
\item{within_ss}{Within-location sum of squared covariate deviations
divided by the template's within-location sample size, with one value
per covariate. Calculated automatically from \code{X} and \code{loc}.
When only \code{xbar} is supplied, the default assumes centered
covariates with population variance one.}
\item{max_search}{Maximum integer sample size evaluated in the
exact search.}
}

\description{
Computes replication thresholds for agreement in Gaussian full
conditional coefficient variances as balanced within-location
replication increases.
}

\details{
Covariance parameters and the remaining regression coefficients are
held fixed, with random intercepts integrated out and a flat prior
on the coefficient being compared. Each covariate is evaluated
separately.

The location means and within-location variation per unit of
replication are held fixed as the sample size varies. Covariates
supplied through \code{X} need not be centered or standardized.

The approximate threshold uses the leading inverse-sample-size term
of the relative variance difference. The exact method checks every
integer sample size from 2 through an upper limit established by a
finite-sample bound. It returns one more than the last tolerance
violation, or 2 if there are none.

If the required search limit exceeds \code{max_search},
\code{m_exact} is \code{NA} and the status is \code{search_limit}.
The returned \code{sufficient_m} remains a sufficient replication
level established by the bound.

When within-location variation is zero, the function reports
\code{m_approx = Inf}, \code{m_exact = NA}, and status
\code{no_within_information}. If the two conditional variance
expressions are identical for the covariate, both thresholds are 2.
}

\value{
A list containing \code{m_star} for the selected method,
\code{m_approx}, \code{m_exact}, \code{gamma}, \code{method},
and a \code{thresholds} data frame with one row per covariate.

The table contains \code{coefficient}, \code{m_approx},
\code{m_exact}, \code{sufficient_m}, \code{status},
\code{searched_through}, \code{last_violation},
\code{leading_constant}, and \code{tail_constant}.

Additional components are \code{g}, \code{tau_ns_2},
\code{tau_sp_2}, the squared covariate projections \code{d_sq},
and the graph Laplacian \code{eigenvalues}.
}