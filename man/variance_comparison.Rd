\name{variance_comparison}
\alias{variance_comparison}
\title{Conditional Gaussian Variance Comparison for an Observed Design}
\usage{
variance_comparison(W, sigma2, rho, X, loc,
  tau0_2 = NULL, tau1_2 = NULL, tau2 = NULL,
  calibration = c("full", "common"))
}
\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive common residual variance.}
\item{rho}{Leroux parameter in \eqn{[0,1)}.}
\item{X}{Numeric covariate vector or matrix; each column specifies one coefficient's conditional target. No intercept is added and no standardization is applied.}
\item{loc}{Integer location index per observation. Every location \code{1:nrow(W)} must be represented; counts may be unequal.}
\item{tau0_2,tau1_2,tau2,calibration}{See \code{\link{variance_calibration}}.}
}
\description{Computes exact conditional coefficient variances at the observed Gaussian design, including unequal replication and multiple covariates.}
\details{
For each covariate \eqn{x}, the conditional precision is
\deqn{I=\frac{\sum_{ij}(x_{ij}-\bar{x}_i)^2}{\sigma^2}
+\bar{x}^T\{C+\operatorname{diag}(\sigma^2/m_i)\}^{-1}\bar{x},}
where \eqn{C=\tau_1^2 Q(\rho)^{-1}} for the spatial model and
\eqn{C=\tau_0^2I} for the iid model. This equals \eqn{x^TV^{-1}x}, with
\eqn{V=\sigma^2 I+ZCZ^T}, while requiring only a location-sized inverse.

All other regression coefficients and covariance parameters are held fixed and random effects are integrated out. The conditional variance is \eqn{1/I} under a flat coefficient prior. For multiple covariates, the output is the reciprocal of each diagonal precision entry, not the diagonal of the inverse joint precision. It does not claim a marginal posterior variance result.

This function evaluates the observed unequal design exactly. It does not replace its counts by a minimum or mean and does not assert a replication threshold for an unspecified future unequal design.
}
\value{A data frame with \code{coefficient}, \code{variance_spatial}, \code{variance_iid}, \code{relative_signed}, \code{relative_difference}, and \code{within_ss}. Relative differences use the iid variance as denominator. Attributes contain \code{calibration}, observed \code{counts}, and the \code{target}.}
\examples{
W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
loc <- rep(1:3, c(2, 3, 4))
X <- cbind(x = seq_along(loc), z = sin(seq_along(loc)))
variance_comparison(W, sigma2 = 0.5, rho = 0.8, X = X,
                    loc = loc, tau0_2 = 0.5)
}
\seealso{\code{\link{variance_curve}}, \code{\link{m_star}}}
