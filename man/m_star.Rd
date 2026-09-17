\name{m_star}
\alias{m_star}
\title{Leading-Order and Exact Integer Replication Thresholds}
\usage{
m_star(W, sigma2, tau2 = NULL, rho, xbar = NULL, X = NULL,
  loc = NULL, gamma = 0.05, tau0_2 = NULL, tau1_2 = NULL,
  method = c("approx", "exact"), within_ss = NULL,
  calibration = c("full", "common"), max_search = 1e6)
}
\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive common residual variance.}
\item{tau2}{Deprecated spatial-scale alias. Prefer named \code{tau0_2} or \code{tau1_2}. The third positional argument retains its old spatial-scale meaning.}
\item{rho}{Leroux parameter in \eqn{[0,1)}.}
\item{xbar,X,loc,within_ss}{Covariate template inputs; see \code{\link{variance_curve}}. Balanced replication is required for the \code{X}, \code{loc} interface.}
\item{gamma}{Positive relative variance tolerance. A value equal to the tolerance is acceptable.}
\item{tau0_2,tau1_2,calibration}{Variance scales and conversion; see \code{\link{variance_calibration}}.}
\item{method}{\code{"approx"} returns the leading-order rule; \code{"exact"} seeks the smallest integer at least 2 after which the exact conditional curve stays within tolerance at every integer replication.}
\item{max_search}{Maximum integer replication evaluated in the exact search. If the rigorous tail bound lies beyond this limit, the exact result is unresolved, not a guessed threshold.}
}
\description{Computes replication guidance for fixed Gaussian covariance parameters and a balanced covariate template using separately specified or calibrated spatial and iid variances.}
\details{
The curve and target are defined in \code{\link{variance_curve}}. Let
\eqn{D_i=(U^T\bar{x})_i^2}, \eqn{c_i=1-\rho+\rho\lambda_i}, and
\eqn{s=W_{\mathrm{within}}/m>0}. The approximate threshold is
\deqn{m_{\mathrm{approx}}=\max\left\{2,\left\lceil
\frac{\sigma^2|\sum_iD_i(1/\tau_0^2-c_i/\tau_1^2)|}{\gamma s}
\right\rceil\right\}.}
With full-field calibration this becomes the whiteboard expression with
\eqn{(\sigma^2/\tau_0^2)\sum_iD_i(1-gc_i)} in the signed numerator.
It is a leading-order approximation, not a finite-sample guarantee.

A rigorous bound is \eqn{R(m)\leq A/m}, with
\deqn{A=\frac{\sigma^2\sum_iD_i|1/\tau_0^2-c_i/\tau_1^2|}{s}.}
Therefore the exact search need only inspect integers through
\eqn{\max\{2,\lceil A/\gamma\rceil\}}. It checks every integer up to this bound, including possible later violations after an initial crossing, and returns one more than the last violation (or 2 if none). All larger integers satisfy the tolerance by the bound. This definition concerns integer replication; crossings at fractional replication are not evaluated.

If \code{max_search} is insufficient, \code{m_exact} is \code{NA}, with status \code{search_limit}, the evaluated endpoint, and a certified sufficient \code{sufficient_m}. Increase the search limit to resolve the exact minimum. When \eqn{s=0}, the \eqn{1/m} bound does not apply: \code{m_approx=Inf}, \code{m_exact=NA}, and status \code{no_within_information}. This does not assert that the exact curve never falls below tolerance. Identical covariance action on the covariate, including calibrated \eqn{\rho=0}, instead gives exact threshold 2.

These are conditional coefficient variance comparisons under flat regression priors: covariance parameters and every other coefficient are fixed, with random effects integrated out. They are not guarantees for marginal posterior inference or coefficient bias. Multiple covariate columns receive separate conditional thresholds.
}
\value{A list with \code{m_star} (the selected method), \code{m_approx}, \code{m_exact}, \code{gamma}, \code{method}, and a \code{thresholds} data frame. The table reports the coefficient, both thresholds, \code{sufficient_m}, \code{status}, \code{searched_through}, \code{last_violation}, \code{leading_constant}, \code{tail_constant}, and \code{within_ss_per_m}. Other fields are \code{g}, \code{tau0_2}, \code{tau1_2}, \code{calibration}, \code{numerator_raw} (legacy unscaled spectral sum), \code{d_sq}, \code{eigenvalues}, \code{within_ss}, and \code{target}. \code{m_exact} is normally unevaluated for \code{method="approx"}.}
\examples{
W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
m_star(W, sigma2 = 0.5, rho = 0.8, tau0_2 = 0.5,
       xbar = c(-0.5, 0, 0.5), method = "exact")
}
\seealso{\code{\link{variance_calibration}}, \code{\link{variance_curve}}, \code{\link{variance_comparison}}}
