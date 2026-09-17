\name{variance_calibration}
\alias{variance_calibration}
\title{Calibrate iid and Spatial Random-Effect Variances}
\usage{
variance_calibration(W, rho, tau0_2 = NULL, tau1_2 = NULL,
  tau2 = NULL, calibration = c("full", "common"))
}
\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal. The proper Leroux covariance also permits disconnected graphs.}
\item{rho}{Leroux parameter, between zero inclusive and one exclusive.}
\item{tau0_2}{Positive iid random-effect variance. Recommended input scale.}
\item{tau1_2}{Positive spatial random-effect scale.}
\item{tau2}{Deprecated alias for \code{tau1_2}; always the spatial scale.}
\item{calibration}{\code{"full"} uses the average marginal variance over all locations. \code{"common"} explicitly reproduces the legacy assumption of a common numerical variance scale.}
}
\description{Converts between variance scales using the full-field trace of the proper Leroux covariance matrix.}
\details{
Let \eqn{L=\operatorname{diag}(W\boldsymbol{1})-W}, with eigenvalues \eqn{\lambda_i}, and set
\deqn{c_i=1-\rho+\rho\lambda_i,\qquad
 g=\frac{1}{n}\sum_{i=1}^n c_i^{-1}.}
For \code{calibration="full"}, specifying only the iid scale sets
\eqn{\tau_1^2=\tau_0^2/g}; specifying only the spatial scale sets
\eqn{\tau_0^2=g\tau_1^2}. All \eqn{n} eigenvalues are used, including zero modes.
The average marginal variance under the specified spatial model is exactly
\eqn{\tau_1^2 g}. Using this value to approximate a fitted iid variance is a calibration assumption, not an exact identity for its posterior mean.
If both variances are supplied, they are used as specified without enforcing the full-field conversion. Their provenance is reported as \code{specified_separate_variances}.
The deprecated \code{tau2} argument retains its spatial-scale meaning; it no longer silently assumes that the iid variance is identical. Request \code{calibration="common"} to reproduce that comparison.
}
\value{A list with \code{g}, \code{tau0_2}, \code{tau1_2}, \code{calibration}, and graph quantities \code{n}, \code{W}, \code{rho}, \code{eigenvalues}, \code{eigenvectors}, \code{c_i}. Eigenvalues are sorted ascending. Eigenvectors are columns.}
\examples{
W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
variance_calibration(W, rho = 0.8, tau0_2 = 0.5)
}
\seealso{\code{\link{m_star}}, \code{\link{variance_curve}}, \code{\link{variance_comparison}}}
