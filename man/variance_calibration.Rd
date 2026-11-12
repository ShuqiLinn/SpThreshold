\name{variance_calibration}
\alias{variance_calibration}
\title{Calibrate Spatial and Nonspatial Variances}

\usage{
variance_calibration(W, rho, tau_ns_2 = NULL, tau_sp_2 = NULL)
}

\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.
Disconnected graphs are allowed.}
\item{rho}{Leroux spatial dependence parameter satisfying
\eqn{0 \leq \rho < 1}.}
\item{tau_ns_2}{Positive nonspatial random-intercept variance.}
\item{tau_sp_2}{Positive spatial covariance scale.}
}

\description{
Converts between the nonspatial random-intercept variance and the spatial
covariance scale by matching average marginal variances.
}

\details{
Let \eqn{\lambda_i} denote the eigenvalues of the graph Laplacian. Define
\deqn{
c_i = 1 - \rho + \rho\lambda_i,
\qquad
g = \frac{1}{n}\sum_{i=1}^{n}c_i^{-1}.
}
Variance matching gives
\deqn{
\tau_{\mathrm{ns}}^2 = g\tau_{\mathrm{sp}}^2.
}
Supply at least one variance parameter. If only one is supplied, the other
is calculated using this relationship. If both are supplied, they are
used as specified without enforcing variance matching.
}

\value{
A list containing \code{tau_ns_2}, \code{tau_sp_2}, the calibration factor
\code{g}, and graph quantities \code{n}, \code{W}, \code{rho},
\code{eigenvalues}, \code{eigenvectors}, and \code{c_i}.
}