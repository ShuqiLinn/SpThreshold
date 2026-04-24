\name{m_star}
\alias{m_star}
\title{Replication Threshold for Spatial vs Nonspatial Inference on a Regression Coefficient}
\description{
Computes the asymptotic replication threshold \eqn{m^*} beyond which a nonspatial model yields effectively equivalent posterior inference for a regression coefficient compared to a Leroux CAR spatial model, based on the closed-form bound derived in the manuscript.
}
\usage{
m_star(W, sigma2, tau2, rho, xbar = NULL, X = NULL, loc = NULL, gamma = 0.05)
}
\arguments{
\item{W}{Adjacency matrix, an \eqn{n} by \eqn{n} symmetric matrix of 0/1 entries, where \eqn{w_{ii'} = 1} indicates that areas \eqn{i} and \eqn{i'} are neighbors. Must correspond to a connected graph.}

\item{sigma2}{Within-area residual variance.}

\item{tau2}{Spatial random-effect variance.}

\item{rho}{Spatial correlation parameter in \eqn{[0, 1)}.}

\item{xbar}{Optional length-\eqn{n} vector of area-level covariate means. If supplied, \code{X} and \code{loc} are ignored.}

\item{X}{Optional covariate vector of length \eqn{N = nm}. Used together with \code{loc} to compute \code{xbar} internally.}

\item{loc}{Optional integer vector of length \eqn{N} identifying the spatial unit for each observation in \code{X}.}

\item{gamma}{Tolerance for the maximum acceptable relative difference in posterior variance. Default is \code{0.05}.}
}
\value{
A list with components
\item{\code{m_star}}{Integer replication threshold.}
\item{\code{gamma}}{Tolerance used in the computation.}
\item{\code{numerator_raw}}{Signed value of \eqn{\sum_i d_i^2 (1 - \lambda_i)} before taking the absolute value.}
\item{\code{d_sq}}{Length-\eqn{n} vector of squared projections \eqn{d_i^2}.}
\item{\code{eigenvalues}}{Eigenvalues of the graph Laplacian, sorted in ascending order.}
}
\details{
The function implements
\deqn{m^* = \max\left\{2, \left\lceil \frac{\sigma^2 \rho \left|\sum_{i=1}^{n} d_i^2 (1 - \lambda_i)\right|}{\gamma \tau^2 \left(n - \sum_{i=1}^{n} d_i^2\right)} \right\rceil \right\},}
where \eqn{\lambda_i} are the eigenvalues of the graph Laplacian of \code{W} and \eqn{d_i = \sum_{j} \bar{x}_{j \cdot} u_{ij}} is the projection of the area-level covariate means onto the \eqn{i}th eigenvector.

When \code{X} and \code{loc} are supplied instead of \code{xbar}, the covariate is expected to be pre-standardized to population mean zero and variance one. A warning is issued otherwise.
}
\examples{
set.seed(1)
n <- 25
W <- create_random_W(n)

xbar <- rnorm(n)
xbar <- xbar - mean(xbar)

m_star(W, sigma2 = 0.5, tau2 = 0.5, rho = 0.8, xbar = xbar)
}

\seealso{
\code{\link{SpThreshold}} for fitting the underlying model; \code{\link{create_random_W}} and \code{\link{create_W_from_shapefile}} for adjacency matrix construction.
}
