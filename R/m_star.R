#' Replication threshold for spatial vs non-spatial inference on beta_1
#'
#' Computes the asymptotic replication threshold \eqn{m^*} beyond which a
#' non-spatial model yields effectively equivalent posterior inference for a
#' regression coefficient compared to a Leroux CAR spatial model, based on the
#' closed-form bound derived in the manuscript.
#'
#' @param W Adjacency matrix, an \eqn{n \times n} symmetric matrix of 0/1
#'   entries, where \eqn{w_{ii'} = 1} indicates that areas \eqn{i} and
#'   \eqn{i'} are neighbors. Must correspond to a connected graph.
#' @param sigma2 Within-area residual variance.
#' @param tau2 Spatial random-effect variance.
#' @param rho Spatial correlation parameter in \eqn{[0, 1)}.
#' @param xbar Optional length-\eqn{n} vector of location-level covariate
#'   means. If supplied, \code{X} and \code{loc} are ignored.
#' @param X Optional covariate vector (length \eqn{N}, where \eqn{N = nm}).
#'   Used together with \code{loc} to compute \code{xbar} internally.
#' @param loc Optional integer vector (length \eqn{N}) identifying the spatial
#'   unit for each observation in \code{X}.
#' @param gamma Tolerance for the maximum acceptable relative difference in
#'   posterior variance. Default is \code{0.05}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{m_star}}{Integer replication threshold.}
#'     \item{\code{gamma}}{Tolerance used in the computation.}
#'     \item{\code{numerator_raw}}{Signed value of
#'       \eqn{\sum_i d_i^2 (1 - \lambda_i)} before taking the absolute value;
#'       useful for diagnosing the sign of the untransformed bound.}
#'     \item{\code{d_sq}}{Length-\eqn{n} vector of squared projections
#'       \eqn{d_i^2}.}
#'     \item{\code{eigenvalues}}{Eigenvalues of the graph Laplacian, sorted in
#'       ascending order.}
#'   }
#'
#' @details The function implements
#' \deqn{m^* = \max\left\{2, \left\lceil \frac{\sigma^2 \rho \left|\sum_{i=1}^{n} d_i^2 (1 - \lambda_i)\right|}{\gamma \tau^2 \left(n - \sum_{i=1}^{n} d_i^2\right)} \right\rceil \right\},}
#' where \eqn{\lambda_i} are the eigenvalues of the graph Laplacian of
#' \code{W} and \eqn{d_i = \sum_{j} \bar{x}_{j \cdot} u_{ij}} is the projection
#' of the location-level covariate means onto the \eqn{i}th eigenvector. The
#' derivation is given in Section 3.2 of the paper.
#'
#' When \code{X} and \code{loc} are supplied instead of \code{xbar}, the
#' covariate is expected to be pre-standardized to population mean zero and
#' variance one; a warning is issued otherwise.
#'
#' @examples
#' # A small random connected graph
#' set.seed(1)
#' n <- 25
#' W <- matrix(0, n, n)
#' perm <- sample(n)
#' for (i in 2:n) {
#'    a <- perm[i]; b <- perm[sample.int(i - 1, 1)]
#'    W[a, b] <- 1; W[b, a] <- 1
#' }
#'
#' # Covariate varying across locations
#' xbar <- rnorm(n)
#' xbar <- xbar - mean(xbar)
#'
#' m_star(W, sigma2 = 0.5, tau2 = 0.5, rho = 0.8, xbar = xbar)
#'
#'
#' @export
m_star <- function(W, sigma2, tau2, rho,
                   xbar = NULL,
                   X = NULL,
                   loc = NULL,
                   gamma = 0.05) {

   ## -- Input dispatch ------------------------------------------------------
   n <- nrow(W)

   if (is.null(xbar) && (is.null(X) || is.null(loc))) {
      stop("Supply either xbar (location-level means) or both X and loc.")
   }

   if (is.null(xbar)) {
      ## Option B: compute xbar from raw covariate and location indicator
      if (length(X) != length(loc)) {
         stop("X and loc must have the same length.")
      }
      if (!all(loc %in% seq_len(n))) {
         stop("loc values must be integers in 1, ..., nrow(W).")
      }

      ## Check standardization: population mean ~ 0 and population variance ~ 1
      x_mean <- mean(X)
      x_var  <- mean((X - x_mean)^2)

      if (abs(x_mean) > 1e-6 || abs(x_var - 1) > 1e-3) {
         warning("X does not appear to be standardized to population mean 0 and variance 1. ",
                 "The bound assumes this standardization; results may be misleading otherwise.")
      }

      xbar <- tapply(X, loc, mean)
      xbar <- as.numeric(xbar)
   }

   if (length(xbar) != n) {
      stop(sprintf("length(xbar) = %d does not match nrow(W) = %d.",
                   length(xbar), n))
   }

   ## -- Validate variance and correlation arguments -------------------------
   if (sigma2 <= 0 || tau2 <= 0) {
      stop("sigma2 and tau2 must be positive.")
   }
   if (rho < 0 || rho >= 1) {
      stop("rho must lie in [0, 1).")
   }
   if (gamma <= 0 ) {
      stop("gamma must be positive")
   }

   ## -- Eigendecomposition of the graph Laplacian ---------------------------
   L <- diag(rowSums(W)) - W
   eig <- eigen(L, symmetric = TRUE)

   ## Order ascending so lambda_1 = 0
   ord <- order(eig$values)
   lambda <- eig$values[ord]
   U      <- eig$vectors[, ord]

   ## -- Projections d_i and summations --------------------------------------
   ## d_i = sum_j xbar_j u_{ij} = (U^T xbar)_i
   d <- as.numeric(crossprod(U, xbar))
   d_sq <- d^2
   sum_d_sq <- sum(d_sq)

   num_signed <- sum(d_sq * (1 - lambda))
   den <- n - sum_d_sq

   ## -- Asymptotic bound ----------------------------------------------------
   ## When the covariate has no within-location variation (C3), standardization 
   ## yields sum(d_sq) ~ n up to floating-point noise, so the denominator n - sum(d_sq) collapses to (numerically) zero
   ## and the bound is infinite.  A tolerance of sqrt(.Machine$double.eps) * n
   ## on the denominator is used to catch this case robustly.
   if (abs(den) < sqrt(.Machine$double.eps) * n) {
      return(list(m_star        = Inf,
                  gamma       = gamma,
                  numerator_raw = num_signed,
                  d_sq          = d_sq,
                  eigenvalues   = lambda))
   }

   m_raw <- (sigma2 * rho * abs(num_signed)) / (gamma * tau2 * den)
   m_out <- max(2, ceiling(m_raw))

   list(m_star        = m_out,
        gamma       = gamma,
        numerator_raw = num_signed,
        d_sq          = d_sq,
        eigenvalues   = lambda)

}
