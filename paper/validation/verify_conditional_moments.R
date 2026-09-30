# Numerical verification of the four expressions in Proposition 1.

create_W <- function(n_locs) {
  W <- matrix(0, n_locs, n_locs)
  perm <- sample.int(n_locs)
  for (i in 2:n_locs) {
    a <- perm[i]
    b <- perm[sample.int(i - 1L, 1L)]
    W[a, b] <- W[b, a] <- 1
  }
  W
}

run_one_check <- function(n, m) {
  stopifnot(n >= 2, m >= 1, n == as.integer(n), m == as.integer(m))
  N <- n * m
  rho <- runif(1, 0, 0.99)
  sigma2 <- runif(1, 0.05, 0.95)
  tau0_2 <- 1 - sigma2
  beta0 <- runif(1, 0.5, 2)  # Nonzero intercept checks its cancellation.
  beta1 <- rnorm(1)

  W <- create_W(n)
  L <- diag(rowSums(W)) - W
  eig <- eigen(L, symmetric = TRUE)
  ord <- order(eig$values)
  lam <- pmax(eig$values[ord], 0)  # Remove roundoff at the zero eigenvalue.
  U <- eig$vectors[, ord, drop = FALSE]
  q <- rho * lam + 1 - rho
  tau1_2 <- n * tau0_2 / sum(1 / q)

  x_loc_var <- runif(1, 0, 1)
  loc_means <- rnorm(n, sd = sqrt(x_loc_var))
  x <- rnorm(N, mean = rep(loc_means, each = m), sd = 1)
  x <- x - mean(x)
  x <- x / sqrt(mean(x^2))
  x_bar <- colMeans(matrix(x, nrow = m, ncol = n))

  Z <- kronecker(diag(n), matrix(1, nrow = m, ncol = 1))
  Q <- rho * L + (1 - rho) * diag(n)
  Q_inv <- chol2inv(chol(Q))
  theta <- drop(t(chol(tau1_2 * Q_inv)) %*% rnorm(n))
  y <- beta0 + beta1 * x + rep(theta, each = m) +
    rnorm(N, sd = sqrt(sigma2))
  y_bar <- colMeans(matrix(y, nrow = m, ncol = n))

  # Closed forms in Proposition 1.
  x_proj <- drop(crossprod(U, x_bar))
  y_proj <- drop(crossprod(U, y_bar))
  d <- x_proj^2
  xy <- sum(x * y)
  denom_sp <- sigma2 * q + m * tau1_2
  denom_ns <- sigma2 + m * tau0_2

  prec_sp_closed <- N / sigma2 -
    m^2 * tau1_2 / sigma2 * sum(d / denom_sp)
  prec_ns_closed <- N / sigma2 -
    m^2 * tau0_2 / (sigma2 * denom_ns) * sum(d)

  mean_sp_closed <- (xy - m^2 * tau1_2 *
    sum(x_proj * y_proj / denom_sp)) / (sigma2 * prec_sp_closed)
  mean_ns_closed <- (xy - m^2 * tau0_2 / denom_ns *
    sum(x_proj * y_proj)) / (sigma2 * prec_ns_closed)

  # Independent matrix calculations, retaining beta0 in the mean formula.
  Omega_sp <- tau1_2 * Z %*% Q_inv %*% t(Z) + sigma2 * diag(N)
  Omega_ns <- tau0_2 * tcrossprod(Z) + sigma2 * diag(N)

  matrix_moments <- function(Omega) {
    R <- chol(Omega)
    solved <- backsolve(R, forwardsolve(t(R), cbind(x, y - beta0)))
    precision <- sum(x * solved[, 1])
    c(precision = precision, mean = sum(x * solved[, 2]) / precision)
  }
  sp <- matrix_moments(Omega_sp)
  ns <- matrix_moments(Omega_ns)

  mean_sp_error <- abs(mean_sp_closed - unname(sp["mean"]))
  mean_ns_error <- abs(mean_ns_closed - unname(ns["mean"]))
  c(
    precision_sp_relative = abs(prec_sp_closed / unname(sp["precision"]) - 1),
    precision_ns_relative = abs(prec_ns_closed / unname(ns["precision"]) - 1),
    mean_sp_scaled = mean_sp_error / max(1, abs(unname(sp["mean"]))),
    mean_ns_scaled = mean_ns_error / max(1, abs(unname(ns["mean"]))),
    mean_sp_absolute = mean_sp_error,
    mean_ns_absolute = mean_ns_error
  )
}

set.seed(2026)
n <- 50
m <- 20
n_sims <- 100
results <- replicate(n_sims, run_one_check(n, m))
max_errors <- apply(results, 1, max)
print(max_errors, digits = 4)

