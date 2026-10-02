# Regression tests for the Gaussian sum-to-zero posterior.
# Run after installing this source: Rscript tests/gaussian_centering.R
# Base R only; tests call the installed compiled package, not a copied sampler.
library(SpThreshold)

near <- function(actual, expected, label, tolerance = 1e-10) {
  error <- max(abs(actual - expected))
  scale <- max(1, abs(expected))
  if (!is.finite(error) || error > tolerance * scale)
    stop(label, ": error ", signif(error, 5), call. = FALSE)
}

theta_update <- getFromNamespace("theta_update", "SpThreshold")
tau2_update <- getFromNamespace("tau2_update", "SpThreshold")
rho_update <- getFromNamespace("rho_update", "SpThreshold")

n <- 5L
W <- matrix(0, n, n)
W[cbind(seq_len(n - 1L), 2:n)] <- 1
W <- W + t(W)
L <- diag(rowSums(W)) - W
# An independent parameterization theta = H * eta, with no singular covariance.
H <- qr.Q(qr(cbind(rep(1, n), diag(n)[, seq_len(n - 1L)])))[, -1L]
near(crossprod(H), diag(n - 1L), "Contrast basis orthonormality")
near(colSums(H), rep(0, n - 1L), "Contrast basis removes the intercept")

make_case <- function(balanced, spatial) {
  counts <- if (balanced) rep(5L, n) else c(1L, 2L, 5L, 11L, 23L)
  loc <- rep(seq_len(n), counts)
  x <- sin(seq_along(loc) * 0.7) + c(-1, -0.6, 0, 0.4, 1.2)[loc]
  X <- cbind("(Intercept)" = 1, x = x - mean(x))
  y <- 0.3 + 0.8 * X[, 2L] + c(-0.9, 0.4, 0.8, -0.5, 0.2)[loc] +
    cos(seq_along(loc) * 1.3)
  Z <- diag(n)[loc, , drop = FALSE]
  rho <- if (spatial) 0.73 else 0
  Q <- rho * L + (1 - rho) * diag(n)
  sigma2 <- 0.7
  tau2 <- 1.1
  precision_eta <- crossprod(H, (diag(counts / sigma2) + Q / tau2) %*% H)
  list(counts = counts, loc = loc, X = X, y = y, Z = Z, rho = rho, Q = Q,
       sigma2 = sigma2, tau2 = tau2, precision_eta = precision_eta,
       root_eta = chol(precision_eta),
       label = paste(if (balanced) "balanced" else "unequal",
                     if (spatial) "spatial" else "iid"),
       spatial = spatial)
}

conditional_mean_eta <- function(case, beta) {
  solve(case$precision_eta,
        crossprod(H, crossprod(case$Z, case$y - drop(case$X %*% beta))) /
          case$sigma2)
}

fit_once <- function(case, theta_init = rep(0, n), beta_init = c(0, 0)) {
  spfit(y = case$y, X = case$X, loc = case$loc, W = W,
        model_indicator = as.integer(case$spatial), mcmc_samples = 2L,
        burnin = 0L, adapt_rho = FALSE, verbose = FALSE,
        beta_init = beta_init, theta_init = theta_init,
        sigma2_init = case$sigma2, tau2_init = case$tau2,
        rho_init = 0.73, proposal_sd_init = 1.4,
        a_sigma2_prior = 2.1, b_sigma2_prior = 0.8,
        a_tau2_prior = 1.7, b_tau2_prior = 0.9,
        a_rho_prior = 1.2, b_rho_prior = 1.4)
}

contrast_logdet <- function(Q) {
  2 * sum(log(diag(chol(crossprod(H, Q %*% H)))))
}

# Check every optimized branch against the n-1-dimensional Gaussian posterior.
# Standardizing conditional draws in this independent basis tests both mean and
# covariance. In particular, subtracting the arithmetic mean under unequal
# replication fails these checks even though the output sums to zero.
cases <- list(make_case(TRUE, FALSE), make_case(TRUE, TRUE),
              make_case(FALSE, FALSE), make_case(FALSE, TRUE))
repetitions <- 4000L
for (k in seq_along(cases)) {
  case <- cases[[k]]
  set.seed(18310L + k)
  standardized <- matrix(NA_real_, repetitions, n - 1L)
  for (b in seq_len(repetitions)) {
    fit <- fit_once(case)
    stopifnot(identical(fit$sampler_version, "gaussian_sum_zero_v2"),
              fit$random_effect_rank == n - 1L)
    theta <- fit$theta[, 2L]
    near(sum(theta), 0, paste(case$label, "sum-zero draw"))
    standardized[b, ] <- drop(case$root_eta %*%
      (crossprod(H, theta) - conditional_mean_eta(case, fit$beta[, 2L])))
  }
  mean_error <- max(abs(colMeans(standardized)))
  covariance_error <- abs(cov(standardized) - diag(n - 1L))
  # Eight standard errors, with separate normal-theory variances for diagonal
  # and off-diagonal sample covariances; deterministic seed and ample margin.
  limits <- matrix(8 / sqrt(repetitions - 1L), n - 1L, n - 1L)
  diag(limits) <- 8 * sqrt(2 / (repetitions - 1L))
  if (mean_error > 8 / sqrt(repetitions) || any(covariance_error > limits))
    stop(case$label, ": draws do not match the contrast-space conditional")

  # Reconstruct the RNG stream for one full update. This checks the actual
  # tau and rho updates in spfit, including branches that bypass the helpers.
  different_legacy_rho_decisions <- 0L
  for (seed in 19001:19100) {
    set.seed(seed)
    fit <- fit_once(case)
    beta <- fit$beta[, 2L]
    theta <- fit$theta[, 2L]
    set.seed(seed)
    invisible(rnorm(ncol(case$X) + n))
    tau_rate <- 0.9 + drop(crossprod(theta, case$Q %*% theta)) / 2
    expected_tau <- 1 / rgamma(1, shape = 1.7 + (n - 1L) / 2,
                               rate = tau_rate)
    near(fit$tau2[2L], expected_tau, paste(case$label, "tau conditional rank"))
    residual <- case$y - drop(case$X %*% beta) - theta[case$loc]
    expected_sigma <- 1 / rgamma(1, shape = 2.1 + length(residual) / 2,
                                 rate = 0.8 + sum(residual^2) / 2)
    near(fit$sigma2[2L], expected_sigma, paste(case$label, "sigma conditional"))
    if (case$spatial) {
      rho_old <- 0.73
      rho_proposal <- plogis(rnorm(1, qlogis(rho_old), 1.4))
      Q_proposal <- rho_proposal * L + (1 - rho_proposal) * diag(n)
      log_ratio <- (contrast_logdet(Q_proposal) - contrast_logdet(case$Q)) / 2 -
        drop(crossprod(theta, (Q_proposal - case$Q) %*% theta)) / (2 * expected_tau) +
        1.2 * log(rho_proposal / rho_old) +
        1.4 * (log1p(-rho_proposal) - log1p(-rho_old))
      log_uniform <- log(runif(1))
      accepted <- log_uniform < log_ratio
      expected_rho <- if (accepted) rho_proposal else rho_old
      near(fit$rho[2L], expected_rho, paste(case$label, "rho contrast determinant"))
      legacy_ratio <- log_ratio + (log1p(-rho_proposal) - log1p(-rho_old)) / 2
      different_legacy_rho_decisions <- different_legacy_rho_decisions +
        as.integer(accepted != (log_uniform < legacy_ratio))
    }
  }
  if (case$spatial && different_legacy_rho_decisions == 0L)
    stop("Seeds did not distinguish the full and contrast rho determinants")
}

# The standalone compiled helpers must implement the same target. For theta,
# compare exact seeded draws to conditioning an unrestricted Gaussian, and
# verify that its induced contrast covariance is the independent oracle above.
for (case in cases) {
  beta <- c(0.2, -0.4)
  precision <- diag(case$counts / case$sigma2) + case$Q / case$tau2
  covariance <- solve(precision)
  covariance_one <- drop(covariance %*% rep(1, n))
  projection <- diag(n) - tcrossprod(covariance_one, rep(1, n)) / sum(covariance_one)
  near(crossprod(H, projection %*% covariance %*% t(projection) %*% H),
       solve(case$precision_eta), paste(case$label, "conditional covariance geometry"))
  seed <- 71911L
  set.seed(seed)
  observed <- theta_update(length(case$y), n, case$y, case$X,
    as.integer(case$loc - 1L), as.numeric(case$counts), beta,
    case$sigma2, case$tau2, case$Q)
  set.seed(seed)
  mean_unrestricted <- solve(precision,
    crossprod(case$Z, case$y - drop(case$X %*% beta)) / case$sigma2)
  unrestricted <- mean_unrestricted + backsolve(chol(precision), rnorm(n))
  expected <- drop(projection %*% unrestricted)
  near(observed, expected, paste(case$label, "standalone theta helper"))
  set.seed(seed)
  observed_tau <- tau2_update(n, expected, case$Q, 1.7, 0.9)
  set.seed(seed)
  expected_tau <- 1 / rgamma(1, shape = 1.7 + (n - 1L) / 2,
    rate = 0.9 + drop(crossprod(expected, case$Q %*% expected)) / 2)
  near(observed_tau, expected_tau, paste(case$label, "standalone tau helper"))
}

theta <- c(-1.0, 0.6, 0.1, -0.4, 0.7)
rho_old <- 0.83
Q_old <- rho_old * L + (1 - rho_old) * diag(n)
full_logdet <- function(Q) as.numeric(determinant(Q, logarithm = TRUE)$modulus)
for (seed in 40101:40200) {
  set.seed(seed)
  observed <- rho_update(n, W, theta, 0.8, rho_old, Q_old,
                        full_logdet(Q_old), 1.8, 1.2, 1.4)
  set.seed(seed)
  proposal <- plogis(rnorm(1, qlogis(rho_old), 1.8))
  Q_proposal <- proposal * L + (1 - proposal) * diag(n)
  log_ratio <- (contrast_logdet(Q_proposal) - contrast_logdet(Q_old)) / 2 -
    drop(crossprod(theta, (Q_proposal - Q_old) %*% theta)) / (2 * 0.8) +
    1.2 * log(proposal / rho_old) +
    1.4 * (log1p(-proposal) - log1p(-rho_old))
  accepted <- log(runif(1)) < log_ratio
  expected_rho <- if (accepted) proposal else rho_old
  expected_Q <- if (accepted) Q_proposal else Q_old
  near(observed$rho, expected_rho, "Standalone rho helper target")
  near(observed$Q, expected_Q, "Standalone rho helper returned Q")
  # Compatibility contract: this helper still accepts/returns the full Q logdet.
  near(observed$Q_log_det, full_logdet(expected_Q), "Standalone rho helper full logdet")
}

# A nonzero initial random-effect mean is absorbed into the intercept without
# changing any observation's initial linear predictor, including unequal sizes.
case <- cases[[4L]]
beta_initial <- c(0.8, -0.2)
theta_initial <- c(1, 2, -1, 0, 3)
fit <- fit_once(case, theta_initial, beta_initial)
near(sum(fit$theta[, 1L]), 0, "Centered initial state")
near(drop(case$X %*% fit$beta[, 1L]) + fit$theta[case$loc, 1L],
     drop(case$X %*% beta_initial) + theta_initial[case$loc],
     "Initial centering preserves the predictor")

# A centered theta prior has the same slope posterior as the original proper
# Leroux prior after integrating the unrestricted flat intercept. Check this
# separately for both designs and covariance models by dense Gaussian algebra.
for (case in cases) {
  original <- case$sigma2 * diag(length(case$y)) +
    case$tau2 * case$Z %*% solve(case$Q) %*% t(case$Z)
  centered <- case$sigma2 * diag(length(case$y)) +
    case$tau2 * case$Z %*% H %*% solve(crossprod(H, case$Q %*% H)) %*%
      t(H) %*% t(case$Z)
  beta_moments <- function(Omega) {
    OX <- solve(Omega, case$X)
    variance <- solve(crossprod(case$X, OX))
    mean <- variance %*% crossprod(OX, case$y)
    c(mean = mean[2L], variance = variance[2L, 2L])
  }
  near(beta_moments(centered), beta_moments(original),
       paste(case$label, "proper and centered slope posterior equivalence"))
}

cat("Gaussian centering regression checks passed: all optimized kernels, helpers,\n",
    "contrast Gaussian geometry, hyperparameter updates and slope equivalence.\n", sep = "")
