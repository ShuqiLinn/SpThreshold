# End-to-end check of the Gaussian sampler against deterministic integration
# of the ORIGINAL, unrestricted proper-Leroux observation covariance.
# Run after installing the updated package:
#   Rscript paper/validation/validate_gaussian_posterior.R
# Optional: SPTHRESHOLD_R_LIB=/path/to/library Rscript ...
# The stored intercept is beta0 + mean(theta); only slopes and covariance
# parameters are compared with the unrestricted model's posterior.

lib <- Sys.getenv("SPTHRESHOLD_R_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
stopifnot(requireNamespace("SpThreshold", quietly = TRUE))

gauss_legendre <- function(k, lower, upper) {
  off <- seq_len(k - 1L) / sqrt(4 * seq_len(k - 1L)^2 - 1)
  J <- matrix(0, k, k)
  J[cbind(seq_len(k - 1L), 2:k)] <- off
  J <- J + t(J)
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = (lower + upper) / 2 + (upper - lower) * e$values[ord] / 2,
       w = (upper - lower) * e$vectors[1, ord]^2)
}

make_data <- function(counts, seed) {
  set.seed(seed)
  n <- length(counts)
  W <- matrix(0, n, n)
  W[cbind(seq_len(n - 1L), 2:n)] <- 1
  W <- W + t(W)
  W[1, n] <- W[n, 1] <- 1
  W[2, 5] <- W[5, 2] <- 1
  loc <- rep(seq_len(n), counts)
  Z <- diag(n)[loc, , drop = FALSE]
  x_area <- c(-1, 0.4, 1.1, -0.6, 0.5, -0.4)
  x_within <- rnorm(length(loc))
  x_within <- x_within - ave(x_within, loc)
  x_within <- x_within + rep(c(0.3, -0.2, 0.4, -0.1, -0.4, 0), counts)
  X <- cbind(`(Intercept)` = 1, within = x_within, area = x_area[loc])
  for (j in 2:3) {
    X[, j] <- X[, j] - mean(X[, j])
    X[, j] <- X[, j] / sqrt(mean(X[, j]^2))
  }
  L <- diag(rowSums(W)) - W
  Q <- 0.65 * L + 0.35 * diag(n)
  theta <- drop(solve(chol(Q), rnorm(n))) * sqrt(0.5)
  y <- drop(X %*% c(0.4, -0.7, 0.9) + Z %*% theta + rnorm(length(loc), sd = sqrt(0.8)))
  list(y = y, X = X, loc = loc, W = W, Z = Z, L = L)
}

# Integrate beta analytically in N(X beta, Omega), retaining the FULL prior
# covariance of theta. No constrained precision, projection, or n-1 prior
# factor from the production sampler is used in this calculation.
reference_moments <- function(dat, spatial, nodes, rho_nodes,
                              log_bounds = c(-5, 3)) {
  sv <- gauss_legendre(nodes, log_bounds[1], log_bounds[2])
  rr <- if (spatial) gauss_legendre(rho_nodes, 0, 1) else list(x = 0, w = 1)
  variances <- exp(sv$x)
  a_sigma <- 6; b_sigma <- 4
  a_tau <- 6; b_tau <- 2.5
  # IG prior times Jacobian for integration with respect to log variance.
  log_prior_s <- -a_sigma * sv$x - b_sigma / variances
  log_prior_v <- -a_tau * sv$x - b_tau / variances
  K <- nodes^2 * length(rr$x)
  values <- matrix(NA_real_, K, 8L)
  colnames(values) <- c("log_weight", "within", "within_second", "area",
                       "area_second", "sigma2", "tau2", "rho")
  N <- length(dat$y)
  p <- ncol(dat$X)
  n <- nrow(dat$W)
  q <- 0L
  for (ir in seq_along(rr$x)) {
    rho <- rr$x[ir]
    Q <- rho * dat$L + (1 - rho) * diag(n)
    B <- dat$Z %*% chol2inv(chol(Q)) %*% t(dat$Z)
    ee <- eigen(B, symmetric = TRUE)
    # Zero eigenvalues arise from the within-area contrast space.
    ev <- pmax(ee$values, 0)
    XR <- crossprod(ee$vectors, dat$X)
    yR <- drop(crossprod(ee$vectors, dat$y))
    rho_log_weight <- log(rr$w[ir]) + if (spatial) dbeta(rho, 2, 2, log = TRUE) else 0
    for (is in seq_len(nodes)) {
      s <- variances[is]
      for (iv in seq_len(nodes)) {
        v <- variances[iv]
        omega_eigen <- s + v * ev
        precision <- 1 / omega_eigen
        info <- crossprod(XR, XR * precision)
        ch <- chol(info)
        beta_var <- chol2inv(ch)
        score <- drop(crossprod(XR, yR * precision))
        beta_mean <- drop(beta_var %*% score)
        residual <- yR - drop(XR %*% beta_mean)
        log_likelihood <- -0.5 * (sum(log(omega_eigen)) +
          2 * sum(log(diag(ch))) + sum(precision * residual^2))
        q <- q + 1L
        values[q, ] <- c(log_likelihood + log_prior_s[is] + log_prior_v[iv] +
          log(sv$w[is]) + log(sv$w[iv]) + rho_log_weight,
          beta_mean[2], beta_var[2, 2] + beta_mean[2]^2,
          beta_mean[3], beta_var[3, 3] + beta_mean[3]^2, s, v, rho)
      }
    }
  }
  weight <- exp(values[, 1] - max(values[, 1]))
  weight <- weight / sum(weight)
  moment <- drop(crossprod(weight, values[, -1, drop = FALSE]))
  names(moment) <- colnames(values)[-1]
  c(within_mean = moment[["within"]],
    within_variance = moment[["within_second"]] - moment[["within"]]^2,
    area_mean = moment[["area"]],
    area_variance = moment[["area_second"]] - moment[["area"]]^2,
    sigma2_mean = moment[["sigma2"]], tau2_mean = moment[["tau2"]],
    if (spatial) c(rho_mean = moment[["rho"]]))
}

batch_mcse <- function(x, batches = 100L) {
  size <- floor(length(x) / batches)
  stopifnot(size > 1)
  batch_means <- colMeans(matrix(x[seq_len(size * batches)], size, batches))
  sd(batch_means) / sqrt(batches)
}

sample_moments <- function(dat, spatial, seed, chains = 2L,
                           burnin = 10000L, retained = 60000L) {
  samples <- vector("list", chains)
  for (chain in seq_len(chains)) {
    set.seed(seed + chain)
    fit <- SpThreshold::spfit(y = dat$y, X = dat$X, loc = dat$loc,
      W = dat$W, model_indicator = as.integer(spatial), family = "gaussian",
      mcmc_samples = burnin + retained + 1L, burnin = burnin,
      a_sigma2_prior = 6, b_sigma2_prior = 4,
      a_tau2_prior = 6, b_tau2_prior = 2.5,
      a_rho_prior = 2, b_rho_prior = 2,
      sigma2_init = c(0.3, 2)[chain], tau2_init = c(0.2, 1.5)[chain],
      rho_init = c(0.15, 0.85)[chain], verbose = FALSE)
    stopifnot(identical(fit$sampler_version, "gaussian_sum_zero_v2"),
              fit$random_effect_rank == nrow(dat$W) - 1L,
              max(abs(colSums(fit$theta))) < 1e-9)
    keep <- seq.int(burnin + 2L, burnin + retained + 1L)
    samples[[chain]] <- cbind(within = fit$beta[2, keep], area = fit$beta[3, keep],
      sigma2 = fit$sigma2[keep], tau2 = fit$tau2[keep],
      if (spatial) cbind(rho = fit$rho[keep]))
  }
  all <- do.call(rbind, samples)
  # Influence functions for posterior variances include uncertainty in the
  # posterior center, and batch means retain serial dependence.
  transforms <- function(z) cbind(
    within_mean = z[, "within"],
    within_variance = (z[, "within"] - mean(all[, "within"]))^2,
    area_mean = z[, "area"],
    area_variance = (z[, "area"] - mean(all[, "area"]))^2,
    sigma2_mean = z[, "sigma2"], tau2_mean = z[, "tau2"],
    if (spatial) cbind(rho_mean = z[, "rho"]))
  moments <- colMeans(transforms(all))
  chain_mcse <- vapply(samples, function(z) apply(transforms(z), 2, batch_mcse),
                       numeric(length(moments)))
  mcse <- sqrt(rowSums(chain_mcse^2)) / chains
  list(estimate = moments, mcse = mcse)
}

validate_gaussian_posterior <- function() {
  designs <- list(balanced = make_data(rep(4L, 6), 8401L),
                  unbalanced = make_data(c(2L, 3L, 5L, 2L, 4L, 6L), 8402L))
  reports <- list()
  for (design in names(designs)) for (spatial in c(FALSE, TRUE)) {
    tag <- paste(design, if (spatial) "spatial" else "iid", sep = "/")
    cat("Checking", tag, "against the unrestricted-covariance posterior ...\n")
    dat <- designs[[design]]
    coarse <- reference_moments(dat, spatial, nodes = 48L, rho_nodes = 32L)
    fine <- reference_moments(dat, spatial, nodes = 64L, rho_nodes = 40L)
    expanded <- reference_moments(dat, spatial, nodes = 80L, rho_nodes = 48L,
                                  log_bounds = c(-6, 4))
    quadrature_error <- pmax(abs(fine - coarse), abs(expanded - fine))
    if (any(quadrature_error > 2e-5 * pmax(1, abs(expanded)))) {
      print(cbind(coarse, fine, expanded, quadrature_error))
      stop("Quadrature has not converged for ", tag, "; refine it before judging the sampler.")
    }
    sampled <- sample_moments(dat, spatial,
      seed = 8410L + 10L * match(design, names(designs)) + as.integer(spatial))
    error <- sampled$estimate - expanded
    allowance <- 5 * sampled$mcse + quadrature_error
    report <- data.frame(design = design, model = if (spatial) "spatial" else "iid",
      quantity = names(expanded), reference = unname(expanded),
      estimate = unname(sampled$estimate), MCSE = unname(sampled$mcse),
      quadrature_error = unname(quadrature_error),
      standardized_error = unname(error / sampled$mcse),
      pass = unname(abs(error) <= allowance), row.names = NULL)
    print(report, digits = 5, row.names = FALSE)
    reports[[tag]] <- report
  }
  report <- do.call(rbind, reports)
  rownames(report) <- NULL
  stopifnot(all(report$pass))
  cat("PASS: all slope and covariance moments agree within 5 batch-MCSEs",
      "plus the checked quadrature error.\n")
  invisible(report)
}

if (sys.nframe() == 0L) validate_gaussian_posterior()
