# Two-stage Gaussian covariance estimation; see man/two_stage_estimate.Rd.
#' @export
two_stage_estimate <- function(y, X, loc, W, rho_cap = 0.9999) {
  input <- .sp_inputs(NULL, NULL, NULL, y, X, loc, W)
  y <- input$y; X <- input$X; loc <- input$loc; W <- input$W
  n <- nrow(W); N <- length(y)
  
  .threshold_scalar(rho_cap, "rho_cap", 0, 1)
  if (ncol(X) < 2L || !all(X[, 1L] == 1) || qr(X)$rank < ncol(X))
    stop("X must have a first column of ones and full-rank predictor columns.")
  if (N <= n || N <= ncol(X))
    stop("More observations than areas and fixed-effect coefficients are required.")
  if (!.is_connected(W)) stop("W must describe a connected graph.")
  
  m <- tabulate(loc, nbins = n)
  Xb <- sweep(rowsum(X, loc, reorder = TRUE), 1L, m, "/")
  Sxx <- crossprod(X - Xb[loc, , drop = FALSE])
  
  graph <- .threshold_graph(W)
  U <- graph$eigenvectors
  lambda <- graph$eigenvalues
  lambda[1L] <- 0
  
  solve_chol <- function(R, b) backsolve(R, forwardsolve(t(R), b))
  
  # Stage 1: nonspatial REML using every column of X.
  dat <- data.frame(y = y, area = factor(loc, levels = seq_len(n)))
  dat$X <- I(X)
  
  fit <- lme4::lmer(
    y ~ 0 + X + (1 | area), data = dat, REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 100000)
    )
  )
  
  messages <- fit@optinfo$conv$lme4$messages
  messages <- messages[!grepl("boundary \\(singular\\) fit", messages)]
  if (any(fit@optinfo$conv$opt != 0) || length(messages))
    stop("Nonspatial REML did not converge: ",
         paste(messages, collapse = "; "))
  
  tau0_2 <- as.numeric(lme4::VarCorr(fit)$area[1L, 1L])
  sigma2 <- stats::sigma(fit)^2
  .threshold_scalar(sigma2, "fitted sigma2")
  
  # Spatial covariance eigenvalues, with average variance tau0_2.
  spatial_variances <- function(rho) {
    if (tau0_2 == 0) return(rep(0, n))
    inverse_q <- 1 / (1 - rho + rho * lambda)
    tau0_2 * inverse_q / mean(inverse_q)
  }
  
  rho <- NA_real_
  tau1_2 <- 0
  contrast_rank <- NA_integer_
  upper_boundary <- FALSE
  status <- "zero_variance_rho_unidentified"
  
  if (tau0_2 > 0) {
    # Stage 2: express the estimated random intercepts as B %*% y.
    Z <- matrix(0, N, n)
    Z[cbind(seq_len(N), loc)] <- 1
    
    shrink <- tau0_2 / (tau0_2 + sigma2 / m)
    VinvX <- (X - sweep(
      Xb[loc, , drop = FALSE], 1L, shrink[loc], "*"
    )) / sigma2
    
    H <- solve(crossprod(X, VinvX), t(VinvX))
    B <- sweep(
      sweep(t(Z), 1L, m, "/") - Xb %*% H,
      1L, shrink, "*"
    )
    
    effects <- lme4::ranef(fit)$area
    theta <- effects[
      match(as.character(seq_len(n)), rownames(effects)), 1L
    ]
    if (max(abs(as.numeric(B %*% y) - theta)) >
        1e-6 * max(1, max(abs(theta))))
      stop("Reconstructed random intercepts do not match the REML fit.")
    
    # Use nonsingular contrasts and standardize their observation noise.
    s <- svd(B, nu = n, nv = 0)
    keep <- which(s$d > max(s$d) * sqrt(.Machine$double.eps))
    if (!length(keep)) stop("No random-intercept contrasts remain.")
    contrast_rank <- length(keep)
    
    A <- sweep(
      crossprod(s$u[, keep, drop = FALSE], B),
      1L, s$d[keep], "/"
    )
    response <- as.numeric(A %*% y)
    K <- A %*% Z %*% U
    
    nll <- function(rho) {
      C <- tcrossprod(sweep(
        K, 2L, sqrt(spatial_variances(rho)), "*"
      ))
      diag(C) <- diag(C) + sigma2
      R <- chol(C)
      0.5 * (
        2 * sum(log(diag(R))) +
          sum(response * solve_chol(R, response))
      )
    }
    
    # Search all grid-local minima, including both endpoints.
    grid <- seq(0, rho_cap, length.out = 1001L)
    values <- vapply(grid, nll, numeric(1))
    if (any(!is.finite(values)))
      stop("Nonfinite second-stage likelihood.")
    
    mid <- seq.int(2L, length(grid) - 1L)
    local <- mid[
      values[mid] <= values[mid - 1L] &
        values[mid] <= values[mid + 1L]
    ]
    
    candidates <- c(0, rho_cap, grid[which.min(values)])
    for (i in local) {
      candidates <- c(
        candidates,
        stats::optimize(
          nll, grid[c(i - 1L, i + 1L)], tol = 1e-9
        )$minimum
      )
    }
    candidates <- unique(candidates)
    rho <- candidates[which.min(vapply(candidates, nll, numeric(1)))]
    
    tau1_2 <- tau0_2 / mean(1 / (1 - rho + rho * lambda))
    upper_boundary <- rho >= rho_cap - 1e-7
    status <- if (upper_boundary) "rho_upper_cap" else "ok"
    
    if (upper_boundary)
      warning("rho reached rho_cap; inspect the search boundary.")
  }
  
  # Scalar full-conditional variances: other coefficients are held fixed.
  variances <- function(r) {
    G <- tcrossprod(sweep(U, 2L, sqrt(spatial_variances(r)), "*"))
    diag(G) <- diag(G) + sigma2 / m
    information <- Sxx / sigma2 +
      crossprod(Xb, solve_chol(chol(G), Xb))
    as.numeric(1 / diag(information)[-1L])
  }
  
  vi <- variances(0)
  vs <- variances(rho)
  
  list(
    estimates = data.frame(
      rho = rho, tau0_2 = tau0_2, tau1_2 = tau1_2,
      sigma2 = sigma2, upper_boundary = upper_boundary
    ),
    conditional_variances = data.frame(
      coefficient = colnames(X)[-1L],
      spatial_variance = vs,
      iid_variance = vi,
      conditional_gap = abs(vs / vi - 1)
    ),
    iid_fit = fit,
    stage2 = list(
      rho = rho, upper_boundary = upper_boundary,
      contrast_rank = contrast_rank, status = status
    )
  )
}