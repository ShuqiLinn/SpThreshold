# Variance calibration for the proper Leroux model. Public documentation lives
# in man/variance_calibration.Rd.
.threshold_scalar <- function(x, name, lower = 0, upper = Inf,
                              include_lower = FALSE, include_upper = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (include_lower) x < lower else x <= lower) ||
      (if (include_upper) x > upper else x >= upper))
    stop(name, " is outside its permitted range.", call. = FALSE)
  invisible(x)
}

.threshold_graph <- function(W) {
  W <- as.matrix(W)
  if (!is.numeric(W) || nrow(W) < 2L || ncol(W) != nrow(W) ||
      any(!is.finite(W)) || any(W < 0) || any(diag(W) != 0) ||
      !isTRUE(all.equal(W, t(W), tolerance = 1e-12, check.attributes = FALSE)))
    stop("W must be a finite, symmetric, nonnegative adjacency matrix with a zero diagonal.",
         call. = FALSE)
  ee <- eigen(diag(rowSums(W)) - W, symmetric = TRUE)
  ord <- order(ee$values)
  lambda <- pmax(0, ee$values[ord])
  # Proper Leroux covariance is also defined for disconnected graphs. All zero
  # modes are included in the full-field trace calibration.
  list(n = nrow(W), W = W, eigenvalues = lambda,
       eigenvectors = ee$vectors[, ord, drop = FALSE])
}

variance_calibration <- function(W, rho, tau0_2 = NULL, tau1_2 = NULL,
                                tau2 = NULL, calibration = c("full", "common")) {
  calibration <- match.arg(calibration)
  .threshold_scalar(rho, "rho", 0, 1, include_lower = TRUE)
  gg <- .threshold_graph(W)
  c_i <- 1 - rho + rho * gg$eigenvalues
  g <- mean(1 / c_i)
  if (!is.null(tau2)) {
    if (!is.null(tau1_2)) stop("Supply tau1_2 or its legacy alias tau2, not both.", call. = FALSE)
    warning("tau2 is a deprecated spatial-scale alias; prefer tau0_2 (iid scale) or tau1_2.",
            call. = FALSE)
    tau1_2 <- tau2
  }
  if (is.null(tau0_2) && is.null(tau1_2))
    stop("Supply tau0_2 (iid variance) or tau1_2 (spatial scale).", call. = FALSE)
  if (!is.null(tau0_2)) .threshold_scalar(tau0_2, "tau0_2")
  if (!is.null(tau1_2)) .threshold_scalar(tau1_2, "tau1_2")
  supplied_both <- !is.null(tau0_2) && !is.null(tau1_2)
  if (calibration == "common") {
    if (supplied_both && !isTRUE(all.equal(tau0_2, tau1_2)))
      stop("calibration='common' requires equal variances.", call. = FALSE)
    if (is.null(tau0_2)) tau0_2 <- tau1_2
    if (is.null(tau1_2)) tau1_2 <- tau0_2
    source <- "legacy_common_variance"
  } else {
    if (is.null(tau0_2)) tau0_2 <- tau1_2 * g
    if (is.null(tau1_2)) tau1_2 <- tau0_2 / g
    source <- if (supplied_both) "specified_separate_variances" else "full_field_average"
  }
  c(gg, list(rho = rho, c_i = c_i, g = g, tau0_2 = tau0_2,
             tau1_2 = tau1_2, calibration = source))
}
