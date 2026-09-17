.threshold_design <- function(X, loc, n) {
  if (is.null(X) || is.null(loc)) stop("Supply both X and loc.", call. = FALSE)
  if (is.null(dim(X))) X <- matrix(X, ncol = 1L)
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 1L || any(!is.finite(X)))
    stop("X must be a finite numeric vector or matrix.", call. = FALSE)
  if (!is.numeric(loc) || length(loc) != nrow(X) || any(!is.finite(loc)) ||
      any(loc != as.integer(loc)) || any(!loc %in% seq_len(n)))
    stop("loc must contain integer indices 1,...,nrow(W), one per row of X.", call. = FALSE)
  counts <- tabulate(loc, nbins = n)
  if (any(counts == 0L)) stop("Every location in W must occur in loc.", call. = FALSE)
  xsum <- rowsum(X, group = factor(loc, levels = seq_len(n)), reorder = TRUE)
  xbar <- xsum / counts
  ss <- colSums((X - xbar[loc, , drop = FALSE])^2)
  nm <- colnames(X)
  if (is.null(nm)) nm <- paste0("beta", seq_len(ncol(X)))
  colnames(X) <- colnames(xbar) <- nm
  list(X = X, loc = loc, counts = counts, xbar = xbar,
       within_ss = ss, coefficient = nm)
}

.threshold_balanced <- function(cal, xbar, X, loc, within_ss) {
  n <- cal$n
  if (is.null(xbar)) {
    dd <- .threshold_design(X, loc, n)
    if (length(unique(dd$counts)) != 1L)
      stop("Varying-m curves require a balanced template. Use variance_comparison() for an unequal observed design.",
           call. = FALSE)
    xbar <- dd$xbar
    if (is.null(within_ss)) within_ss <- dd$within_ss / dd$counts[1L]
  } else {
    if (!is.null(X) || !is.null(loc))
      stop("Supply xbar or X with loc, not both.", call. = FALSE)
    if (is.null(dim(xbar))) xbar <- matrix(xbar, ncol = 1L)
    xbar <- as.matrix(xbar)
    if (!is.numeric(xbar) || nrow(xbar) != n || ncol(xbar) < 1L || any(!is.finite(xbar)))
      stop("xbar must have nrow(W) finite numeric entries or rows.", call. = FALSE)
    if (is.null(within_ss)) {
      if (any(abs(colMeans(xbar)) > 1e-7))
        stop("For xbar alone, provide centered means under population standardization, or specify within_ss.", call. = FALSE)
      within_ss <- n - colSums(xbar^2)
    }
  }
  p <- ncol(xbar)
  if (!is.numeric(within_ss) || length(within_ss) != p || any(!is.finite(within_ss)))
    stop("within_ss must contain one finite value per coefficient.", call. = FALSE)
  tol <- 100 * .Machine$double.eps * pmax(1, colSums(xbar^2), abs(within_ss))
  if (any(within_ss < -tol))
    stop("within_ss is negative: xbar is inconsistent with the stated population standardization.", call. = FALSE)
  within_ss[abs(within_ss) <= tol] <- 0
  d_sq <- crossprod(cal$eigenvectors, xbar)^2
  if (any(within_ss + colSums(d_sq) == 0))
    stop("A zero covariate has no information for its coefficient.", call. = FALSE)
  nm <- colnames(xbar)
  if (is.null(nm)) nm <- paste0("beta", seq_len(p))
  colnames(d_sq) <- nm
  list(xbar = xbar, d_sq = d_sq, within_ss = within_ss,
       coefficient = nm)
}

.threshold_components <- function(cal, design, sigma2) {
  delta <- 1 / cal$tau0_2 - cal$c_i / cal$tau1_2
  ds <- design$d_sq
  # C0-C1 is a sum of mode differences, each bounded by |delta_i| d_i^2.
  B <- sigma2 * colSums(ds * delta)
  A <- sigma2 * colSums(ds * abs(delta))
  s <- design$within_ss
  identical <- colSums(ds * abs(delta)) == 0
  leading <- ifelse(s > 0, B / s, NA_real_)
  tail <- ifelse(s > 0, A / s, Inf)
  leading[identical] <- 0
  tail[identical] <- 0
  list(delta = delta, leading_signed = leading,
       tail_constant = tail, identical = identical)
}

.threshold_curve <- function(m, cal, design, sigma2) {
  p <- ncol(design$d_sq)
  cc <- .threshold_components(cal, design, sigma2)
  a0 <- m / (sigma2 + m * cal$tau0_2)
  a1 <- outer(m, cal$c_i, "*") /
    (outer(rep(sigma2, length(m)), cal$c_i, "*") +
       outer(m * cal$tau1_2, rep(1, cal$n), "*"))
  within <- outer(m / sigma2, design$within_ss, "*")
  I0 <- within + outer(a0, colSums(design$d_sq), "*")
  I1 <- within + a1 %*% design$d_sq
  signed <- ((a0 - a1) %*% design$d_sq) / I1
  leading <- outer(1 / m, cc$leading_signed, "*")
  leading[, cc$identical] <- 0
  bound <- outer(1 / m, cc$tail_constant, "*")
  bound[, cc$identical] <- 0
  flatten <- function(z) as.vector(t(z))
  data.frame(m = rep(m, each = p), coefficient = rep(design$coefficient, length(m)),
    variance_spatial = flatten(1 / I1), variance_iid = flatten(1 / I0),
    relative_signed = flatten(signed), relative_difference = flatten(abs(signed)),
    leading_signed = flatten(leading), leading_difference = flatten(abs(leading)),
    tail_bound = flatten(bound),
    within_ss_per_m = rep(design$within_ss, length(m)), stringsAsFactors = FALSE)
}

variance_curve <- function(W, sigma2, rho, m, xbar = NULL, X = NULL,
                           loc = NULL, tau0_2 = NULL, tau1_2 = NULL,
                           tau2 = NULL, within_ss = NULL,
                           calibration = c("full", "common")) {
  .threshold_scalar(sigma2, "sigma2")
  if (!is.numeric(m) || !length(m) || any(!is.finite(m)) || any(m <= 0))
    stop("m must contain positive finite replication values.", call. = FALSE)
  cal <- variance_calibration(W, rho, tau0_2, tau1_2, tau2, calibration)
  dd <- .threshold_balanced(cal, xbar, X, loc, within_ss)
  out <- .threshold_curve(m, cal, dd, sigma2)
  attr(out, "calibration") <- cal[c("g", "tau0_2", "tau1_2", "calibration")]
  attr(out, "target") <- "Gaussian variance conditional on covariance parameters and all other regression coefficients; random effects integrated out"
  out
}

variance_comparison <- function(W, sigma2, rho, X, loc,
                                tau0_2 = NULL, tau1_2 = NULL, tau2 = NULL,
                                calibration = c("full", "common")) {
  .threshold_scalar(sigma2, "sigma2")
  cal <- variance_calibration(W, rho, tau0_2, tau1_2, tau2, calibration)
  dd <- .threshold_design(X, loc, cal$n)
  # Orthogonally split individual data into within-location contrasts and the
  # location means. This is equivalent to diag(X' V^-1 X), but avoids an N x N
  # inverse and the cancellation in a direct Woodbury subtraction.
  C1 <- cal$tau1_2 * tcrossprod(sweep(cal$eigenvectors, 2L,
                                                 sqrt(cal$c_i), "/"))
  B1 <- C1 + diag(sigma2 / dd$counts, cal$n)
  I1 <- dd$within_ss / sigma2 + colSums(dd$xbar * solve(B1, dd$xbar))
  I0 <- dd$within_ss / sigma2 + colSums(dd$xbar^2 /
                                     (cal$tau0_2 + sigma2 / dd$counts))
  if (any(I0 <= 0) || any(I1 <= 0))
    stop("Every covariate must have positive conditional information.", call. = FALSE)
  out <- data.frame(coefficient = dd$coefficient,
    variance_spatial = 1 / I1, variance_iid = 1 / I0,
    relative_signed = I0 / I1 - 1, relative_difference = abs(I0 / I1 - 1),
    within_ss = dd$within_ss, stringsAsFactors = FALSE)
  attr(out, "calibration") <- cal[c("g", "tau0_2", "tau1_2", "calibration")]
  attr(out, "counts") <- dd$counts
  attr(out, "target") <- "Gaussian variance conditional on covariance parameters and all other regression coefficients; random effects integrated out"
  out
}
