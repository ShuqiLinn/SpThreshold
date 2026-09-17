# See man/m_star.Rd for the exact and leading-order definitions.
m_star <- function(W, sigma2, tau2 = NULL, rho, xbar = NULL, X = NULL,
                   loc = NULL, gamma = 0.05, tau0_2 = NULL, tau1_2 = NULL,
                   method = c("approx", "exact"), within_ss = NULL,
                   calibration = c("full", "common"), max_search = 1e6) {
  method <- match.arg(method)
  .threshold_scalar(sigma2, "sigma2")
  .threshold_scalar(gamma, "gamma")
  .threshold_scalar(max_search, "max_search", 1)
  if (max_search != floor(max_search)) stop("max_search must be an integer.", call. = FALSE)
  cal <- variance_calibration(W, rho, tau0_2, tau1_2, tau2, calibration)
  dd <- .threshold_balanced(cal, xbar, X, loc, within_ss)
  cc <- .threshold_components(cal, dd, sigma2)
  p <- length(dd$within_ss)
  approx <- pmax(2, ceiling(abs(cc$leading_signed) / gamma))
  approx[dd$within_ss == 0] <- Inf
  approx[cc$identical] <- 2
  sufficient <- pmax(2, ceiling(cc$tail_constant / gamma))
  sufficient[cc$identical] <- 2
  exact <- rep(NA_real_, p)
  searched <- rep(NA_real_, p)
  last_bad <- rep(NA_real_, p)
  status <- ifelse(dd$within_ss > 0, "approximation_only", "no_within_information")
  exact[cc$identical] <- 2
  status[cc$identical] <- "identical_for_this_covariate"
  if (method == "exact") {
    for (j in seq_len(p)) {
      if (cc$identical[j] || dd$within_ss[j] == 0) next
      upper <- sufficient[j]
      stop_at <- min(upper, max_search)
      dj <- list(d_sq = dd$d_sq[, j, drop = FALSE],
                 within_ss = dd$within_ss[j], coefficient = dd$coefficient[j])
      bad <- 1
      # Chunked evaluation caps memory while checking every integer. No
      # monotonicity or first-crossing assumption is made.
      for (start in seq(2, stop_at, by = 10000)) {
        mm <- seq(start, min(stop_at, start + 9999), by = 1)
        rel <- .threshold_curve(mm, cal, dj, sigma2)$relative_difference
        hit <- mm[rel > gamma]
        if (length(hit)) bad <- max(bad, hit)
      }
      searched[j] <- stop_at
      last_bad[j] <- if (bad == 1) NA_real_ else bad
      if (stop_at >= upper) {
        exact[j] <- max(2, bad + 1)
        status[j] <- "certified_integer_threshold"
      } else status[j] <- "search_limit"
    }
  }
  selected <- if (method == "exact") exact else approx
  names(selected) <- names(approx) <- names(exact) <- dd$coefficient
  tab <- data.frame(coefficient = dd$coefficient, m_approx = unname(approx),
    m_exact = unname(exact), sufficient_m = sufficient,
    status = status, searched_through = searched, last_violation = last_bad,
    leading_constant = abs(cc$leading_signed), tail_constant = cc$tail_constant,
    within_ss_per_m = dd$within_ss, stringsAsFactors = FALSE)
  list(m_star = selected, m_approx = approx, m_exact = exact, gamma = gamma,
       method = method, thresholds = tab, g = cal$g, tau0_2 = cal$tau0_2,
       tau1_2 = cal$tau1_2, calibration = cal$calibration,
       numerator_raw = colSums(dd$d_sq * (1 - cal$eigenvalues)),
       d_sq = if (p == 1L) as.vector(dd$d_sq) else dd$d_sq,
       eigenvalues = cal$eigenvalues, within_ss = dd$within_ss,
       target = "Conditional Gaussian coefficient variance; not marginal posterior variance")
}
