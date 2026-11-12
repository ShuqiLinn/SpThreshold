#' Extract posterior draws
#'
#' @param fit A fit returned by \code{spfit} or \code{spfit_glmm}.
#' @param burnin Number of sampling updates to discard. Defaults to the
#'   burn-in stored in \code{fit}, or zero if unavailable.
#' @param thin Retain every \code{thin}th update after burn-in.
#' @param include_theta Include random-intercept draws.
#' @return A data frame with one row per retained draw.
#' @export
posterior_draws <- function(fit, burnin = NULL, thin = 1L,
                            include_theta = FALSE) {
   nstate <- ncol(fit$beta)
   
   if (is.null(burnin))
      burnin <- if (is.null(fit$burnin)) 0L else fit$burnin
   
   .sp_integer(burnin, "burnin", 0L)
   .sp_integer(thin, "thin", 1L)
   
   if (burnin + thin > nstate - 1L)
      stop("burnin and thin leave no retained sampling updates.")
   
   updates <- seq.int(burnin + thin, nstate - 1L, by = thin)
   cols <- updates + 1L
   
   bn <- rownames(fit$beta)
   if (is.null(bn))
      bn <- paste0("beta[", seq_len(nrow(fit$beta)), "]")
   
   out <- as.data.frame(t(fit$beta[, cols, drop = FALSE]))
   names(out) <- make.unique(bn)
   
   for (nm in c("sigma2", "tau2", "rho")) {
      if (!is.null(fit[[nm]])) {
         if (nm %in% names(out))
            stop("A coefficient name conflicts with parameter name: ", nm)
         
         out[[nm]] <- as.numeric(fit[[nm]])[cols]
      }
   }
   
   if (include_theta) {
      tn <- paste0("theta[", seq_len(nrow(fit$theta)), "]")
      
      if (any(tn %in% names(out)))
         stop("Coefficient names conflict with theta column names.")
      
      th <- as.data.frame(t(fit$theta[, cols, drop = FALSE]))
      names(th) <- tn
      out <- cbind(out, th)
   }
   
   rownames(out) <- NULL
   attr(out, "updates") <- updates
   attr(out, "burnin") <- burnin
   attr(out, "thin") <- thin
   
   out
}

#' Summarize posterior draws
#'
#' @inheritParams posterior_draws
#' @param probs Probabilities for the lower and upper credible interval limits.
#' @return A data frame with posterior means, medians, variances, standard
#'   deviations, interval limits, draw counts, effective sample sizes,
#'   and Monte Carlo standard errors of the means.
#' @export
posterior_summary <- function(fit, burnin = NULL, thin = 1L,
                              include_theta = FALSE,
                              probs = c(0.025, 0.975)) {
   if (!is.numeric(probs) || length(probs) != 2L ||
       any(!is.finite(probs)) || any(probs < 0 | probs > 1) ||
       probs[1] >= probs[2])
      stop("probs must be two increasing probabilities in [0,1].")
   
   draws <- posterior_draws(fit, burnin, thin, include_theta)
   
   if (nrow(draws) < 2L)
      stop("At least two retained draws are required for a variance summary.")
   
   vals <- lapply(draws, function(x) {
      v <- stats::var(x)
      ess <- .sp_ess(x)
      q <- stats::quantile(x, probs = probs, names = FALSE)
      
      c(
         mean = mean(x),
         median = stats::median(x),
         variance = v,
         sd = sqrt(v),
         lower = q[1],
         upper = q[2],
         n_draws = length(x),
         ess = ess,
         mcse_mean = if (is.finite(ess)) sqrt(v / ess) else NA_real_
      )
   })
   
   out <- data.frame(
      parameter = names(draws),
      do.call(rbind, vals),
      row.names = NULL,
      check.names = FALSE
   )
   
   attr(out, "burnin") <- attr(draws, "burnin")
   attr(out, "thin") <- attr(draws, "thin")
   attr(out, "probs") <- probs
   
   out
}

.sp_ess <- function(x) {
   n <- length(x)
   v <- stats::var(x)
   
   if (n < 4L || !is.finite(v) || v == 0)
      return(NA_real_)
   
   ac <- as.numeric(stats::acf(
      x,
      lag.max = min(n - 1L, max(100L, ceiling(10 * sqrt(n)))),
      plot = FALSE,
      demean = TRUE
   )$acf)
   
   k <- length(ac) %/% 2L
   pairs <- ac[seq.int(1L, 2L * k, by = 2L)] +
      ac[seq.int(2L, 2L * k, by = 2L)]
   
   bad <- which(pairs <= 0)
   if (length(bad))
      pairs <- pairs[seq_len(bad[1] - 1L)]
   
   if (!length(pairs))
      return(as.numeric(n))
   
   tau <- -1 + 2 * sum(cummin(pairs))
   min(n, n / max(tau, 1))
}

#' @rdname posterior_summary
#' @param object A fitted model.
#' @param ... Arguments passed to \code{posterior_summary}.
#' @export
summary.spfit <- function(object, ...) {
   posterior_summary(object, ...)
}