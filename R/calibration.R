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
   
   list(
      n = nrow(W),
      W = W,
      eigenvalues = lambda,
      eigenvectors = ee$vectors[, ord, drop = FALSE]
   )
}

variance_calibration <- function(W, rho, tau_ns_2 = NULL, tau_sp_2 = NULL) {
   .threshold_scalar(rho, "rho", 0, 1, include_lower = TRUE)
   
   if (is.null(tau_ns_2) && is.null(tau_sp_2))
      stop("Supply tau_ns_2 or tau_sp_2.", call. = FALSE)
   
   if (!is.null(tau_ns_2))
      .threshold_scalar(tau_ns_2, "tau_ns_2")
   
   if (!is.null(tau_sp_2))
      .threshold_scalar(tau_sp_2, "tau_sp_2")
   
   gg <- .threshold_graph(W)
   c_i <- 1 - rho + rho * gg$eigenvalues
   g <- mean(1 / c_i)
   
   if (is.null(tau_ns_2))
      tau_ns_2 <- tau_sp_2 * g
   
   if (is.null(tau_sp_2))
      tau_sp_2 <- tau_ns_2 / g
   
   c(gg, list(
      rho = rho,
      c_i = c_i,
      g = g,
      tau_ns_2 = tau_ns_2,
      tau_sp_2 = tau_sp_2
   ))
}