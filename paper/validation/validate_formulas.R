create_W <- function(n_locs) {
   W <- matrix(0, nrow = n_locs, ncol = n_locs)
   perm <- sample(n_locs)
   for (i in 2:n_locs) {
      a <- perm[i]
      b <- perm[sample.int(i - 1, 1)]
      W[a, b] <- 1
      W[b, a] <- 1
   }
   W
}

run_one_check <- function(n, m) {
   rho <- runif(1, 0, 0.99)
   sigma2 <- runif(1, 0.05, 0.95)
   tau2 <- 1 - sigma2
   x_loc_var <- runif(1, 0, 1)
   
   W <- create_W(n)
   L <- diag(rowSums(W)) - W
   
   loc_means <- rnorm(n, mean = 0, sd = sqrt(x_loc_var))
   x <- as.vector(sapply(loc_means, function(mu) rnorm(m, mean = mu, sd = 1)))
   x <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
   x_bar <- tapply(x, rep(1:n, each = m), mean)
   
   e <- eigen(L)
   ord <- order(e$values)
   lam <- e$values[ord]
   U <- e$vectors[, ord]
   
   d <- sapply(1:n, function(i) sum(x_bar * U[, i])^2)
   
   prec_sp_closed <- (n * m / sigma2) -
      (m^2 * tau2 / sigma2) * sum(d / (sigma2 * (rho * lam + 1 - rho) + m * tau2))
   prec_ns_closed <- (n * m / sigma2) -
      (m^2 * tau2 / (sigma2 * (sigma2 + m * tau2))) * sum(d)
   
   Z <- diag(n) %x% rep(1, m)
   Q_inv_sp <- chol2inv(chol(rho * L + (1 - rho) * diag(n)))
   Omega_sp <- tau2 * Z %*% Q_inv_sp %*% t(Z) + sigma2 * diag(n * m)
   Omega_ns <- tau2 * tcrossprod(Z) + sigma2 * diag(n * m)
   
   prec_sp_matrix <- as.numeric(crossprod(x, chol2inv(chol(Omega_sp)) %*% x))
   prec_ns_matrix <- as.numeric(crossprod(x, chol2inv(chol(Omega_ns)) %*% x))
   
   out <- c(abs(prec_sp_closed - prec_sp_matrix) / abs(prec_sp_matrix),
            abs(prec_ns_closed - prec_ns_matrix) / abs(prec_ns_matrix))
   print(out)
   out
}

set.seed(2026)
n <- 50
m <- 20
n_sims <- 100
results <- replicate(n_sims, run_one_check(n, m))
max(results[1, ])
max(results[2, ])