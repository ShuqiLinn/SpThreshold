library(lme4)
library(SpThreshold)

homes <- read.csv("srrs2.dat", strip.white = TRUE)
county <- read.csv("cty.dat", strip.white = TRUE)
names(homes) <- tolower(trimws(names(homes)))
names(county) <- tolower(trimws(names(county)))
homes <- homes[trimws(homes$state) == "MN", ]
homes$fips <- 1000L * homes$stfips + homes$cntyfips
excluded <- homes[homes$activity == 0, c("idnum", "fips", "activity", "floor")]
homes <- homes[homes$activity > 0, ]
county$fips <- 1000L * county$stfips + county$ctfips
county <- unique(county[c("fips", "uppm")])
uranium <- county$uppm[match(homes$fips, county$fips)]

y <- log(homes$activity)
X <- cbind("(Intercept)" = 1, floor = homes$floor, log_uranium = log(uranium))
covariate_centers <- colMeans(X[, -1, drop = FALSE])
X[, -1] <- scale(X[, -1, drop = FALSE], center = TRUE, scale = FALSE)
fips <- sort(unique(homes$fips))
loc <- match(homes$fips, fips)
n <- length(fips)
N <- length(y)
m <- tabulate(loc, nbins = n)
coefficients <- colnames(X)[-1]

W <- readRDS("radon_adjacency.rds")
W <- W[as.character(fips), as.character(fips)]

spectral <- eigen(diag(rowSums(W)) - W, symmetric = TRUE)
order_eigen <- order(spectral$values)
lambda <- spectral$values[order_eigen]
lambda[1] <- 0
U <- spectral$vectors[, order_eigen]
Xb <- sweep(rowsum(X, loc, reorder = TRUE), 1, m, "/")
Xw <- X - Xb[loc, ]
Sxx <- crossprod(Xw)
Z <- matrix(0, N, n)
Z[cbind(seq_len(N), loc)] <- 1

solve_chol <- function(R, b) backsolve(R, forwardsolve(t(R), b))

dat <- data.frame(y = y, floor = X[, "floor"],
                  log_uranium = X[, "log_uranium"], county = factor(loc))
iid_fit <- lmer(y ~ floor + log_uranium + (1 | county), data = dat,
                REML = TRUE, control = lmerControl(
                    optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
tau0_2 <- as.numeric(VarCorr(iid_fit)$county[1, 1])
sigma2 <- sigma(iid_fit)^2

shrink <- tau0_2 / (tau0_2 + sigma2 / m)
VinvX <- (X - sweep(Xb[loc, ], 1, shrink[loc], "*")) / sigma2
H <- solve(crossprod(X, VinvX), t(VinvX))
B <- sweep(sweep(t(Z), 1, m, "/") - Xb %*% H, 1, shrink, "*")
s <- svd(B, nu = n, nv = 0)
keep <- which(s$d > max(s$d) * sqrt(.Machine$double.eps))
A <- sweep(crossprod(s$u[, keep, drop = FALSE], B), 1, s$d[keep], "/")
z <- as.numeric(A %*% y)
K <- A %*% Z %*% U

spatial_variances <- function(rho) {
    inverse_q <- 1 / (1 - rho + rho * lambda)
    tau0_2 * inverse_q / mean(inverse_q)
}

nll <- function(rho) {
    C <- tcrossprod(sweep(K, 2, sqrt(spatial_variances(rho)), "*"))
    diag(C) <- diag(C) + sigma2
    R <- chol(C)
    0.5 * (2 * sum(log(diag(R))) + sum(z * solve_chol(R, z)))
}

rho_grid <- seq(0, 0.9999, length.out = 1001)
values <- vapply(rho_grid, nll, numeric(1))
mid <- 2:(length(rho_grid) - 1)
minima <- mid[values[mid] <= values[mid - 1] & values[mid] <= values[mid + 1]]
candidates <- c(0, 0.9999, rho_grid[which.min(values)])
for (i in minima) {
    candidates <- c(candidates,
        optimize(nll, rho_grid[c(i - 1, i + 1)], tol = 1e-9)$minimum)
}
candidates <- unique(candidates)
rho <- candidates[which.min(vapply(candidates, nll, numeric(1)))]
tau1_2 <- tau0_2 / mean(1 / (1 - rho + rho * lambda))
estimates <- data.frame(rho, tau0_2, tau1_2, sigma2,
                        upper_boundary = rho >= 0.9999 - 1e-7)

conditional_variance <- function(rho) {
    G <- tcrossprod(sweep(U, 2, sqrt(spatial_variances(rho)), "*"))
    diag(G) <- diag(G) + sigma2 / m
    information <- Sxx / sigma2 + crossprod(Xb, solve_chol(chol(G), Xb))
    setNames(1 / diag(information)[-1], coefficients)
}
v_spatial <- conditional_variance(rho)
v_iid <- conditional_variance(0)
conditional <- data.frame(coefficient = coefficients,
    spatial_variance = unname(v_spatial), iid_variance = unname(v_iid),
    conditional_gap = unname(abs(v_spatial / v_iid - 1)))

sample_summary <- data.frame(homes_before = N + nrow(excluded),
    zero_readings_excluded = nrow(excluded), homes_retained = N, counties = n,
    minimum_county_size = min(m), median_county_size = median(m),
    maximum_county_size = max(m))
print(sample_summary, row.names = FALSE)
print(estimates, digits = 7, row.names = FALSE)

draws <- 20000L
burnin <- 15000L
thin <- 5L
chains <- 4L
seed <- 20261001L
tolerance <- 0.05
settings <- list(draws_per_chain = draws, burnin_updates = burnin, thin = thin,
    chains = chains, seed = seed, tolerance = tolerance,
    updates_per_chain = burnin + draws * thin,
    mcmc_samples = burnin + draws * thin + 1L,
    a_sigma2 = 0.01, b_sigma2 = 0.01, a_tau2 = 0.01, b_tau2 = 0.01,
    a_rho = 1, b_rho = 1)
description <- packageDescription("SpThreshold")
package <- list(version = description$Version,
                RemoteSha = description$RemoteSha, Built = description$Built)
sampled <- list(spatial = vector("list", chains), iid = vector("list", chains))
metadata <- list()
beta_start <- as.numeric(lm.fit(X, y)$coefficients)
start_scale <- c(0.6, 1, 1.6, 2)

for (model in names(sampled)) {
    for (id in seq_len(chains)) {
        cat("Fitting", model, "chain", id, "of", chains, "\n")
        spatial <- model == "spatial"
        chain_seed <- seed + 101L * id + if (spatial) 0L else 10000L
        tau_start <- max(if (spatial) tau1_2 else tau0_2, 0.01) * start_scale[id]
        set.seed(chain_seed)
        fit <- spfit(
            y = y, X = X, loc = loc, W = W,
            model_indicator = as.integer(spatial), family = "gaussian",
            mcmc_samples = settings$mcmc_samples, burnin = burnin, adapt_rho = TRUE,
            a_sigma2_prior = 0.01, b_sigma2_prior = 0.01,
            a_tau2_prior = 0.01, b_tau2_prior = 0.01,
            a_rho_prior = 1, b_rho_prior = 1,
            beta_init = beta_start + c(-0.3, -0.1, 0.1, 0.3)[id],
            theta_init = rep(0, n), sigma2_init = sigma2 * start_scale[id],
            tau2_init = tau_start, rho_init = c(0.15, 0.4, 0.7, 0.95)[id],
            verbose = FALSE
        )
        sampled[[model]][[id]] <- posterior_draws(fit, burnin = burnin,
                                                 thin = thin, include_theta = FALSE)
        metadata[[paste(model, id, sep = "_")]] <- list(
            model = model, chain = id, seed = chain_seed, kernel = fit$kernel,
            acceptance_count = fit$accept_rho, final_proposal_sd = fit$final_proposal_sd
        )
        rm(fit)
    }
}

rank_rhat <- function(v) {
    half <- nrow(v) %/% 2L
    s <- cbind(v[seq_len(half), ], v[nrow(v) - half + seq_len(half), ])
    basic <- function(a) {
        w <- mean(apply(a, 2, var))
        b <- nrow(a) * var(colMeans(a))
        if (w == 0) return(if (b == 0) 1 else Inf)
        sqrt(((nrow(a) - 1) * w + b) / (nrow(a) * w))
    }
    normalize <- function(a) matrix(
        qnorm((rank(a) - 0.375) / (length(a) + 0.25)), nrow(a))
    max(basic(normalize(s)), basic(normalize(abs(s - median(s)))))
}

fit_summary <- lapply(sampled, function(a) {
    beta <- lapply(a, function(x) as.matrix(x[, coefficients, drop = FALSE]))
    all_beta <- do.call(rbind, beta)
    center <- colMeans(all_beta)
    influence <- lapply(beta, function(x)
        nrow(all_beta) / (nrow(all_beta) - 1) * sweep(x, 2, center)^2)
    diagnostic <- Map(function(x, z)
        cbind(x, setNames(as.data.frame(z), paste0("variance_influence_", coefficients))),
        a, influence)
    rhat <- vapply(names(diagnostic[[1]]), function(nm)
        rank_rhat(do.call(cbind, lapply(diagnostic, function(x) x[[nm]]))), numeric(1))
    list(mean = center, variance = apply(all_beta, 2, var),
         ci = apply(all_beta, 2, quantile, probs = c(0.025, 0.975), names = FALSE),
         beta = beta, influence = influence, rhat = rhat)
})

vs <- fit_summary$spatial$variance
vi <- fit_summary$iid$variance
batch_sizes <- c(100L, 250L, 500L, 1000L)
batch_se <- function(f, size) {
    n_batch <- draws %/% size
    vv <- vapply(f$influence, function(x) {
        bm <- rowsum(x[seq_len(n_batch * size), ],
                     rep(seq_len(n_batch), each = size)) / size
        apply(bm, 2, var)
    }, numeric(length(coefficients)))
    sqrt(rowMeans(vv) / (chains * n_batch))
}
errors <- vapply(batch_sizes, function(size) {
    es <- batch_se(fit_summary$spatial, size)
    ei <- batch_se(fit_summary$iid, size)
    sqrt((es / vi)^2 + (vs * ei / vi^2)^2)
}, numeric(length(coefficients)))
colnames(errors) <- paste0("batch_", batch_sizes)
se <- apply(errors, 1, max)
signed <- vs / vi - 1
lo <- signed - qnorm(0.975) * se
hi <- signed + qnorm(0.975) * se

comparison <- data.frame(coefficient = coefficients,
    conditional_gap = conditional$conditional_gap,
    spatial_variance = unname(vs), iid_variance = unname(vi),
    marginal_signed = unname(signed), marginal_gap = unname(abs(signed)),
    MCSE = unname(se), MC_lower = ifelse(lo * hi > 0, pmin(abs(lo), abs(hi)), 0),
    MC_upper = pmax(abs(lo), abs(hi)), row.names = NULL)
comparison$conditional_assessment <- ifelse(comparison$conditional_gap <= tolerance,
                                           "Within", "Above")
comparison$marginal_assessment <- ifelse(comparison$MC_upper <= tolerance, "Within",
    ifelse(comparison$MC_lower > tolerance, "Above", "Unresolved by Monte Carlo interval"))
comparison_percent <- comparison[c("coefficient", "conditional_gap", "marginal_gap",
    "MCSE", "MC_lower", "MC_upper", "conditional_assessment", "marginal_assessment")]
percent_columns <- c("conditional_gap", "marginal_gap", "MCSE", "MC_lower", "MC_upper")
comparison_percent[percent_columns] <- 100 * comparison_percent[percent_columns]
names(comparison_percent)[match(percent_columns, names(comparison_percent))] <-
    paste0(percent_columns, "_pct")

coefficient_summary <- do.call(rbind, lapply(names(fit_summary), function(model) {
    f <- fit_summary[[model]]
    data.frame(model = model, coefficient = coefficients, mean = unname(f$mean),
        variance = unname(f$variance), sd = unname(sqrt(f$variance)),
        CrI_lower = unname(f$ci[1, ]), CrI_upper = unname(f$ci[2, ]))
}))
diagnostics <- data.frame(model = names(fit_summary), chains = chains,
    draws_per_chain = draws, max_rank_split_Rhat = vapply(fit_summary,
        function(f) max(f$rhat), numeric(1)))
detailed_rhat <- do.call(rbind, lapply(names(fit_summary), function(model)
    data.frame(model = model, parameter = names(fit_summary[[model]]$rhat),
               Rhat = unname(fit_summary[[model]]$rhat))))

design <- list(y = y, X = X, loc = loc, W = W, n = n, N = N, m = m,
    Xb = Xb, Sxx = Sxx, U = U, lambda = lambda, coefficient = coefficients,
    excluded = excluded, home_ids = homes$idnum, county_fips = fips,
    covariate_centers = covariate_centers)
two_stage <- list(design = design, estimates = estimates,
    conditional_variances = conditional, iid_fit = iid_fit,
    stage2 = list(rho = rho, upper_boundary = estimates$upper_boundary,
                  contrast_rank = length(keep)), sample_summary = sample_summary,
    county_sizes = data.frame(fips = fips, n_homes = m))
radon_results <- list(design = design, estimates = estimates,
    conditional_variances = conditional, comparison = comparison,
    comparison_percent = comparison_percent, coefficient_summary = coefficient_summary,
    diagnostics = diagnostics, detailed_rhat = detailed_rhat, batch_MCSE = errors,
    draws = sampled, chain_metadata = metadata, settings = settings, package = package,
    sample_summary = sample_summary, two_stage = two_stage, session_info = sessionInfo())
saveRDS(radon_results, "radon_results.rds")

print(comparison_percent, digits = 5, row.names = FALSE)
print(diagnostics, digits = 5, row.names = FALSE)
