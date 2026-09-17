#' Fit Bernoulli, binomial, or negative-binomial spatial mixed models
#'
#' Uses Polya-Gamma augmentation with compiled Gaussian updates. See
#' \code{\link{spfit_glmm}} for the explicit centering constraint and priors.
#' @export
spfit_glmm <- function(y, X, loc, W, model_indicator, mcmc_samples,
                       family = c("bernoulli", "binomial", "negative_binomial"),
                       trials = 1, nb_size = NULL, offset = 0,
                       burnin = 0L, adapt_rho = TRUE,
                       beta_prior_mean = 0, beta_prior_sd = 10,
                       center = TRUE, pg_method = "devroye",
                       a_tau2_prior = 0.01, b_tau2_prior = 0.01,
                       a_rho_prior = 1, b_rho_prior = 1,
                       beta_init = NULL, theta_init = NULL,
                       tau2_init = 1, rho_init = 0.5, proposal_sd_init = 0.30) {
    cl <- match.call()
    family <- match.arg(family)
    pg_method <- match.arg(pg_method, "devroye")
    if (!requireNamespace("pgdraw", quietly = TRUE))
        stop("Install pgdraw to fit the non-Gaussian families.", call. = FALSE)
    scalar <- function(z, name, positive = FALSE, integer = FALSE) {
        if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
            (positive && z <= 0) || (integer && z != floor(z)))
            stop(name, " must be a finite ", if (positive) "positive ",
                 if (integer) "integer." else "number.", call. = FALSE)
    }
    expand <- function(z, len, name, infinite = FALSE) {
        if (!is.numeric(z) || !(length(z) %in% c(1L, len)) ||
            anyNA(z) || (!infinite && any(!is.finite(z))))
            stop(name, " must be numeric, with length 1 or ", len,
                 if (!infinite) ", and finite." else ".", call. = FALSE)
        rep(z, length.out = len)
    }
    if (!is.numeric(y) || !is.null(dim(y)) || length(y) < 1L ||
        any(!is.finite(y)) || any(y < 0 | y != floor(y)))
        stop("y must be a nonempty vector of finite nonnegative counts.", call. = FALSE)
    N <- length(y)
    if (!is.matrix(X) || !is.numeric(X) || nrow(X) != N || ncol(X) < 1L ||
        any(!is.finite(X)))
        stop("X must be a finite numeric matrix with length(y) rows.", call. = FALSE)
    p <- ncol(X)
    if (inherits(W, "Matrix")) W <- as.matrix(W)
    if (!is.matrix(W) || !is.numeric(W) || nrow(W) < 1L ||
        nrow(W) != ncol(W) || any(!is.finite(W)) || any(W < 0) ||
        any(diag(W) != 0) || !isSymmetric(W, tol = 1e-10))
        stop("W must be a finite, symmetric, nonnegative square matrix with zero diagonal.",
             call. = FALSE)
    n <- nrow(W)
    if (!is.numeric(loc) || length(loc) != N || any(!is.finite(loc)) ||
        any(loc != floor(loc)) || any(loc < 1 | loc > n))
        stop("loc must identify each observation with an integer in 1:nrow(W).", call. = FALSE)
    scalar(model_indicator, "model_indicator", integer = TRUE)
    if (!(model_indicator %in% 0:1)) stop("model_indicator must be 0 or 1.", call. = FALSE)
    if (model_indicator == 1 && !.is_connected(W))
        stop("The spatial model requires a connected graph.", call. = FALSE)
    scalar(mcmc_samples, "mcmc_samples", positive = TRUE, integer = TRUE)
    scalar(burnin, "burnin", integer = TRUE)
    if (mcmc_samples > .Machine$integer.max)
        stop("mcmc_samples exceeds the supported integer range.", call. = FALSE)
    if (burnin < 0 || burnin >= mcmc_samples)
        stop("burnin must be between 0 and mcmc_samples - 1.", call. = FALSE)
    if (!is.logical(adapt_rho) || length(adapt_rho) != 1L || is.na(adapt_rho))
        stop("adapt_rho must be TRUE or FALSE.", call. = FALSE)
    if (!is.logical(center) || length(center) != 1L || is.na(center))
        stop("center must be TRUE or FALSE.", call. = FALSE)
    if (center && !any(vapply(seq_len(p), function(j) all(X[, j] == 1), logical(1))))
        stop("center = TRUE requires an intercept column of 1s in X.", call. = FALSE)
    beta_prior_mean <- expand(beta_prior_mean, p, "beta_prior_mean")
    beta_prior_sd <- expand(beta_prior_sd, p, "beta_prior_sd", infinite = TRUE)
    if (any(beta_prior_sd <= 0)) stop("beta_prior_sd must be positive (or Inf).", call. = FALSE)
    prior_precision <- 1 / beta_prior_sd^2
    if (any(is.infinite(beta_prior_sd)) && qr(X)$rank < p)
        stop("A flat coefficient prior requires full column rank in X.", call. = FALSE)
    offset <- expand(offset, N, "offset")
    trials <- expand(trials, N, "trials")
    if (family %in% c("bernoulli", "binomial")) {
        if (!is.null(nb_size)) stop("nb_size applies only to negative_binomial.", call. = FALSE)
        if (any(trials < 1 | trials != floor(trials)))
            stop("trials must be positive integers.", call. = FALSE)
        if (family == "bernoulli" && any(trials != 1))
            stop("Bernoulli trials must equal 1.", call. = FALSE)
        if (any(y > trials)) stop("Binomial counts cannot exceed trials.", call. = FALSE)
        shape <- trials
        kappa <- y - trials / 2
        offset_pg <- offset
    } else {
        if (is.null(nb_size))
            stop("Supply the fixed negative-binomial dispersion as nb_size.", call. = FALSE)
        scalar(nb_size, "nb_size", positive = TRUE, integer = TRUE)
        if (any(trials != 1)) stop("trials does not apply to negative_binomial.", call. = FALSE)
        shape <- y + nb_size
        kappa <- (y - nb_size) / 2
        offset_pg <- offset - log(nb_size)
    }
    if (any(!is.finite(shape)) || any(shape > .Machine$integer.max))
        stop("Polya-Gamma shapes exceed the supported numeric range.", call. = FALSE)
    priors <- list(a_tau2_prior = a_tau2_prior, b_tau2_prior = b_tau2_prior,
                   a_rho_prior = a_rho_prior, b_rho_prior = b_rho_prior,
                   tau2_init = tau2_init, proposal_sd_init = proposal_sd_init)
    for (nm in names(priors)) scalar(priors[[nm]], nm, positive = TRUE)
    scalar(rho_init, "rho_init")
    if (rho_init <= 0 || rho_init >= 1) stop("rho_init must lie in (0,1).", call. = FALSE)
    if (is.null(beta_init)) beta_init <- rep(0, p)
    if (is.null(theta_init)) theta_init <- rep(0, n)
    if (!is.numeric(beta_init) || length(beta_init) != p || any(!is.finite(beta_init)))
        stop("beta_init must have ncol(X) finite entries.", call. = FALSE)
    if (!is.numeric(theta_init) || length(theta_init) != n || any(!is.finite(theta_init)))
        stop("theta_init must have nrow(W) finite entries.", call. = FALSE)
    if (center && abs(sum(theta_init)) > 1e-10 * max(1, sum(abs(theta_init))))
        stop("theta_init must sum to zero when center = TRUE.", call. = FALSE)
    fit <- sp_glmm_pg_cpp(as.numeric(y), X, as.integer(loc - 1L), W,
                          as.integer(model_indicator), as.integer(mcmc_samples),
                          as.integer(burnin), adapt_rho, shape, kappa, offset_pg,
                          beta_prior_mean, prior_precision, as.numeric(beta_init),
                          as.numeric(theta_init), tau2_init, rho_init, proposal_sd_init,
                          a_tau2_prior, b_tau2_prior, a_rho_prior, b_rho_prior,
                          center)
    rownames(fit$beta) <- colnames(X)
    fit$call <- cl
    fit$model_indicator <- as.integer(model_indicator)
    fit$mcmc_samples <- as.integer(mcmc_samples)
    fit$burnin <- as.integer(burnin)
    fit$family <- family
    fit$link <- if (family == "negative_binomial") "log" else "logit"
    fit$trials <- if (family == "negative_binomial") NULL else trials
    fit$nb_size <- nb_size
    fit$offset <- offset
    fit$center <- center
    fit$pg_method <- pg_method
    fit$pg_package <- "pgdraw"
    fit$beta_prior_mean <- beta_prior_mean
    fit$beta_prior_sd <- beta_prior_sd
    fit$engine <- "polya_gamma"
    fit$initial_state_included <- TRUE
    fit$design <- list(y=as.numeric(y),X=X,loc=as.integer(loc),W=W,offset=offset)
    class(fit) <- c("spfit","list")
    fit
}
