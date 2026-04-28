#' Fit a Bayesian hierarchical multilevel model for areal spatial data
#'
#' Fits the model
#' \deqn{Y_{ij} = \beta_0 + \beta_1 x_{ij} + \theta_i + \varepsilon_{ij}}
#' under either a Leroux conditional autoregressive prior on \eqn{\theta}
#' (\code{model_indicator = 1}, the spatial model) or an iid Gaussian prior on
#' \eqn{\theta} (\code{model_indicator = 0}, the nonspatial model). Posterior
#' sampling is performed by a Markov chain Monte Carlo algorithm implemented
#' in C++, with adaptive Metropolis tuning for \eqn{\rho} during burn-in.
#'
#' Two input interfaces are supported and are mutually exclusive:
#' \itemize{
#'    \item \strong{Data frame interface:} pass \code{data}, \code{formula}, and
#'      \code{area}.
#'    \item \strong{Vector interface:} pass \code{y}, \code{X}, and \code{loc}.
#' }
#' Mixing the two interfaces is an error.
#'
#' @param data A data frame containing the response, covariate(s), and the
#'   area-identifier column. Required if using the data frame interface.
#' @param formula A two-sided formula of the form \code{y ~ x}, where \code{y}
#'   is the response column in \code{data} and \code{x} is the covariate
#'   column. Required if using the data frame interface.
#' @param area A length-1 character string giving the column name in
#'   \code{data} that identifies each observation's spatial unit. Values must
#'   be integers in \code{1:n}, where \code{n = nrow(W)}. Required if using
#'   the data frame interface.
#'
#' @param y Numeric response vector of length \eqn{N = nm}. Required if using
#'   the vector interface.
#' @param X Numeric design matrix of dimension \eqn{N \times p}, including the
#'   intercept column if needed. Required if using the vector interface.
#' @param loc Integer vector of length \eqn{N} identifying the spatial unit of
#'   each observation. Values must be integers in \code{1:n}. Required if
#'   using the vector interface.
#'
#' @param W An \eqn{n \times n} symmetric adjacency matrix of 0/1 entries.
#'   Required for both spatial and nonspatial models (used to determine
#'   \code{n}). For the spatial model, \code{W} must correspond to a
#'   connected graph.
#' @param model_indicator \code{1} for the spatial (Leroux CAR) model;
#'   \code{0} for the nonspatial (iid random effect) model.
#' @param mcmc_samples Total number of MCMC iterations to draw, including the
#'   burn-in period. The user is responsible for discarding burn-in samples
#'   and applying any thinning post-hoc.
#' @param burnin Number of iterations during which the Metropolis proposal
#'   standard deviation for \eqn{\rho} is adapted (spatial model only). Has
#'   no effect on the returned samples.
#' @param adapt_rho Logical; whether to apply Robbins-Monro adaptation of the
#'   Metropolis proposal SD during \code{burnin}. Spatial model only.
#'
#' @param a_sigma2_prior,b_sigma2_prior Shape and rate parameters of the
#'   inverse-gamma prior on \eqn{\sigma^2}. Defaults to \code{0.01}.
#' @param a_tau2_prior,b_tau2_prior Shape and rate parameters of the
#'   inverse-gamma prior on \eqn{\tau^2}. Defaults to \code{0.01}.
#' @param a_rho_prior,b_rho_prior Shape parameters of the Beta prior on
#'   \eqn{\rho}. Defaults to \code{1.0} (\eqn{\rho \sim \text{Uniform}(0, 1)}).
#'
#' @param beta_init,theta_init,sigma2_init,tau2_init,rho_init Optional
#'   starting values. If \code{NULL}, defaults are used: \code{beta = 0},
#'   \code{theta = 0}, \code{sigma2 = 1}, \code{tau2 = 1}, \code{rho = 0.5}.
#' @param proposal_sd_init Initial Metropolis proposal SD on the logit-\eqn{\rho}
#'   scale. Spatial model only. Default is \code{0.30}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{beta}}{A \eqn{p \times} \code{mcmc_samples} matrix; each
#'       column is one posterior draw of \eqn{\boldsymbol{\beta}}.}
#'     \item{\code{theta}}{An \eqn{n \times} \code{mcmc_samples} matrix.}
#'     \item{\code{sigma2}, \code{tau2}}{Vectors of length \code{mcmc_samples}.}
#'     \item{\code{rho}}{Vector of length \code{mcmc_samples} (spatial model
#'       only).}
#'     \item{\code{accept_rho}}{Total acceptance count of the Metropolis
#'       step (spatial model only).}
#'     \item{\code{final_proposal_sd}}{Final adapted proposal SD on the
#'       logit-\eqn{\rho} scale (spatial model only).}
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{model_indicator}}{The value of \code{model_indicator} used.}
#'   }
#'
#' @details Burn-in and thinning are not applied internally. The returned
#'   sample matrices and vectors all have \code{mcmc_samples} columns/entries,
#'   including the initial state. Apply burn-in by discarding the first
#'   \code{burnin} columns/entries:
#' \preformatted{
#' fit <- spfit(..., mcmc_samples = 11000, burnin = 1000)
#' beta_kept <- fit$beta[, (1001:11000)]
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 25
#' m <- 10
#' W <- create_random_W(n)
#' loc <- rep(1:n, each = m)
#' x <- rnorm(n*m)
#' y <- 1 + 2*x + rnorm(n*m, sd = 0.5)
#'
#' # Vector interface
#' fit <- spfit(y = y, X = cbind(1, x), loc = loc, W = W,
#'              model_indicator = 1, mcmc_samples = 5000, burnin = 1000)
#'
#' # Data frame interface (equivalent)
#' df <- data.frame(y = y, x = x, area = loc)
#' fit2 <- spfit(data = df, formula = y ~ x, area = "area", W = W,
#'               model_indicator = 1, mcmc_samples = 5000, burnin = 1000)
#' }
#'
#' @export
spfit <- function(data = NULL,
                  formula = NULL,
                  area = NULL,
                  y = NULL,
                  X = NULL,
                  loc = NULL,
                  W,
                  model_indicator,
                  mcmc_samples,
                  burnin = 0L,
                  adapt_rho = TRUE,
                  a_sigma2_prior = NULL,
                  b_sigma2_prior = NULL,
                  a_tau2_prior = NULL,
                  b_tau2_prior = NULL,
                  a_rho_prior = NULL,
                  b_rho_prior = NULL,
                  beta_init = NULL,
                  theta_init = NULL,
                  sigma2_init = NULL,
                  tau2_init = NULL,
                  rho_init = NULL,
                  proposal_sd_init = NULL){

   call <- match.call()

   ## -- Required arguments --------------------------------------------------
   if(missing(W)){
      stop("W (adjacency matrix) is required.")
      }
   if(missing(model_indicator)){
      stop("model_indicator is required (0 = nonspatial, 1 = spatial).")
      }
   if(missing(mcmc_samples)){
      stop("mcmc_samples is required.")
      }

   ## -- Determine which interface is being used -----------------------------
   df_args <- list(data = data, formula = formula, area = area)
   vec_args <- list(y = y, X = X, loc = loc)

   df_supplied <- vapply(df_args, function(a) !is.null(a), logical(1))
   vec_supplied <- vapply(vec_args, function(a) !is.null(a), logical(1))

   any_df <- any(df_supplied)
   any_vec <- any(vec_supplied)
   all_df <- all(df_supplied)
   all_vec <- all(vec_supplied)

   if(any_df && any_vec){
      stop("Mix of data frame and vector inputs detected. Use ONE interface ",
           "exclusively: either (data, formula, area) OR (y, X, loc).")
      }

   if(any_df && !all_df){
      missing_df <- names(df_args)[!df_supplied]
      stop(sprintf("Data frame interface incomplete. Missing: %s. ",
                   paste(missing_df, collapse = ", ")),
           "All of (data, formula, area) must be supplied.")
      }

   if(any_vec && !all_vec){
      missing_vec <- names(vec_args)[!vec_supplied]
      stop(sprintf("Vector interface incomplete. Missing: %s. ",
                   paste(missing_vec, collapse = ", ")),
           "All of (y, X, loc) must be supplied.")
      }

   if(!any_df && !any_vec){
      stop("No data supplied. Use either the data frame interface ",
           "(data, formula, area) or the vector interface (y, X, loc).")
      }

   ## -- Translate data frame interface into vectors -------------------------
   if(all_df){

      if(!is.data.frame(data)){
         stop("data must be a data frame.")
         }
      if(!inherits(formula, "formula")){
         stop("formula must be a formula object, e.g. y ~ x.")
         }
      if(!is.character(area) || length(area) != 1L){
         stop("area must be a single character string naming a column in data.")
         }
      if(!(area %in% names(data))){
         stop(sprintf("Column '%s' not found in data.", area))
         }

      ## Extract response and design matrix from formula
      mf <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
      y <- stats::model.response(mf)
      X <- stats::model.matrix(formula, data = mf)
      loc <- data[[area]]

      ## Check for missing values in area column (model.frame already handled
      ## the formula columns)
      if(any(is.na(loc))){
         stop(sprintf("Column '%s' contains NA values.", area))
         }
      }

   ## -- Validate vector inputs ----------------------------------------------
   if(!is.numeric(y)){
      stop("y must be numeric.")
      }
   if(any(is.na(y))){
      stop("y contains NA values.")
      }

   if(!is.numeric(X) || !is.matrix(X)){
      stop("X must be a numeric matrix. ",
           "If you have a single covariate, use cbind(1, x) for an intercept ",
           "and slope, or matrix(x, ncol = 1) for slope only.")
      }
   if(any(is.na(X))){
      stop("X contains NA values.")
      }

   N <- length(y)
   if(nrow(X) != N){
      stop(sprintf("nrow(X) = %d does not match length(y) = %d.",
                   nrow(X), N))
      }

   if(!is.numeric(loc) || any(loc != floor(loc))){
      stop("loc must contain integer values identifying spatial units.")
      }
   loc <- as.integer(loc)
   if(length(loc) != N){
      stop(sprintf("length(loc) = %d does not match length(y) = %d.",
                   length(loc), N))
      }

   ## -- Validate W ----------------------------------------------------------
   if(!is.numeric(W)){
      stop("W must be a numeric matrix.")
      }
   ## Coerce sparse Matrix to dense for the C++ engine, but warn loudly only
   ## if the matrix isn't square or symmetric
   if(inherits(W, "Matrix")){
      W <- as.matrix(W)
      }
   if(!is.matrix(W)){
      stop("W must be a matrix.")
      }
   if(nrow(W) != ncol(W)){
      stop("W must be square.")
      }
   if(!isSymmetric(W, tol = .Machine$double.eps^0.5)){
      stop("W must be symmetric.")
      }
   if(any(W < 0)){
      stop("W must contain non-negative entries (typically 0/1).")
      }
   if(any(diag(W) != 0)){
      stop("Diagonal of W must be zero (no self-loops).")
      }
   if(any(is.na(W))){
      stop("W contains NA values.")
      }

   n <- nrow(W)

   ## Check loc values match W dimensions
   loc_unique <- sort(unique(loc))
   if(min(loc_unique) < 1L || max(loc_unique) > n){
      stop(sprintf("loc values must be integers in 1:%d (= nrow(W)). ",
                   n),
           sprintf("Got range [%d, %d].",
                   min(loc), max(loc)))
      }

   ## Check connectivity if spatial
   if(model_indicator == 1L){
      if(!.is_connected(W)){
         stop("W does not correspond to a connected graph. The spatial ",
              "model requires a single connected component.")
         }
      }

   ## -- Validate model_indicator and mcmc_samples ---------------------------
   if(!(model_indicator %in% c(0L, 1L))){
      stop("model_indicator must be 0 (nonspatial) or 1 (spatial). ",
           sprintf("Got: %s.", model_indicator))
      }
   model_indicator <- as.integer(model_indicator)

   if(!is.numeric(mcmc_samples) || mcmc_samples != floor(mcmc_samples) ||
      mcmc_samples < 1L){
      stop("mcmc_samples must be a positive integer.")
      }
   mcmc_samples <- as.integer(mcmc_samples)

   if(!is.numeric(burnin) || burnin != floor(burnin) || burnin < 0L){
      stop("burnin must be a non-negative integer.")
      }
   burnin <- as.integer(burnin)
   if(burnin >= mcmc_samples){
      stop(sprintf("burnin (%d) must be less than mcmc_samples (%d).",
                   burnin, mcmc_samples))
      }

   if(!is.logical(adapt_rho) || length(adapt_rho) != 1L){
      stop("adapt_rho must be a single logical value.")
      }

   ## -- Validate optional priors --------------------------------------------
   .check_positive_scalar <- function(x, name){
      if(is.null(x)) return(invisible())
      if(!is.numeric(x) || length(x) != 1L || x <= 0){
         stop(sprintf("%s must be a single positive number.", name))
         }
      }
   .check_positive_scalar(a_sigma2_prior, "a_sigma2_prior")
   .check_positive_scalar(b_sigma2_prior, "b_sigma2_prior")
   .check_positive_scalar(a_tau2_prior,   "a_tau2_prior")
   .check_positive_scalar(b_tau2_prior,   "b_tau2_prior")
   .check_positive_scalar(a_rho_prior,    "a_rho_prior")
   .check_positive_scalar(b_rho_prior,    "b_rho_prior")
   .check_positive_scalar(proposal_sd_init, "proposal_sd_init")

   ## -- Validate optional initial values ------------------------------------
   p <- ncol(X)
   if(!is.null(beta_init)){
      if(!is.numeric(beta_init) || length(beta_init) != p){
         stop(sprintf("beta_init must be a numeric vector of length %d (= ncol(X)).",
                      p))
         }
      }
   if(!is.null(theta_init)){
      if(!is.numeric(theta_init) || length(theta_init) != n){
         stop(sprintf("theta_init must be a numeric vector of length %d (= nrow(W)).",
                      n))
         }
      }
   .check_positive_scalar(sigma2_init, "sigma2_init")
   .check_positive_scalar(tau2_init,   "tau2_init")
   if(!is.null(rho_init)){
      if(!is.numeric(rho_init) || length(rho_init) != 1L ||
         rho_init <= 0 || rho_init >= 1){
         stop("rho_init must be a single number in (0, 1).")
         }
      }

   ## -- Convert loc from 1-indexed (R) to 0-indexed (C++) -------------------
   loc_zero <- as.integer(loc - 1L)

   ## -- Call the C++ engine -------------------------------------------------
   fit <- SpThreshold(mcmc_samples       = mcmc_samples,
                      y                  = as.numeric(y),
                      X                  = X,
                      loc                = loc_zero,
                      model_indicator    = model_indicator,
                      W                  = W,
                      proposal_sd_init   = proposal_sd_init,
                      a_sigma2_prior     = a_sigma2_prior,
                      b_sigma2_prior     = b_sigma2_prior,
                      a_tau2_prior       = a_tau2_prior,
                      b_tau2_prior       = b_tau2_prior,
                      a_rho_prior        = a_rho_prior,
                      b_rho_prior        = b_rho_prior,
                      beta_init          = beta_init,
                      theta_init         = theta_init,
                      sigma2_init        = sigma2_init,
                      tau2_init          = tau2_init,
                      rho_init           = rho_init,
                      burnin             = burnin,
                      adapt_rho          = adapt_rho)

   ## -- Attach metadata for the user ----------------------------------------
   fit$call            <- call
   fit$model_indicator <- model_indicator
   fit$mcmc_samples    <- mcmc_samples
   fit$burnin          <- burnin

   fit
   }


## ----------------------------------------------------------------------------
## Internal helpers
## ----------------------------------------------------------------------------

#' Check whether an adjacency matrix corresponds to a connected graph
#'
#' BFS from node 1; returns TRUE if every node is reachable.
#'
#' @param W Symmetric 0/1 adjacency matrix.
#' @return Logical.
#' @keywords internal
#' @noRd
.is_connected <- function(W){

   n <- nrow(W)
   if(n <= 1L) return(TRUE)

   visited <- rep(FALSE, n)
   visited[1L] <- TRUE
   queue <- 1L

   while(length(queue) > 0L){
      current <- queue[1L]
      queue <- queue[-1L]
      neighbors <- which(W[current, ] != 0)
      new <- neighbors[!visited[neighbors]]
      visited[new] <- TRUE
      queue <- c(queue, new)
      }

   all(visited)
   }
