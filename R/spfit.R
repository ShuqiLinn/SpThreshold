#' Fit a Bayesian multilevel model with spatial or iid random effects
#'
#' Gaussian models use a Leroux CAR or iid Gaussian random-effect prior.
#' Bernoulli, binomial and negative-binomial models use Polya--Gamma augmentation.
#' Multiple covariates and unequal observation counts are supported.
#'
#' @param data,formula,area Data-frame interface: a data frame, a two-sided
#'   formula, and the name of its integer location column. Formula offsets are
#'   supported. Use this interface or the vector interface, not both.
#' @param y,X,loc Vector interface: numeric response, numeric design matrix,
#'   and integer location identifiers in \code{1:nrow(W)}. Include an intercept
#'   in \code{X} if required. For binomial responses, use success counts in
#'   \code{y} and specify \code{trials} separately.
#' @param W Symmetric nonnegative adjacency matrix with zero diagonal.
#' @param model_indicator One for the spatial model; zero for iid random effects.
#' @param mcmc_samples Number of saved states, including the initial state.
#'   There are \code{mcmc_samples - 1} sampling updates.
#' @param burnin Number of updates during which rho proposals are adapted.
#'   Draws are not discarded internally. See \code{\link{posterior_draws}}.
#' @param adapt_rho Adapt the logit-scale rho proposal during burn-in.
#' @param a_sigma2_prior,b_sigma2_prior Gaussian residual inverse-gamma shape
#'   and rate parameters. NULL uses 0.01.
#' @param a_tau2_prior,b_tau2_prior Random-effect inverse-gamma shape and rate.
#'   NULL uses 0.01.
#' @param a_rho_prior,b_rho_prior Beta prior parameters for rho. NULL uses one.
#' @param beta_init,theta_init,sigma2_init,tau2_init,rho_init Starting values.
#'   NULL uses zero coefficients/random effects, unit variances, and rho=0.5.
#' @param proposal_sd_init Initial logit-rho proposal SD. NULL uses 0.3.
#' @param family One of \code{"gaussian"}, \code{"bernoulli"},
#'   \code{"binomial"}, or \code{"negative_binomial"}; \code{"negbin"} is an alias.
#' @param trials Binomial trial counts; scalar or one per observation.
#' @param nb_size Fixed positive negative-binomial size. Required for that
#'   family; dispersion is not estimated. See \code{\link{spfit_glmm}}.
#' @param offset Known linear-predictor offset, scalar or one per observation.
#' @param beta_prior_mean,beta_prior_sd Regression prior for count/binary
#'   families: independent normal means and SDs, default zero and ten. Infinite
#'   SD requests a flat prior. Gaussian fits use the existing flat prior.
#' @param center Sum-to-zero constraint for the count/binary families, using
#'   a constrained Gaussian draw. The Gaussian sampler keeps its existing
#'   centering convention; \code{center=FALSE} is not supported for Gaussian.
#' @param pg_method \code{"devroye"} uses exact integer-shape Polya--Gamma
#'   sampling. \code{"hybrid"} permits noninteger NB sizes using BayesLogit's
#'   hybrid algorithm, which can use numerical approximations.
#' @param verbose Print Gaussian sampling progress.
#' @return A \code{spfit} list containing draws of beta (parameters by states),
#'   theta (locations by states), tau2 and, for spatial models, rho. Gaussian
#'   fits also return sigma2. Metadata identify the family, initial-state
#'   convention, burn-in, and fitted design. The Gaussian \code{kernel} records
#'   which numerical implementation was used.
#' @details Gaussian sampling preserves the original update order, priors and
#'   centering. Binary/count families have explicitly specified regression
#'   priors and a proper sum-to-zero implementation when \code{center=TRUE}.
#'   Their thresholds are not given by the Gaussian replication formulas.
#'   Use \code{posterior_draws()} or \code{posterior_summary()} to discard the
#'   initial state and the requested burn-in updates consistently.
#' @examples
#' set.seed(1)
#' W <- create_random_W(5)
#' loc <- rep(1:5, each = 4)
#' x <- rnorm(20)
#' y <- 1 + x + rnorm(20)
#' fit <- spfit(y = y, X = cbind(1, x), loc = loc, W = W,
#'              model_indicator = 0, mcmc_samples = 101, burnin = 20)
#' posterior_summary(fit)
#' @export
spfit <- function(data=NULL, formula=NULL, area=NULL, y=NULL, X=NULL, loc=NULL,
                  W, model_indicator, mcmc_samples, burnin=0L, adapt_rho=TRUE,
                  a_sigma2_prior=NULL, b_sigma2_prior=NULL,
                  a_tau2_prior=NULL, b_tau2_prior=NULL,
                  a_rho_prior=NULL, b_rho_prior=NULL,
                  beta_init=NULL, theta_init=NULL, sigma2_init=NULL,
                  tau2_init=NULL, rho_init=NULL, proposal_sd_init=NULL,
                  family="gaussian", trials=1, nb_size=NULL, offset=0,
                  beta_prior_mean=0, beta_prior_sd=10, center=TRUE,
                  pg_method=c("devroye","hybrid"), verbose=FALSE) {
  cl <- match.call()
  if (identical(family,"negbin")) family <- "negative_binomial"
  family <- match.arg(family,c("gaussian","bernoulli","binomial","negative_binomial"))
  z <- .sp_inputs(data,formula,area,y,X,loc,W)
  y <- z$y; X <- z$X; loc <- z$loc; W <- z$W
  .sp_integer(model_indicator,"model_indicator",0L)
  if (!model_indicator %in% 0:1) stop("model_indicator must be 0 or 1.")
  .sp_integer(mcmc_samples,"mcmc_samples",2L)
  .sp_integer(burnin,"burnin",0L)
  if (burnin >= mcmc_samples-1L) stop("burnin must leave at least one sampling update.")
  for (nm in c("adapt_rho","center","verbose")) {
    v <- get(nm)
    if (!is.logical(v)||length(v)!=1L||is.na(v)) stop(nm," must be TRUE or FALSE.")
  }
  if (model_indicator==1L && !.is_connected(W)) stop("The spatial model requires a connected graph.")
  if (!is.numeric(offset)||!length(offset)%in%c(1L,length(y))||any(!is.finite(offset)))
    stop("offset must be a finite scalar or have length(y) entries.")
  offset <- rep_len(offset,length(y)) + z$offset
  if (family!="gaussian") {
    if(any(!vapply(list(a_sigma2_prior,b_sigma2_prior,sigma2_init),is.null,logical(1))))
      stop("sigma2 controls apply only to the Gaussian family.")
    args <- list(y=y,X=X,loc=loc,W=W,model_indicator=model_indicator,
      mcmc_samples=mcmc_samples,family=family,trials=trials,nb_size=nb_size,
      offset=offset,burnin=burnin,adapt_rho=adapt_rho,
      beta_prior_mean=beta_prior_mean,beta_prior_sd=beta_prior_sd,center=center,
      pg_method=match.arg(pg_method),
      a_tau2_prior=a_tau2_prior,b_tau2_prior=b_tau2_prior,
      a_rho_prior=a_rho_prior,b_rho_prior=b_rho_prior,
      beta_init=beta_init,theta_init=theta_init,tau2_init=tau2_init,
      rho_init=rho_init,proposal_sd_init=proposal_sd_init)
    args <- args[!vapply(args,is.null,logical(1))]
    fit <- do.call(spfit_glmm,args)
  } else {
    if (!center) stop("Gaussian fits use the existing centering convention.")
    if (qr(X)$rank<ncol(X)) stop("X must have full column rank for the Gaussian flat prior.")
    vals <- list(a_sigma2_prior=a_sigma2_prior,b_sigma2_prior=b_sigma2_prior,
      a_tau2_prior=a_tau2_prior,b_tau2_prior=b_tau2_prior,
      a_rho_prior=a_rho_prior,b_rho_prior=b_rho_prior,
      sigma2_init=sigma2_init,tau2_init=tau2_init,proposal_sd_init=proposal_sd_init)
    for (nm in names(vals)) if(!is.null(vals[[nm]])) .sp_positive(vals[[nm]],nm)
    if (!is.null(rho_init) && (length(rho_init)!=1L||!is.finite(rho_init)||rho_init<=0||rho_init>=1))
      stop("rho_init must be in (0,1).")
    if (!is.null(beta_init) && (length(beta_init)!=ncol(X)||any(!is.finite(beta_init))))
      stop("beta_init must have one finite entry per design column.")
    if (!is.null(theta_init) && (length(theta_init)!=nrow(W)||any(!is.finite(theta_init))))
      stop("theta_init must have one finite entry per location.")
    fit <- SpThreshold(mcmc_samples=as.integer(mcmc_samples),y=y-offset,X=X,
      loc=as.integer(loc-1L),model_indicator=as.integer(model_indicator),W=W,
      proposal_sd_init=proposal_sd_init,a_sigma2_prior=a_sigma2_prior,
      b_sigma2_prior=b_sigma2_prior,a_tau2_prior=a_tau2_prior,
      b_tau2_prior=b_tau2_prior,a_rho_prior=a_rho_prior,b_rho_prior=b_rho_prior,
      beta_init=beta_init,theta_init=theta_init,sigma2_init=sigma2_init,
      tau2_init=tau2_init,rho_init=rho_init,burnin=as.integer(burnin),
      adapt_rho=adapt_rho,verbose=verbose)
  }
  rownames(fit$beta) <- colnames(X)
  fit$call <- cl; fit$model_indicator <- as.integer(model_indicator)
  fit$mcmc_samples <- as.integer(mcmc_samples); fit$burnin <- as.integer(burnin)
  fit$family <- family; fit$initial_state_included <- TRUE
  fit$design <- list(y=y,X=X,loc=loc,W=W,offset=offset)
  class(fit) <- unique(c("spfit",class(fit)))
  fit
}

.sp_inputs <- function(data,formula,area,y,X,loc,W) {
  df <- !vapply(list(data,formula,area),is.null,logical(1))
  ve <- !vapply(list(y,X,loc),is.null,logical(1))
  if (any(df)&&any(ve)) stop("Use either (data, formula, area) or (y, X, loc).")
  if (!(all(df)||all(ve))) stop("Supply all of (data, formula, area) or all of (y, X, loc).")
  off <- 0
  if(all(df)) {
    if(!is.data.frame(data)||!inherits(formula,"formula")||length(area)!=1L||!is.character(area)||!area%in%names(data))
      stop("Invalid data, formula or area column.")
    mf <- stats::model.frame(formula,data=data,na.action=stats::na.fail)
    y <- stats::model.response(mf); X <- stats::model.matrix(formula,mf)
    loc <- data[[area]]; fo <- stats::model.offset(mf)
    if(!is.null(fo)) off <- fo
  }
  if(!is.numeric(y)||!is.null(dim(y))||!length(y)||any(!is.finite(y)))
    stop("y must be a finite numeric response vector; binomial counts use trials separately.")
  if(!is.matrix(X)||!is.numeric(X)||nrow(X)!=length(y)||ncol(X)<1L||any(!is.finite(X)))
    stop("X must be a finite numeric matrix with length(y) rows.")
  if(is.null(colnames(X))) colnames(X) <- ifelse(vapply(seq_len(ncol(X)),function(j) all(X[,j]==1),logical(1)),"(Intercept)",paste0("beta",seq_len(ncol(X))))
  if(anyDuplicated(colnames(X))) colnames(X) <- make.unique(colnames(X))
  if(inherits(W,"Matrix")) W <- as.matrix(W)
  if(!is.matrix(W)||!is.numeric(W)||nrow(W)!=ncol(W)||nrow(W)<1L||any(!is.finite(W))||
     any(W<0)||any(diag(W)!=0)||!isSymmetric(W,tol=1e-10))
    stop("W must be a finite symmetric nonnegative square matrix with zero diagonal.")
  if(!is.numeric(loc)||length(loc)!=length(y)||any(!is.finite(loc))||any(loc!=floor(loc))||any(!loc%in%seq_len(nrow(W))))
    stop("loc must contain integer location IDs in 1:nrow(W), one per response.")
  if(!all(seq_len(nrow(W))%in%loc)) stop("Every row of W must have at least one observed location ID.")
  list(y=as.numeric(y),X=X,loc=as.integer(loc),W=W,offset=off)
}
.sp_integer <- function(x,name,lower) {
  if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x!=floor(x)||x<lower||x>.Machine$integer.max)
    stop(name," must be an integer >= ",lower,".")
}
.sp_positive <- function(x,name) {
  if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<=0) stop(name," must be finite and positive.")
}
.is_connected <- function(W) {
  n <- nrow(W); seen <- rep(FALSE,n); seen[1L] <- TRUE; queue <- 1L
  while(length(queue)) {
    cur <- queue[1L]; queue <- queue[-1L]
    add <- which(W[cur,]!=0 & !seen); seen[add] <- TRUE; queue <- c(queue,add)
  }
  all(seen)
}
