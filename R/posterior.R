#' Extract retained posterior draws
#'
#' The package stores the initial state in column one. This function discards
#' that state and the requested number of sampling updates, then thins by update
#' number. For example, burnin=20 and thin=2 retains updates 22, 24, and so on.
#' @param fit A fit returned by \code{spfit} or \code{spfit_glmm}.
#' @param burnin Number of initial updates to discard. NULL uses the fit metadata,
#'   or zero for an older fit without metadata.
#' @param thin Retain every thin-th update after burn-in.
#' @param include_theta Include columns for the random effects.
#' @return A data frame with one row per retained draw. Attributes \code{updates},
#'   \code{burnin}, and \code{thin} describe the retained update numbers.
#' @export
posterior_draws <- function(fit, burnin=NULL, thin=1L, include_theta=FALSE) {
  if (!is.list(fit)||is.null(fit$beta)||!is.matrix(fit$beta))
    stop("fit must contain a beta matrix with parameters in rows and states in columns.")
  nstate <- ncol(fit$beta)
  if (nstate<2L) stop("fit contains no sampling updates after its initial state.")
  if (is.null(burnin)) burnin <- if(is.null(fit$burnin)) 0L else fit$burnin
  .sp_integer(burnin,"burnin",0L); .sp_integer(thin,"thin",1L)
  if(!is.logical(include_theta)||length(include_theta)!=1L||is.na(include_theta))
    stop("include_theta must be TRUE or FALSE.")
  if(burnin+thin>nstate-1L) stop("burnin and thin leave no retained sampling updates.")
  updates <- seq.int(burnin+thin,nstate-1L,by=thin)
  cols <- updates+1L
  bn <- rownames(fit$beta)
  if(is.null(bn)) bn <- paste0("beta[",seq_len(nrow(fit$beta)),"]")
  out <- as.data.frame(t(fit$beta[,cols,drop=FALSE]))
  names(out) <- make.unique(bn)
  for(nm in c("sigma2","tau2","rho")) if(!is.null(fit[[nm]])) {
    v <- as.numeric(fit[[nm]])
    if(length(v)!=nstate) stop(nm," has an inconsistent number of states.")
    if(nm%in%names(out)) stop("A coefficient name conflicts with parameter name: ",nm)
    out[[nm]] <- v[cols]
  }
  if(include_theta) {
    if(is.null(fit$theta)||!is.matrix(fit$theta)||ncol(fit$theta)!=nstate)
      stop("fit must contain a theta matrix with the same number of states as beta.")
    tn <- paste0("theta[",seq_len(nrow(fit$theta)),"]")
    if(any(tn%in%names(out))) stop("Coefficient names conflict with theta column names.")
    th <- as.data.frame(t(fit$theta[,cols,drop=FALSE])); names(th) <- tn
    out <- cbind(out,th)
  }
  rownames(out) <- NULL
  attr(out,"updates") <- updates; attr(out,"burnin") <- burnin; attr(out,"thin") <- thin
  out
}

#' Summarize the marginal posterior from MCMC draws
#'
#' Computes sample means, sample variances and quantiles of the retained draws.
#' No conditional plug-in approximation or Rao--Blackwellization is applied.
#' @param fit,burnin,thin,include_theta See \code{\link{posterior_draws}}.
#' @param probs Two quantile probabilities for lower and upper summaries.
#' @return A data frame with parameter, mean, variance, sd, lower, upper,
#'   n_draws, ess, and mcse_mean. ESS uses an initial-positive, monotone paired
#'   autocorrelation estimate within this chain; it is not a multi-chain
#'   convergence diagnostic. Undefined ESS is NA for constant draws.
#' @export
posterior_summary <- function(fit, burnin=NULL, thin=1L, include_theta=FALSE,
                              probs=c(.025,.975)) {
  if(!is.numeric(probs)||length(probs)!=2L||any(!is.finite(probs))||
     any(probs<0|probs>1)||probs[1]>=probs[2]) stop("probs must be two increasing probabilities in [0,1].")
  draws <- posterior_draws(fit,burnin,thin,include_theta)
  if(nrow(draws)<2L) stop("At least two retained draws are required for a variance summary.")
  vals <- lapply(draws,function(x) {
    v <- stats::var(x); ess <- .sp_ess(x)
    q <- stats::quantile(x,probs=probs,names=FALSE)
    c(mean=mean(x),variance=v,sd=sqrt(v),lower=q[1],upper=q[2],
      n_draws=length(x),ess=ess,mcse_mean=if(is.finite(ess)) sqrt(v/ess) else NA_real_)
  })
  out <- data.frame(parameter=names(draws),do.call(rbind,vals),row.names=NULL,check.names=FALSE)
  attr(out,"burnin") <- attr(draws,"burnin"); attr(out,"thin") <- attr(draws,"thin")
  attr(out,"probs") <- probs
  out
}

.sp_ess <- function(x) {
  n <- length(x)
  if(n<4L||!is.finite(stats::var(x))||stats::var(x)==0) return(NA_real_)
  ac <- as.numeric(stats::acf(x,lag.max=min(n-1L,max(100L,ceiling(10*sqrt(n)))),
                            plot=FALSE,demean=TRUE)$acf)
  k <- length(ac)%/%2L
  pairs <- ac[seq.int(1L,2L*k,by=2L)]+ac[seq.int(2L,2L*k,by=2L)]
  bad <- which(pairs<=0)
  if(length(bad)) pairs <- pairs[seq_len(bad[1]-1L)]
  if(!length(pairs)) return(as.numeric(n))
  tau <- -1+2*sum(cummin(pairs))
  min(n,n/max(tau,1))
}

#' @rdname posterior_summary
#' @param object A fitted model.
#' @param ... Passed to \code{posterior_summary}.
#' @export
summary.spfit <- function(object, ...) posterior_summary(object,...)
