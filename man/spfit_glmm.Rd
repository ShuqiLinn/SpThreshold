\name{spfit_glmm}
\alias{spfit_glmm}
\title{Fit Bernoulli, binomial, or negative-binomial spatial mixed models}
\usage{
spfit_glmm(y, X, loc, W, model_indicator, mcmc_samples,
  family = c("bernoulli", "binomial", "negative_binomial"),
  trials = 1, nb_size = NULL, offset = 0, burnin = 0L,
  adapt_rho = TRUE, beta_prior_mean = 0, beta_prior_sd = 10,
  center = TRUE, pg_method = "devroye",
  a_tau2_prior = 0.01, b_tau2_prior = 0.01,
  a_rho_prior = 1, b_rho_prior = 1,
  beta_init = NULL, theta_init = NULL, tau2_init = 1,
  rho_init = 0.5, proposal_sd_init = 0.30)
}
\arguments{
\item{y}{Nonnegative integer response vector. For Bernoulli, entries are 0 or 1; for binomial, entries are success counts.}
\item{X}{Numeric design matrix with one row per observation. Include a column of ones for an intercept. Multiple covariates are supported.}
\item{loc}{Integer location indices in \code{1:nrow(W)}, one per observation. Replication may be unequal and locations may have no observations.}
\item{W}{Symmetric, nonnegative adjacency matrix with zero diagonal. A connected graph is required for the spatial model.}
\item{model_indicator}{1 for Leroux CAR random effects or 0 for iid Gaussian random effects.}
\item{mcmc_samples}{Number of returned states, including the initial state in column 1. There are \code{mcmc_samples - 1} transitions.}
\item{family}{\code{"bernoulli"}, \code{"binomial"}, or \code{"negative_binomial"}.}
\item{trials}{Positive integer number of binomial trials, scalar or one per observation. Must be 1 for Bernoulli; unused for negative binomial.}
\item{nb_size}{Fixed positive integer negative-binomial dispersion parameter, required for negative binomial. It is not estimated.}
\item{offset}{Known additive offset in the linear predictor, scalar or one per observation.}
\item{burnin}{Number of transitions during which the spatial rho proposal may adapt. No samples are removed internally.}
\item{adapt_rho}{Whether to adapt the rho proposal during burn-in.}
\item{beta_prior_mean,beta_prior_sd}{Independent normal coefficient prior means and standard deviations, scalar or one per coefficient. Defaults to \eqn{N(0,10^2)}. A standard deviation of \code{Inf} requests an improper flat prior; full column rank is then required, and the user must establish posterior propriety, especially under separation.}
\item{center}{If \code{TRUE}, fit the explicit constraint \eqn{\sum_i\theta_i=0}. Requires an intercept column. If \code{FALSE}, use the full, proper Gaussian random-effect prior. This is a model choice, not a display transformation.}
\item{pg_method}{\code{"devroye"}, the only supported value, uses the \pkg{pgdraw} package's exact accept-reject sampler for integer Polya-Gamma shapes.}
\item{a_tau2_prior,b_tau2_prior}{Inverse-gamma shape and scale for the random-effect variance: density proportional to \eqn{t^{-a-1}\exp(-b/t)}.}
\item{a_rho_prior,b_rho_prior}{Beta prior shape parameters for rho.}
\item{beta_init,theta_init}{Initial coefficient and random-effect vectors, default zero. With \code{center = TRUE}, theta must sum to zero.}
\item{tau2_init,rho_init,proposal_sd_init}{Initial variance, rho in \eqn{(0,1)}, and standard deviation of the logit-rho proposal.}
}
\description{
Fits non-Gaussian multilevel models using Polya-Gamma augmentation, with
\pkg{pgdraw} random-number generation and compiled Gaussian conditional updates.
}
\details{
The linear predictor is \eqn{\eta_j=X_j\beta+\theta_{loc_j}+offset_j}.
Bernoulli and binomial use \eqn{\operatorname{logit}(p_j)=\eta_j}.
Negative binomial uses \eqn{\log(\mu_j)=\eta_j}, with variance
\eqn{\mu_j+\mu_j^2/r} for \eqn{r=nb\_size}; its Polya-Gamma tilt is
\eqn{\eta_j-\log(r)}, shape \eqn{y_j+r}, and linear coefficient
\eqn{(y_j-r)/2}. Binomial shape and linear coefficient are \eqn{trials_j}
and \eqn{y_j-trials_j/2}.

The proper coefficient prior is specific to these new families and avoids
assuming posterior propriety under separated binary data. It does not modify
the Gaussian model's coefficient prior. A sufficiently justified flat prior
can be requested explicitly through \code{beta_prior_sd = Inf}.

For the centered model, the Gaussian conditional draw is conditioned exactly
on \eqn{1^T\theta=0}. Ordinary mean subtraction is not used, because the
Polya-Gamma weights are generally unequal. The tau2 update uses dimension
\eqn{n-1}, and the rho density uses the determinant on that same contrast
space. For \code{center = FALSE}, both use dimension \eqn{n} and the full
Leroux precision. These internal prior normalizations do not define or alter
the full-field Gaussian variance calibration used in threshold calculations.

No residual Gaussian variance parameter is present. Threshold formulas for
the Gaussian response model do not become marginal bounds for these families.
Use the posterior draws to assess their behavior empirically. The negative
binomial dispersion is fixed: a study that requires estimating it needs an
additional model and sampler, which this function does not supply.

The Polya-Gamma draws use \pkg{pgdraw} (version 1.1 or later), which requires
positive integer shapes. Binomial trial counts must be positive integers;
negative-binomial size must be a fixed positive integer. Shapes must not exceed
\code{.Machine$integer.max}. The sampling work can increase for large counts.

No burn-in removal or thinning is performed. After \code{b} transitions of
burn-in, the first eligible update is stored in column \code{b + 2}. For example,
with 1000 burn-in transitions, use columns \code{1002:mcmc_samples}, optionally
thinning from those updates.
}
\value{
A list containing \code{beta} and \code{theta} matrices (one state per column),
\code{tau2} and \code{rho} vectors, acceptance and proposal information,
and family, link, offset, binomial trials or fixed negative-binomial size,
coefficient prior, centering, Polya-Gamma method, and iteration metadata.
The initial values are the first state. The iid fit has rho identically zero.
There is no \code{sigma2} component.
}
\references{
Polson, N. G., Scott, J. G., and Windle, J. (2013).
Bayesian inference for logistic models using Polya-Gamma latent variables.
\emph{Journal of the American Statistical Association}, 108, 1339--1349.
\doi{10.1080/01621459.2013.829001}.

pgdraw package documentation:
\url{https://cran.r-project.org/package=pgdraw}.
}
\examples{
set.seed(18)
W <- matrix(c(0,1,1,0), 2, 2)
loc <- rep(1:2, c(8,12))
x <- rnorm(length(loc))
y <- rbinom(length(loc), 1, plogis(0.3 + x))
fit <- spfit_glmm(y, cbind(1,x), loc, W, 0, 30,
                  family = "bernoulli", burnin = 10)
stopifnot(max(abs(colSums(fit$theta))) < 1e-10)
}
