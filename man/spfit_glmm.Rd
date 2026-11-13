\name{spfit_glmm}
\alias{spfit_glmm}
\title{Fit Non-Gaussian Multilevel Models}

\usage{
spfit_glmm(
  y, X, loc, W, model_indicator, mcmc_samples,
  family = c("bernoulli", "binomial", "negative_binomial"),
  trials = 1, nb_size = NULL, offset = 0, burnin = 0L,
  adapt_rho = TRUE, beta_prior_mean = 0, beta_prior_sd = 10,
  center = TRUE, pg_method = "devroye",
  a_tau2_prior = 0.01, b_tau2_prior = 0.01,
  a_rho_prior = 1, b_rho_prior = 1,
  beta_init = NULL, theta_init = NULL, tau2_init = 1,
  rho_init = 0.5, proposal_sd_init = 0.30
)
}

\arguments{
\item{y}{Nonnegative integer response vector. Bernoulli responses
must be zero or one; binomial responses cannot exceed \code{trials}.}

\item{X}{Numeric design matrix with one row per observation.
Include a column of ones for an intercept.}

\item{loc}{Integer location indices from \code{1} to \code{nrow(W)},
one per observation. Replication may be unequal, and locations
may have no observations.}

\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.
A connected graph is required for the spatial model.}

\item{model_indicator}{1 for the Leroux spatial model or 0 for
the nonspatial Gaussian random-intercept model.}

\item{mcmc_samples}{Number of stored states, including the initial
state. Must be at least 2.}

\item{family}{Response family: \code{"bernoulli"},
\code{"binomial"}, or \code{"negative_binomial"}.}

\item{trials}{Positive integer binomial trial counts, supplied as
a scalar or one value per observation. Must equal 1 for Bernoulli
and negative binomial responses.}

\item{nb_size}{Fixed positive integer negative binomial size parameter.
Required for \code{family = "negative_binomial"}; otherwise
\code{NULL}. This parameter is not estimated.}

\item{offset}{Known additive offset in the linear predictor,
supplied as a scalar or one value per observation.}

\item{burnin}{Number of initial sampling updates treated as burn-in.
The spatial correlation proposal may adapt during these updates.
No states are removed from the returned fit.}

\item{adapt_rho}{Whether to adapt the spatial correlation proposal
during burn-in.}

\item{beta_prior_mean, beta_prior_sd}{Means and standard deviations
of independent normal coefficient priors, supplied as scalars or
one value per coefficient. Defaults to \eqn{N(0,10^2)}.
A standard deviation of \code{Inf} specifies a flat prior.
Flat priors require full column rank in \code{X} and a proper posterior.}

\item{center}{Whether to constrain the random intercepts to sum to zero.
Defaults to \code{TRUE}, which requires an intercept column in \code{X}.
With \code{FALSE}, the Gaussian random-intercept prior is unconstrained.}

\item{pg_method}{Polya-Gamma sampling method. The only supported
value is \code{"devroye"}, using \pkg{pgdraw}.}

\item{a_tau2_prior, b_tau2_prior}{Inverse-gamma shape and scale
parameters for the random-intercept variance, with density
proportional to \eqn{t^{-a-1}\exp(-b/t)}.}

\item{a_rho_prior, b_rho_prior}{Beta prior shape parameters for
the spatial correlation parameter.}

\item{beta_init, theta_init}{Initial coefficient and random-intercept
vectors. Defaults to zero vectors. When \code{center = TRUE},
\code{theta_init} must sum to zero.}

\item{tau2_init}{Initial positive random-intercept variance parameter.}

\item{rho_init}{Initial spatial correlation parameter,
strictly between zero and one.}

\item{proposal_sd_init}{Initial positive standard deviation
of the proposal on the logit scale of spatial correlation.}
}

\description{
Fits Bernoulli, binomial, and negative binomial multilevel models
using Polya-Gamma augmentation and compiled Gaussian conditional
updates, with spatial or nonspatial Gaussian random intercepts.
}

\details{
The linear predictor is
\eqn{\eta_j = X_j\beta + \theta_{\mathrm{loc}_j} + \mathrm{offset}_j}.
Bernoulli and binomial models use
\eqn{\operatorname{logit}(p_j) = \eta_j}.
The negative binomial model uses \eqn{\log(\mu_j) = \eta_j},
with variance \eqn{\mu_j + \mu_j^2/r}, where
\eqn{r} is \code{nb_size}.

The returned fit includes the initial state and all sampling updates.
Use \code{\link{posterior_draws}} or
\code{\link{posterior_summary}} to apply burn-in and thinning.
}

\value{
An object of class \code{spfit} containing:
\describe{
\item{\code{beta}}{Coefficient draws, with coefficients in rows
and stored states in columns.}

\item{\code{theta}}{Random-intercept draws, with locations in rows
and stored states in columns.}

\item{\code{tau2}}{Draws of the model-specific random-intercept
variance parameter.}

\item{\code{rho}}{Spatial correlation draws, included only
for the spatial model.}

\item{\code{accept_rho}}{Number of accepted spatial correlation
proposals. Zero for the nonspatial model.}

\item{\code{final_proposal_sd}}{Final proposal standard deviation
on the logit scale of spatial correlation.}
}

The object also records the response family, link, offsets,
trial counts or fixed negative binomial size, coefficient priors,
centering setting, design, and sampling settings.
There is no \code{sigma2} component.
}