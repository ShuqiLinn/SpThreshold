\name{SpThreshold}
\alias{SpThreshold}
\title{Fit Gaussian Multilevel Models}

\description{
Fits Gaussian multilevel models with Leroux spatial or independent
Gaussian random intercepts using MCMC. Multiple covariates and
unequal within-location sample sizes are supported.
}

\usage{
SpThreshold(
  mcmc_samples, y, X, loc, model_indicator,
  W = NULL, proposal_sd_init = NULL,
  a_sigma2_prior = NULL, b_sigma2_prior = NULL,
  a_tau2_prior = NULL, b_tau2_prior = NULL,
  a_rho_prior = NULL, b_rho_prior = NULL,
  beta_init = NULL, theta_init = NULL,
  sigma2_init = NULL, tau2_init = NULL, rho_init = NULL,
  burnin = 0L, adapt_rho = TRUE, verbose = TRUE
)
}

\arguments{
\item{mcmc_samples}{Number of MCMC iterations, counting initialization
as the first iteration.}

\item{y}{Numeric response vector.}

\item{X}{Full-rank numeric design matrix with one row per observation
and a column of ones for the intercept.}

\item{loc}{Integer location indices starting at zero, one per
observation. At least two locations are required.}

\item{model_indicator}{1 for spatial or 0 for nonspatial
random intercepts.}

\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.
Required for the spatial model.}

\item{proposal_sd_init}{Initial standard deviation of the
logit-scale proposal for \eqn{\rho}. \code{NULL} uses 0.30.}

\item{a_sigma2_prior, b_sigma2_prior}{Inverse-gamma shape and scale
parameters for the residual variance.
Each defaults to 0.01 when \code{NULL}.}

\item{a_tau2_prior, b_tau2_prior}{Inverse-gamma shape and scale
parameters for the random-intercept variance.
Each defaults to 0.01 when \code{NULL}.}

\item{a_rho_prior, b_rho_prior}{Beta prior shape parameters for
\eqn{\rho}. Each defaults to 1 when \code{NULL}.}

\item{beta_init, theta_init}{Initial coefficient and random-intercept
vectors. \code{NULL} uses zero vectors.}

\item{sigma2_init, tau2_init}{Initial variance parameters.
\code{NULL} uses 1 for each.}

\item{rho_init}{Initial spatial correlation in \eqn{(0,1)}.
\code{NULL} uses 0.5.}

\item{burnin}{Number of initial updates during which the
\eqn{\rho} proposal may adapt. Draws are not discarded internally.}

\item{adapt_rho}{Whether to adapt the \eqn{\rho} proposal during
burn-in, targeting an acceptance rate of 0.234.}

\item{verbose}{Whether to print sampling progress.}
}

\details{
Regression coefficients have flat priors. Inverse-gamma priors
use the density proportional to \eqn{t^{-a-1}\exp(-b/t)}.

Random intercepts are represented with zero sum, with their mean
absorbed into the intercept. Under the flat intercept prior,
this preserves the slope and covariance parameter posteriors
of the unrestricted model.

For the standard interface with location indices starting at one,
use \code{\link{spfit}}.
}

\value{
A list containing coefficient and random-intercept matrices
\code{beta} and \code{theta}, with one stored state per column,
and variance parameter vectors \code{sigma2} and \code{tau2}.

Spatial fits also return \code{rho}, the accepted proposal count
\code{accept_rho}, and the final proposal standard deviation
\code{final_proposal_sd}.

The first stored state is the initialization. To summarize draws,
use \code{\link{posterior_summary}} with \code{burnin}
specified explicitly.
}