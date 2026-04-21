\name{SpThreshold}
\alias{SpThreshold}
\title{Bayesian Multilevel Model for Areal Data with or without a Spatial Random Effect}
\description{
Fits a Bayesian multilevel model for repeatedly measured areal spatial data under either a Leroux conditional autoregressive (CAR) prior on the area-level random effects (spatial model) or an independent and identically distributed prior (nonspatial model). Posterior sampling proceeds via Markov chain Monte Carlo.
}
\usage{
SpThreshold(mcmc_samples, y, X, loc, model_indicator,
            W = NULL, proposal_sd_init = NULL,
            a_sigma2_prior = NULL, b_sigma2_prior = NULL,
            a_tau2_prior = NULL, b_tau2_prior = NULL,
            a_rho_prior = NULL, b_rho_prior = NULL,
            beta_init = NULL, theta_init = NULL,
            sigma2_init = NULL, tau2_init = NULL, rho_init = NULL,
            burnin = 0, adapt_rho = TRUE)
}
\arguments{
\item{mcmc_samples}{Integer number of MCMC iterations (including burn-in).}

\item{y}{Length-\eqn{N} numeric vector of outcomes, where \eqn{N = nm}.}

\item{X}{\eqn{N} by \eqn{p} design matrix for the fixed effects, typically containing an intercept column and one or more covariate columns.}

\item{loc}{Length-\eqn{N} integer vector identifying the spatial area for each observation, with values in \eqn{0, 1, \ldots, n-1}.}

\item{model_indicator}{Integer, either \code{0} (nonspatial) or \code{1} (spatial Leroux CAR).}

\item{W}{Required when \code{model_indicator = 1}. An \eqn{n} by \eqn{n} symmetric binary adjacency matrix with zero diagonal.}

\item{proposal_sd_init}{Initial standard deviation of the logit-scale random-walk proposal for \eqn{\rho}. Default: \code{0.30}.}

\item{a_sigma2_prior, b_sigma2_prior}{Shape and rate of the inverse-gamma prior on \eqn{\sigma^2}. Defaults: \code{0.01, 0.01}.}

\item{a_tau2_prior, b_tau2_prior}{Shape and rate of the inverse-gamma prior on \eqn{\tau^2}. Defaults: \code{0.01, 0.01}.}

\item{a_rho_prior, b_rho_prior}{Shape parameters of the beta prior on \eqn{\rho}. Defaults: \code{1, 1} (Uniform\eqn{(0, 1)}).}

\item{beta_init, theta_init, sigma2_init, tau2_init, rho_init}{Initial values for MCMC. Defaults: \code{0} vector for \code{beta_init} and \code{theta_init}, \code{1.00} for \code{sigma2_init} and \code{tau2_init}, \code{0.50} for \code{rho_init}.}

\item{burnin}{Integer number of initial iterations during which the \eqn{\rho} proposal standard deviation is adapted (if \code{adapt_rho = TRUE}). Should be less than \code{mcmc_samples}.}

\item{adapt_rho}{Logical. If \code{TRUE} (default), the \eqn{\rho} proposal standard deviation is adapted during burn-in using a Robbins-Monro scheme targeting an acceptance rate of \code{0.234}.}
}
\value{
A list with components
\item{\code{beta}}{\eqn{p} by \code{mcmc_samples} matrix of \eqn{\boldsymbol{\beta}} posterior samples.}
\item{\code{theta}}{\eqn{n} by \code{mcmc_samples} matrix of \eqn{\boldsymbol{\theta}} posterior samples.}
\item{\code{sigma2}}{Length-\code{mcmc_samples} vector of \eqn{\sigma^2} posterior samples.}
\item{\code{tau2}}{Length-\code{mcmc_samples} vector of \eqn{\tau^2} posterior samples.}
\item{\code{rho}}{(Spatial model only.) Length-\code{mcmc_samples} vector of \eqn{\rho} posterior samples.}
\item{\code{accept_rho}}{(Spatial model only.) Total number of accepted \eqn{\rho} proposals.}
\item{\code{final_proposal_sd}}{(Spatial model only.) Final proposal standard deviation after adaptation.}
}
\details{
Full details of the model specification, prior distributions, and MCMC sampling algorithm are given in the \code{SpThreshold_Model_Details} vignette.
}
\seealso{
\code{\link{m_star}} for the closed-form replication threshold; \code{\link{create_random_W}} and \code{\link{create_W_from_shapefile}} for adjacency matrix construction.
}
