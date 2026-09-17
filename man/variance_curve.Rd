\name{variance_curve}
\alias{variance_curve}
\title{Exact and Leading-Order Conditional Variance Curves}
\usage{
variance_curve(W, sigma2, rho, m, xbar = NULL, X = NULL,
  loc = NULL, tau0_2 = NULL, tau1_2 = NULL, tau2 = NULL,
  within_ss = NULL, calibration = c("full", "common"))
}
\arguments{
\item{W}{Symmetric nonnegative adjacency matrix with zero diagonal.}
\item{sigma2}{Positive residual variance, shared between the two specified models.}
\item{rho}{Leroux parameter in \eqn{[0,1)}.}
\item{m}{Positive replication values at which to evaluate the balanced-design curve.}
\item{xbar}{Location means for one covariate, or a matrix with one covariate per column. Supply this or \code{X} with \code{loc}.}
\item{X}{A covariate vector or matrix for a balanced observed design. No intercept is added. Variables need not be centered or standardized.}
\item{loc}{Integer location indices \code{1:nrow(W)}, one per row of \code{X}. Every location must be represented equally often.}
\item{tau0_2,tau1_2,tau2,calibration}{Variance-scale inputs and calibration; see \code{\link{variance_calibration}}. \code{tau2} is a deprecated spatial-scale alias.}
\item{within_ss}{Within-location sum of squared covariate deviations per unit of balanced replication, one value per covariate. See Details.}
}
\description{Evaluates the exact Gaussian conditional variance difference at specified parameters and its leading-order approximation as balanced replication varies.}
\details{
The target is the variance of one coefficient conditional on all other regression coefficients and the covariance parameters, with Gaussian random effects integrated out and a flat coefficient prior. For covariance \eqn{V}, this is
\eqn{1/(\boldsymbol{x}_j^T V^{-1}\boldsymbol{x}_j)}.
It is not the \eqn{j}th diagonal of \eqn{(X^T V^{-1}X)^{-1}}, nor a marginal posterior variance integrating over covariance parameters. Each matrix column is evaluated separately with this same conditional target.

Write \eqn{D_i=(U^T\bar{x})_i^2}, \eqn{s=W_{\mathrm{within}}/m}, and \eqn{c_i=1-\rho+\rho\lambda_i}. The exact precisions are
\deqn{I_1(m)=\frac{ms}{\sigma^2}+\sum_i\frac{m c_i D_i}{\sigma^2c_i+m\tau_1^2},
\qquad I_0(m)=\frac{ms}{\sigma^2}+\frac{m\sum_iD_i}{\sigma^2+m\tau_0^2}.}
The signed relative difference is \eqn{I_0/I_1-1}, corresponding to
\eqn{(\operatorname{Var}_1-\operatorname{Var}_0)/\operatorname{Var}_0}.
For \eqn{s>0}, its leading signed term is
\deqn{\frac{\sigma^2\sum_iD_i(1/\tau_0^2-c_i/\tau_1^2)}{ms}.}
The absolute exact difference is bounded by \eqn{A/m}, where
\deqn{A=\frac{\sigma^2}{s}\sum_iD_i\left|1/\tau_0^2-c_i/\tau_1^2\right|.}
This bound follows by bounding each between-location precision difference by its absolute limit. It does not assume monotonicity of their signed sum.

For raw \code{X}, the location means and \eqn{s} are computed from the balanced input design; those quantities are held fixed as \code{m} varies. For \code{xbar} without \code{within_ss}, the population-standardized template is assumed: centered covariates with population variance one, giving \eqn{s=n-\sum_iD_i}. To use any other scale, specify \code{within_ss} explicitly. A varying-m curve describes this fixed template, not newly generated covariates at each m.

When \eqn{s=0}, the exact curve remains evaluable, but no decay bound or leading \eqn{1/m} approximation is supplied. Identical covariance action on the covariate gives zero difference even in this case, including calibrated \eqn{\rho=0}.
}
\value{A long data frame with columns \code{m}, \code{coefficient}, \code{variance_spatial}, \code{variance_iid}, \code{relative_signed}, \code{relative_difference}, \code{leading_signed}, \code{leading_difference}, \code{tail_bound}, and \code{within_ss_per_m}. Attributes \code{calibration} and \code{target} document the comparison.}
\examples{
W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
variance_curve(W, sigma2 = 0.5, rho = 0.8, m = c(2, 10, 100),
               xbar = c(-0.5, 0, 0.5), tau0_2 = 0.5)
}
\seealso{\code{\link{variance_calibration}}, \code{\link{m_star}}, \code{\link{variance_comparison}}}
