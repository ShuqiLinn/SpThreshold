#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
arma::vec glmm_precision_draw(const arma::mat& A, const arma::vec& h,
                              bool constrain = false) {
    arma::mat U;
    if (!arma::chol(U, A)) Rcpp::stop("Non-positive Gaussian conditional precision.");
    arma::vec mu = arma::solve(arma::trimatu(U),
                              arma::solve(arma::trimatl(U.t()), h));
    arma::vec z(h.n_elem);
    for (arma::uword j = 0; j < z.n_elem; ++j) z(j) = R::rnorm(0.0, 1.0);
    arma::vec draw = mu + arma::solve(arma::trimatu(U), z);
    if (constrain) {
        // Exact Gaussian conditioning on 1' theta = 0. Subtracting the
        // ordinary mean is incorrect with unequal Polya-Gamma weights.
        arma::vec one(h.n_elem, arma::fill::ones);
        arma::vec a1 = arma::solve(arma::trimatu(U),
                                  arma::solve(arma::trimatl(U.t()), one));
        draw -= a1 * (arma::sum(draw) / arma::sum(a1));
    }
    return draw;
}

double glmm_log_rho(double rho, const arma::vec& lambda,
                    const arma::vec& theta, const arma::mat& L,
                    double tau2, bool center, double a, double b) {
    if (rho <= 0.0 || rho >= 1.0) return R_NegInf;
    double logdet = 0.0;
    for (arma::uword k = center ? 1 : 0; k < lambda.n_elem; ++k)
        logdet += std::log(1.0 - rho + rho * lambda(k));
    double quad = (1.0-rho)*arma::dot(theta, theta) +
                   rho*arma::dot(theta, L*theta);
    // Includes the Jacobian for the logit-rho random walk.
    return 0.5*logdet - 0.5*quad/tau2 + a*std::log(rho) + b*std::log1p(-rho);
}
}

// [[Rcpp::export]]
Rcpp::List sp_glmm_pg_cpp(arma::vec y, arma::mat X, arma::uvec loc,
                         arma::mat W, int model_indicator, int mcmc_samples,
                         int burnin, bool adapt_rho, arma::vec pg_shape,
                         arma::vec kappa, arma::vec offset_pg,
                         arma::vec beta_mean, arma::vec beta_precision,
                         arma::vec beta, arma::vec theta, double tau2,
                         double rho, double proposal_sd,
                         double a_tau, double b_tau, double a_rho,
                         double b_rho, bool center) {
    const arma::uword N = y.n_elem, n = W.n_rows, p = X.n_cols;
    arma::mat L = arma::diagmat(arma::sum(W, 1)) - W;
    arma::vec lambda;
    if (!arma::eig_sym(lambda, L)) Rcpp::stop("Graph eigendecomposition failed.");
    lambda.transform([](double x) { return std::max(0.0, x); });
    if (model_indicator == 0) rho = 0.0;
    arma::mat Q = (1.0-rho)*arma::eye(n,n) + rho*L;
    Rcpp::Environment pg_ns = Rcpp::Environment::namespace_env("pgdraw");
    Rcpp::Function pg = pg_ns["pgdraw"];
    arma::mat beta_samples(p,mcmc_samples), theta_samples(n,mcmc_samples);
    arma::vec tau_samples(mcmc_samples), rho_samples(mcmc_samples);
    beta_samples.col(0) = beta; theta_samples.col(0) = theta;
    tau_samples(0) = tau2; rho_samples(0) = rho;
    int accepted = 0, batch_accepted = 0;
    const double rank = center ? n-1.0 : n;
    for (int it = 1; it < mcmc_samples; ++it) {
        arma::vec tilt = X*beta + theta.elem(loc) + offset_pg;
        // Synchronize R's seed around the R-level PG routine. The remainder
        // of this loop draws through R's RNG in C++; preserve one RNG stream.
        PutRNGstate();
        Rcpp::NumericVector om_r = pg(Rcpp::Named("b") = pg_shape,
                                     Rcpp::Named("c") = tilt);
        GetRNGstate();
        arma::vec omega = Rcpp::as<arma::vec>(om_r);
        if (omega.n_elem != N || !omega.is_finite() || arma::any(omega <= 0.0))
            Rcpp::stop("Polya-Gamma sampler returned invalid draws.");
        arma::mat weighted_X = X.each_col() % omega;
        arma::mat A_beta = X.t()*weighted_X + arma::diagmat(beta_precision);
        arma::vec h_beta = X.t()*(kappa - omega % (theta.elem(loc)+offset_pg)) +
                           beta_precision % beta_mean;
        beta = glmm_precision_draw(A_beta, h_beta);
        arma::vec w_sum(n,arma::fill::zeros), h_theta(n,arma::fill::zeros);
        arma::vec residual = kappa - omega % (X*beta + offset_pg);
        for (arma::uword j = 0; j < N; ++j) {
            w_sum(loc(j)) += omega(j);
            h_theta(loc(j)) += residual(j);
        }
        theta = glmm_precision_draw(Q/tau2 + arma::diagmat(w_sum), h_theta, center);
        tau2 = 1.0 / R::rgamma(a_tau + 0.5*rank,
                              1.0/(b_tau + 0.5*arma::dot(theta,Q*theta)));
        if (!std::isfinite(tau2) || tau2 <= 0.0)
            Rcpp::stop("Invalid tau2 draw; check the prior and data scale.");
        if (model_indicator == 1) {
            double proposal = R::rnorm(std::log(rho)-std::log1p(-rho), proposal_sd);
            double rho_prop = proposal >= 0 ? 1.0/(1.0+std::exp(-proposal)) :
                               std::exp(proposal)/(1.0+std::exp(proposal));
            double delta = glmm_log_rho(rho_prop,lambda,theta,L,tau2,center,a_rho,b_rho) -
                           glmm_log_rho(rho,lambda,theta,L,tau2,center,a_rho,b_rho);
            if (std::log(R::runif(0.0,1.0)) < delta) {
                rho = rho_prop;
                Q = (1.0-rho)*arma::eye(n,n) + rho*L;
                ++accepted; ++batch_accepted;
            }
            if (adapt_rho && it <= burnin && it % 100 == 0) {
                double rate = batch_accepted/100.0;
                double factor = std::exp(std::min(0.5,1.0/std::sqrt(it/100.0)));
                if (rate > 0.284) proposal_sd *= factor;
                if (rate < 0.184) proposal_sd /= factor;
                proposal_sd = std::max(0.01,std::min(15.0,proposal_sd));
                batch_accepted = 0;
            }
        }
        beta_samples.col(it) = beta; theta_samples.col(it) = theta;
        tau_samples(it) = tau2; rho_samples(it) = rho;
        if (it % 25 == 0) Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(Rcpp::Named("beta")=beta_samples,
                              Rcpp::Named("theta")=theta_samples,
                              Rcpp::Named("tau2")=tau_samples,
                              Rcpp::Named("rho")=rho_samples,
                              Rcpp::Named("accept_rho")=accepted,
                              Rcpp::Named("final_proposal_sd")=proposal_sd);
}
