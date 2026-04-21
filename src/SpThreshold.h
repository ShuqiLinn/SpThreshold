#ifndef __SpThreshold__
#define __SpThreshold__

double sigma2_update(int N,
                     arma::vec y,
                     arma::mat X,
                     arma::uvec loc,
                     arma::vec beta_old,
                     arma::vec theta_old,
                     double a_sigma2,
                     double b_sigma2);

double tau2_update(int n,
                   arma::vec theta,
                   arma::mat Q,
                   double a_tau2,
                   double b_tau2);

arma::vec beta_update(int N,
                      int p,
                      arma::vec y,
                      arma::mat X,
                      arma::uvec loc,
                      arma::vec theta_old,
                      double sigma2_old);

arma::vec theta_update(int N,
                       int n,
                       arma::vec y,
                       arma::mat X,
                       arma::uvec loc,
                       arma::vec m_vec,
                       arma::vec beta_old,
                       double sigma2_old,
                       double tau2_old,
                       arma::mat Q);

Rcpp::List rho_update(int n,
                      arma::mat W,
                      arma::vec theta,
                      double tau2,
                      double rho_old,
                      arma::mat Q_old,
                      double Q_log_det_old,
                      double proposal_sd,
                      double a_rho,
                      double b_rho);

Rcpp::List SpThreshold(int mcmc_samples,
                       arma::vec y,
                       arma::mat X,
                       arma::uvec loc,
                       int model_indicator,
                       Rcpp::Nullable<Rcpp::NumericMatrix> W,
                       Rcpp::Nullable<double> proposal_sd_init,
                       Rcpp::Nullable<double> a_sigma2_prior,
                       Rcpp::Nullable<double> b_sigma2_prior,
                       Rcpp::Nullable<double> a_tau2_prior,
                       Rcpp::Nullable<double> b_tau2_prior,
                       Rcpp::Nullable<double> a_rho_prior,
                       Rcpp::Nullable<double> b_rho_prior,
                       Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                       Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                       Rcpp::Nullable<double> sigma2_init,
                       Rcpp::Nullable<double> tau2_init,
                       Rcpp::Nullable<double> rho_init,
                       int burnin,
                       bool adapt_rho);

#endif // __SpThreshold__
