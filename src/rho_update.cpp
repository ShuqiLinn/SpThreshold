#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List rho_update(int n,
                      arma::mat W,
                      arma::vec theta,
                      double tau2,
                      double rho_old,
                      arma::mat Q_old,
                      double Q_log_det_old,
                      double proposal_sd,
                      double a_rho,
                      double b_rho){

//Current state:  log-posterior evaluated at rho_old.
//log f(rho | theta, tau2) = 0.5 * log|Q(rho)|
//                           - 0.5/tau2 * theta^T Q(rho) theta
//                           + (a_rho - 1)*log(rho) + (b_rho - 1)*log(1 - rho)
double quad_old = dot(theta, (Q_old*theta));
double lp_curr = 0.50*Q_log_det_old +
                 -0.50*quad_old/tau2 +
                 (a_rho - 1.00)*log(rho_old) +
                 (b_rho - 1.00)*log(1.00 - rho_old);

//Propose on logit scale:  logit_prop ~ Normal(logit(rho_old), proposal_sd^2)
double logit_rho = log(rho_old/(1.00 - rho_old));
double logit_prop = R::rnorm(logit_rho, proposal_sd);
double rho_prop = 1.00/(1.00 + exp(-logit_prop));

//Build Q(rho_prop) = rho_prop*(D - W) + (1 - rho_prop)*I
arma::mat D = arma::diagmat(arma::sum(W, 1));
arma::mat Q_prop = rho_prop*(D - W) +
                   (1.00 - rho_prop)*arma::eye(n, n);

//Log-determinant of Q(rho_prop)
double Q_log_det_prop = 0.00;
double sign = 0.00;
arma::log_det(Q_log_det_prop, sign, Q_prop);

//Proposed-state log-posterior
double quad_prop = dot(theta, (Q_prop*theta));
double lp_prop = 0.50*Q_log_det_prop +
                 -0.50*quad_prop/tau2 +
                 (a_rho - 1.00)*log(rho_prop) +
                 (b_rho - 1.00)*log(1.00 - rho_prop);

//Jacobian of logit transformation:  |drho/dlogit| = rho*(1 - rho)
//log-Jacobian = log(rho) + log(1 - rho)
double log_jac_curr = log(rho_old) +
                      log(1.00 - rho_old);
double log_jac_prop = log(rho_prop) +
                      log(1.00 - rho_prop);

double log_acc = lp_prop +
                 -lp_curr +
                 log_jac_prop +
                 -log_jac_curr;

//Accept/reject
int accept = 0;
double rho = rho_old;
arma::mat Q = Q_old;
double Q_log_det = Q_log_det_old;

if(log(R::runif(0.00, 1.00)) < log_acc){
  
   rho = rho_prop;
   Q = Q_prop;
   Q_log_det = Q_log_det_prop;
   accept = 1;
  
   }

return Rcpp::List::create(Rcpp::Named("rho")       = rho,
                          Rcpp::Named("Q")         = Q,
                          Rcpp::Named("Q_log_det") = Q_log_det,
                          Rcpp::Named("accept")    = accept);

}
