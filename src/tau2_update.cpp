#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double tau2_update(int n,
                   arma::vec theta,
                   arma::mat Q,
                   double a_tau2,
                   double b_tau2){

double quad = dot(theta, (Q*theta));

double a_tau2_update = n/2.00 +
                       a_tau2;
double b_tau2_update = 0.50*quad +
                       b_tau2;

double tau2 = 1.00/R::rgamma(a_tau2_update,
                             (1.00/b_tau2_update));

return(tau2);

}
