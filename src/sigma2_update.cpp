#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_update(int N,
                     arma::vec y,
                     arma::mat X,
                     arma::uvec loc,
                     arma::vec beta_old,
                     arma::vec theta_old,
                     double a_sigma2,
                     double b_sigma2){

arma::vec theta_expanded(N); theta_expanded.fill(0.00);
for(int i = 0; i < N; ++i){
   theta_expanded(i) = theta_old(loc(i));
   }

arma::vec resid = y +
                  -X*beta_old +
                  -theta_expanded;

double a_sigma2_update = N/2.00 +
                         a_sigma2;
double b_sigma2_update = 0.50*dot(resid, resid) +
                         b_sigma2;

double sigma2 = 1.00/R::rgamma(a_sigma2_update,
                               (1.00/b_sigma2_update));

return(sigma2);

}
