#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(int N,
                      int p,
                      arma::vec y,
                      arma::mat X,
                      arma::uvec loc,
                      arma::vec theta_old,
                      double sigma2_old){

arma::vec theta_expanded(N); theta_expanded.fill(0.00);
for(int i = 0; i < N; ++i){
   theta_expanded(i) = theta_old(loc(i));
   }

arma::vec resid = y +
                  -theta_expanded;

arma::mat cov_beta = inv_sympd(trans(X)*X/sigma2_old);
arma::vec mean_beta = cov_beta*(trans(X)*resid/sigma2_old);

arma::vec z(p); z.fill(0.00);
for(int k = 0; k < p; ++k){
   z(k) = R::rnorm(0.00, 1.00);
   }

arma::vec beta = mean_beta +
                 trans(arma::chol(cov_beta))*z;

return(beta);

}
