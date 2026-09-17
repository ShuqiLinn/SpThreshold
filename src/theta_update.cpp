#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec theta_update(int N,
                       int n,
                       arma::vec y,
                       arma::mat X,
                       arma::uvec loc,
                       arma::vec m_vec,
                       arma::vec beta_old,
                       double sigma2_old,
                       double tau2_old,
                       arma::mat Q){

//Residual after fixed effects:  y - X*beta
arma::vec resid(N); resid.fill(0.00);
resid = y +
        -X*beta_old;

//Z^T (y - X*beta):  per-location sum of residuals
arma::vec Zt_resid(n); Zt_resid.fill(0.00);
for(int i = 0; i < N; ++i){
   Zt_resid(loc(i)) = Zt_resid(loc(i)) +
                      resid(i);
   }

//Precision of the full conditional:  diag(m_vec)/sigma2 + Q/tau2
arma::mat Sigma_theta_inv = arma::diagmat(m_vec)/sigma2_old +
                            Q/tau2_old;

// Precision = U' U.  Solve triangular systems without explicitly forming
// the inverse covariance; dense Cholesky does not use a permutation.
arma::mat upper = arma::chol(Sigma_theta_inv);
arma::vec lower_solution = arma::solve(arma::trimatl(upper.t()), Zt_resid/sigma2_old);
arma::vec mean_theta = arma::solve(arma::trimatu(upper), lower_solution);

//Draw z ~ N(0, I) using R's RNG stream
arma::vec z(n); z.fill(0.00);
for(int k = 0; k < n; ++k){
   z(k) = R::rnorm(0.00, 1.00);
   }

//Cov(U^{-1} z) = (U' U)^{-1}.
arma::vec theta = mean_theta + arma::solve(arma::trimatu(upper), z);

//Sum-to-zero constraint
theta = theta +
        -arma::mean(theta);

return(theta);

}
