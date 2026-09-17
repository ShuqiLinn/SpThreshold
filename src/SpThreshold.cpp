#include "RcppArmadillo.h"
#include "SpThreshold.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List SpThreshold(int mcmc_samples,
                       arma::vec y,
                       arma::mat X,
                       arma::uvec loc,
                       int model_indicator,
                       Rcpp::Nullable<Rcpp::NumericMatrix> W = R_NilValue,
                       Rcpp::Nullable<double> proposal_sd_init = R_NilValue,
                       Rcpp::Nullable<double> a_sigma2_prior = R_NilValue,
                       Rcpp::Nullable<double> b_sigma2_prior = R_NilValue,
                       Rcpp::Nullable<double> a_tau2_prior = R_NilValue,
                       Rcpp::Nullable<double> b_tau2_prior = R_NilValue,
                       Rcpp::Nullable<double> a_rho_prior = R_NilValue,
                       Rcpp::Nullable<double> b_rho_prior = R_NilValue,
                       Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                       Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                       Rcpp::Nullable<double> sigma2_init = R_NilValue,
                       Rcpp::Nullable<double> tau2_init = R_NilValue,
                       Rcpp::Nullable<double> rho_init = R_NilValue,
                       int burnin = 0,
                       bool adapt_rho = true,
                       bool verbose = true){

//Data dimensions
if(mcmc_samples < 1) Rcpp::stop("mcmc_samples must be positive.");
if(y.n_elem == 0 || X.n_rows != y.n_elem || loc.n_elem != y.n_elem || X.n_cols == 0)
  Rcpp::stop("Incompatible or empty y, X, and loc.");
if(model_indicator != 0 && model_indicator != 1) Rcpp::stop("model_indicator must be 0 or 1.");
int N = y.n_elem;
int p = X.n_cols;

//Number of spatial units = max(loc) + 1 (0-indexed)
int n = loc.max() + 1;

//Per-location replication counts (for theta_update)
arma::vec m_vec(n); m_vec.fill(0.00);
for(int i = 0; i < N; ++i){
   m_vec(loc(i)) = m_vec(loc(i)) +
                   1.00;
   }

//Adjacency matrix (only needed for spatial model)
arma::mat W_mat(n, n); W_mat.fill(0.00);
if(model_indicator == 1){
  
   if(W.isNull()){
     Rcpp::stop("W (adjacency matrix) required for spatial model (model_indicator = 1).");
     }
   W_mat = Rcpp::as<arma::mat>(W);
  
   }

//Priors
double a_sigma2 = 0.01;
if(a_sigma2_prior.isNotNull()){
   a_sigma2 = Rcpp::as<double>(a_sigma2_prior);
   }

double b_sigma2 = 0.01;
if(b_sigma2_prior.isNotNull()){
   b_sigma2 = Rcpp::as<double>(b_sigma2_prior);
   }

double a_tau2 = 0.01;
if(a_tau2_prior.isNotNull()){
   a_tau2 = Rcpp::as<double>(a_tau2_prior);
   }

double b_tau2 = 0.01;
if(b_tau2_prior.isNotNull()){
   b_tau2 = Rcpp::as<double>(b_tau2_prior);
   }

double a_rho = 1.00;
if(a_rho_prior.isNotNull()){
   a_rho = Rcpp::as<double>(a_rho_prior);
   }

double b_rho = 1.00;
if(b_rho_prior.isNotNull()){
   b_rho = Rcpp::as<double>(b_rho_prior);
   }

//Adaptive tuning
double proposal_sd = 0.30;
if(proposal_sd_init.isNotNull()){
   proposal_sd = Rcpp::as<double>(proposal_sd_init);
   }
int adapt_interval = 100;
double target_accept = 0.234;

//Storage (each iteration is a column)
arma::mat beta_samples(p, mcmc_samples); beta_samples.fill(0.00);
arma::mat theta_samples(n, mcmc_samples); theta_samples.fill(0.00);
arma::vec sigma2_samples(mcmc_samples); sigma2_samples.fill(0.00);
arma::vec tau2_samples(mcmc_samples); tau2_samples.fill(0.00);
arma::vec rho_samples(mcmc_samples); rho_samples.fill(0.00);

//Initial values
arma::vec beta(p); beta.fill(0.00);
if(beta_init.isNotNull()){
   beta = Rcpp::as<arma::vec>(beta_init);
   }

arma::vec theta(n); theta.fill(0.00);
if(theta_init.isNotNull()){
   theta = Rcpp::as<arma::vec>(theta_init);
   }

double sigma2 = 1.00;
if(sigma2_init.isNotNull()){
   sigma2 = Rcpp::as<double>(sigma2_init);
   }

double tau2 = 1.00;
if(tau2_init.isNotNull()){
   tau2 = Rcpp::as<double>(tau2_init);
   }

double rho = 0.50;
if(rho_init.isNotNull()){
   rho = Rcpp::as<double>(rho_init);
   }

// Cache all quantities that do not depend on the MCMC state.  The
// Laplacian eigenvectors are also eigenvectors of the balanced conditional
// precision, so no matrix inverse or Cholesky is required inside that loop.
if(beta.n_elem != (arma::uword)p || theta.n_elem != (arma::uword)n)
  Rcpp::stop("Initial beta or theta has the wrong length.");
if(!(sigma2 > 0 && tau2 > 0 && a_sigma2 > 0 && b_sigma2 > 0 &&
     a_tau2 > 0 && b_tau2 > 0 && a_rho > 0 && b_rho > 0 && proposal_sd > 0))
  Rcpp::stop("Variance, prior, and proposal parameters must be positive.");
if(model_indicator == 1 && !(rho > 0 && rho < 1))
  Rcpp::stop("Initial rho must lie strictly between zero and one.");
arma::mat XtX_inverse;
if(!arma::inv_sympd(XtX_inverse, X.t()*X))
  Rcpp::stop("X must have full column rank under the flat regression prior.");
const arma::mat beta_root = arma::chol(XtX_inverse).t();
const arma::vec Xty = X.t()*y;
arma::vec Zty(n, arma::fill::zeros);
arma::mat ZtX(n, p, arma::fill::zeros);
for(int i = 0; i < N; ++i){
  Zty(loc(i)) += y(i);
  ZtX.row(loc(i)) += X.row(i);
}
const bool balanced = arma::all(m_vec == m_vec(0));
arma::mat laplacian, eigenvectors;
arma::vec eigenvalues;
if(model_indicator == 1){
  if(W_mat.n_rows != (arma::uword)n || W_mat.n_cols != (arma::uword)n ||
     !W_mat.is_finite() || !W_mat.is_symmetric(1e-10) || W_mat.min() < 0)
    Rcpp::stop("W must be a finite symmetric nonnegative n by n matrix.");
  laplacian = arma::diagmat(arma::sum(W_mat, 1)) - W_mat;
  if(!arma::eig_sym(eigenvalues, eigenvectors, laplacian))
    Rcpp::stop("Laplacian eigendecomposition failed.");
  // A nonnegative symmetric adjacency has a positive semidefinite Laplacian.
  // Remove only tiny negative roundoff at its zero eigenvalues.
  eigenvalues.transform([](double x){ return std::max(0.0, x); });
}
arma::vec eigen_Zty, eigen_ones;
arma::mat eigen_ZtX;
if(model_indicator == 1 && balanced){
  eigen_Zty = eigenvectors.t()*Zty;
  eigen_ZtX = eigenvectors.t()*ZtX;
  eigen_ones = eigenvectors.t()*arma::ones(n);
}
const std::string kernel = model_indicator == 0 ? "iid_diagonal" :
  (balanced ? "balanced_spectral" : "unbalanced_precision_cholesky");

//Store first iteration
beta_samples.col(0) = beta;
theta_samples.col(0) = theta;
sigma2_samples(0) = sigma2;
tau2_samples(0) = tau2;
rho_samples(0) = rho;

//Metropolis tracking
int acc_rho_total = 0;
int acc_rho_batch = 0;

//Main loop
for(int iter = 1; iter < mcmc_samples; ++iter){
  
   //1) beta update: cached X'X, X'y, and X'Z; same flat-prior conditional.
   arma::vec z_beta(p);
   for(int k = 0; k < p; ++k) z_beta(k) = R::rnorm(0.0, 1.0);
   beta = XtX_inverse*(Xty - ZtX.t()*theta) + sqrt(sigma2)*beta_root*z_beta;

   //2) theta update.  Each branch draws from the same uncentered Gaussian
   // conditional as before, then applies the existing centering operation.
   const arma::vec rhs = (Zty - ZtX*beta)/sigma2;
   arma::vec z_theta(n);
   for(int k = 0; k < n; ++k) z_theta(k) = R::rnorm(0.0, 1.0);
   arma::vec theta_spectral;
   if(model_indicator == 0){
     const arma::vec precision = m_vec/sigma2 + 1.0/tau2;
     theta = rhs/precision + z_theta/arma::sqrt(precision);
   } else if(balanced){
     const arma::vec precision = m_vec(0)/sigma2 +
       (1.0 - rho + rho*eigenvalues)/tau2;
     theta_spectral = (eigen_Zty - eigen_ZtX*beta)/(sigma2*precision) +
       z_theta/arma::sqrt(precision);
     theta = eigenvectors*theta_spectral;
   } else {
     arma::mat precision = rho*laplacian/tau2;
     precision.diag() += m_vec/sigma2 + (1.0-rho)/tau2;
     // precision = U' U; U^{-1} z has covariance precision^{-1}.
     const arma::mat upper = arma::chol(precision);
     const arma::vec lower_solution = arma::solve(arma::trimatl(upper.t()), rhs);
     theta = arma::solve(arma::trimatu(upper), lower_solution + z_theta);
   }
   const double theta_mean = arma::mean(theta);
   theta -= theta_mean;
   if(model_indicator == 1 && balanced) theta_spectral -= theta_mean*eigen_ones;

   //3) tau2 update: preserve n/2 (no degrees-of-freedom recalibration).
   const double quad_i = arma::dot(theta, theta);
   const double quad_l = model_indicator == 1 ?
     arma::dot(eigenvalues, arma::square(balanced ? theta_spectral : arma::vec(eigenvectors.t()*theta))) : 0.0;
   const double quad_q = model_indicator == 1 ?
     rho*quad_l + (1.0-rho)*quad_i : quad_i;
   tau2 = 1.0/R::rgamma(n/2.0 + a_tau2, 1.0/(0.5*quad_q + b_tau2));

   //4) sigma2 update
   sigma2 = sigma2_update(N, y, X, loc, beta, theta, a_sigma2, b_sigma2);

   //5) rho update: the same logit random walk and Jacobian, using cached
   // eigenvalues for log|Q| and the two quadratic forms for theta'Q theta.
   if(model_indicator == 1){
      const double logit_rho = log(rho/(1.0-rho));
      const double logit_prop = R::rnorm(logit_rho, proposal_sd);
      const double rho_prop = 1.0/(1.0 + exp(-logit_prop));
      double log_acc = R_NegInf;
      if(rho_prop > 0.0 && rho_prop < 1.0){
        const double logdet_old = arma::accu(arma::log(1.0-rho + rho*eigenvalues));
        const double logdet_prop = arma::accu(arma::log(1.0-rho_prop + rho_prop*eigenvalues));
        log_acc = 0.5*(logdet_prop-logdet_old) -
          0.5*(rho_prop-rho)*(quad_l-quad_i)/tau2 +
          a_rho*(log(rho_prop)-log(rho)) +
          b_rho*(log1p(-rho_prop)-log1p(-rho));
      }
      const int accept = log(R::runif(0.0, 1.0)) < log_acc ? 1 : 0;
      if(accept) rho = rho_prop;
      acc_rho_total = acc_rho_total +
                      accept;
      acc_rho_batch = acc_rho_batch +
                      accept;
     
      //Robbins-Monro adaptation during burnin.  The adaptation factor tapers
      //with the iteration count (not the batch count), so early batches take
      //smaller steps as iter grows, matching the Roberts-Rosenthal scheme.
      if(adapt_rho && iter <= burnin && iter % adapt_interval == 0){
        
         double batch_rate = (double)acc_rho_batch/(double)adapt_interval;
         double adapt_factor = exp(std::min(0.50, 1.00/sqrt((double)iter/(double)adapt_interval)));
        
         if(batch_rate > target_accept + 0.05){
           proposal_sd = proposal_sd*adapt_factor;
           }
         if(batch_rate < target_accept - 0.05){
           proposal_sd = proposal_sd/adapt_factor;
           }
        
         //Clamp
         if(proposal_sd < 0.01){
           proposal_sd = 0.01;
           }
         if(proposal_sd > 15.00){
           proposal_sd = 15.00;
           }
        
         acc_rho_batch = 0;
        
         }
     
      }
  
   //Save samples
   beta_samples.col(iter) = beta;
   theta_samples.col(iter) = theta;
   sigma2_samples(iter) = sigma2;
   tau2_samples(iter) = tau2;
   rho_samples(iter) = rho;
  
   //Progress and user interrupt
   if((iter + 1) % 10 == 0){
     Rcpp::checkUserInterrupt();
     }
  
   if(verbose && (iter + 1) % std::max(1, int(round(mcmc_samples*0.10))) == 0){
     
      double completion = round(100.00*(iter + 1)/(double)mcmc_samples);
      Rcpp::Rcout << "Progress: " << completion << "%";
     
      if(model_indicator == 1){
         double accrate = round(100.00*(double)acc_rho_total/(double)iter);
         Rcpp::Rcout << " | rho acceptance: " << accrate << "%"
                     << " | proposal_sd: " << proposal_sd;
         }
     
      Rcpp::Rcout << std::endl;
     
      }
  
   }

//Build return list
if(model_indicator == 1){
  
   return Rcpp::List::create(Rcpp::Named("beta")              = beta_samples,
                             Rcpp::Named("theta")             = theta_samples,
                             Rcpp::Named("sigma2")            = sigma2_samples,
                             Rcpp::Named("tau2")              = tau2_samples,
                             Rcpp::Named("rho")               = rho_samples,
                             Rcpp::Named("accept_rho")        = acc_rho_total,
                             Rcpp::Named("final_proposal_sd") = proposal_sd,
                             Rcpp::Named("kernel") = kernel);
  
   }

return Rcpp::List::create(Rcpp::Named("beta")   = beta_samples,
                          Rcpp::Named("theta")  = theta_samples,
                          Rcpp::Named("sigma2") = sigma2_samples,
                          Rcpp::Named("tau2")   = tau2_samples,
                          Rcpp::Named("kernel") = kernel);

}
