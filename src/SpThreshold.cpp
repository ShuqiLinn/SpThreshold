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
                       bool adapt_rho = true){

//Data dimensions
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

//Initialize Q and log|Q|
//Spatial:     Q = rho*(D - W) + (1 - rho)*I
//Nonspatial:  Q = I  (effectively rho = 0)
arma::mat Q(n, n); Q.fill(0.00);
double Q_log_det = 0.00;
double sign_det = 0.00;

if(model_indicator == 1){
  
   arma::mat D = arma::diagmat(arma::sum(W_mat, 1));
   Q = rho*(D - W_mat) +
       (1.00 - rho)*arma::eye(n, n);
   arma::log_det(Q_log_det, sign_det, Q);
  
   }
if(model_indicator == 0){
  
   Q = arma::eye(n, n);
   Q_log_det = 0.00;
  
   }

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
  
   //1) beta update
   beta = beta_update(N, p, y, X, loc, theta, sigma2);
  
   //2) theta update
   theta = theta_update(N, n, y, X, loc, m_vec, beta, sigma2, tau2, Q);
  
   //3) tau2 update
   tau2 = tau2_update(n, theta, Q, a_tau2, b_tau2);
  
   //4) sigma2 update
   sigma2 = sigma2_update(N, y, X, loc, beta, theta, a_sigma2, b_sigma2);
  
   //5) rho update (spatial only)
   if(model_indicator == 1){
     
      Rcpp::List rho_out = rho_update(n, W_mat, theta, tau2, rho, Q,
                                       Q_log_det, proposal_sd, a_rho, b_rho);
      rho       = Rcpp::as<double>(rho_out["rho"]);
      Q         = Rcpp::as<arma::mat>(rho_out["Q"]);
      Q_log_det = Rcpp::as<double>(rho_out["Q_log_det"]);
      int accept = Rcpp::as<int>(rho_out["accept"]);
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
  
   if((iter + 1) % int(round(mcmc_samples*0.10)) == 0){
     
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
                             Rcpp::Named("final_proposal_sd") = proposal_sd);
  
   }

return Rcpp::List::create(Rcpp::Named("beta")   = beta_samples,
                          Rcpp::Named("theta")  = theta_samples,
                          Rcpp::Named("sigma2") = sigma2_samples,
                          Rcpp::Named("tau2")   = tau2_samples);

}
