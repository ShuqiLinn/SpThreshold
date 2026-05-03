################
#Generating data
################

#W: Spatial adjacency matrix
create_W_fun<-function(n_locs){

              W<-matrix(0, 
                        nrow = n_locs, 
                        ncol = n_locs)
              perm<-sample(n_locs)
              for(i in 2:n_locs){
 
                 a<-perm[i]
                 b<-perm[sample.int(i - 1, 1)]
                 W[a,b]<-1
                 W[b,a]<-1
 
                 }

              return(W)

              }

#set.seed(849)
n<-50
W<-create_W_fun(n_locs = n)

sigma2<-0.90
tau2<-1.00 - sigma2
rho<-0.80
gamma<-0.05

#Generate x averages within each location
m_big<-1000000
x<-rep(NA,
       times = (n*m_big))
for(j in 1:n){ 
  
   loc_mean<-rnorm(n = 1)
   x[(1 + (j-1)*m_big):(j*m_big)]<-rnorm(n = m_big,
                                         mean = loc_mean,
                                         sd = 0.50)

   }
x<-(x - mean(x))/sqrt(mean((x - mean(x))^2))
x_bar_vec<-sapply(c(1:n), 
                  function(j) mean(x[(1 + (j-1)*m_big):(j*m_big)]))

#Eigendecomposition: Reorder to get \lambda_1 = 0 first
e<-eigen(diag(rowSums(W)) - W)
ord<-order(e$values) 
e_val<-e$values[ord]
e_vec<-e$vectors[,ord]

#############################################
#Validating the spectral decomposition: Q_inv
#############################################
Q_inv_true<-chol2inv(chol(rho*(diag(rowSums(W)) - W) + (1.00 - rho)*diag(n)))
Q_inv_spectral<-tcrossprod(e_vec[,1])/(rho*e_val[1] + 1.00 - rho)
for(j in 2:n){
   Q_inv_spectral<-Q_inv_spectral +
                   tcrossprod(e_vec[,j])/(rho*e_val[j] + 1.00 - rho)
   }
mean(round(Q_inv_true, 10) == round(Q_inv_spectral, 10))

#################################################################
#Validating the spectral decomposition: Full conditional variance
#################################################################
m_test<-20
x_test<-rnorm(n = (n*m_test))
x_test<-(x_test - mean(x_test))/sqrt(mean((x_test - mean(x_test))^2))
x_bar_vec_test<-sapply(c(1:n), 
                       function(j) mean(x_test[(1 + (j-1)*m_test):(j*m_test)]))

d_num_test<-
d_den_test<-rep(NA, 
                times = n)
for(i in 1:n){

   d_den_test[i]<-sum(x_bar_vec_test*e_vec[,i])^2 
   d_num_test[i]<-d_den_test[i]*(1.00 - e_val[i])

   }

Z<-diag(n)%x%rep(1, times = m_test)
Omega_true<-tau2*tcrossprod((Z%*%Q_inv_true), Z) +
            sigma2*diag(n*m_test)
beta1_var_true<-1.00/crossprod(x_test, (chol2inv(chol(Omega_true))%*%x_test))
beta1_var_spectral<-1.00/(n*m_test/sigma2 - ((m_test^2)*tau2/sigma2)*sum(d_den_test/(sigma2*rho*e_val + sigma2 - sigma2*rho + m_test*tau2)))
as.numeric(round(beta1_var_true, 10) == round(beta1_var_spectral, 10))

#########################################################################################
#Validating the spectral decomposition: Relative difference in full conditional variances
#########################################################################################
Q_zero_inv_true<-diag(n)
Omega_zero_true<-tau2*tcrossprod((Z%*%Q_zero_inv_true), Z) +
                 sigma2*diag(n*m_test)
beta1_var_zero_true<-1.00/crossprod(x_test, (chol2inv(chol(Omega_zero_true))%*%x_test))
rel_diff_true<-(beta1_var_true - beta1_var_zero_true)/beta1_var_zero_true
rel_diff_spectral<-((m_test^2)*tau2*sigma2*rho*sum(d_num_test/((sigma2 + m_test*tau2)*(sigma2*rho*e_val + sigma2 - sigma2*rho + m_test*tau2))))/(n*m_test - (m_test^2)*tau2*sum(d_den_test/(sigma2*rho*e_val + sigma2 - sigma2*rho + m_test*tau2)))
as.numeric(round(rel_diff_true, 10) == round(rel_diff_spectral, 10))

##################
#Approximate bound
##################
d_num<-
d_den<-rep(NA, 
           times = n)
for(i in 1:n){

   d_den[i]<-sum(x_bar_vec*e_vec[,i])^2  
   d_num[i]<-d_den[i]*(1.00 - e_val[i])

   }

m_approx<-max(2, 
              ceiling(abs((rho*sigma2/tau2)*sum(d_num))/(gamma*(n - sum(d_den)))))
m_approx

############
#Exact bound
############
m_exact_fun<-function(m_test,
                      x_bar_vec){

                      d_num<-
                      d_den<-rep(NA, 
                                 times = n)
                      for(i in 1:n){

                         d_den[i]<-(sum(x_bar_vec*e_vec[,i])^2)/(sigma2*rho*e_val[i] + sigma2 - sigma2*rho + m_test*tau2)  
                         d_num[i]<-d_den[i]*(1.00 - e_val[i])/(sigma2 + m_test*tau2)

                         }

                      return(abs((m_test^2)*tau2*sigma2*rho*sum(d_num)/(n*m_test - (m_test^2)*tau2*sum(d_den))))

                      }

abs_diff<-rep(NA,
              times = max(100, ceiling(m_approx)))
for(j in 1:max(100, ceiling(m_approx))){
   abs_diff[j]<-m_exact_fun(m_test = j,
                            x_bar_vec = x_bar_vec)
   }
plot(abs_diff,
     type = "l")
abline(h = gamma,
       col = "red")
abline(v = m_approx,
       col = "blue")
