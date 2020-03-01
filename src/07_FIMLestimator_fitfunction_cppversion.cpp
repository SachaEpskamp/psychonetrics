// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// mat test(
//   mat X,
//   mat Y,
//   uvec ind
// ){
//   X(ind,ind) += Y;
//   return X;
// }



// Inner function
// [[Rcpp::export]]
double fimlEstimator_Gauss_group_cpp_inner(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List dat,
    double epsilon,
    double n) {
  
  double logdet = 0;
  double n_part;
  
  // Observed indices:
  vec obs = dat["obs"];
  uvec inds = find(obs == true);
  
  // sample size:
  n_part = dat["n"];
  
  // Subset matrices:
  arma::mat sigma_p = sigma(inds,inds);
  arma::mat kappa_p(sigma_p);
  arma::vec mu_p = mu(inds);
  
  // Observed values:
  arma::mat S = dat["S"];
  arma::vec means = dat["means"];
  
  // log det:
  arma::vec ev = arma::eig_sym(sigma_p);
  bool ispos = true;
  for (int j = 0; j < ev.size(); j++){
    if (ev[j] < sqrt(epsilon)){
      ispos = false;
      break;
    }
  }
  if (ispos){
    kappa_p = inv(sigma_p);
    logdet = log(det(kappa_p));
  } else {
    kappa_p = pinv(sigma_p);
    logdet = log(epsilon);
  }
  
  // Likelihood:
  double result = n_part * (trace(S * kappa_p) + 
                            dot((means - mu_p).t(), kappa_p * (means - mu_p)) - 
                            logdet);
  
  // Return
  return result;
}


// Outer normal function:
// [[Rcpp::export]] 
double fimlEstimator_Gauss_group_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon,
    double n) {
  // Rf_PrintValue(wrap("USED"));
  
  double result = 0;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
    
    result += fimlEstimator_Gauss_group_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
      epsilon,
      n);
  }
  
  
  // Return
  return (1/n) * result;
}



// Outer function PER GROUP:
// [[Rcpp::export]] 
double fimlEstimator_Gauss_group_cpp_perGroup(
    Rcpp::List sigma, 
    Rcpp::List kappa,
    Rcpp::List mu,
    Rcpp::List fimldata,
    double epsilon,
    double n) {
  // Rf_PrintValue(wrap("USED"));
  
  double result = 0;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
    
    result += fimlEstimator_Gauss_group_cpp_inner(
      sigma[i], 
      kappa[i],
      mu[i],
      fimldata[i],
              epsilon,
              n);
  }
  
  
  // Return
  return (1/n) * result;
}
