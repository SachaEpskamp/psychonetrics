// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
double fimlEstimator_Gauss_group_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon,
    double n) {
  
  double result = 0;
  double logdet = 0;
  double n_part;
  
  // Loop over groups
  for (int i = 0; i < fimldata.size(); i++){
  
    // Obtain list:
    Rcpp::List dat = fimldata[i];
    
    // Observed indices:
    vec obs = dat["obs"];
    uvec inds = find(obs == true);
  
  // sample size:
  n_part = dat["n"];
    
    // Subset matrices:
    arma::mat sigma_p = sigma(inds,inds);
    arma::mat kappa_p = inv(sigma_p);
    arma::vec mu_p = mu(inds);
    
    // Observed values:
    arma::mat S = dat["S"];
    arma::vec means = dat["means"];
    
    // log det:
    arma::vec ev = arma::eig_sym(kappa_p);
    bool ispos = true;
    for (int j = 0; j < ev.size(); j++){
      if (ev[j] < sqrt(epsilon)){
        ispos = false;
        break;
      }
    }
    if (ispos){
      logdet = log(det(kappa_p));
    } else {
      logdet = log(epsilon);
    }
    
    // Likelihood:
    result += n_part * (trace(S * kappa_p) + 
      dot((means - mu_p).t(), kappa_p * (means - mu_p)) - 
      logdet);
    
  }
  
  
  // Return
  return (1/n) * result;
}
