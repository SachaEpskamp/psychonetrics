// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

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


// [[Rcpp::export]]
double logLikelihood_gaussian_subgroup_fiml_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  double result = 0;
  double logdet = 0;
  double n_part;
  double log2pi = log(2*M_PI);
  
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
    arma::mat kappa_p(sigma_p);
    arma::vec mu_p = mu(inds);
    
    int nvar = mu_p.size();
    
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
    result += n_part * (
      logdet - nvar * log2pi - trace(S * kappa_p) -
        dot((means - mu_p).t(), kappa_p * (means - mu_p))
    );
    
    // nvar <- ncol(kappa)
      // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  }
  
  
  // Return
  return result;
}
