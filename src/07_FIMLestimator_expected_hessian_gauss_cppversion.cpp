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
arma::mat expected_hessian_fiml_Gaussian_group_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // double result = 0;
  // double logdet = 0;
  double n_part;
  // double log2pi = log(2*M_PI);
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Hes = zeros(nvars + nmeans, nvars + nmeans);
  
  // Integers
  int i;
  int j;
  // int ii;
  // int jj;
  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){
    
    // Obtain list:
    Rcpp::List dat = fimldata[i];
    
    // Observed indices:
    vec obs = dat["obs"];
    
    uvec inds = find(obs == true);
    
    // sample size:
    n_part = dat["n"];
    
    // Duplication matrix:
    arma::sp_mat D = dat["D"];
    arma::sp_mat L = dat["L"];
    
    // Subset matrices:
    arma::mat sigma_p = sigma(inds,inds);
    arma::mat kappa_p(sigma_p);
    arma::vec mu_p = mu(inds);
    
    int nvar = mu_p.size();
    
    // Observed values:
    arma::mat S = dat["S"];
    arma::vec means = dat["means"];
    
    // inverse:
    arma::vec ev = arma::eig_sym(sigma_p);
    bool ispos = true;
    for ( j = 0; j < ev.size(); j++){
      if (ev[j] < sqrt(epsilon)){
        ispos = false;
        break;
      }
    }
    if (ispos){
      kappa_p = inv(sigma_p);
      // logdet = log(det(kappa_p));
    } else {
      kappa_p = pinv(sigma_p);
      // logdet = log(epsilon);
    }
    
    // Mean part:
    arma::mat meanpart = 2.0 * kappa_p;
    // Sigma part:
    arma::mat sigmapart = D.t() * kron(kappa_p, kappa_p) * D;
    
    // Zeroes:
    arma::mat O = zeros(nvar, nvar * (nvar + 1) / 2);
    
    // 
    // Rf_PrintValue(wrap(meanpart));
    // Rf_PrintValue(wrap(sigmapart));
    
    // Join:
    Hes += n_part * L.t() * join_cols(
      join_rows(meanpart, O),
      join_rows(O.t(), sigmapart)
    ) * L;
    // nvar <- ncol(kappa)
    // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  }
  
  
  // Return
  return Hes;
}
