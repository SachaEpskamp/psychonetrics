// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"

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

// Inner function:
// [[Rcpp::export]]
arma::mat expected_hessian_fiml_Gaussian_group_cpp_inner(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    Rcpp::List dat,
    double epsilon) {
  

  // double result = 0;
  // double logdet = 0;
  double n_part;
  // double log2pi = log(2*M_PI);
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  // arma::mat Hes = zeros(nvars + nmeans, nvars + nmeans);
  
  // Integers
  // int i;
  int j;
  // int ii;
  // int jj;
  
  // Loop over groups
  // for ( i = 0; i < fimldata.size(); i++){
  
  // Obtain list:
  // Rcpp::List dat = fimldata[i];
  
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
  // arma::vec ev = arma::eig_sym(sigma_p);
  // bool ispos = true;
  // for ( j = 0; j < ev.size(); j++){
  //   if (ev[j] < sqrt(epsilon)){
  //     ispos = false;
  //     break;
  //   }
  // }
  // bool ispos = sigma_p.is_sympd();
  // 
  // if (ispos){
  //   kappa_p = inv_sympd(sigma_p);
  //   // logdet = log(det(kappa_p));
  // } else {
  //   kappa_p = pinv(sigma_p);
  //   // logdet = log(epsilon);
  // }
  kappa_p = solve_symmetric_cpp_matrixonly(sigma_p, epsilon);
  
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
  arma::mat Hes = n_part * L.t() * join_cols(
    join_rows(meanpart, O),
    join_rows(O.t(), sigmapart)
  ) * L;
  // nvar <- ncol(kappa)
  // res <-  attr(kappa, "logdet") - nvar * log((2*pi)) - sum(diag(SK)) - t(means - mu) %*% kappa %*% (means - mu)
  // }
  
  
  // Return
  return Hes;
}

// Outer function:
// [[Rcpp::export]]
arma::mat expected_hessian_fiml_Gaussian_group_cpp(
    const arma::mat& sigma, 
    const arma::mat& kappa,
    const arma::vec& mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // double result = 0;
  // double logdet = 0;
  // double n_part;
  // double log2pi = log(2*M_PI);
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Hes = zeros(nvars + nmeans, nvars + nmeans);
  
  // Integers
  int i;
  // int j;
  // int ii;
  // int jj;
  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){

    // Join:
    Hes += expected_hessian_fiml_Gaussian_group_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
      epsilon);
    
  }
  
  
  // Return
  return Hes;
}

// Outer function:
// [[Rcpp::export]]
arma::mat expected_hessian_fiml_Gaussian_group_cpp_fullFIML(
    Rcpp::List sigma, 
    Rcpp::List kappa,
    Rcpp::List mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // double result = 0;
  // double logdet = 0;
  // double n_part;
  // double log2pi = log(2*M_PI);
  
  // Number of parameters
  arma::vec firstmu = mu[0];
  int nmeans = firstmu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Hes = zeros(nvars + nmeans, nvars + nmeans);
  
  // Integers
  int i;
  // int j;
  // int ii;
  // int jj;
  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){
    
    // Join:
    Hes += expected_hessian_fiml_Gaussian_group_cpp_inner(
      sigma[i], 
      kappa[i],
      mu[i],
      fimldata[i],
              epsilon);
    
  }
  
  
  // Return
  return Hes;
}

