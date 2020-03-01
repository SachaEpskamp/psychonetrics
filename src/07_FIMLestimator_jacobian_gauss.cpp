#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Inner function
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List dat,
    double epsilon) {
  // Rf_PrintValue(wrap("USED"));
  
  double n_part;
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  // arma::mat Jac = zeros(1, nvars + nmeans);
  
  // Integers
  // int i;
  int j;
  
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
  } else {
    kappa_p = pinv(sigma_p);
  }
  
  // Mean part:
  arma::mat meanpart = -2 * (means - mu_p).t() * kappa_p;
  // Sigma part:
  arma::mat matpart = kappa_p * (S + (means - mu_p) * (means - mu_p).t() - sigma_p) * kappa_p;
  
  //Vectorize:
  arma::vec vecpart = vectorise(matpart);
  arma::mat sigmapart = -1.0 * (D.t() * vecpart).t();
  
  // Join:
  arma::mat Jac =  n_part * join_rows(meanpart, sigmapart) * L;
  
  // }
  
  // Return
  return Jac;
}


// Outer function
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp(
    arma::mat sigma, 
    arma::mat kappa,
    arma::vec mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Jac = zeros(1, nvars + nmeans);
  
  // Integers
  int i;

  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){

    // Join:
    Jac += jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
      sigma, 
      kappa,
      mu,
      fimldata[i],
       epsilon);
    
  }
  
  // Return
  return Jac;
}


// Outer function
// [[Rcpp::export]]
arma::mat jacobian_fiml_gaussian_subgroup_sigma_cpp_perGroup(
    Rcpp::List sigma, 
    Rcpp::List kappa,
    Rcpp::List mu,
    Rcpp::List fimldata,
    double epsilon) {
  
  // Number of parameters
  int nmeans = mu.size();
  int nvars = nmeans * (nmeans + 1) / 2;
  
  // Empty Jacobian:
  arma::mat Jac = zeros(1, nvars + nmeans);
  
  // Integers
  int i;
  
  
  // Loop over groups
  for ( i = 0; i < fimldata.size(); i++){
    
    // Join:
    Jac += jacobian_fiml_gaussian_subgroup_sigma_cpp_inner(
      sigma[i], 
      kappa[i],
      mu[i],
      fimldata[i],
              epsilon);
    
  }
  
  // Return
  return Jac;
}



