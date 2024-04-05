// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Eigenvalues of symmetric matrix:
// [[Rcpp::export]]
arma::sp_mat Mmatrix_cpp(
    Rcpp::DataFrame parDF
){
  /// Obtain pars:
  arma::vec par = parDF["par"];
  int allPar = par.n_rows;
  int maxPar = max(par);
  
  // Create matrix:
  sp_mat M(allPar, maxPar);
  
  // Loop over the data:
  int i;
  for (i=0;i<allPar;i++){
    if (par(i) > 0){
      M(i, par(i)-1) = 1;
    }
    
  }
  
  return(M);
}


// [[Rcpp::export]]
arma::sp_mat Mmatrix_cpp_list(
    Rcpp::List parDF
){
  /// Obtain pars:
  arma::vec par = parDF["par"];
  int allPar = par.n_rows;
  int maxPar = max(par);
  
  // Create matrix:
  sp_mat M(allPar, maxPar);
  
  // Loop over the data:
  int i;
  for (i=0;i<allPar;i++){
    if (par(i) > 0){
      M(i, par(i)-1) = 1;
    }
    
  }
  
  return(M);
}
