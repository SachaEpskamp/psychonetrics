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

// I kron X
// [[Rcpp::export]]
arma::sp_mat kronecker_I_X(
  const arma::mat& X,
  int n
){
  int i,j,k;
  
  // Size of X:
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  
  // return matrix:
  arma::sp_mat res = sp_mat(n*nrow,n*ncol);
  
  // Fill the matrices:
  for (i=0; i<nrow; i++){
    for (j=0;j<ncol; j++){
      for (k=0;k<n;k++){
        res(i + k*nrow, j + k*ncol) = X(i,j);
        
        
      }
    }
  }
  
  return(res);
}

// X kron I
// [[Rcpp::export]]
arma::sp_mat kronecker_X_I(
    const arma::mat& X,
    int n
){
  int i,j,k;
  
  // Size of X:
  int nrow = X.n_rows;
  int ncol = X.n_cols;
  
  // return matrix:
  arma::sp_mat res = sp_mat(n*nrow,n*ncol);
  
  // Fill the matrices:
  for (i=0; i<nrow; i++){
    for (j=0;j<ncol; j++){
      for (k=0;k<n;k++){
        
        
        res(k + i*n, k + j*n) = X(i,j);
        
        
      }
    }
  }
  
  return(res);
}

// Kron of diagonal matrix with itself:
// [[Rcpp::export]]
arma::sp_mat kronecker_diag_sparse(
  arma::sp_mat X
){
  int i,j;
  
  // Size of X:
  int n = X.n_rows;

  // return matrix:
  arma::sp_mat res = sp_mat(n*n,n*n);
  
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      res(i*n + j, i*n + j) = X(i,i) * X(j,j);
    }
  }
  
  return(res);
}

// [[Rcpp::export]]
arma::sp_mat kronecker_diag(
    arma::mat X
){
  int i,j;
  
  // Size of X:
  int n = X.n_rows;
  
  // return matrix:
  arma::sp_mat res = sp_mat(n*n,n*n);
  
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      res(i*n + j, i*n + j) = X(i,i) * X(j,j);
    }
  }
  
  return(res);
}


