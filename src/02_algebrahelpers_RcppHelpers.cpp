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
arma::vec eig_sym_cpp(
  arma::mat X
){
  return(eig_sym(X));
}

// Check symmetric pd:
// [[Rcpp::export]]
bool sympd_cpp(
  arma::mat X
){
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }

  
  // Check if posdef:
  double epsilon = 1.490116e-08;
  bool posdef = eig_sym(X)[0] > -epsilon; //X.is_sympd();
  
  // return:
  return(posdef);
}

// Symmetric solve:
// [[Rcpp::export]]
Rcpp::List solve_symmetric_cpp(
  arma::mat X,
  bool logdet,
  double epsilon
){
  double logdetval = R_NegInf;
  Rcpp::List res;
  int i;
  int nvar = X.n_cols;
  
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // // Check if posdef:
  // bool posdef = X.is_sympd();
  // Check if posdef:
  double lowestEV = eig_sym(X)[0];
  bool posdef = lowestEV > -epsilon; //X.is_sympd();
  
  // If not, pseudoinverse:
  if (!posdef){
    arma::mat inv = pinv(X);
    res["inv"] = inv;
    
    if (logdet){
      logdetval = log(epsilon);
      res["logdet"] = logdetval;
    }
    
  } else {
    // Smal spectral shift:
    if (lowestEV < 0){
      for (i=0;i<nvar;i++){
        X(i,i) -= lowestEV;
      }
    }
    // Rf_PrintValue(wrap(lowestEV));
    // invert:
    arma::mat inv = inv_sympd(X); // FIXME
    // arma::mat inv = inv(X);
    res["inv"] = inv;
    
    if (logdet){
      logdetval = real(log_det(inv));
      res["logdet"] = logdetval;
    }
  }
  
  return(res);
}

