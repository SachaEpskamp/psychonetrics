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
  bool posdef = X.is_sympd();
  
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
  
  // Check if symmetric:
  bool issym = X.is_symmetric();
  
  // if not, make symmetric:
  if (!issym){
    X = 0.5* (X + X.t());
  }
  
  // Check if posdef:
  bool posdef = X.is_sympd();
  
  // If not, pseudoinverse:
  if (!posdef){
    arma::mat inv = pinv(X);
    res["inv"] = inv;
    
    if (logdet){
      logdetval = log(epsilon);
      res["logdet"] = logdetval;
    }
    
  } else {
    // invert:
    arma::mat inv = inv_sympd(X);
    res["inv"] = inv;
    
    if (logdet){
      logdetval = real(log_det(inv));
      res["logdet"] = logdetval;
    }
  }
  
  return(res);
}

