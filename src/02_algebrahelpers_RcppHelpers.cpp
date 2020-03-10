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
    // Small spectral shift:
    if (lowestEV < epsilon){
      for (i=0;i<nvar;i++){
        //X(i,i) -= lowestEV;
        X(i,i) = X(i,i) -lowestEV + epsilon;
      }
    }
    // Rf_PrintValue(wrap(lowestEV));
    // invert:
    arma::mat inv = inv_sympd(X); // FIXME
    // arma::mat inv = inv(X);
    res["inv"] = inv;
    
    if (logdet){
      double logepsilon = log(epsilon);
      logdetval =  log(det(inv));
      if (logdetval < logepsilon){
        logdetval = logepsilon;
      }
      
      res["logdet"] = logdetval;
    }
  }
  
  return(res);
}


// [[Rcpp::export]]
arma::mat solve_symmetric_cpp_matrixonly(
    arma::mat X,
    double epsilon
){
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
    arma::mat res = pinv(X);
    return(res);
    
  } else {
    // Small spectral shift:
    if (lowestEV < epsilon){
      for (i=0;i<nvar;i++){
        X(i,i) = X(i,i) -lowestEV + epsilon;
      }
    }

    // invert:
    arma::mat res = inv_sympd(X); // FIXME
    return(res);
  }
  

}
