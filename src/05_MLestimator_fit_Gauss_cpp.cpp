// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// GROUP FIT FUNCTION //
// [[Rcpp::export]]
double maxLikEstimator_Gauss_group_cpp(
    const Rcpp::List& grouplist
){
  arma::mat S = grouplist["S"];
  arma::mat kappa = grouplist["kappa"];
  arma::vec means = grouplist["means"];
  arma::vec mu = grouplist["mu"];
  

  double logdet;
  if (sympd_cpp(kappa) == false){
    // kappa is not positive definite (e.g. a pseudo-inverse of an indefinite
    // implied covariance matrix). The trace and Mahalanobis terms are then
    // unbounded from below, so clamping the log determinant alone (as was
    // done previously) still lets optimizers diverge to -Inf. Return a large
    // penalty instead so the optimizer backs away from the improper region:
    return(1e20);
  } else {
    logdet = log(det(kappa));
    // Plugin for any non-finite log determinant: det(kappa) may be a tiny
    // *negative* number for a near-singular matrix, in which case log()
    // returns NaN (which the old == -Inf check let through):
    if (!std::isfinite(logdet)){
      logdet = log(1.490116e-08);
    }
  }



  arma::mat resvec = trace(S * kappa) + (means - mu).t() * kappa * (means - mu) - logdet;
  double res = resvec(0,0);

  
  return(res);
}


// full fit function 
// [[Rcpp::export]]
double maxLikEstimator_Gauss_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];

  // Result:
  double fit = 0;
  
  for (int i=0; i<nGroup;i++){
    fit += (nPerGroup(i) / nTotal) * maxLikEstimator_Gauss_group_cpp(groupmodels[i]);
  }

  return(fit);
}

