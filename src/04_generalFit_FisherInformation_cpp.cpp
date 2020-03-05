// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Dense - sparse - sparse
// [[Rcpp::export]]
arma::mat FisherInformation_inner_cpp_DSS(
    const arma::mat& estimator,
    const arma::sp_mat& model,
    const arma::sp_mat& manual
) {
  // Sparse part first:
  arma::sp_mat sparse = model * manual;
  
  // Now second part:
  arma::mat Fis = 0.5 * sparse.t() * estimator * sparse;
  
  // Return
  return Fis;
}
// 0.5 * t(manualPart) %*% t(modelPart) %*% estimatorPartHessian %*% modelPart %*% manualPart

// Dense - dense - sparse
// [[Rcpp::export]]
arma::mat FisherInformation_inner_cpp_DDS(
    const arma::mat& estimator,
    const arma::mat& model,
    const arma::sp_mat& manual
) {
  // Dense part first:
  arma::mat dense = model.t() * estimator * model;
  
  // Now full:
  arma::mat Fis = 0.5 * manual.t() * dense * manual;
  
  return Fis;
}

