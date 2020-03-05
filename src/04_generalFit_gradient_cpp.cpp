// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Dense - sparse - sparse
// [[Rcpp::export]]
arma::mat gradient_inner_cpp_DSS(
    const arma::mat& estimator,
    const arma::sp_mat& model,
    const arma::sp_mat& manual
) {
  // Sparse part first:
  arma::sp_mat sparse = model * manual;
  
  // Now second part:
  arma::vec grad = vectorise(estimator * sparse);
  
  // Return
  return grad;
}

// Dense - dense - sparse
// [[Rcpp::export]]
arma::mat gradient_inner_cpp_DDS(
    const arma::mat& estimator,
    const arma::mat& model,
    const arma::sp_mat& manual
) {
  arma::mat part1 = estimator * model;
  // Return
  arma::vec grad = vectorise(part1 * manual);
  
  return grad;
}

