// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::mat fimlEstimator_Gauss_subgroup_cpp(arma::mat sigma, arma::mat kappa) {
  return sigma * kappa;
}
