// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Penalized ML fit function for Gaussian distribution
// Wraps the standard ML fit and adds L1 + L2 penalty terms

#include <RcppArmadillo.h>
#include "05_MLestimator_fit_Gauss_cpp.h"
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double penMaxLikEstimator_Gauss_cpp(
    const Rcpp::List& prep,
    const arma::vec& x,
    const OptimWorkspace& ws
){
  // Base ML fit:
  double fit = maxLikEstimator_Gauss_cpp(prep);

  // Add penalty terms:
  double alpha = ws.penalty_alpha;
  const arma::vec& lam = ws.penalty_lambda_vec;

  // L2 (ridge): 0.5 * sum(lambda_i * (1-alpha) * x_i^2)
  fit += 0.5 * arma::dot(lam * (1.0 - alpha), x % x);

  // L1 (LASSO): sum(lambda_i * alpha * |x_i|)
  fit += arma::dot(lam * alpha, arma::abs(x));

  return fit;
}
