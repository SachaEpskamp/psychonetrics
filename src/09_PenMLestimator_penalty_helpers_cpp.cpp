// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Penalty gradient helper for Penalized Maximum Likelihood (PML) estimation
// Adds L2 gradient + L1 subgradient to the already-computed ML gradient

#include <RcppArmadillo.h>
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

void addPenaltyGradient_cpp(
    arma::vec& grad,
    const arma::vec& x,
    const OptimWorkspace& ws
){
  // L2 gradient: lambda * (1-alpha) * x
  grad += (ws.penalty_lambda_vec * (1.0 - ws.penalty_alpha)) % x;

  // L1 subgradient: lambda * alpha * sign(x) for x != 0, 0 for x == 0
  grad += ws.penalty_lambda_vec * ws.penalty_alpha % arma::sign(x);
}
