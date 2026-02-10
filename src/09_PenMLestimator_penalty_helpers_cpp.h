#ifndef PENML_PENALTY_HELPERS_H
#define PENML_PENALTY_HELPERS_H

#include <RcppArmadillo.h>
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Add L2 gradient + L1 subgradient to the gradient vector (modifies in-place)
void addPenaltyGradient_cpp(
    arma::vec& grad,
    const arma::vec& x,
    const OptimWorkspace& ws
);

#endif
