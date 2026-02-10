#ifndef PENMLFITISING_H
#define PENMLFITISING_H

#include <RcppArmadillo.h>
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double penMaxLikEstimator_Ising_cpp(
    const Rcpp::List& prep,
    const arma::vec& x,
    const OptimWorkspace& ws
);

#endif
