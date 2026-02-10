#ifndef PENMLFITGAUSS_H
#define PENMLFITGAUSS_H

#include <RcppArmadillo.h>
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double penMaxLikEstimator_Gauss_cpp(
    const Rcpp::List& prep,
    const arma::vec& x,
    const OptimWorkspace& ws
);

#endif
