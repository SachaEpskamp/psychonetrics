#ifndef MLFITGAUSS_H
#define MLFITGAUSS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double maxLikEstimator_Gauss_cpp(
    const Rcpp::List& prep
);

#endif
