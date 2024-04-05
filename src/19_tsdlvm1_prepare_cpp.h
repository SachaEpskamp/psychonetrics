#ifndef TSDLVM1_PREPARE_H
#define TSDLVM1_PREPARE_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List prepare_tsdlvm1_cpp(
    arma::vec x,
    const S4& model
);

#endif
