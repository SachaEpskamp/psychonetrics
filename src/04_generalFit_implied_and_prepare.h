#ifndef IMPPREP_H
#define IMPPREP_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List impliedModel_cpp(
    const S4& model,
    bool all = false
);

Rcpp::List prepareModel_cpp(
    arma::vec x,
    const S4& model
);

#endif
