#ifndef RCPP_H
#define RCPP_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

bool sparse_or_dense_cpp(
    const arma::mat& X
);

int diag_sparse_dense_cpp(
    const arma::mat& X
);

#endif
