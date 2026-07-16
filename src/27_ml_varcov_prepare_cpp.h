#ifndef MLVARCOVPREPARECPP_H
#define MLVARCOVPREPARECPP_H

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List prepare_ml_varcov_cpp(
    arma::vec x,
    const S4& model
);

#endif
