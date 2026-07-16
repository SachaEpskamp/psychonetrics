#ifndef MLVARCOVIMPLIEDCPP_H
#define MLVARCOVIMPLIEDCPP_H

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List implied_ml_varcov_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all
);

Rcpp::List implied_ml_varcov_cpp(
    const S4& model,
    bool all
);

#endif
