#ifndef VAR1_IMPLIED_H
#define VAR1_IMPLIED_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List implied_var1_cpp(
    const S4& model,
    bool all = false
);

Rcpp::List implied_var1_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
);

#endif
