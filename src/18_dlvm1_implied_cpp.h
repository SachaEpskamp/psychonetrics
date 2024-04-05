#ifndef DLVM1_IMPLIED_H
#define DLVM1_IMPLIED_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List implied_dlvm1_cpp(
    const S4& model,
    bool all = false
);

#endif
