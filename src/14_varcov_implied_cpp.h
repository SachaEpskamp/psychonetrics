#ifndef VARCOV_IMPLIED_H
#define VARCOV_IMPLIED_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List implied_varcov_cpp(
    const S4& model,
    bool all = false
);

#endif
