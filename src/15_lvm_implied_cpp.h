#ifndef LVM_IMPLIED_H
#define LVM_IMPLIED_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List implied_lvm_cpp(
    const S4& model,
    bool all = false
);

// Core version that takes pre-formed model matrices (skips formModelMatrices_cpp)
Rcpp::List implied_lvm_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
);

#endif
