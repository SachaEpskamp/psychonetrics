#ifndef GENERALGRAD_H
#define GENERALGRAD_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::vec psychonetrics_gradient_cpp(
    arma::vec x,
    const S4& model
);

arma::vec psychonetrics_gradient_cpp_prepared(
    Rcpp::List prep,
    arma::sp_mat manualPart
);

#endif