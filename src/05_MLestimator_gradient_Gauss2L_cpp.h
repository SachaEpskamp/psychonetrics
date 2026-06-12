#ifndef MLGRADGAUSS2L_H
#define MLGRADGAUSS2L_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat jacobian_gaussian2L_sigma_group_cpp(
    const Rcpp::List& grouplist
);

arma::mat jacobian_gaussian2L_sigma_cpp(
    const Rcpp::List& prep
);

#endif
