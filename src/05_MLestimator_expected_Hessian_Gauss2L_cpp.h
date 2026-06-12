#ifndef MLEXPHESSGAUSS2L_H
#define MLEXPHESSGAUSS2L_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat expected_hessian_Gauss2L_group_cpp(
    const Rcpp::List& grouplist
);

arma::mat expected_hessian_Gauss2L_cpp(
    const Rcpp::List& prep
);

#endif
