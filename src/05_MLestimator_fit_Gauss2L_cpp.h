#ifndef MLFITGAUSS2L_H
#define MLFITGAUSS2L_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double minustwo_logl_Gauss2L_noconstant_cpp(
    const arma::vec& mu,
    const arma::mat& SW,
    const arma::mat& SB,
    const Rcpp::List& twolevel
);

double maxLikEstimator_Gauss2L_group_cpp(
    const Rcpp::List& grouplist
);

double maxLikEstimator_Gauss2L_cpp(
    const Rcpp::List& prep
);

#endif
