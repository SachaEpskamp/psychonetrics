#ifndef DMLGAUSS_H
#define DMLGAUSS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat jacobian_gaussian_sigma_cpp(
    const Rcpp::List& prep
);


#endif
