#ifndef DWLS_H
#define DWLS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat ULS_gradient_Gauss_cpp(
    const Rcpp::List& prep
);


#endif
