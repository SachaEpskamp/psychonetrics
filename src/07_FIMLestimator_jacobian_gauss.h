#ifndef DFIML_H
#define DFIML_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat jacobian_fiml_gaussian_sigma_cpp(
    const Rcpp::List& prep
);


#endif
