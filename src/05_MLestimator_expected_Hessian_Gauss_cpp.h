#ifndef HESMLGAUSS_H
#define HESMLGAUSS_H

#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat expected_hessian_Gaussian_cpp(
    const Rcpp::List& prep
);


#endif
