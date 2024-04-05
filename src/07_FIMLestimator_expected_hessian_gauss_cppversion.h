#ifndef HESFIML_H
#define HESFIML_H

#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



arma::mat expected_hessian_fiml_Gaussian_cppVersion(
    const Rcpp::List& prep
);


#endif
