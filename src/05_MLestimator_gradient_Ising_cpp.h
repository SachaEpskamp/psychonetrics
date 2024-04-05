#ifndef DMLISING_H
#define DMLISING_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat jacobian_Ising_cpp(
    const Rcpp::List& prep
);


#endif
