#ifndef DVAR1_H
#define DVAR1_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_phi_theta_var1_cpp(
        const Rcpp::List& grouplist
);


#endif
