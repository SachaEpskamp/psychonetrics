#ifndef MODMAT_H
#define MODMAT_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::sp_mat Mmatrix_cpp(
    Rcpp::DataFrame parDF
);

arma::sp_mat Mmatrix_cpp_list(
    Rcpp::List parDF
);

#endif
