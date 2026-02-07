#ifndef RCPP_H
#define RCPP_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::sp_mat kronecker_I_X(
    const arma::mat& X,
    int n
);

arma::sp_mat kronecker_X_I(
    const arma::mat& X,
    int n
);

// Dense versions (return arma::mat directly, faster when result is used as dense):
arma::mat kronecker_I_X_dense(
    const arma::mat& X,
    int n
);

arma::mat kronecker_X_I_dense(
    const arma::mat& X,
    int n
);

arma::sp_mat kronecker_diag_sparse(
    arma::sp_mat X
);

arma::sp_mat kronecker_diag(
        arma::mat X
);

#endif
