#ifndef PDC_HELPERS_H
#define PDC_HELPERS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat PDC_to_beta_cpp(
    const arma::mat& PDC,
    const arma::mat& sigma
);

// Fills Tmat (d vec beta / d vec PDC) and Xmat (d vec beta / d theta_c):
void PDC_reparam_cpp(
    const arma::mat& PDC,
    const arma::mat& beta,
    const arma::mat& sigma,
    const arma::mat& aug,
    const arma::sp_mat& D,
    const arma::sp_mat& C,
    arma::mat& Tmat,
    arma::mat& Xmat
);

#endif
