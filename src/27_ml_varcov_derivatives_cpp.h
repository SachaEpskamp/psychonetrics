#ifndef MLVARCOVDERIVCPP_H
#define MLVARCOVDERIVCPP_H

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_phi_theta_ml_varcov_group_cpp(
    const Rcpp::List& grouplist
);

arma::mat d_phi_theta_ml_varcov_wide_group_cpp(
    const Rcpp::List& grouplist
);

arma::mat d_phi_theta_ml_varcov_cpp(
    const Rcpp::List& prep
);

#endif
