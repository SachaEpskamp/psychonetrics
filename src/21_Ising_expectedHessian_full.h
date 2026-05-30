#ifndef ISING_EXPECTEDHESSIAN_FULL_H
#define ISING_EXPECTEDHESSIAN_FULL_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat expected_hessian_Ising_group_full_cpp(
    const arma::mat& omega,
    const arma::vec& tau,
    const arma::vec& delta,
    double beta,
    const arma::vec& responses,
    double min_sum
);

arma::mat expected_hessian_Ising_full_cpp(const Rcpp::List& prep);

#endif
