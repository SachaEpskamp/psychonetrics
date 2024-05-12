#ifndef ISINGHELPERS_H
#define ISINGHELPERS_H

#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


double computeZ_cpp(
                arma::mat graph,
                arma::vec tau,
                double beta,
                arma::vec responses,
                double min_sum
);


Rcpp::List isingExpectation(
    arma::mat graph,
    arma::vec tau,
    double beta,
    arma::vec responses,
    double min_sum
);


#endif
