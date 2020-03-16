#ifndef RCPPHELPERS_H
#define RCPPHELPERS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::vec eig_sym_cpp(
        arma::mat X
);

Rcpp::List solve_symmetric_cpp(
        arma::mat X,
        bool logdet,
        double epsilon
);

arma::mat solve_symmetric_cpp_matrixonly(
        arma::mat X,
        double epsilon
);


arma::mat bdiag_psychonetrics(
                Rcpp::List mats
);

arma::vec vech(
                arma::mat& X,
                bool diag = true
);


#endif