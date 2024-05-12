#ifndef GENERALFIT_H
#define GENERALFIT_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double psychonetrics_fitfunction_cpp(
    const arma::vec& x,
    const S4& model
);

#endif
