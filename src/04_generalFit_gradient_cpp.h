#ifndef GENERALGRAD_H
#define GENERALGRAD_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::vec psychonetrics_gradient_cpp(
    arma::vec x,
    const S4& model,
    bool sparsemodel = false,
    bool useM = false
);

void psychonetrics_gradient_cpp_inner(
    const arma::vec& x,
    arma::vec& grad,
    const S4& model,
    bool sparsemodel = false,
    bool useM = false
);

#endif
