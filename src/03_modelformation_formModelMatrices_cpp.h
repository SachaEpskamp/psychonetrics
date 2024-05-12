#ifndef MODELMATS_H
#define MODELMATS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List formModelMatrices_cpp(
                const S4& model
);


#endif
