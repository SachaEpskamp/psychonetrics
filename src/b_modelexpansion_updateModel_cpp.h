#ifndef UDPATEMODEL_H
#define UDPATEMODEL_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

S4 updateModel_cpp(
                arma::vec x,
                const S4& model,
                bool updateMatrices
);


#endif
