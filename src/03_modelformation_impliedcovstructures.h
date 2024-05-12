#ifndef IMPCOVS_H
#define IMPCOVS_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

Rcpp::List impliedcovstructures_cpp(
                Rcpp::List x,
                std::string name = "",
                std::string type = "cov",
                bool all = false
);


#endif
