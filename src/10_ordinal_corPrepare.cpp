// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// General function to do the following:
// 1. Compute means for continuous variables, thresholds for ordinal
// 2. Compute Pearson/polychoric/polyserial correlations
// 3. Compute WLS weights matrix
// [[Rcpp::export]]
List corPrepare_cpp(
    List Data, // Data as data frame
    LogicalVector isOrdered
) { 
  
  
  
  List Result;
  return(Result);
}
