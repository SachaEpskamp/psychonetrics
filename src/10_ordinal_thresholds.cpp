// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Compute thresholds for one variable
// [[Rcpp::export]]
NumericVector computeThresholds(
    IntegerVector y // Data
  ) { 
  // Remove NA:
  y = y[!is_na(y)];
  
  int nSample = y.length();
  int maxLevel = max(y);
  int i, j;

  // Vector to store the occurances:
  IntegerVector table(maxLevel + 1, 0);
  for (i = 0; i < nSample; i++){
    table[y[i]]++;
  }
  
  // Empirical cumulative distribution (ingoring the last, not needed):
  NumericVector ECD(maxLevel, 0.0);
  for (i=0;i<maxLevel;i++){
    for (j=0;j<=i;j++){
      ECD[i] += (1.0/(double)nSample) * table[j];
    }
  }
  
  // Transform to quantiles for the thresholds:
  NumericVector thresholds = qnorm(ECD,0.0,1.0,1,0);
 
 return(thresholds);
}
