#ifndef ORDINALCCPP_H
#define ORDINALCCPP_H

#include <Rcpp.h>
#include <math.h>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double computeMean(Rcpp::NumericVector y);

Rcpp::NumericVector computeThresholds(Rcpp::IntegerVector y);

double pearsonCov(Rcpp::NumericVector y1, Rcpp::NumericVector y2, double mean1, double mean2, bool unbiased = false);

IntegerVector toOrdinal(NumericVector var);

int maxInt(IntegerVector y);

IntegerMatrix cpp_table(
    IntegerVector y1,
    IntegerVector y2
);

double polychoric_fit_summary(double rho, IntegerMatrix tab, NumericVector t1, NumericVector t2);

double binormal_density(double x1, double x2, double rho, double sigma1 = 1.0, double sigma2 = 1.0, double mu1 = 0.0, double mu2 = 0.0);


#endif