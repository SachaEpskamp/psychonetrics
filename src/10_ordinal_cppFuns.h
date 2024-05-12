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

double estimate_polychoric(IntegerVector y1, IntegerVector y2, NumericVector t1, NumericVector t2,
                    double tol = 0.00000001, double stepsize = 1, int maxIt = 1000,
                    double zeroAdd = 0.5);

double threshold_grad_singlesubject(int y, int j, NumericVector t);

double polychor_grad_singlesubject(int y1, int y2, double rho, NumericVector t1, NumericVector t2, double pi);


double bthreshold_grad_singlesubject(int y1, int y2, double rho, int tIndex, NumericVector t1, NumericVector t2, double pi);

double ordered_bivariate_likelihood(int y1, int y2, double rho, NumericVector t_aug1, NumericVector t_aug2);


#endif
