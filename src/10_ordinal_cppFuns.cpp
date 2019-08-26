// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
using namespace Rcpp;
using namespace arma;
#include <RcppNumerical.h>
using namespace Numer;

// FIXME: I want to put all these functions nicely in seperate files, but can't figure out how and get weird errors.
// So my solution is that I put all helper functions in this file. 

// Silly function to compute the mean as well:
// [[Rcpp::export]]
double computeMean(
    NumericVector y
){
  // Iterators
  int p;
  
  // Remove NA:
  y = y[!is_na(y)];
  
  // sample size:
  int n = y.length();
  
  // Mean:
  double mean = 0;
  for (p=0;p<n;p++){
    mean += y[p];
  }
  mean *= (1.0/(double)n);
  return(mean);
}


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



// Function to compute the ML pearon covariance (or variance):
// [[Rcpp::export]]
double pearsonCov(
    Rcpp::NumericVector y1, Rcpp::NumericVector y2, 
    double mean1, double mean2, bool unbiased = false
){
  
  // Iterators
  int p;
  
  // Remove NA:
  LogicalVector noMis = !is_na(y1) & !is_na(y2);
  y1 = y1[noMis];
  y2 = y2[noMis];
  
  // sample size:
  int n = y1.length();
  
  if (y2.length() != n){
    Rf_error("Length of y1 is not equal to length of y2.");
  }
  
  // FIXME: Instead of this, let's use the supplied means already, as just like with thresholds,
  // we treat the first order statistics as fixed in computing the ML estimate of the second order
  // statistics. <- all this stuff is due to pairwise.
  
  // Compute means (FIXME: Not sure if this is needed, as means are computed before,
  // but those are based on all observations of each vector. e.g., if y1 has more missings
  // than y2, then mean of y1 will be based on something else).
  // double mean1 = computeMean(y1);
  // double mean2 = computeMean(y2);
  
  
  
  // Cov:
  double cov = 0;
  
  // Compute:
  for (p=0;p<n;p++){
    cov += (y1[p] - mean1) * (y2[p] - mean2);
  }
  
  // Divide by n:
  if (unbiased){
    cov *= (1.0/(double)(n-1));  
  } else {
    cov *= (1.0/(double)n);
  }
  
  
  // Return:
  return(cov);
}

// Silly maximum function:
int maxInt(IntegerVector y){
  int n = y.length();
  int max = 0;
  for (int i=0; i<n; i++){
    if (y[i] > max){
      max = y[i];
    }
  }
  return(max);
}

// Function to transform numeric varianle in ordinal levels (starting at 0):
// [[Rcpp::export]]
IntegerVector toOrdinal(
    NumericVector var
){
  // N subjects:
  int nCase = var.length();
  
  // Iterators:
  int p, j;
  
  // Missings:
  LogicalVector curNA = is_na(var);
  
  // All unique values sorted
  NumericVector uniqueVals = unique(na_omit(var));
  
  // Sort these:
  std::sort(uniqueVals.begin(), uniqueVals.end());
  
  // Get the size:
  int nLevelCur = uniqueVals.size();
  
  // Integervector to store data:
  IntegerVector ordinalVar(nCase);
  
  // Fill per person:
  for (p=0;p<nCase;p++){
    if (curNA[p] == true){
      ordinalVar[p] = NA_INTEGER;
    } else {
      for (j=0;j<nLevelCur;j++){
        if (var[p] == uniqueVals[j]){
          ordinalVar[p] = j;
        }
      } 
    }
  }
  
  return(ordinalVar);
}
// 
// // // POLYCHORIC CORRELATION ESTIMATOR FUNCTIONS //
// [[Rcpp::export]]
IntegerMatrix cpp_table(
    IntegerVector y1,
    IntegerVector y2
){
  // y1 and y2 must be correctly scored, integers starting at 0.
  // Remove NAs:
  LogicalVector noMis = !is_na(y1) & !is_na(y2);
  y1 = y1[noMis];
  y2 = y2[noMis];
  
  // sample size:
  int n = y1.length();
  
  if (y2.length() != n){
    Rf_error("Length of y1 is not equal to length of y2.");
  }
  
  // iterators:
  int p;
  
  // Number of levels:
  int nLevels1 = maxInt(y1) + 1;
  int nLevels2 = maxInt(y2) + 1;
  
  // Empty table (starting at zeroes):
  IntegerMatrix table(nLevels1, nLevels2);
  
  // Fill the table:
  for (p=0;p<n;p++){
    table(y1[p],y2[p])++;
  }
  
  // Return:
  return(table);
}

// polychoric correlation fit function, as function of rho, based on SUMMARY STATISTICS (used in optimization)
// [[Rcpp::export]]
double polychoric_fit_summary(double rho, IntegerMatrix tab, NumericVector t1, NumericVector t2){
  // Iterators:
  int i, j;
  
  // Double used later:
  double pi;
  
  // Levels:
  int nLevels1 = tab.nrow();
  int nLevels2 = tab.ncol();
  
  // Check size of t1 and t2:
  if (t1.length() != nLevels1 - 1){
    Rf_error("Length of thresholds of variable 1 is not equal to number of levels - 1");
  }
  if (t2.length() != nLevels2 - 1){
    Rf_error("Length of thresholds of variable 2 is not equal to number of levels - 1");
  }
  
  // Make augmented thresholds:
  NumericVector t1_aug(nLevels1 + 1);
  NumericVector t2_aug(nLevels2 + 1);
  
  // Fill first with -Inf:
  t1_aug[0] = -9999; // R_NegInf;
  t2_aug[0] = -9999; // R_NegInf;
  
  // Fill last with Inf:
  t1_aug[nLevels1] = 9999; // R_PosInf;
  t2_aug[nLevels2] = 9999; // R_PosInf;
  
  // Fill rest with thresholds:
  for (i=1;i<nLevels1;i++){
    t1_aug[i] = t1[i-1];
  }
  for (i=1;i<nLevels2;i++){
    t2_aug[i] = t2[i-1];
  }
  
  
  // Compute sample size:
  int n=0;
  for (i=0;i<nLevels1;i++){
    for (j=0;j<nLevels2;j++){
      n += tab(i,j);
    }
  }
  
  
  // Current sum:
  double curL = 0.0;
  
  // Loop over all levels:
  for (i=0;i<nLevels1;i++){
    for (j=0;j<nLevels2;j++){
      
      
      // Probabiliy:
      pi = pbv::pbv_rcpp_pbvnorm0(t1_aug[i+1],t2_aug[j+1], rho) - 
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i],t2_aug[j+1], rho) -
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i+1],t2_aug[j], rho) + 
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i],t2_aug[j], rho);
      
      
      curL +=  tab(i,j) * log(pi);
    }
  }
  
  // Multiply with -2/n:
  curL *= -2.0 / (double)n;
  
  return curL;
}

// Simple bivariate normal function that does not use arma or matrix multiplication (because I am lazy),
// Taken from https://www.statisticshowto.datasciencecentral.com/bivariate-normal-distribution/
// [[Rcpp::export]]
double binormal_density(double x1, double x2, double rho, double sigma1 = 1.0, double sigma2 = 1.0, double mu1 = 0.0, double mu2 = 0.0){
  
  double z = pow(x1 - mu1,2.0) / pow(sigma1, 2.0) - 2.0 * rho * (x1 - mu1) * (x2 - mu2) / (sigma1 * sigma2) + pow(x2 - mu2, 2.0) / pow(sigma2, 2.0);
  double res = 1.0 / (2.0 * M_PI * sigma1 * sigma2 * sqrt(1.0 - pow(rho, 2.0))) * exp(-1.0 * z / (2.0 * (1.0 - pow(rho, 2.0))));
  
  return res;
}

// polychoric correlation gradient, as function of rho, based on SUMMARY STATISTICS (used in optimization)
// [[Rcpp::export]]
double polychoric_grad_summary(double rho, IntegerMatrix tab, NumericVector t1, NumericVector t2){
  // FIXME: Copypaste from fit function, could be nicer...
  // Iterators:
  int i, j;
  
  // Double used later:
  double pi, num;
  
  // Levels:
  int nLevels1 = tab.nrow();
  int nLevels2 = tab.ncol();
  
  // Check size of t1 and t2:
  if (t1.length() != nLevels1 - 1){
    Rf_error("Length of thresholds of variable 1 is not equal to number of levels - 1");
  }
  if (t2.length() != nLevels2 - 1){
    Rf_error("Length of thresholds of variable 2 is not equal to number of levels - 1");
  }
  
  // Make augmented thresholds:
  NumericVector t1_aug(nLevels1 + 1);
  NumericVector t2_aug(nLevels2 + 1);
  
  // Fill first with -Inf:
  t1_aug[0] = -9999; // R_NegInf;
  t2_aug[0] = -9999; // R_NegInf;
  
  // Fill last with Inf:
  t1_aug[nLevels1] = 9999; // R_PosInf;
  t2_aug[nLevels2] = 9999; // R_PosInf;
  
  // Fill rest with thresholds:
  for (i=1;i<nLevels1;i++){
    t1_aug[i] = t1[i-1];
  }
  for (i=1;i<nLevels2;i++){
    t2_aug[i] = t2[i-1];
  }
  
  
  // Compute sample size:
  int n=0;
  for (i=0;i<nLevels1;i++){
    for (j=0;j<nLevels2;j++){
      n += tab(i,j);
    }
  }
  
  
  // Current sum:
  double curGrad = 0.0;
  
  // Loop over all levels:
  for (i=0;i<nLevels1;i++){
    for (j=0;j<nLevels2;j++){
      
      // Numerator:
      num = binormal_density(t1_aug[i+1],t2_aug[j+1], rho) - 
        binormal_density(t1_aug[i],t2_aug[j+1], rho) -
        binormal_density(t1_aug[i+1],t2_aug[j], rho) + 
        binormal_density(t1_aug[i],t2_aug[j], rho);
      
      // Denominator:
      pi = pbv::pbv_rcpp_pbvnorm0(t1_aug[i+1],t2_aug[j+1], rho) - 
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i],t2_aug[j+1], rho) -
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i+1],t2_aug[j], rho) + 
        pbv::pbv_rcpp_pbvnorm0(t1_aug[i],t2_aug[j], rho);
      
      
      curGrad +=  tab(i,j) * num / pi;
    }
  }
  
  // Multiply with -2/n:
  curGrad *= -2.0 / (double)n;
  
  return curGrad;
}


////Test from https://cran.rstudio.com/web/packages/RcppNumerical/vignettes/introduction.html
// P(0.3 < X < 0.8), X ~ Beta(a, b)
class BetaPDF: public Func
{
private:
  double a;
  double b;
public:
  BetaPDF(double a_, double b_) : a(a_), b(b_) {}
  
  double operator()(const double& x) const
  {
    return R::dbeta(x, a, b, 0);
  }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
{
  const double a = 3, b = 10;
  const double lower = 0.3, upper = 0.8;
  const double true_val = R::pbeta(upper, a, b, 1, 0) -
    R::pbeta(lower, a, b, 1, 0);
  
  BetaPDF f(a, b);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


