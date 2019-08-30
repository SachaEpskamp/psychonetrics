// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// FIXME: I want to put all these functions nicely in seperate files, but can't figure out how and get weird errors.
// So my solution is that I put all helper functions in this file. 


// Silly wrapper for bivariate distribution to return 0 / 1 when needed:
double mypbinorm(double u1,double u2,double rho){
  double res = 0.0;
  
  // Special cases:
  // One is < -50 (basically -Inf), just return 0
  if (u1 < -50 || u2 < -50){
    return(0.0);
  }
  
  // If both are above 50, return 1:
  if (u1 > 50 && u2 > 50){
    return(1.0);
  }
  
  // If u1 > 50, just do univariate
  if (u1 > 50){
    return(R::pnorm(u2,0.0,1.0,1,0));
  }
  
  // If u2 > 50, just do univariate
  if (u2 > 50){
    return(R::pnorm(u1,0.0,1.0,1,0));
  }
  
  // Finally, return the bivariate one if needed:
  res = pbv::pbv_rcpp_pbvnorm0(u1,u2,rho);      
  return(res);
}

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
double polychoric_fit_summary(double rho, NumericMatrix tab, NumericVector t1, NumericVector t2){
  // Iterators:
  int i, j;
  
  // Double used later:
  double pi, tabmult;
  
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
  t1_aug[0] = -99; // R_NegInf;
  t2_aug[0] = -99; // R_NegInf;
  
  // Fill last with Inf:
  t1_aug[nLevels1] = 99; // R_PosInf;
  t2_aug[nLevels2] = 99; // R_PosInf;
  
  // Fill rest with thresholds:
  for (i=1;i<nLevels1;i++){
    t1_aug[i] = t1[i-1];
  }
  for (i=1;i<nLevels2;i++){
    t2_aug[i] = t2[i-1];
  }
  
  
  // Compute sample size:
  double n=0.0;
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
      pi = mypbinorm(t1_aug[i+1],t2_aug[j+1], rho) - 
        mypbinorm(t1_aug[i],t2_aug[j+1], rho) -
        mypbinorm(t1_aug[i+1],t2_aug[j], rho) + 
        mypbinorm(t1_aug[i],t2_aug[j], rho);
      
      
      // Adjust table to 0.5 if it is zero:
      // if (tab(i,j) == 0 && nLevels1 == 2 && nLevels2 == 2){
      //   tabmult = 1.0;
      // } else {
      tabmult = (double)tab(i,j);
      // }
      
      curL +=  tabmult * log(pi);
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
double polychoric_grad_summary(double rho, NumericMatrix tab, NumericVector t1, NumericVector t2){
  // FIXME: Copypaste from fit function, could be nicer...
  // Iterators:
  int i, j;
  
  // Double used later:
  double pi, num, tabmult;
  
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
  t1_aug[0] = -99; // R_NegInf;
  t2_aug[0] = -99; // R_NegInf;
  
  // Fill last with Inf:
  t1_aug[nLevels1] = 99; // R_PosInf;
  t2_aug[nLevels2] = 99; // R_PosInf;
  
  // Fill rest with thresholds:
  for (i=1;i<nLevels1;i++){
    t1_aug[i] = t1[i-1];
  }
  for (i=1;i<nLevels2;i++){
    t2_aug[i] = t2[i-1];
  }
  
  
  // Compute sample size:
  double n=0.0;
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
      pi = mypbinorm(t1_aug[i+1],t2_aug[j+1], rho) - 
        mypbinorm(t1_aug[i],t2_aug[j+1], rho) -
        mypbinorm(t1_aug[i+1],t2_aug[j], rho) + 
        mypbinorm(t1_aug[i],t2_aug[j], rho);
      
      // Adjust table to 0.5 if it is zero:
      // if (tab(i,j) == 0 && nLevels1 == 2 && nLevels2 == 2){
      //   tabmult = 1.0;
      // } else {
      tabmult = (double)tab(i,j);
      
      // To prevent numerical issues, set pi to be at minimum a very small number:
      pi = std::max(pi, 0.000001);
      
      curGrad +=  tabmult * num / pi;
    }
  }
  
  // Multiply with -2/n:
  curGrad *= -2.0 / (double)n;
  
  return curGrad;
}

// FUCK THIS!
// 
// ////Test from https://cran.rstudio.com/web/packages/RcppNumerical/vignettes/introduction.html
// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::depends(RcppNumerical)]]
// 
// #include <RcppNumerical.h>
// 
// using namespace Numer;
// 
// // f = 100 * (x2 - x1^2)^2 + (1 - x1)^2
// // True minimum: x1 = x2 = 1
// class Rosenbrock: public MFuncGrad
// {
// public:
//   double f_grad(Constvec& x, Refvec grad)
//   {
//     double t1 = x[1] - x[0] * x[0];
//     double t2 = 1 - x[0];
//     grad[0] = -400 * x[0] * t1 - 2 * t2;
//     grad[1] = 200 * t1;
//     return 100 * t1 * t1 + t2 * t2;
//   }
// };
// 
// // [[Rcpp::export]]
// Rcpp::List optim_test()
// {
//   Eigen::VectorXd x(2);
//   x[0] = -1.2;
//   x[1] = 1;
//   double fopt;
//   Rosenbrock f;
//   int res = optim_lbfgs(f, x, fopt);
//   return Rcpp::List::create(
//     Rcpp::Named("xopt") = x,
//     Rcpp::Named("fopt") = fopt,
//     Rcpp::Named("status") = res
//   );
// }
// 
// 
// //class for a polychoric correlation
// class Polychoric: public MFuncGrad
// {
// public:
//   double f_grad(Constvec& rho, Refvec grad)
//   {
//     // Fill in gradient here:
//     grad[0] = polychoric_grad_summary(rho[0], tab, t1, t2);
//     
//     // Fill in fit function here:
//     double fit = polychoric_fit_summary(rho[0], tab, t1, t2);
//     
//     // return:
//     return fit;
//   }
// };
// 
// // [[Rcpp::export]]
// Rcpp::List optim_test2(IntegerMatrix tab, NumericVector t1, NumericVector t2)
// {
//   Eigen::VectorXd rho(1);
//   rho[0] = 0.0;
//   double fopt;
//   Polychoric f;
//   int res = optim_lbfgs(f, rho, fopt);
//   return Rcpp::List::create(
//     Rcpp::Named("xopt") = rho,
//     Rcpp::Named("fopt") = fopt,
//     Rcpp::Named("status") = res
//   );
// }

// 
// My own Gradient Descent optimizer:
// [[Rcpp::export]]
double estimate_polychoric(IntegerVector y1, IntegerVector y2, NumericVector t1, NumericVector t2,
                           double tol = 0.000001, double stepsize = 1, int maxIt = 1000,
                           double zeroAdd = 0.5){
  // Table:
  IntegerMatrix tab = cpp_table(y1,y2);
  
  
  int i,j;
  int n1 = tab.nrow();
  int n2 = tab.ncol();
  bool zeroMargin = false;
  
  for (i=0;i<n1;i++){
    for (j=0;j<n2;j++){
      if (tab(i,j) == 0){
        Rf_warning("Zero frequency cell found in cross-table! Polychoric estimate might be unstable!");
        zeroMargin = true;
      }
    }  
  }
  
  // Make this a numeric matrix instead:
  NumericMatrix tabNumeric(n1, n2);
  for (i=0;i<n1;i++){
    for (j=0;j<n2;j++){
      tabNumeric(i,j) = (double)tab(i,j);
    }  
  }
  
  
  // For a 2x2 table with zero margins, adjust:
  // From https://github.com/yrosseel/lavaan/blob/master/R/lav_polychor.R#L327-L336
  if (zeroMargin == true && n1 == 2 && n2 == 2){
    if (tab(0,0) == 0 || tab(1,1) == 0){
      tabNumeric(0,0) += zeroAdd;
      tabNumeric(1,1) += zeroAdd;
      tabNumeric(1,0) -= zeroAdd;
      tabNumeric(0,1) -= zeroAdd;
    } else {
      tabNumeric(0,0) -= zeroAdd;
      tabNumeric(1,1) -= zeroAdd;
      tabNumeric(1,0) += zeroAdd;
      tabNumeric(0,1) += zeroAdd;
    }
  }
  
  
  // Current iteration:
  int curIt = 0;
  // Start value:
  double rho = 0.0;
  
  // Some doubles I'll need:
  double curFit, newFit, curGrad, newGrad, delta;
  
  double gamma = stepsize;
  
  curFit = newFit = polychoric_fit_summary(rho, tabNumeric, t1, t2);
  curGrad = newGrad = polychoric_grad_summary(rho, tabNumeric, t1, t2);
  
  // Start iterating:
  do {
    curIt++;
    curFit = newFit;
    curGrad = newGrad;
    
    // Gradient descent step:
    delta = -1.0 *  gamma * curGrad;
    
    // Check if out of bounds, and make smaller step if needed:
    while (rho + delta < -1 || rho + delta > 1){
      delta *= 0.5;
    }

    rho += delta;
    
    // Update fit:
    newFit = polychoric_fit_summary(rho, tabNumeric, t1, t2);
    
    // Update gradient:
    newGrad = polychoric_grad_summary(rho, tabNumeric, t1, t2);
    
    // Update step size:
    gamma = std::abs(delta * (newGrad - curGrad)) / std::pow(newGrad - curGrad, 2.0);
    
  } while (curIt < maxIt && std::abs(newGrad) > tol && rho > -1.0 && rho < 1.0);
  
  if (curIt >= maxIt){
    Rf_error("Polychoric correlation estimator did not converge.");
  }
  
  return(rho);
}

// Gradient of fit to threshold for a single subject
// [[Rcpp::export]]
double threshold_grad_singlesubject(int y, int j, NumericVector t_aug){
  
  double res = 0;
  // int nThresh = t.length();
  // 
  // // Make augmented thresholds:
  // NumericVector t_aug(nThresh + 2);
  // 
  // // Fill first with -Inf:
  // t_aug[0] = -99; // R_NegInf;
  // 
  // // Fill last with Inf:
  // t_aug[nThresh + 1] = 99; // R_PosInf;
  // 
  // // FILL THE REST TOO
  // for (int i = 0; i < nThresh; i++){
  //   t_aug[i+1] = t[i];
  // }
  // 
  double t_0 = t_aug[y];
  double t_1 = t_aug[y+1];
  
  
  if (y == j){
    res = R::dnorm(t_1,0.0,1.0,0) / (R::pnorm(t_1,0.0,1.0,1,0) - R::pnorm(t_0,0.0,1.0,1,0));
  } else if (y-1 == j){
    res =  -1.0 * R::dnorm(t_0,0.0,1.0,0) / (R::pnorm(t_1,0.0,1.0,1,0) - R::pnorm(t_0,0.0,1.0,1,0));
  }
  return res;
}

// Gradient of fit to polychoric correlation for a single subject
// [[Rcpp::export]]
double polychor_grad_singlesubject(int y1, int y2, double rho, NumericVector t_aug1, NumericVector t_aug2, double pi){
  
  double res = 0;
  // int nThresh1 = t1.length();
  // int nThresh2 = t2.length();
  // 
  // // Make augmented thresholds:
  // NumericVector t_aug1(nThresh1 + 2);
  // NumericVector t_aug2(nThresh2 + 2);
  // 
  // 
  // // Fill first with -Inf:
  // t_aug1[0] = -99; // R_NegInf;
  // t_aug2[0] = -99; // R_NegInf;
  // 
  // // Fill last with Inf:
  // t_aug1[nThresh1 + 1] = 99; // R_PosInf;
  // t_aug2[nThresh2 + 1] = 99; // R_PosInf;
  // 
  // // FILL THE REST TOO
  // for (int i = 0; i < nThresh1; i++){
  //   t_aug1[i+1] = t1[i];
  // }
  // for (int i = 0; i < nThresh2; i++){
  //   t_aug2[i+1] = t2[i];
  // }
  // 
  double t_01 = t_aug1[y1];
  double t_11 = t_aug1[y1+1];
  double t_02 = t_aug2[y2];
  double t_12 = t_aug2[y2+1];
  
  // Compute:
  // Numerator:
  double num = binormal_density(t_11,t_12, rho) - 
    binormal_density(t_01,t_12, rho) -
    binormal_density(t_11,t_02, rho) + 
    binormal_density(t_01,t_02, rho);
  
  // // Denominator:
  // double pi = mypbinorm(t_11,t_12, rho) - 
  //   mypbinorm(t_01,t_12, rho) -
  //   mypbinorm(t_11,t_02, rho) + 
  //   mypbinorm(t_01,t_02, rho);
  
  
  res = num / pi;
  
  return res;
}


// Gradient of BIVARIATE fit to threshold (t1) for a single subject
// [[Rcpp::export]]
double bthreshold_grad_singlesubject(int y1, int y2, double rho, int tIndex, NumericVector t_aug1, NumericVector t_aug2, double pi){
  
  // Empty Gradient:
  double grad = 0.0;
  
  // The gradient will only be nonzero if y1 == tIndex (y is below) or y1 == tIndex + 1 (above):
  if (y1 != tIndex && y1 != tIndex + 1){
    return(grad);
  } else {
    
    double t_01 = t_aug1[y1];
    double t_11 = t_aug1[y1+1];
    double t_02 = t_aug2[y2];
    double t_12 = t_aug2[y2+1];

    // double pi = mypbinorm(t_11,t_12, rho) - 
    //   mypbinorm(t_01,t_12, rho) -
    //   mypbinorm(t_11,t_02, rho) + 
    //   mypbinorm(t_01,t_02, rho);
    // 
    double rhopart = std::pow(1.0 - std::pow(rho, 2.0), 0.5);
    
    // If tIndex == y1 (y is below), the gradient becomes:
    if (y1 == tIndex){
      grad = 1.0/pi * R::dnorm(t_11,0.0,1.0,0) * (
        R::pnorm((t_12 - rho * t_11) / rhopart,0.0,1.0,1,0) - 
          R::pnorm((t_02 - rho * t_11) / rhopart,0.0,1.0,1,0)
      );      
    }
    
    if (y1 == tIndex+1){
      grad = 1.0/pi * R::dnorm(t_01,0.0,1.0,0) * ( -1.0 *
        R::pnorm((t_12 - rho * t_01) / rhopart,0.0,1.0,1,0) +
          R::pnorm((t_02 - rho * t_01) / rhopart,0.0,1.0,1,0)
      );
    }
    return(grad);
  }
  
}


// Function to compute the bivariate likelihood for an ORDINAL variable:
double ordered_bivariate_likelihood(int y1, int y2, double rho, NumericVector t_aug1, NumericVector t_aug2){
  
  double t_01 = t_aug1[y1];
  double t_11 = t_aug1[y1+1];
  double t_02 = t_aug2[y2];
  double t_12 = t_aug2[y2+1];
  
  double pi = mypbinorm(t_11,t_12, rho) - 
    mypbinorm(t_01,t_12, rho) -
    mypbinorm(t_11,t_02, rho) + 
    mypbinorm(t_01,t_02, rho);
  

  // To prevent numerical issues, set pi to be at minimum a very small number:
  pi = std::max(pi, 0.000001);
  
  return pi;
}




