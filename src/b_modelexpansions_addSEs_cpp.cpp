// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
S4 addSEs_cpp(
    const S4& xOld
){
  
  // Setup:
  int i;
  
  // Copy model:
  S4 x = clone(xOld);
  
  
  
  // If not computed, warn user!
  bool computed = x.slot("computed");
  if (!computed){
    Rf_warning("Model was not computed, interpret standard errors and p-values with care!");
  }
  
  // Parameter table:
  Rcpp::List parameters = x.slot("parameters");
  
  // Parameter indices and estimates:
  NumericVector est = parameters["est"];
  NumericVector absest = abs(est);
  NumericVector se = parameters["se"];
  NumericVector p = parameters["p"];
  IntegerVector par = parameters["par"];
  

  
  int nparTotal = par.length();
  int nparFree = max(par);
  
  
  // If no constrained parameters, nothing to do!
  if (nparFree == 0){
    
    return(x);
    
  }
  
  
  
  // Clear old:
  for (i=0; i<nparTotal;i++){
    se[i] = NA_REAL;
    p[i] = NA_REAL;
  }

    // Sample size:
    S4 sample = x.slot("sample");
    Rcpp::List groups = sample.slot("groups");
    NumericVector nobs = groups["nobs"];
    double n = sum(nobs);
    
    // Assuming information has been computed ....
    arma::mat Fisher = x.slot("information");
    arma::mat Hinv = 1.0/n * solve_symmetric_cpp_matrixonly(Fisher);
      
      
    // Obtain SEs
    arma::vec SEs = sqrt(abs(Hinv.diag()));
    
    // 
    // 
    // x =sqrt(abs(diag(solve_symmetric(H))))
    // sqrt(2/n) * x -
    // SEs

    
    // Add standard errors:
    for (i=0; i<nparTotal;i++){
      if (par[i] > 0){
        se[i] = SEs(par[i]-1);
        p[i] = 2.0 * R::pnorm(absest[i], 0.0, se[i], 0, 0);
      }
    }
    
    
    // Write back:
    parameters["se"] = se;
    parameters["p"] = p;
    x.slot("parameters") = parameters;
    
      // Return:
      return(x);
      
}
