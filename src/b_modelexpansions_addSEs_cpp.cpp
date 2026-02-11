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
    const S4& xOld,
    bool verbose = true,
    bool approximate_SEs = false
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

/*
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
    */

    arma::mat Hinv;

    if (approximate_SEs){
      Hinv = 1.0/n * solve_symmetric_cpp_matrixonly(Fisher, 1.490116e-08, true);
    } else {
      Hinv = 1.0/n * solve_symmetric_cpp_matrixonly(Fisher, 1.490116e-08, false);
    }

    // Check for NA:
    double first_element = (double)Hinv(0,0);

    if (R_IsNA(first_element)){

      if (!approximate_SEs){
        // Auto-fallback: try approximate SEs before giving up
        Hinv = 1.0/n * solve_symmetric_cpp_matrixonly(Fisher, 1.490116e-08, true);
        first_element = (double)Hinv(0,0);

        if (R_IsNA(first_element)){
          // Even approximate inversion failed:
          Rf_warning("Standard errors could not be obtained because the Fischer information matrix could not be inverted. This may be a symptom of a non-identified model or due to convergence issues.");

          // Write NAs:
          for (i=0; i<nparTotal;i++){
            if (par[i] > 0){
              se[i] = NA_REAL;
              p[i] = NA_REAL;
            }
          }
        } else {
          // Approximate inversion succeeded:
          Rf_warning("Exact standard errors could not be obtained because the Fischer information matrix could not be inverted. Falling back to approximate standard errors. This can occur with zero cells in crosstables (common in multi-group Ising models) or near-boundary estimates. Interpret with care.");

          // Obtain SEs
          arma::vec SEs = sqrt(abs(Hinv.diag()));

          // Add standard errors:
          for (i=0; i<nparTotal;i++){
            if (par[i] > 0){
              se[i] = SEs(par[i]-1);
              p[i] = 2.0 * R::pnorm(absest[i], 0.0, se[i], 0, 0);
            }
          }
        }

      } else {
        // User explicitly requested approximate SEs and it still failed:
        Rf_warning("Standard errors could not be obtained because the Fischer information matrix could not be inverted. This may be a symptom of a non-identified model or due to convergence issues.");

        // Write NAs:
        for (i=0; i<nparTotal;i++){
          if (par[i] > 0){
            se[i] = NA_REAL;
            p[i] = NA_REAL;
          }
        }
      }

    } else {

      // Obtain SEs
      arma::vec SEs = sqrt(abs(Hinv.diag()));

      // Add standard errors:
      for (i=0; i<nparTotal;i++){
        if (par[i] > 0){
          se[i] = SEs(par[i]-1);
          p[i] = 2.0 * R::pnorm(absest[i], 0.0, se[i], 0, 0);
        }
      }

    }


    // Write back:
    parameters["se"] = se;
    parameters["p"] = p;
    x.slot("parameters") = parameters;
    
      // Return:
      return(x);
      
}
