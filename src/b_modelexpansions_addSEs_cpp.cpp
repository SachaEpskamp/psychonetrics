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

    // Estimator:
    std::string estimator = Rcpp::as<std::string>(x.slot("estimator"));

    // For the least-squares estimators ULS and DWLS the weight matrix used in
    // estimation is not equal to Gamma^-1, so the naive (1/n) Info^-1 standard
    // errors are not robust. The robust sandwich
    //   (Delta'WDelta)^-1 Delta'W Gamma W Delta (Delta'WDelta)^-1 / n
    // is implemented in R (getVCOV). Delegate to it so the C++ and R SE paths
    // agree. For all other estimators (ML/FIML/WLS) the naive form below is
    // exact and unchanged.
    bool delegated = (estimator == "ULS" || estimator == "DWLS");

    arma::mat Hinv;

    if (delegated){

      // Robust VCOV computed in R:
      Rcpp::Environment pkg = Rcpp::Environment::namespace_env("psychonetrics");
      Rcpp::Function getVCOVfun = pkg["getVCOV"];
      Hinv = Rcpp::as<arma::mat>(getVCOVfun(x, Rcpp::Named("approximate_SEs") = approximate_SEs));

      // Check for NA:
      double first_element = (double)Hinv(0,0);

      if (R_IsNA(first_element)){
        // Inversion failed in R:
        Rf_warning("Standard errors could not be obtained because the information matrix could not be inverted. This may be a symptom of a non-identified model or due to convergence issues.");

        for (i=0; i<nparTotal;i++){
          if (par[i] > 0){
            se[i] = NA_REAL;
            p[i] = NA_REAL;
          }
        }
      } else {
        // Obtain SEs
        arma::vec SEs = sqrt(abs(Hinv.diag()));

        for (i=0; i<nparTotal;i++){
          if (par[i] > 0){
            se[i] = SEs(par[i]-1);
            p[i] = 2.0 * R::pnorm(absest[i], 0.0, se[i], 0, 0);
          }
        }
      }

    } else {

    // Assuming information has been computed ....
    arma::mat Fisher = x.slot("information");

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

    }


    // Write back:
    parameters["se"] = se;
    parameters["p"] = p;
    x.slot("parameters") = parameters;
    
      // Return:
      return(x);
      
}
