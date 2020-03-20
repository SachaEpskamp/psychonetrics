// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "05_MLestimator_fit_Gauss_cpp.h"
#include "05_MLestimator_fit_Ising.h"
#include "06_ULS_fitfunction_cpp.h"
#include "07_FIMLestimator_fitfunction_cppversion.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double psychonetrics_fitfunction_cpp_prepared(
  Rcpp::List prep
){
  double fit;
  
  
  // What estimator:
  std::string estimator= prep["estimator"];
  
  
  // What distribution:
  std::string distribution = prep["distribution"];
  
  
  
  if (estimator == "ML"){
    
    if (distribution == "Gaussian"){
      fit = maxLikEstimator_Gauss_cpp(prep);
    } else if (distribution == "Ising"){
      
      fit = maxLikEstimator_Ising_cpp(prep);
      
    } else {
      Rf_error("Distribution not supported for ML estimator.");
    }
    
    
  } else if (estimator == "ULS" || estimator == "WLS" || estimator == "DWLS"){
    
    fit = ULS_Gauss_cpp(prep);
    
  } else if (estimator == "FIML"){
    
    fit = fimlestimator_Gauss_cpp(prep);
    
  } else {
    Rf_error("Estimator not supported.");
  }

  
  
  return(fit);
}


// [[Rcpp::export]]
double psychonetrics_fitfunction_cpp(
    arma::vec x,
    const S4& model
){
  // Prepare model:
  Rcpp::List prep = prepareModel_cpp(x, model);
  
  double fit = psychonetrics_fitfunction_cpp_prepared(prep);
  
  return(fit);
}