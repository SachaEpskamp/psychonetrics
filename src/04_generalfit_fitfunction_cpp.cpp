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
double psychonetrics_fitfunction_cpp(
    const arma::vec& x,
    const S4& model
){
  double fit;
  
  // Prepare model:
  Rcpp::List prep = prepareModel_cpp(x, model);
  
  // What estimator:
  std::string estimator= prep["estimator"];
  
  
  // What distribution:
  std::string distribution = prep["distribution"];
  
  // Group models:
  Rcpp::List groupmodels = prep["groupModels"];
  
  // Loop over:
  int nGroup = groupmodels.length();
  for (int g=0; g<nGroup; g++){
    Rcpp::List grouplist = groupmodels[g];
    
    // check if element proper is there:
    if (grouplist.containsElementNamed("proper")){
      bool proper = grouplist["proper"];
      if (!proper){
        // Rf_PrintValue(wrap("a matrix was not positive definite in optimization... fit function may not be accurate."));
        // FIXME:
        // Rf_warning("A matrix was not positive definite... fit function may not be accurate. Check your results for consistency using other optimizers with setoptimizer.");
        // return(1e20);
        // return(NA_REAL);
      }
    }
    
  }
  
  
  
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
