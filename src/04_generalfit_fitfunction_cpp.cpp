// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "05_MLestimator_fit_Gauss_cpp.h"
#include "05_MLestimator_fit_Gauss2L_cpp.h"
#include "05_MLestimator_fit_Ising.h"
#include "06_ULS_fitfunction_cpp.h"
#include "07_FIMLestimator_fitfunction_cppversion.h"
#include "04_generalfit_optimWorkspace.h"
#include "09_PenMLestimator_fit_Gauss_cpp.h"
#include "09_PenMLestimator_fit_Ising_cpp.h"
#include "09_PenMLestimator_fit_FIML_Gauss_cpp.h"

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
        // Implied matrices are improper (not positive definite): the
        // estimators below are not guaranteed to be bounded from below in
        // that case (e.g. the ML trace/Mahalanobis terms with an indefinite
        // pseudo-inverse), which lets optimizers diverge. Return a large
        // penalty so the optimizer backs away from the improper region:
        return(1e20);
      }
    }
    
  }
  
  
  
  if (estimator == "ML"){

    if (distribution == "Gaussian"){
      fit = maxLikEstimator_Gauss_cpp(prep);
    } else if (distribution == "Spin"){
      fit = maxLikEstimator_Ising_cpp(prep);
    } else if (distribution == "TwoLevelGaussian"){
      fit = maxLikEstimator_Gauss2L_cpp(prep);
    } else {
      Rcpp::stop("Distribution not supported for ML estimator.");
    }

  } else if (estimator == "PML"){

    const OptimWorkspace& ws = getOrBuildWorkspace(model);
    if (distribution == "Gaussian"){
      fit = penMaxLikEstimator_Gauss_cpp(prep, x, ws);
    } else if (distribution == "Spin"){
      fit = penMaxLikEstimator_Ising_cpp(prep, x, ws);
    } else {
      Rcpp::stop("Distribution not supported for PML estimator.");
    }

  } else if (estimator == "ULS" || estimator == "WLS" || estimator == "DWLS"){

    fit = ULS_Gauss_cpp(prep);

  } else if (estimator == "PFIML"){

    const OptimWorkspace& ws = getOrBuildWorkspace(model);
    fit = penFIMLEstimator_Gauss_cpp(prep, x, ws);

  } else if (estimator == "FIML"){

    fit = fimlestimator_Gauss_cpp(prep);

  } else {
    Rcpp::stop("Estimator not supported.");
  }

  return(fit);
}
