// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "04_generalfit_optimWorkspace.h"
#include "03_modelformation_impliedcovstructures.h"
#include "22_ml_lvm_implied_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List prepare_ml_lvm_cpp(
    arma::vec x,
    const S4& model
){
  int g, i;

  // Form model matrices using cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  Rcpp::List mats = formModelMatrices_direct(x, ws.mapping);

  // Two-level sufficient-statistics ML estimator (twin of the twolevel_ML
  // flag in prepare_ml_lvm)?
  bool twolevel_ML = (ws.estimator == "ML" && ws.distribution == "TwoLevelGaussian");

  // Compute implied matrices using core function:
  Rcpp::List imp = implied_ml_lvm_cpp_core(mats, model, false, twolevel_ML);

  // Read constant data from cached workspace (no S4 slot reads):
  bool corinput = ws.corinput;
  bool meanstructure = ws.meanstructure;
  int nVar = ws.nVar;
  int nGroup = ws.nGroup;
  arma::vec nPerGroup = ws.nPerGroup;
  const Rcpp::List& extramatrices = ws.extramatrices;
  const Rcpp::List& types = ws.types;
  double nTotal = ws.nTotal;
  const Rcpp::List& S = ws.sampleCovs;
  const Rcpp::List& means = ws.sampleMeans;
  const Rcpp::List& thresholds = ws.sampleThresholds;

  // Two-level sufficient statistics (per group), guarded for objects saved
  // before the 'twolevel' slot existed (twin of the guard in prepare_ml_lvm):
  const Rcpp::List& twolevelStats = ws.sampleTwolevel;
  if (twolevel_ML && twolevelStats.length() != nGroup){
    Rcpp::stop("estimator = 'ML' for ml_lvm requires two-level sufficient statistics, which are missing from the model. Rebuild the model with ml_lvm(..., estimator = 'ML').");
  }

  // Group models:
  Rcpp::List groupModels(nGroup);
  
  for (g=0; g<nGroup; g++){
    Rcpp::List grouplist = imp[g];
    growlist(grouplist, extramatrices);
    growlist(grouplist, types);
    
    grouplist["S"] = S[g];
    grouplist["means"] = means[g];
    grouplist["corinput"] = corinput;
    grouplist["meanstructure"] = meanstructure;

    if (twolevel_ML){
      grouplist["twolevel"] = twolevelStats[g];
    }
    
    // Tau:
    if (!grouplist.containsElementNamed("tau")){
      arma::mat tau(1,nVar);
      for (i=0; i<nVar; i++){
        tau(0,i) = NA_REAL;
      }
      grouplist["tau"] = tau;
    }
    
    if (thresholds.length() > 0){
      grouplist["thresholds"] = thresholds[g];
    }
    
    groupModels[g] = grouplist;
  }
  
  
  
  Rcpp::List result;
  
  result["nPerGroup"] = nPerGroup;
  result["nTotal"] = nTotal;
  result["nGroup"] = nGroup;
  result["groupModels"] = groupModels;
  // Flag so that the model Jacobian dispatch (d_phi_theta_ml_lvm_cpp) picks
  // the two-level variant with rows [mu; vech Sigma_W; vech Sigma_B]:
  result["twolevel_ML"] = twolevel_ML;

  return(result);
}



