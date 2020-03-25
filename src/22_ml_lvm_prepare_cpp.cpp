// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "b_modelexpansion_updateModel_cpp.h"
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
  
  // New model:
  S4 newMod = updateModel_cpp(x,model,false);
  
  // sample slot:
  S4 sample = newMod.slot("sample");
  
  // Things needed:
  bool corinput = sample.slot("corinput");
  bool meanstructure = newMod.slot("meanstructure");
  
  // Groups table:
  Rcpp::List group = sample.slot("groups");
  
  // Variables table:
  Rcpp::List variables = sample.slot("variables");
  arma::vec vars = variables["id"];
  int nVar = vars.n_elem;
  
  // Number of groups:
  arma::vec id = group["id"];
  int nGroup = id.n_elem;
  
  // Number of observations per group:
  arma::vec nPerGroup = group["nobs"];
  
  // Compute implied matrices:
  Rcpp::List imp = implied_ml_lvm_cpp(newMod, false);
  
  // Extra matrices:
  Rcpp::List extramatrices = model.slot("extramatrices");
  
  // Types:
  Rcpp::List types = model.slot("types");
  
  // Total sample:
  double nTotal = sum(nPerGroup);
        
  // Sample stats:
  Rcpp::List S = sample.slot("covs");
  Rcpp::List means =  sample.slot("means");
  Rcpp::List thresholds = sample.slot("thresholds");

  
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
  
  return(result);
}



