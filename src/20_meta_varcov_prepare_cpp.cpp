// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_formModelMatrices_direct.h"
#include "03_modelformation_impliedcovstructures.h"
#include "20_meta_varcov_implied_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List prepare_meta_varcov_cpp(
    arma::vec x,
    const S4& model
){
  int g;
  
  // Form model matrices directly (no S4 clone needed):
  Rcpp::List mats = formModelMatrices_cpp_direct(x, model);

  // Compute implied matrices using core function:
  Rcpp::List imp = implied_meta_varcov_cpp_core(mats, model, false);

  // sample slot (read from original model - unchanged by parameter updates):
  S4 sample = model.slot("sample");

  // Things needed:
  bool corinput = sample.slot("corinput");
  bool meanstructure = model.slot("meanstructure");

  // Groups table:
  Rcpp::List group = sample.slot("groups");

  // Variables table:
  Rcpp::List variables = sample.slot("variables");
  arma::vec vars = variables["id"];
  // int nVar = vars.n_elem;

  // Number of groups:
  arma::vec id = group["id"];
  int nGroup = id.n_elem;

  // Number of observations per group:
  arma::vec nPerGroup = group["nobs"];
  
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
    grouplist["corinput"] = false;
    grouplist["meanstructure"] = meanstructure;
    grouplist["metacor"] = corinput;
    
   
    groupModels[g] = grouplist;
  }
  
  
  
  Rcpp::List result;
  
  result["nPerGroup"] = nPerGroup;
  result["nTotal"] = nTotal;
  result["nGroup"] = nGroup;
  result["groupModels"] = groupModels;
  
  return(result);
}



