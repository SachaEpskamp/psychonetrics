// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "04_generalfit_optimWorkspace.h"
#include "03_modelformation_impliedcovstructures.h"
#include "21_Ising_implied_cpp.h"
#include "21_Ising_helperfunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List prepare_Ising_cpp(
    arma::vec x,
    const S4& model
){
  int g;
  
  // Form model matrices using cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  Rcpp::List mats = formModelMatrices_direct(x, ws.mapping);

  // Compute implied matrices using core function:
  Rcpp::List imp = implied_Ising_cpp_core(mats, model, false);

  // sample slot (read from original model - unchanged by parameter updates):
  S4 sample = model.slot("sample");

  // Things needed:
  // bool corinput = sample.slot("corinput");
  // bool meanstructure = model.slot("meanstructure");

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
  Rcpp::List squares = sample.slot("squares");
  Rcpp::List means =  sample.slot("means");

  
  
  // Group models:
  Rcpp::List groupModels(nGroup);
  
  for (g=0; g<nGroup; g++){
    Rcpp::List grouplist = imp[g];
    growlist(grouplist, extramatrices);
    growlist(grouplist, types);
    
    grouplist["squares"] = squares[g];
    grouplist["means"] = means[g];
    grouplist["nobs"] = nPerGroup(g);
    
    
    // Compute expectation:
    Rcpp::List exp = isingExpectation(
      grouplist["omega"], grouplist["tau"],  grouplist["beta"],  grouplist["responses"], extramatrices["min_sum"]
    );
    
    growlist(grouplist, exp);
    
    groupModels[g] = grouplist;
  }
  
  
  
  Rcpp::List result;
  
  result["nPerGroup"] = nPerGroup;
  result["nTotal"] = nTotal;
  result["nGroup"] = nGroup;
  result["groupModels"] = groupModels;
  
  return(result);
}



