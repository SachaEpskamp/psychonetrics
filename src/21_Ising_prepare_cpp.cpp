// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "b_modelexpansion_updateModel_cpp.h"
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
  
  // New model:
  S4 newMod = updateModel_cpp(x,model,false);
  
  // sample slot:
  S4 sample = newMod.slot("sample");
  
  // Things needed:
  // bool corinput = sample.slot("corinput");
  // bool meanstructure = newMod.slot("meanstructure");
  
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
  
  // Compute implied matrices:
  Rcpp::List imp = implied_Ising_cpp(newMod, false);
  
  // Extra matrices:
  Rcpp::List extramatrices = model.slot("extramatrices");
  
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
    
    grouplist["squares"] = squares[g];
    grouplist["means"] = means[g];
    grouplist["nobs"] = nPerGroup(g);
    
    // Compute expectation:
    Rcpp::List exp = isingExpectation(
      grouplist["omega"], grouplist["tau"],  grouplist["beta"],  grouplist["responses"]
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



