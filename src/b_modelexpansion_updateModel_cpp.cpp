// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
S4 updateModel_cpp(
    arma::vec x,
    const S4& model,
    bool updateMatrices
){
  S4 newMod = clone(model); // FIXME: creates a copy, but it avoids a ton of weird stuff happening otherwise... Could do better.
  
  int i;
    
  // Extract pars:
  Rcpp::List parsList = newMod.slot("parameters");

  // Extract ests:
  arma::vec est = parsList["est"];
  arma::vec par = parsList["par"];
  
  // max par:
  int maxPar = max(par);
  
  int totPar = est.n_elem;
  
  // Need to update?
  if (maxPar > 0){
    for (i=0; i<totPar;i++){
      if (par(i) > 0){
        est(i) = x(par(i)-1);
      }
    }
  }
  
  // Add model:
  if (updateMatrices){
    Rf_error("updateMatrices not yet implemented in C++");
    
  }
  
  // Write back:
  parsList["est"] = est;
  newMod.slot("parameters") = parsList;
  
  return(newMod);
}

