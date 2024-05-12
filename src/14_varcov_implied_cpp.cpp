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
Rcpp::List implied_varcov_cpp(
    const S4& model,
    bool all = false
){
  S4 sample = model.slot("sample");
  Rcpp::List means = sample.slot("means");
  
  
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);
  
  // Add implied cov structure:
  Rcpp::List types = model.slot("types");
  std::string type = types["y"];
  x = impliedcovstructures_cpp(x, "", type, all);
  
  int nGroup = x.length();
  int g;
  
  for (g=0; g<nGroup; g++){
    bool proper = true;
    
    Rcpp::List grouplist = x[g];
    
    // If mu is not there, obtain from the sample means:
    if (!grouplist.containsElementNamed("mu")){
      arma::vec groupmean = means[g];
      grouplist["mu"] = groupmean;
    }
    
    
    // Check kappa:
    if (!grouplist.containsElementNamed("kappa")){
      arma::mat sigma = grouplist["sigma"];
      arma::mat kappa = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
      grouplist["kappa"] = kappa;
    }
    
    
    // Return properness:
    grouplist["proper"] = proper;
    
    x[g] = grouplist;
  }
  
  return(x);
}
