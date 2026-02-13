// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "04_generalfit_optimWorkspace.h"
#include "03_modelformation_impliedcovstructures.h"
#include "23_meta_lvm_implied_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List prepare_meta_lvm_cpp(
    arma::vec x,
    const S4& model
){
  int g;

  // Form model matrices using cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  Rcpp::List mats = formModelMatrices_direct(x, ws.mapping);

  // Compute implied matrices using core function:
  Rcpp::List imp = implied_meta_lvm_cpp_core(mats, model, false);

  // Read constant data from cached workspace (no S4 slot reads):
  bool corinput = ws.corinput;
  bool meanstructure = ws.meanstructure;
  int nGroup = ws.nGroup;
  arma::vec nPerGroup = ws.nPerGroup;
  const Rcpp::List& extramatrices = ws.extramatrices;
  const Rcpp::List& types = ws.types;
  double nTotal = ws.nTotal;
  const Rcpp::List& S = ws.sampleCovs;
  const Rcpp::List& means = ws.sampleMeans;


  // Group models:
  Rcpp::List groupModels(nGroup);

  for (g=0; g<nGroup; g++){
    Rcpp::List grouplist = imp[g];
    growlist(grouplist, extramatrices);
    growlist(grouplist, types);

    grouplist["S"] = S[g];
    grouplist["means"] = means[g];
    grouplist["corinput"] = false;       // Meta-level data is not corinput
    grouplist["meanstructure"] = meanstructure;
    grouplist["metacor"] = corinput;     // Original study-level corinput flag

    // Override D with D_c: the ML gradient uses D for the meta-level sigma (nCov x nCov),
    // not the LVM-level D (nNode x nNode):
    grouplist["D"] = extramatrices["D_c"];

    groupModels[g] = grouplist;
  }



  Rcpp::List result;

  result["nPerGroup"] = nPerGroup;
  result["nTotal"] = nTotal;
  result["nGroup"] = nGroup;
  result["groupModels"] = groupModels;

  return(result);
}
