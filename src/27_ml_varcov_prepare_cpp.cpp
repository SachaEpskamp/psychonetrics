// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "04_generalfit_optimWorkspace.h"
#include "27_ml_varcov_implied_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// C++ twin of prepare_ml_varcov (R). Builds the per-group model lists consumed
// by the two-level Gauss2L estimator and the ml_varcov model Jacobian.
// [[Rcpp::export]]
Rcpp::List prepare_ml_varcov_cpp(
    arma::vec x,
    const S4& model
){
  int g;

  // Form model matrices using the cached workspace mapping:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  Rcpp::List mats = formModelMatrices_direct(x, ws.mapping);

  // Implied structures (two-level: mu, sigma_within, sigma_between + helpers):
  Rcpp::List imp = implied_ml_varcov_cpp_core(mats, model, false);

  // Constant data from the workspace:
  int nGroup = ws.nGroup;
  arma::vec nPerGroup = ws.nPerGroup;
  const Rcpp::List& extramatrices = ws.extramatrices;
  const Rcpp::List& types = ws.types;
  double nTotal = ws.nTotal;
  const Rcpp::List& S = ws.sampleCovs;
  const Rcpp::List& means = ws.sampleMeans;

  // Two-level sufficient-statistics ML estimator? (The FIML estimator uses the
  // wide-format matrices from the implied function and the per-pattern FIML
  // data added by the generic prepare wrapper instead.)
  const bool twolevel_ML = (ws.estimator == "ML");

  // Two-level sufficient statistics per group (ML estimator only):
  const Rcpp::List& twolevelStats = ws.sampleTwolevel;
  if (twolevel_ML && twolevelStats.length() != nGroup){
    Rcpp::stop("estimator = 'ML' for ml_varcov requires two-level sufficient statistics, which are missing from the model. Rebuild the model with ml_varcov(...).");
  }

  // Group models:
  Rcpp::List groupModels(nGroup);
  for (g = 0; g < nGroup; g++){
    Rcpp::List grouplist = imp[g];
    growlist(grouplist, extramatrices);
    growlist(grouplist, types);

    grouplist["S"] = S[g];
    grouplist["means"] = means[g];
    if (twolevel_ML){
      grouplist["twolevel"] = twolevelStats[g];
    }

    groupModels[g] = grouplist;
  }

  Rcpp::List result;
  result["nPerGroup"] = nPerGroup;
  result["nTotal"] = nTotal;
  result["nGroup"] = nGroup;
  result["groupModels"] = groupModels;
  // Flag so that the model Jacobian dispatch picks the two-level variant
  // rather than the wide FIML variant:
  result["twolevel_ML"] = twolevel_ML;

  return(result);
}
