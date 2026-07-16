// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Core implementation that takes pre-formed model matrices. C++ twin of
// R/27_ml_varcov_implied.R: ml_varcov is always the two-level sufficient-
// statistics model, so this simply forms the within/between covariance
// structures and stores mu, sigma_within and sigma_between (symmetrized). There
// is no wide-format / FIML branch and no latent layer.
Rcpp::List implied_ml_varcov_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all
){
  // Read the two types from the cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;

  std::string within = types["within"];
  std::string between = types["between"];

  // Implied covariance structures for the two blocks (writes sigma_within /
  // sigma_between plus the derivative helpers delta_IminOinv_<block> etc. when
  // all = false), exactly as impliedcovstructures does in the R twin:
  x = impliedcovstructures_cpp(x, "within", within, all);
  x = impliedcovstructures_cpp(x, "between", between, all);

  int nGroup = x.length();

  for (int g = 0; g < nGroup; g++){
    Rcpp::List grouplist = x[g];

    arma::mat mu = grouplist["mu"];
    arma::mat sigma_within = grouplist["sigma_within"];
    arma::mat sigma_between = grouplist["sigma_between"];

    // Force symmetric (as in the R twin):
    arma::mat sw = 0.5 * (sigma_within + sigma_within.t());
    arma::mat sb = 0.5 * (sigma_between + sigma_between.t());

    grouplist["mu"] = mu;
    grouplist["sigma_within"] = sw;
    grouplist["sigma_between"] = sb;

    if (all){
      grouplist["sigma_crosssection"] = sw + sb;
    }

    x[g] = grouplist;
  }

  return(x);
}

// Forms matrices from the S4 model, then computes the implied structures:
// [[Rcpp::export]]
Rcpp::List implied_ml_varcov_cpp(
    const S4& model,
    bool all = false
){
  Rcpp::List x = formModelMatrices_cpp(model);
  return implied_ml_varcov_cpp_core(x, model, all);
}
