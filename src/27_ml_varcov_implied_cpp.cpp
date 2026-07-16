// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Core implementation that takes pre-formed model matrices. C++ twin of
// R/27_ml_varcov_implied.R: forms the within/between covariance structures and
// stores mu, sigma_within and sigma_between (symmetrized). What else is
// produced depends on the estimator:
// - "ML" (two-level sufficient statistics): nothing else is needed.
// - "FIML": additionally the wide-format (one row per cluster) mean vector and
//   covariance matrix fullSigma = I_nMax (x) sigma_within +
//   J_nMax (x) sigma_between, subset by the design pattern, plus its inverse,
//   stored under mu / sigma / kappa (mirrors the ml_lvm wide format).
Rcpp::List implied_ml_varcov_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all
){
  // Read the types and estimator from the cached workspace:
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

  // Wide format needed for the FIML estimator:
  const bool wide = (ws.estimator == "FIML");
  arma::uvec obsInds;
  bool completeDesign = true;
  int nVar = 0, nMaxInCluster = 0;
  if (wide){
    const Rcpp::List& extramats = ws.extramatrices;
    arma::mat design = extramats["designPattern"];
    nVar = design.n_rows;
    nMaxInCluster = design.n_cols;

    // Indices (column-major over the design matrix: position t outer, variable
    // i inner) of the observed wide-format positions, as in the R twin's
    // as.vector(designPattern) == 1 subsetting:
    arma::uvec allInds(nVar * nMaxInCluster);
    int nObs = 0;
    for (int t = 0; t < nMaxInCluster; t++){
      for (int i = 0; i < nVar; i++){
        if (design(i, t) == 1){
          allInds(nObs) = t * nVar + i;
          nObs++;
        }
      }
    }
    obsInds = allInds.head(nObs);
    completeDesign = (nObs == nVar * nMaxInCluster);
  }

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

    if (wide){
      // Wide-format mean vector and covariance matrix, subset by the design
      // (also under all = true, so that getmatrix() returns the wide sigma /
      // kappa / mu exactly as it does for ml_lvm FIML models):
      arma::vec mu_p = arma::vectorise(mu);
      arma::vec fullMu(nMaxInCluster * nVar);
      for (int t = 0; t < nMaxInCluster; t++){
        fullMu.subvec(t * nVar, (t + 1) * nVar - 1) = mu_p;
      }

      // fullSigma = I (x) SW + J (x) SB:
      arma::mat fullSigma = (arma::mat)kronecker_I_X(sw, nMaxInCluster);
      for (int t = 0; t < nMaxInCluster; t++){
        for (int tt = 0; tt < nMaxInCluster; tt++){
          fullSigma.submat(t * nVar, tt * nVar, (t + 1) * nVar - 1, (tt + 1) * nVar - 1) += sb;
        }
      }

      if (!completeDesign){
        arma::vec subMu = fullMu(obsInds);
        arma::mat subSigma = fullSigma.submat(obsInds, obsInds);
        fullMu = subMu;
        fullSigma = subSigma;
      }
      fullSigma = 0.5 * (fullSigma + fullSigma.t());

      grouplist["mu"] = fullMu;
      grouplist["sigma"] = fullSigma;

      // Precision (with positive-definiteness flag, as the ml_lvm wide format;
      // the generic C++ fit/gradient return the penalty value when !proper):
      bool proper = true;
      grouplist["kappa"] = solve_symmetric_cpp_matrixonly_withcheck(fullSigma, proper);
      grouplist["proper"] = proper;
    }

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
