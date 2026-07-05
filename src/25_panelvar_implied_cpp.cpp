// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "04_generalfit_optimWorkspace.h"
#include "03_modelformation_PDC_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Implied model for panelvar. This is the dlvm1 implied structure with
// lambda = I: no factor loadings, no residual variances, observed means.
// Core implementation that takes pre-formed model matrices:
Rcpp::List implied_panelvar_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  int i, t, tt;

  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;

  std::string within_latent = types["within_latent"];
  std::string between_latent = types["between_latent"];
  std::string temporal = "raw";
  if (types.containsElementNamed("temporal")) temporal = as<std::string>(types["temporal"]);

  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "zeta_within", within_latent, all);
  x = impliedcovstructures_cpp(x, "zeta_between", between_latent, all);

  int nGroup = x.length();
  int g;

  // General stuff:
  const Rcpp::List& extramats = ws.extramatrices;
  arma::mat design = extramats["design"];

  int nVar = design.n_rows;
  int nTime = design.n_cols;

  for (g=0; g<nGroup; g++){
    bool proper = true;

    Rcpp::List grouplist = x[g];

    // PDC temporal parameterization: compute beta from PDC and the
    // innovation covariance:
    if (temporal == "PDC"){
      arma::mat PDCmat = grouplist["PDC"];
      arma::mat sigma_zw_tmp = grouplist["sigma_zeta_within"];
      grouplist["beta"] = PDC_to_beta_cpp(PDCmat, sigma_zw_tmp);
    }

    // Model matrices:
    arma::mat mu = grouplist["mu"];
    arma::mat beta = grouplist["beta"];

    arma::mat sigma_zeta_within = grouplist["sigma_zeta_within"];
    arma::mat sigma_zeta_between = grouplist["sigma_zeta_between"];

    // Beta star:
    arma::mat I2 = eye(nVar*nVar, nVar*nVar);
    arma::mat BetaStar = inv(I2 - kron(beta, beta));

    // Implied mean vector (stationary means, repeated at every wave):
    arma::vec impMu = mu;

    arma::vec fullMu(nTime * nVar);
    for (i=0;i<nTime;i++){
      fullMu.submat(i*nVar,0,(i+1)*nVar-1,0) = impMu;
    }

    Rcpp::List allSigmas_within(nTime);
    arma::mat stationary_sigma_within = matrixform(BetaStar * vectorise(sigma_zeta_within));

    allSigmas_within[0] = stationary_sigma_within;

    // Fill implied:
    if (nTime > 1){
      for (t=1;t<nTime;t++){
        arma::mat lastsigma = allSigmas_within[t-1];
        allSigmas_within[t] = beta * lastsigma;
      }
    }

    // Create the block Toeplitz (full within-person cov matrix):
    arma::mat fullSigma_within  = blockToeplitz_cpp(allSigmas_within);

    // Full between-person cov matrix:
    arma::mat fullSigma_between(nTime * nVar, nTime * nVar);

    for (t=0;t<nTime;t++){
      for (tt=0;tt<nTime;tt++){
        fullSigma_between.submat(t*nVar,tt*nVar,(t+1)*nVar-1,(tt+1)*nVar-1) = sigma_zeta_between;
      }
    }

    // Full implied covmat:
    arma::mat fullSigma = fullSigma_within + fullSigma_between;

    // Subset:
    uvec inds = find(vectorise(design) == 0);

    // Subset and add to the list:
    arma::mat finalsigma = fullSigma;
    finalsigma.shed_rows(inds);
    finalsigma.shed_cols(inds);
    fullMu.shed_rows(inds);

    grouplist["mu"] = fullMu;
    grouplist["sigma"] = finalsigma;

    // Precision:
    grouplist["kappa"] =  solve_symmetric_cpp_matrixonly_withcheck(finalsigma, proper);

    // Extra matrices needed in optimization:
    if (!all){
      grouplist["BetaStar"] = BetaStar;
      grouplist["allSigmas_within"] = allSigmas_within;
      grouplist["IkronBeta"] = kronecker_I_X(beta, nVar);
    } else {
      // Named as in dlvm1 (with lambda = I the observed-level and
      // "latent"-level structures coincide):
      grouplist["sigma_within"] = stationary_sigma_within;
      grouplist["sigma_between"] = sigma_zeta_between;
      grouplist["sigma_within_full"] = fullSigma_within;
      grouplist["sigma_eta_within"] = stationary_sigma_within;
      grouplist["sigma_eta_within_lag1"] = allSigmas_within[1];
      grouplist["sigma_crosssection"] = stationary_sigma_within + sigma_zeta_between;

      // Add PDC:
      arma::mat kappa_zeta_within;

      if (grouplist.containsElementNamed("kappa_zeta_within")){
        kappa_zeta_within = as<arma::mat>(grouplist["kappa_zeta_within"]);
      } else {
        kappa_zeta_within = solve_symmetric_cpp_matrixonly_withcheck(sigma_zeta_within, proper);
      }

      grouplist["PDC"] = computePDC_cpp(grouplist["beta"], kappa_zeta_within, grouplist["sigma_zeta_within"]);
    }

    // Return properness:
    grouplist["proper"] = proper;

    x[g] = grouplist;
  }

  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_panelvar_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_panelvar_cpp_core(x, model, all);
}
