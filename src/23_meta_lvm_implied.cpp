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

// Core implementation that takes pre-formed model matrices
// Combines LVM implied (sigma_y from lambda, beta, sigma_zeta, sigma_epsilon)
// with meta-analytic structure (mu = vech(sigma_y), sigma = sigma_randomEffects + V)
Rcpp::List implied_meta_lvm_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  int s;

  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;

  std::string latent = types["latent"];
  std::string residual = types["residual"];
  std::string randomEffects = types["randomEffects"];

  // Add implied covariance structures:
  x = impliedcovstructures_cpp(x, "zeta", latent, all);
  x = impliedcovstructures_cpp(x, "epsilon", residual, all);
  x = impliedcovstructures_cpp(x, "randomEffects", randomEffects, all);

  int g = 0; // meta_lvm is single-group

  bool proper = true;

  // General stuff from workspace:
  const Rcpp::List& extramats = ws.extramatrices;
  const arma::vec& nPerGroupVec = ws.nPerGroup;

  std::string est = Rcpp::as<std::string>(extramats["Vestimation"]);

  arma::mat V = Rcpp::as<arma::mat>(extramats["V"]);
  Rcpp::List Vall = extramats["Vall"];

  Rcpp::List grouplist = x[g];

  // ---- LVM implied covariance (same as lvm_implied_cpp_core) ----
  // Model matrices:
  arma::mat beta = grouplist["beta"];
  arma::mat lambda = grouplist["lambda"];
  arma::mat sigma_zeta = grouplist["sigma_zeta"];
  arma::mat sigma_epsilon = grouplist["sigma_epsilon"];

  // BetaStar = (I - beta)^{-1}
  int n = beta.n_rows;
  arma::mat I = eye(n, n);
  arma::mat BetaStar = inv(I - beta);
  arma::mat Lambda_BetaStar = lambda * BetaStar;

  // Store intermediates for derivatives (when !all):
  if (!all){
    grouplist["BetaStar"] = BetaStar;
    grouplist["Lambda_BetaStar"] = Lambda_BetaStar;

    arma::mat Betasta_sigmaZeta = BetaStar * sigma_zeta;
    arma::mat tBetakronBeta = kron(BetaStar.t(), BetaStar);

    grouplist["Betasta_sigmaZeta"] = Betasta_sigmaZeta;
    grouplist["tBetakronBeta"] = tBetakronBeta;
  }

  // Factor part:
  arma::mat factorPart = Lambda_BetaStar * sigma_zeta * Lambda_BetaStar.t();

  // Implied sigma_y (sigma_epsilon diagonal is free):
  arma::mat sigma_y = factorPart + sigma_epsilon;

  // Force symmetric:
  sigma_y = 0.5 * (sigma_y + sigma_y.t());

  // Store sigma_y:
  grouplist["sigma_y"] = sigma_y;

  // ---- Meta-analytic structure (same as meta_varcov_implied_cpp_core) ----
  arma::mat sigma_randomEffects = grouplist["sigma_randomEffects"];

  if (est == "averaged"){
    // mu = vech(sigma_y) with diagonal (covariance input):
    arma::vec mu = vech(sigma_y, true);
    grouplist["mu"] = mu;

    // sigma = sigma_randomEffects + V:
    arma::mat sigma = sigma_randomEffects + V;
    grouplist["sigma"] = sigma;
    grouplist["kappa"] = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);

  } else {
    int nStudy = (int)nPerGroupVec(g);

    // Shared mu across studies:
    arma::vec mu = vech(sigma_y, true);

    Rcpp::List mulist(nStudy);
    Rcpp::List sigmalist(nStudy);
    Rcpp::List kappalist(nStudy);

    for (s = 0; s < nStudy; s++){
      mulist[s] = mu;

      arma::mat Vcur = Vall[s];
      arma::mat cursigma = sigma_randomEffects + Vcur;
      sigmalist[s] = cursigma;
      kappalist[s] = solve_symmetric_cpp_matrixonly_withcheck(cursigma, proper);
    }

    grouplist["mu"] = mulist;
    grouplist["sigma"] = sigmalist;
    grouplist["kappa"] = kappalist;
  }

  grouplist["proper"] = proper;

  x[g] = grouplist;

  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_meta_lvm_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_meta_lvm_cpp_core(x, model, all);
}
