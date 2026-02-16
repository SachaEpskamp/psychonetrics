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
// Combines VAR(1) implied (Sigma0, Sigma1 from beta, sigma_zeta)
// with meta-analytic structure (mu = [vech(Sigma0), vec(Sigma1)], sigma = sigma_randomEffects + V)
Rcpp::List implied_meta_var1_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  int i, j, s;

  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;

  std::string zeta = types["zeta"];
  std::string randomEffects = types["randomEffects"];

  // Add implied covariance structures:
  x = impliedcovstructures_cpp(x, "zeta", zeta, all);
  x = impliedcovstructures_cpp(x, "randomEffects", randomEffects, all);

  int g = 0; // meta_var1 is single-group

  bool proper = true;

  // General stuff from workspace:
  const Rcpp::List& extramats = ws.extramatrices;
  const arma::vec& nPerGroupVec = ws.nPerGroup;

  std::string est = Rcpp::as<std::string>(extramats["Vestimation"]);

  arma::mat V = Rcpp::as<arma::mat>(extramats["V"]);
  Rcpp::List Vall = extramats["Vall"];

  Rcpp::List grouplist = x[g];

  // ---- VAR(1) implied covariance ----
  // Model matrices:
  arma::mat beta = grouplist["beta"];
  arma::mat sigma_zeta = grouplist["sigma_zeta"];

  // Number of nodes:
  int nNode = beta.n_rows;

  // BetaStar = (I_{n^2} - beta kron beta)^{-1}
  arma::mat I2 = eye(nNode * nNode, nNode * nNode);
  arma::mat BetaStar = inv(I2 - kron(beta, beta));

  // Implied stationary distribution (vectorized):
  arma::vec sigmaZetaVec = vectorise(sigma_zeta);
  arma::vec vecSigma0 = BetaStar * sigmaZetaVec;

  // Reshape to matrix:
  arma::mat Sigma0(nNode, nNode);
  int curind = 0;
  for (j = 0; j < nNode; j++){
    for (i = 0; i < nNode; i++){
      Sigma0(i, j) = vecSigma0(curind);
      curind++;
    }
  }

  // Force symmetric:
  Sigma0 = 0.5 * (Sigma0 + Sigma0.t());

  // Implied lag-1 cross-covariance:
  arma::mat Sigma1 = beta * Sigma0;

  // Store intermediates for derivatives (when !all):
  if (!all){
    grouplist["BetaStar"] = BetaStar;
    grouplist["sigmaZetaVec"] = sigmaZetaVec;
    grouplist["IkronBeta"] = kronecker_I_X(beta, nNode);
    grouplist["Sigma0"] = Sigma0;
    grouplist["Sigma1"] = Sigma1;
  }

  // ---- Meta-analytic structure ----
  arma::mat sigma_randomEffects = grouplist["sigma_randomEffects"];

  // Meta-analytic mean = [vech(Sigma0), vec(Sigma1)]:
  arma::vec mu_vech = vech(Sigma0, true);
  arma::vec mu_vec = vectorise(Sigma1);
  arma::vec mu = join_cols(mu_vech, mu_vec);

  if (est == "averaged"){
    grouplist["mu"] = mu;

    // sigma = sigma_randomEffects + V:
    arma::mat sigma = sigma_randomEffects + V;
    grouplist["sigma"] = sigma;
    grouplist["kappa"] = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);

  } else {
    int nStudy = (int)nPerGroupVec(g);

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
Rcpp::List implied_meta_var1_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_meta_var1_cpp_core(x, model, all);
}
