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
Rcpp::List implied_lvm_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;
  std::string latent = types["latent"];
  std::string residual = types["residual"];

  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "zeta", latent, all);
  x = impliedcovstructures_cpp(x, "epsilon", residual, all);

  int nGroup = x.length();
  int g;

  for (g=0; g<nGroup; g++){
    bool proper = true;

    Rcpp::List grouplist = x[g];

    // Model matrices:
    arma::mat beta = grouplist["beta"];
    arma::mat lambda = grouplist["lambda"];
    arma::mat sigma_zeta = grouplist["sigma_zeta"];
    arma::mat sigma_epsilon = grouplist["sigma_epsilon"];

    // Matrices I need in every model framework when estimating:
    int n = beta.n_rows;
    arma::mat I = eye(n,n);
    arma::mat BetaStar = inv(I - beta);
    arma::mat Lambda_BetaStar = lambda * BetaStar;

    // If NOT all (I know it is confusing), store these matrices (and some more):
    if (!all){
      grouplist["BetaStar"] = BetaStar;
      grouplist["Lambda_BetaStar"] = Lambda_BetaStar;

      arma::mat Betasta_sigmaZeta = BetaStar * sigma_zeta;
      arma::mat tBetakronBeta = kron(BetaStar.t(), BetaStar);

      grouplist["Betasta_sigmaZeta"] = Betasta_sigmaZeta;
      grouplist["tBetakronBeta"] = tBetakronBeta;

    }

    // Implied means (only if mean structure is modeled):
    if (grouplist.containsElementNamed("nu")){
      arma::mat nu = grouplist["nu"];
      arma::mat nu_eta = grouplist["nu_eta"];
      arma::vec mu = nu + lambda * BetaStar * nu_eta;
      grouplist["mu"] = mu;
    }

    // Factor part:
    arma::mat factorPart = Lambda_BetaStar * sigma_zeta * Lambda_BetaStar.t();

    // When corinput=TRUE, derive diagonal of sigma_epsilon from diag(sigma)=1 constraint:
    bool corinput_flag = ws.corinput;
    if (grouplist.containsElementNamed("corinput")){
      corinput_flag = Rcpp::as<bool>(grouplist["corinput"]);
    }
    if (corinput_flag){
      arma::vec derivedDiag = 1.0 - diagvec(factorPart);
      sigma_epsilon.diag() = derivedDiag;
      grouplist["sigma_epsilon"] = sigma_epsilon;
    }

    arma::mat sigma = factorPart + sigma_epsilon;
    grouplist["sigma"] = sigma;

    // FIXME: Make symmetric?

    arma::mat kappa = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
    grouplist["kappa"] = kappa;

    // Return properness:
    grouplist["proper"] = proper;

    x[g] = grouplist;
  }

  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_lvm_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_lvm_cpp_core(x, model, all);
}

