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

// Core implementation that takes pre-formed model matrices
Rcpp::List implied_meta_varcov_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  int s;

  S4 sample = model.slot("sample");
  Rcpp::List means = sample.slot("means");
  bool corinput = sample.slot("corinput");
  Rcpp::List groupslist = sample.slot("groups");
  NumericVector nPerGroup = groupslist["nobs"];

  // Types:
  Rcpp::List types = model.slot("types");

  std::string y = types["y"];
  std::string randomEffects = types["randomEffects"];

  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "y", y, all);
  x = impliedcovstructures_cpp(x, "randomEffects", randomEffects, all);

  // int nGroup = x.length();
  int g = 0;

  bool proper = true;

  // General stuff:
  Rcpp::List extramats = model.slot("extramatrices");

  std::string est = extramats["Vestimation"];

  arma::mat V = extramats["V"];
  Rcpp::List Vall = extramats["Vall"];

  Rcpp::List grouplist = x[g];

  // Model matrices:
  arma::mat sigma_y = grouplist["sigma_y"];
  arma::mat sigma_randomEffects = grouplist["sigma_randomEffects"];

  if (est == "averaged"){
    arma::vec mu;

    if (corinput){
      mu = vech(cov2cor_cpp(sigma_y), false);
      grouplist["mu"] = mu;
    } else {
      mu = vech(sigma_y, true);
      grouplist["mu"] = mu;
    }

    arma::mat sigma = sigma_randomEffects + V;
    grouplist["sigma"] = sigma;
    grouplist["kappa"] = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);

  } else {
    int nStudy = nPerGroup[g];

    arma::mat mu;
    if (corinput){
      mu = vech(cov2cor_cpp(sigma_y), false);
    } else {
      mu = vech(sigma_y, true);
    }

    Rcpp::List mulist(nStudy);
    Rcpp::List sigmalist(nStudy);
    Rcpp::List kappalist(nStudy);

    for (s=0;s<nStudy;s++){
      mulist[s] = mu;

      arma::mat Vcur = Vall[s];
      arma::mat cursigma = sigma_randomEffects + Vcur;;
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
Rcpp::List implied_meta_varcov_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_meta_varcov_cpp_core(x, model, all);
}

