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
Rcpp::List implied_lvm_cpp(
    const S4& model,
    bool all = false
){
  S4 sample = model.slot("sample");
  Rcpp::List means = sample.slot("means");
  
  
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);
  
  
  // Types:
  Rcpp::List types = model.slot("types");
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
    arma::mat nu = grouplist["nu"];
    arma::mat nu_eta = grouplist["nu_eta"];
    
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
    
    // Implied means:
    arma::vec mu = nu + lambda * BetaStar * nu_eta;
    grouplist["mu"] = mu;
    
    arma::mat sigma = Lambda_BetaStar * sigma_zeta * Lambda_BetaStar.t() + sigma_epsilon;
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

