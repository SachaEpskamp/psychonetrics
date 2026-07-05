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
Rcpp::List implied_dlvm1_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  int i, t, tt;

  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;

  std::string within_latent = types["within_latent"];
  std::string within_residual = types["within_residual"];
  std::string between_latent = types["between_latent"];
  std::string between_residual = types["between_residual"];

  // Add implied cov structure:
  x = impliedcovstructures_cpp(x, "zeta_within", within_latent, all);
  x = impliedcovstructures_cpp(x, "epsilon_within", within_residual, all);
  x = impliedcovstructures_cpp(x, "zeta_between", between_latent, all);
  x = impliedcovstructures_cpp(x, "epsilon_between", between_residual, all);

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
    
    // Model matrices:
    arma::mat nu = grouplist["nu"];
    arma::mat mu_eta = grouplist["mu_eta"];
    arma::mat lambda = grouplist["lambda"];
    arma::mat beta = grouplist["beta"];
    
    arma::mat sigma_zeta_within = grouplist["sigma_zeta_within"];
    arma::mat sigma_zeta_between = grouplist["sigma_zeta_between"];
    arma::mat sigma_epsilon_within = grouplist["sigma_epsilon_within"];
    arma::mat sigma_epsilon_between = grouplist["sigma_epsilon_between"];

    // Residual temporal effects (may be absent in models created before
    // beta_epsilon existed; those keep the white-noise residual structure):
    bool has_beta_eps = grouplist.containsElementNamed("beta_epsilon");
    arma::mat beta_epsilon;
    bool eps_ar = false;
    if (has_beta_eps){
      beta_epsilon = as<arma::mat>(grouplist["beta_epsilon"]);
      eps_ar = arma::any(arma::vectorise(beta_epsilon) != 0);
    }
    
    // Matrices I need in every model framework when estimating:
    // Beta star:
    int n_lat = beta.n_rows;
    arma::mat I2 = eye(n_lat*n_lat, n_lat*n_lat);
    
    // Some stuff needed now:
    arma::mat BetaStar = inv(I2 - kron(beta, beta));
    
    // Implied mean vector:
    arma::vec impMu = nu + lambda * mu_eta;
    
    arma::vec fullMu(nTime * nVar);
    for (i=0;i<nTime;i++){
      fullMu.submat(i*nVar,0,(i+1)*nVar-1,0) = impMu;
    }
    
    
    
    Rcpp::List allSigmas_within(nTime);
    arma::mat stationary_sigma_within_latent = matrixform(BetaStar * vectorise(sigma_zeta_within));
    
    allSigmas_within[0] = stationary_sigma_within_latent;
    
    // Fill implied:
    if (nTime > 1){
      for (t=1;t<nTime;t++){
        arma::mat lastsigma = allSigmas_within[t-1];
        allSigmas_within[t] = beta * lastsigma;
      }
    }
    
    
    
    // Create the block Toeplitz:
    arma::mat fullSigma_within_latent  = blockToeplitz_cpp(allSigmas_within);

    // Residual part: white noise (I_T kron sigma_epsilon_within), or, when a
    // residual temporal process is modeled (beta_epsilon non-zero), the block
    // Toeplitz of the stationary residual AR(1) process with innovation
    // covariance sigma_epsilon_within:
    Rcpp::List allSigmas_epsilon(nTime);
    arma::mat BetaStarEps;
    arma::mat fullSigma_epsilon;
    if (eps_ar){
      arma::mat I2eps = eye(nVar*nVar, nVar*nVar);
      BetaStarEps = inv(I2eps - kron(beta_epsilon, beta_epsilon));
      arma::mat stationary_sigma_epsilon = matrixform(BetaStarEps * vectorise(sigma_epsilon_within));
      allSigmas_epsilon[0] = stationary_sigma_epsilon;
      if (nTime > 1){
        for (t=1;t<nTime;t++){
          arma::mat lastsigma_eps = allSigmas_epsilon[t-1];
          allSigmas_epsilon[t] = beta_epsilon * lastsigma_eps;
        }
      }
      fullSigma_epsilon = blockToeplitz_cpp(allSigmas_epsilon);
    } else {
      fullSigma_epsilon = arma::mat(kronecker_I_X(sigma_epsilon_within, nTime));
    }

    // Full within-subject cov matrix:
    arma::mat fullSigma_within = kronecker_I_X(lambda, nTime) * fullSigma_within_latent *  kronecker_I_X(lambda.t(), nTime) +  fullSigma_epsilon;
    
    // Full between-subject cov matrix:
    arma::mat fullSigma_between(nTime * nVar, nTime * nVar);
    arma::mat subSigma_between = lambda * sigma_zeta_between * lambda.t() + sigma_epsilon_between;
    
    
    for (t=0;t<nTime;t++){
      for (tt=0;tt<nTime;tt++){
        
      fullSigma_between.submat(t*nVar,tt*nVar,(t+1)*nVar-1,(tt+1)*nVar-1) = subSigma_between;
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
    
    // FIXME: forcing symmetric, but not sure why this is needed...
    // grouplist["kappa = 0.5*(grouplist["kappa + t(grouplist["kappa))
    
    // Let's round to make sparse if possible:
    // grouplist["kappa = as(round(grouplist["kappa,14),"Matrix")
    
    
    // Extra matrices needed in optimization:
    if (!all){
      grouplist["BetaStar"] = BetaStar;
      // grouplist["E = Emat(nrow(grouplist["beta),grouplist["beta)
      grouplist["allSigmas_within"] = allSigmas_within;
      grouplist["IkronBeta"] = kronecker_I_X(beta, n_lat);
      grouplist["lamWkronlamW"] = kron(lambda, lambda);
      if (eps_ar){
        grouplist["BetaStarEps"] = BetaStarEps;
        grouplist["allSigmas_epsilon"] = allSigmas_epsilon;
        grouplist["IkronBetaEps"] = kronecker_I_X(beta_epsilon, nVar);
      }
    } else {
      // Stationary within-person residual covariance (equals the innovation
      // covariance sigma_epsilon_within when beta_epsilon is zero/absent):
      arma::mat sigma_epsilon_stationary;
      if (eps_ar){
        sigma_epsilon_stationary = as<arma::mat>(allSigmas_epsilon[0]);
      } else {
        sigma_epsilon_stationary = sigma_epsilon_within;
      }

      arma::mat sigma_within = lambda * stationary_sigma_within_latent * lambda.t() + sigma_epsilon_stationary;
      
      
      
      grouplist["sigma_within"] = sigma_within;

      grouplist["sigma_between"] = subSigma_between;
      grouplist["sigma_within_full"] = fullSigma_within;
      grouplist["sigma_eta_within"] = stationary_sigma_within_latent;
      grouplist["sigma_eta_within_lag1"] = allSigmas_within[1];
      grouplist["sigma_crosssection"] = sigma_within + subSigma_between;
      grouplist["sigma_epsilon_within_stationary"] = sigma_epsilon_stationary;
      
    }
    
    // Add PDC:
    // FIXME: This should not be needed?
    // if (!is.null(grouplist["kappa_zeta_within)){
    arma::mat kappa_zeta_within;
    
    if (grouplist.containsElementNamed("kappa_zeta_within")){
      kappa_zeta_within = as<arma::mat>(grouplist["kappa_zeta_within"]);
    } else {
      kappa_zeta_within = solve_symmetric_cpp_matrixonly_withcheck(sigma_zeta_within, proper);
    }
    
    grouplist["PDC"] = computePDC_cpp(grouplist["beta"], kappa_zeta_within, grouplist["sigma_zeta_within"]);
    // } else {
    //   grouplist["PDC = t(grouplist["beta)
    //   grouplist["PDC[] = 0
    // }
    
    // Return properness:
    grouplist["proper"] = proper;
    
    
    x[g] = grouplist;
  }
  
  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_dlvm1_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_dlvm1_cpp_core(x, model, all);
}

