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
Rcpp::List implied_ml_lvm_cpp_core(
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
  arma::mat design = extramats["designPattern"];
  
  int nVar = design.n_rows;
  int nMaxInCluster = design.n_cols;

  // Indices (column-major over the design matrix, i.e. measurement t outer,
  // variable i inner) of the observed variables. With an incomplete design
  // (some variables missing at some measurements) the implied mu/sigma must
  // be subset to the observed positions, exactly as in the R twin
  // implied_ml_lvm (R/22_ml_lvm_implied.R):
  arma::uvec obsInds(nVar * nMaxInCluster);
  int nObs = 0;
  for (int tt2 = 0; tt2 < nMaxInCluster; tt2++){
    for (int i2 = 0; i2 < nVar; i2++){
      if (design(i2, tt2) == 1){
        obsInds(nObs) = tt2 * nVar + i2;
        nObs++;
      }
    }
  }
  obsInds = obsInds.head(nObs);
  bool completeDesign = (nObs == nVar * nMaxInCluster);


  for (g=0; g<nGroup; g++){
    bool proper = true;
    
    Rcpp::List grouplist = x[g];
    
    // Model matrices:
    arma::mat nu = grouplist["nu"];
    arma::mat nu_eta = grouplist["nu_eta"];
    arma::mat lambda = grouplist["lambda"];
    arma::mat beta_within = grouplist["beta_within"];
    arma::mat beta_between = grouplist["beta_between"];
    arma::mat sigma_zeta_within = grouplist["sigma_zeta_within"];
    arma::mat sigma_zeta_between = grouplist["sigma_zeta_between"];
    arma::mat sigma_epsilon_within = grouplist["sigma_epsilon_within"];
    arma::mat sigma_epsilon_between = grouplist["sigma_epsilon_between"];
    
    // Matrices I need in every model framework when estimating:
    // Beta star:
    int n_lat = lambda.n_cols;
    arma::mat I = eye(n_lat, n_lat);
    // arma::mat I2 = eye(n_lat*n_lat, n_lat*n_lat);
    
    // Some stuff needed now:
    // Beta star within/between. When beta == 0 (the default) the inverse is
    // exactly I, so it (and the Kronecker products below) can be skipped
    // (bit-identical shortcut):
    const bool beta_within_zero = beta_within.is_zero();
    const bool beta_between_zero = beta_between.is_zero();
    arma::mat BetaStar_within;
    if (beta_within_zero){
      BetaStar_within = I;
    } else {
      BetaStar_within = inv(I - beta_within);
    }
    arma::mat BetaStar_between;
    if (beta_between_zero){
      BetaStar_between = I;
    } else {
      BetaStar_between = inv(I - beta_between);
    }
    
    arma::mat Betasta_sigmaZeta_within = BetaStar_within * sigma_zeta_within;
    arma::mat Betasta_sigmaZeta_between = BetaStar_between * sigma_zeta_between;
    
    // Implied mean vector:
    arma::vec impMu = nu + lambda * BetaStar_between * nu_eta;
    
    arma::vec fullMu(nMaxInCluster * nVar);
    for (i=0;i<nMaxInCluster;i++){
      fullMu.submat(i*nVar,0,(i+1)*nVar-1,0) = impMu;
    }
    
    
    // Implied within covariance:
    arma::mat sigma_eta_within = Betasta_sigmaZeta_within * BetaStar_within.t();
    arma::mat sigma_eta_between = Betasta_sigmaZeta_between * BetaStar_between.t();
    
    // Implied structures:
    arma::mat sigma_within = lambda * sigma_eta_within * lambda.t() + sigma_epsilon_within;
    arma::mat sigma_between = lambda * sigma_eta_between * lambda.t() + sigma_epsilon_between;
    
    

    
    // FIXME: fiew things here I do sparse. Smart or not?
    // Create the block Toeplitz:
    arma::mat fullSigma_within  = (arma::mat)kronecker_I_X(sigma_within, nMaxInCluster);
    
    // Full between-subject cov matrix:
    arma::mat fullSigma_between(nMaxInCluster * nVar, nMaxInCluster * nVar);

    for (t=0;t<nMaxInCluster;t++){
      for (tt=0;tt<nMaxInCluster;tt++){
        
        fullSigma_between.submat(t*nVar,tt*nVar,(t+1)*nVar-1,(tt+1)*nVar-1) = sigma_between;
      }
    }
    

    
    // Full implied covmat:
    arma::mat fullSigma = fullSigma_within + fullSigma_between;

    // Subset to the observed positions when the design is incomplete
    // (matches the R twin, which subsets by as.vector(design) == 1):
    if (!completeDesign){
      arma::vec subMu = fullMu(obsInds);
      arma::mat subSigma = fullSigma.submat(obsInds, obsInds);
      // Force symmetry as in the R twin:
      subSigma = 0.5 * (subSigma + subSigma.t());
      fullMu = subMu;
      fullSigma = subSigma;
    }

    // Add to list:
    grouplist["mu"] = fullMu;
    grouplist["sigma"] = fullSigma;



    // Precision:
    grouplist["kappa"] =  solve_symmetric_cpp_matrixonly_withcheck(fullSigma, proper);
    
    // FIXME: forcing symmetric, but not sure why this is needed...
    // grouplist["kappa = 0.5*(grouplist["kappa + t(grouplist["kappa))
    
    // Let's round to make sparse if possible:
    // grouplist["kappa = as(round(grouplist["kappa,14),"Matrix")
    
    
    // Extra matrices needed in optimization:
    if (!all){
      arma::mat Lambda_BetaStar_within = lambda *  BetaStar_within;
      arma::mat Lambda_BetaStar_between = lambda *  BetaStar_between;
      arma::mat tBetakronBeta_within;
      if (beta_within_zero){
        tBetakronBeta_within = eye(n_lat*n_lat, n_lat*n_lat);
      } else {
        tBetakronBeta_within = kron( BetaStar_within.t() , BetaStar_within );
      }
      arma::mat tBetakronBeta_between;
      if (beta_between_zero){
        tBetakronBeta_between = eye(n_lat*n_lat, n_lat*n_lat);
      } else {
        tBetakronBeta_between = kron( BetaStar_between.t()  , BetaStar_between  );
      }
      
      grouplist["BetaStar_within"] = BetaStar_within;
      grouplist["BetaStar_between"] = BetaStar_between;
      grouplist["Betasta_sigmaZeta_within"] = Betasta_sigmaZeta_within;
      grouplist["Betasta_sigmaZeta_between"] = Betasta_sigmaZeta_between;
      grouplist["Lambda_BetaStar_within"] = Lambda_BetaStar_within;
      grouplist["Lambda_BetaStar_between"] = Lambda_BetaStar_between;
      grouplist["tBetakronBeta_within"] = tBetakronBeta_within;
      grouplist["tBetakronBeta_between"] = tBetakronBeta_between;
      
      grouplist["sigma_eta_within"] = sigma_eta_within;
        grouplist["sigma_eta_between"] = sigma_eta_between;
        grouplist["lamWkronlamW"] = kron(lambda, lambda);
      
      
    } else {
    
    
    grouplist["sigma_within"] = sigma_within;
      grouplist["sigma_between"] = sigma_between;
      grouplist["sigma_eta_within"] = sigma_eta_within;
      grouplist["sigma_eta_between"] = sigma_eta_between;
      grouplist["sigma_eta"] = sigma_eta_within + sigma_eta_between;
      grouplist["sigma_crosssection"] = sigma_within + sigma_between;
      
      
    }
    
    
    // Return properness:
    grouplist["proper"] = proper;
    
    
    x[g] = grouplist;
  }
  
  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_ml_lvm_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_ml_lvm_cpp_core(x, model, all);
}

