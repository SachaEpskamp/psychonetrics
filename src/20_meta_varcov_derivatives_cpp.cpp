// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "14_varcov_derivatives_cpp.h"
#include "02_algebrahelpers_RcppHelpers.h"



// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Group grad:
// [[Rcpp::export]]
arma::mat d_phi_theta_meta_varcov_group_cpp(
    const Rcpp::List& grouplist
){
  // objects needed now:
  arma::mat sigma_y = grouplist["sigma_y"];
  std::string y = grouplist["y"];
  std::string randomEffects = grouplist["randomEffects"];
  bool metacor = grouplist["metacor"];
  
  
  // Number of variables:
  int nvar = sigma_y.n_rows;
  
  // Number of correlations:
  int ncor = nvar * (nvar-1) / 2;
  
  // Number of modeled elements:
  int nmod = ncor + (metacor == false) * nvar;
  
  // Number of observations:
  int nobs = nmod + // correlations
    nmod * (nmod+1) / 2; // Random effects
  
  // Mean part:
  int meanPart_start = 0;
  int meanPart_end = nmod - 1;
  
  // Variance part:
  int varPart_start = nmod;
  int varPart_end = nmod + nmod*(nmod+1)/2 - 1;
  
  // Empty Jacobian:
  arma::mat Jac = zeros(nobs,nobs);
  
  // Fill the mean part with model for cors/covs:
  arma::sp_mat Lmat;
  if (metacor){
    arma::sp_mat Lstar = grouplist["Lstar"];
    Lmat = Lstar;
  } else {
    arma::sp_mat L = grouplist["L"];
    Lmat = L;
  }
  
  // Derivatives for y:
  // Fill the simga part:
  if (y == "cov"){
    Jac.submat(meanPart_start,meanPart_start,meanPart_end,meanPart_end ) =  eye(nmod,nmod);
    
    
  } else if (y == "chol"){
    
    Jac.submat(meanPart_start,meanPart_start,meanPart_end,meanPart_end ) =  d_sigma_cholesky_cpp( 
      grouplist["lowertri_y"], Lmat, grouplist["C"], grouplist["In"]
    );
    
    
  } else if (y == "ggm"){
    int netPart_start = 0;
    int netPart_end = nvar * (nvar - 1) / 2 - 1;
    int scalingPart_start = netPart_end + 1;
    int scalingPart_end = scalingPart_start + nvar - 1;
    
    
    if (metacor){
      
      Jac.submat(meanPart_start,netPart_start,meanPart_end,netPart_end) =  d_sigma_omega_corinput_cpp(
        Lmat, grouplist["delta_IminOinv_y"], grouplist["A"], grouplist["delta_y"], grouplist["Dstar"], grouplist["IminOinv_y"], grouplist["In"]
      );
      
      
    } else {
      
      Jac.submat(meanPart_start,netPart_start,meanPart_end,netPart_end ) =  d_sigma_omega_cpp( 
        Lmat, grouplist["delta_IminOinv_y"], grouplist["A"], grouplist["delta_y"], grouplist["Dstar"]
      );
      
      Jac.submat(meanPart_start,scalingPart_start,meanPart_end,scalingPart_end ) =  d_sigma_delta_cpp( 
        Lmat, grouplist["delta_IminOinv_y"], grouplist["In"],  grouplist["A"]
      );
      
      
    }
    
    
  } else if (y == "prec"){
    
    
    Jac.submat(meanPart_start,meanPart_start,meanPart_end,meanPart_end ) =  d_sigma_kappa_cpp( 
      Lmat, grouplist["D"], grouplist["sigma_y"]
    );
    
  } else if (y == "cor"){
    
    int corPart_start = 0;
    int corPart_end = nvar * (nvar-1) / 2 - 1;
    int sdPart_start = corPart_end + 1;
    int sdPart_end = sdPart_start + nvar - 1;
    
    Jac.submat(meanPart_start,corPart_start,meanPart_end,corPart_end ) =  d_sigma_rho_cpp( 
      Lmat, grouplist["SD_y"], grouplist["A"], grouplist["Dstar"]
    );
    
    if (metacor == false){
      
      Jac.submat(meanPart_start,sdPart_start,meanPart_end,sdPart_end ) =  d_sigma_SD_cpp( 
        Lmat, grouplist["SD_IplusRho_y"], grouplist["In"], grouplist["A"]
      );
      
    }
    
  }
  
  
  
  
  // Derivatives for random effects:
  int nEl = nmod * (nmod+1) / 2;
  // Fill the simga part:
  if (randomEffects == "cov"){
    Jac.submat(varPart_start,varPart_start,varPart_end,varPart_end ) =  eye(nEl,nEl);
    
    
  } else if (randomEffects == "chol"){
    
    Jac.submat(varPart_start,varPart_start,varPart_end,varPart_end ) =  d_sigma_cholesky_cpp( 
      grouplist["lowertri_randomEffects"], grouplist["L_c"], grouplist["C_c"], grouplist["In_c"]
    );
    
    
  } else if (randomEffects == "ggm"){
    int netPart_start = nmod;
    int netPart_end = nmod + nmod * (nmod-1) / 2 - 1;
    int scalingPart_start = netPart_end + 1;
    int scalingPart_end = scalingPart_start + nmod - 1;
    
    
    // if (metacor){
    //   
    //   Jac.submat(varPart_start,netPart_start,varPart_end,netPart_end ) =  d_sigma_omega_corinput_cpp(
    //     grouplist["L_c"], grouplist["delta_IminOinv_randomEffects"], grouplist["A_c"], grouplist["delta_randomEffects"], grouplist["Dstar_c"], grouplist["IminOinv_randomEffects"], grouplist["In_c"]
    //   );
    //   
    //   
    // } else {
      
      Jac.submat(varPart_start,netPart_start,varPart_end,netPart_end ) =  d_sigma_omega_cpp( 
        grouplist["L_c"], grouplist["delta_IminOinv_randomEffects"], grouplist["A_c"], grouplist["delta_randomEffects"], grouplist["Dstar_c"]
      );
      
      Jac.submat(varPart_start,scalingPart_start,varPart_end,scalingPart_end ) =  d_sigma_delta_cpp( 
        grouplist["L_c"], grouplist["delta_IminOinv_randomEffects"], grouplist["In_c"],  grouplist["A_c"]
      );
      
      
    // }
    
    
  } else if (randomEffects == "prec"){
    
    
    Jac.submat(varPart_start,varPart_start,varPart_end,varPart_end ) =  d_sigma_kappa_cpp( 
      grouplist["L_c"], grouplist["D_c"], grouplist["sigma_randomEffects"]
    );
    
  } else if (randomEffects == "cor"){
    
    int corPart_start = nmod;
    int corPart_end = nmod + nmod * (nmod-1) / 2 - 1;
    int sdPart_start = corPart_end + 1;
    int sdPart_end = sdPart_start + nmod - 1;
    
  
    Jac.submat(varPart_start,corPart_start,varPart_end,corPart_end ) =  d_sigma_rho_cpp( 
      grouplist["L_c"], grouplist["SD_randomEffects"], grouplist["A_c"], grouplist["Dstar_c"]
    );
    
    // if (metacor == false){
      
      Jac.submat(varPart_start,sdPart_start,varPart_end,sdPart_end ) =  d_sigma_SD_cpp( 
        grouplist["L_c"], grouplist["SD_IplusRho_randomEffects"], grouplist["In_c"], grouplist["A_c"]
      );
      
    // }
    
  }
  
  
  // Return:
  return(Jac);
}

// Full grad:
// [[Rcpp::export]]
arma::mat d_phi_theta_meta_varcov_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_meta_varcov_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}
