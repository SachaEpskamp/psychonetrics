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


// Inner functions to use:
// [[Rcpp::export]]
arma::mat d_mu_lambda_tsdlvm1_cpp(
    const arma::mat& mu_eta,
    const arma::sp_mat& I_y){
  
  arma::mat res = (arma::mat)kronecker_X_I(mu_eta.t(), I_y.n_rows);
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigmak_lambda_tsdlvm1_cpp(
    const arma::mat& lambda,
    int k, 
    const arma::mat& Sigma_eta_0, 
    const arma::mat& Sigma_eta_1, 
    const arma::sp_mat& C_y_eta, 
    const arma::sp_mat& I_y, 
    const arma::sp_mat& L_y){
  
  arma::mat sigEta;
  if (k == 0){
    sigEta = Sigma_eta_0;
  } else {
    sigEta= Sigma_eta_1;
  }
  
  sp_mat within = (
    kronecker_I_X(lambda * sigEta,I_y.n_rows) * C_y_eta + 
      kronecker_X_I( lambda * sigEta.t(), I_y.n_rows)  
  );
  
  arma::mat res = (arma::mat)(within);
  
  if (k == 0){
    return(L_y * res);
  } else {
    return(res);
  }
}

// [[Rcpp::export]]
arma::mat d_sigma0_sigma_zeta_tsdlvm1_cpp(
    const arma::mat& BetaStar, 
    const arma::sp_mat& D_eta){
  arma::mat res = BetaStar * D_eta;
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigma0_beta_tsdlvm1_cpp(
    const arma::mat& BetaStar,
    const arma::sp_mat& I_eta,
    const arma::mat& Sigma_eta_1,
    const arma::sp_mat& C_eta_eta){
  
  int n = I_eta.n_rows;
  int nn = n * n;
  
  arma::mat res = (eye(nn, nn) + C_eta_eta) * BetaStar * kronecker_X_I(Sigma_eta_1, I_eta.n_rows); 
  
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigma1_beta_tsdlvm1_cpp(
  const arma::mat& J_sigma_beta,
  const arma::sp_mat& IkronBeta,
  const arma::mat& Sigma_eta_0,
  const arma::sp_mat& I_eta){
  
  arma::mat res = IkronBeta * J_sigma_beta + kronecker_X_I(Sigma_eta_0.t(), I_eta.n_rows);
  return(res);
}



// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_tsdlvm1_group_cpp(
    const Rcpp::List& grouplist
){
  // Stuff I need now:
  std::string zeta = grouplist["zeta"];
  std::string epsilon = grouplist["epsilon"];
  
  arma::mat lambda = grouplist["lambda"];
  arma::sp_mat P = grouplist["P"];
  
  arma::sp_mat L_eta = grouplist["L_eta"];
  arma::sp_mat L_y = grouplist["L_y"];
  
  
  // Number of variables:
  int  nvar = lambda.n_rows * 2;
  
  // Number of nodes:
  int nNode = lambda.n_rows;
  
  // Number of latents:
  int nLat = lambda.n_cols;
  
  // Number of observations:
  int nobs = nvar + // Means
    (nvar * (nvar+1))/2; // Variances
  
  // total number of elements:
  int nelement = nNode + // Exogenous means
    nNode + // intercepts
    nLat + // Latent means
    nNode*(nNode+1)/2 + // Exogenous variances
    nNode * nLat + // Factor loadings
    nLat * (nLat + 1) / 2 + // Contemporaneous
    nLat * nLat + // Beta
    nNode * (nNode+1) / 2; // Residual
  
  
  arma::mat Jac = zeros(nobs,nelement);
  
  // Indices:
  arma::vec meanInds_exo = seq_len_inds(0,nNode);
  arma::vec meanInds_endo = seq_len_inds(meanInds_exo(1) + 1, nNode);
  arma::vec sigmaStarInds =seq_len_inds(meanInds_endo(1) + 1, nNode*(nNode+1)/2); 
  arma::vec sigma0Inds =  seq_len_inds(sigmaStarInds(1) + 1, nNode*(nNode+1)/2);
  arma::vec sigma1Inds = seq_len_inds(sigma0Inds(1) + 1, nNode * nNode);
  
  // Indices model:
  arma::vec exomeanInds = seq_len_inds(0,nNode);
  arma::vec interceptInds = seq_len_inds(exomeanInds(1) + 1, nNode);
  arma::vec latmeanInds = seq_len_inds(interceptInds(1) + 1, nLat); 
  arma::vec exovarInds =  seq_len_inds(latmeanInds(1) + 1, nNode*(nNode+1)/2);
  arma::vec lambdaInds =  seq_len_inds(exovarInds(1) + 1, nNode * nLat);
  arma::vec contInds =  seq_len_inds(lambdaInds(1) + 1, nLat * (nLat + 1) / 2);
  arma::vec betaInds =  seq_len_inds(contInds(1) + 1, nLat * nLat);
  arma::vec residInds = seq_len_inds(betaInds(1) + 1, nNode*(nNode+1)/2);
  
  
  //
  
  
  ////// Augmentation parts //////
  arma::mat aug_zeta;
  arma::mat aug_epsilon;
  
  // Contemporaneous
  if (zeta == "chol"){
    aug_zeta =   d_sigma_cholesky_cpp(grouplist["lowertri_zeta"],L_eta,grouplist["C_eta_eta"],grouplist["I_eta"]);
  } else if (zeta == "prec"){
    aug_zeta = d_sigma_kappa_cpp(L_eta, grouplist["D_eta"], grouplist["sigma_zeta"]);
  } else if (zeta == "ggm"){
    aug_zeta = join_rows(
      d_sigma_omega_cpp(L_eta, grouplist["delta_IminOinv_zeta"], grouplist["A_eta"], grouplist["delta_zeta"], grouplist["Dstar_eta"]),
      d_sigma_delta_cpp(L_eta, grouplist["delta_IminOinv_zeta"], grouplist["I_eta"], grouplist["A_eta"])
    );
  } else if (zeta == "cov"){
    aug_zeta = eye(nLat*(nLat+1)/2, nLat*(nLat+1)/2);
  }
  
  // Residual
  if (epsilon == "chol"){
    aug_epsilon =   d_sigma_cholesky_cpp(grouplist["lowertri_epsilon"],L_y,grouplist["C_y_y"],grouplist["I_y"]);
  } else if (epsilon == "prec"){
    aug_epsilon = d_sigma_kappa_cpp(L_y, grouplist["D_y"], grouplist["sigma_epsilon"]);
  } else if (epsilon == "ggm"){
    aug_epsilon = join_rows(
      d_sigma_omega_cpp(L_y, grouplist["delta_IminOinv_epsilon"], grouplist["A_y"], grouplist["delta_epsilon"], grouplist["Dstar_y"]),
      d_sigma_delta_cpp(L_y, grouplist["delta_IminOinv_epsilon"], grouplist["I_y"], grouplist["A_y"])
    );
  } else if (epsilon == "cov"){
    aug_epsilon = eye(nNode*(nNode+1)/2, nNode*(nNode+1)/2);
  }

  // Exogenous mean part:
  Jac.submat(meanInds_exo(0),exomeanInds(0),meanInds_exo(1),exomeanInds(1)) =  eye(nNode,nNode);


  // Intercept part:
  Jac.submat(meanInds_endo(0),interceptInds(0),meanInds_endo(1),interceptInds(1)) =  eye(nNode,nNode);

  // Fill latent mean part:
  Jac.submat(meanInds_endo(0),latmeanInds(0),meanInds_endo(1),latmeanInds(1)) =  lambda;


  // Fill mean to factor loading part:
  Jac.submat(meanInds_endo(0),lambdaInds(0),meanInds_endo(1),lambdaInds(1)) =   d_mu_lambda_tsdlvm1_cpp(grouplist["mu_eta"], grouplist["I_y"]);

  // Exogenous block variances:
  Jac.submat(sigmaStarInds(0),exovarInds(0),sigmaStarInds(1),exovarInds(1)) =  d_sigma_cholesky_cpp(grouplist["exo_cholesky"],L_y,grouplist["C_y_y"],grouplist["I_y"]);
  
    ////// Sigma 0 part:

    // Sigma 0 to lambda:
    Jac.submat(sigma0Inds(0),lambdaInds(0),sigma0Inds(1),lambdaInds(1)) =  d_sigmak_lambda_tsdlvm1_cpp(grouplist["lambda"], 0, grouplist["Sigma_eta_0"], grouplist["Sigma_eta_1"], grouplist["C_y_eta"], grouplist["I_y"], grouplist["L_y"]);

    // Sigma 0 to contemporaneous
    arma::mat J_sigma_zeta = d_sigma0_sigma_zeta_tsdlvm1_cpp(grouplist["BetaStar"], grouplist["D_eta"]) * aug_zeta;
    arma::mat lamWkronlamW = grouplist["lamWkronlamW"];
    
    Jac.submat(sigma0Inds(0),contInds(0),sigma0Inds(1),contInds(1)) =  L_y * lamWkronlamW * J_sigma_zeta;
    
      
    // Fill s0 to beta part (and store for later use):
    arma::mat J_sigma_beta = d_sigma0_beta_tsdlvm1_cpp( grouplist["BetaStar"],  grouplist["I_eta"],  grouplist["Sigma_eta_1"],  grouplist["C_eta_eta"]);
    Jac.submat(sigma0Inds(0),betaInds(0),sigma0Inds(1),betaInds(1)) =  L_y  * lamWkronlamW * J_sigma_beta;
    

    // Fill s0 to sigma_epsilon_within:
    Jac.submat(sigma0Inds(0),residInds(0),sigma0Inds(1),residInds(1)) = aug_epsilon;
      

    //// Sigma 1 parts

    // Fill s1 to lambda part:
    Jac.submat(sigma1Inds(0),lambdaInds(0),sigma1Inds(1),lambdaInds(1)) = d_sigmak_lambda_tsdlvm1_cpp(grouplist["lambda"], 1, grouplist["Sigma_eta_0"], grouplist["Sigma_eta_1"], grouplist["C_y_eta"], grouplist["I_y"], grouplist["L_y"]);;
    
    

    // Fill s1 to sigma_zeta part:
    arma::sp_mat IkronBeta = grouplist["IkronBeta"];
    J_sigma_zeta = IkronBeta * J_sigma_zeta;
    Jac.submat(sigma1Inds(0),contInds(0),sigma1Inds(1),contInds(1)) = lamWkronlamW * J_sigma_zeta;
    
    // Fill s1 to beta part (and store for later use):
    J_sigma_beta =  d_sigma1_beta_tsdlvm1_cpp(J_sigma_beta, grouplist["IkronBeta"], grouplist["Sigma_eta_0"], grouplist["I_eta"]);
    Jac.submat(sigma1Inds(0),betaInds(0),sigma1Inds(1),betaInds(1)) = lamWkronlamW * J_sigma_beta;
    
    
    // Permute:
    Jac = P * Jac;
  

  // Return:
  return(Jac);
}



// 
// [[Rcpp::export]]
arma::mat d_phi_theta_tsdlvm1_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_tsdlvm1_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}




