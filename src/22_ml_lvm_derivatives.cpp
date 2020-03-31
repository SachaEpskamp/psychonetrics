// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "14_varcov_derivatives_cpp.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "15_lvm_derivatives_cpp.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_ml_lvm_group_cpp(
    const Rcpp::List& grouplist
){
  // Setup:
  // int t;
  
  
  // Things I need now:
  std::string within_latent  = grouplist["within_latent"];
  std::string within_residual  = grouplist["within_residual"];
  std::string between_latent  = grouplist["between_latent"];
  std::string between_residual  = grouplist["between_residual"];
  
  arma::mat lambda  = grouplist["lambda"];
  
  
  arma::mat design = grouplist["designPattern"];
  arma::sp_mat P = grouplist["P"];
  arma::sp_mat L_eta = grouplist["L_eta"];
  arma::sp_mat L_y = grouplist["L_y"];
  arma::sp_mat D_y = grouplist["D_y"];
  
  
  // Number of variables:
  int nVar = design.n_rows;
  
  // Number of times points:
  // int nMaxInCluster = design.n_cols;
  
  // Number of latents:
  int nLat = lambda.n_cols;
  
  // Number of within latents:
  // nLat <- ncol(lambda)
  
  // Number of between latents:
  // nLat <- ncol(lambda)
  
  // I need to construct the Jacobian for the following "observations:"
  int nobs = nVar + // Means
    (nVar * (nVar+1))/2 + // Variances
    (nVar * nVar); // covariances between members of cluster
  
  // total number of elements:
  int nelement = nVar + // intercepts in nu
    nLat + // Latent means
    nVar * nLat + // Factor loadings
    nLat * nLat + // Beta within
    nLat * (nLat+1) / 2 + // within factor variances
    nVar * (nVar + 1) / 2 + // Within residuals
    // nVar * nLat + // Between-subject factor loadings
    nLat * nLat + // Beta between
    nLat * (nLat + 1) / 2 + // Between-subject factor variances
    nVar * (nVar + 1) / 2; // between residuals
  
  // Form empty matrix:
  arma::mat Jac = zeros(nobs,nelement);
  
  
  
  // Indices model:
  // arma::vec interceptInds = seq_len_inds(0, nvar);
  //
  // arma::vec nuetaInds = seq_len_inds(interceptInds(1) + 1, nlat);
  //
  // Indices:
  arma::vec meanInds = seq_len_inds(0, nVar);
  
  arma::vec sigInds_inPerson = seq_len_inds(meanInds(1) + 1, nVar*(nVar+1)/2);
  arma::vec sigInds_betPersons = seq_len_inds(sigInds_inPerson(1) + 1, nVar*nVar);
  
  
  // Indices model:
  arma::vec nu_inds = seq_len_inds(0, nVar);
  arma::vec nu_eta_inds = seq_len_inds(nu_inds(1) + 1, nLat);
  arma::vec lambda_inds = seq_len_inds(nu_eta_inds(1) + 1, nVar * nLat);
  
  
  
  // lambda_inds <- max(nu_eta_inds) + seq_len(nVar * nLat)
  arma::vec beta_within_inds = seq_len_inds(lambda_inds(1) + 1, nLat * nLat);
  arma::vec sigma_zeta_within_inds = seq_len_inds(beta_within_inds(1) + 1,  nLat * (nLat+1) / 2);
  arma::vec sigma_epsilon_within_inds =  seq_len_inds(sigma_zeta_within_inds(1) + 1,  nVar * (nVar + 1) / 2);
  
  arma::vec beta_between_inds = seq_len_inds(sigma_epsilon_within_inds(1) + 1, nLat * nLat);
  arma::vec sigma_zeta_between_inds = seq_len_inds(beta_between_inds(1) + 1, nLat * (nLat + 1) / 2);
  arma::vec sigma_epsilon_between_inds = seq_len_inds(sigma_zeta_between_inds(1) + 1, nVar * (nVar + 1) / 2 );
  
  // 
  //
  // AUGMENTATION PARTS //
  arma::mat aug_within_latent;
  arma::mat aug_within_residual;
  arma::mat aug_between_latent;
  arma::mat aug_between_residual;
  
  // Within latent
  
  // Within latent
  if (within_latent == "chol"){
    
    aug_within_latent = d_sigma_cholesky_cpp(
      grouplist["lowertri_zeta_within"],
               grouplist["L_eta"],
                        grouplist["C_eta_eta"],
                                 grouplist["I_eta"]);
    
  } else if (within_latent == "prec"){
    
    aug_within_latent = d_sigma_kappa_cpp(
      grouplist["L_eta"],grouplist["D_eta"],grouplist["sigma_zeta_within"]
    );
    
  } else if (within_latent == "ggm"){
    
    aug_within_latent = join_rows(
      d_sigma_omega_cpp(grouplist["L_eta"],grouplist["delta_IminOinv_zeta_within"],grouplist["A_eta"],grouplist["delta_zeta_within"],grouplist["Dstar_eta"]),
                                  d_sigma_delta_cpp(grouplist["L_eta"], grouplist["delta_IminOinv_zeta_within"], grouplist["I_eta"], grouplist["A_eta"])
    );
    
  } else if (within_latent == "cov"){
    
    aug_within_latent = eye(nLat*(nLat+1)/2, nLat*(nLat+1)/2);
    
  }
  
  // Between latent
  if (between_latent == "chol"){
    
    aug_between_latent = d_sigma_cholesky_cpp(
      grouplist["lowertri_zeta_between"],
               grouplist["L_eta"],
                        grouplist["C_eta_eta"],
                                 grouplist["I_eta"]);
    
  } else if (between_latent == "prec"){
    
    aug_between_latent = d_sigma_kappa_cpp(
      grouplist["L_eta"],grouplist["D_eta"],grouplist["sigma_zeta_between"]
    );
    
  } else if (between_latent == "ggm"){
    
    aug_between_latent = join_rows(
      d_sigma_omega_cpp(grouplist["L_eta"],grouplist["delta_IminOinv_zeta_between"],grouplist["A_eta"],grouplist["delta_zeta_between"],grouplist["Dstar_eta"]),
                                  d_sigma_delta_cpp(grouplist["L_eta"], grouplist["delta_IminOinv_zeta_between"], grouplist["I_eta"], grouplist["A_eta"])
    );
    
  } else if (between_latent == "cov"){
    
    aug_between_latent = eye(nLat*(nLat+1)/2, nLat*(nLat+1)/2);
    
  }
  
  // Within residual
  if (within_residual == "chol"){
    
    aug_within_residual = d_sigma_cholesky_cpp(
      grouplist["lowertri_epsilon_within"],
               grouplist["L_y"],
                        grouplist["C_y_y"],
                                 grouplist["I_y"]);
    
  } else if (within_residual == "prec"){
    
    aug_within_residual = d_sigma_kappa_cpp(
      grouplist["L_y"],grouplist["D_y"],grouplist["sigma_epsilon_within"]
    );
    
  } else if (within_residual == "ggm"){
    
    aug_within_residual = join_rows(
      d_sigma_omega_cpp(grouplist["L_y"],grouplist["delta_IminOinv_epsilon_within"],grouplist["A_y"],grouplist["delta_epsilon_within"],grouplist["Dstar_y"]),
                                  d_sigma_delta_cpp(grouplist["L_y"], grouplist["delta_IminOinv_epsilon_within"], grouplist["I_y"], grouplist["A_y"])
    );
    
  } else if (within_residual == "cov"){
    
    aug_within_residual = eye(nVar*(nVar+1)/2, nVar*(nVar+1)/2);
    
  }
  
  
  // Between residual
  if (between_residual == "chol"){
    
    aug_between_residual = d_sigma_cholesky_cpp(
      grouplist["lowertri_epsilon_between"],
               grouplist["L_y"],
                        grouplist["C_y_y"],
                                 grouplist["I_y"]);
    
  } else if (between_residual == "prec"){
    
    aug_between_residual = d_sigma_kappa_cpp(
      grouplist["L_y"],grouplist["D_y"],grouplist["sigma_epsilon_between"]
    );
    
  } else if (between_residual == "ggm"){
    
    aug_between_residual = join_rows(
      d_sigma_omega_cpp(grouplist["L_y"],grouplist["delta_IminOinv_epsilon_between"],grouplist["A_y"],grouplist["delta_epsilon_between"],grouplist["Dstar_y"]),
                                  d_sigma_delta_cpp(grouplist["L_y"], grouplist["delta_IminOinv_epsilon_between"], grouplist["I_y"], grouplist["A_y"])
    );
    
  } else if (between_residual == "cov"){
    
    aug_between_residual = eye(nVar*(nVar+1)/2, nVar*(nVar+1)/2);
    
  }
  
  
  // fill intercept part:
  Jac.submat(meanInds(0),nu_inds(0),meanInds(1),nu_inds(1)) =  eye(nVar, nVar);
  
  
  // Fill latent mean part:
  Jac.submat(meanInds(0),nu_eta_inds(0),meanInds(1),nu_eta_inds(1)) =  d_mu_nu_eta_lvm_cpp(grouplist["Lambda_BetaStar_between"]);
  
  // Fill mean to factor loading part:
  Jac.submat(meanInds(0),lambda_inds(0),meanInds(1),lambda_inds(1)) = d_mu_lambda_lvm_cpp(
    grouplist["nu_eta"], grouplist["BetaStar_between"], grouplist["I_y"]
  );
  
  // Fill s_in to lambda part:
  Jac.submat(sigInds_inPerson(0),lambda_inds(0),sigInds_inPerson(1),lambda_inds(1)) = 
    d_sigma_lambda_lvm_cpp(grouplist["L_y"],grouplist["Lambda_BetaStar_within"],grouplist["Betasta_sigmaZeta_within"],grouplist["I_y"], grouplist["C_y_eta"]) + 
    d_sigma_lambda_lvm_cpp(grouplist["L_y"],grouplist["Lambda_BetaStar_between"],grouplist["Betasta_sigmaZeta_between"],grouplist["I_y"],grouplist["C_y_eta"]);
  
  // beta
  Jac.submat(sigInds_inPerson(0),beta_within_inds(0),sigInds_inPerson(1),beta_within_inds(1)) = d_sigma_beta_lvm_cpp(
    grouplist["L_y"], grouplist["lambda"], grouplist["Betasta_sigmaZeta_within"], grouplist["C_eta_eta"], grouplist["I_eta"], grouplist["tBetakronBeta_within"]);
  
  
  // Fill s_in  to sigma_zeta_within part:
  Jac.submat(sigInds_inPerson(0),sigma_zeta_within_inds(0),sigInds_inPerson(1),sigma_zeta_within_inds(1)) = d_sigma_sigma_zeta_lvm_cpp(
    grouplist["L_y"],grouplist["Lambda_BetaStar_within"],grouplist["D_eta"]) * aug_within_latent;
  
    
  // Fill s_in to sigma_epsilon_within:
  Jac.submat(sigInds_inPerson(0),sigma_epsilon_within_inds(0),sigInds_inPerson(1),sigma_epsilon_within_inds(1)) = 
    aug_within_residual;
  
  // Between cluster model:
  
  
  // Beta part:
  arma::mat beta_between_part = d_sigma_beta_lvm_cpp(grouplist["L_y"], grouplist["lambda"], grouplist["Betasta_sigmaZeta_between"], 
                                                     grouplist["C_eta_eta"], grouplist["I_eta"], grouplist["tBetakronBeta_between"]);
  
  Jac.submat(sigInds_inPerson(0),beta_between_inds(0),sigInds_inPerson(1),beta_between_inds(1)) =   beta_between_part;
  
  
  // Fill s_in  to sigma_zeta_between part (and store for later use):
  arma::mat J_sigma_zeta_between = d_sigma_sigma_zeta_lvm_cpp(grouplist["L_y"],grouplist["Lambda_BetaStar_between"], grouplist["D_eta"]) * aug_between_latent;
  Jac.submat(sigInds_inPerson(0),sigma_zeta_between_inds(0),sigInds_inPerson(1),sigma_zeta_between_inds(1)) = J_sigma_zeta_between;
  
  
  // Fill s_in to sigma_epsilon_within:
  arma::mat J_sigma_epsilon_between = aug_between_residual;
  Jac.submat(sigInds_inPerson(0),sigma_epsilon_between_inds(0),sigInds_inPerson(1),sigma_epsilon_between_inds(1)) = J_sigma_epsilon_between;
  
  
  // Gradients for the between person blocks
  // Lambda
  arma::sp_mat I = arma::speye(nVar * nVar, nVar * nVar);
    Jac.submat(sigInds_betPersons(0),lambda_inds(0),sigInds_betPersons(1),lambda_inds(1))  =  
      d_sigma_lambda_lvm_cpp(I,grouplist["Lambda_BetaStar_between"], grouplist["Betasta_sigmaZeta_between"], grouplist["I_y"], grouplist["C_y_eta"]);
    
  
  // Beta part:
  Jac.submat(sigInds_betPersons(0),beta_between_inds(0),sigInds_betPersons(1),beta_between_inds(1)) = D_y * beta_between_part;

  // Fill s_in  to sigma_zeta_between part (and store for later use):
  Jac.submat(sigInds_betPersons(0),sigma_zeta_between_inds(0),sigInds_betPersons(1),sigma_zeta_between_inds(1)) = D_y * J_sigma_zeta_between;
  
  // Fill s_in to sigma_epsilon_within:
  Jac.submat(sigInds_betPersons(0),sigma_epsilon_between_inds(0),sigInds_betPersons(1),sigma_epsilon_between_inds(1)) = D_y * J_sigma_epsilon_between;
  
  
  
  // Permute:
  Jac = P * Jac;
  
  // Return:
  return(Jac);
}



// [[Rcpp::export]]
arma::mat d_phi_theta_ml_lvm_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_ml_lvm_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}




