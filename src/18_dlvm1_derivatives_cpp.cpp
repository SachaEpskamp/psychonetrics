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
arma::mat d_mu_lambda_dlvm1_cpp(
    const arma::mat& mu_eta,
    const arma::sp_mat& I_y){
  arma::mat res = (arma::mat)kronecker_X_I(mu_eta.t(), I_y.n_rows);
  return(res);
}


// [[Rcpp::export]]
arma::mat d_sigmak_lambda_dlvm1_cpp(
    const arma::mat& lambda, 
    int k,  
    const Rcpp::List allSigmas_within,
    const arma::sp_mat& C_y_eta, 
    const arma::sp_mat& I_y, 
    const arma::sp_mat& L_y,
    const arma::mat& sigma_zeta_between){
  
  arma::mat sigmak = allSigmas_within[k];
  
  arma::sp_mat within = (
    kronecker_I_X(lambda * sigmak, I_y.n_rows) * C_y_eta + 
      kronecker_X_I(lambda * sigmak.t(),I_y.n_rows)  
  );
  
  arma::sp_mat between = kronecker_I_X(lambda * sigma_zeta_between,I_y.n_rows) * C_y_eta + 
    kronecker_X_I(lambda * sigma_zeta_between, I_y.n_rows);  
  
  arma::mat res = (arma::mat)(within + between);
  
  if (k == 0){
    return(L_y * res);
  } else {
    return(res);
  }
}


// [[Rcpp::export]]
arma::mat d_sigma0_sigma_zeta_within_dlvm1_cpp(
    const arma::mat& BetaStar, 
    const arma::sp_mat& D_eta){
  
  arma::mat res = BetaStar * D_eta;
  
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigma0_beta_dlvm1_cpp(
    const arma::mat& BetaStar,
    const arma::sp_mat& I_eta,
    const Rcpp::List allSigmas_within,
    const arma::sp_mat& C_eta_eta){
  
  arma::mat sigma1 = allSigmas_within[1];
  
  int n = I_eta.n_rows;
  
  arma::mat res = (eye(n * n, n * n) + C_eta_eta) * BetaStar * kronecker_X_I(sigma1, I_eta.n_rows); 
  
  
  
  return(res);
}


// [[Rcpp::export]]
arma::mat d_sigmak_beta_dlvm1_cpp(
    const arma::mat& J_sigma_beta,
    const arma::sp_mat& IkronBeta, 
    int k,  
    Rcpp::List allSigmas_within,
    const arma::sp_mat& I_eta){
  
  arma::mat sigma = allSigmas_within[k-1];
  
  arma::mat res = IkronBeta * J_sigma_beta + kronecker_X_I(sigma.t(), I_eta.n_rows);
  
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigmak_sigma_zeta_between_dlvm1_cpp(
    const arma::mat& lambda,
    const arma::sp_mat& D_eta){
  
  arma::mat res = kron(lambda, lambda) * D_eta;
  
  return(res);
}



// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_dlvm1_group_cpp(
    const Rcpp::List& grouplist
){
  // Setup:
  int t;
  
  
  // Things I need now:
  std::string within_latent  = grouplist["within_latent"];
  std::string within_residual  = grouplist["within_residual"];
  std::string between_latent  = grouplist["between_latent"];
  std::string between_residual  = grouplist["between_residual"];
  
  arma::mat lambda  = grouplist["lambda"];
  
  
  arma::mat design = grouplist["design"];
  arma::sp_mat P = grouplist["P"];
  arma::sp_mat L_eta = grouplist["L_eta"];
  arma::sp_mat L_y = grouplist["L_y"];
  arma::sp_mat D_y = grouplist["D_y"];
  
  
  // Number of variables:
  int nVar = design.n_rows;
  
  // Number of times points:
  int nTime = design.n_cols;
  
  // Number of latents:
  int nLat = lambda.n_cols;
  
  // Number of within latents:
  // nLat <- ncol(lambda)
  
  // Number of between latents:
  // nLat <- ncol(lambda)
  
  // I need to construct the Jacobian for the following "observations:"
  int nobs = nVar + // Means
    (nVar * (nVar+1))/2 + // Variances
    (nVar * nVar) * (nTime - 1); // lagged variances
  
  // total number of elements:
  int nelement = nVar + // intercepts in nu
    nLat + // Latent means
    nVar * nLat + // Factor loadings
    nLat * (nLat+1) / 2 + // within factor variances
    (nLat * nLat) + // temporal effects
    nVar * (nVar + 1) / 2 + // Within residuals
    // nVar * nLat + // Between-subject factor loadings
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
  
  Rcpp::List sigInds(nTime);
  sigInds[0] = seq_len_inds(meanInds(1) + 1, nVar*(nVar+1)/2);
  
  // For each lag:
  for (t=1; t<nTime;t++){
    arma::vec prev = sigInds[t-1];
    sigInds[t]  = seq_len_inds(prev(1) + 1, nVar*nVar);
  }
  
  
  // Indices model:
  arma::vec nu_inds = seq_len_inds(0, nVar);
  arma::vec mu_eta_inds = seq_len_inds(nu_inds(1) + 1, nLat);
  arma::vec lambda_inds = seq_len_inds(mu_eta_inds(1) + 1, nVar * nLat);
  

  
  // lambda_inds <- max(mu_eta_inds) + seq_len(nVar * nLat)
  arma::vec sigma_zeta_within_inds = seq_len_inds(lambda_inds(1) + 1,  nLat * (nLat+1) / 2);
  arma::vec beta_inds = seq_len_inds(sigma_zeta_within_inds(1) + 1, nLat * nLat);
  arma::vec sigma_epsilon_within_inds =  seq_len_inds(beta_inds(1) + 1,  nVar * (nVar + 1) / 2);
  
  arma::vec sigma_zeta_between_inds = seq_len_inds(sigma_epsilon_within_inds(1) + 1, nLat * (nLat + 1) / 2);
  arma::vec sigma_epsilon_between_inds = seq_len_inds(sigma_zeta_between_inds(1) + 1, nVar * (nVar + 1) / 2 );
  
  
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
  Jac.submat(meanInds(0),mu_eta_inds(0),meanInds(1),mu_eta_inds(1)) =  lambda;
  
  // Fill mean to factor loading part:
  Jac.submat(meanInds(0),lambda_inds(0),meanInds(1),lambda_inds(1)) = d_mu_lambda_dlvm1_cpp(grouplist["mu_eta"], grouplist["I_y"]);
  
  // Fill s0 to lambda part:
  arma::vec curSigInds = sigInds[0];
  Jac.submat(curSigInds(0),lambda_inds(0),curSigInds(1),lambda_inds(1)) = d_sigmak_lambda_dlvm1_cpp(lambda, 0, grouplist["allSigmas_within"], grouplist["C_y_eta"], grouplist["I_y"], grouplist["L_y"], grouplist["sigma_zeta_between"]);
  
  
  // Fill s0  to sigma_zeta_within part (and store for later use):
  arma::mat lamWkronlamW = grouplist["lamWkronlamW"];
  
  
  
  arma::mat J_sigma_zeta_within = d_sigma0_sigma_zeta_within_dlvm1_cpp( 
    grouplist["BetaStar"], grouplist["D_eta"]
  ) * aug_within_latent;
  
  
  
  Jac.submat(curSigInds(0),sigma_zeta_within_inds(0),curSigInds(1),sigma_zeta_within_inds(1)) = L_y * lamWkronlamW * J_sigma_zeta_within;
  
  // Fill s0 to beta part (and store for later use):
  arma::mat J_sigma_beta = d_sigma0_beta_dlvm1_cpp( grouplist["BetaStar"],  grouplist["I_eta"],  grouplist["allSigmas_within"],  grouplist["C_eta_eta"]);
  Jac.submat(curSigInds(0),beta_inds(0),curSigInds(1),beta_inds(1)) =  L_y  * lamWkronlamW * J_sigma_beta;
  
  // Fill s0 to sigma_epsilon_within:
  Jac.submat(curSigInds(0),sigma_epsilon_within_inds(0),curSigInds(1),sigma_epsilon_within_inds(1)) = aug_within_residual;
  
  
  // Fill s0 to sigma_zeta_between, and store for later use:
  arma::mat J_sigmak_sigma_zeta_between = d_sigmak_sigma_zeta_between_dlvm1_cpp(
    lambda,  grouplist["D_eta"]
  )  * aug_between_latent;
  
  Jac.submat(curSigInds(0),sigma_zeta_between_inds(0),curSigInds(1),sigma_zeta_between_inds(1)) =  L_y * J_sigmak_sigma_zeta_between;
  
  
  // Fill sigma_epsilon_between inds:
  Jac.submat(curSigInds(0),sigma_epsilon_between_inds(0),curSigInds(1),sigma_epsilon_between_inds(1))  =  aug_between_residual;
  
  // Obtain IkronBeta:
  arma::sp_mat IkronBeta = grouplist["IkronBeta"];
  
  // Now for every further lag ...
  for (t=1; t<nTime;t++){
    arma::vec curSigInds = sigInds[t];
    
    // Fill sk to lambda part:
    Jac.submat(curSigInds(0),lambda_inds(0),curSigInds(1),lambda_inds(1))  = d_sigmak_lambda_dlvm1_cpp(lambda, t, grouplist["allSigmas_within"], grouplist["C_y_eta"], grouplist["I_y"], grouplist["L_y"], grouplist["sigma_zeta_between"]);

    // Fill sk  to sigma_zeta_within part (and store for later use):
    J_sigma_zeta_within = IkronBeta * J_sigma_zeta_within;
    Jac.submat(curSigInds(0),sigma_zeta_within_inds(0),curSigInds(1),sigma_zeta_within_inds(1))  = lamWkronlamW * J_sigma_zeta_within;
    
    // Fill sk to beta part (and store for later use):
    J_sigma_beta = d_sigmak_beta_dlvm1_cpp(J_sigma_beta,IkronBeta,t, grouplist["allSigmas_within"],  grouplist["I_eta"]);
    Jac.submat(curSigInds(0),beta_inds(0),curSigInds(1),beta_inds(1)) = lamWkronlamW * J_sigma_beta;
    
    // Fill s0 to sigma_zeta_between
    Jac.submat(curSigInds(0),sigma_zeta_between_inds(0),curSigInds(1),sigma_zeta_between_inds(1)) = J_sigmak_sigma_zeta_between;

    // Fill s0 to sigma_zeta_between, and store for later use:
    Jac.submat(curSigInds(0),sigma_zeta_between_inds(0),curSigInds(1),sigma_zeta_between_inds(1)) = J_sigmak_sigma_zeta_between;
      
    // Fill sigma_epsilon_between inds:
    Jac.submat(curSigInds(0),sigma_epsilon_between_inds(0),curSigInds(1),sigma_epsilon_between_inds(1)) = D_y * aug_between_residual;
    
  }
  
  // Permute:
  Jac = P * Jac;
  
  // Return:
  return(Jac);
}



// [[Rcpp::export]]
arma::mat d_phi_theta_dlvm1_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_dlvm1_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}




