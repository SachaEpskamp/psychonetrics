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
arma::mat d_mu_nu_lvm_cpp(
    const arma::mat& nu){
  arma::mat res = eye(nu.n_rows, nu.n_rows);
  return(res);
}


// [[Rcpp::export]]
arma::mat d_mu_nu_eta_lvm_cpp(
    arma::mat Lambda_BetaStar
){
  arma::mat res(Lambda_BetaStar);
  return(res);
}

// [[Rcpp::export]]
arma::mat d_mu_lambda_lvm_cpp(
    const arma::mat& nu_eta,
    const arma::mat& BetaStar,
    const arma::sp_mat& In){
  
  arma::sp_mat res = kronecker_X_I(
    nu_eta.t() * BetaStar.t()
  , In.n_rows
  );
  
  return((arma::mat)res);
}



// [[Rcpp::export]]
arma::mat d_mu_beta_lvm_cpp(
    const arma::mat& nu_eta,
    const arma::mat& lambda,
    const arma::mat& tBetakronBeta){
  
  arma::mat res = kron(nu_eta.t(), lambda) * tBetakronBeta;
  
  return(res);
}


// Derivative of factor loadings to vars:
// [[Rcpp::export]]
arma::mat d_sigma_lambda_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& Lambda_BetaStar,
    const arma::mat& Betasta_sigmaZeta,
    const arma::sp_mat& In,
    const arma::sp_mat& C){
  

  arma::mat inner = Lambda_BetaStar * (Betasta_sigmaZeta.t());
  
  arma::sp_mat res = L * (
    kronecker_X_I(inner, In.n_rows ) + 
      (kronecker_I_X(inner, In.n_rows) * C)
  );
  
  
  return((arma::mat)res);
}




// Derivative of beta to vars:
// [[Rcpp::export]]
arma::mat d_sigma_beta_lvm_cpp(
    const arma::sp_mat& L, 
    const arma::mat& lambda, 
    const arma::mat& Betasta_sigmaZeta, 
    const arma::sp_mat& Cbeta,
    const arma::sp_mat& Inlatent,
    const arma::mat& tBetakronBeta){
  
  arma::mat res = L * (kron(lambda, lambda) * (arma::mat)(
    kronecker_X_I(Betasta_sigmaZeta,Inlatent.n_rows) + 
      kronecker_I_X(Betasta_sigmaZeta,Inlatent.n_rows) * Cbeta
  ) * tBetakronBeta);
  
  
  return(res);
}


// Derivative of latent variance-covariance matrix:
// [[Rcpp::export]]
arma::mat d_sigma_sigma_zeta_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& Lambda_BetaStar,
    const arma::sp_mat& Deta){
  arma::mat res = L * kron(Lambda_BetaStar, Lambda_BetaStar) * Deta;
  
  return(res);
}

// Derivative of sigma_zeta to cholesky:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_cholesky_lvm_cpp(
    const arma::mat& lowertri_zeta,
    const arma::sp_mat& L_eta,
    const arma::sp_mat& Cbeta,
    const arma::sp_mat& Inlatent){
  arma::mat res = d_sigma_cholesky_cpp(lowertri_zeta,L_eta,Cbeta,Inlatent);
  return(res);
}

// Derivative of sigma_zeta to precision:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_kappa_lvm_cpp(
    const arma::sp_mat& L_eta,
    const arma::sp_mat& Deta,
    const arma::mat& sigma_zeta){
  arma::mat res = d_sigma_kappa_cpp(L_eta,Deta,sigma_zeta);
  
  return(res);
}


// Derivative of sigma_zeta to ggm:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_ggm_lvm_cpp(
    const arma::sp_mat& L_eta,
    const arma::mat& delta_IminOinv_zeta,
    const arma::sp_mat& Aeta,
    const arma::mat& delta_zeta,
    const arma::sp_mat& Dstar_eta,
    const arma::sp_mat& Inlatent){
  arma::mat res = join_rows(
    d_sigma_omega_cpp(L_eta, delta_IminOinv_zeta, Aeta, delta_zeta, Dstar_eta),
    d_sigma_delta_cpp(L_eta, delta_IminOinv_zeta, Inlatent, Aeta)
  );
  
  return(res);
}


// Residual vars
// [[Rcpp::export]]
arma::mat d_sigma_epsilon_cholesky_lvm_cpp(
    const arma::mat& lowertri_epsilon,
    const arma::sp_mat& L,
    const arma::sp_mat& C_chol,
    const arma::sp_mat& In){
  arma::mat res = d_sigma_cholesky_cpp(lowertri_epsilon,L,C_chol,In);
  return(res);
}

// Derivative of sigma_epsilon to precision:
// [[Rcpp::export]]
arma::mat d_sigma_epsilon_kappa_lvm_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& D,
    const arma::mat& sigma_epsilon){
  arma::mat res = d_sigma_kappa_cpp(L, D, sigma_epsilon);
  
  return(res);
}

// Derivative of sigma_epsilon to ggm:
// [[Rcpp::export]]
arma::mat d_sigma_epsilon_ggm_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv_epsilon,
    const arma::sp_mat& A,
    const arma::mat& delta_epsilon,
    const arma::sp_mat& Dstar,
    const arma::sp_mat& In){
  arma::mat res = join_rows(
    d_sigma_omega_cpp(L, delta_IminOinv_epsilon, A, delta_epsilon, Dstar),
    d_sigma_delta_cpp(L, delta_IminOinv_epsilon, In, A)
  );
  
  return(res);
}


// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_lvm_group_cpp(
    const Rcpp::List& grouplist
){
  // Some things I need now:
  arma::mat lambda = grouplist["lambda"];
  std::string latent = grouplist["latent"];
  std::string residual = grouplist["residual"];
  
  // Number of variables:
  int nvar = lambda.n_rows;
  
  // Number of latents:
  int nlat = lambda.n_cols;
  
  // Number of observations:
  int nobs = nvar + // Means
    (nvar * (nvar+1))/2; // Variances
  
  // total number of elements:
  int nelement = nvar + // Means
    nlat + // Latent intercpets
    nvar * nlat + // factor loadings
    nlat * nlat + // beta elements
    nlat * (nlat + 1)/2 + // Latent variances and covariances
    nvar * ( nvar + 1)/2; // Residual network and scaling
  
  // Empty Jacobian:
  arma::mat Jac = zeros(nobs,nelement);
  
  // Indices:
  arma::vec meanInds = seq_len_inds(0, nvar);
  arma::vec sigmaInds = seq_len_inds(nvar, nvar*(nvar+1)/2); 
  
  
  // Indices model:
  arma::vec interceptInds = seq_len_inds(0, nvar);
  
  arma::vec nuetaInds = seq_len_inds(interceptInds(1) + 1, nlat);
  
  arma::vec lambdaInds = seq_len_inds(nuetaInds(1) + 1, nlat*nvar);  
  
  arma::vec betaInds = seq_len_inds(lambdaInds(1) + 1, nlat*nlat);    
  
  arma::vec sigmazetaInds = seq_len_inds(betaInds(1) + 1,nlat*(nlat+1)/2);
  
  arma::vec sigmaepsilonInds = seq_len_inds(sigmazetaInds(1) + 1, nvar*(nvar+1)/2);
  
  // fill intercept part:
  Jac.submat(meanInds(0),interceptInds(0),meanInds(1),interceptInds(1)) =  d_mu_nu_lvm_cpp( 
    grouplist["nu"]
  );
  
  
  // Fill latent intercept part:
  Jac.submat(meanInds(0),nuetaInds(0),meanInds(1),nuetaInds(1)) =  d_mu_nu_eta_lvm_cpp(
    grouplist["Lambda_BetaStar"]
  );
  
  // Fill factor loading parts:
  Jac.submat(meanInds(0),lambdaInds(0),meanInds(1),lambdaInds(1)) =  d_mu_lambda_lvm_cpp( 
    grouplist["nu_eta"], grouplist["BetaStar"], grouplist["In"]
  );
  

  Jac.submat(sigmaInds(0),lambdaInds(0),sigmaInds(1),lambdaInds(1)) =  d_sigma_lambda_lvm_cpp( 
    grouplist["L"],  grouplist["Lambda_BetaStar"],  grouplist["Betasta_sigmaZeta"],  grouplist["In"],  grouplist["C"]
  );
  
  
  
  // Fill the beta parts:
  Jac.submat(meanInds(0),betaInds(0),meanInds(1),betaInds(1)) =  d_mu_beta_lvm_cpp( 
    grouplist["nu_eta"], grouplist["lambda"], grouplist["tBetakronBeta"]
  );
  Jac.submat(sigmaInds(0),betaInds(0),sigmaInds(1),betaInds(1)) =  d_sigma_beta_lvm_cpp( 
    grouplist["L"],  grouplist["lambda"],  grouplist["Betasta_sigmaZeta"], grouplist["Cbeta"],  grouplist["Inlatent"],  grouplist["tBetakronBeta"]                                                
  );
  // 
  // 
  
  // Fill latent variances part:
  Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) =  d_sigma_sigma_zeta_lvm_cpp( 
    grouplist["L"],  grouplist["Lambda_BetaStar"],  grouplist["Deta"]
  );


  if (latent == "chol"){
    Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) =  
      Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) * 
      d_sigma_zeta_cholesky_lvm_cpp( 
      grouplist["lowertri_zeta"],  grouplist["L_eta"],  grouplist["Cbeta"],  grouplist["Inlatent"]
    );
  
  } else if (latent == "prec"){
    
    Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) =  
      Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) * 
      d_sigma_zeta_kappa_lvm_cpp( 
        grouplist["L_eta"],  grouplist["Deta"],  grouplist["sigma_zeta"]
      );
    
  } else if (latent == "ggm"){
    
    Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) =  
      Jac.submat(sigmaInds(0),sigmazetaInds(0),sigmaInds(1),sigmazetaInds(1)) * 
      d_sigma_zeta_ggm_lvm_cpp( 
        grouplist["L_eta"],  grouplist["delta_IminOinv_zeta"],  grouplist["Aeta"],
                  grouplist["delta_zeta"],grouplist["Dstar_eta"],grouplist["Inlatent"]
      );
  }

  // Residual variances:
  if (residual == "cov"){
    
    Jac.submat(sigmaInds(0),sigmaepsilonInds(0),sigmaInds(1),sigmaepsilonInds(1)) = eye(nvar*(nvar+1)/2, nvar*(nvar+1)/2);
    
    
  } else if (residual == "chol"){
    
    Jac.submat(sigmaInds(0),sigmaepsilonInds(0),sigmaInds(1),sigmaepsilonInds(1)) = d_sigma_epsilon_cholesky_lvm_cpp(
      grouplist["lowertri_epsilon"], grouplist["L"], grouplist["C_chol"], grouplist["In"]
    );
    
  } else if (residual == "prec"){
    
    Jac.submat(sigmaInds(0),sigmaepsilonInds(0),sigmaInds(1),sigmaepsilonInds(1)) = d_sigma_epsilon_kappa_lvm_cpp(
      grouplist["L"], grouplist["D"], grouplist["sigma_epsilon"]
    );
    
    
  } else if (residual == "ggm"){
    
    Jac.submat(sigmaInds(0),sigmaepsilonInds(0),sigmaInds(1),sigmaepsilonInds(1)) = d_sigma_epsilon_ggm_lvm_cpp(
      grouplist["L"], grouplist["delta_IminOinv_epsilon"], grouplist["A"], 
      grouplist["delta_epsilon"], grouplist["Dstar"], grouplist["In"]
    );
    

  }
  // 
  
  // Return:
  return(Jac);
}
// 
// [[Rcpp::export]]
arma::mat d_phi_theta_lvm_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_lvm_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}




