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
  bool corinput = false;
  bool meanstructure = true;

  if (grouplist.containsElementNamed("corinput")){
    corinput = grouplist["corinput"];
  }
  if (grouplist.containsElementNamed("meanstructure")){
    meanstructure = grouplist["meanstructure"];
  }

  // Number of variables:
  int nvar = lambda.n_rows;

  // Number of latents:
  int nlat = lambda.n_cols;

  // Number of thresholds:
  int nThresh = 0;
  arma::mat tau;
  if (grouplist.containsElementNamed("tau")){
    tau = Rcpp::as<arma::mat>(grouplist["tau"]);
    for (int j = 0; j < (int)tau.n_cols; j++){
      for (int i = 0; i < (int)tau.n_rows; i++){
        if (!R_IsNA(tau(i,j))){
          nThresh++;
        }
      }
    }
  }

  // Number of means:
  int nMeans = meanstructure ? nvar : 0;

  // Number of observation rows (means/thresh + vech(sigma)):
  int nMean_Thresh = nThresh + nMeans;
  int nobs = nMean_Thresh + (nvar * (nvar+1))/2;

  // total number of parameter columns:
  int nelement = nThresh + // Thresholds
    (meanstructure ? nvar : 0) + // Means (nu)
    (meanstructure ? nlat : 0) + // Latent intercepts (nu_eta)
    nvar * nlat + // factor loadings
    nlat * nlat + // beta elements
    nlat * (nlat + 1)/2 + // Latent variances and covariances
    nvar * ( nvar + 1)/2; // Residual network and scaling

  // Empty Jacobian:
  arma::mat Jac = zeros(nobs,nelement);

  // Observation indices:
  arma::vec meanInds = seq_len_inds(0, nMean_Thresh);
  arma::vec sigmaInds = seq_len_inds(nMean_Thresh, nvar*(nvar+1)/2);

  // Parameter indices:
  int curInd = 0;

  // Threshold indices:
  arma::vec tauInds;
  if (nThresh > 0){
    tauInds = seq_len_inds(curInd, nThresh);
    curInd = tauInds(1) + 1;
  }

  // Mean structure indices:
  arma::vec interceptInds;
  arma::vec nuetaInds;
  if (meanstructure){
    interceptInds = seq_len_inds(curInd, nvar);
    curInd = interceptInds(1) + 1;
    nuetaInds = seq_len_inds(curInd, nlat);
    curInd = nuetaInds(1) + 1;
  }

  arma::vec lambdaInds = seq_len_inds(curInd, nlat*nvar);
  curInd = lambdaInds(1) + 1;

  arma::vec betaInds = seq_len_inds(curInd, nlat*nlat);
  curInd = betaInds(1) + 1;

  arma::vec sigmazetaInds = seq_len_inds(curInd, nlat*(nlat+1)/2);
  curInd = sigmazetaInds(1) + 1;

  arma::vec sigmaepsilonInds = seq_len_inds(curInd, nvar*(nvar+1)/2);

  // Fill threshold part (identity):
  if (nThresh > 0){
    Jac.submat(meanInds(0), tauInds(0), meanInds(0) + nThresh - 1, tauInds(1)) = eye(nThresh, nThresh);
  }

  // Fill mean structure parts:
  if (meanstructure){
    // fill intercept part:
    Jac.submat(meanInds(0),interceptInds(0),meanInds(1),interceptInds(1)) =  d_mu_nu_lvm_cpp(
      grouplist["nu"]
    );

    // Fill latent intercept part:
    Jac.submat(meanInds(0),nuetaInds(0),meanInds(1),nuetaInds(1)) =  d_mu_nu_eta_lvm_cpp(
      grouplist["Lambda_BetaStar"]
    );

    // Fill factor loading parts (mean part):
    Jac.submat(meanInds(0),lambdaInds(0),meanInds(1),lambdaInds(1)) =  d_mu_lambda_lvm_cpp(
      grouplist["nu_eta"], grouplist["BetaStar"], grouplist["In"]
    );

    // Fill the beta parts (mean part):
    Jac.submat(meanInds(0),betaInds(0),meanInds(1),betaInds(1)) =  d_mu_beta_lvm_cpp(
      grouplist["nu_eta"], grouplist["lambda"], grouplist["tBetakronBeta"]
    );
  }

  // Fill factor loading parts (sigma part):
  Jac.submat(sigmaInds(0),lambdaInds(0),sigmaInds(1),lambdaInds(1)) =  d_sigma_lambda_lvm_cpp(
    grouplist["L"],  grouplist["Lambda_BetaStar"],  grouplist["Betasta_sigmaZeta"],  grouplist["In"],  grouplist["C"]
  );

  // Fill the beta parts (sigma part):
  Jac.submat(sigmaInds(0),betaInds(0),sigmaInds(1),betaInds(1)) =  d_sigma_beta_lvm_cpp(
    grouplist["L"],  grouplist["lambda"],  grouplist["Betasta_sigmaZeta"], grouplist["Cbeta"],  grouplist["Inlatent"],  grouplist["tBetakronBeta"]
  );

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

  // Cut out the rows not needed:
  if (corinput){
    // Remove diagonal variance rows from sigma part
    arma::uvec keepRows;
    // Keep all mean/threshold rows
    arma::uvec meanRows = arma::linspace<arma::uvec>(0, nMean_Thresh - 1, nMean_Thresh);
    // For sigma rows, keep only off-diagonal (lower triangular, not diagonal)
    arma::uvec offDiagRows;
    int idx = 0;
    for (int j = 0; j < nvar; j++){
      for (int i = j; i < nvar; i++){
        if (i != j){
          offDiagRows.resize(offDiagRows.n_elem + 1);
          offDiagRows(offDiagRows.n_elem - 1) = nMean_Thresh + idx;
        }
        idx++;
      }
    }
    if (nMean_Thresh > 0){
      keepRows = join_cols(meanRows, offDiagRows);
    } else {
      keepRows = offDiagRows;
    }
    Jac = Jac.rows(keepRows);
  }
  if (!meanstructure && nThresh == 0 && nMean_Thresh > 0){
    // Remove the first nvar rows (means) if no mean structure and no thresholds
    // Only when there are actual mean rows to remove (nMean_Thresh > 0)
    if (Jac.n_rows > (unsigned)nvar){
      Jac = Jac.rows(nvar, Jac.n_rows - 1);
    }
  }

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




