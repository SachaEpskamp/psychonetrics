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

// kronecker_diag
// kronecker_X_I
// kronecker_I_X


// arma::mat d_sigma_cholesky_cpp(
//     const arma::mat& lowertri,
//     const arma::sp_mat& L,
//     const arma::sp_mat& C,
//     const arma::sp_mat In
// ){
//   // return(L * (kron(In,In) + C) * (kron(lowertri, (arma::mat)In) * L.t()));
//   
//   arma::sp_mat res = L * (kronecker_diag(In) + C) * (kronecker_X_I(lowertri, In.n_rows) * L.t());
//   return((arma::mat)res);
//   
// }


// [[Rcpp::export]]
arma::mat d_mu_mu_var1_cpp(
    const arma::mat& beta)
{
  
  return eye(beta.n_rows, beta.n_cols);
}


// [[Rcpp::export]]
arma::mat d_sigmastar_exo_cholesky_var1_cpp(
    const arma::sp_mat& In, 
    const arma::sp_mat& L, 
    const arma::sp_mat& C, 
    const arma::mat& exo_cholesky){
  
  arma::mat res(L * (
      kronecker_I_X(exo_cholesky,  In.n_rows) * C * L.t() + 
        kronecker_X_I(exo_cholesky,  In.n_rows) * L.t()
  ));
  
  return(res);
}

// [[Rcpp::export]]
arma::mat d_sigma0_beta_var1_cpp(
    const arma::mat& BetaStar,
    const arma::sp_mat& In, 
    const arma::mat& sigma,
    const arma::sp_mat& C, 
    const arma::sp_mat& L){
  
  int n = In.n_rows;
  arma::mat sigma1 = sigma.submat(n,0,2*n-1,(n-1));
  
  arma::sp_mat InXIn = speye(n*n,n*n);
  
  arma::mat res = (L * (InXIn + C)) * BetaStar * kronecker_X_I(sigma1,n); 
  return(res);
}




// Derivative of sigma_zeta to cholesky:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_cholesky_var1_cpp(
    const arma::mat& lowertri_zeta,
    const arma::sp_mat& L,
    const arma::sp_mat& C,
    const arma::sp_mat& In){
  
  return(d_sigma_cholesky_cpp(lowertri_zeta,L,C,In));
}

// Derivative of sigma_zeta to precision:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_kappa_var1_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& D2, 
    const arma::mat& sigma_zeta){
  
  return(d_sigma_kappa_cpp(L, D2, sigma_zeta));
}

// Derivative of sigma_zeta to ggm:
// [[Rcpp::export]]
arma::mat d_sigma_zeta_ggm_var1_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv_zeta,
    const arma::sp_mat& A,
    const arma::mat& delta_zeta,
    const arma::sp_mat& Dstar,
    const arma::sp_mat& In){
  
  
  arma::mat matomega = d_sigma_omega_cpp(L, delta_IminOinv_zeta, A, delta_zeta, Dstar);
  
  arma::mat matdelta = d_sigma_delta_cpp(L,  delta_IminOinv_zeta, In, A);
  
  
  return(join_rows(matomega, matdelta));
}



// Derivative of sigma1 with respect to beta:
// [[Rcpp::export]]
arma::mat d_sigma1_beta_var1_cpp(
    const arma::sp_mat& IkronBeta,
    const arma::sp_mat& D2,
    const arma::mat& Jb,
    const arma::mat& sigma,
    const arma::mat& beta,
    const arma::sp_mat& In){
  int n = beta.n_rows;
  arma::mat sigma0 = sigma.submat(n,n,2*n-1,2*n-1);
  return( 
    IkronBeta * D2 * Jb + kronecker_X_I(sigma0,In.n_rows)
  );
}


// Derivative of sigma1 with respect to omega:
// [[Rcpp::export]]
arma::mat d_sigma1_sigma_zeta_var1_cpp(
    const arma::sp_mat& IkronBeta,
    const arma::sp_mat& D2,
    const arma::mat& Js){
  
  return(IkronBeta * D2 * Js);
}



// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_var1_group_cpp(
    const Rcpp::List& grouplist
){
  // objects needed now:
  arma::sp_mat P = grouplist["P"];
  arma::mat beta = grouplist["beta"];
  std::string zeta = grouplist["zeta"];
  
  
  // Number of variables:
  int nvar = beta.n_rows * 2;
  
  // Number of nodes:
  int nNode = nvar / 2;
  
  // Number of observations:
  int nobs = nvar + // Means
    (nvar * (nvar+1))/2; // Variances
  
  // total number of elements:
  int nelement = nvar + (nNode*(nNode+1)/2) +  (nNode * nNode) + (nNode * (nNode+1) / 2); // Contemporaneous network and var-cov
  
  
  // Empty Jacobian:
  arma::mat Jac = zeros(nobs,nelement);
  
  
  // Indices:
  int meanInds_start = 0; 
  int meanInds_end = nvar - 1;
  int sigmaStarInds_start = nvar;
  int sigmaStarInds_end = nvar + nNode*(nNode+1)/2 - 1;
  int sigma0Inds_start = sigmaStarInds_end + 1;
  int sigma0Inds_end = sigma0Inds_start + nNode*(nNode+1)/2 - 1;
  int sigma1Inds_start =  sigma0Inds_end + 1;
  int sigma1Inds_end = sigma1Inds_start + nNode * nNode - 1;
  
  // Indices model:
  int interceptInds_start = 0;
  int interceptInds_end = nvar - 1;
  int exovarInds_start = nvar;
  int exovarInds_end = nvar + nNode * (nNode + 1)/2 - 1;
  int betaInds_start = exovarInds_end + 1;
  int betaInds_end = betaInds_start + nNode * nNode - 1;
  int sigmazetaInds_start = betaInds_end + 1;
  int sigmazetaInds_end = sigmazetaInds_start + nNode * (nNode + 1)/2 - 1;
  

  
  // fill intercept part:
  Rcpp::List dummylist = List::create(d_mu_mu_var1_cpp(grouplist["beta"]), d_mu_mu_var1_cpp(grouplist["beta"]));
  
  Jac.submat(meanInds_start,interceptInds_start,meanInds_end,interceptInds_end ) =  bdiag_psychonetrics(dummylist);

  
  // Fill exo var part:
  Jac.submat(sigmaStarInds_start,exovarInds_start,sigmaStarInds_end,exovarInds_end ) =   d_sigmastar_exo_cholesky_var1_cpp(
    grouplist["In"], grouplist["L"], grouplist["C"], grouplist["exo_cholesky"]
  );
  
  // Fill sigma0 to beta part:
  arma::mat Jb = d_sigma0_beta_var1_cpp(
    grouplist["BetaStar"], grouplist["In"], grouplist["sigma"], grouplist["C"], grouplist["L"]
  );
  
  Jac.submat(sigma0Inds_start,betaInds_start,sigma0Inds_end,betaInds_end) =  Jb;
  
  
  // Fill sigma0 to sigma_zeta part:
  Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) =  d_sigma0_sigma_zeta_var1_cpp(
    grouplist["L"], grouplist["BetaStar"], grouplist["D2"]
  );
  
  
  // Augment:
  
  if (zeta == "chol"){
    
    Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) =  
      Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) * 
      d_sigma_zeta_cholesky_var1_cpp(
        grouplist["lowertri_zeta"], grouplist["L"], grouplist["C"],  grouplist["In"]
      );
    
  } else if (zeta == "prec"){
    
    Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) =  
      Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) * 
      d_sigma_zeta_kappa_var1_cpp(
        grouplist["L"], grouplist["D2"], grouplist["sigma_zeta"]
      );
    
    
  } else if (zeta == "ggm"){
    
    Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) =  
      Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end) * 
      d_sigma_zeta_ggm_var1_cpp(
        grouplist["L"], grouplist["delta_IminOinv_zeta"], grouplist["A"],grouplist["delta_zeta"], grouplist["Dstar"], grouplist["In"]
    
      );
    
  }
  
  // Store:
  arma::mat Js = Jac.submat(sigma0Inds_start,sigmazetaInds_start,sigma0Inds_end,sigmazetaInds_end);
  
  // Fill sigma1 to beta part:
  Jac.submat(sigma1Inds_start,betaInds_start,sigma1Inds_end,betaInds_end) = d_sigma1_beta_var1_cpp(
    grouplist["IkronBeta"],  grouplist["D2"],  Jb,  grouplist["sigma"],  grouplist["beta"],  grouplist["In"]
  );
  
  // Fil to zeta part:
  Jac.submat(sigma1Inds_start,sigmazetaInds_start,sigma1Inds_end,sigmazetaInds_end) = d_sigma1_sigma_zeta_var1_cpp(
    grouplist["IkronBeta"],  grouplist["D2"],  Js
  );
  
  
  
  // Permute:
  Jac = P * Jac;
  
  // Return:
  return(Jac);
}
// 
// [[Rcpp::export]]
arma::mat d_phi_theta_var1_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_var1_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}




