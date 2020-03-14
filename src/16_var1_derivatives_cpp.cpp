// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "14_varcov_derivatives_cpp.h"


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
  
  arma::mat res = L * (InXIn + C) * BetaStar * kronecker_X_I(sigma1,n); 
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
    const arma::sp_mat& delta_zeta,
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






