#ifndef LVMCPP_H
#define LVMCPP_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;




arma::mat d_mu_nu_lvm_cpp(
    const arma::mat& nu);



arma::mat d_mu_nu_eta_lvm_cpp(
    arma::mat Lambda_BetaStar
);


arma::mat d_mu_lambda_lvm_cpp(
    const arma::mat& nu_eta,
    const arma::mat& BetaStar,
    const arma::sp_mat& In);


arma::mat d_mu_beta_lvm_cpp(
    const arma::mat& nu_eta,
    const arma::mat& lambda,
    const arma::mat& tBetakronBeta);



arma::mat d_sigma_lambda_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& Lambda_BetaStar,
    const arma::mat& Betasta_sigmaZeta,
    const arma::sp_mat& In,
    const arma::sp_mat& C);




arma::mat d_sigma_beta_lvm_cpp(
    const arma::sp_mat& L, 
    const arma::mat& lambda, 
    const arma::mat& Betasta_sigmaZeta, 
    const arma::sp_mat& Cbeta,
    const arma::sp_mat& Inlatent,
    const arma::mat& tBetakronBeta);



arma::mat d_sigma_sigma_zeta_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& Lambda_BetaStar,
    const arma::sp_mat& Deta);



arma::mat d_sigma_zeta_cholesky_lvm_cpp(
    const arma::mat& lowertri_zeta,
    const arma::sp_mat& L_eta,
    const arma::sp_mat& Cbeta,
    const arma::sp_mat& Inlatent);


arma::mat d_sigma_zeta_kappa_lvm_cpp(
    const arma::sp_mat& L_eta,
    const arma::sp_mat& Deta,
    const arma::mat& sigma_zeta);


arma::mat d_sigma_zeta_ggm_lvm_cpp(
    const arma::sp_mat& L_eta,
    const arma::mat& delta_IminOinv_zeta,
    const arma::sp_mat& Aeta,
    const arma::mat& delta_zeta,
    const arma::sp_mat& Dstar_eta,
    const arma::sp_mat& Inlatent);


arma::mat d_sigma_epsilon_cholesky_lvm_cpp(
    const arma::mat& lowertri_epsilon,
    const arma::sp_mat& L,
    const arma::sp_mat& C_chol,
    const arma::sp_mat& In);


arma::mat d_sigma_epsilon_kappa_lvm_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& D,
    const arma::mat& sigma_epsilon);


arma::mat d_sigma_epsilon_ggm_lvm_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv_epsilon,
    const arma::sp_mat& A,
    const arma::mat& delta_epsilon,
    const arma::sp_mat& Dstar,
    const arma::sp_mat& In);



arma::mat d_phi_theta_lvm_group_cpp(
    const Rcpp::List& grouplist
);



arma::mat d_phi_theta_lvm_cpp(
    const Rcpp::List& prep
);



#endif
