#ifndef DVAR1_H
#define DVAR1_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_phi_theta_var1_cpp(
        const Rcpp::List& grouplist
);

// Individual VAR(1) derivative helpers (used by meta_var1):
arma::mat d_sigma0_beta_var1_cpp(
    const arma::mat& BetaStar,
    const arma::sp_mat& In,
    const arma::mat& sigma,
    const arma::sp_mat& C,
    const arma::sp_mat& L);

arma::mat d_sigma_zeta_cholesky_var1_cpp(
    const arma::mat& lowertri_zeta,
    const arma::sp_mat& L,
    const arma::sp_mat& C,
    const arma::sp_mat& In);

arma::mat d_sigma_zeta_kappa_var1_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& D2,
    const arma::mat& sigma_zeta);

arma::mat d_sigma_zeta_ggm_var1_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv_zeta,
    const arma::sp_mat& A,
    const arma::mat& delta_zeta,
    const arma::sp_mat& Dstar,
    const arma::sp_mat& In);

arma::mat d_sigma1_beta_var1_cpp(
    const arma::sp_mat& IkronBeta,
    const arma::sp_mat& D2,
    const arma::mat& Jb,
    const arma::mat& sigma,
    const arma::mat& beta,
    const arma::sp_mat& In);

arma::mat d_sigma1_sigma_zeta_var1_cpp(
    const arma::sp_mat& IkronBeta,
    const arma::sp_mat& D2,
    const arma::mat& Js);

#endif
