#ifndef VARCOVCPP_H
#define VARCOVCPP_H

#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
#include <pbv.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat d_sigma_cholesky_cpp(
                const arma::mat& lowertri,
                const arma::sp_mat& L,
                const arma::sp_mat& C,
                const arma::sp_mat In
);

arma::mat d_sigma_delta_cpp(
                const arma::sp_mat& L,
                const arma::mat& delta_IminOinv,
                const arma::sp_mat& In,
                const arma::sp_mat& A
);

arma::mat d_sigma_omega_cpp(
                const arma::sp_mat& L,
                const arma::mat& delta_IminOinv,
                const arma::sp_mat& A,
                const arma::mat& delta,
                const arma::sp_mat& Dstar
);

arma::mat d_sigma_kappa_cpp(
                const arma::sp_mat& L,
                const arma::sp_mat& D,
                const arma::mat& sigma);

arma::mat d_sigma_rho_cpp(
                const arma::sp_mat& L,
                const arma::mat& SD,
                const arma::sp_mat& A,
                const arma::sp_mat& Dstar);

arma::mat d_sigma_SD_cpp(
                const arma::sp_mat& L,
                const arma::mat& SD_IplusRho,
                const arma::sp_mat& In,
                const arma::sp_mat& A);

arma::mat d_sigma_omega_corinput_cpp(
                const arma::sp_mat& L,
                const arma::mat& delta_IminOinv,
                const arma::sp_mat& A,
                const arma::mat&delta,
                const arma::sp_mat& Dstar,
                const arma::mat& IminOinv,
                const arma::sp_mat& In);

arma::mat d_sigma0_sigma_zeta_var1_cpp(
                const arma::sp_mat& L,
                const arma::mat& BetaStar,
                const arma::sp_mat& D2);


arma::mat d_phi_theta_varcov_group_cpp(
    const Rcpp::List& grouplist
);


arma::mat d_phi_theta_varcov_group_cpp(
    const Rcpp::List& grouplist
);

arma::mat d_phi_theta_varcov_cpp(
    const Rcpp::List& prep
);

#endif
