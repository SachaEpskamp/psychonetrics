// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// C++ twin of the analytic expected Hessian of the two-level Gaussian ML fit
// function with respect to phi = [mu; vech(Sigma_W); vech(Sigma_B)]. Mirrors
// R/05_MLestimator_expected_hessian_Gauss2L.R exactly; see that file for the
// formulas and the scaling convention (TWICE the per-cluster-unit expected
// information per group, so the generic 0.5 * M'J'HJM assembly yields the
// unit information).

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "05_MLestimator_expected_Hessian_Gauss2L_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// GROUP EXPECTED HESSIAN //
// [[Rcpp::export]]
arma::mat expected_hessian_Gauss2L_group_cpp(
        const Rcpp::List& grouplist
){
    arma::vec mu = grouplist["mu"];
    arma::mat SW = grouplist["sigma_within"];
    arma::mat SB = grouplist["sigma_between"];
    Rcpp::List twolevel = grouplist["twolevel"];

    // Duplication matrix for the p observed variables:
    arma::sp_mat D = grouplist["D_y"];

    int p = mu.n_elem;
    int k = p * (p + 1) / 2;

    double N = Rcpp::as<double>(twolevel["N"]);
    double J = Rcpp::as<double>(twolevel["J"]);
    Rcpp::List sizes = twolevel["sizes"];
    arma::vec nj_vec = Rcpp::as<arma::vec>(sizes["nj"]);
    arma::vec m_vec = Rcpp::as<arma::vec>(sizes["m"]);

    // Mean part and the three covariance blocks (WW, WB, BB):
    arma::mat I_mu = zeros(p, p);
    arma::mat I_WW = zeros(k, k);
    arma::mat I_WB = zeros(k, k);
    arma::mat I_BB = zeros(k, k);

    int nSizes = nj_vec.n_elem;
    for (int s = 0; s < nSizes; s++){
        double nj = nj_vec(s);
        double m = m_vec(s);

        // Omega_s^-1 = n_s * (Sigma_W + n_s Sigma_B)^-1:
        arma::mat iSj = solve_symmetric_cpp_matrixonly(SW + nj * SB);

        // Mean part: m_s * Omega_s^-1:
        I_mu += m * nj * iSj;

        // G_s = 0.5 D'(Sigma_s^-1 (x) Sigma_s^-1)D; with D_s = [(1/n_s) I, I]
        // and H_s = 0.5 D'(Omega_s^-1 (x) Omega_s^-1)D = n_s^2 G_s, the blocks
        // of D_s' H_s D_s are: WW = G_s, WB = n_s G_s, BB = n_s^2 G_s:
        arma::mat G = 0.5 * D.t() * kron(iSj, iSj) * D;

        I_WW += m * G;
        I_WB += m * nj * G;
        I_BB += m * nj * nj * G;
    }

    // Pooled-within part (zero weight when all clusters have size 1):
    if (N > J){
        arma::mat iSW = solve_symmetric_cpp_matrixonly(SW);
        I_WW += (N - J) * 0.5 * D.t() * kron(iSW, iSW) * D;
    }

    // Assemble [mean; (WW, WB; WB, BB)] block-diagonally and return TWICE the
    // per-cluster unit information (see header):
    int npar = p + 2 * k;
    arma::mat Out = zeros(npar, npar);
    Out.submat(0, 0, p - 1, p - 1) = (2.0 / J) * I_mu;
    Out.submat(p, p, p + k - 1, p + k - 1) = (2.0 / J) * I_WW;
    Out.submat(p, p + k, p + k - 1, p + 2 * k - 1) = (2.0 / J) * I_WB;
    Out.submat(p + k, p, p + 2 * k - 1, p + k - 1) = (2.0 / J) * I_WB;
    Out.submat(p + k, p + k, p + 2 * k - 1, p + 2 * k - 1) = (2.0 / J) * I_BB;

    return(Out);
}


// full expected Hessian (mirrors expected_hessian_Gauss2L)
// [[Rcpp::export]]
arma::mat expected_hessian_Gauss2L_cpp(
        const Rcpp::List& prep
){
    Rcpp::List groupmodels = prep["groupModels"];
    int nGroup = groupmodels.length();
    arma::vec nPerGroup = prep["nPerGroup"];
    double nTotal = prep["nTotal"];

    Rcpp::List grouphessians(nGroup);

    for (int i = 0; i < nGroup; i++){
        arma::mat grouphes = (nPerGroup(i) / nTotal) * expected_hessian_Gauss2L_group_cpp(groupmodels[i]);
        grouphessians[i] = grouphes;
    }

    arma::mat res = bdiag_psychonetrics(grouphessians);

    return(res);
}
