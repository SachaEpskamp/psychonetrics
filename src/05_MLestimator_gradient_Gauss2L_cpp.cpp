// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// C++ twin of the gradient of the two-level Gaussian ML fit function with
// respect to the distribution parameters phi = [mu; vech(Sigma_W);
// vech(Sigma_B)]. Mirrors R/05_MLestimator_gradient_Gauss2L.R exactly; see
// that file for the formulas.

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "05_MLestimator_gradient_Gauss2L_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Helper: matrix derivative -> vech gradient (double off-diagonal entries,
// keep the diagonal; twin of vech_grad_2L):
static arma::vec vech_grad_2L_cpp(const arma::mat& G){
    int p = G.n_rows;
    arma::vec out(p * (p + 1) / 2);
    int ind = 0;
    // Column-major lower triangle including diagonal, matching
    // G[lower.tri(G, diag = TRUE)] in R:
    for (int j = 0; j < p; j++){
        for (int i = j; i < p; i++){
            out(ind) = (i == j) ? G(i, j) : 2.0 * G(i, j);
            ind++;
        }
    }
    return(out);
}


// GROUP JACOBIAN FUNCTION: (1/J) d(-2 l*)/d phi as 1 x npar row vector //
// [[Rcpp::export]]
arma::mat jacobian_gaussian2L_sigma_group_cpp(
        const Rcpp::List& grouplist
){
    arma::vec mu = grouplist["mu"];
    arma::mat SW = grouplist["sigma_within"];
    arma::mat SB = grouplist["sigma_between"];
    Rcpp::List twolevel = grouplist["twolevel"];

    int p = mu.n_elem;

    double N = Rcpp::as<double>(twolevel["N"]);
    double J = Rcpp::as<double>(twolevel["J"]);
    arma::mat S_PW = twolevel["S_PW"];
    Rcpp::List sizes = twolevel["sizes"];
    arma::vec nj_vec = Rcpp::as<arma::vec>(sizes["nj"]);
    arma::vec m_vec = Rcpp::as<arma::vec>(sizes["m"]);
    Rcpp::List mean_d = twolevel["mean_d"];
    Rcpp::List cov_d = twolevel["cov_d"];

    arma::vec g_mu = zeros(p);
    arma::mat G_SW = zeros(p, p);
    arma::mat G_SB = zeros(p, p);

    int nSizes = nj_vec.n_elem;
    for (int s = 0; s < nSizes; s++){
        double nj = nj_vec(s);
        double m = m_vec(s);
        arma::mat iSj = solve_symmetric_cpp_matrixonly(SW + nj * SB);
        arma::vec yc = Rcpp::as<arma::vec>(mean_d[s]) - mu;
        arma::mat A = Rcpp::as<arma::mat>(cov_d[s]) + yc * yc.t();
        arma::mat K = iSj - nj * iSj * A * iSj;
        g_mu -= m * 2.0 * nj * (iSj * yc);
        G_SW += m * K;
        G_SB += m * nj * K;
    }

    if (N > J){
        arma::mat iSW = solve_symmetric_cpp_matrixonly(SW);
        G_SW += (N - J) * (iSW - iSW * S_PW * iSW);
    }

    // Bind [mu; vech SW; vech SB] and divide by J:
    arma::vec grad = join_cols(g_mu, vech_grad_2L_cpp(G_SW), vech_grad_2L_cpp(G_SB)) / J;

    // Return as 1 x npar row matrix:
    return(grad.t());
}


// full Jacobian function (mirrors jacobian_gaussian2L_sigma)
// [[Rcpp::export]]
arma::mat jacobian_gaussian2L_sigma_cpp(
        const Rcpp::List& prep
){
    Rcpp::List groupmodels = prep["groupModels"];
    int nGroup = groupmodels.length();
    arma::vec nPerGroup = prep["nPerGroup"];
    double nTotal = prep["nTotal"];

    // Jacobian:
    Rcpp::List groupgradients(nGroup);

    for (int i = 0; i < nGroup; i++){
        arma::mat groupgrad = (nPerGroup(i) / nTotal) * jacobian_gaussian2L_sigma_group_cpp(groupmodels[i]);
        groupgradients[i] = groupgrad;
    }

    arma::mat res = cbind_psychonetrics(groupgradients);

    return(res);
}
