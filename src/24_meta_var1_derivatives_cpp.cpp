// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "14_varcov_derivatives_cpp.h"
#include "16_var1_derivatives_cpp.h"
#include "02_algebrahelpers_RcppHelpers.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Group-level Jacobian for meta_var1
// Combines:
//   - Mean part: VAR(1) Jacobian of [vech(Sigma0), vec(Sigma1)] w.r.t. VAR(1) parameters
//   - Variance part: random effects Jacobian (same as meta_lvm/meta_varcov)
// [[Rcpp::export]]
arma::mat d_phi_theta_meta_var1_group_cpp(
    const Rcpp::List& grouplist
){
  // Read model matrices:
  arma::mat beta = grouplist["beta"];
  std::string zeta = grouplist["zeta"];
  std::string randomEffects = grouplist["randomEffects"];

  // Number of nodes:
  int nNode = beta.n_rows;

  // Number of modeled elements in the mean part:
  // vech(Sigma0) = nNode*(nNode+1)/2 + vec(Sigma1) = nNode^2
  int nSigma0 = nNode * (nNode + 1) / 2;
  int nSigma1 = nNode * nNode;
  int nmod = nSigma0 + nSigma1;

  // Number of observations: mean part + variance part
  int nobs = nmod +                  // mean part (vech(Sigma0), vec(Sigma1))
    nmod * (nmod + 1) / 2;           // variance part (vech of random effects covariance)

  // Mean part and variance part row indices:
  int meanPart_start = 0;
  int meanPart_end = nmod - 1;
  int varPart_start = nmod;
  int varPart_end = nobs - 1;

  // Sigma0 and Sigma1 row indices within mean part:
  int sigma0_start = 0;
  int sigma0_end = nSigma0 - 1;
  int sigma1_start = nSigma0;
  int sigma1_end = nmod - 1;

  // ---- VAR(1) parameter count ----
  // Parameters: beta (nNode^2) + sigma_zeta (depends on parameterization)
  int nBeta = nNode * nNode;
  int nSigmaZeta;

  if (zeta == "ggm"){
    nSigmaZeta = nNode * (nNode - 1) / 2 + nNode;  // omega (off-diag) + delta (diag)
  } else if (zeta == "prec"){
    nSigmaZeta = nNode * (nNode + 1) / 2;
  } else {
    // "cov" or "chol"
    nSigmaZeta = nNode * (nNode + 1) / 2;
  }
  int nVAR1 = nBeta + nSigmaZeta;

  // ---- Random effects parameter count ----
  int nRan;
  if (randomEffects == "cov" || randomEffects == "chol" || randomEffects == "prec"){
    nRan = nmod * (nmod + 1) / 2;
  } else {
    // "ggm" or "cor": omega/rho (off-diag) + delta/SD (diag)
    nRan = nmod * (nmod - 1) / 2 + nmod;
  }

  // Total parameters:
  int nTotal = nVAR1 + nRan;

  // Parameter column indices:
  int betaInds_start = 0;
  int betaInds_end = nBeta - 1;
  int sigmazetaInds_start = nBeta;
  int sigmazetaInds_end = nBeta + nSigmaZeta - 1;
  int ranInds_start = nVAR1;
  int ranInds_end = nTotal - 1;

  // Empty Jacobian:
  arma::mat Jac = zeros(nobs, nTotal);

  // ============================================================
  // MEAN PART: VAR(1) Jacobian of [vech(Sigma0), vec(Sigma1)]
  //            w.r.t. [vec(beta), vech/vec(sigma_zeta params)]
  // ============================================================

  // Read extramatrices needed for VAR(1) derivatives:
  arma::mat BetaStar = grouplist["BetaStar"];
  arma::sp_mat In = grouplist["In"];
  arma::sp_mat C = grouplist["C"];
  arma::sp_mat L = grouplist["L"];
  arma::sp_mat D2 = grouplist["D2"];
  arma::sp_mat IkronBeta = grouplist["IkronBeta"];
  arma::mat Sigma0 = grouplist["Sigma0"];
  arma::mat Sigma1 = grouplist["Sigma1"];

  // Construct dummy 2n x 2n sigma for var1 derivative helpers:
  // sigma = [0,         Sigma1']
  //         [Sigma1,    Sigma0 ]
  int n = nNode;
  arma::mat dummySigma = zeros(2*n, 2*n);
  dummySigma.submat(n, n, 2*n-1, 2*n-1) = Sigma0;
  dummySigma.submat(n, 0, 2*n-1, n-1) = Sigma1;
  dummySigma.submat(0, n, n-1, 2*n-1) = Sigma1.t();

  // d(vech(Sigma0))/d(beta):
  arma::mat Jb = d_sigma0_beta_var1_cpp(BetaStar, In, dummySigma, C, L);
  Jac.submat(sigma0_start, betaInds_start, sigma0_end, betaInds_end) = Jb;

  // d(vech(Sigma0))/d(sigma_zeta):
  Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) =
    d_sigma0_sigma_zeta_var1_cpp(L, BetaStar, D2);

  // Chain rule for sigma_zeta parameterization:
  if (zeta == "chol"){
    Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) =
      Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) *
      d_sigma_zeta_cholesky_var1_cpp(
        grouplist["lowertri_zeta"], grouplist["L"], grouplist["C"], grouplist["In"]
      );

  } else if (zeta == "prec"){
    Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) =
      Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) *
      d_sigma_zeta_kappa_var1_cpp(
        grouplist["L"], grouplist["D2"], grouplist["sigma_zeta"]
      );

  } else if (zeta == "ggm"){
    Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) =
      Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end) *
      d_sigma_zeta_ggm_var1_cpp(
        grouplist["L"], grouplist["delta_IminOinv_zeta"], grouplist["A"],
        grouplist["delta_zeta"], grouplist["Dstar"], grouplist["In"]
      );
  }

  // Store for sigma1 derivatives:
  arma::mat Js = Jac.submat(sigma0_start, sigmazetaInds_start, sigma0_end, sigmazetaInds_end);

  // d(vec(Sigma1))/d(beta):
  Jac.submat(sigma1_start, betaInds_start, sigma1_end, betaInds_end) =
    d_sigma1_beta_var1_cpp(IkronBeta, D2, Jb, dummySigma, beta, In);

  // d(vec(Sigma1))/d(sigma_zeta):
  Jac.submat(sigma1_start, sigmazetaInds_start, sigma1_end, sigmazetaInds_end) =
    d_sigma1_sigma_zeta_var1_cpp(IkronBeta, D2, Js);

  // ============================================================
  // VARIANCE PART: random effects Jacobian
  // (same structure as meta_lvm derivatives)
  // ============================================================
  int nEl = nmod * (nmod + 1) / 2;

  if (randomEffects == "cov"){
    Jac.submat(varPart_start, ranInds_start, varPart_end, ranInds_end) = eye(nEl, nEl);

  } else if (randomEffects == "chol"){
    Jac.submat(varPart_start, ranInds_start, varPart_end, ranInds_end) = d_sigma_cholesky_cpp(
      grouplist["lowertri_randomEffects"], grouplist["L_c"], grouplist["C_c"], grouplist["In_c"]
    );

  } else if (randomEffects == "ggm"){
    int netPart_start = nVAR1;
    int netPart_end = nVAR1 + nmod * (nmod - 1) / 2 - 1;
    int scalingPart_start = netPart_end + 1;
    int scalingPart_end = scalingPart_start + nmod - 1;

    Jac.submat(varPart_start, netPart_start, varPart_end, netPart_end) = d_sigma_omega_cpp(
      grouplist["L_c"], grouplist["delta_IminOinv_randomEffects"], grouplist["A_c"],
      grouplist["delta_randomEffects"], grouplist["Dstar_c"]
    );

    Jac.submat(varPart_start, scalingPart_start, varPart_end, scalingPart_end) = d_sigma_delta_cpp(
      grouplist["L_c"], grouplist["delta_IminOinv_randomEffects"], grouplist["In_c"], grouplist["A_c"]
    );

  } else if (randomEffects == "prec"){
    Jac.submat(varPart_start, ranInds_start, varPart_end, ranInds_end) = d_sigma_kappa_cpp(
      grouplist["L_c"], grouplist["D_c"], grouplist["sigma_randomEffects"]
    );

  } else if (randomEffects == "cor"){
    int corPart_start = nVAR1;
    int corPart_end = nVAR1 + nmod * (nmod - 1) / 2 - 1;
    int sdPart_start = corPart_end + 1;
    int sdPart_end = sdPart_start + nmod - 1;

    Jac.submat(varPart_start, corPart_start, varPart_end, corPart_end) = d_sigma_rho_cpp(
      grouplist["L_c"], grouplist["SD_randomEffects"], grouplist["A_c"], grouplist["Dstar_c"]
    );

    Jac.submat(varPart_start, sdPart_start, varPart_end, sdPart_end) = d_sigma_SD_cpp(
      grouplist["L_c"], grouplist["SD_IplusRho_randomEffects"], grouplist["In_c"], grouplist["A_c"]
    );
  }

  // Return:
  return(Jac);
}

// Full Jacobian across all groups:
// [[Rcpp::export]]
arma::mat d_phi_theta_meta_var1_cpp(
    const Rcpp::List& prep
){

  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  Rcpp::List groupgradients(nGroup);

  for (int i=0; i<nGroup; i++){
    arma::mat groupgrad = d_phi_theta_meta_var1_group_cpp(groupmodels[i]);
    groupgradients[i] = groupgrad;
  }

  arma::mat res = bdiag_psychonetrics(groupgradients);

  return(res);
}
