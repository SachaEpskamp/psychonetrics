// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "14_varcov_derivatives_cpp.h"
#include "15_lvm_derivatives_cpp.h"
#include "02_algebrahelpers_RcppHelpers.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Group-level Jacobian for meta_lvm
// Combines:
//   - Mean part: LVM Jacobian of vech(sigma_y) w.r.t. LVM parameters
//   - Variance part: random effects Jacobian (same as meta_varcov)
// [[Rcpp::export]]
arma::mat d_phi_theta_meta_lvm_group_cpp(
    const Rcpp::List& grouplist
){
  // Read model matrices:
  arma::mat lambda = grouplist["lambda"];
  std::string latent = grouplist["latent"];
  std::string residual = grouplist["residual"];
  std::string randomEffects = grouplist["randomEffects"];

  // Number of variables:
  int nvar = lambda.n_rows;

  // Number of latents:
  int nlat = lambda.n_cols;

  // Number of modeled elements in the mean (vech with diagonal for covariance input):
  int nmod = nvar * (nvar + 1) / 2;

  // Number of observations: mean part + variance part
  int nobs = nmod + nmod * (nmod + 1) / 2;

  // LVM parameter count (no mean structure, no thresholds):
  // lambda (nvar*nlat) + beta (nlat^2) + sigma_zeta (nlat*(nlat+1)/2) + sigma_epsilon (nvar*(nvar+1)/2)
  int nLVM = nvar * nlat + nlat * nlat + nlat * (nlat + 1) / 2 + nvar * (nvar + 1) / 2;

  // Random effects parameter count:
  int nRan;
  if (randomEffects == "cov" || randomEffects == "chol" || randomEffects == "prec"){
    nRan = nmod * (nmod + 1) / 2;
  } else {
    // "ggm" or "cor": omega/rho (off-diag) + delta/SD (diag)
    nRan = nmod * (nmod - 1) / 2 + nmod;
  }

  // Total parameters:
  int nTotal = nLVM + nRan;

  // Index boundaries:
  int meanPart_start = 0;
  int meanPart_end = nmod - 1;
  int varPart_start = nmod;
  int varPart_end = nobs - 1;

  int ranInds_start = nLVM;
  int ranInds_end = nTotal - 1;

  // Empty Jacobian:
  arma::mat Jac = zeros(nobs, nTotal);

  // ============================================================
  // MEAN PART: LVM Jacobian of vech(sigma_y) w.r.t. LVM params
  // ============================================================
  // Built inline using individual LVM derivative helper functions
  // (avoids naming conflict with meta-level sigma)
  // Structure mirrors d_phi_theta_lvm_group_cpp sigma part (lines 318-385)
  // but without mean structure and without corinput row-stripping.

  // Parameter column indices within LVM block:
  int curInd = 0;

  int lambdaInds_start = curInd;
  int lambdaInds_end = curInd + nlat * nvar - 1;
  curInd = lambdaInds_end + 1;

  int betaInds_start = curInd;
  int betaInds_end = curInd + nlat * nlat - 1;
  curInd = betaInds_end + 1;

  int sigmazetaInds_start = curInd;
  int sigmazetaInds_end = curInd + nlat * (nlat + 1) / 2 - 1;
  curInd = sigmazetaInds_end + 1;

  int sigmaepsilonInds_start = curInd;
  int sigmaepsilonInds_end = curInd + nvar * (nvar + 1) / 2 - 1;

  // Fill factor loading part (sigma part):
  Jac.submat(meanPart_start, lambdaInds_start, meanPart_end, lambdaInds_end) = d_sigma_lambda_lvm_cpp(
    grouplist["L"], grouplist["Lambda_BetaStar"], grouplist["Betasta_sigmaZeta"], grouplist["In"], grouplist["C"]
  );

  // Fill beta part (sigma part):
  Jac.submat(meanPart_start, betaInds_start, meanPart_end, betaInds_end) = d_sigma_beta_lvm_cpp(
    grouplist["L"], grouplist["lambda"], grouplist["Betasta_sigmaZeta"], grouplist["Cbeta"], grouplist["Inlatent"], grouplist["tBetakronBeta"]
  );

  // Fill latent variances part:
  Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) = d_sigma_sigma_zeta_lvm_cpp(
    grouplist["L"], grouplist["Lambda_BetaStar"], grouplist["Deta"]
  );

  // Chain rule for latent variance parameterization:
  if (latent == "chol"){
    Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) =
      Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) *
      d_sigma_zeta_cholesky_lvm_cpp(
        grouplist["lowertri_zeta"], grouplist["L_eta"], grouplist["Cbeta"], grouplist["Inlatent"]
      );

  } else if (latent == "prec"){
    Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) =
      Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) *
      d_sigma_zeta_kappa_lvm_cpp(
        grouplist["L_eta"], grouplist["Deta"], grouplist["sigma_zeta"]
      );

  } else if (latent == "ggm"){
    Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) =
      Jac.submat(meanPart_start, sigmazetaInds_start, meanPart_end, sigmazetaInds_end) *
      d_sigma_zeta_ggm_lvm_cpp(
        grouplist["L_eta"], grouplist["delta_IminOinv_zeta"], grouplist["Aeta"],
        grouplist["delta_zeta"], grouplist["Dstar_eta"], grouplist["Inlatent"]
      );
  }

  // Fill residual variances part:
  if (residual == "cov"){
    Jac.submat(meanPart_start, sigmaepsilonInds_start, meanPart_end, sigmaepsilonInds_end) =
      eye(nmod, nvar * (nvar + 1) / 2);

  } else if (residual == "chol"){
    Jac.submat(meanPart_start, sigmaepsilonInds_start, meanPart_end, sigmaepsilonInds_end) =
      d_sigma_epsilon_cholesky_lvm_cpp(
        grouplist["lowertri_epsilon"], grouplist["L"], grouplist["C_chol"], grouplist["In"]
      );

  } else if (residual == "prec"){
    Jac.submat(meanPart_start, sigmaepsilonInds_start, meanPart_end, sigmaepsilonInds_end) =
      d_sigma_epsilon_kappa_lvm_cpp(
        grouplist["L"], grouplist["D"], grouplist["sigma_epsilon"]
      );

  } else if (residual == "ggm"){
    Jac.submat(meanPart_start, sigmaepsilonInds_start, meanPart_end, sigmaepsilonInds_end) =
      d_sigma_epsilon_ggm_lvm_cpp(
        grouplist["L"], grouplist["delta_IminOinv_epsilon"], grouplist["A"],
        grouplist["delta_epsilon"], grouplist["Dstar"], grouplist["In"]
      );
  }

  // ============================================================
  // VARIANCE PART: random effects Jacobian
  // (same structure as meta_varcov derivatives)
  // ============================================================
  int nEl = nmod * (nmod + 1) / 2;

  if (randomEffects == "cov"){
    Jac.submat(varPart_start, ranInds_start, varPart_end, ranInds_end) = eye(nEl, nEl);

  } else if (randomEffects == "chol"){
    Jac.submat(varPart_start, ranInds_start, varPart_end, ranInds_end) = d_sigma_cholesky_cpp(
      grouplist["lowertri_randomEffects"], grouplist["L_c"], grouplist["C_c"], grouplist["In_c"]
    );

  } else if (randomEffects == "ggm"){
    int netPart_start = nLVM;
    int netPart_end = nLVM + nmod * (nmod - 1) / 2 - 1;
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
    int corPart_start = nLVM;
    int corPart_end = nLVM + nmod * (nmod - 1) / 2 - 1;
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
arma::mat d_phi_theta_meta_lvm_cpp(
    const Rcpp::List& prep
){

  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  Rcpp::List groupgradients(nGroup);

  for (int i=0; i<nGroup; i++){
    arma::mat groupgrad = d_phi_theta_meta_lvm_group_cpp(groupmodels[i]);
    groupgradients[i] = groupgrad;
  }

  arma::mat res = bdiag_psychonetrics(groupgradients);

  return(res);
}
