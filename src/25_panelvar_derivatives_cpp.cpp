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

// Derivatives for the panelvar model. These are the dlvm1 derivatives with
// lambda = I (see 18_dlvm1_derivatives_cpp.cpp): no factor loadings, no
// residual (co)variances, and observed means in place of latent means.

// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_panelvar_group_cpp(
    const Rcpp::List& grouplist
){
  // Setup:
  int t;

  // Things I need now:
  std::string within_latent  = grouplist["within_latent"];
  std::string between_latent  = grouplist["between_latent"];

  arma::mat design = grouplist["design"];
  arma::sp_mat P = grouplist["P"];
  arma::sp_mat L = grouplist["L"];
  arma::sp_mat In = grouplist["In"];
  arma::sp_mat D2 = grouplist["D2"];
  arma::sp_mat C = grouplist["C"];

  arma::mat BetaStar = grouplist["BetaStar"];
  Rcpp::List allSigmas_within = grouplist["allSigmas_within"];
  arma::sp_mat IkronBeta = grouplist["IkronBeta"];

  // Number of variables:
  int nVar = design.n_rows;

  // Number of time points:
  int nTime = design.n_cols;

  // I need to construct the Jacobian for the following "observations:"
  int nobs = nVar + // Means
    (nVar * (nVar+1))/2 + // Variances
    (nVar * nVar) * (nTime - 1); // lagged variances

  // total number of elements:
  int nelement = nVar + // means
    nVar * (nVar+1) / 2 + // within-person contemporaneous block
    (nVar * nVar) + // temporal effects
    nVar * (nVar + 1) / 2; // between-person block

  // Form empty matrix:
  arma::mat Jac = zeros(nobs,nelement);

  // Indices (observations):
  arma::vec meanInds = seq_len_inds(0, nVar);

  Rcpp::List sigInds(nTime);
  sigInds[0] = seq_len_inds(meanInds(1) + 1, nVar*(nVar+1)/2);

  // For each lag:
  for (t=1; t<nTime;t++){
    arma::vec prev = sigInds[t-1];
    sigInds[t]  = seq_len_inds(prev(1) + 1, nVar*nVar);
  }

  // Indices model:
  arma::vec mu_inds = seq_len_inds(0, nVar);
  arma::vec sigma_zeta_within_inds = seq_len_inds(mu_inds(1) + 1,  nVar * (nVar+1) / 2);
  arma::vec beta_inds = seq_len_inds(sigma_zeta_within_inds(1) + 1, nVar * nVar);
  arma::vec sigma_zeta_between_inds = seq_len_inds(beta_inds(1) + 1, nVar * (nVar + 1) / 2);

  // AUGMENTATION PARTS //
  arma::mat aug_within_latent;
  arma::mat aug_between_latent;

  // Within-person block:
  if (within_latent == "chol"){

    aug_within_latent = d_sigma_cholesky_cpp(
      grouplist["lowertri_zeta_within"], L, C, In);

  } else if (within_latent == "prec"){

    aug_within_latent = d_sigma_kappa_cpp(
      L, D2, grouplist["sigma_zeta_within"]);

  } else if (within_latent == "ggm"){

    aug_within_latent = join_rows(
      d_sigma_omega_cpp(L, grouplist["delta_IminOinv_zeta_within"], grouplist["A"], grouplist["delta_zeta_within"], grouplist["Dstar"]),
      d_sigma_delta_cpp(L, grouplist["delta_IminOinv_zeta_within"], In, grouplist["A"])
    );

  } else if (within_latent == "cov"){

    aug_within_latent = eye(nVar*(nVar+1)/2, nVar*(nVar+1)/2);

  } else if (within_latent == "cor"){

    aug_within_latent = join_rows(
      d_sigma_rho_cpp(L, grouplist["SD_zeta_within"], grouplist["A"], grouplist["Dstar"]),
      d_sigma_SD_cpp(L, grouplist["SD_IplusRho_zeta_within"], In, grouplist["A"])
    );

  }

  // Between-person block:
  if (between_latent == "chol"){

    aug_between_latent = d_sigma_cholesky_cpp(
      grouplist["lowertri_zeta_between"], L, C, In);

  } else if (between_latent == "prec"){

    aug_between_latent = d_sigma_kappa_cpp(
      L, D2, grouplist["sigma_zeta_between"]);

  } else if (between_latent == "ggm"){

    aug_between_latent = join_rows(
      d_sigma_omega_cpp(L, grouplist["delta_IminOinv_zeta_between"], grouplist["A"], grouplist["delta_zeta_between"], grouplist["Dstar"]),
      d_sigma_delta_cpp(L, grouplist["delta_IminOinv_zeta_between"], In, grouplist["A"])
    );

  } else if (between_latent == "cov"){

    aug_between_latent = eye(nVar*(nVar+1)/2, nVar*(nVar+1)/2);

  } else if (between_latent == "cor"){

    aug_between_latent = join_rows(
      d_sigma_rho_cpp(L, grouplist["SD_zeta_between"], grouplist["A"], grouplist["Dstar"]),
      d_sigma_SD_cpp(L, grouplist["SD_IplusRho_zeta_between"], In, grouplist["A"])
    );

  }

  // fill mean part:
  Jac.submat(meanInds(0),mu_inds(0),meanInds(1),mu_inds(1)) =  eye(nVar, nVar);

  arma::vec curSigInds = sigInds[0];

  // Fill s0 to sigma_zeta_within part (and store for later use):
  arma::mat J_sigma_zeta_within = (BetaStar * D2) * aug_within_latent;

  Jac.submat(curSigInds(0),sigma_zeta_within_inds(0),curSigInds(1),sigma_zeta_within_inds(1)) = L * J_sigma_zeta_within;

  // Fill s0 to beta part (and store for later use):
  arma::mat sigma1 = allSigmas_within[1];
  arma::mat J_sigma_beta = (eye(nVar * nVar, nVar * nVar) + C) * BetaStar * kronecker_X_I(sigma1, nVar);
  Jac.submat(curSigInds(0),beta_inds(0),curSigInds(1),beta_inds(1)) =  L * J_sigma_beta;

  // Fill s0 to sigma_zeta_between part (and store for later use):
  arma::mat J_sigma_zeta_between = D2 * aug_between_latent;

  Jac.submat(curSigInds(0),sigma_zeta_between_inds(0),curSigInds(1),sigma_zeta_between_inds(1)) =  L * J_sigma_zeta_between;

  // Now for every further lag ...
  for (t=1; t<nTime;t++){
    arma::vec curSigInds = sigInds[t];

    // Fill sk to sigma_zeta_within part (and store for later use):
    J_sigma_zeta_within = IkronBeta * J_sigma_zeta_within;
    Jac.submat(curSigInds(0),sigma_zeta_within_inds(0),curSigInds(1),sigma_zeta_within_inds(1))  = J_sigma_zeta_within;

    // Fill sk to beta part (and store for later use):
    arma::mat sigma_prev = allSigmas_within[t-1];
    J_sigma_beta = IkronBeta * J_sigma_beta + kronecker_X_I(sigma_prev.t(), nVar);
    Jac.submat(curSigInds(0),beta_inds(0),curSigInds(1),beta_inds(1)) = J_sigma_beta;

    // Fill sk to sigma_zeta_between part:
    Jac.submat(curSigInds(0),sigma_zeta_between_inds(0),curSigInds(1),sigma_zeta_between_inds(1)) = J_sigma_zeta_between;
  }

  // Permute:
  Jac = P * Jac;

  // Return:
  return(Jac);
}


// [[Rcpp::export]]
arma::mat d_phi_theta_panelvar_cpp(
    const Rcpp::List& prep
){

  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  Rcpp::List groupgradients(nGroup);

  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_panelvar_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }

  arma::mat res =  bdiag_psychonetrics(groupgradients);

  return(res);
}
