// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"
#include "14_varcov_derivatives_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Covariance-structure augmentation for one block (twin of the aug() closure in
// R/27_ml_varcov_derivatives.R). Returns a (kp x kp) matrix mapping the block's
// covariance parameters to vech(Sigma_<block>), reusing the shared single-level
// varcov d_sigma_* helpers. 'suffix' is "within" or "between"; the p-dimensional
// kernels (L, In, C, D, Dstar, A) are the same for both blocks.
static arma::mat aug_ml_varcov(
    const std::string& type,
    const std::string& suffix,
    const Rcpp::List& gl,
    int kp
){
  if (type == "cov"){
    return eye(kp, kp);
  } else if (type == "chol"){
    return d_sigma_cholesky_cpp(
      gl[std::string("lowertri_") + suffix], gl["L"], gl["C"], gl["In"]);
  } else if (type == "prec"){
    return d_sigma_kappa_cpp(
      gl["L"], gl["D"], gl[std::string("sigma_") + suffix]);
  } else if (type == "ggm"){
    return join_rows(
      d_sigma_omega_cpp(gl["L"], gl[std::string("delta_IminOinv_") + suffix],
                        gl["A"], gl[std::string("delta_") + suffix], gl["Dstar"]),
      d_sigma_delta_cpp(gl["L"], gl[std::string("delta_IminOinv_") + suffix],
                        gl["In"], gl["A"])
    );
  } else if (type == "cor"){
    return join_rows(
      d_sigma_rho_cpp(gl["L"], gl[std::string("SD_") + suffix], gl["A"], gl["Dstar"]),
      d_sigma_SD_cpp(gl["L"], gl[std::string("SD_IplusRho_") + suffix], gl["In"], gl["A"])
    );
  }
  Rcpp::stop("Unsupported ml_varcov type: " + type);
  return arma::mat();
}

// Group Jacobian of phi = [mu; vech(Sigma_within); vech(Sigma_between)] w.r.t.
// the model parameters, for the two-level sufficient-statistics ML estimator.
// C++ twin of d_phi_theta_ml_varcov_group (R): identity mean block, then the
// varcov covariance-structure Jacobian for each level. No lambda/beta/residual.
// [[Rcpp::export]]
arma::mat d_phi_theta_ml_varcov_group_cpp(
    const Rcpp::List& grouplist
){
  std::string within = grouplist["within"];
  std::string between = grouplist["between"];

  arma::mat sigma_within = grouplist["sigma_within"];
  int nVar = sigma_within.n_rows;
  int kp = nVar * (nVar + 1) / 2;

  arma::mat aug_within = aug_ml_varcov(within, "within", grouplist, kp);
  arma::mat aug_between = aug_ml_varcov(between, "between", grouplist, kp);

  int nobs = nVar + kp + kp;
  int nelement = nVar + kp + kp;
  arma::mat Jac = zeros(nobs, nelement);

  // Row/column blocks (mirroring the R index construction):
  //   rows    [ mu (p) ; vech Sigma_within (kp) ; vech Sigma_between (kp) ]
  //   columns [ mu (p) ; within params (kp)      ; between params (kp)     ]

  // Mean part: identity (mu enters phi's mean block directly):
  Jac.submat(0, 0, nVar - 1, nVar - 1) = eye(nVar, nVar);

  // Within covariance block:
  Jac.submat(nVar, nVar, nVar + kp - 1, nVar + kp - 1) = aug_within;

  // Between covariance block:
  Jac.submat(nVar + kp, nVar + kp, nVar + 2 * kp - 1, nVar + 2 * kp - 1) = aug_between;

  return(Jac);
}

// Wide-format (estimator "FIML") group Jacobian, twin of
// d_phi_theta_ml_varcov_wide_group (R). The compact rows are
//   [ mu (p) ; vech(within-position block) (kp) ; vec(between-positions block)
//     (p^2) ]
// where the within-position (diagonal) block of the wide covariance matrix is
// Sigma_within + Sigma_between and every between-positions (off-diagonal)
// block is Sigma_between. P (from the constructor's extra matrices) maps these
// compact rows to the wide distribution parameters [wide mu; vech(wide Sigma)]:
// [[Rcpp::export]]
arma::mat d_phi_theta_ml_varcov_wide_group_cpp(
    const Rcpp::List& grouplist
){
  std::string within = grouplist["within"];
  std::string between = grouplist["between"];

  arma::mat sigma_within = grouplist["sigma_within"];
  int nVar = sigma_within.n_rows;
  int kp = nVar * (nVar + 1) / 2;

  arma::mat aug_within = aug_ml_varcov(within, "within", grouplist, kp);
  arma::mat aug_between = aug_ml_varcov(between, "between", grouplist, kp);

  arma::sp_mat P = grouplist["P"];
  arma::sp_mat D = grouplist["D"];   // p-dimensional duplication matrix

  int nobs = nVar + kp + nVar * nVar;
  int nelement = nVar + kp + kp;
  arma::mat Jac = zeros(nobs, nelement);

  // Mean part:
  Jac.submat(0, 0, nVar - 1, nVar - 1) = eye(nVar, nVar);

  // Within-position (diagonal) blocks equal Sigma_within + Sigma_between, so
  // both parameter blocks contribute:
  Jac.submat(nVar, nVar, nVar + kp - 1, nVar + kp - 1) = aug_within;
  Jac.submat(nVar, nVar + kp, nVar + kp - 1, nVar + 2 * kp - 1) = aug_between;

  // Between-positions (off-diagonal) blocks equal Sigma_between; the compact
  // rows are the full vec, so map vech -> vec with the duplication matrix:
  Jac.submat(nVar + kp, nVar + kp, nVar + kp + nVar * nVar - 1, nVar + 2 * kp - 1) = D * aug_between;

  // Map the compact rows to the wide distribution parameters:
  arma::mat res = P * Jac;
  return(res);
}

// Full Jacobian across groups (block-diagonal), twin of d_phi_theta_ml_varcov.
// The prepare function flags which estimator variant applies; if the flag is
// absent, the wide variant is recognized by the P permutation matrix, which
// only FIML models carry:
// [[Rcpp::export]]
arma::mat d_phi_theta_ml_varcov_cpp(
    const Rcpp::List& prep
){
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  bool twolevel_ML;
  if (prep.containsElementNamed("twolevel_ML")){
    twolevel_ML = Rcpp::as<bool>(prep["twolevel_ML"]);
  } else {
    Rcpp::List g0 = groupmodels[0];
    twolevel_ML = !g0.containsElementNamed("P");
  }

  Rcpp::List groupgradients(nGroup);
  for (int i = 0; i < nGroup; i++){
    if (twolevel_ML){
      groupgradients[i] = d_phi_theta_ml_varcov_group_cpp(groupmodels[i]);
    } else {
      groupgradients[i] = d_phi_theta_ml_varcov_wide_group_cpp(groupmodels[i]);
    }
  }

  arma::mat res = bdiag_psychonetrics(groupgradients);
  return(res);
}
