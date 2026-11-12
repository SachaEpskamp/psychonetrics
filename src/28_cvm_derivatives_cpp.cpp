// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Derivatives for the composite variable model (cvm). C++ twin of
// R/28_cvm_derivatives.R; see that file for the derivation.

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

// Chain from vech(sigma) of a typed covariance block to its parameters:
static arma::mat cvm_typechain_cpp(
    const std::string& type,
    const std::string& suffix,
    const Rcpp::List& grouplist,
    const arma::sp_mat& L,
    const arma::sp_mat& D,
    const arma::sp_mat& Dstar,
    const arma::sp_mat& In,
    const arma::sp_mat& C_chol,
    const arma::sp_mat& A
){
  int n = In.n_rows;
  arma::mat res;

  if (type == "cov"){

    res = eye(n*(n+1)/2, n*(n+1)/2);

  } else if (type == "chol"){

    arma::mat lowertri = grouplist["lowertri_" + suffix];
    res = d_sigma_cholesky_cpp(lowertri, L, C_chol, In);

  } else if (type == "prec"){

    arma::mat sigma = grouplist["sigma_" + suffix];
    res = d_sigma_kappa_cpp(L, D, sigma);

  } else if (type == "ggm"){

    arma::mat delta_IminOinv = grouplist["delta_IminOinv_" + suffix];
    arma::mat delta = grouplist["delta_" + suffix];
    res = join_rows(
      d_sigma_omega_cpp(L, delta_IminOinv, A, delta, Dstar),
      d_sigma_delta_cpp(L, delta_IminOinv, In, A)
    );

  } else if (type == "cor"){

    arma::mat SD = grouplist["SD_" + suffix];
    arma::mat SD_IplusRho = grouplist["SD_IplusRho_" + suffix];
    res = join_rows(
      d_sigma_rho_cpp(L, SD, A, Dstar),
      d_sigma_SD_cpp(L, SD_IplusRho, In, A)
    );

  } else {
    Rcpp::stop("Unknown type.");
  }

  return(res);
}

// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_cvm_group_cpp(
    const Rcpp::List& grouplist
){
  // Model matrices and types:
  arma::mat lambda = grouplist["lambda"];
  arma::mat lambda_f = grouplist["lambda_f"];
  arma::mat weights = grouplist["weights"];
  std::string latent = grouplist["latent"];
  std::string residual = grouplist["residual"];
  std::string composite = grouplist["composite"];
  bool meanstructure = true;
  if (grouplist.containsElementNamed("meanstructure")){
    meanstructure = grouplist["meanstructure"];
  }

  arma::mat BetaStar = grouplist["BetaStar"];
  arma::mat Lambda_BetaStar = grouplist["Lambda_BetaStar"];
  arma::mat Betasta_sigmaZeta = grouplist["Betasta_sigmaZeta"];
  arma::mat tBetakronBeta = grouplist["tBetakronBeta"];

  // Extra matrices:
  arma::sp_mat L = grouplist["L"];
  arma::sp_mat In = grouplist["In"];
  arma::sp_mat C = grouplist["C"];
  arma::sp_mat Cbeta = grouplist["Cbeta"];
  arma::sp_mat Inlatent = grouplist["Inlatent"];
  arma::sp_mat Deta = grouplist["Deta"];

  bool has_res = grouplist.containsElementNamed("E_res");
  bool has_comp = grouplist.containsElementNamed("E_comp");

  // Dimensions:
  int nvar = lambda.n_rows;
  int nlat = lambda.n_cols;
  int nres = 0;
  int ncmp = 0;
  arma::mat Ec;
  if (has_res){
    arma::sp_mat E_res = grouplist["E_res"];
    nres = E_res.n_cols;
  }
  if (has_comp){
    arma::sp_mat E_comp = grouplist["E_comp"];
    ncmp = E_comp.n_cols;
    Ec = arma::mat(E_comp);
  }
  int nvech = nvar * (nvar + 1) / 2;
  int nvechLat = nlat * (nlat + 1) / 2;

  // Composites:
  arma::vec compvec = arma::vectorise(Rcpp::as<arma::mat>(grouplist["compvec"]));
  arma::uvec compIdx = arma::find(compvec == 1);
  int nCompLat = compIdx.n_elem;

  arma::mat Ilat = eye(nlat, nlat);

  // Number of rows:
  int nMean = meanstructure ? nvar : 0;
  int nobs = nMean + nvech;

  // Number of columns:
  int nelement = (meanstructure ? nvar : 0) + // nu
    (meanstructure ? nlat : 0) + // nu_eta
    nvar * nlat + // lambda
    nvar * nlat + // weights
    nlat * nlat + // beta
    nvechLat + // latent block
    nres * (nres + 1) / 2 + // residual block
    ncmp * (ncmp + 1) / 2; // composite indicator block

  // Empty Jacobian:
  arma::mat Jac = zeros(nobs, nelement);

  // Row indices (first and last):
  int sigmaFirst = nMean;
  int sigmaLast = nMean + nvech - 1;

  // Column indices (first):
  int curInd = 0;
  int nuFirst = 0, nuetaFirst = 0;
  if (meanstructure){
    nuFirst = curInd;
    curInd += nvar;
    nuetaFirst = curInd;
    curInd += nlat;
  }
  int lambdaFirst = curInd;
  curInd += nvar * nlat;
  int weightsFirst = curInd;
  curInd += nvar * nlat;
  int betaFirst = curInd;
  curInd += nlat * nlat;
  int latentFirst = curInd;
  curInd += nvechLat;
  int resFirst = curInd;
  curInd += nres * (nres + 1) / 2;
  int compFirst = curInd;

  // Mean part:
  if (meanstructure){
    arma::mat nu = grouplist["nu"];
    arma::mat mu_eta = grouplist["mu_eta"];
    arma::mat Af = grouplist["Af"];
    arma::mat Mc = grouplist["Mc"];

    Jac.submat(0, nuFirst, nMean - 1, nuFirst + nvar - 1) = eye(nvar, nvar) + lambda_f * Mc * weights.t();
    Jac.submat(0, nuetaFirst, nMean - 1, nuetaFirst + nlat - 1) = lambda_f * Af;
    Jac.submat(0, lambdaFirst, nMean - 1, lambdaFirst + nvar*nlat - 1) = kron(mu_eta.t(), eye(nvar, nvar));
    Jac.submat(0, weightsFirst, nMean - 1, weightsFirst + nvar*nlat - 1) = lambda_f * Mc * kron(Ilat, nu.t());
    Jac.submat(0, betaFirst, nMean - 1, betaFirst + nlat*nlat - 1) = lambda_f * kron(mu_eta.t(), Af);
  }

  // Variance part. lvm derivatives, evaluated at the completed matrices:
  arma::mat J_lambda = d_sigma_lambda_lvm_cpp(L, Lambda_BetaStar, Betasta_sigmaZeta, In, C);
  arma::mat J_beta = d_sigma_beta_lvm_cpp(L, lambda, Betasta_sigmaZeta, Cbeta, Inlatent, tBetakronBeta);
  arma::mat J_psi = d_sigma_sigma_zeta_lvm_cpp(L, Lambda_BetaStar, Deta);

  // Chain to the latent parameterization:
  arma::mat Jtype;
  if (latent == "cov"){
    Jtype = eye(nvechLat, nvechLat);
  } else if (latent == "ggm"){
    Jtype = d_sigma_zeta_ggm_lvm_cpp(
      grouplist["L_eta"], grouplist["delta_IminOinv_zeta"], grouplist["Aeta"],
               grouplist["delta_zeta"], grouplist["Dstar_eta"], grouplist["Inlatent"]
    );
  } else {
    Rcpp::stop("Latent type not supported in the cvm framework.");
  }
  arma::mat J_latent = J_psi * Jtype;

  if (nCompLat > 0){
    arma::mat V = grouplist["sigma_eta"];
    arma::mat TW = grouplist["TW"];
    arma::mat Tfull = grouplist["Tfull"];
    arma::vec dv = arma::vectorise(Rcpp::as<arma::mat>(grouplist["dinv"]));
    arma::sp_mat Aeta = grouplist["Aeta"];
    arma::mat Q = lambda - lambda_f;
    arma::mat X = lambda * V - TW;
    int p;

    // U: derivative of vech(sigma_zeta) with respect to the computed quantities
    // (the variance for latent = "cov" and for composites that are isolated in
    // the latent network, the scaling delta otherwise):
    arma::mat U(nvechLat, nCompLat, arma::fill::zeros);
    if (latent == "cov"){
      for (p=0; p<nCompLat; p++){
        int k = compIdx(p);
        // Position of element (k,k) in vech:
        U(k * nlat - (k * (k - 1)) / 2, p) = 1;
      }
    } else {
      arma::mat omega = grouplist["omega_zeta"];
      for (p=0; p<nCompLat; p++){
        int k = compIdx(p);
        bool isolated = !arma::any(omega.row(k) != 0);
        if (isolated){
          U(k * nlat - (k * (k - 1)) / 2, p) = 1;
        } else {
          U.col(p) = Jtype.col(nlat * (nlat - 1) / 2 + k);
        }
      }
    }

    // Rv: derivative of diag_C(V) with respect to vech(sigma_zeta):
    arma::mat Rv(nCompLat, nvechLat);
    // RB: derivative of diag_C(V) with respect to vec(beta):
    arma::mat RB(nCompLat, nlat * nlat);
    arma::sp_mat tDeta = Deta.t();
    for (p=0; p<nCompLat; p++){
      int k = compIdx(p);
      arma::vec a = BetaStar.row(k).t();
      arma::mat aa = a * a.t();
      arma::vec rvrow = tDeta * arma::vectorise(aa);
      Rv.row(p) = rvrow.t();
      arma::mat av = a * V.col(k).t();
      RB.row(p) = 2.0 * arma::vectorise(av).t();
    }

    // Implicit derivative of the computed quantities:
    arma::mat HC = Rv * U;
    arma::mat HCinv;
    bool invok = arma::inv(HCinv, HC);
    if (!invok){
      HCinv = arma::pinv(HC);
    }
    arma::mat K = (J_psi * U) * HCinv;

    J_latent -= K * (Rv * Jtype);
    J_beta -= K * RB;

    // Derivative with respect to the composite loadings Q (lvm form at X):
    arma::mat J_Q = d_sigma_lambda_lvm_cpp(L, X, Ilat, In, C);

    // Derivative with respect to the composite variances d:
    arma::mat J_d(nvech, nlat, arma::fill::zeros);
    for (p=0; p<nCompLat; p++){
      int k = compIdx(p);
      arma::vec q = Q.col(k);
      arma::vec xk = X.col(k);
      arma::mat M1 = q * q.t() + dv(k) * (q * xk.t() + xk * q.t());
      J_d.col(k) = K.col(p) - vech(M1);
    }

    // Weights:
    arma::mat KdW = 2.0 * arma::mat(Aeta.t() * kron(Ilat, TW.t()));
    arma::mat J_W = J_Q * kron(arma::diagmat(dv), Tfull) + J_d * KdW;

    // Composite indicator covariances (T):
    arma::sp_mat LEE_comp = grouplist["LEE_comp"];
    arma::sp_mat D_comp = grouplist["D_comp"];
    arma::mat Wc = Ec.t() * weights;
    arma::mat WcD = Wc.t();
    WcD.each_col() %= dv;
    arma::mat J_Tvec = arma::mat(LEE_comp) + J_Q * kron(WcD, Ec) + J_d * arma::mat(Aeta.t() * kron(Wc.t(), Wc.t()));
    arma::mat J_T = J_Tvec * D_comp;

    // Chain to the composite parameterization:
    J_T = J_T * cvm_typechain_cpp(composite, "c", grouplist,
                                  grouplist["L_comp"], D_comp, grouplist["Dstar_comp"],
                                  grouplist["In_comp"], grouplist["C_chol_comp"], grouplist["A_comp"]);

    Jac.submat(sigmaFirst, weightsFirst, sigmaLast, weightsFirst + nvar*nlat - 1) = J_W;
    Jac.submat(sigmaFirst, compFirst, sigmaLast, compFirst + ncmp*(ncmp+1)/2 - 1) = J_T;
  }

  Jac.submat(sigmaFirst, lambdaFirst, sigmaLast, lambdaFirst + nvar*nlat - 1) = J_lambda;
  Jac.submat(sigmaFirst, betaFirst, sigmaLast, betaFirst + nlat*nlat - 1) = J_beta;
  Jac.submat(sigmaFirst, latentFirst, sigmaLast, latentFirst + nvechLat - 1) = J_latent;

  // Residuals of the factor indicators:
  if (nres > 0){
    arma::sp_mat Sel_res = grouplist["Sel_res"];
    arma::mat chain = cvm_typechain_cpp(residual, "epsilon", grouplist,
                                        grouplist["L_res"], grouplist["D_res"], grouplist["Dstar_res"],
                                        grouplist["In_res"], grouplist["C_chol_res"], grouplist["A_res"]);
    Jac.submat(sigmaFirst, resFirst, sigmaLast, resFirst + nres*(nres+1)/2 - 1) = arma::mat(Sel_res * chain);
  }

  // Return:
  return(Jac);
}

// [[Rcpp::export]]
arma::mat d_phi_theta_cvm_cpp(
    const Rcpp::List& prep
){

  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  Rcpp::List groupgradients(nGroup);

  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_cvm_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }

  arma::mat res =  bdiag_psychonetrics(groupgradients);

  return(res);
}
