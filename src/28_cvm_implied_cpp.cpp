// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Implied model for the composite variable model (cvm). C++ twin of
// R/28_cvm_implied.R; see that file for a description of the model.

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"
#include "04_generalfit_optimWorkspace.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Solve a (small) linear system, falling back to the pseudo-inverse:
static arma::vec cvm_linsolve(const arma::mat& G, const arma::vec& rhs){
  arma::vec res;
  bool ok = arma::solve(res, G, rhs, arma::solve_opts::no_approx);
  if (!ok){
    res = arma::pinv(G) * rhs;
  }
  return res;
}

// Core implementation that takes pre-formed model matrices
Rcpp::List implied_cvm_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  // Read constant data from cached workspace:
  const OptimWorkspace& ws = getOrBuildWorkspace(model);
  const Rcpp::List& types = ws.types;
  const Rcpp::List& extramatrices = ws.extramatrices;
  std::string latent = types["latent"];
  std::string residual = types["residual"];
  std::string composite = types["composite"];

  // Which blocks are in the model?
  bool has_res = extramatrices.containsElementNamed("E_res");
  bool has_comp = extramatrices.containsElementNamed("E_comp");

  // Which constructs are composites?
  arma::vec compvec = arma::vectorise(Rcpp::as<arma::mat>(extramatrices["compvec"]));
  arma::uvec compIdx = arma::find(compvec == 1);
  arma::uvec facIdx = arma::find(compvec != 1);
  int nCompLat = compIdx.n_elem;

  // Selection matrices:
  arma::mat Ec;
  arma::mat Ef;
  if (has_comp){
    Ec = arma::mat(Rcpp::as<arma::sp_mat>(extramatrices["E_comp"]));
  }
  if (has_res){
    Ef = arma::mat(Rcpp::as<arma::sp_mat>(extramatrices["E_res"]));
  }

  // Implied covariance structures of the residual and composite indicator blocks:
  if (has_res){
    x = impliedcovstructures_cpp(x, "epsilon", residual, all);
  }
  if (has_comp){
    x = impliedcovstructures_cpp(x, "c", composite, all);
  }

  int nGroup = x.length();
  int g, p;

  // Whether the composite part is proper, per group:
  std::vector<bool> properComp(nGroup, true);

  // Step 1: composite variances and the computed latent entries
  for (g=0; g<nGroup; g++){
    Rcpp::List grouplist = x[g];

    arma::mat W = grouplist["weights"];
    arma::mat B = grouplist["beta"];
    int nLat = W.n_cols;
    int nVar = W.n_rows;

    // Composite indicator covariances (T), embedded in a nVar x nVar matrix:
    arma::mat Tfull(nVar, nVar, arma::fill::zeros);
    arma::mat TW(nVar, nLat, arma::fill::zeros);
    arma::vec d(nLat, arma::fill::zeros);
    arma::vec dinv(nLat, arma::fill::zeros);

    if (has_comp){
      arma::mat sigma_c = grouplist["sigma_c"];
      Tfull = Ec * sigma_c * Ec.t();
      TW = Tfull * W;
      for (p=0; p<nCompLat; p++){
        int k = compIdx(p);
        d(k) = arma::accu(W.col(k) % TW.col(k));
        double val = 1.0 / d(k);
        if (std::isfinite(val)){
          dinv(k) = val;
        }
      }
    }

    // (I - B)^-1:
    arma::mat I = eye(nLat, nLat);
    arma::mat BetaStar;
    if (B.is_zero()){
      BetaStar = I;
    } else {
      BetaStar = inv(I - B);
    }

    // Computed variances / scalings:
    bool proper = true;

    if (latent == "cov"){
      arma::mat S0 = grouplist["sigma_zeta"];
      if (nCompLat > 0){
        for (p=0; p<nCompLat; p++){
          S0(compIdx(p), compIdx(p)) = 0;
        }
        arma::mat G = arma::square(BetaStar.submat(compIdx, compIdx));
        arma::mat BSB = BetaStar * S0 * BetaStar.t();
        arma::vec BSBdiag = BSB.diag();
        arma::vec rhs = d.elem(compIdx) - BSBdiag.elem(compIdx);
        arma::vec psi = cvm_linsolve(G, rhs);
        // The computed residual variances may become negative (explained
        // variance larger than the composite variance): the model is continued
        // smoothly, as for a free variance parameter (see R/28_cvm_implied.R).
        for (p=0; p<nCompLat; p++){
          S0(compIdx(p), compIdx(p)) = psi(p);
        }
      }
      grouplist["sigma_zeta"] = S0;

    } else if (latent == "ggm"){
      arma::mat omega = grouplist["omega_zeta"];
      arma::mat deltamat = grouplist["delta_zeta"];
      arma::vec delta = deltamat.diag();

      // Residual variances of composites that are isolated in the network;
      // these are placed directly in sigma_zeta after it is formed (NA: not used):
      arma::mat varpatch(nLat, 1);
      varpatch.fill(NA_REAL);

      if (nCompLat > 0){
        arma::mat IminOinv = solve_symmetric_cpp_matrixonly(I - omega);

        // Composites that are isolated in the network get their residual
        // variance from the linear system; all other composites must be
        // exogenous, in which case the scaling follows in closed form:
        std::vector<arma::uword> isoVec;
        double eps = std::sqrt(std::numeric_limits<double>::epsilon());

        for (p=0; p<nCompLat; p++){
          int k = compIdx(p);
          delta(k) = 0;
          bool isolated = !arma::any(omega.row(k) != 0);
          if (isolated){
            isoVec.push_back(k);
          } else {
            // Must be exogenous: no incoming effects (tested on beta itself,
            // which is exact):
            if (arma::any(B.row(k) != 0)){
              Rcpp::stop("A composite cannot have both incoming effects (beta) and edges in the latent network (omega_zeta).");
            }
          }
        }
        for (p=0; p<nCompLat; p++){
          int k = compIdx(p);
          bool isolated = !arma::any(omega.row(k) != 0);
          if (!isolated){
            double val = d(k) / IminOinv(k,k);
            if (!std::isfinite(val) || val <= 0){
              // Only possible when a covariance matrix in the model is not
              // positive definite:
              proper = false;
              val = eps;
            }
            delta(k) = std::sqrt(val);
          }
        }
        if (isoVec.size() > 0){
          arma::uvec isoIdx(isoVec);
          arma::mat S0 = arma::diagmat(delta) * IminOinv * arma::diagmat(delta);
          arma::mat G = arma::square(BetaStar.submat(isoIdx, isoIdx));
          arma::mat BSB = BetaStar * S0 * BetaStar.t();
          arma::vec BSBdiag = BSB.diag();
          arma::vec rhs = d.elem(isoIdx) - BSBdiag.elem(isoIdx);
          arma::vec psi = cvm_linsolve(G, rhs);
          for (p=0; p<(int)isoIdx.n_elem; p++){
            double val = psi(p);
            varpatch(isoIdx(p), 0) = val;
            // The scaling is only defined for a positive variance:
            if (!std::isfinite(val) || val <= 0){
              val = eps;
            }
            delta(isoIdx(p)) = std::sqrt(val);
          }
        }
      }

      arma::mat newdelta = arma::diagmat(delta);
      grouplist["delta_zeta"] = newdelta;
      grouplist["compvar_patch"] = varpatch;

    } else {
      Rcpp::stop("Latent type is not supported in the cvm framework.");
    }

    properComp[g] = proper;

    // Store:
    grouplist["Tfull"] = Tfull;
    grouplist["TW"] = TW;
    arma::mat dmat = d;
    arma::mat dinvmat = dinv;
    grouplist["dvec"] = dmat;
    grouplist["dinv"] = dinvmat;
    grouplist["BetaStar"] = BetaStar;

    x[g] = grouplist;
  }

  // Implied latent covariance structure (now that the computed entries are in):
  x = impliedcovstructures_cpp(x, "zeta", latent, all);

  // Step 2: implied means and covariances
  for (g=0; g<nGroup; g++){
    bool proper = true;

    Rcpp::List grouplist = x[g];

    arma::mat W = grouplist["weights"];
    arma::mat B = grouplist["beta"];
    int nLat = W.n_cols;
    arma::mat Tfull = grouplist["Tfull"];
    arma::mat TW = grouplist["TW"];
    arma::vec dinv = arma::vectorise(Rcpp::as<arma::mat>(grouplist["dinv"]));
    arma::mat BetaStar = grouplist["BetaStar"];
    arma::mat sigma_zeta = grouplist["sigma_zeta"];
    arma::mat I = eye(nLat, nLat);

    // latent = "ggm": the residual variance of a composite that is isolated in
    // the network is placed in sigma_zeta directly (identical to the squared
    // scaling when positive, and a smooth continuation otherwise):
    if (grouplist.containsElementNamed("compvar_patch")){
      arma::mat patch = grouplist["compvar_patch"];
      bool patched = false;
      for (p=0; p<(int)patch.n_rows; p++){
        if (!R_IsNA(patch(p,0)) && !std::isnan(patch(p,0))){
          sigma_zeta(p,p) = patch(p,0);
          patched = true;
        }
      }
      if (patched){
        grouplist["sigma_zeta"] = sigma_zeta;
      }
    }

    // Composite loadings:
    arma::mat Q = TW;
    Q.each_row() %= dinv.t();

    // Factor loadings only (composite columns are computed):
    arma::mat lambda_f = grouplist["lambda"];
    for (p=0; p<nCompLat; p++){
      lambda_f.col(compIdx(p)).zeros();
    }

    // Full loading matrix:
    arma::mat lambda = lambda_f + Q;

    arma::mat Lambda_BetaStar = lambda * BetaStar;
    arma::mat Betasta_sigmaZeta = BetaStar * sigma_zeta;
    arma::mat sigma_eta = Betasta_sigmaZeta * BetaStar.t();

    // Residual part:
    arma::mat residComp = Tfull - Q * TW.t();
    arma::mat residPart = residComp;
    if (has_res){
      arma::mat sigma_epsilon = grouplist["sigma_epsilon"];
      residPart += Ef * sigma_epsilon * Ef.t();
    }

    // Implied covariances:
    arma::mat sigma = Lambda_BetaStar * sigma_zeta * Lambda_BetaStar.t() + residPart;
    sigma = 0.5 * (sigma + sigma.t());

    // Mean structure:
    bool hasmeans = grouplist.containsElementNamed("nu");
    arma::mat Af(nLat, nLat, arma::fill::zeros);
    arma::mat Mc;
    arma::mat mu_eta;

    if (hasmeans){
      arma::mat nu = grouplist["nu"];
      arma::mat nu_eta = grouplist["nu_eta"];
      if (facIdx.n_elem > 0){
        arma::mat If = eye(facIdx.n_elem, facIdx.n_elem);
        Af.submat(facIdx, facIdx) = inv(If - B.submat(facIdx, facIdx));
      }
      Mc = I + Af * B;
      for (p=0; p<(int)facIdx.n_elem; p++){
        Mc.col(facIdx(p)).zeros();
      }
      mu_eta = Af * nu_eta + Mc * (W.t() * nu);
      arma::vec mu = arma::vectorise(nu + lambda_f * mu_eta);
      grouplist["mu"] = mu;
      grouplist["mu_eta"] = mu_eta;
    }

    // Store:
    grouplist["lambda"] = lambda;
    grouplist["sigma"] = sigma;
    arma::mat kappa = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
    grouplist["kappa"] = kappa;
    grouplist["sigma_eta"] = sigma_eta;

    if (!all){
      // Extra matrices needed for the derivatives:
      grouplist["lambda_f"] = lambda_f;
      grouplist["Lambda_BetaStar"] = Lambda_BetaStar;
      grouplist["Betasta_sigmaZeta"] = Betasta_sigmaZeta;
      arma::mat tBetakronBeta;
      if (B.is_zero()){
        tBetakronBeta = eye(nLat*nLat, nLat*nLat);
      } else {
        tBetakronBeta = kron(BetaStar.t(), BetaStar);
      }
      grouplist["tBetakronBeta"] = tBetakronBeta;
      if (hasmeans){
        grouplist["Af"] = Af;
        grouplist["Mc"] = Mc;
      }
    } else {
      // Computed composite intercepts:
      if (hasmeans){
        arma::mat nu_eta_full = (I - B) * mu_eta;
        grouplist["nu_eta"] = nu_eta_full;
      }
      // Implied (singular) residual covariances of the composite indicators:
      if (has_comp){
        arma::mat sigma_epsilon_c = Ec.t() * residComp * Ec;
        grouplist["sigma_epsilon_c"] = sigma_epsilon_c;
      }
      // Remove the internal matrices:
      Rcpp::CharacterVector nms = grouplist.names();
      Rcpp::List cleaned;
      for (int i=0; i<grouplist.length(); i++){
        std::string nm = Rcpp::as<std::string>(nms[i]);
        if (nm != "Tfull" && nm != "TW" && nm != "dvec" && nm != "dinv" && nm != "BetaStar" && nm != "compvar_patch"){
          cleaned[nm] = grouplist[i];
        }
      }
      grouplist = cleaned;
    }

    // Return properness:
    grouplist["proper"] = proper && properComp[g];

    x[g] = grouplist;
  }

  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_cvm_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_cvm_cpp_core(x, model, all);
}
