// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Helpers for the PDC (partial directed correlations) parameterization of
// temporal effects. See R/03_modelformation_PDC.R for the math; the PDC
// matrix encodes from = row, to = column (transposed orientation relative
// to beta):
//   beta[j,i] = PDC[i,j] * sqrt(kappa_ii * sigma_jj) / sqrt(1 - PDC[i,j]^2)
// with sigma the innovation covariance and kappa = sigma^{-1}.

#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::mat PDC_to_beta_cpp(
    const arma::mat& PDC,
    const arma::mat& sigma
){
  int n = sigma.n_rows;
  bool proper = true;
  arma::mat K = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
  arma::vec sdiag = sigma.diag();
  arma::vec kdiag = K.diag();

  arma::mat B(n, n, fill::zeros);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      double p = PDC(i, j);
      double denom = 1.0 - p * p;
      if (denom < 1e-14) denom = 1e-14;
      // beta[j,i]:
      B(j, i) = p * std::sqrt(kdiag(i) * sdiag(j)) / std::sqrt(denom);
    }
  }
  return B;
}

void PDC_reparam_cpp(
    const arma::mat& PDC,
    const arma::mat& beta,
    const arma::mat& sigma,
    const arma::mat& aug,
    const arma::sp_mat& D,
    const arma::sp_mat& C,
    arma::mat& Tmat,
    arma::mat& Xmat
){
  int n = sigma.n_rows;
  bool proper = true;
  arma::mat K = solve_symmetric_cpp_matrixonly_withcheck(sigma, proper);
  arma::vec sdiag = sigma.diag();
  arma::vec kdiag = K.diag();
  arma::mat Dm = (arma::mat)D;
  int nvech = Dm.n_cols;

  // T = C * diag( sqrt(kappa_ii sigma_jj) * (1 - PDC[i,j]^2)^(-3/2) ),
  // diagonal in vec(PDC) (column-major [i,j]) order:
  arma::vec tdiag(n * n);
  for (int j = 0; j < n; j++){
    for (int i = 0; i < n; i++){
      double p = PDC(i, j);
      double denom = 1.0 - p * p;
      if (denom < 1e-14) denom = 1e-14;
      tdiag(j * n + i) = std::sqrt(kdiag(i) * sdiag(j)) * std::pow(denom, -1.5);
    }
  }
  Tmat = (arma::mat)C * arma::diagmat(tdiag);

  // d sigma_jj / d vech(sigma): duplication-matrix rows of diagonal elements:
  arma::mat dsig(n, nvech);
  for (int j = 0; j < n; j++){
    dsig.row(j) = Dm.row(j * n + j);
  }
  // d kappa_ii / d vech(sigma) = -(K[i,] kron K[i,]) * D:
  arma::mat dkap(n, nvech);
  for (int i = 0; i < n; i++){
    arma::rowvec kk(n * n);
    for (int a = 0; a < n; a++){
      for (int b = 0; b < n; b++){
        kk(a * n + b) = K(i, a) * K(i, b);
      }
    }
    dkap.row(i) = -kk * Dm;
  }

  arma::mat X(n * n, nvech, fill::zeros);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      double b = beta(j, i);
      if (b != 0.0){
        // vec index of beta[j,i]:
        X.row(i * n + j) = b / (2.0 * sdiag(j)) * dsig.row(j) + b / (2.0 * kdiag(i)) * dkap.row(i);
      }
    }
  }
  Xmat = X * aug;
}
