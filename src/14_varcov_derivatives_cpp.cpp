// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat d_sigma_cholesky_cpp(
    const arma::mat& lowertri,
    const arma::sp_mat& L,
    const arma::sp_mat& C,
    const arma::sp_mat In
){
  // return(L * (kron(In,In) + C) * (kron(lowertri, (arma::mat)In) * L.t()));
  
  arma::sp_mat res = L * (kronecker_diag(In) + C) * (kronecker_X_I(lowertri, In.n_rows) * L.t());
  return((arma::mat)res);
  
}

// d_sigma_cholesky <- function(lowertri,L,C,In,...){
//   
//   res <- L %*% ((In %x% In) + C) %*% ((lowertri %x% In) %*% t(L))
//   
//   as.matrix(res)
// }


// [[Rcpp::export]]
arma::mat d_sigma_delta_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv,
    const arma::sp_mat& In,
    const arma::sp_mat& A,
    const arma::sp_mat& delta
){
  // return(
  //   L * (
  //       kron(delta_IminOinv, (arma::mat)In) + 
  //         kron((arma::mat)In, delta_IminOinv) // FIXME: Kronecker product with identity matrix is much easier...
  //   ) * A
  // );
  arma::sp_mat res = L * (
    kronecker_X_I(delta_IminOinv, In.n_rows) + 
      kronecker_I_X(delta_IminOinv, In.n_rows) 
  ) * A;
  
  return((arma::mat)res);
}

// d_sigma_delta <- function(L,delta_IminOinv,In,A,delta,...){
//   res <- L %*% (
//       (delta_IminOinv%x% In) + 
//         (In %x% delta_IminOinv)
//   ) %*% A
//   
//   as.matrix(res)
// }
// 

// [[Rcpp::export]]
arma::mat d_sigma_omega_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv,
    const arma::sp_mat& A,
    const arma::sp_mat& delta,
    const arma::sp_mat& Dstar
){
  return(L * kron(delta_IminOinv, delta_IminOinv) * Dstar);
}

// 
// d_sigma_omega <- function(L,delta_IminOinv,A,delta,Dstar,...){
// # L %*% (delta %x% delta) %*% (IminOinv %x% IminOinv) %*% Dstar
//   
// # delta_IminOinv <- delta %*% IminOinv
//   res <- L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
//     
// # all(a == b)
//     as.matrix(res)
// }

// [[Rcpp::export]]
arma::mat d_sigma_kappa_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& D,
    const arma::mat& sigma){
  return(- L * kron(sigma, sigma) * D);
}

// [[Rcpp::export]]
arma::mat d_sigma_rho_cpp(
    const arma::sp_mat& L,
    const arma::sp_mat& SD,
    const arma::sp_mat& A,
    const arma::sp_mat& Dstar){
  
  // return(L * kron((arma::mat)SD, (arma::mat)SD) * Dstar); // FIXME: This kronecker prduct can be MUCH better...
  arma::sp_mat res = L * kronecker_diag(SD) * Dstar;
  return((arma::mat)res); 
}


// [[Rcpp::export]]
arma::mat d_sigma_SD_cpp(
    const arma::sp_mat& L,
    const arma::mat& SD_IplusRho,
    const arma::sp_mat& In,
    const arma::sp_mat& A){
  
  // return(L * (
  //     kron(SD_IplusRho, (arma::mat)In) + 
  //       kron((arma::mat)In, SD_IplusRho)
  // ) * A);
  
  arma::sp_mat res = L * (
    kronecker_X_I(SD_IplusRho, In.n_rows) + 
      kronecker_I_X(SD_IplusRho, In.n_rows)
  ) * A;
  
  return((arma::mat)res);
}

// [[Rcpp::export]]
arma::mat d_sigma_omega_corinput_cpp(
    const arma::sp_mat& L,
    const arma::mat& delta_IminOinv,
    const arma::sp_mat& A,
    const arma::sp_mat&delta,
    const arma::sp_mat& Dstar,
    const arma::mat& IminOinv,
    const arma::sp_mat& In){
  

  // First make the diagonal matrix needed:
  arma::vec d = diagvec(IminOinv);
  // sp_mat dmat(d.size(), d.size());
  for (int i=0; i < d.size(); i++){
    d[i] = pow(d[i], -1.5);
  }
  mat dmat = diagmat(d);

  // sp_mat dmat(d.size(), d.size());
  // for (int i=0; i < d.size(); i++){
  //   dmat(i,i) = pow(d[i], -1.5);
  // }
  
  // Dense kronecker products:
  arma::mat kron2 = kron(IminOinv, IminOinv);
  arma::mat kron1 =  kronecker_diag(delta) * kron2;
    // kron(delta_IminOinv, delta_IminOinv);

    
  // Sparse inner part:
  arma::sp_mat sparse = 0.5 *  (kronecker_X_I(delta_IminOinv, In.n_rows) + kronecker_I_X(delta_IminOinv, In.n_rows));
  arma::mat dense = A * dmat * A.t();
  
  // Return value:
  arma::mat res = 
    L * (kron1 - sparse * dense * kron2) * Dstar;

  return(res);
  
  // return(
  //   L * (
  //       kron(delta_IminOinv, delta_IminOinv) -
  //         0.5 *  (kronecker_X_I(delta_IminOinv, In.n_rows) + kronecker_I_X(delta_IminOinv, In.n_rows)) * A * 
  //         dmat * A.t() * kron(IminOinv, IminOinv)
  //   ) * Dstar
  // );
  
}