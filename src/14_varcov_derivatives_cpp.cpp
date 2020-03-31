// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

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
  arma::sp_mat res = L * (kronecker_diag_sparse(In) + C) * (kronecker_X_I(lowertri, In.n_rows) * L.t());
  
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
    const arma::sp_mat& A
){
  arma::sp_mat inner =  kronecker_X_I(delta_IminOinv, In.n_rows) + kronecker_I_X(delta_IminOinv, In.n_rows);
  
  
  arma::sp_mat res = L * inner * A;
  
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
    const arma::mat& delta,
    const arma::sp_mat& Dstar
){
  
  arma::mat kronmat = kron(delta_IminOinv, delta_IminOinv);
  
  
  
  return(L * kronmat * Dstar);
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
    const arma::mat& SD,
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
    const arma::mat& delta,
    const arma::sp_mat& Dstar,
    const arma::mat& IminOinv,
    const arma::sp_mat& In){
  
  
  // First make the diagonal matrix needed:
  arma::vec d = diagvec(IminOinv);
  // sp_mat dmat(d.size(), d.size());
  int dsize = d.size();
  for (int i=0; i < dsize; i++){
    d[i] = pow(d[i], -1.5);
  }
  mat dmat = diagmat(d);
  
  // sp_mat dmat(d.size(), d.size());
  // for (int i=0; i < d.size(); i++){
  //   dmat(i,i) = pow(d[i], -1.5);
  // }
  
  // Dense kronecker products:
  arma::mat kron2 = kron(IminOinv, IminOinv);
  arma::mat kron1 =  (arma::mat)kronecker_diag(delta) * kron2;
  // kron(delta_IminOinv, delta_IminOinv);
  
  
  // Sparse inner part:
  arma::sp_mat sparse = 0.5 *  (kronecker_X_I(delta_IminOinv, In.n_rows) + kronecker_I_X(delta_IminOinv, In.n_rows));
  arma::mat dense = A * dmat * A.t();
  
  // Return value:
  arma::mat res = 
    L * (kron1 - (arma::mat)sparse * dense * kron2) * Dstar;
  
  return(res);
  
  // return(
  //   L * (
  //       kron(delta_IminOinv, delta_IminOinv) -
  //         0.5 *  (kronecker_X_I(delta_IminOinv, In.n_rows) + kronecker_I_X(delta_IminOinv, In.n_rows)) * A * 
  //         dmat * A.t() * kron(IminOinv, IminOinv)
  //   ) * Dstar
  // );
  
}


// [[Rcpp::export]]
arma::mat d_sigma0_sigma_zeta_var1_cpp(
    const arma::sp_mat& L,
    const arma::mat& BetaStar,
    const arma::sp_mat& D2){
  return(L * BetaStar * D2);
}


// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_varcov_group_cpp(
    const Rcpp::List& grouplist
){
  
  // objects needed now:
  arma::mat sigma = grouplist["sigma"];
  std::string y = grouplist["y"];
  bool corinput = grouplist["corinput"];
  bool meanstructure = grouplist["meanstructure"];
  arma::mat mu = grouplist["mu"];
  
  
  
  // Number of variables:
  int nvar = sigma.n_rows;
  
  // ints needed:
  int i, j;
  
  // If not missing tau:
  arma::mat tau = grouplist["tau"];
  
  // Check if all NA:
  bool noThresholds = true;
  
  int taunrows = tau.n_rows;
  int tauncols = tau.n_cols;
  for (i = 0; i < taunrows && noThresholds; i++){
    
    for (j = 0; j < tauncols && noThresholds; j++){
      
      if (is_finite(tau(i,j))){
        noThresholds = false;
      }
    }  
  }
  
  // if (grouplist.containsElementNamed("tau")){
  //   tau = (arma::mat)grouplist["tau"];
  //   } else {
  //     tau = zeros(1, nvar);
  //   for (i=0;i<nvar; i++){
  //     tau(1,i) = NA_REAL;
  //   }
  // }
  
  // Count how many means and thresholds:
  int nThresh = 0;
  int nMean = 0;
  int nrow_tau = tau.n_rows;
  
  
  for (i=0;i<nvar;i++){
    
    if (arma::is_finite(mu(i,0))){
      nMean++;
    }
    
    for (j = 0; j < nrow_tau; j++){
      if (arma::is_finite(tau(j,i))){
        nThresh++;
      }
    }
  }
  
  int nMean_Thresh = nMean + nThresh;
  
  // number of vars:
  int nvarcov = nvar * (nvar + 1) / 2;
  
  // Number of observations:
  int nobs = nMean_Thresh + (nvar * (nvar + 1)) / 2;
  
  // Count number of parametesr:
  int npars = 0;
  if (meanstructure){
    npars = nobs - corinput * nvar;
  } else {
    npars = nobs - corinput * nvar - nMean;
  }
  
  // Empty JAcobian:
  arma::mat Jac = zeros(nobs, npars);
  
  
  // Indixes for start:
  int meanPart_start = 0;
  int meanPart_end = nMean_Thresh - 1;
  int varPart_start = nMean_Thresh;
  int varPart_end = nMean_Thresh + nvarcov- 1;
  
  // For parameters:
  int varPartPars_start = meanstructure * nMean_Thresh + nThresh;
  int varPartPars_end = meanstructure * nMean_Thresh + nThresh+ nvarcov - 1;
  
  
  // Fill mean part with diagonal:
  if (meanstructure || nThresh > 0){
    Jac.submat(meanPart_start,meanPart_start,meanPart_end,meanPart_end ) =  eye(nMean_Thresh,nMean_Thresh);
  }
  
  // Fill the simga part:
  if (y == "cov"){
    Jac.submat(varPart_start,varPartPars_start,varPart_end,varPartPars_end ) =  eye(nvarcov,nvarcov);
    
    
  } else if (y == "chol"){
    
    Jac.submat(varPart_start,varPartPars_start,varPart_end,varPartPars_end ) =  d_sigma_cholesky_cpp( 
      grouplist["lowertri"], grouplist["L"], grouplist["C"], grouplist["In"]
    );
    
    
  } else if (y == "ggm"){
    int netPart_start = meanstructure * nMean_Thresh + nThresh;
    int netPart_end = meanstructure * nMean_Thresh + nThresh + nvar * (nvar-1) / 2 - 1;
    int scalingPart_start = netPart_end + 1;
    int scalingPart_end = scalingPart_start + nvar - 1;
    
    
    if (corinput){
      
      Jac.submat(varPart_start,netPart_start,varPart_end,netPart_end ) =  d_sigma_omega_corinput_cpp(
        grouplist["L"], grouplist["delta_IminOinv"], grouplist["A"], grouplist["delta"], grouplist["Dstar"], grouplist["IminOinv"], grouplist["In"]
      );
      
      
    } else {
      
      Jac.submat(varPart_start,netPart_start,varPart_end,netPart_end ) =  d_sigma_omega_cpp( 
        grouplist["L"], grouplist["delta_IminOinv"], grouplist["A"], grouplist["delta"], grouplist["Dstar"]
      );
      
      Jac.submat(varPart_start,scalingPart_start,varPart_end,scalingPart_end ) =  d_sigma_delta_cpp( 
        grouplist["L"], grouplist["delta_IminOinv"], grouplist["In"],  grouplist["A"]
      );
      
      
    }
    
    
  } else if (y == "prec"){
    
    
    Jac.submat(varPart_start,varPart_start,varPart_end,varPart_end ) =  d_sigma_kappa_cpp( 
      grouplist["L"], grouplist["D"], grouplist["sigma"]
    );
    
  } else if (y == "cor"){
    
    int corPart_start = meanstructure * nMean_Thresh + nThresh;
    int corPart_end = meanstructure * nMean_Thresh + nThresh + nvar * (nvar-1) / 2 - 1;
    int sdPart_start = corPart_end + 1;
    int sdPart_end = sdPart_start + nvar - 1;
    
    Jac.submat(varPart_start,corPart_start,varPart_end,corPart_end ) =  d_sigma_rho_cpp( 
      grouplist["L"], grouplist["SD"], grouplist["A"], grouplist["Dstar"]
    );
    
    if (corinput == false){
      
      Jac.submat(varPart_start,sdPart_start,varPart_end,sdPart_end ) =  d_sigma_SD_cpp( 
        grouplist["L"], grouplist["SD_IplusRho"], grouplist["In"], grouplist["A"]
      );
      
    }
    
  }
  
  
  // If corinput, cut out the diagonal rows:
  if (corinput){
    
    
    arma::mat I = eye(nvar, nvar );
    arma::vec dummyvec = join_cols(
      zeros<vec>(nMean_Thresh),
      vech(I)
    );
    
    
    
    
    uvec remove = find(dummyvec > 0);
    // Rf_PrintValue(wrap(remove));
    Jac.shed_rows(remove);
    
    
    
  }
  
  
  if (meanstructure == false){
    if (noThresholds){
      Jac.shed_rows(0, nvar - 1);  
    }
    
    // FIXME: Add error for thresholds
  }
  
  // Return:
  return(Jac);
}
// 
// [[Rcpp::export]]
arma::mat d_phi_theta_varcov_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_varcov_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);
  
  return(res);
}

