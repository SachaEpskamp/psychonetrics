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

// [[Rcpp::export]]
Rcpp::List impliedcovstructures_cpp(
    Rcpp::List x,
    std::string name = "",
    std::string type = "cov",
    bool all = false
){
  
  int nGroup = x.length();
  int g;

  // Strings:
  std::string sigma;
  std::string omega;
  std::string delta;
  std::string kappa;
  std::string lowertri;
  std::string IminOinv;
  std::string delta_IminOinv;
  
  std::string rho;
  std::string SD;
  
  std::string IplusRho;
  std::string SD_IplusRho;
  
  if (name != ""){
    sigma = "sigma_" + name;
    omega = "omega_" + name;
    delta = "delta_" + name;
    kappa = "kappa_" + name;
    lowertri = "lowertri_" + name;
    IminOinv = "IminOinv_" + name;
    delta_IminOinv = "delta_IminOinv_" + name;
    
    rho = "rho_" + name;
    SD = "SD_" + name;
    
    IplusRho = "IplusRho_" + name;
    SD_IplusRho = "SD_IplusRho_" + name;
    
  } else {
    sigma = "sigma";
    omega = "omega";
    delta = "delta";
    kappa = "kappa";
    lowertri = "lowertri";
    IminOinv = "IminOinv";
    delta_IminOinv = "delta_IminOinv";
    
    rho = "rho";
    SD = "SD";
    IplusRho = "IplusRho";
    SD_IplusRho = "SD_IplusRho";
  }
  
  // For each group:
  
  for (g=0; g < nGroup; g++){
    
    bool proper = true;
    
    Rcpp::List grouplist = x[g]; // FIXME: This will copy the list.
    
    
    // Form the models:
    // Contemporaneous:
    if (type == "cov"){
      
      // Only need to do things if all = TRUE:
      if (all){
        arma::mat mat_sigma = grouplist[sigma];
        // bool anyNon0 = false;
        // for (int i=0;i<mat_sigma.n_rows;i++){
        //   for (int j=0; j<mat_sigma.n_cols;j++){
        //     if (mat_sigma(i,j) != 0){
        //       anyNon0 = true;
        //     }
        //   }
        // }
        
        if (anyNon0(mat_sigma)){
          arma::mat mat_kappa = solve_symmetric_cpp_matrixonly_withcheck(mat_sigma, proper);
          grouplist[kappa] = mat_kappa;
          grouplist[omega] = wi2net_cpp(mat_kappa);
          grouplist[rho] = cov2cor_cpp(mat_sigma);
          grouplist[SD] = SDmat(mat_sigma);
        }
      }
      
      
    }
    else if(type == "chol"){
      
      // form cov matrix:
      arma::mat mat_lowertri = grouplist[lowertri];
      arma::mat mat_sigma = mat_lowertri * mat_lowertri.t();
      grouplist[sigma] = mat_sigma;
      
      // Return precision and network if all = TRUE:
      if (all){
        if (anyNon0(mat_sigma)){
          arma::mat mat_kappa = solve_symmetric_cpp_matrixonly_withcheck(mat_sigma, proper);
          grouplist[kappa] = mat_kappa;
          grouplist[omega] = wi2net_cpp(mat_kappa);
          grouplist[rho] = cov2cor_cpp(mat_sigma);
          grouplist[SD] = SDmat(mat_sigma);
          
        }
        
      }
    } else if (type == "ggm"){
      
      arma::mat mat_omega = grouplist[omega];
      arma::mat I = eye(mat_omega.n_rows, mat_omega.n_cols);
      
      // First check if the delta Matrix is present (it is ignored when corinput = TRUE only, so don't need to know that that argument was used):
      arma::mat IminO = I - mat_omega;
      arma::mat IminO_inv = solve_symmetric_cpp_matrixonly_withcheck(IminO, proper);
      
      if (grouplist.containsElementNamed(delta.c_str()) == false){
        // FIXME: non positive definite matrices are even worse here... So I am trying to solve this with a spectral shift for now:
        // FIXME: use spectral shift in R code here ..
        grouplist[delta] = (arma::mat)invSDmat(IminO_inv);
      }
      
      arma::mat mat_delta = grouplist[delta];
      mat_delta = diagmat(mat_delta);
      
      // Compute sigma:
      arma::mat mat_sigma = mat_delta * IminO_inv * mat_delta;
      grouplist[sigma] = mat_sigma;
      
      // Stuff needed if all = TRUE:
      if (all){
        if (anyNon0(mat_sigma)){
          arma::mat mat_kappa = solve_symmetric_cpp_matrixonly_withcheck(mat_sigma, proper);
          grouplist[kappa] = mat_kappa;
          grouplist[rho] = cov2cor_cpp(mat_sigma);
          grouplist[SD] = SDmat(mat_sigma);
        }
        
      }
      
      // Extra matrix needed:
      if (!all){
        grouplist[IminOinv] = IminO_inv;
        grouplist[delta_IminOinv] =  mat_delta * IminO_inv;
      }
    } else if (type == "prec"){
      
      arma::mat mat_kappa = grouplist[kappa];
      arma::mat mat_sigma = solve_symmetric_cpp_matrixonly_withcheck(mat_kappa, proper);
      // Precision matrix
      grouplist[sigma] = mat_sigma;
      
      if (all) {
        grouplist[omega] = wi2net_cpp(mat_kappa);
        grouplist[rho] = cov2cor_cpp(mat_sigma);
        grouplist[SD] = SDmat(mat_sigma);
      }
    }  else if (type == "cor"){
      
      arma::mat mat_rho = grouplist[rho];
      int nvar = mat_rho.n_rows;
      
      // First check if the SD Matrix is present (it is ignored when corinput = TRUE only, so don't need to know that that argument was used):
      if (grouplist.containsElementNamed(SD.c_str()) == false){
        grouplist[SD] = eye(nvar, nvar);
      }
      
      arma::mat mat_SD = grouplist[SD];
      mat_SD = diagmat(mat_SD);
      arma::mat mat_sigma = mat_SD * (eye(nvar, nvar) + mat_rho) * mat_SD;
      grouplist[sigma] = mat_sigma;
      
      // Stuff needed if all = TRUE:
      if (all){
        if (anyNon0(mat_sigma)){
          
          arma::mat mat_kappa = solve_symmetric_cpp_matrixonly_withcheck(mat_sigma, proper);
          grouplist[kappa] = mat_kappa;
          grouplist[omega] = wi2net_cpp(mat_kappa);
        }
        
      }
      
      // Extra matrix needed:
      if (!all){
        arma::mat mat_IplusRho = eye(nvar,nvar) + mat_rho;;
        grouplist[IplusRho] = mat_IplusRho;
        grouplist[SD_IplusRho] = mat_SD * mat_IplusRho;
      }
      
    }
    
    
    // Return properness:
    grouplist["proper"] = proper;
    
    // Return to list:
    x[g] = grouplist;
  }
  
  
  
  return(x);
}


