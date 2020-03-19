// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "21_Ising_helperfunctions.h"
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// GROUP JACOBIAN FUNCTION //
// [[Rcpp::export]]
arma::mat ULS_Gauss_exphes_pergroup_cpp(
    const Rcpp::List& grouplist
){
  // FIXME: This is a copy of the fit function. Move to the prepare function perhaps.
  
  arma::mat WLS_W = grouplist["WLS.W"];
  
  std::string estimator = grouplist["estimator"];
  
  // If DWLS, only use the diagonal:
  if (estimator == "DWLS"){
    WLS_W = diagmat(WLS_W);
  }
  
  arma::mat Hes = 2 * WLS_W;
  
  return(Hes);
}


// full Jacobian function 
// [[Rcpp::export]]
arma::mat expected_hessian_ULS_Gaussian_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // JAcobian:
  Rcpp::List grouphessians(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat grouphes =  ((nPerGroup(i)+1) / nTotal) * ULS_Gauss_exphes_pergroup_cpp(groupmodels[i]);
    grouphessians[i]  = grouphes;
  }
  
  
  arma::mat res =  bdiag_psychonetrics(grouphessians);
  
  
  return(res);
}

