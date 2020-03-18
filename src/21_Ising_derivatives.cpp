// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// FULL GROUP JACOBIAN ///
// [[Rcpp::export]]
arma::mat d_phi_theta_Ising_group_cpp(
    const Rcpp::List& grouplist
){
  arma::mat omega = grouplist["omega"];
  int n = omega.n_rows;
  int ntotal = n +  n * (n-1) / 2 + 1;

  arma::mat Jac = eye(ntotal,ntotal);

  // Return:
  return(Jac);
}


// [[Rcpp::export]]
arma::mat d_phi_theta_Ising_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();

  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad = d_phi_theta_Ising_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  arma::mat res =  bdiag_psychonetrics(groupgradients);

  return(res);
}

