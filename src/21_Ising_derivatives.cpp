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
  // Distribution parameters phi = (tau, omega[lower.tri], delta, beta).
  // The map from model parameters theta to phi is the identity (the Ising and
  // BlumeCapel models map their parameters directly onto the Spin distribution);
  // for the Ising model the delta entries are fixed to zero in the parameter
  // table, so the corresponding identity columns are simply never used.
  int ntotal = n +  n * (n-1) / 2 + n + 1;

  arma::mat Jac = eye(ntotal,ntotal);
  
  
  // Types:
  std::string beta_model = grouplist["beta_model"];
  bool log_beta  = beta_model == "log_beta";

  if (log_beta){
    // Fix last derivative if needed:
    arma::mat log_beta = grouplist["log_beta"];
    Jac(ntotal-1, ntotal-1) = exp(log_beta(0,0));
  }
  
  
  
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

