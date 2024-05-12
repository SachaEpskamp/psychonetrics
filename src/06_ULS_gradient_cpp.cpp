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
arma::mat ULS_Gauss_gradient_pergroup_cpp(
    const Rcpp::List& grouplist
){
  // FIXME: This is a copy of the fit function. Move to the prepare function perhaps.
  
  arma::mat S = grouplist["S"];
  arma::vec means = grouplist["means"];
  
  arma::mat sigma = grouplist["sigma"];
  arma::vec mu = grouplist["mu"];
  
  arma::mat WLS_W = grouplist["WLS.W"];
  
  std::string estimator = grouplist["estimator"];
  
  bool meanstructure =  grouplist["meanstructure"];
  bool corinput =  grouplist["corinput"];
  
  // Empty vectors:
  arma::vec obs(1); // FIXME: Start with 1 element
  arma::vec imp(1); // FIXME: Start with 1 element
  
  // Integers I need
  int i;
  int nvar = mu.n_elem;
  
  // If DWLS, only use the diagonal:
  if (estimator == "DWLS"){
    WLS_W = diagmat(WLS_W);
  }
  
  // If no tau, do normal stuff:
  if (!grouplist.containsElementNamed("tau") || !grouplist.containsElementNamed("thresholds")){
    if (meanstructure){
      
      obs = join_cols(obs,means);
      imp = join_cols(imp,mu);
      
    } 
    
    
  } else {
    Rcpp::List threhsolds = grouplist["thresholds"];
    arma::mat tau = grouplist["tau"];
    
    for (i = 0; i < nvar; i++){
      if (is_finite(means(i))){
        // FIXME:  This is just silly ...
        arma::vec obselem(1);
        obselem(0) = means(i);
        arma::vec impelem(1);
        impelem(0) = mu(i);
        
        obs = join_cols(obs,obselem);
        imp = join_cols(imp,impelem);
      } else {
        arma::vec sampthresh = threhsolds[i];
        arma::vec modthresh = tau.submat(0,i,sampthresh.n_elem-1,i);
        
        obs = join_cols(obs,sampthresh);
        imp = join_cols(imp,modthresh);
      }
      
      
    }
    
  }
  
  
  // Add variances (if needed) and covariances:
  if (corinput){
    
    obs = join_cols(obs,vech(S, false));
    imp = join_cols(imp,vech(sigma, false));
    
  } else {
    
    obs = join_cols(obs,vech(S, true));
    imp = join_cols(imp,vech(sigma, true));
    
  }
  
  
  
  
  // FIXME: Cut out the first element
  obs = obs.subvec(1, obs.n_elem - 1);
  imp = imp.subvec(1, imp.n_elem - 1);
  
  // Now compute the gradient:
  arma::mat Jac = -2 * (obs - imp).t() * WLS_W;
  
  return(Jac);
}


// full Jacobian function 
// [[Rcpp::export]]
arma::mat ULS_gradient_Gauss_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // JAcobian:
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad =  ((nPerGroup(i)+1) / nTotal) * ULS_Gauss_gradient_pergroup_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  
  arma::mat res =  cbind_psychonetrics(groupgradients);
  
  
  return(res);
}

