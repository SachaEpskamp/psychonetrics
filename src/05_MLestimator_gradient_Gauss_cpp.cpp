// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



arma::mat jacobian_gaussian_group_sigmaVersion_meanPart_cpp(
    const arma::mat& sigma,
    const arma::vec& mu,
    const arma::vec& means,
    const arma::mat& kappa){
  
  arma::mat grad_mean = -2 * (means - mu).t() * kappa;
  return(grad_mean);
}


arma::mat jacobian_gaussian_group_sigmaVersion_sigmaPart_cpp(
    const arma::mat& S,
    const arma::vec& means,
    const arma::vec& mu,
    const arma::mat& sigma,
    const arma::sp_mat& D,
    const arma::mat& kappa){
  
  // int n = S.n_cols;
  
  arma::mat innermat = S + (means - mu) * (means - mu).t() - sigma;
  arma::mat grad_sigma = (-D.t() * vectorise(kappa * innermat * kappa)).t();
  
  return(grad_sigma);
}

// GROUP JACOBIAN FUNCTION //
// [[Rcpp::export]]
arma::mat jacobian_gaussian_sigma_group_cpp(
    const Rcpp::List& grouplist
){
  // Stuff I need now:
  bool corinput = false;
  if (grouplist.containsElementNamed("corinput")){
    corinput = grouplist["corinput"];
  }
  
  bool meanstructure = true;
  if (grouplist.containsElementNamed("meanstructure")){
    meanstructure = grouplist["meanstructure"];
  }

  
  
  arma::mat sigma = grouplist["sigma"];
  int n = sigma.n_rows;
  
  // Sigma part:
  arma::mat grad_sigma = jacobian_gaussian_group_sigmaVersion_sigmaPart_cpp(grouplist["S"],grouplist["means"],grouplist["mu"],grouplist["sigma"],grouplist["D"],grouplist["kappa"]);
  
  arma::mat Out = grad_sigma;
  
  // FIXME: Nicer to not bother with computing the parts not needed
  // Cut out variances if needed:
  if (corinput){
    if (meanstructure){
      arma::mat I = eye(n, n );
      arma::vec dummyvec = join_cols(
        zeros<vec>(n),
        vech(I)
      );
      
      uvec remove = find(dummyvec > 0);
      
      grad_sigma.shed_cols(remove);      
    } else {
      arma::mat I = eye(n, n );

      uvec remove = find(vech(I) > 0);

      grad_sigma.shed_cols(remove);
    }

  }
  
  if (meanstructure){
    arma::mat grad_mean = jacobian_gaussian_group_sigmaVersion_meanPart_cpp(grouplist["sigma"],grouplist["mu"],grouplist["means"],grouplist["kappa"]);
    Out = join_rows(grad_mean, grad_sigma);
  } else {
    Out = grad_sigma;
  }
  

  return(Out);
}


// full Jacobian function 
// [[Rcpp::export]]
arma::mat jacobian_gaussian_sigma_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // JAcobian:
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad =  (nPerGroup(i) / nTotal) * jacobian_gaussian_sigma_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  
  arma::mat res =  cbind_psychonetrics(groupgradients);
  
  
  return(res);
}

