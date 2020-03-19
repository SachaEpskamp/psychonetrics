// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



arma::mat expected_hessian_Gaussian_group_meanPart_cpp(
  const arma::mat& kappa){
  
  arma::mat hes = 2 * kappa;
  
  return(hes);
}


arma::mat expected_hessian_Gaussian_group_varPart_cpp(
    const arma::sp_mat& D,
    const arma::mat& kappa){
  
  arma::mat hes = D.t() * kron(kappa, kappa) * D;
  
  return(hes);
}

// GROUP JACOBIAN FUNCTION //
// [[Rcpp::export]]
arma::mat expected_hessian_Gaussian_group_cpp(
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
  arma::mat Hes_sigma = expected_hessian_Gaussian_group_varPart_cpp(grouplist["D"],grouplist["kappa"]);
  
  arma::mat Out = Hes_sigma;
  
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
      
      Hes_sigma.shed_cols(remove);      
      Hes_sigma.shed_rows(remove);      
      
    } else {
      arma::mat I = eye(n, n );

      uvec remove = find(vech(I) > 0);

      Hes_sigma.shed_cols(remove);
      Hes_sigma.shed_rows(remove);
    }

  }
  
  if (meanstructure){
    arma::mat Hes_mean = expected_hessian_Gaussian_group_meanPart_cpp(grouplist["kappa"]);
    
    Rcpp::List Hessians(2);
    Hessians[0] = Hes_mean;
    Hessians[1] = Hes_sigma;
    Out = bdiag_psychonetrics(Hessians);
    
  } else {
    
    Out = Hes_sigma;
    
  }
  

  return(Out);
}


// full Jacobian function 
// [[Rcpp::export]]
arma::mat expected_hessian_Gaussian_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // JAcobian:
  Rcpp::List grouphessians(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat grouphes =  (nPerGroup(i) / nTotal) * expected_hessian_Gaussian_group_cpp(groupmodels[i]);
    grouphessians[i]  = grouphes;
  }
  
  
  arma::mat res =  bdiag_psychonetrics(grouphessians);
  
  
  return(res);
}

