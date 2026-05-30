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
arma::mat jacobian_Ising_group_cpp(
    const Rcpp::List& grouplist
){
  arma::vec thresholds = grouplist["tau"];
  arma::vec delta = grouplist["delta"];
  arma::vec responses = grouplist["responses"];
  arma::mat squares = grouplist["squares"];
  arma::vec means = grouplist["means"];
  arma::mat graph = grouplist["omega"];
  double beta = grouplist["beta"];
  double nobs = grouplist["nobs"];

  // double Z = grouplist["Z"];
  double exp_H = grouplist["exp_H"];
  arma::vec exp_v1  = grouplist["exp_v1"];
  arma::mat exp_v2  = grouplist["exp_v2"];

  arma::vec v1 = means * nobs;
  arma::mat v2 = squares;

  // Hamiltonian (includes the Blume-Capel quadratic term):
  arma::vec part1 = thresholds % v1;
  arma::mat part2 = graph % v2;


  double H =  (
    - sum(part1) +
      sum(delta % v2.diag()) -
      sum(vech(part2, false))
  );


  // Thresholds gradient:
  arma::vec threshGrad = (
    2 * beta * exp_v1 - 2 * beta * v1 / nobs
  );


  // Network: gradient
  arma::mat graphmatrixgrad = 2 * beta * exp_v2 - 2 * beta * v2 / nobs;
  arma::vec graphGrad = vech(graphmatrixgrad, false);


  // delta gradient: d Fml / d delta_i = 2 beta ( obs E[x_i^2] - model E[x_i^2] ).
  // delta enters the exponent with a minus sign, so this has the opposite sign
  // pattern to the threshold gradient. Diagonal of 'squares' holds sum_obs x_i^2.
  arma::vec deltaGrad = (
    2 * beta * v2.diag() / nobs - 2 * beta * exp_v2.diag()
  );


  // beta gradient
  arma::vec betaGrad(1);
  betaGrad(0) = (
    2*(H/nobs - exp_H)
  );

  // Final gradient (order: tau, omega[lower.tri], delta, beta):
  arma::mat grad = join_rows(join_rows(threshGrad.t(), graphGrad.t()), join_rows(deltaGrad.t(), betaGrad));

    // Return:
    return(grad);
}


// full Jacobian function 
// [[Rcpp::export]]
arma::mat jacobian_Ising_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // JAcobian:
  Rcpp::List groupgradients(nGroup);
  
  for (int i=0; i<nGroup;i++){
    arma::mat groupgrad =  (nPerGroup(i) / nTotal) * jacobian_Ising_group_cpp(groupmodels[i]);
    groupgradients[i]  = groupgrad;
  }
  
  
  arma::mat res =  cbind_psychonetrics(groupgradients);
  
  
  return(res);
}

