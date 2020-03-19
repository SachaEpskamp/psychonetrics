// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "21_Ising_helperfunctions.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// GROUP FIT FUNCTION //
// [[Rcpp::export]]
double maxLikEstimator_Ising_group_cpp(
    const Rcpp::List& grouplist
){
  arma::vec thresholds = grouplist["tau"];
  arma::vec responses = grouplist["responses"];
  arma::mat squares = grouplist["squares"];
  arma::vec means = grouplist["means"];
  arma::mat graph = grouplist["omega"];
  double beta = grouplist["beta"];
  double nobs = grouplist["nobs"];
  
  
  
  // Compute Z:
  double Z = computeZ_cpp(graph, thresholds, beta, responses);
  
  // Compute summary statistics:
  // FIXME: Not nice, will make things double
  arma::vec v1 = means * nobs;
  arma::mat v2 = squares;
  
  arma::vec part1 = thresholds % v1;
  arma::mat part2 = graph % v2;
  
  
  double H =  (
    - sum(part1) - 
      sum(vech(part2, false))
  );
  
  // Fml
  double Fml = 2 * log(Z) + 2 * beta * H / nobs;
  
  
  // Return:
  return(Fml);
}


// full fit function 
// [[Rcpp::export]]
double maxLikEstimator_Ising_cpp(
    const Rcpp::List& prep
){
  
  Rcpp::List groupmodels = prep["groupModels"];
  int nGroup = groupmodels.length();
  arma::vec nPerGroup = prep["nPerGroup"];
  double nTotal = prep["nTotal"];
  
  // Result:
  double fit = 0;
  
  for (int i=0; i<nGroup;i++){
    fit += (nPerGroup(i) / nTotal) * maxLikEstimator_Ising_group_cpp(groupmodels[i]);
  }
  
  return(fit);
}

