// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Hammiltonian (copied from IsingSampler):
// [[Rcpp::export]]
double H(
    arma::vec state,
    arma::mat graph,
    arma::vec tau
){
  
  
  double Res = 0;
  int N = graph.n_rows;
  for (int i=0;i<N;i++)
  {
    Res -=  tau(i) * state(i);
    for (int j=i;j<N;j++)
    {
      if (j!=i) Res -= graph(i,j) * state(i) * state(j);
    }
  }
  return(Res);
  
}

// Potential (copied from IsingSampler):
// [[Rcpp::export]]
double Pot(
    arma::vec state,
    arma::mat graph,
    arma::vec tau,
    double beta
){
  exp(-1.0 * beta * H(state,graph,tau));
}


// [[Rcpp::export]]
double computeZ_cpp(
   arma::mat graph,
   arma::vec tau,
   double beta,
   arma::vec responses 
){
  int nVar = graph.n_cols;
  double r1 = responses(0);
  double r2 = responses(1);
  int i, j;
  
  // Current score pattern:
  arma::vec curstate(nVar);
  
  for (i = 0; i < nVar; i++){
    curstate[i] = r1;
  }
  
  
  // Current Z:
  double Z = 0;
  
  bool all_r2 = false;
  
  // Repeat:
  do{
    // Update Z:
    Z += Pot(curstate, graph, tau, beta);
    
    // Check if leftmost state is not r2, if not, make r2:
    if (curstate(0) == r2){
      
      // Search for the first element that is not r2:
      all_r2 = true;
      for (i=0; i<nVar && all_r2; i++){
        
        if (curstate(i) == r1){
          // Set all_r2 to false:
          all_r2 = false;
          
          // Make this element r2:
          curstate(i) = r2;
          
          // Make all elements left of this r1:
          for (j=0;j<i;j++){
            curstate(j) = r1;
          }
          
          break;
        }
        
      }
      
    } else {
      curstate(0) = r2;
    }
    
    

  } while (!all_r2);

  
  
  return(Z);
}