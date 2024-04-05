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
  return exp(-1.0 * beta * H(state,graph,tau));
}


// Z and expected values:
// [[Rcpp::export]]
Rcpp::List isingExpectation(
    arma::mat graph,
    arma::vec tau,
    double beta,
    arma::vec responses,
    double min_sum
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
  
  // Current exp_v1
  arma::vec exp_v1 = zeros(nVar, 1);
  
  // Current exp_v2"
  arma::mat exp_v2 = zeros(nVar, nVar);
  
  // Current expected Hamiltonian:
  double exp_H = 0;
  
  // Dummy for potential:
  double curpot, curH;
  
  // Current bool to stop:
  bool all_r2 = false;
  
  // Do we need to update? for min_sum cutoff model
  bool update = true;
  bool alwaysUpdate = min_sum == R_NegInf;
  
  // Repeat:
  do{
    
    // First check if we need to update:
    if (!alwaysUpdate){
      update = sum(curstate) >= min_sum; 
    }
    
    
    // Update the expected values:
    if (update){
      curH = H(curstate, graph, tau);
      
      curpot = exp(-1.0 * beta * curH);
      
      // Update Z:
      Z += curpot;
      
      // Update exp_v1:
      exp_v1 += curpot * curstate;
      
      // Update exp_v2:
      exp_v2 += curpot * curstate * curstate.t();
      
      // Update exp_H:
      exp_H += curpot * curH;
    }
    
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
  
  // Normalize expectations:
  exp_v1 = exp_v1 / Z;
  exp_v2 = exp_v2 / Z;
  exp_H = exp_H / Z;
  
  // Make return list:
  List L = List::create(Named("Z") = Z , Named("exp_v1") = exp_v1, Named("exp_v2") = exp_v2, Named("exp_H") = exp_H);
  
  
  return(L);
}


// Only Z:
// [[Rcpp::export]]
double computeZ_cpp(
    arma::mat graph,
    arma::vec tau,
    double beta,
    arma::vec responses,
    double min_sum
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
  
  // Do we need to update? for min_sum cutoff model
  bool update = true;
  bool alwaysUpdate = min_sum == R_NegInf;
  
  
  // Repeat:
  do{
    if (!alwaysUpdate){
      update = sum(curstate) >= min_sum; 
    }
    
    if (update){
      // Update Z:
      Z += Pot(curstate, graph, tau, beta); 
    }
    
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
