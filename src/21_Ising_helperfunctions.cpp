// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
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
//
// The state space is enumerated with a mixed-radix counter (digit 0 fastest)
// over the ordered 'responses' vector. The same ordered response set is used
// for every variable. When 'responses' has length 2 this visits exactly the
// same states, in the same order, as the original binary r1/r2 counter, so
// Z and all expectations are bit-identical to the two-outcome implementation.
// [[Rcpp::export]]
Rcpp::List isingExpectation(
    arma::mat graph,
    arma::vec tau,
    double beta,
    arma::vec responses,
    double min_sum
){
  int nVar = graph.n_cols;
  int nResp = responses.n_elem;
  int i;

  // Mixed-radix index of the current state into 'responses', and the
  // corresponding state values:
  std::vector<int> idx(nVar, 0);
  arma::vec curstate(nVar);
  for (i = 0; i < nVar; i++){
    curstate[i] = responses(0);
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
  bool done = false;

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

    // Increment the mixed-radix counter (digit 0 fastest). Carry over any
    // digits that have reached the last response level, resetting them to
    // the first level:
    i = 0;
    while (i < nVar && idx[i] == nResp - 1){
      idx[i] = 0;
      curstate[i] = responses(0);
      i++;
    }
    if (i == nVar){
      // All digits at the last level: every state has been visited.
      done = true;
    } else {
      idx[i]++;
      curstate[i] = responses(idx[i]);
    }



  } while (!done);

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
  int nResp = responses.n_elem;
  int i;

  // Mixed-radix index of the current state into 'responses', and the
  // corresponding state values:
  std::vector<int> idx(nVar, 0);
  arma::vec curstate(nVar);
  for (i = 0; i < nVar; i++){
    curstate[i] = responses(0);
  }


  // Current Z:
  double Z = 0;

  bool done = false;

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

    // Increment the mixed-radix counter (digit 0 fastest):
    i = 0;
    while (i < nVar && idx[i] == nResp - 1){
      idx[i] = 0;
      curstate[i] = responses(0);
      i++;
    }
    if (i == nVar){
      done = true;
    } else {
      idx[i]++;
      curstate[i] = responses(idx[i]);
    }



  } while (!done);



  return(Z);
}
