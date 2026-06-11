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

// Hamiltonian for the Spin distribution. Generalizes the classical Ising
// Hamiltonian with a per-node quadratic (Blume-Capel) term delta_i:
//
//   H(x) = - sum_i tau_i x_i + sum_i delta_i x_i^2 - sum_{i<j} omega_ij x_i x_j
//
// so that P(x) propto exp(-beta H(x)) =
//   exp( beta [ sum_i tau_i x_i - sum_i delta_i x_i^2 + sum_{i<j} omega_ij x_i x_j ] ).
// Setting delta = 0 recovers the classical Ising Hamiltonian exactly.
// [[Rcpp::export]]
double H(
    arma::vec state,
    arma::mat graph,
    arma::vec tau,
    arma::vec delta
){


  double Res = 0;
  int N = graph.n_rows;
  for (int i=0;i<N;i++)
  {
    Res -=  tau(i) * state(i);
    Res +=  delta(i) * state(i) * state(i);
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
    arma::vec delta,
    double beta
){
  return exp(-1.0 * beta * H(state,graph,tau,delta));
}


// Maximum of the log-potential -beta * H(state) over all (included) states.
// Used by the log-sum-exp accumulation below to keep the partition function
// and moment computations finite even for extreme parameter values (e.g.
// very large |tau|), where raw exp(-beta*H) would overflow to Inf or
// underflow to 0. Cheap first pass: only Hamiltonian evaluations.
double maxLogPot_Ising(
    const arma::mat& graph,
    const arma::vec& tau,
    const arma::vec& delta,
    double beta,
    const arma::vec& responses,
    double min_sum
){
  int nVar = graph.n_cols;
  int nResp = responses.n_elem;
  int i;

  std::vector<int> idx(nVar, 0);
  arma::vec curstate(nVar);
  for (i = 0; i < nVar; i++){
    curstate[i] = responses(0);
  }

  double M = R_NegInf;
  bool done = false;
  bool update = true;
  bool alwaysUpdate = min_sum == R_NegInf;

  do{
    if (!alwaysUpdate){
      update = sum(curstate) >= min_sum;
    }

    if (update){
      double logpot = -1.0 * beta * H(curstate, graph, tau, delta);
      if (logpot > M) M = logpot;
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

  return(M);
}

// Z and expected values:
//
// The state space is enumerated with a mixed-radix counter (digit 0 fastest)
// over the ordered 'responses' vector. The same ordered response set is used
// for every variable. When 'responses' has length 2 this visits exactly the
// same states, in the same order, as the original binary r1/r2 counter, so
// Z and all expectations are bit-identical to the two-outcome implementation.
//
// Accumulation is performed in the log domain (log-sum-exp): weights are
// w = exp(-beta*H - M) with M the maximum log-potential, so that logZ =
// M + log(sum(w)) stays finite even when raw exp(-beta*H) would overflow
// or underflow, and the moments are normalized by sum(w). For backward
// compatibility the raw Z = exp(logZ) is still returned (it may be Inf/0
// in extreme regimes; consumers should prefer logZ).
// [[Rcpp::export]]
Rcpp::List isingExpectation(
    arma::mat graph,
    arma::vec tau,
    arma::vec delta,
    double beta,
    arma::vec responses,
    double min_sum
){
  int nVar = graph.n_cols;
  int nResp = responses.n_elem;
  int i;

  // First pass: maximum log-potential over all included states:
  double M = maxLogPot_Ising(graph, tau, delta, beta, responses, min_sum);
  // Degenerate case (no included states): keep M = 0 so the loop below
  // accumulates zeros and logZ becomes -Inf, as before:
  bool anyState = (M != R_NegInf);
  if (!anyState) M = 0.0;

  // Mixed-radix index of the current state into 'responses', and the
  // corresponding state values:
  std::vector<int> idx(nVar, 0);
  arma::vec curstate(nVar);
  for (i = 0; i < nVar; i++){
    curstate[i] = responses(0);
  }


  // Current sum of shifted weights w = exp(-beta*H - M):
  double sumw = 0;

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
      curH = H(curstate, graph, tau, delta);

      curpot = exp(-1.0 * beta * curH - M);

      // Update sum of weights:
      sumw += curpot;

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

  // Log partition function:
  double logZ = anyState ? (M + log(sumw)) : R_NegInf;

  // Raw partition function (backward compatibility; may overflow to Inf or
  // underflow to 0 in extreme regimes):
  double Z = exp(logZ);

  // Normalize expectations by the sum of shifted weights:
  exp_v1 = exp_v1 / sumw;
  exp_v2 = exp_v2 / sumw;
  exp_H = exp_H / sumw;

  // Make return list:
  List L = List::create(Named("Z") = Z, Named("logZ") = logZ, Named("exp_v1") = exp_v1, Named("exp_v2") = exp_v2, Named("exp_H") = exp_H);


  return(L);
}


// Only logZ (log partition function), computed with log-sum-exp so the
// result stays finite even when raw exp(-beta*H) would overflow/underflow:
// [[Rcpp::export]]
double computeLogZ_cpp(
    arma::mat graph,
    arma::vec tau,
    arma::vec delta,
    double beta,
    arma::vec responses,
    double min_sum
){
  int nVar = graph.n_cols;
  int nResp = responses.n_elem;
  int i;

  // First pass: maximum log-potential over all included states:
  double M = maxLogPot_Ising(graph, tau, delta, beta, responses, min_sum);
  if (M == R_NegInf){
    // No included states:
    return(R_NegInf);
  }

  // Mixed-radix index of the current state into 'responses', and the
  // corresponding state values:
  std::vector<int> idx(nVar, 0);
  arma::vec curstate(nVar);
  for (i = 0; i < nVar; i++){
    curstate[i] = responses(0);
  }


  // Current sum of shifted weights:
  double sumw = 0;

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
      // Update sum of shifted weights:
      sumw += exp(-1.0 * beta * H(curstate, graph, tau, delta) - M);
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



  return(M + log(sumw));
}

// Only Z (kept for backward compatibility; may overflow to Inf or underflow
// to 0 in extreme regimes — prefer computeLogZ_cpp):
// [[Rcpp::export]]
double computeZ_cpp(
    arma::mat graph,
    arma::vec tau,
    arma::vec delta,
    double beta,
    arma::vec responses,
    double min_sum
){
  return(exp(computeLogZ_cpp(graph, tau, delta, beta, responses, min_sum)));
}
