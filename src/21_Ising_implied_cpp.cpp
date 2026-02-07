// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "02_algebragelpers_kronecker.h"
#include "02_algebrahelpers_RcppHelpers.h"
#include "03_modelformation_formModelMatrices_cpp.h"
#include "03_modelformation_impliedcovstructures.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Core implementation that takes pre-formed model matrices
Rcpp::List implied_Ising_cpp_core(
    Rcpp::List x,
    const S4& model,
    bool all = false
){
  // Types:
  Rcpp::List types = model.slot("types");
  std::string beta_model = types["beta_model"];
  bool log_beta  = beta_model == "log_beta";

  int nGroup = x.length();
  int g, i, j;

  for (g=0; g<nGroup; g++){
    // Grouplist
    Rcpp::List grouplist = x[g];

   if (log_beta){
     arma::mat beta = grouplist["log_beta"];

     for (i=0;i<beta.n_rows;i++){
       for (j=0;j<beta.n_cols;j++){
           beta(i,j) = exp(beta(i,j));
       }
     }

     grouplist["beta"] = beta;
   } else {

     arma::mat log_beta = grouplist["beta"];

     for (i=0;i<log_beta.n_rows;i++){
       for (j=0;j<log_beta.n_cols;j++){
         log_beta(i,j) = log(log_beta(i,j));
       }
     }

     grouplist["log_beta"] = log_beta;


   }

   x[g] = grouplist;
  }


  return(x);
}

// Original version: forms matrices from the S4 model
// [[Rcpp::export]]
Rcpp::List implied_Ising_cpp(
    const S4& model,
    bool all = false
){
  // Form basic model matrices:
  Rcpp::List x = formModelMatrices_cpp(model);

  return implied_Ising_cpp_core(x, model, all);
}

