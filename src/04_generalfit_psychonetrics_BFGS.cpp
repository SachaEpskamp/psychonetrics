// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <math.h>
#include "04_generalFit_implied_and_prepare.h"
#include "04_generalfit_fitfunction_cpp.h"
#include "04_generalFit_gradient_cpp.h"
#include "b_modelexpansion_updateModel_cpp.h"
#include "02_algebrahelpers_modelMatrix_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
S4 psychonetrics_BFGS(
    const S4& model,
    arma::mat Hstart
){
  int i;
  
  // Copy model:
  S4 newMod(model);
  
  // Parameters:
  Rcpp::List pars = newMod.slot("parameters");
  
  // Total pars:
  arma::vec est = pars["est"];
  arma::vec parnum = pars["par"];
  int nPar = max(parnum);
  
  // If 0, nothing to do!
  if (nPar == 0){
    newMod.slot("computed") = true;
    return(newMod);    
  }
  
  // Form the start vector:
  arma::vec x(nPar);
  
  for (i=0;i<nPar;i++){
    if (parnum(i) > 0){
      x(parnum(i)-1) = est(i);
    }
  }
  
  
  // Manual part:
  arma::sp_mat manualPart = Mmatrix_cpp_list(newMod.slot("parameters"));
  
  // Initial Hessian:
  arma::mat Hes = Hstart; //eye(nPar, nPar);
  
  // Setup things I need:
  arma::vec curGrad;
  arma::vec newGrad;
  arma::vec p;
  arma::vec y;
  double alpha_k ;
  double dummyfit ;
  double propfit;
  double propalpha;
  bool alphasearchcont;
  arma::vec denominator;
  
  // Converged:
  bool converged = false;
  
  int nIter = 0;
  // START OPTIMIZER:
  do{
    nIter++;
    // Current gradient:
    curGrad = psychonetrics_gradient_cpp(x, newMod);
    
    // Step 1: Find initial direction:
    p = solve(Hes, -curGrad);
    
    // STep 2:Find acceptable stepsize:
    // FIXME: This is a silly optimizer ....
    alpha_k = 1;
    
    // dummyfit = psychonetrics_fitfunction_cpp(x + alpha_k * p, newMod);
    // propfit;
    // propalpha;
    // 
    // alphasearchcont = true;
    // 
    // do{
    //   propalpha = alpha_k / 2;
    //   propfit = psychonetrics_fitfunction_cpp(x + propalpha * p, newMod);
    // 
    //   if (dummyfit < propfit){
    //     alphasearchcont = false;
    //   } else{
    //     alpha_k = propalpha;
    //     dummyfit = propfit;
    //   }
    // } while (alphasearchcont);
    
    
    // Step 3: update x:
    arma::vec s = alpha_k * p;
    x = x + s;
    
    // Step 4, set y:
    newGrad = psychonetrics_gradient_cpp(x, newMod);
    y = newGrad - curGrad;
    curGrad = newGrad;
    
    // Step 5: Update Hes estimate:
    denominator = (s.t() * Hes * s);
      Hes = Hes + (y * y.t()) / dot(y,s) - ( Hes * s * s.t() * Hes.t() )/ denominator(0);
      
    // Check if converged:
    converged = sum(abs(curGrad)) < 0.01;
    
  } while (!converged);
  
  
  // Update model:
  newMod = updateModel_cpp(x,newMod,false);
  
  // Set computed:
  newMod.slot("computed") = true;
  
  // Rf_PrintValue(wrap(nIter));
  
  // Return:
  return(newMod);
}